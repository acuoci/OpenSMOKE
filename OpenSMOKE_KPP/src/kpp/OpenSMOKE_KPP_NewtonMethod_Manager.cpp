/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "OpenSMOKE_KPP_NewtonMethod_Manager.h"
#include "OpenSMOKE_KPP_ReactorNetwork.h"
#include <mpi.h>
#include <petscvec.h>

const KPP_NonLinearSystem_Method OpenSMOKE_KPP_NewtonMethod_Manager::default_method_ = KPP_NLS_KELLEY;
const int		OpenSMOKE_KPP_NewtonMethod_Manager::default_maxit_	= 40;
const int		OpenSMOKE_KPP_NewtonMethod_Manager::default_maxarm_ = 20;
const int		OpenSMOKE_KPP_NewtonMethod_Manager::default_shamanskii = 3;	// TODO
const double	OpenSMOKE_KPP_NewtonMethod_Manager::default_alpha_	= 1.e-4;
const double	OpenSMOKE_KPP_NewtonMethod_Manager::default_relativeTolerance_ = 1.e-8;
const double	OpenSMOKE_KPP_NewtonMethod_Manager::default_absoluteTolerance_ = 1.e-12;
const double	OpenSMOKE_KPP_NewtonMethod_Manager::default_safetyReductionCoefficient_ = 0.99;

void OpenSMOKE_KPP_NewtonMethod_Manager::ErrorMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class: OpenSMOKE_KPP_NewtonMethod_Manager"	<< std::endl;
    std::cout << "Error: " << message_				<< std::endl;
    std::cout << "Press enter to continue... "		<< std::endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_NewtonMethod_Manager::WarningMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class:	  OpenSMOKE_KPP_NewtonMethod_Manager"		<< std::endl;
    std::cout << "Warning: " << message_					<< std::endl;
	std::cout << std::endl;
}

OpenSMOKE_KPP_NewtonMethod_Manager::OpenSMOKE_KPP_NewtonMethod_Manager()
{
}

OpenSMOKE_KPP_NewtonMethod_Manager::OpenSMOKE_KPP_NewtonMethod_Manager(BzzVector& x_, void (*ptResiduals_)(BzzVector& x, Vec &f), int (*ptSolution_)(BzzVector& x, BzzVector& dir, const bool jacobianFlag), void (*ptAnalysis_)(), OpenSMOKE_KPP_ReactorNetwork& network)
{
	Setup(x_, ptResiduals_, ptSolution_, ptAnalysis_, network);
}

void OpenSMOKE_KPP_NewtonMethod_Manager::Setup(BzzVector& x_, void (*ptResiduals_)(BzzVector& x, Vec& f), int (*ptSolution_)(BzzVector& x, BzzVector& dir, const bool jacobianFlag), void (*ptAnalysis_)(), OpenSMOKE_KPP_ReactorNetwork& network)
{
	procrank_ = MPI::COMM_WORLD.Get_rank();
	nprocs_ = MPI::COMM_WORLD.Get_size();
	numworkers_ = nprocs_ - 1;

	MASTER = 0;
	FROM_MASTER = 1;
	FROM_WORKER = 2;

	n = network.NumberOfEquations();

	xt_ptr = new double[n + 1];
	ptResiduals = ptResiduals_;
	ptSolution  = ptSolution_;
	ptAnalysis  = ptAnalysis_;

	ChangeDimensions(n, &x);
	ChangeDimensions(n, &xt);
	
	if(procrank_ == 0)
	{
	    ChangeDimensions(n, &xOld);
	    ChangeDimensions(n, &step);
	    ChangeDimensions(n, &direction);
	}


	VecCreate(PETSC_COMM_WORLD, &f0);
	f0size_ = network.LocalNumberOfReactors()[procrank_] * network.NumberOfSpecies();
	VecSetSizes(f0, f0size_, PETSC_DECIDE);
	VecSetFromOptions(f0);

	VecCreate(PETSC_COMM_WORLD, &fOld);
	fOldsize_ = network.LocalNumberOfReactors()[procrank_] * network.NumberOfSpecies();
	VecSetSizes(fOld, fOldsize_, PETSC_DECIDE);
	VecSetFromOptions(fOld);

	VecCreate(PETSC_COMM_WORLD, &ft);
	ftsize_ = network.LocalNumberOfReactors()[procrank_] * network.NumberOfSpecies();
	VecSetSizes(ft, ftsize_, PETSC_DECIDE);
	VecSetFromOptions(ft);

	ratioNormInf	= new RingVector<double>(5);
	ratioNorm2		= new RingVector<double>(5);

	// Parameters
	absoluteTolerance_ = default_absoluteTolerance_;
	relativeTolerance_ = default_relativeTolerance_;
	maxarm_ =  default_maxarm_;
	maxit_  =  default_maxit_;
	alpha_  =  default_alpha_;	

	SetMethod(default_method_);

	// Counters
	itc    =  0;
	itsham = isham_;

	// Evaluate residuals at first iterate
	// Evaluate residuals at first iterate
	if(procrank_ == 0)
	    x = x_;

	double *x_ptr = new double[network.NumberOfEquations() + 1];
	if(procrank_ == 0)
	    for(int i = 1; i <= network.NumberOfEquations(); i++)
		x_ptr[i] = x[i];

	MPI::COMM_WORLD.Bcast(&x_ptr[1], network.NumberOfEquations(), MPI::DOUBLE, MASTER);

	if(procrank_ > 0)
	    for(int i = 1; i <= network.NumberOfEquations(); i++)
		x[i] = x_ptr[i];


	ptResiduals(x, f0);
	VecNorm(f0, NORM_2, &fnrm);

	fnrm0 = 1.;

	// Print data on video
	verbose = true;

	delete [] x_ptr;
}

void OpenSMOKE_KPP_NewtonMethod_Manager::SetMethod(const KPP_NonLinearSystem_Method method, const int m)
{
	method_ = method;

	if (method_ == KPP_NLS_KELLEY)
	{
		isham_  = -1;
		rsham_  =  0.50;
	}
	else if (method_ == KPP_NLS_NEWTON)
	{
		isham_  = 1;
		rsham_  = 0.;
	}
	else if (method_ == KPP_NLS_CHORD)
	{
		isham_  = -1;
		rsham_  =  1.;
	}
	else if (method_ == KPP_NLS_SHAMANSKII)
	{
		isham_  = m;
		rsham_  = 1.;
	}
}

void OpenSMOKE_KPP_NewtonMethod_Manager::SetRelativeTolerance(const double relativeTolerance)
{
	relativeTolerance_ = relativeTolerance;
}

void OpenSMOKE_KPP_NewtonMethod_Manager::SetAbsoluteTolerance(const double absoluteTolerance)
{
	absoluteTolerance_ = absoluteTolerance;
}

void OpenSMOKE_KPP_NewtonMethod_Manager::SetMaximumNumberOfIterations(const int maximumIterations)
{
	maxit_ = maximumIterations;
}

void OpenSMOKE_KPP_NewtonMethod_Manager::SetMaximumNumberOfArmijioIterations(const int maximumArmijioIterations)
{
	maxarm_ = maximumArmijioIterations;
}

void OpenSMOKE_KPP_NewtonMethod_Manager::SetAlpha(const double alpha)
{
	alpha_ = alpha;
}

// -2: Complete Armijio failure
// -1: Error in solving the linear system
//  0: Solution succesfully reached
//  1: Maximum number of iteration (but the solution is feasible)

int OpenSMOKE_KPP_NewtonMethod_Manager::Solve(BzzVector& solution)
{
	// Compute the stop tolerance
	stop_tol = absoluteTolerance_ + relativeTolerance_*fnrm;

	int flag = 0;
	int arm_flag=0;
	int jac_age=0;

	while (itc<maxit_ && fnrm>stop_tol)
	{
		if(procrank_ == 0)
		    startTime = MPI::Wtime();

		double ratio = fnrm/fnrm0;
		fnrm0 = fnrm;
		itc++;

		if (verbose == true && procrank_ == 0)
		{
			std::cout << " * Iteration:           " << itc		<< std::endl;
			std::cout << "   Current norm ||ft||: " << fnrm		<< std::endl;
			std::cout << "   Requested norm:      " << stop_tol	<< std::endl;
			std::cout << "   Current ratio:       " << ratio		<< std::endl;
		}

		// Updating Jacobian and Evaluating Residuals
		if (itc==1 || ratio>rsham_ || itsham==0 || arm_flag==1)
		{
			if (verbose == true && procrank_ == 0)
			{
					 if (itc == 1)		std::cout << "   Jacobian is updated (First Iteration)"		<< std::endl;
				else if (ratio>rsham_)	std::cout << "   Jacobian is updated (Ratio too large)"		<< std::endl;
				else if (itsham==0)		std::cout << "   Jacobian is updated (Requested by the user)"	<< std::endl;
				else if (arm_flag==1)	std::cout << "   Jacobian is updated (After Armijio Failure)"	<< std::endl;
			}

			itsham = isham_;
			jac_age = -1;

			
			if(procrank_ == 0)
			    startTime = MPI::Wtime();
		

			flag = ptSolution(x, direction, 1);

			if (flag<0) return -1;

			if(procrank_ == 0)
			    endTime = MPI::Wtime();
			
			if (verbose == true && procrank_ == 0) 
				std::cout << "   Elapsed time: " << endTime-startTime << std::endl;
		}
		// Evaluating Residuals (without updating Jacobian)
		else 
		{
			if (verbose == true && procrank_ == 0)
				std::cout << "   Jacobian is not updated" << std::endl;

			if(procrank_ == 0)
			    startTime = MPI::Wtime();

			flag = ptSolution(x, direction, 0);

			if (flag<0) return -1;

			double endTime = MPI::Wtime();
			
			if (verbose == true && procrank_ == 0)
				std::cout << "   Elapsed time: " << endTime-startTime << std::endl;
		}

		itsham--;
		jac_age++;

		if(procrank_ == 0)
		    xOld = x;

		VecCopy(f0, fOld);


		// Armijio Search
		arm_flag = ArmijoStep(x,f0,direction);

		if(arm_flag == 2)
		{
		    flag = 2;
		    convergence_ = false;
		    break;
		}

		// Failure of Armjio rule
		if (arm_flag == 1)
		{
			// If the Jacobian was old, the procedure can be repeated after updating the jacobian
			if (jac_age > 0)
			{
				if (verbose == true && procrank_ == 0)
					std::cout << "Armjio failure: Jacobian is recalculated" << std::endl;

				double *xOld_ptr = new double[n + 1];
				if(procrank_ == 0)
				{
				    x  = xOld;
				    for(int i = 1; i <= n; i++)
				        xOld_ptr[i] = xOld[i];
				}

				MPI::COMM_WORLD.Bcast(&xOld_ptr[1], n, MPI::DOUBLE, MASTER);

				if(procrank_ > 0)
				{
				    for(int i = 1; i <= n; i++)
					x[i] = xOld_ptr[i];
				}

				delete [] xOld_ptr;


				VecCopy(fOld, f0);
			}
			// Else the Newton method failed
			else
			{
				if (verbose == true && procrank_ == 0)
					std::cout << "Complete Armjio failure" << std::endl;

				if(procrank_ == 0)
				    solution = xOld;

				convergence_ = false;


				VecCopy(fOld, f0);
				VecNorm(f0, NORM_2, &fnrm);
				return -2;
			}
		}

		VecNorm(f0, NORM_2, &fnrm);		
		VecNorm(f0, NORM_INFINITY, &f0MaxAbs);		

		ratio = fnrm/fnrm0;

		ratioNormInf->Append(f0MaxAbs);
		ratioNorm2->Append(fnrm);

		double endTime = MPI::Wtime();

		if (verbose == true && procrank_ == 0)
		    std::cout << "CPU Time: " << endTime - startTime << std::endl;

		if (arm_flag == 0)
		    ptAnalysis();

		if ( ((ratioNormInf->MeanRatios() > 0.99) && (ratioNorm2->MeanRatios() > 0.99) && itc>=3 ))
		{
			flag = 2;
			break;
		}

	//	if ( (ratioNormInf->MeanRatios() < 0.70) && (ratioNorm2->MeanRatios() < 0.70) && itc>=3 )
	//	{
	//		flag = 3;
	//		break;
	//	}
	}

	if(procrank_ == 0)
	    solution = x;

	if (itc>maxit_)	flag=1;

	if (fnrm < stop_tol)    convergence_ =  true;

	else convergence_ = false;

	return flag;
}

int OpenSMOKE_KPP_NewtonMethod_Manager::ArmijoStep(BzzVector& x, Vec& f0, const BzzVector& direction)
{
	
	if (verbose == true && procrank_ == 0)
	{
		std::cout << " ******************************************************************** " << std::endl;
		std::cout << "                           Armijio Search                             " << std::endl;
		std::cout << " ******************************************************************** " << std::endl;
	}

	for(int i = 0; i <= n; i++)
	    xt_ptr[i] = 0.;

	int arm_flag = 0;
	int iarm	 = 0;

	double lambda		= 1.;
	double lamm		= 1.;
	double lamc   	  = lambda;
	
	if(procrank_ == 0)
	{
	    Product(lambda, direction, &step);
	    Sum(x,step, &xt);
	    for(int i = 1; i <= n; i++)
		xt_ptr[i] = xt[i];
	}

	MPI::COMM_WORLD.Bcast(&xt_ptr[1], n, MPI::DOUBLE, MASTER);

	if(procrank_ > 0)
	{
	    for(int i = 1; i <= n; i++)
		xt[i] = xt_ptr[i];
	}
	    

	ptResiduals(xt, ft);

	double nft, nf0, ff0, ffc, ffm = 0.;


	VecNorm(ft, NORM_2, &nft);

	VecNorm(f0, NORM_2, &nf0);
	ff0 = nf0*nf0; 
	ffc = nft*nft; 
	ffm = nft*nft;

	if (verbose == true && procrank_ == 0)
	{
		std::cout << " * Starting norm      ||f0||: " << nf0 << std::endl;
		std::cout << " * Current norm       ||ft||: " << nft << std::endl;
		std::cout << " * Objective function   F0:   " << ff0 << std::endl;
		std::cout << " * Objective function   Fc:   " << ffc << std::endl;
		std::cout << " * Objective function   Fm:   " << ffm << std::endl;
	}

	while ( nft >= ((1.-alpha_*lambda)*nf0) )
	{
		if (iarm == 0)	lambda *= 0.50;
		else		lambda  = parabolicModel(lamc, lamm, ff0, ffc, ffm);

		if(procrank_ == 0)
		{
	    	    Product(lambda, direction, &step);
	    	    Sum(x,step, &xt);
	    	    for(int i = 1; i <= n; i++)
			xt_ptr[i] = xt[i];
		}

		MPI::COMM_WORLD.Barrier();

//		if(procrank_ == 0)	cout << "Qui arrivo" << endl;	getchar();

		MPI::COMM_WORLD.Bcast(&xt_ptr[1], n, MPI::DOUBLE, MASTER);

//		if(procrank_ == 0)	cout << "Qui arrivo" << endl;	getchar();

		if(procrank_ > 0)
		{
	    	    for(int i = 1; i <= n; i++)
			xt[i] = xt_ptr[i];
		}

		if(procrank_ == 0)
		{
        	    lamm = lamc;
        	    lamc = lambda;
		}

//		if(procrank_ == 0)	cout << "Qui arrivo" << endl;	getchar();

		MPI::COMM_WORLD.Bcast(&lamm, 1, MPI::DOUBLE, MASTER);
		MPI::COMM_WORLD.Bcast(&lamc, 1, MPI::DOUBLE, MASTER);
		MPI::COMM_WORLD.Bcast(&lambda, 1, MPI::DOUBLE, MASTER);


		ptResiduals(xt, ft);

		VecNorm(ft, NORM_2, &nft);
        	ffm = ffc;
        	ffc = nft*nft;

        	iarm++;

		if (verbose == true && procrank_ == 0)
		{
			std::cout << " * Iteration:               " << iarm << std::endl;
			std::cout << "   Lambda:                  " << lambda << std::endl;
			std::cout << "   Current norm     ||ft||: " << nft << std::endl;
			std::cout << "   Objective function Fc:   " << ffc << std::endl;
			std::cout << "   Objective function Fm:   " << ffm << std::endl;
			std::cout << "   Required reduction:      " << (1.-alpha_*lambda) << std::endl;
		}

		if (iarm > maxarm_)
		{
		    if (verbose == true && procrank_ == 0)
			std::cout << "Armijo failure, too many reductions" << std::endl;

            	    arm_flag = 1;
            	    return arm_flag;
		}
	}

	if (verbose == true && procrank_ == 0)
	{
		std::cout << std::endl;
		std::cout << " ********************** Armijio Search Succeeded! ******************* " << std::endl;
		std::cout << "  Norm: " << nf0 <<" ==> " << nft << "  Reduction: " << nft/nf0 << std::endl;
		std::cout << std::endl;
	}

	if(nft/nf0 > 0.999)
	{
	    arm_flag = 2;
	    return arm_flag;
	}

	x = xt;
	VecCopy(ft, f0);

	return arm_flag;
}
