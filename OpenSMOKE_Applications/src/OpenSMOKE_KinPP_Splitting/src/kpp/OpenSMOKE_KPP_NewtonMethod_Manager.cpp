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

const KPP_NonLinearSystem_Method OpenSMOKE_KPP_NewtonMethod_Manager::default_method_ = KPP_NLS_KELLEY;
const int		OpenSMOKE_KPP_NewtonMethod_Manager::default_maxit_	= 40;
const int		OpenSMOKE_KPP_NewtonMethod_Manager::default_maxarm_ = 20;
const int		OpenSMOKE_KPP_NewtonMethod_Manager::default_shamanskii = 3;	// TODO
const double	OpenSMOKE_KPP_NewtonMethod_Manager::default_alpha_	= 1.e-4;
const double	OpenSMOKE_KPP_NewtonMethod_Manager::default_relativeTolerance_ = 1.e-8;
const double	OpenSMOKE_KPP_NewtonMethod_Manager::default_absoluteTolerance_ = 1.e-12;
const double	OpenSMOKE_KPP_NewtonMethod_Manager::default_safetyReductionCoefficient_ = 0.99;

void OpenSMOKE_KPP_NewtonMethod_Manager::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_NewtonMethod_Manager"	<< endl;
    cout << "Error: " << message_				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_NewtonMethod_Manager::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_NewtonMethod_Manager"		<< endl;
    cout << "Warning: " << message_					<< endl;
	cout << endl;
}

OpenSMOKE_KPP_NewtonMethod_Manager::OpenSMOKE_KPP_NewtonMethod_Manager()
{
}

OpenSMOKE_KPP_NewtonMethod_Manager::OpenSMOKE_KPP_NewtonMethod_Manager(BzzVector& x_, void (*ptResiduals_)(BzzVector& x, BzzVector &f), int (*ptSolution_)(BzzVector& x, BzzVector& dir, const bool jacobianFlag), void (*ptAnalysis_)())
{
	Setup(x_, ptResiduals_, ptSolution_, ptAnalysis_);
}

void OpenSMOKE_KPP_NewtonMethod_Manager::Setup(BzzVector& x_, void (*ptResiduals_)(BzzVector& x, BzzVector &f), int (*ptSolution_)(BzzVector& x, BzzVector& dir, const bool jacobianFlag), void (*ptAnalysis_)())
{
	n = x_.Size();
	ptResiduals = ptResiduals_;
	ptSolution  = ptSolution_;
	ptAnalysis  = ptAnalysis_;

	ChangeDimensions(n, &xt);
	ChangeDimensions(n, &ft);
	ChangeDimensions(n, &xOld);
	ChangeDimensions(n, &fOld);
	ChangeDimensions(n, &x);
	ChangeDimensions(n, &f0);
	ChangeDimensions(n, &step);
	ChangeDimensions(n, &direction);

	ratioNormInf	= new RingVector<double>(3);
	ratioNorm2		= new RingVector<double>(3);

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
	x = x_;
	ptResiduals(x, f0);
	fnrm  = f0.Norm2();
	fnrm0 = 1.;

	// Print data on video
	verbose = true;
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
	double stop_tol = absoluteTolerance_ + relativeTolerance_*fnrm;

	int flag = 0;
	int arm_flag=0;
	int jac_age=0;

	while (itc<maxit_ && fnrm>stop_tol)
	{
		double startTime = BzzGetCpuTime();

		double ratio = fnrm/fnrm0;
		fnrm0 = fnrm;
		itc++;

		if (verbose == true)
		{
			cout << " * Iteration:           " << itc		<< endl;
			cout << "   Current norm ||ft||: " << fnrm		<< endl;
			cout << "   Requested norm:      " << stop_tol	<< endl;
			cout << "   Current ratio:       " << ratio		<< endl;
		}
		// Updating Jacobian and Evaluating Residuals
		if (itc==1 || ratio>rsham_ || itsham==0 || arm_flag==1)
		{
			if (verbose == true)
			{
					 if (itc == 1)		cout << "   Jacobian is updated (First Iteration)"		<< endl;
				else if (ratio>rsham_)	cout << "   Jacobian is updated (Ratio too large)"		<< endl;
				else if (itsham==0)		cout << "   Jacobian is updated (Requested by the user)"	<< endl;
				else if (arm_flag==1)	cout << "   Jacobian is updated (After Armijio Failure)"	<< endl;
			}

			itsham = isham_;
			jac_age = -1;

			double startTime = BzzGetCpuTime();
			flag = ptSolution(x, direction, 1);
			if (flag<0) return -1;
			double endTime = BzzGetCpuTime();
			
			if (verbose == true) 
				cout << "   Elapsed time: " << endTime-startTime << endl;
		}
		// Evaluating Residuals (without updating Jacobian)
		else 
		{
			if (verbose == true)
				cout << "   Jacobian is not updated" << endl;
			
			double startTime = BzzGetCpuTime();
			flag = ptSolution(x, direction, 0);
			if (flag<0) return -1;
			double endTime = BzzGetCpuTime();
			
			if (verbose == true)
				cout << "   Elapsed time: " << endTime-startTime << endl;
		}

		itsham--;
		jac_age++;
		
		xOld = x;
		fOld = f0;
		
		// Armijio Search
		arm_flag = ArmijoStep(x,f0,direction);

		// Failure of Armjio rule
		if (arm_flag == 1)
		{
			// If the Jacobian was old, the procedure can be repeated after updating the jacobian
			if (jac_age > 0)
			{
				if (verbose == true)
					cout << "Armjio failure: Jacobian is recalculated" << endl;
				x  = xOld;
				f0 = fOld;
			}
			// Else the Newton method failed
			else
			{
				if (verbose == true)
					cout << "Complete Armjio failure" << endl;
				solution = xOld;
				return -2;
			}
		}

		fnrm = f0.Norm2();
		ratio = fnrm/fnrm0;

		ratioNormInf->Append(f0.MaxAbs());
		ratioNorm2->Append(fnrm);

		double endTime = BzzGetCpuTime();

		if (verbose == true)
			cout << "CPU Time: " << endTime - startTime << endl;

		if (arm_flag == 0)
			ptAnalysis();

		if ( (ratioNormInf->MeanRatios() > 0.999) && (ratioNorm2->MeanRatios() > 0.999) && itc>=4 )
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

	solution = x;

	if (itc>maxit_)	flag=1;

	return flag;
}

int OpenSMOKE_KPP_NewtonMethod_Manager::ArmijoStep(BzzVector& x, BzzVector& f0, const BzzVector& direction)
{
	if (verbose == true)
	{
		cout << " ******************************************************************** " << endl;
		cout << "                           Armijio Search                             " << endl;
		cout << " ******************************************************************** " << endl;
	}

	int arm_flag = 0;
	int iarm	 = 0;

	double lambda	= 1.;
	double lamm		= 1.;
	double lamc     = lambda;

	Product(lambda, direction, &step);
	Sum(x,step, &xt);
	ptResiduals(xt, ft);

	double nft = ft.Norm2(); 
	double nf0 = f0.Norm2(); 
	double ff0 = nf0*nf0; 
	double ffc = nft*nft; 
	double ffm = nft*nft;

	if (verbose == true)
	{
		cout << " * Starting norm      ||f0||: " << nf0 << endl;
		cout << " * Current norm       ||ft||: " << nft << endl;
		cout << " * Objective function   F0:   " << ff0 << endl;
		cout << " * Objective function   Fc:   " << ffc << endl;
		cout << " * Objective function   Fm:   " << ffm << endl;
	}

	while ( nft >= ((1.-alpha_*lambda)*nf0) )
	{
		if (iarm == 0)	lambda *= 0.50;
		else			lambda  = parabolicModel(lamc, lamm, ff0, ffc, ffm);
		
		Product(lambda, direction, &step);
		Sum(x,step, &xt);
        lamm = lamc;
        lamc = lambda;

		ptResiduals(xt, ft);
        nft = ft.Norm2();
        ffm = ffc;
        ffc = nft*nft;
        iarm++;

		if (verbose == true)
		{
			cout << " * Iteration:               " << iarm << endl;
			cout << "   Lambda:                  " << lambda << endl;
			cout << "   Current norm     ||ft||: " << nft << endl;
			cout << "   Objective function Fc:   " << ffc << endl;
			cout << "   Objective function Fm:   " << ffm << endl;
			cout << "   Required reduction:      " << (1.-alpha_*lambda) << endl;
		}

		if (iarm > maxarm_)
		{
			if (verbose == true)
				cout << "Armijo failure, too many reductions" << endl;

            arm_flag = 1;
            return arm_flag;
		}
	}

	if (verbose == true)
	{
		cout << endl;
		cout << " ********************** Armijio Search Succeeded! ******************* " << endl;
		cout << "  Norm: " << nf0 <<" ==> " << nft << "  Reduction: " << nft/nf0 << endl;
		cout << endl;
	}

	x  = xt;
	f0 = ft;

	return arm_flag;
}
