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

#include "OpenSMOKE_KPP_ODE_Manager.h"
#include "OpenSMOKE_KPP_ReactorNetwork.h"
#include <mpi.h>
#include <petscvec.h>

const double	OpenSMOKE_KPP_ODE_Manager::default_absoluteTolerance_			= 1.e-10;
const double	OpenSMOKE_KPP_ODE_Manager::default_relativeTolerance_			= 1.e-6;
const double	OpenSMOKE_KPP_ODE_Manager::default_initialTimeStep_				= 1.e-4;
const double	OpenSMOKE_KPP_ODE_Manager::default_maxTimeStep_					= 1.e2;
const double	OpenSMOKE_KPP_ODE_Manager::default_timeStepIncrementFactor_		= 2.;
const int		OpenSMOKE_KPP_ODE_Manager::default_maximumIterations_			= 10000;
const int		OpenSMOKE_KPP_ODE_Manager::default_updatingFrequencyTimeStep_	= 3;
const double	OpenSMOKE_KPP_ODE_Manager::default_safetyReductionCoefficient_  = 0.99;

					
void OpenSMOKE_KPP_ODE_Manager::ErrorMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class: OpenSMOKE_KPP_ODE_Manager"	<< std::endl;
    std::cout << "Error: " << message_				<< std::endl;
    std::cout << "Press enter to continue... "		<< std::endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_ODE_Manager::WarningMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class:	  OpenSMOKE_KPP_ODE_Manager"		<< std::endl;
    std::cout << "Warning: " << message_					<< std::endl;
	std::cout << std::endl;
}

OpenSMOKE_KPP_ODE_Manager::OpenSMOKE_KPP_ODE_Manager(BzzVector& x_, void (*ptResiduals_)(BzzVector& x, Vec &f), int (*ptSolution_)(BzzVector& x, BzzVector& dir, double& dt, const bool jacobianFlag), void (*ptAnalysis_)(), OpenSMOKE_KPP_ReactorNetwork& network)
{

	procrank_ = MPI::COMM_WORLD.Get_rank();
	nprocs_ = MPI::COMM_WORLD.Get_size();
	numworkers_ = nprocs_ - 1;

	MASTER = 0;
	FROM_MASTER = 1;
	FROM_WORKER = 2;

	n = network.NumberOfEquations();

	ptResiduals = ptResiduals_;
	ptSolution  = ptSolution_;
	ptAnalysis  = ptAnalysis_;

	ChangeDimensions(n, &x);

	if(procrank_ == 0)
	    ChangeDimensions(n, &step);

	VecCreate(PETSC_COMM_WORLD, &f);
	nResiduals_ = network.LocalNumberOfReactors()[procrank_] * network.NumberOfSpecies();
	VecSetSizes(f, nResiduals_, PETSC_DECIDE);
	VecSetFromOptions(f);

	ratioNormInf		= new RingVector<double>(8);
	ratioNorm1		= new RingVector<double>(8);
	ratioNorm2		= new RingVector<double>(8);

	// Parameters
	absoluteTolerance_					= default_absoluteTolerance_;
	relativeTolerance_					= default_relativeTolerance_;
	updatingFrequencyTimeStep_			= default_updatingFrequencyTimeStep_;
	maximumIterations_					= default_maximumIterations_;
	startingDeltaTime_					= default_initialTimeStep_;
	maximumDeltaTime_					= default_maxTimeStep_;
	increasingFactorDeltaTime_			= default_timeStepIncrementFactor_;

	requestedRatio_						= 0.;
	requestedIterationsJacobian_		= 100;

	requestedDeltaTime_					= startingDeltaTime_;
	numberIterationsWithoutCorrections_ = updatingFrequencyTimeStep_;
	numberIterationsWithCorrections_	= updatingFrequencyTimeStep_;

	// Counters
	iTimeStep_						= true;
	iteration_						= 0;
	iterationsWithoutCorrections_	= 0;
	iterationsWithCorrections_		= 0;
	iterationsJacobian_				= requestedIterationsJacobian_;

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



	ptResiduals(x, f);

	VecNorm(f, NORM_2, &norm2);

	norm2Starting  = norm2 ;
	norm2Previous  = 1.;

	delete [] x_ptr;
}

void OpenSMOKE_KPP_ODE_Manager::SetRelativeTolerance(const double relativeTolerance)
{
	relativeTolerance_ = relativeTolerance;
}

void OpenSMOKE_KPP_ODE_Manager::SetAbsoluteTolerance(const double absoluteTolerance)
{
	absoluteTolerance_ = absoluteTolerance;
}

void OpenSMOKE_KPP_ODE_Manager::SetMaximumNumberOfIterations(const int maximumIterations)
{
	maximumIterations_ = maximumIterations;
}

void OpenSMOKE_KPP_ODE_Manager::SetStartingDeltaTime(const double startingDeltaTime)
{
	startingDeltaTime_	= startingDeltaTime;
	requestedDeltaTime_	= startingDeltaTime_;
}

void OpenSMOKE_KPP_ODE_Manager::SetMaximumDeltaTime(const double maximumDeltaTime)
{
	maximumDeltaTime_ = maximumDeltaTime;
}

void OpenSMOKE_KPP_ODE_Manager::SetIncreasingFactorDeltaTime(const double increasingFactorDeltaTime)
{
	increasingFactorDeltaTime_ = increasingFactorDeltaTime;
}

void OpenSMOKE_KPP_ODE_Manager::SetUpdatingFrequencyTimeStep(const int updatingFrequencyTimeStep)
{
	updatingFrequencyTimeStep_ = updatingFrequencyTimeStep;
	numberIterationsWithoutCorrections_ = updatingFrequencyTimeStep_;
	numberIterationsWithCorrections_	= updatingFrequencyTimeStep_;
}


int OpenSMOKE_KPP_ODE_Manager::Solve(BzzVector& solution)
{
	// Compute the stop tolerance
	int flag=0;
//	MPI::COMM_WORLD.Bcast(&norm2, 1, MPI::DOUBLE, MASTER);
	double stop_tol = absoluteTolerance_ + relativeTolerance_*norm2;


	while (iteration_ < maximumIterations_ && norm2 > stop_tol)
	{
		if(procrank_ == 0)
		    double startTime = MPI::Wtime();

		double deltat = requestedDeltaTime_;

		double ratio  = norm2/norm2Previous;
		norm2Previous = norm2;

		iteration_++;
		
		if(procrank_ == 0)
		{
		    std::cout << " * Iteration:            " << iteration_				<< std::endl;
		    std::cout << "   Requested time:       " << requestedDeltaTime_		<< std::endl;
		    std::cout << "   Current norm  ||f||:  " << norm2					<< std::endl;
		    std::cout << "   Previous step ||s||:  " << step.Norm2()				<< std::endl;
		    std::cout << "   Current ratio:        " << ratio					<< std::endl;
		    std::cout << "   Requested norm:       " << stop_tol					<< std::endl;
		    std::cout << "   Global ratio:         " << norm2/norm2Starting		<< std::endl;
		}
		  
		// Updating Jacobian and Evaluating Residuals
		if (iteration_==1 || ratio>requestedRatio_ || iterationsJacobian_==0 || iTimeStep_==true)
		{
			     if (iteration_ == 1 && procrank_ == 0)			std::cout << "   Jacobian is updated (First Iteration)"			<< std::endl;
			else if (ratio>requestedRatio_ && procrank_ == 0)		std::cout << "   Jacobian is updated (Ratio too large)"			<< std::endl;
			else if (iterationsJacobian_ == 0 && procrank_ == 0)	std::cout << "   Jacobian is updated (Requested by the user)"	<< std::endl;
			else if (iTimeStep_ == true && procrank_ == 0)		std::cout << "   Jacobian is updated (Time step changed)"		<< std::endl;
			
			iterationsJacobian_ = requestedIterationsJacobian_;

			if(procrank_ == 0)
			startTime = MPI::Wtime();

			flag = ptSolution(x, step, deltat, true);

			if (flag<0) return -1;

			if(procrank_ == 0)
			endTime = MPI::Wtime();
			
			if(procrank_ == 0)
			    std::cout << "   Elapsed time: " << endTime-startTime << std::endl;
		}
		// Evaluating Residuals (without updating Jacobian)
		else 
		{
			if(procrank_ == 0)
			std::cout << "   Jacobian is not updated" << std::endl;

			if(procrank_ == 0)
			startTime = MPI::Wtime();

			flag = ptSolution(x, step, deltat, false);

			if (flag<0) return -1;
			if(procrank_ == 0)
			    endTime = MPI::Wtime();

			if(procrank_ == 0)
			std::cout << "   Elapsed time: " << endTime-startTime << std::endl;
		}

		// Update Jacobian counter
		iterationsJacobian_--;

		// Evaluating if correction of time step occurred
		if (deltat < requestedDeltaTime_)
		{
			iterationsWithCorrections_++;
			iterationsWithoutCorrections_ = 0;
		}
		else
		{
			iterationsWithoutCorrections_++;
			iterationsWithCorrections_ = 0;
		}
		TimeStepPolicy();

		// Evaluating new norm of residuals
		ptResiduals(x,f);

		VecNorm(f, NORM_1, &norm1);
		VecNorm(f, NORM_2, &norm2);
		VecNorm(f, NORM_INFINITY, &fMaxAbs);		

		ratioNormInf->Append(fMaxAbs);
		ratioNorm2->Append(norm2);
		ratioNorm1->Append(norm1);
		
		double endTime;
		if(procrank_ == 0)
		    endTime = MPI::Wtime();
		
		if(procrank_ == 0)
		    std::cout << "CPU Time: " << endTime - startTime << std::endl;

		ptAnalysis();
		
		// Slow-Convergence: it is better to increase the time step
		if ( (ratioNormInf->MeanRatios() > 0.99) && (ratioNorm1->MeanRatios() > 0.99)  && iteration_>= 20)
		{
			flag = 2;
			break;
		}

/*		// Very-Slow-Convergence: it is better to increase the time step
		if ( ( (ratioNormInf->MeanRatios() > 0.999) && (ratioNorm1->MeanRatios() > 0.999) ) && ( (ratioNormInf->MeanRatios() < 1.) && (ratioNorm1->MeanRatios() < 1. )  && iteration_>= 10) )
		{
			flag = 3;
			break;
		}*/

		// Super-convergence: switch to Global NLS
		if ( (ratioNormInf->MeanRatios() < 0.750) && (ratioNorm1->MeanRatios() < 0.750)  && iteration_>= 10)
		{
			flag = 4;
			break;
		}

	}

	solution = x;

	if (iteration_ > maximumIterations_)	
		flag=1;

	return flag;
}

void OpenSMOKE_KPP_ODE_Manager::TimeStepPolicy()
{
	iTimeStep_ = false;

	if (iterationsWithoutCorrections_ >= numberIterationsWithoutCorrections_)
	{
		requestedDeltaTime_ *= increasingFactorDeltaTime_;
		requestedDeltaTime_  = std::min(requestedDeltaTime_, maximumDeltaTime_);
		iterationsWithoutCorrections_ = 0;
		iterationsWithCorrections_    = 0;
		
		iTimeStep_ = true;
	}
	
	if (iterationsWithCorrections_ >= numberIterationsWithCorrections_)
	{
		requestedDeltaTime_ /= increasingFactorDeltaTime_;
		iterationsWithoutCorrections_ = 0;
		iterationsWithCorrections_    = 0;

		iTimeStep_ = true;
	}
}
