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

#include <algorithm> 
#include "OpenSMOKE_KPP_ODE_Manager.h"
using namespace std;

const double	OpenSMOKE_KPP_ODE_Manager::default_absoluteTolerance_			= 1.e-10;
const double	OpenSMOKE_KPP_ODE_Manager::default_relativeTolerance_			= 1.e-6;
const double	OpenSMOKE_KPP_ODE_Manager::default_initialTimeStep_				= 1.e-4;
const double	OpenSMOKE_KPP_ODE_Manager::default_maxTimeStep_					= 1.e2;
const double	OpenSMOKE_KPP_ODE_Manager::default_timeStepIncrementFactor_		= 2.;
const int		OpenSMOKE_KPP_ODE_Manager::default_maximumIterations_			= 10000;
const int		OpenSMOKE_KPP_ODE_Manager::default_updatingFrequencyTimeStep_	= 3;
const double	OpenSMOKE_KPP_ODE_Manager::default_safetyReductionCoefficient_  = 0.99;

					
void OpenSMOKE_KPP_ODE_Manager::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_ODE_Manager"	<< endl;
    cout << "Error: " << message_				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_ODE_Manager::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_ODE_Manager"		<< endl;
    cout << "Warning: " << message_					<< endl;
	cout << endl;
}

OpenSMOKE_KPP_ODE_Manager::OpenSMOKE_KPP_ODE_Manager(BzzVector& x_, void (*ptResiduals_)(BzzVector& x, BzzVector &f), int (*ptSolution_)(BzzVector& x, BzzVector& dir, double& dt, const bool jacobianFlag), void (*ptAnalysis_)())
{
	n = x_.Size();

	ptResiduals = ptResiduals_;
	ptSolution  = ptSolution_;
	ptAnalysis  = ptAnalysis_;

	ChangeDimensions(n, &x);
	ChangeDimensions(n, &step);
	ChangeDimensions(n, &f);

	ratioNormInf	= new RingVector<double>(6);
	ratioNorm2		= new RingVector<double>(6);

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
	x = x_;
	ptResiduals(x, f);
	norm2          = f.Norm2();
	norm2Starting  = norm2 ;
	norm2Previous  = 1.;
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
	double stop_tol = absoluteTolerance_ + relativeTolerance_*norm2;

	while (iteration_ < maximumIterations_ && norm2 > stop_tol)
	{
		double startTime = BzzGetCpuTime();
		double deltat = requestedDeltaTime_;

		double ratio  = norm2/norm2Previous;
		norm2Previous = norm2;
		
		iteration_++;

		cout << " * Iteration:            " << iteration_				<< endl;
		cout << "   Requested time:       " << requestedDeltaTime_		<< endl;
		cout << "   Current norm  ||f||:  " << norm2					<< endl;
		cout << "   Previous step ||s||:  " << step.Norm2()				<< endl;
		cout << "   Current ratio:        " << ratio					<< endl;
		cout << "   Requested norm:       " << stop_tol					<< endl;
		cout << "   Global ratio:         " << norm2/norm2Starting		<< endl;
		  
		// Updating Jacobian and Evaluating Residuals
		if (iteration_==1 || ratio>requestedRatio_ || iterationsJacobian_==0 || iTimeStep_==true)
		{
			     if (iteration_ == 1)		cout << "   Jacobian is updated (First Iteration)"			<< endl;
			else if (ratio>requestedRatio_)		cout << "   Jacobian is updated (Ratio too large)"			<< endl;
			else if (iterationsJacobian_ == 0)	cout << "   Jacobian is updated (Requested by the user)"	<< endl;
			else if (iTimeStep_ == true)		cout << "   Jacobian is updated (Time step changed)"		<< endl;
			
			iterationsJacobian_ = requestedIterationsJacobian_;

			double startTime = BzzGetCpuTime();
			flag = ptSolution(x, step, deltat, true);
			if (flag<0) return -1;
			double endTime = BzzGetCpuTime();
			cout << "   Elapsed time: " << endTime-startTime << endl;
		}
		// Evaluating Residuals (without updating Jacobian)
		else 
		{
			cout << "   Jacobian is not updated" << endl;
			double startTime = BzzGetCpuTime();
			flag = ptSolution(x, step, deltat, false);
			if (flag<0) return -1;
			double endTime = BzzGetCpuTime();
			cout << "   Elapsed time: " << endTime-startTime << endl;
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
		norm2 = f.Norm2();

		ratioNormInf->Append(f.MaxAbs());
		ratioNorm2->Append(norm2);
		
		double endTime = BzzGetCpuTime();

		cout << "CPU Time: " << endTime - startTime << endl;

		ptAnalysis();

		// Slow-Convergence: it is better to increase the time step
		if ( (ratioNormInf->MeanRatios() > 0.99) && (ratioNorm2->MeanRatios() > 0.99)  && iteration_>=50)
		{
			flag = 2;
			break;
		}

		// Very-Slow-Convergence: it is better to increase the time step
	//	if ( ( (ratioNormInf->MeanRatios() > 0.999) && (ratioNorm2->MeanRatios() > 0.999) ) && iteration_>=10)
	//	{
	//		flag = 3;
	//		break;
	//	}

		// Super-convergence: switch to Global NLS
	//	if ( (ratioNormInf->MeanRatios() < 0.700) && (ratioNorm2->MeanRatios() < 0.700)  && iteration_>=10)
	//	{
	//		flag = 4;
	//		break;
	//	}
	}

	solution = x;

	if (iteration_ > maximumIterations_)	
		flag=1;;

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
