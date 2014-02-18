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

#ifndef OpenSMOKE_KPP_ODE_Manager_H
#define OpenSMOKE_KPP_ODE_Manager_H

#include "BzzMath.hpp"
#include "OpenSMOKE_KPP_Definitions.h"

class OpenSMOKE_KPP_ODE_Manager
{
public:
	
	OpenSMOKE_KPP_ODE_Manager(BzzVector& xStarting_, void (*ptResiduals_)(BzzVector& x, BzzVector &f), int (*ptSolution_)(BzzVector& x, BzzVector &f, double& dt, const bool jacobianFlag), void (*ptAnalysis)());
	
	int Solve(BzzVector& solution);

	void SetRelativeTolerance(const double relativeTolerance);
	void SetAbsoluteTolerance(const double absoluteTolerance);
	void SetMaximumNumberOfIterations(const int maximumIterations);
	void SetStartingDeltaTime(const double startingDeltaTime);
	void SetMaximumDeltaTime(const double maximumDeltaTime);
	void SetIncreasingFactorDeltaTime(const double increasingFactorDeltaTime);
	void SetUpdatingFrequencyTimeStep(const int updatingFrequencyTimeStep);

	inline double LastRequestedDeltaTime() const { return requestedDeltaTime_; }
	
private:

	int	   n;
	bool   iTimeStep_;
	double norm2;
	double norm2Previous;
	double norm2Starting;

	// Vectors
	BzzVector x;
	BzzVector step;
	BzzVector f;

	// Parameters
	double startingDeltaTime_;
	double maximumDeltaTime_;
	double increasingFactorDeltaTime_;
	double requestedDeltaTime_;
	int numberIterationsWithoutCorrections_;
	int numberIterationsWithCorrections_;
	int updatingFrequencyTimeStep_;

	double absoluteTolerance_;
	double relativeTolerance_;
	int maximumIterations_;
	double requestedRatio_;
	int requestedIterationsJacobian_;

	// Counters
	int iteration_;
	int iterationsJacobian_;
	int iterationsWithoutCorrections_;
	int iterationsWithCorrections_;

	RingVector<double>* ratioNormInf;
	RingVector<double>* ratioNorm2;

private:

	void (*ptResiduals)(BzzVector& x, BzzVector& f);
	int  (*ptSolution)(BzzVector& x, BzzVector& dir, double& dt, const bool jacobianFlag);
	void (*ptAnalysis)();

	void TimeStepPolicy();

private:

	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);

public:

	static const double	default_absoluteTolerance_;
	static const double	default_relativeTolerance_;
	static const double	default_initialTimeStep_;
	static const double	default_maxTimeStep_;
	static const double	default_timeStepIncrementFactor_;
	static const int	default_maximumIterations_;
	static const int	default_updatingFrequencyTimeStep_;
	static const double default_safetyReductionCoefficient_;

};

#endif // OpenSMOKE_KPP_ODE_Manager_H