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

#ifndef OpenSMOKE_KPP_NewtonMethod_Manager_H
#define OpenSMOKE_KPP_NewtonMethod_Manager_H

#include "BzzMath.hpp"
#include "OpenSMOKE_KPP_Definitions.h"
#include <petscvec.h>

class OpenSMOKE_KPP_ReactorNetwork;

class OpenSMOKE_KPP_NewtonMethod_Manager
{
public:

	OpenSMOKE_KPP_NewtonMethod_Manager();
	OpenSMOKE_KPP_NewtonMethod_Manager(BzzVector& xStarting_, void (*ptResiduals_)(BzzVector& x, Vec &f), int (*ptSolution_)(BzzVector& x, BzzVector &f, const bool jacobianFlag), void (*ptAnalysis)(), OpenSMOKE_KPP_ReactorNetwork& network);
	void Setup(BzzVector& xStarting_, void (*ptResiduals_)(BzzVector& x, Vec &f), int (*ptSolution_)(BzzVector& x, BzzVector &f, const bool jacobianFlag), void (*ptAnalysis)(), OpenSMOKE_KPP_ReactorNetwork& network);

	inline bool convergence() const {return convergence_;}

	int Solve(BzzVector& solution);

	void SetMethod(const KPP_NonLinearSystem_Method method, const int m = default_shamanskii);
	void SetRelativeTolerance(const double relativeTolerance);
	void SetAbsoluteTolerance(const double absoluteTolerance);
	void SetMaximumNumberOfIterations(const int maximumIterations);
	void SetMaximumNumberOfArmijioIterations(const int maximumArmijioIterations);
	void SetAlpha(const double alpha);
	
private:

	int	n;
	double startTime, endTime;
	double stop_tol;
	
	// Iterators
	int itc;
	int itsham;

	bool convergence_;

	// Parameters
	KPP_NonLinearSystem_Method method_;
	int    maxit_;
	int    maxarm_;
	double alpha_;
	int    isham_;
	double rsham_;
	double relativeTolerance_;
	double absoluteTolerance_;
	double f0MaxAbs;
	double *xt_ptr;

	int f0size_;
	int fOldsize_;
	int ftsize_;

	// Norms
	double fnrm;
	double fnrm0;

	// Vectors
	BzzVector xt;
	Vec ft;
	BzzVector xOld;
	Vec fOld;
	BzzVector x;
	Vec f0;
	BzzVector step;
	BzzVector direction;

	bool verbose;

	//Parallel utilities
	int procrank_, numworkers_, nprocs_, MASTER, FROM_MASTER, FROM_WORKER;

	RingVector<double>* ratioNormInf;
	RingVector<double>* ratioNorm2;

private:

	void (*ptResiduals)(BzzVector& x, Vec& f);
	int  (*ptSolution)(BzzVector& x, BzzVector& dir, const bool jacobianFlag);
	void (*ptAnalysis)();
	int ArmijoStep(BzzVector& x, Vec& f0, const BzzVector& direction);

private:

	void ErrorMessage(const std::string message_);
	void WarningMessage(const std::string message_);

public:

	static const KPP_NonLinearSystem_Method  default_method_;
	static const int		default_maxit_;
	static const int		default_maxarm_;
	static const int		default_shamanskii;
	static const double		default_alpha_;
	static const double		default_relativeTolerance_;
	static const double		default_absoluteTolerance_;
	static const double		default_safetyReductionCoefficient_;
};

#endif // OpenSMOKE_KPP_NewtonMethod_Manager_H
