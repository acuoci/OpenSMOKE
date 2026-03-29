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

#ifndef OpenSMOKE_KPP_SingleReactorStatistics_H
#define OpenSMOKE_KPP_SingleReactorStatistics_H

#include "BzzMath.hpp"
#include "OpenSMOKE_KPP_Definitions.h"
#include "OpenSMOKE_KPP_Communicator.h"

class OpenSMOKE_KPP_Communicator;

class OpenSMOKE_KPP_SingleReactorStatistics
{
public:
	
	OpenSMOKE_KPP_SingleReactorStatistics(const int numberOfSpecies, const int numberOfReactors, const std::string fileName, OpenSMOKE_KPP_Communicator* communicator);
	void Reset();
	void Analysis(BzzOdeStiffObject &o);
	void PrintOnFile();
	
private:

	std::ofstream fOutput;
	
	int numberOfSpecies_;
	int numberOfReactors_;

	int iteration_;
	int numberOfSteps_;
	int numberOfFunctions_;
	int numberOfFunctionsForJacobian_;
	int numberOfAnalyticalJacobians_;
	int numberOfNumericalJacobians_;
	int numberOfFactorizations_;

	int *numberOfStepsGlob_;
	int *numberOfFunctionsGlob_;
	int *numberOfFunctionsForJacobianGlob_;
	int *numberOfAnalyticalJacobiansGlob_;
	int *numberOfNumericalJacobiansGlob_;
	int *numberOfFactorizationsGlob_;

	BzzVector residuals_;
	double	maxAbsResidual_;
	double	meanAbsResidual_;

	double	*maxAbsResidualGlob_;
	double	*meanAbsResidualGlob_;

	//Parallel utilities
	int nprocs_, procrank_, numworkers_;
	int MASTER, source, mtype, FROM_MASTER, FROM_WORKER;

private:

	void ErrorMessage(const std::string message_);
	void WarningMessage(const std::string message_);
	OpenSMOKE_KPP_Communicator* communicator_;
};


#endif // OpenSMOKE_KPP_SingleReactorStatistics_H
