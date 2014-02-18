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

#include "OpenSMOKE_KPP_SingleReactorStatistics.h"
#include <iostream>
#include <iomanip>
#include <omp.h>

void OpenSMOKE_KPP_SingleReactorStatistics::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_SingleReactorStatistics"	<< endl;
    cout << "Error: " << message_				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_SingleReactorStatistics::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_SingleReactorStatistics"		<< endl;
    cout << "Warning: " << message_					<< endl;
	cout << endl;
}

OpenSMOKE_KPP_SingleReactorStatistics::OpenSMOKE_KPP_SingleReactorStatistics(const int numberOfSpecies, const int numberOfReactors, const string fileName)
{
	iteration_ = 0;

	numberOfSpecies_  = numberOfSpecies;
	numberOfReactors_ = numberOfReactors;

	numberOfSteps_ = 0;
	numberOfFunctions_ = 0;
	numberOfFunctionsForJacobian_ = 0;
	numberOfAnalyticalJacobians_ = 0;
	numberOfNumericalJacobians_ = 0;
	numberOfFactorizations_ = 0;

	maxAbsResidual_  = 0.;
	meanAbsResidual_ = 0.;

	ChangeDimensions( numberOfSpecies_, &residuals_);

	fOutput.open(fileName.c_str(), ios::out);
	fOutput.setf(ios::scientific);
}

void OpenSMOKE_KPP_SingleReactorStatistics::Reset()
{
	numberOfSteps_ = 0;
	numberOfFunctions_ = 0;
	numberOfFunctionsForJacobian_ = 0;
	numberOfAnalyticalJacobians_ = 0;
	numberOfNumericalJacobians_ = 0;
	numberOfFactorizations_ = 0;
	
	maxAbsResidual_  = 0.;
	meanAbsResidual_ = 0.;
}

void OpenSMOKE_KPP_SingleReactorStatistics::Analysis(BzzOdeStiffObject &o)
{
	residuals_ = o.GetY1InMeshPoint();
	
	numberOfSteps_ += o.GetNumStep();
	numberOfFunctions_ += o.GetNumFunction();
	numberOfFunctionsForJacobian_ += o.GetNumFunctionForJacobian();
	numberOfAnalyticalJacobians_ += o.GetNumAnalyticalJacobian();
	numberOfNumericalJacobians_ += o.GetNumNumericalJacobian();
	numberOfFactorizations_ += o.GetNumFactorization();

	maxAbsResidual_   = max(maxAbsResidual_, residuals_.MaxAbs());
	meanAbsResidual_ += residuals_.GetSumAbsElements()/double(numberOfSpecies_);
}

void OpenSMOKE_KPP_SingleReactorStatistics::PrintOnFile()
{
//	double nEquations = double(numberOfReactors_*numberOfSpecies_);
	iteration_++;

	fOutput << setw(16) << left << iteration_;
	fOutput << setw(16) << left << meanAbsResidual_/double(numberOfReactors_);
	fOutput << setw(16) << left << maxAbsResidual_;
	fOutput << setw(16) << left << numberOfSteps_/double(numberOfReactors_);
	fOutput << setw(16) << left << numberOfFunctions_/double(numberOfReactors_);
	fOutput << setw(16) << left << numberOfFunctionsForJacobian_/double(numberOfReactors_);
	fOutput << setw(16) << left << numberOfAnalyticalJacobians_/double(numberOfReactors_);
	fOutput << setw(16) << left << numberOfNumericalJacobians_/double(numberOfReactors_);
	fOutput << setw(16) << left << numberOfFactorizations_/double(numberOfReactors_);
	fOutput << endl;
}
