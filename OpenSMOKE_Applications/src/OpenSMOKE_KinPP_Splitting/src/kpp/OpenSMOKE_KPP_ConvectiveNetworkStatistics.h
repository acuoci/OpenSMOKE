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

#ifndef OpenSMOKE_KPP_ConvectiveNetworkStatistics_H
#define OpenSMOKE_KPP_ConvectiveNetworkStatistics_H

#include "BzzMath.hpp"

class OpenSMOKE_KPP_ConvectiveNetworkStatistics
{
friend class OpenSMOKE_KPP_ReactorNetwork;

public:
	
	OpenSMOKE_KPP_ConvectiveNetworkStatistics(const int numberOfSpecies, const int numberOfReactors, const string fileName);
	double ReductionOfTimeStep(BzzVector& vNew, BzzVector& vOld);
	void Reset();
	void Analysis();
	void Analysis(const double deltat, const double a, BzzVector& omegaNew);
	void PrintOnFile();

	void SetRobust()	{ iRobust_ = true;}
	void UnsetRobust()	{ iRobust_ = false;}

private:

	ofstream fOutput;
	
	int numberOfSpecies_;
	int numberOfReactors_;

	int iteration_;
	int iterationsWithoutCorrections_;
	int iterationsWithCorrections_;

	BzzVector  sumOmega_;

	double currentTimeStep_;
	double requestedTimeStep_;
	double iTimeCorrectedMin_;
	double iTimeCorrectedMax_;
	double uncorrectedOmegaMin_;
	double uncorrectedOmegaMax_;
	double sumCorrectedOmegaMin_;
	double sumCorrectedOmegaMax_;

	bool iRobust_;

private:

	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
};

#endif // OpenSMOKE_KPP_ConvectiveNetworkStatistics_H