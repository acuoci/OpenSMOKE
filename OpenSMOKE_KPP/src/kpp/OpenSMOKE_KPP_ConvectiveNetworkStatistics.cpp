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

#include "OpenSMOKE_KPP_ConvectiveNetworkStatistics.h"
#include <iostream>
#include <iomanip>
#include <omp.h>

void OpenSMOKE_KPP_ConvectiveNetworkStatistics::ErrorMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class: OpenSMOKE_KPP_ConvectiveNetworkStatistics"	<< std::endl;
    std::cout << "Error: " << message_				<< std::endl;
    std::cout << "Press enter to continue... "		<< std::endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_ConvectiveNetworkStatistics::WarningMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class:	  OpenSMOKE_KPP_ConvectiveNetworkStatistics"		<< std::endl;
    std::cout << "Warning: " << message_					<< std::endl;
	std::cout << std::endl;
}

OpenSMOKE_KPP_ConvectiveNetworkStatistics::OpenSMOKE_KPP_ConvectiveNetworkStatistics(const int numberOfSpecies, const int numberOfReactors, const std::string fileName)
{
	iteration_ = 0;
	iterationsWithoutCorrections_ = 0;
	iterationsWithCorrections_    = 0;

	numberOfSpecies_  = numberOfSpecies;
	numberOfReactors_ = numberOfReactors;

	currentTimeStep_   = 0.;
	requestedTimeStep_ = 0.;
	iTimeCorrectedMin_ = 0;
	iTimeCorrectedMax_ = 0;
	uncorrectedOmegaMin_ =  1e16;
	uncorrectedOmegaMax_ = -1e16;
	sumCorrectedOmegaMin_ = 1.;
	sumCorrectedOmegaMax_ = 1.;

	iRobust_ = true;

	ChangeDimensions( numberOfReactors_, &sumOmega_);

	fOutput.open(fileName.c_str(), std::ios::out);
	fOutput.setf(std::ios::scientific);
}

void OpenSMOKE_KPP_ConvectiveNetworkStatistics::Reset()
{
	iTimeCorrectedMin_ = 0;
	iTimeCorrectedMax_ = 0;
	uncorrectedOmegaMin_ =  1e16;
	uncorrectedOmegaMax_ = -1e16;
}

void OpenSMOKE_KPP_ConvectiveNetworkStatistics::Analysis()
{
	if (iTimeCorrectedMin_ == 0 && iTimeCorrectedMax_ == 0)
	{
		iterationsWithoutCorrections_++;
		iterationsWithCorrections_ = 0;
	}
	else
	{
		iterationsWithCorrections_++;
		iterationsWithoutCorrections_ = 0;
	}
}

void OpenSMOKE_KPP_ConvectiveNetworkStatistics::PrintOnFile()
{
//	double nEquations = double(numberOfReactors_*numberOfSpecies_);
	iteration_++;

	fOutput << std::setw(16) << std::left << iteration_;
	fOutput << std::setw(16) << std::left << iterationsWithoutCorrections_;
	fOutput << std::setw(16) << std::left << iterationsWithCorrections_;
	fOutput << std::setw(16) << std::left << currentTimeStep_;
	fOutput << std::setw(16) << std::left << requestedTimeStep_;
	fOutput << std::setw(16) << std::left << currentTimeStep_/requestedTimeStep_;
	fOutput << std::setw(16) << std::left << iTimeCorrectedMin_;
	fOutput << std::setw(16) << std::left << iTimeCorrectedMax_;
	fOutput << std::setw(16) << std::left << uncorrectedOmegaMin_;
	fOutput << std::setw(16) << std::left << uncorrectedOmegaMax_;
	fOutput << std::setw(16) << std::left << sumCorrectedOmegaMin_;
	fOutput << std::setw(16) << std::left << sumCorrectedOmegaMax_;
	fOutput << std::endl;
}

double OpenSMOKE_KPP_ConvectiveNetworkStatistics::ReductionOfTimeStep(BzzVector& vNew, BzzVector& vOld)
{
	double a = 1.;

	// Values lower than 0
	{
		if (iRobust_ == true)
		{
			for (int i=1;i<=vNew.Size();i++)
			{
				if(vNew[i]<-1.e-16)
				{
					iTimeCorrectedMin_++;

					a = std::min(a,vOld[i]/(vOld[i]-vNew[i]));
					if (a<0. || a>1.)
						ErrorMessage("Safety reduction factor outside the boundaries...");
				}
			}
		}

		else
		{
			int kMin;
			uncorrectedOmegaMin_ = vNew.Min(&kMin);
			if (uncorrectedOmegaMin_ > 1.)
			{
				iTimeCorrectedMin_++;
				a = std::min(a, vOld[kMin]/(vOld[kMin]-vNew[kMin]));
			}
		}
	}
		
	// Values higher than 1
	{
		int kMax;
		uncorrectedOmegaMax_ = vNew.Max(&kMax);
		if (uncorrectedOmegaMax_ > 1.)
		{
			iTimeCorrectedMax_++;
			a = std::min(a, (1.-vOld[kMax])/(vNew[kMax]-vOld[kMax]));
		}
	}

	return a;
}

void OpenSMOKE_KPP_ConvectiveNetworkStatistics::Analysis(const double deltat, const double a, BzzVector& omegaNew)
{
	sumOmega_ = 0.;
	for(int k=1;k<=numberOfReactors_;k++)
	{
		int position = (k-1)*numberOfSpecies_;
		for(int j=1;j<=numberOfSpecies_;j++)
			sumOmega_[k] += omegaNew[position+j];
	}

	sumCorrectedOmegaMin_ = sumOmega_.Min();
	sumCorrectedOmegaMax_ = sumOmega_.Max();

	currentTimeStep_    = deltat*a;
	requestedTimeStep_  = deltat;

	Analysis();
	PrintOnFile();

	// Data on video
	std::cout << " * Current correction time step:      " << a << std::endl;
	std::cout << " * Requested time step:               " << deltat << std::endl;
	std::cout << " * Current  time step:                " << deltat*a << std::endl;
}
