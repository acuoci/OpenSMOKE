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

#include <iostream>
#include <iomanip>
#include <sstream>
#include "OpenSMOKE_KPP_ReactorNetwork.h"
#include "OpenSMOKE_KPP_DataManager.h"
#include "OpenSMOKE_KPP_ReactorNetwork_Residuals.h"

const int OpenSMOKE_KPP_ReactorNetwork_Residuals::nRanges_		= 12;
const double OpenSMOKE_KPP_ReactorNetwork_Residuals::endingRange_	= 1.e-10;
const double OpenSMOKE_KPP_ReactorNetwork_Residuals::startingRange_	= 1.e-4;

void OpenSMOKE_KPP_ReactorNetwork_Residuals::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_ReactorNetwork_Residuals"	<< endl;
    cout << "Error: " << message_				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_ReactorNetwork_Residuals"		<< endl;
    cout << "Warning: " << message_					<< endl;
	cout << endl;
}

OpenSMOKE_KPP_ReactorNetwork_Residuals::OpenSMOKE_KPP_ReactorNetwork_Residuals(OpenSMOKE_KPP_ReactorNetwork& network) :
network_(network)
{
	nTimeSteps_ = 10;
	nFrequencyRefinedAnalysis_ = 10;

	ChangeDimensions(network_.NumberOfEquations(), &residuals_);

	ChangeDimensions(network_.NumberOfReactors(), &speciesResiduals_);
	ChangeDimensions(network_.NumberOfSpecies(),  &speciesNormInf_);
	ChangeDimensions(network_.NumberOfSpecies(),  &speciesNorm1_);
	ChangeDimensions(network_.NumberOfSpecies(),  &speciesNorm2_);

	ChangeDimensions(network_.NumberOfSpecies(),  &reactorResiduals_);
	ChangeDimensions(network_.NumberOfReactors(), &reactorNormInf_);
	ChangeDimensions(network_.NumberOfReactors(), &reactorNorm1_);
	ChangeDimensions(network_.NumberOfReactors(), &reactorNorm2_);

	ChangeDimensions(nTimeSteps_, &historyNormInf_);
	ChangeDimensions(nTimeSteps_, &historyNorm1_);
	ChangeDimensions(nTimeSteps_, &historyNorm2_);

	openOutputFileAndControl(fResiduals, network_.data().nameResidualFile());
	fResiduals.setf(ios::scientific);
	WriteResidualsLabels();

	openOutputFileAndControl(fResidualsSpecies, network_.data().nameResidualSpeciesFile());
	fResidualsSpecies.setf(ios::scientific);
	WriteResidualsSpeciesLabels();

	openOutputFileAndControl(fResidualsReactorStatistics, network_.data().nameResidualReactorStatisticsFile());
	fResidualsReactorStatistics.setf(ios::scientific);
	openOutputFileAndControl(fResidualsSpeciesStatistics, network_.data().nameResidualSpeciesStatisticsFile());
	fResidualsSpeciesStatistics.setf(ios::scientific);
	WriteResidualsStatisticsLabels();

	// Statistical data
	{
		double alfa = pow(startingRange_/endingRange_, 1./(nRanges_-1));
		ChangeDimensions( nRanges_, &ranges_);
		ranges_[1] = startingRange_;
		for (int i=2;i<=ranges_.Size();i++)
			ranges_[i] = ranges_[i-1]/alfa;
	}

	currentStatus_   = network_.status();
	localIteration_  = 0;

	ratioNormInf	= new RingVector<double>(10);
	ratioNorm1		= new RingVector<double>(10);
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::SetNumberOfTimeSteps(const int nTimeSteps)
{
	nTimeSteps_ = nTimeSteps;
	ChangeDimensions(nTimeSteps_, &historyNormInf_);
	ChangeDimensions(nTimeSteps_, &historyNorm1_);
	ChangeDimensions(nTimeSteps_, &historyNorm2_);
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::Analysis()
{	
	// Convenctional Analysis
	{
		maxSpecies = 0.;
		minSpecies = 1.;

		// Build residuals vector
		for (int k=1;k<=network_.NumberOfReactors();k++)
		{
			int position = (k-1)*network_.NumberOfSpecies() + 1;
			network_.reactors(k).Residuals(position, residuals_, network_);

			maxSpecies = max( network_.reactors(k).omegaMax(), maxSpecies);
			minSpecies = min( network_.reactors(k).omegaMin(), minSpecies);
		}

		normInf_	= residuals_.MaxAbs();
		norm1_		= residuals_.GetSumAbsElements();
		norm2_		= residuals_.Norm2();

		ratioNormInf->Append(normInf_);
		ratioNorm1->Append(norm1_);

		// Update current status
		if ( currentStatus_ != network_.status() || network_.status()==KPP_NETWORK_STATUS_START )
		{
			currentStatus_		= network_.status();
			localIteration_		= 0;
			startingNormInf_	= normInf_;
			startingNorm1_		= norm1_;
			startingNorm2_		= norm2_;
		}
		localIteration_++;

		// Write on file
		WriteResidualsOnFile();
	}

	// Refined Analysis
	if ( (network_.iteration()-1)%nFrequencyRefinedAnalysis_ == 0 )
	{
		AnalysisRefined();
		WriteResidualsSpeciesOnFile();
		Statistics();
	}
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::AnalysisRefined()
{
	// Species residuals
	for (int j=1;j<=network_.NumberOfSpecies();j++)
	{
		int count = j;
		for (int k=1;k<=network_.NumberOfReactors();k++)
			speciesResiduals_[k] = residuals_[count+=network_.NumberOfSpecies()];
		
		speciesNormInf_[j]	= speciesResiduals_.MaxAbs();
		speciesNorm1_[j]	= speciesResiduals_.GetSumAbsElements()/double(network_.NumberOfReactors());
		speciesNorm2_[j]	= speciesResiduals_.Norm2()/double(network_.NumberOfReactors());
	}

	// Reactor residuals
	int count=1;
	for (int k=1;k<=network_.NumberOfReactors();k++)
	{
		for (int j=1;j<=network_.NumberOfSpecies();j++)
			reactorResiduals_[j] = residuals_[count++];

		reactorNormInf_[k]	= reactorResiduals_.MaxAbs();
		reactorNorm1_[k]	= reactorResiduals_.GetSumAbsElements()/double(network_.NumberOfSpecies());
		reactorNorm2_[k]	= reactorResiduals_.Norm2()/double(network_.NumberOfSpecies());
	}
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsOnVideo()
{
	cout << endl;
	cout << "// ********************************************************************************* // " << endl;
	cout << "//                               Residual Analysis                                   // " << endl;
	cout << "// ********************************************************************************* // " << endl;
	cout << "   #   normInf        norm1          norm2          minOmega       maxOmega             " << endl;

	cout << "       ";
	cout << setw(15) << left << scientific << normInf_;
	cout << setw(15) << left << scientific << norm1_/double( network_.NumberOfEquations() );
	cout << setw(15) << left << scientific << norm2_/double( network_.NumberOfEquations() );
	cout << setw(15) << left << scientific << minSpecies;
	cout << setw(15) << left << scientific << maxSpecies;
	cout << setw(15) << left << scientific << ratioNormInf->MeanRatios();
	cout << setw(15) << left << scientific << ratioNorm1->MeanRatios();
	cout << endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsOnFile()
{
	fResiduals << setw(9)  << left << network_.iteration();
	fResiduals << setw(9)  << left << localIteration_;
	fResiduals << setw(9)  << left << currentStatus_;
	fResiduals << setw(16) << left << normInf_;
	fResiduals << setw(16) << left << norm1_/double( network_.NumberOfEquations() );
	fResiduals << setw(16) << left << norm2_/double( network_.NumberOfEquations() );
	fResiduals << setw(16) << left << normInf_/startingNormInf_;
	fResiduals << setw(16) << left << norm1_/startingNorm1_;
	fResiduals << setw(16) << left << norm2_/startingNorm2_;
	fResiduals << setw(16) << left << norm1_;
	fResiduals << setw(16) << left << norm2_;
	fResiduals << endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsSpeciesOnFile()
{
	fResidualsSpecies << setw(9)  << left << network_.iteration();
	for (int j=1;j<=network_.NumberOfSpecies();j++)
		fResidualsSpecies << setw(16)  << left << speciesNorm2_[j];
	fResidualsSpecies << endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsLabels()
{
	fResiduals << setw(9)   << left << "Glob. #";
	fResiduals << setw(9)   << left << "Local #";
	fResiduals << setw(9)   << left << "Status";
	fResiduals << setw(16)  << left << "NormInf";
	fResiduals << setw(16)  << left << "MeanNorm1";
	fResiduals << setw(16)  << left << "MeanNorm2";
	fResiduals << setw(16)  << left << "NormInf(ratio)";
	fResiduals << setw(16)  << left << "Norm1(ratio)";
	fResiduals << setw(16)  << left << "Norm2(ratio)";
	fResiduals << setw(16)  << left << "Norm1";
	fResiduals << setw(16)  << left << "Norm2";
	fResiduals << endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsStatisticsLabels()
{
	fResidualsReactorStatistics << setw(9)  << left << "Iter.";
	for (int j=1;j<=nRanges_;j++)
		fResidualsReactorStatistics << setw(16)  << left << "range";
	fResidualsReactorStatistics << endl;

	fResidualsSpeciesStatistics << setw(9)  << left << "Iter.";
	for (int j=1;j<=nRanges_;j++)
		fResidualsSpeciesStatistics << setw(16)  << left << "range";
	fResidualsSpeciesStatistics << endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsSpeciesLabels()
{
	fResidualsSpecies << setw(9)  << left << "Iter.";
	
	int count = 2;
	for (int j=1;j<=network_.NumberOfSpecies();j++)
	{
		stringstream index; index << count++;
		string label = network_.mix().names[j] + "(" + index.str() + ")";
		
		fResidualsSpecies << setw(16)  << left << label;
	}
	fResidualsSpecies << endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::Statistics()
{
	// Reactor statistics
	{
		BzzVectorInt countRanges(ranges_.Size()+1);

		for(int k=1;k<=network_.NumberOfReactors();k++)
		{
			bool iFound=false;
			for(int i=ranges_.Size();i>=1;i--)
				if (reactorNorm2_[k] <= ranges_[i])
				{
					iFound=true;
					countRanges[i+1]++;
					break;
				}
	
			if (iFound == false)
				countRanges[1]++;
		}

		fResidualsReactorStatistics << setw(9)  << left << network_.iteration();
		for (int i=1;i<=ranges_.Size()+1;i++)
			fResidualsReactorStatistics << setw(16)  << left << 100.*double(countRanges[i])/double(network_.NumberOfReactors());
		fResidualsReactorStatistics << endl;
	}

	// Species statistics
	{
		BzzVectorInt countRanges(ranges_.Size()+1);

		for(int k=1;k<=network_.NumberOfSpecies();k++)
		{
			bool iFound=false;
			for(int i=ranges_.Size();i>=1;i--)
				if (speciesNorm2_[k] <= ranges_[i])
				{
					iFound=true;
					countRanges[i+1]++;
					break;
				}
	
			if (iFound == false)
				countRanges[1]++;
		}

		fResidualsSpeciesStatistics << setw(9)  << left << network_.iteration();
		for (int i=1;i<=ranges_.Size()+1;i++)
			fResidualsSpeciesStatistics << setw(16)  << left << 100.*double(countRanges[i])/double(network_.NumberOfSpecies());
		fResidualsSpeciesStatistics << endl;
	}
}
