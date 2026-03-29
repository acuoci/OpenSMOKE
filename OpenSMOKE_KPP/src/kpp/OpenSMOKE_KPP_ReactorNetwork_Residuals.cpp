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
#include <mpi.h>
#include <petscvec.h>
#include "OpenSMOKE_KPP_ReactorNetwork.h"
#include "OpenSMOKE_KPP_DataManager.h"
#include "OpenSMOKE_KPP_ReactorNetwork_Residuals.h"
#include "OpenSMOKE_KPP_Communicator.h"

const int OpenSMOKE_KPP_ReactorNetwork_Residuals::nRanges_		= 12;
const double OpenSMOKE_KPP_ReactorNetwork_Residuals::endingRange_	= 1.e-10;
const double OpenSMOKE_KPP_ReactorNetwork_Residuals::startingRange_	= 1.e-4;

void OpenSMOKE_KPP_ReactorNetwork_Residuals::ErrorMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class: OpenSMOKE_KPP_ReactorNetwork_Residuals"	<< std::endl;
    std::cout << "Error: " << message_				<< std::endl;
    std::cout << "Press enter to continue... "		<< std::endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WarningMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class:	  OpenSMOKE_KPP_ReactorNetwork_Residuals"		<< std::endl;
    std::cout << "Warning: " << message_					<< std::endl;
	std::cout << std::endl;
}

OpenSMOKE_KPP_ReactorNetwork_Residuals::OpenSMOKE_KPP_ReactorNetwork_Residuals(OpenSMOKE_KPP_ReactorNetwork& network, OpenSMOKE_KPP_Communicator* communicator) :
network_(network), communicator_(communicator)
{
	nTimeSteps_ = 10;
	nFrequencyRefinedAnalysis_ = 10;

	nprocs_ = MPI::COMM_WORLD.Get_size();
	numworkers_ = nprocs_ - 1;
	procrank_ = MPI::COMM_WORLD.Get_rank();

	NR_P = network_.LocalNumberOfReactors();
	offset = network_.offs();

	MASTER = 0;
	FROM_MASTER = 1;
	FROM_WORKER = 2;

	InitializeResidualsVector();
	InitializeSpeciesResiduals();
	communicator_->InitializeArray(r_residualsMaxAbs);
	communicator_->InitializeArray(r_residualsNorm1);
	communicator_->InitializeArray(r_residualsNorm2);

	if(procrank_ == 0)
	{
//	    ChangeDimensions(network_.NumberOfEquations(), &residuals_);

	    ChangeDimensions(network_.NumberOfReactors(), &speciesResiduals_);
	    ChangeDimensions(network_.NumberOfSpecies(),  &speciesNormInf_);
	    ChangeDimensions(network_.NumberOfSpecies(),  &speciesNorm1_);
	    ChangeDimensions(network_.NumberOfSpecies(),  &speciesNorm2_);

	    ChangeDimensions(network_.NumberOfReactors(), &reactorNormInf_);
	    ChangeDimensions(network_.NumberOfReactors(), &reactorNorm1_);
	    ChangeDimensions(network_.NumberOfReactors(), &reactorNorm2_);

	    ChangeDimensions(nTimeSteps_, &historyNormInf_);
	    ChangeDimensions(nTimeSteps_, &historyNorm1_);
	    ChangeDimensions(nTimeSteps_, &historyNorm2_);
	
	    openOutputFileAndControl(fResiduals, network_.data().nameResidualFile());
	    fResiduals.setf(std::ios::scientific);
	    WriteResidualsLabels();

	    openOutputFileAndControl(fResidualsSpecies, network_.data().nameResidualSpeciesFile());
	    fResidualsSpecies.setf(std::ios::scientific);
	    WriteResidualsSpeciesLabels();

	    openOutputFileAndControl(fResidualsReactorStatistics, network_.data().nameResidualReactorStatisticsFile());
	    fResidualsReactorStatistics.setf(std::ios::scientific);
	    openOutputFileAndControl(fResidualsSpeciesStatistics, network_.data().nameResidualSpeciesStatisticsFile());
	    fResidualsSpeciesStatistics.setf(std::ios::scientific);
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

	    ratioNormInf	= new RingVector<double>(10);
	    ratioNorm1		= new RingVector<double>(10);
	    ratioNorm2		= new RingVector<double>(10);
	}
	    localIteration_  = 0;
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
	MPI::Status status;
	// Convenctional Analysis

    	React_omegamax = 0;
    	React_omegamin = 1;

	maxSpecies = 0;
	minSpecies = 1;

	{
	    int counter = 0;
	    double *value = new double[nResiduals_];
	    for (int k=1;k<=NR_P[procrank_];k++)
	    {
	    	network_.reactors(k).Residuals(network_);
		for(int i = 1; i <= network_.NumberOfSpecies(); i++)
		{
		    value[counter] = network_.reactors(k).residuals()[i];
		    counter++;

	    	    React_omegamax = std::max( network_.reactors(k).omegaMax(), React_omegamax);
	    	    React_omegamin = std::min( network_.reactors(k).omegaMin(), React_omegamin);
		}
	    }

	    VecSetValues(EqResiduals_, nResiduals_, res_place, value, INSERT_VALUES);
	    delete [] value;
	}

	VecAssemblyBegin(EqResiduals_);

    	for(int p = 0; p <= numworkers_; p++)
    	{
            MPI::COMM_WORLD.Barrier();
            if(p == procrank_)
            {
                mtype = FROM_WORKER;
                MPI::COMM_WORLD.Send(&React_omegamax, 1, MPI::DOUBLE, MASTER, mtype);
                MPI::COMM_WORLD.Send(&React_omegamin, 1, MPI::DOUBLE, MASTER, mtype);
            }
            if(procrank_ == 0)
            {
                source = p;
                mtype = FROM_WORKER;
                MPI::COMM_WORLD.Recv(&React_omegamax, 1, MPI::DOUBLE, source, mtype, status);
                MPI::COMM_WORLD.Recv(&React_omegamin, 1, MPI::DOUBLE, source, mtype, status);
                maxSpecies = std::max( React_omegamax, maxSpecies);
                minSpecies = std::min( React_omegamin, minSpecies);
            }
        }

	VecAssemblyEnd(EqResiduals_);
	
	VecNorm(EqResiduals_, NORM_INFINITY, &normInf_);
	VecNorm(EqResiduals_, NORM_1, &norm1_);	
	VecNorm(EqResiduals_, NORM_2, &norm2_);

	if(procrank_ == 0)
	{
	    ratioNormInf->Append(normInf_);
	    ratioNorm1->Append(norm1_);
	    ratioNorm2->Append(norm2_);

	    // Write on file
	    if(procrank_ == 0)	WriteResidualsOnFile();
	}

	// Update current status
	if ( currentStatus_ != network_.status() || network_.status()==KPP_NETWORK_STATUS_START )
	{
	    currentStatus_		= network_.status();
	    localIteration_		= 0;
	    startingNormInf_	= normInf_;
	    startingNorm1_		= norm1_;
	    startingNorm2_		= norm2_;
	}

	// Refined Analysis
	if ( (network_.iteration()-1)%nFrequencyRefinedAnalysis_ == 0 )
	{
	    AnalysisRefined();
	    if(procrank_ == 0) WriteResidualsSpeciesOnFile();
	    if(procrank_ == 0) Statistics();
	}

	    localIteration_++;	
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::AnalysisRefined()
{
	VecGetArray(EqResiduals_, &get_residual);
	
	//Species Residuals
	for(int j = 0; j < network_.NumberOfSpecies(); j++)
	{
	    double local_normInf, local_norm1, local_norm2;
	    int count = 0;
	    double *value = new double[network_.LocalNumberOfReactors()[procrank_]];
	    for(int k = 0; k < network_.LocalNumberOfReactors()[procrank_]; k++)
	    {
		value[k] = get_residual[j + count * network_.NumberOfSpecies()];
		count++;
	    }

	    VecSetValues(SpeciesResiduals_, nSpecies_, species_place, value, INSERT_VALUES);

	    VecAssemblyBegin(SpeciesResiduals_);
	    VecAssemblyEnd(SpeciesResiduals_);


	    VecNorm(SpeciesResiduals_, NORM_INFINITY, &local_normInf);
	    VecNorm(SpeciesResiduals_, NORM_1, &local_norm1);
	    VecNorm(SpeciesResiduals_, NORM_2, &local_norm2);

	    if(procrank_ == 0)
	    {
		speciesNormInf_[j+1] = local_normInf;
		speciesNorm1_[j+1] = local_norm1 / double(network_.NumberOfReactors());
		speciesNorm2_[j+1] = local_norm2 / double(network_.NumberOfReactors());
	    }
	    delete [] value;
	}

	// Reactor residuals
	int count = 0;
	for(int k = 1; k <= network_.LocalNumberOfReactors()[procrank_]; k++)
	{
	    int index = k+offset[procrank_];
	    BzzVector local_residuals(network_.NumberOfSpecies());
	    for(int i = 0; i < network_.NumberOfSpecies(); i++)
	    {
		local_residuals(i+1) = get_residual[count];
		count++;
	    }
	    
	    r_residualsMaxAbs[index] = local_residuals.MaxAbs();
	    r_residualsNorm1[index] = local_residuals.GetSumAbsElements()/double(network_.NumberOfSpecies());
	    r_residualsNorm2[index] = local_residuals.Norm2()/double(network_.NumberOfSpecies());

/*	    if(procrank_ + 1 + nprocs_ * (k-1) == network_.NumberOfReactors() - 1)
	    {
		cout << procrank_ + 1 + nprocs_ * (k-1) << std::endl;
		for(int i = 1; i <= network_.NumberOfSpecies(); i++)
		    std::cout << network_.reactors(k).residuals()[i] << std::endl;
	    }*/
		
	}

	communicator_->GatherArray(r_residualsMaxAbs, NR_P, offset);
	communicator_->GatherArray(r_residualsNorm1, NR_P, offset);
	communicator_->GatherArray(r_residualsNorm2, NR_P, offset);

	if(procrank_ == 0)
	{
	    for(int k = 1; k <= network_.NumberOfReactors();k++)
	    {
		reactorNormInf_[k] = r_residualsMaxAbs[k];
		reactorNorm1_[k] = r_residualsNorm1[k];
		reactorNorm2_[k] = r_residualsNorm2[k];
	    }
	}
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsOnVideo()
{
	std::cout << std::endl;
	std::cout << "// ********************************************************************************* // " << std::endl;
	std::cout << "//                               Residual Analysis                                   // " << std::endl;
	std::cout << "// ********************************************************************************* // " << std::endl;
	std::cout << "   #   normInf        norm1          norm2          minOmega       maxOmega             " << std::endl;

	std::cout << "       ";
	std::cout << std::setw(15) << std::left << std::scientific << normInf_;
	std::cout << std::setw(15) << std::left << std::scientific << norm1_/double( network_.NumberOfEquations() );
	std::cout << std::setw(15) << std::left << std::scientific << norm2_/double( network_.NumberOfEquations() );
	std::cout << std::setw(15) << std::left << std::scientific << minSpecies;
	std::cout << std::setw(15) << std::left << std::scientific << maxSpecies;
	std::cout << std::setw(15) << std::left << std::scientific << ratioNormInf->MeanRatios();
	std::cout << std::setw(15) << std::left << std::scientific << ratioNorm1->MeanRatios();
	std::cout << std::setw(15) << std::left << std::scientific << ratioNorm2->MeanRatios();
	std::cout << std::endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsOnFile()
{

	fResiduals << std::setw(9)  << std::left << network_.iteration();
	fResiduals << std::setw(9)  << std::left << localIteration_;
	fResiduals << std::setw(9)  << std::left << currentStatus_;
	fResiduals << std::setw(16) << std::left << normInf_;
	fResiduals << std::setw(16) << std::left << norm1_/double( network_.NumberOfEquations() );
	fResiduals << std::setw(16) << std::left << norm2_/double( network_.NumberOfEquations() );
	fResiduals << std::setw(16) << std::left << normInf_/startingNormInf_;
	fResiduals << std::setw(16) << std::left << norm1_/startingNorm1_;
	fResiduals << std::setw(16) << std::left << norm2_/startingNorm2_;
	fResiduals << std::setw(16) << std::left << norm1_;
	fResiduals << std::setw(16) << std::left << norm2_;
	fResiduals << std::endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsSpeciesOnFile()
{
	fResidualsSpecies << std::setw(9)  << std::left << network_.iteration();
	for (int j=1;j<=network_.NumberOfSpecies();j++)
		fResidualsSpecies << std::setw(16)  << std::left << speciesNorm2_[j];
	fResidualsSpecies << std::endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsLabels()
{
	fResiduals << std::setw(9)   << std::left << "Glob. #";
	fResiduals << std::setw(9)   << std::left << "Local #";
	fResiduals << std::setw(9)   << std::left << "Status";
	fResiduals << std::setw(16)  << std::left << "NormInf";
	fResiduals << std::setw(16)  << std::left << "MeanNorm1";
	fResiduals << std::setw(16)  << std::left << "MeanNorm2";
	fResiduals << std::setw(16)  << std::left << "NormInf(ratio)";
	fResiduals << std::setw(16)  << std::left << "Norm1(ratio)";
	fResiduals << std::setw(16)  << std::left << "Norm2(ratio)";
	fResiduals << std::setw(16)  << std::left << "Norm1";
	fResiduals << std::setw(16)  << std::left << "Norm2";
	fResiduals << std::endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsStatisticsLabels()
{
	fResidualsReactorStatistics << std::setw(9)  << std::left << "Iter.";
	for (int j=1;j<=nRanges_;j++)
		fResidualsReactorStatistics << std::setw(16)  << std::left << "range";
	fResidualsReactorStatistics << std::endl;

	fResidualsSpeciesStatistics << std::setw(9)  << std::left << "Iter.";
	for (int j=1;j<=nRanges_;j++)
		fResidualsSpeciesStatistics << std::setw(16)  << std::left << "range";
	fResidualsSpeciesStatistics << std::endl;
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::WriteResidualsSpeciesLabels()
{
	fResidualsSpecies << std::setw(9)  << std::left << "Iter.";
	
	int count = 2;
	for (int j=1;j<=network_.NumberOfSpecies();j++)
	{
		std::stringstream index; index << count++;
		std::string label = network_.mix()[0].names[j] + "(" + index.str() + ")";
		
		fResidualsSpecies << std::setw(16)  << std::left << label;
	}
	fResidualsSpecies << std::endl;
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

		fResidualsReactorStatistics << std::setw(9)  << std::left << network_.iteration();
		for (int i=1;i<=ranges_.Size()+1;i++)
			fResidualsReactorStatistics << std::setw(16)  << std::left << 100.*double(countRanges[i])/double(network_.NumberOfReactors());
		fResidualsReactorStatistics << std::endl;
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

		fResidualsSpeciesStatistics << std::setw(9)  << std::left << network_.iteration();
		for (int i=1;i<=ranges_.Size()+1;i++)
			fResidualsSpeciesStatistics << std::setw(16)  << std::left << 100.*double(countRanges[i])/double(network_.NumberOfSpecies());
		fResidualsSpeciesStatistics << std::endl;
	}
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::InitializeResidualsVector()
{
	nResiduals_ = network_.LocalNumberOfReactors()[procrank_] * network_.NumberOfSpecies();

	communicator_->InitializePetscVector(EqResiduals_, nResiduals_, residuals_low, residuals_high);

	res_place = new Petsc64bitInt[nResiduals_];

	std::vector<int*> res_place_vec;
	res_place_vec.resize(NR_P[procrank_] + 1);

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int size = network_.NumberOfSpecies();
	    res_place_vec[k] = new int[size];
	}


	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = k + offset[procrank_];
	    for(int i = 0; i < network_.NumberOfSpecies(); i++)
	        res_place_vec[k][i] = (index - 1) * network_.NumberOfSpecies() + i;
	}

	int counter = 0;
	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    for(int i = 0; i < network_.NumberOfSpecies(); i++)
	    {
		res_place[counter] = res_place_vec[k][i];
		counter++;
	    }
	}

	for(int k = 1; k <= NR_P[procrank_]; k++)
	    delete [] res_place_vec[k];
}

void OpenSMOKE_KPP_ReactorNetwork_Residuals::InitializeSpeciesResiduals()
{
	nSpecies_ = network_.LocalNumberOfReactors()[procrank_];

	communicator_->InitializePetscVector(SpeciesResiduals_, nSpecies_, speciesresiduals_low, speciesresiduals_high);

	species_place = new Petsc64bitInt[NR_P[procrank_]];

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = k + offset[procrank_];
	    species_place[k-1] = index - 1;
	}
}
