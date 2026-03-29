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
#include "OpenSMOKE_KPP_Communicator.h"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <mpi.h>

void OpenSMOKE_KPP_SingleReactorStatistics::ErrorMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class: OpenSMOKE_KPP_SingleReactorStatistics"	<< std::endl;
    std::cout << "Error: " << message_				<< std::endl;
    std::cout << "Press enter to continue... "		<< std::endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_SingleReactorStatistics::WarningMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class:	  OpenSMOKE_KPP_SingleReactorStatistics"		<< std::endl;
    std::cout << "Warning: " << message_					<< std::endl;
	std::cout << std::endl;
}

OpenSMOKE_KPP_SingleReactorStatistics::OpenSMOKE_KPP_SingleReactorStatistics(const int numberOfSpecies, const int numberOfReactors, const std::string fileName, OpenSMOKE_KPP_Communicator* communicator) :
communicator_(communicator)
{
	iteration_ = 0;

	nprocs_ = MPI::COMM_WORLD.Get_size();
	procrank_ = MPI::COMM_WORLD.Get_rank();
	numworkers_ = nprocs_ - 1;

	MASTER = 0;
	FROM_MASTER = 1;
	FROM_WORKER = 2;

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

	if(procrank_ == 0)
	{
	    communicator_->InitializeArray(numberOfStepsGlob_, nprocs_);
	    communicator_->InitializeArray(numberOfFunctionsGlob_, nprocs_);
	    communicator_->InitializeArray(numberOfFunctionsForJacobianGlob_, nprocs_);
	    communicator_->InitializeArray(numberOfAnalyticalJacobiansGlob_, nprocs_);
	    communicator_->InitializeArray(numberOfNumericalJacobiansGlob_, nprocs_);
	    communicator_->InitializeArray(numberOfFactorizationsGlob_, nprocs_);

	    communicator_->InitializeArray(maxAbsResidualGlob_, nprocs_);
	    communicator_->InitializeArray(meanAbsResidualGlob_, nprocs_);
	}

	ChangeDimensions( numberOfSpecies_, &residuals_);

	if(procrank_ == 0)
	{
	    fOutput.open(fileName.c_str(), std::ios::out);
	    fOutput.setf(std::ios::scientific);
	}
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

	{	
	    residuals_ = o.GetY1InMeshPoint();
	
	    numberOfSteps_ += o.GetNumStep();
	    numberOfFunctions_ += o.GetNumFunction();
	    numberOfFunctionsForJacobian_ += o.GetNumFunctionForJacobian();
	    numberOfAnalyticalJacobians_ += o.GetNumAnalyticalJacobian();
	    numberOfNumericalJacobians_ += o.GetNumNumericalJacobian();
	    numberOfFactorizations_ += o.GetNumFactorization();

	    maxAbsResidual_   = std::max(maxAbsResidual_, residuals_.MaxAbs());
	    meanAbsResidual_ += residuals_.GetSumAbsElements()/double(numberOfSpecies_);
	}
}

void OpenSMOKE_KPP_SingleReactorStatistics::PrintOnFile()
{
	MPI::Status status;

	for(int p = 0; p <= numworkers_; p++)
    	{
            mtype = FROM_WORKER;
            MPI::COMM_WORLD.Barrier();
            if(p == procrank_)
            {
                MPI::COMM_WORLD.Send(&numberOfSteps_, 1, MPI::INT, MASTER, mtype);
                MPI::COMM_WORLD.Send(&numberOfFunctions_, 1, MPI::INT, MASTER, mtype);
                MPI::COMM_WORLD.Send(&numberOfFunctionsForJacobian_, 1, MPI::INT, MASTER, mtype);
                MPI::COMM_WORLD.Send(&numberOfAnalyticalJacobians_, 1, MPI::INT, MASTER, mtype);
                MPI::COMM_WORLD.Send(&numberOfNumericalJacobians_, 1, MPI::INT, MASTER, mtype);
                MPI::COMM_WORLD.Send(&numberOfFactorizations_, 1, MPI::INT, MASTER, mtype);

                MPI::COMM_WORLD.Send(&maxAbsResidual_, 1, MPI::DOUBLE, MASTER, mtype);
                MPI::COMM_WORLD.Send(&meanAbsResidual_, 1, MPI::DOUBLE, MASTER, mtype);
            }
            if(procrank_ == 0)
            {
                source = p;
                MPI::COMM_WORLD.Recv(&numberOfStepsGlob_[source], 1, MPI::INT, source, mtype, status);
                MPI::COMM_WORLD.Recv(&numberOfFunctionsGlob_[source], 1, MPI::INT, source, mtype, status);
                MPI::COMM_WORLD.Recv(&numberOfFunctionsForJacobianGlob_[source], 1, MPI::INT, source, mtype, status);
                MPI::COMM_WORLD.Recv(&numberOfAnalyticalJacobiansGlob_[source], 1, MPI::INT, source, mtype, status);
                MPI::COMM_WORLD.Recv(&numberOfNumericalJacobiansGlob_[source], 1, MPI::INT, source, mtype, status);
                MPI::COMM_WORLD.Recv(&numberOfFactorizationsGlob_[source], 1, MPI::INT, source, mtype, status);

                MPI::COMM_WORLD.Recv(&maxAbsResidualGlob_[source], 1, MPI::DOUBLE, source, mtype, status);
                MPI::COMM_WORLD.Recv(&meanAbsResidualGlob_[source], 1, MPI::DOUBLE, source, mtype, status);
            }
        }

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= numworkers_; i++)
	    {
		numberOfStepsGlob_[0] += numberOfStepsGlob_[i];
		numberOfFunctionsGlob_[0] += numberOfFunctionsGlob_[i];
		numberOfFunctionsForJacobianGlob_[0] += numberOfFunctionsForJacobianGlob_[i];
		numberOfAnalyticalJacobiansGlob_[0] += numberOfAnalyticalJacobiansGlob_[i];
		numberOfNumericalJacobiansGlob_[0] += numberOfNumericalJacobiansGlob_[i];
		numberOfFactorizationsGlob_[0] += numberOfFactorizationsGlob_[i];

		maxAbsResidualGlob_[0] = std::max(maxAbsResidualGlob_[0], maxAbsResidualGlob_[i]);
		meanAbsResidualGlob_[0] += meanAbsResidualGlob_[i];
	    }

	    iteration_++;

	    fOutput << std::setw(16) << std::left << iteration_;
	    fOutput << std::setw(16) << std::left << meanAbsResidualGlob_[0]/double(numberOfReactors_);
	    fOutput << std::setw(16) << std::left << maxAbsResidualGlob_[0];
	    fOutput << std::setw(16) << std::left << numberOfStepsGlob_[0]/double(numberOfReactors_);
	    fOutput << std::setw(16) << std::left << numberOfFunctionsGlob_[0]/double(numberOfReactors_);
	    fOutput << std::setw(16) << std::left << numberOfFunctionsForJacobianGlob_[0]/double(numberOfReactors_);
	    fOutput << std::setw(16) << std::left << numberOfAnalyticalJacobiansGlob_[0]/double(numberOfReactors_);
	    fOutput << std::setw(16) << std::left << numberOfNumericalJacobiansGlob_[0]/double(numberOfReactors_);
	    fOutput << std::setw(16) << std::left << numberOfFactorizationsGlob_[0]/double(numberOfReactors_);
	    fOutput << std::endl;
	}
}
