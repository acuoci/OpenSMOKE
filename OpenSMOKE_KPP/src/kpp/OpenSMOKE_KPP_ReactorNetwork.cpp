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

#include "mkl.h"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <mpi.h>
#include <petscvec.h>
#include <sstream>
#include <math.h>
#include <vector>
#include <fstream>
#include "linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.h"
//#include "linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.h"
#include "linear_solvers/OpenSMOKE_LIS_Unsymmetric.h"
#include "OpenSMOKE_KPP_BlockMatrixNetwork.h"
#include "OpenSMOKE_KPP_ReactorNetwork.h"
#include "OpenSMOKE_KPP_DataManager.h"
#include "OpenSMOKE_KPP_ConvectiveNetworkStatistics.h"
#include "OpenSMOKE_KPP_SingleReactorStatistics.h"
#include "OpenSMOKE_KPP_ReactorNetwork_Residuals.h"
#include "OpenSMOKE_KPP_NewtonMethod_Manager.h"
#include "OpenSMOKE_KPP_ODE_Manager.h"
#include "OpenSMOKE_KPP_Communicator.h"
#include "OpenSMOKE_KPP_SingleReactor_KineticsManager.h"

OpenSMOKE_KPP_ReactorNetwork *ptNetwork;

void externalResiduals(BzzVector &x, Vec &f)
{
	ptNetwork->Residuals(x,f);
}

int externalSolutionNLS(BzzVector &x, BzzVector &direction, const bool jacobianFlag)
{
	return ptNetwork->SolutionNLS(x, direction, jacobianFlag);
}

int externalSolutionODE(BzzVector &x, BzzVector &step, double& dt, const bool jacobianFlag)
{
	return ptNetwork->SolutionODE(x, step, dt, jacobianFlag);
}

void externalResidualsAnalysis()
{
	ptNetwork->ResidualsAnalysis();
}

void OpenSMOKE_KPP_ReactorNetwork::ErrorMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class: OpenSMOKE_KPP_ReactorNetwork"	<< std::endl;
    std::cout << "Error: " << message_				<< std::endl;
    std::cout << "Press enter to continue... "		<< std::endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_ReactorNetwork::WarningMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class:	  OpenSMOKE_KPP_ReactorNetwork"		<< std::endl;
    std::cout << "Warning: " << message_					<< std::endl;
	std::cout << std::endl;
}

OpenSMOKE_KPP_ReactorNetwork::OpenSMOKE_KPP_ReactorNetwork(OpenSMOKE_ReactingGas* mix, OpenSMOKE_KPP_DataManager& data, std::ofstream& fLog, std::ofstream& fWarning) :
mix_(mix), data_(data), fLog_(fLog), fWarning_(fWarning)
{
	ptNetwork	= this;

	nprocs_ = MPI::COMM_WORLD.Get_size();
	procrank_ = MPI::COMM_WORLD.Get_rank();
	numworkers_ = nprocs_ - 1;

	offset = new int[nprocs_ + 1];
	NR_P = new int[nprocs_];
	
	MASTER = 0;
	FROM_MASTER = 1;
	FROM_WORKER = 2;

	iteration_  	 = 0;
	globalTime_ 	 = 0.;
	deltat		 = 0.;
	LoopIndex_	 = 0;

	iAlreadyGlobalPardiso_ = false;
	iAlreadyGlobalMumps_ = false;
	iAlreadyGlobalLis_ = false;
	iAlreadyGlobalGaussSiedel_ = false;
        iBooleanMatrix_ = false;
	
	newtonconvergence_ = false;
	CSTRConvergengeAlert = false;

	nLocalCSTRSequence_  = 0;
	nGlobalCSTRSequence_ = 0;
        
	
	indexGlobalNLS_ = 0;
        globalCounterNLS_ = 0;
        localCounterNLS_ = 0;
        globalCounterNLSJacobians_ = 0;
        
        indexGlobalODE_ = 0;
        globalCounterODE_ = 0;
        localCounterODE_ = 0;
        globalCounterODEJacobians_ = 0;     
        
	lastODEStep_ = data_.GlobalODE_InitialTimeStep();
	currentDeltaTimeGlobalODE_ = data_.GlobalODE_InitialTimeStep()*5.;

}

OpenSMOKE_KPP_ReactorNetwork::~OpenSMOKE_KPP_ReactorNetwork(void)
{
}
void OpenSMOKE_KPP_ReactorNetwork::SetCommunicator(OpenSMOKE_KPP_Communicator* communicator)
{
	communicator_ = communicator;
}

void OpenSMOKE_KPP_ReactorNetwork::PreparingIndicesReactorsThreads()
{
	// Preparing
	indicesReactorsThreads_ = new BzzVectorInt[data_.nThreads()];

	// Policy 1
	int nBlock = NR_P[procrank_]/data_.nThreads();
	for (int k=0;k<data_.nThreads();k++)
            ChangeDimensions(nBlock, &indicesReactorsThreads_[k]);
        
       
        BzzVectorInt icount(data_.nThreads());
        icount = 1;
        int reactor = 1;

        while (reactor<=nBlock*data_.nThreads())
        {
            for (int k=0;k<data_.nThreads();k++)
                indicesReactorsThreads_[k][icount[k+1]++] = reactor++;
        }

        while(reactor<=NR_P[procrank_])
      	{
            for (int k=0;k<data_.nThreads();k++)
            {
                indicesReactorsThreads_[k].Append(reactor++);
                if (reactor>NR_P[procrank_]) break;
            }
        }

	// Policy 2
/*	int avereactors = NR_P[procrank_]/data_.nThreads();
    	int extra = NR_P[procrank_]%data_.nThreads();

	int shared_offset[data_.nThreads() + 1];
	int nBlock[data_.nThreads() + 1];
	shared_offset[0] = 0;
   	shared_offset[nprocs_] = NR_P[procrank_];
	
  	for(int dest = 0; dest < data_.nThreads(); dest++)
  	{
   	    nBlock[dest] = (dest < extra) ? avereactors + 1 : avereactors;
   	    if(dest != data_.nThreads())
   	         shared_offset[dest+1] = shared_offset[dest] + nBlock[dest];
  	}*/
            
        // Check
        int sum = 0;
        for (int k=0;k<data_.nThreads();k++)
            sum += indicesReactorsThreads_[k].Size();
        if (sum != NR_P[procrank_])
            ErrorMessage("Wrong divide and conquer policy!");
}

void OpenSMOKE_KPP_ReactorNetwork::MemoryAllocation()
{
	if(procrank_ == 0)
	{
	    ChangeDimensions(NR_, NR_, &C_);

	    ChangeDimensions(NR_*mix_[0].NumberOfSpecies(), &lastGlobalODESolution_);
	    ChangeDimensions(NR_*mix_[0].NumberOfSpecies(), &lastGlobalNLSSolution_);
	    ChangeDimensions(NR_*mix_[0].NumberOfSpecies(), &auxVector_NRxNC_1);
	}

	    residuals_ = new OpenSMOKE_KPP_ReactorNetwork_Residuals(*this, communicator_);
}

bool IsANumber(const char value)
{
    if (value == '0' || value == '1' || value == '2' || value == '3' || value == '4' || value == '5'
            || value == '6' || value == '7' || value == '8' || value == '9')
        return true;
    else
        return false;
}
void OpenSMOKE_KPP_ReactorNetwork::ReadFirstGuess()
{
	std::ifstream fInput(data_.nameFirstGuessFile().c_str(), std::ios::in);

	unsigned int NCDetailed;
	unsigned int NCReduced;

	fInput >> NCDetailed;
	fInput >> NCReduced;

	ChangeDimensions(NCReduced, &jReduced);
	std::vector<std::string> namesReduced;
        namesReduced.resize(NCReduced+1);

        for (int j=1;j<=jReduced.Size();j++)
            fInput >> namesReduced[j];
        
        if (IsANumber(namesReduced[1].at(0)) == true)
        {
            for (int j=1;j<=jReduced.Size();j++)
                jReduced[j] = atoi(namesReduced[j].c_str());

        }
        else
        {
            for (int j=1;j<=jReduced.Size();j++)
                jReduced[j] = mix_[0].recognize_species(namesReduced[j]);
        }
        
	if(procrank_ == 0)
	{
            for (int j=1;j<=jReduced.Size();j++)
            	std::cout << jReduced[j] << " " << mix_[0].names[jReduced[j]] << std::endl;
	}
        
        
	fInput >> NR_;

//	if(NR_ > 60000)		ErrorMessage("This version supports a std::maximum of 60000 cells");

	DistributeCells();


	// Allocating kinetics
//	if (data_.iSaveKineticConstants() == true)
	{
		kinetics_ = new OpenSMOKE_KPP_SingleReactor_KineticsManager[NR_P[procrank_]+1];
	}
        
        // Indices threads
        PreparingIndicesReactorsThreads();

	// Allocating reactors
	reactors_ = new OpenSMOKE_KPP_SingleReactor[NR_P[procrank_]+1];


	    {
		for(int tid = 0; tid < data_.nThreads(); tid++)
		{
		    for(int kk = 1; kk <= indicesReactorsThreads_[tid].Size(); kk++)
		    {
			int k = indicesReactorsThreads_[tid][kk];
			reactors_[k].Setup(k, k + offset[procrank_], &mix_[tid], jReduced);
		    }
		}
	    }

	for (int k=1;k<=NR_P[procrank_];k++)
	    reactors_[k].SetDataManager(&data_);

		
	ChangeDimensions(NR_P[procrank_], NCDetailed, &residualMatrix_);

	unsigned int dummy;
	fInput >> dummy;
	if (dummy != NCReduced)
		ErrorMessage("The number of reduced species in the FirstGuess file does not fit!");
	if (mix_[0].NumberOfSpecies() != NCDetailed)
		ErrorMessage("The number of detailed species in the FirstGuess file does not fit with the number of species in the detailed kinetic file!");


	MPI::COMM_WORLD.Barrier();

	for(int i = 0; i <= numworkers_; i++)
	{
	    if(i == procrank_)
	    {
		double dummy_;
		for(int p = 1; p <= NCReduced*offset[procrank_]; p++)
		{
                    fInput >> dummy_;
                }
		for (int k = 1; k<=NR_P[procrank_] ;k++)
                {
		    if ( (k + offset[procrank_])%(NR_/10) == 0) std::cout << "Reading mass fractions: " << k + offset[procrank_] << "/" << NR_ << std::endl;
                    reactors_[k].ReadMassFractions(fInput);
		}
	    }
	}

	if (data_.iBackup() == true)
	    ReadBackupFile(data_.nameInputBackupFile());		//Check at the end if it's done well!!!

	fInput.close();

	communicator_->InitializeParameters();

}

void OpenSMOKE_KPP_ReactorNetwork::ReadTopology()
{
	std::ifstream fInput(data_.nameTopologyFile().c_str(), std::ios::in);

	unsigned int dummy;
	fInput >> dummy;
	if (dummy != NR_)
		ErrorMessage("The number of reactors does not fit!");
	
        TopologyFakeReading();

	{
	    double dummy_;

	    for(int i = 1; i <= counter[procrank_]; i++) fInput >> dummy_;

            for (int k=1; k<=NR_P[procrank_]; k++)
            {
                if ( (k + offset[procrank_])%(NR_/10) == 0)	std::cout << "Reading reactor topology: " << k + offset[procrank_] << "/" << NR_ << std::endl;
                reactors_[k].ReadReactorProperties(fInput);
                reactors_[k].ReadTopology(fInput);
            }
	}	for (int k=1;k<=NR_P[procrank_];k++)

	//Min and Max
	Tmin = new double[nprocs_];
	Tmax = new double[nprocs_];
	for(int i = 0; i < nprocs_; i++)
	{
	    Tmin[i] = 10000;
	    Tmax[i] = 0;
	}

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    if(reactors_[k].temperature() > Tmax[procrank_]) Tmax[procrank_] = reactors_[k].temperature();
	    else if(reactors_[k].temperature() < Tmin[procrank_]) Tmin[procrank_] = reactors_[k].temperature();
	}

	for(int p = 0; p < nprocs_; p++)
	{
	    MPI::COMM_WORLD.Bcast(&Tmax[p], 1, MPI::DOUBLE, p);
	    MPI::COMM_WORLD.Bcast(&Tmin[p], 1, MPI::DOUBLE, p);
	}

	Tmin_all = 10000;
	Tmax_all = 0;

	for(int i = 0; i < nprocs_; i++)
	{
	    if(Tmin_all > Tmin[i])  Tmin_all = Tmin[i];
	    if(Tmax_all < Tmax[i])  Tmax_all = Tmax[i];
	}

	if(data_.TemperatureMin() == 0.)	data_.SetMinimumTemperature(Tmin_all - 20., "K");
	if(data_.TemperatureMax() == 0.)	data_.SetMaximumTemperature(Tmax_all * 1.07, "K");

	if(procrank_ == 0)
	{
	    std::cout << "Minimum temperature:			" << Tmin_all << " K" << std::endl;
	    std::cout << "Maximum temperature:			" << Tmax_all << " K" << std::endl;
	    std::cout << "Minimum temperature (effective):	" << data_.TemperatureMin() << " K" << std::endl;
	    std::cout << "Maximum temperature (effective):	" << data_.TemperatureMax() << " K" << std::endl;
	}


	if (data_.correction() != KPP_CORRECTION_NONE)
        {
            if (data_.TemperatureMin() >= Tmin_all-10)
                ErrorMessage("The imposed std::minimum temperature is too large! At least 10 K below than std::minimum!");

            if ( data_.TemperatureMax() <= Tmax_all*1.02)
                ErrorMessage("The imposed std::minimum temperature is too small! At least 2 percent higher than std::maximum!");
        }

	if (data_.iSaveKineticConstants() == true)
	{
		if (procrank_ == 0)	std::cout << "Evaluating kinetic constants for all the reactors in the network..." << std::endl;

		if (data_.iSymbolicKinetics() == false)
		{
		int tid;
		#pragma omp parallel private(tid)
		{
		    tid = omp_get_thread_num();
		    // Fluctuating species
		    std::vector<bool> fluctuatingReactions;

			std::stringstream tid_; tid_ << tid;
			std::stringstream omp_get_num_thread_; omp_get_num_thread_ << omp_get_num_threads();
			std::stringstream indicesReactorsThreads_str; indicesReactorsThreads_str << indicesReactorsThreads_[tid].Size();

//		    string msg =  tid_.str() + " " + omp_get_num_thread_.str() + " " + indicesReactorsThreads_str.str();
//			std::cout << msg << std::endl; 

		    mix_[tid].kinetics.ReactionIndices(data_.fluctuatingSpecies(), fluctuatingReactions);

		    // Kinetic constants
			BzzVector temperatures(indicesReactorsThreads_[tid].Size());
			BzzVector pressures_Pa(indicesReactorsThreads_[tid].Size());
			BzzVector concentrations(indicesReactorsThreads_[tid].Size());

//		std::cout << "Qui arrivo...A" << std::endl;

			for (int kk=1;kk<=indicesReactorsThreads_[tid].Size();kk++)
			{
				int k = indicesReactorsThreads_[tid][kk];
				reactors_[k].index_omp_ = kk;
				temperatures[kk] = reactors_[k].temperature();
				pressures_Pa[kk] = reactors_[k].pressure();
				concentrations[kk] = reactors_[k].pressure()/(Constants::R_J_kmol*reactors_[k].temperature());
			}
//		cout << "Qui arrivo...B" << std::endl;
			mix_[tid].InitializeMap(indicesReactorsThreads_[tid].Size());
			mix_[tid].ComputeKineticParameters_map(temperatures, concentrations, pressures_Pa);
//		cout << "Qui arrivo...B1" << std::endl;
			BzzVector correction_k1(mix_[tid].NumberOfReactions());
			BzzVector correction_k2(mix_[tid].NumberOfReactions());
			BzzVector correction_uKeq(mix_[tid].NumberOfReactions());
//		cout << "Qui arrivo...B2" << std::endl;
			BzzMatrix matrix_correction_k1(indicesReactorsThreads_[tid].Size(),   mix_[tid].NumberOfReactions() );
			BzzMatrix matrix_correction_k2(indicesReactorsThreads_[tid].Size(),   mix_[tid].NumberOfReactions() );
			BzzMatrix matrix_correction_uKeq(indicesReactorsThreads_[tid].Size(), mix_[tid].NumberOfReactions() );
	//	cout << "Qui arrivo...C" << std::endl;
			double v = -1.;
			for (int kk=1;kk<=indicesReactorsThreads_[tid].Size();kk++)
			{
	
				int k = indicesReactorsThreads_[tid][kk];
				if (data_.correction() != KPP_CORRECTION_NONE)
					v = reactors_[k].variance();

				kinetics_[k].SetMinMax(data_.TemperatureMin(), data_.TemperatureMax());
				kinetics_[k].Setup(data_.correction(), &mix_[tid], reactors_[k].temperature(), reactors_[k].pressure(), v,data_, fluctuatingReactions, correction_k1, correction_uKeq, correction_k2 );
				reactors_[k].SetKinetics(&kinetics_[k]);

				matrix_correction_k1.SetRow(kk, correction_k1);
				matrix_correction_k2.SetRow(kk, correction_k2);
				matrix_correction_uKeq.SetRow(kk, correction_uKeq);
			}


//		cout << "Qui arrivo...D" << std::endl;
			//if(procrank_ == 0)	std::cout << "Correcting the kinetic constants for all the reactors in the network..." << std::endl;
			mix_[tid].CorrectKineticParameters_map(matrix_correction_k1, matrix_correction_uKeq, matrix_correction_k2);

		      }
		    }
		    if (data_.iSymbolicKinetics() == true)
		    {
			// Fluctuating species
		        std::vector<bool> fluctuatingReactions;


		        mix_[0].kinetics.ReactionIndices(data_.fluctuatingSpecies(), fluctuatingReactions);

			BzzVector correction_k1(mix_[0].NumberOfReactions());
			BzzVector correction_k2(mix_[0].NumberOfReactions());
			BzzVector correction_uKeq(mix_[0].NumberOfReactions());

			double v = -1.;
			for (int k=1;k<=NR_P[procrank_];k++)
			{
				if (data_.correction() != KPP_CORRECTION_NONE)
					v = reactors_[k].variance();

				kinetics_[k].SetMinMax(data_.TemperatureMin(), data_.TemperatureMax());
				kinetics_[k].Setup(data_.correction(), &mix_[0], reactors_[k].temperature(), reactors_[k].pressure(), v,data_, fluctuatingReactions, correction_k1, correction_uKeq, correction_k2 );
				reactors_[k].SetKinetics(&kinetics_[k]);
			}

			if(procrank_ == 0)	std::cout << "Correcting the kinetic constants for all the reactors in the network..." << std::endl;
		     }
	     }

	   else if (data_.iSaveKineticConstants() == false)
	  {
	     ErrorMessage("Save kinetic constants not implemented in the parallel version");
	  }

	fInput.close();

	MPI::COMM_WORLD.Barrier();
}

void OpenSMOKE_KPP_ReactorNetwork::BuildNetwork()
{
	MPI::COMM_WORLD.Barrier();
	// Memory Allocation
	if(procrank_ == 0)
	    std::cout << " * Memory allocation..." << std::endl;

	MemoryAllocation();

	if(procrank_ == 0)
            std::cout << "Done: Memory Allocation" << std::endl; //getchar();

	MPI::COMM_WORLD.Barrier();

	// Build network
	if(procrank_ == 0)
	    std::cout << " * Assembling mass flow rate topology..." << std::endl;

	TopologyParametersSize();
	TopologyParameters();

        {
            for(int k = 1; k <= NR_; k++)
            {
                for(int j = 1; j <= ReactOutSize[k]; j++)
                {
                    if(ReactID[React_out[k][j]] == procrank_)
                    {
                        reactors_[React_out[k][j] - offset[procrank_]].UpdateInflow(k, React_cOut[k][j], React_diffusion[k][j]);
                    }
                }
            }
    	}

	if(procrank_ == 0)
    	    std::cout << "Done: Assembling mass flow rate topology" << std::endl; // getchar();
        
	// Assembling local reactors
	if(procrank_ == 0)
	    std::cout << " * Assembling local topology..." << std::endl;

	for (int k=1;k<=NR_P[procrank_];k++)
	    reactors_[k].Assembling();
	
	if(procrank_ == 0)
	    std::cout << "Done: Assembling local topology" << std::endl; //getchar();

	DeleteTopologyParametersSize();
	TopologyParametersSize();
	TopologyParameters();



	// Mass flow rate Errors
	if(procrank_ == 0)
	    std::cout << " * Mass flow errors (before corrections)" << std::endl;

	MassFlowErrors();

	if(procrank_ == 0)
	    std::cout << "Done: Mass flow errors" << std::endl; //getchar();

	// Adjust mass flows
	if(procrank_ == 0)
	    std::cout << " * Correcting mass flow rates..." << std::endl;

	CorrectingMassFlowRates();

	if(procrank_ == 0)
    	    std::cout << "Done: Correcting mass flow rates" << std::endl; //getchar();

	// Mass flow rate Errors
	if(procrank_ == 0)
	    std::cout << " * Mass flow errors (after corrections)" << std::endl;

	MassFlowErrors();

	if(procrank_ == 0)
    	    std::cout << "Done: Mass flow errors" << std::endl; //getchar();

	// External feeds and outputs indices
	FeedsAndOutputs();

	if(procrank_ == 0)
	    std::cout << "Done: External feeds and outputs indices" << std::endl; //getchar();

	VariableTopologyParameters();


	// Convection-Diffusion Matrix
	if(procrank_ == 0)
	    std::cout << " * Assembling Convection-Diffusion Matrix..." << std::endl;

	AssemblingConvectionDiffusionMatrix();
	if(procrank_ == 0)
	    std::cout << "Done: Convection-Diffusion Matrix" << std::endl; //getchar();

	// External feeds Right Hand Sides
	if(procrank_ == 0)
	    std::cout << " * Assembling Right Hand Sides..." << std::endl;

	AssemblingRightHandSides();

	if(procrank_ == 0)
	    std::cout << "Done: Assembling Right Hand Sides" << std::endl; //getchar();

	// Exchange info with single reactors
	if(procrank_ == 0)
	    std::cout << " * Convection-Diffusion matrix communication..." << std::endl;

	InitializeNetworkData();


	for(int k = 1; k <= NR_; k++)
	{
	    network_index_ = k;
	    network_ID_ = ReactID[k];

	    PrepareCDMatrixCommunication();

	    if(procrank_ == ReactID[k])
	    {
            reactors_[k - offset[ReactID[k]]].ConvectionDiffusionMatrixCommunication(*this);
	    }
	}

	DeleteTopologyParametersSize();
	TopologyParametersSize();
	TopologyParameters();
	
	if(procrank_ == 0)
	    std::cout << "Done: Convection-Diffusion matrix communication" << std::endl; //getchar();
        
	// Statistics on C matrix
	SummaryConvectionDiffusionMatrix();

	if(procrank_ == 0)
    	    std::cout << "Done: SummaryConvectionDiffusionMatrix" << std::endl; //getchar();

	// Write additional info
	WriteReactorNetworkData();

	if(procrank_ == 0)
    	    std::cout << "Done: WriteReactorNetworkData" << std::endl; //getchar();

	//Organize Convection Diffusion Parameters in a Vector Form
	ConvDiffParameters();

	SendAndBroadcastOmega();

	if(procrank_ == 0)
	    std::cout << "First Residuals Analysis..." << std::endl;
	// Initial residuals

	ResidualsAnalysis();

	if(procrank_ == 0)
	    std::cout << "Done: Residuals Analysis" << std::endl; //getchar();

	// Initialize OdeSystem
	if(procrank_ == 0)
	    std::cout << " * Initialize reactors..." << std::endl;
	

    	    odePool = new ODE_Pool(data_.nThreads(), NumberOfSpecies());
    	    for(int k=0;k<data_.nThreads();k++)
		odePool->Set(k, reactors_[k+1]);
	    odePool->Initialize();

//	for (int k=1;k<=NR_;k++)
//		reactors_[k].SetInitialConditions();

	if(procrank_ == 0)
    	    std::cout << "Done: Initialize reactors" << std::endl; //getchar();
	// Open Files
	if(procrank_ == 0)
	    std::cout << " * Open statistics files..." << std::endl;

	statistics_ = new OpenSMOKE_KPP_SingleReactorStatistics(mix_[0].NumberOfSpecies(), NR_, "SequenceStatistics.out", communicator_);

	if(procrank_ == 0)
	{
	    statisticsConvective_ = new OpenSMOKE_KPP_ConvectiveNetworkStatistics(mix_[0].NumberOfSpecies(), NR_, "SequenceConvection.out");

	    std::cout << "Done: Open statistics files" << std::endl; //getchar();
	}
	
	if (data_.GlobalODE_SparseLinearSolver()  == KPP_SPARSESOLVER_PARDISO || data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
	{
		if(procrank_ == 0)
		    ErrorMessage("Pardiso not implemented in the parallel version!");
		// PARDISO Global Linear Systems 
//		pardisoMatrixGlobal = new OpenSMOKE_PARDISO_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW);
//		pardisoMatrixGlobal->SetFillInReducingOrdering(data_.FillInReducingOrdering());
	}
	
	if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS || data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
	{
		if(procrank_ == 0)
		    ErrorMessage("MUMPS not implemented in the parallel version!");
		// MUMPS Global Linear Systems
		//mumpsMatrixGlobal = new OpenSMOKE_MUMPS_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW);
	}
	
	if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_LIS || data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
	{
		// LIS Global Linear Systems
		lisMatrixGlobal = new OpenSMOKE_LIS_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW);
		lisMatrixGlobal->SetData(&data_);
		lisMatrixGlobal->SetCommunicator(communicator_);
	}

	if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL || data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL)
	{
		if(procrank_ == 0)
		    ErrorMessage("Gauss Seidel not implemented in the parallel version!");
		// OpenSMOKE Gauss Siedel Global Linear Systems
		openSMOKEMatrixGlobal = new OpenSMOKE_KPP_BlockMatrixNetwork();
	}

	if(procrank_ == 0)
	    std::cout << "Done: Global Solver" << std::endl; //getchar();
	
	// Courant number (the std::maximum value is 1)
	if(procrank_ == 0)
	    ChangeDimensions(NR_, &timeSteps);
	
	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	    {
		int index = k + offset[procrank_];
	        timeStepsArray[index] = reactors_[k].mass()/reactors_[k].M() * 1.;
	    }
	}

	communicator_->GatherArray(timeStepsArray, NR_P, offset);

	if(procrank_ == 0)
	{
	    for(int k = 1; k <= NR_; k++)
		timeSteps[k] = timeStepsArray[k];

	    std::cout << "Time step analysis... " << std::endl;
	    std::cout << " * Maximum time step: " << timeSteps.Max() << " s" << std::endl;
	    std::cout << " * Minimum time step: " << timeSteps.Min() << " s" << std::endl;
	    std::cout << " * Mean time step:    " << Mean(timeSteps) << " s" << std::endl;
	    std::cout << std::endl;

	    timeSteps *= data_.PredictorCorrector_CourantCorrectionCoefficient();

	    // Open Output Files
	    std::string sequence_name  = data_.nameFolderOutput() + "/Residuals/Sequence.out";
	    std::string globalnls_name = data_.nameFolderOutput() + "/Residuals/GlobalNLS.out";
	    std::string globalode_name = data_.nameFolderOutput() + "/Residuals/GlobalODE.out";
	

	
	    fGlobalNLS_.open(globalnls_name.c_str());
	    fGlobalNLS_.setf(std::ios::scientific);       
	    fGlobalNLS_ << std::setw(12) << std::left << "Iter.(1)";
	    fGlobalNLS_ << std::setw(12) << std::left << "Glob-NLS(2)";
	    fGlobalNLS_ << std::setw(12) << std::left << "Loc-NLS(3)";
	    fGlobalNLS_ << std::setw(16) << std::left << "Dummy(4)";
	    fGlobalNLS_ << std::setw(16) << std::left << "Reduction(5)";
	    fGlobalNLS_ << std::setw(12) << std::left << "Jacob(6)";        
	    fGlobalNLS_ << std::setw(12) << std::left << "LS-Iter(7)";
	    fGlobalNLS_ << std::setw(16) << std::left << "LS-Norm1(8)";
	    fGlobalNLS_ << std::endl;
        
	    fSequence_.open(sequence_name.c_str());
	    fSequence_.setf(std::ios::scientific);        
            fSequence_ << std::setw(12) << std::left << "Iter.(1)";
            fSequence_ << std::setw(12) << std::left << "Glob-Seq.(2)";
            fSequence_ << std::setw(12) << std::left << "Loc-Seq.(3)";
            fSequence_ << std::setw(16) << std::left << "NormInf(4)";
            fSequence_ << std::setw(16) << std::left << "Norm1-Mean(5)";
            fSequence_ << std::setw(16) << std::left << "Norm2-Mean(6)";
            fSequence_ << std::setw(16) << std::left << "Norm1(7)";
            fSequence_ << std::setw(16) << std::left << "Norm2(8)";        
            fSequence_ << std::setw(12) << std::left << "Direct(9)";
            fSequence_ << std::setw(12) << std::left << "Newton(10)";
            fSequence_ << std::setw(12) << std::left << "ODE-I(11)";
            fSequence_ << std::setw(12) << std::left << "ODE-II(12)";
            fSequence_ << std::setw(12) << std::left << "ODE-III(13)";
            fSequence_ << std::setw(12) << std::left << "ODE-IV(14)";     
            fSequence_ << std::setw(12) << std::left << "Jacobian(15)";        
            fSequence_ << std::setw(12) << std::left << "Newton(16)";        
            fSequence_ << std::setw(12) << std::left << "Failures(17)";        
            fSequence_ << std::endl;
        
	    fGlobalODE_.open(globalode_name.c_str());
	    fGlobalODE_.setf(std::ios::scientific);        
            fGlobalODE_ << std::setw(12) << std::left << "Iter.(1)";
            fGlobalODE_ << std::setw(12) << std::left << "Glob-ODE(2)";
            fGlobalODE_ << std::setw(12) << std::left << "Loc-ODE(3)";
            fGlobalODE_ << std::setw(16) << std::left << "TimeStep(4)";
            fGlobalODE_ << std::setw(16) << std::left << "Reduction(5)";
            fGlobalODE_ << std::setw(12) << std::left << "Jacob(6)";        
            fGlobalODE_ << std::setw(12) << std::left << "LS-Iter(7)";
            fGlobalODE_ << std::setw(16) << std::left << "LS-Norm1(8)";
            fGlobalODE_ << std::endl;
        
        
            cpu_ = new CPUTimer(data_.nameFolderOutput() + "/Additional/CPUTime.out");
            std::cout << "Done: All" << std::endl; //getchar();
	}

	InitializeCSTRVectors();
}

void OpenSMOKE_KPP_ReactorNetwork::AssemblingConvectionDiffusionMatrix()
{
	communicator_->InitializeArray(React_transport);

    	{
            for(int k = 1; k <= NR_P[procrank_]; k++)
            {
		int index = k + offset[procrank_];
                React_transport[index] = reactors_[k].M();
            }
        }

	communicator_->GatherArray(React_transport, NR_P, offset);

	if(procrank_ == 0)
        {
            for (int k=1; k <= NR_; k++)
            {
                C_(k,k) = React_transport[k];

                for(int i = 1; i <= ReactInSize[k]; i++)
                {
                    C_(k, React_in[k][i]) -= React_cIn[k][i];
                }

                for(int i = 1; i <= ReactNeighSize[k]; i++)
                {
                    C_(k, React_neighbours[k][i]) -= React_diffusion[k][i];
                }
            }
        }
        
        delete [] React_transport;
}

void OpenSMOKE_KPP_ReactorNetwork::AssemblingRightHandSides()
{
	ReactfInSize = mix_[0].NumberOfSpecies();
	communicator_->InitializeVector(React_fIn, ReactfInSize);

        {
            for(int k = 1; k <= NR_P[procrank_]; k++)
            {
		int index = k + offset[procrank_];
                for(int j = 1; j <= ReactfInSize; j++)
                    React_fIn[index][j] = reactors_[k].fIn()[j];
            }
        }

	communicator_->GatherVector(React_fIn, ReactfInSize);

}

void OpenSMOKE_KPP_ReactorNetwork::SummaryConvectionDiffusionMatrix()
{
	if(procrank_ == 0)
	{
	    int lower, upper;
	    BzzVectorInt numEquationsForEachVariable;
	    BzzVectorInt numVariablesForEachEquation;
	    C_.GetNumEquationsForEachVariable(&numEquationsForEachVariable);
	    C_.GetNumVariablesForEachEquation(&numVariablesForEachEquation);
	    C_.FindBands(&lower,&upper);

	    long long int non_zero_coeff = double(C_.Rows())*double(C_.Rows());

	    std::cout << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	    std::cout << "//           Convection/Diffusion Matrix Analysis              // " << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	    std::cout << "Number of equations:             " << C_.Rows() << std::endl;
	    std::cout << "Number of non zero coefficients: " << C_.CountElements() << "/" << non_zero_coeff << std::endl;
	    std::cout << "Filling:                         " << C_.CountElements()/double(C_.Rows())/double(C_.Rows())*100. << "%" << std::endl;
	    std::cout << "Mean number of non zero coefficients per equation: " << double(numVariablesForEachEquation.GetSumElements())/double(C_.Rows()) << std::endl;
	    std::cout << "Max number of non zero coefficients per equation:  " << numVariablesForEachEquation.Max() << std::endl;
	    std::cout << "Max number of equations per variables:             " << double(numEquationsForEachVariable.GetSumElements())/double(C_.Rows())<< std::endl;
	    std::cout << "Mean number of equations per variables:            " << numEquationsForEachVariable.Max() << std::endl;
	    std::cout << "Band width: upper(" << upper << ")  lower(" << lower << ")" << std::endl;
	    std::cout << std::endl;
	}
}

void OpenSMOKE_KPP_ReactorNetwork::CorrectingMassFlowRates()
{

    	BzzVector inputMassFlowRate(NR_);
    	BzzVector outputMassFlowRate(NR_);
    	BzzMatrixSparse splitting(NR_, NR_);

    	if(procrank_ > 0)
	{
	    ChangeDimensions(0, &inputMassFlowRate);
	    ChangeDimensions(0, &outputMassFlowRate);
	}

    	splitting.SetDiagonal(0, 1.);

	communicator_->InitializeArray(mOut);
	communicator_->InitializeArray(inputMassFlowRate_aux);
	communicator_->InitializeArray(outputMassFlowRate_aux);
	communicator_->InitializeVector(Split_Out, ReactOutSize);


	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = k + offset[procrank_];
            mOut[index] = reactors_[k].MassFlowOut();
            inputMassFlowRate_aux[index] = reactors_[k].fInTot();
	}

	communicator_->GatherArray(mOut, NR_P, offset);
	communicator_->GatherArray(inputMassFlowRate_aux, NR_P, offset);	

        if(procrank_ == 0)
        {
	    for(int i = 1; i <= NR_; i++)
	    {
	        inputMassFlowRate(i) = inputMassFlowRate_aux[i];
	    }

	    for(int k = 1; k <= NR_; k++)
	    {
	        for(int j = 1; j <= ReactOutSize[k]; j++)
	        {
                splitting(k, React_out[k][j]) = -React_cOut[k][j]/mOut[k];
	        }
	    }   

            Transpose(&splitting);
	
        }

        if (data_.sparseLinearSolverMassFlowRate() == KPP_SPARSESOLVER_BZZ)
        {
	    if(procrank_ == 0)
	    {
		std::cout << "Using BzzSolver..." << std::endl;
        	BzzFactorizedSparseGauss splittingFactorized;
        	splittingFactorized = splitting;
        	outputMassFlowRate = inputMassFlowRate;
        	Solve(&splittingFactorized, &outputMassFlowRate);
	    }
        }

        if (data_.sparseLinearSolverMassFlowRate() == KPP_SPARSESOLVER_PARDISO)
        {
	    if(procrank_ == 0)
	    {
	        std::cout << "Using Pardiso Solver..." << std::endl;
                OpenSMOKE_PARDISO_Unsymmetric pardisoSplitting(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
                pardisoSplitting.SetFillInReducingOrdering(data_.FillInReducingOrdering());
                pardisoSplitting.UnsetMessage();
                pardisoSplitting.SetSparsityPattern(splitting);
                pardisoSplitting.UpdateMatrix(splitting);
                pardisoSplitting.NumericalFactorization();
                pardisoSplitting.Solve(inputMassFlowRate, outputMassFlowRate);
                std::cout << outputMassFlowRate.Max() << " " << outputMassFlowRate.Min() << " " << Mean(outputMassFlowRate) << std::endl;
		pardisoSplitting.Delete();
	    }
        }

	else if (data_.sparseLinearSolverMassFlowRate() == KPP_SPARSESOLVER_LIS)
	{
	    if(procrank_ == 0) std::cout << "Using LIS Solver..." << std::endl;
	    OpenSMOKE_LIS_Unsymmetric lisSplitting(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
	    lisSplitting.SetCommunicator(communicator_);
	    lisSplitting.SetData(&data_);
	    lisSplitting.UnsetMessage();
	    lisSplitting.ResetAllCounters();
	    lisSplitting.SetSparsityPattern(splitting);
	    lisSplitting.CompleteMatrix();
	    lisSplitting.NumericalFactorization();

	    lisSplitting.UpdateVectors(inputMassFlowRate, outputMassFlowRate);
	    lisSplitting.Solve(inputMassFlowRate, outputMassFlowRate);
	    lisSplitting.Delete();	    
	    if(procrank_ == 0)
		std::cout << outputMassFlowRate.Max() << " " << outputMassFlowRate.Min() << " " << Mean(outputMassFlowRate) << std::endl;
	}

    	if (data_.sparseLinearSolverMassFlowRate() == KPP_SPARSESOLVER_MUMPS)
        {
	    ErrorMessage("Mumps sparse linear solver not implemented yet in the parallel form!");
        }
	

	if(procrank_ == 0)
	{
	    for(int k = 1; k <= NR_; k++)
            {
                outputMassFlowRate_aux[k] = outputMassFlowRate[k];
            }
	}

	communicator_->BroadcastArray(outputMassFlowRate_aux, NR_);

        if(procrank_ == 0)
        {
            // Update out flow rates from each reactor
            Transpose(&splitting);

            for(int k = 1; k <= NR_; k++)
            {
        	for(int j = 1; j <= ReactOutSize[k]; j++)
       		{
                    Split_Out[k][j] = splitting(k, React_out[k][j]);
                }
            }
        }

	communicator_->BroadcastVector(Split_Out, ReactOutSize);

	for(int k = 1; k <= NR_P[procrank_]; k++)
    	{
	    int index = k + offset[procrank_];
            for(int j = 1; j <= ReactOutSize[index]; j++)
            	reactors_[k].setcOut(j, -Split_Out[index][j]*outputMassFlowRate_aux[index]);
    	}

	TopologyParameters();

	// Update in flow rates to each reactor
        for(int k = 1; k <= NR_P[procrank_]; k++)
        {
	    int index = k + offset[procrank_];
            for(int j = 1; j <= ReactInSize[index]; j++)
            {
                int fromReactor = React_in[index][j];
                for(int i = 1; i <= ReactOutSize[fromReactor]; i++)
                {
                    int toReactor = React_out[fromReactor][i];
                    if (toReactor == index)
                    {
                        reactors_[k].setcIn(j, React_cOut[fromReactor][i]);
                        break;
                    }
                }
            }
        }

            // Update external output flows
            for (int k=1;k<=NR_P[procrank_];k++)
                reactors_[k].UpdateOutputFlow();

}

void OpenSMOKE_KPP_ReactorNetwork::SetTimeStep(const double deltat_)
{
	// Set Time Step
	deltat = deltat_;
	
	// Set A Matrix: A = I + C*dt
	{
		A_ = C_;

		double* ptrVal;
		int i, j;
		double val;

		A_.BeginScanning();
		if (data_.PredictorCorrector_DeferredConvection() == true)
		{
			if (data_.PredictorCorrector_MultiTimeSplitting() == false)
			{
				while(ptrVal = A_.Scanning(&i,&j,&val))
				{
					*ptrVal *= (deltat/reactors_[i].mass());
					if (i==j)	
						*ptrVal  = 1.;
				}
			}
			else
			{
				while(ptrVal = A_.Scanning(&i,&j,&val))
				{
					*ptrVal *= (timeSteps[i]/reactors_[i].mass());
					if (i==j)	
						*ptrVal  = 1.;
				}
			}
		}
		else
		{
			if (data_.PredictorCorrector_MultiTimeSplitting() == false)
			{
				while(ptrVal = A_.Scanning(&i,&j,&val))
				{
					*ptrVal *= (deltat/reactors_[i].mass()); 
					if (i==j)	
						*ptrVal += 1.;
				}
			}
			else
			{
				while(ptrVal = A_.Scanning(&i,&j,&val))
				{
					*ptrVal *= (timeSteps[i]/reactors_[i].mass());
					if (i==j)	
						*ptrVal += 1.;
				}
			}
		}
	}

	// Set Convection Matrix
	if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
	{
		pardisoMatrixConvection->ResetCounters();
		pardisoMatrixConvection->UpdateMatrix(A_);
		pardisoMatrixConvection->NumericalFactorization();
		pardisoMatrixConvection->UnsetMessage();
	}
/*	else if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
	{
		mumpsMatrixConvection->ResetCounters();
		mumpsMatrixConvection->UpdateMatrix(A_);
		mumpsMatrixConvection->NumericalFactorization();
		mumpsMatrixConvection->UnsetMessage();
	}*/
	else if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
	{
		lisMatrixConvection->ResetCounters();
		lisMatrixConvection->UpdateMatrix(A_);
		lisMatrixConvection->NumericalFactorization();
		lisMatrixConvection->UnsetMessage();
		lisMatrixConvection->SetMessage();
		lisMatrixConvection->EnableInitialSolution();
	}
}

void OpenSMOKE_KPP_ReactorNetwork::PredictorCorrector(const double deltat)
{	
/*	// Statistics
	statisticsConvective_->Reset();

	// Memory allocation (TO MOVE)
	BzzVector b(NR_);
	BzzVector solution(NR_);
	BzzVector omegaHalf(NumberOfEquations());
	BzzVector omegaOld(NumberOfEquations());
	BzzVector aux_omega(NR_);

	// Save Old(k) mass fractions in every computational cell
	{
		int count=1;
		for(int j=1;j<=mix_[0].NumberOfSpecies();j++)
			for(int k=1;k<=NR_;k++)
				omegaOld[count++] = reactors_[k].omega()[j];
	}
	
	// 1. PREDICTOR 
	cout << " * Predictor..." << std::endl;
	for(int j=1;j<=mix_[0].NumberOfSpecies();j++)
	{
		int position = (j-1)*NR_;
		if (data_.PredictorCorrector_DeferredConvection() == true)
		{
			if (data_.PredictorCorrector_MultiTimeSplitting() == false)
			{
				for(int k=1;k<=NR_;k++)
					b[k] = (1.-reactors_[k].M()*deltat/reactors_[k].mass())*omegaOld[position+k] + 
								deltat/reactors_[k].mass()*(RHS_[j][k]+reactors_[k].RV(j));
			}
			else
			{
				for(int k=1;k<=NR_;k++)
					b[k] = (1.-reactors_[k].M()*timeSteps[k]/reactors_[k].mass())*omegaOld[position+k] + 
								timeSteps[k]/reactors_[k].mass()*(RHS_[j][k]+reactors_[k].RV(j));
			}
		}
		else
		{
			if (data_.PredictorCorrector_MultiTimeSplitting() == false)
			{
				for(int k=1;k<=NR_;k++)
					b[k] =	omegaOld[position+k] + 
							deltat/reactors_[k].mass()*(RHS_[j][k]+reactors_[k].RV(j));
			}
			else
			{
				for(int k=1;k<=NR_;k++)
					b[k] =	omegaOld[position+k] + 
							timeSteps[k]/reactors_[k].mass()*(RHS_[j][k]+reactors_[k].RV(j));
			}
		}

		// Linear system solution
		if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
		{
			pardisoMatrixConvection->Solve(b, solution);
		}
		else if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
		{
			mumpsMatrixConvection->Solve(b, solution);
		}
		else if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
		{
			BzzVector aux(NR_);
			for(int i=1;i<=NR_;i++)
				aux[i] = omegaOld[(j-1)*NR_+i];
			lisMatrixConvection->SetInitialSolution(aux);
			lisMatrixConvection->Solve(b, solution);
		}

		// Omega Half
		omegaHalf.SetBzzVector(position+1, solution);
	}

	// Cleaning machine errors
	CleanVectorOnTheBottom(		0., -1.e-16,		omegaHalf);
	CleanVectorOnTheTop(		1.,	 1.+1.e-16,	omegaHalf);
	
	// Correction of time step
	double a = statisticsConvective_->ReductionOfTimeStep(omegaHalf, omegaOld);

	// Correction (time step)
	if (statisticsConvective_->iTimeCorrectedMin_ > 0 ||  statisticsConvective_->iTimeCorrectedMax_ > 0)
	{
		a *= data_.PredictorCorrector_TimeStepReductionFactor();
		for(int j=1;j<=NumberOfEquations();j++)
			omegaHalf[j] = a*omegaHalf[j] + (1.-a)*omegaOld[j];
	}

	// Additional Analysis
	statisticsConvective_->Analysis(deltat, a, omegaHalf);

	// Normalization
	for(int k=1;k<=NR_;k++)
	{
		double sum = 0.;
		for(int j=1;j<=mix_[0].NumberOfSpecies();j++)
			sum += omegaHalf[(j-1)*NR_+k];
		for(int j=1;j<=mix_[0].NumberOfSpecies();j++)
			omegaHalf[(j-1)*NR_+k]/=sum;
	}

	// Communication: Update transport
	{
		cout << " * Communication" << std::endl;
		// a. Updating of mass fractions for every computational cell (k+1/2)
		for(int j=1;j<=mix_[0].NumberOfSpecies();j++)
		{
			int position=(j-1)*NR_;
			for(int k=1;k<=NR_;k++)
				aux_omega[k] = omegaHalf[position+k];
			UpdateMassFractions(j, aux_omega);
		}
		
		// b. Distribution of data between cell (k+1/2)
		for(int k=1;k<=NR_;k++)
			reactors_[k].DistributeMassFlowRates(*this);
	
		// c. Go back to old (k) solution
		for(int j=1;j<=mix_[0].NumberOfSpecies();j++)
		{
			int position=(j-1)*NR_;
			for(int k=1;k<=NR_;k++)
				aux_omega[k] = omegaOld[position+k];
			UpdateMassFractions(j, aux_omega);
		}
	}

	// 2. CORRECTOR (Serial)
	if (data_.iOpenMP() == false)
	{
                BzzMatrix tmpMatrix(NumberOfSpecies(), NumberOfSpecies());
            
		if (data_.SequenceCSTR_UpdatingNetwork() == true)	// Continous
		{
			if (data_.PredictorCorrector_AlgebraicConstraint() == true)
			{
				cout << " * Corrector (continous, serial, dae)..." << std::endl;
				for(int k=1;k<=NR_;k++)
                                {
                                    odePool->Set(0, reactors_[k]);
                                    reactors_[k].SolveCSTR_CorrectorContinous_Smart(data_.SingleReactor_IntegrationTime(), *this, tmpMatrix, odePool->o(0));
                                }
			}
			else
			{
				cout << " * Corrector (continous, serial, ode)..." << std::endl;
				for(int k=1;k<=NR_;k++)
                                {
                                        odePool->Set(0, reactors_[k]);
                                        reactors_[k].SolveCSTR_CorrectorContinous(a*deltat, *this, odePool->o(0));
                                }
			}
		}
		else	// Discrete
		{
			if (data_.PredictorCorrector_AlgebraicConstraint() == true)
			{
				cout << "* Corrector (discrete, serial, dae)..." << std::endl;		
				for(int k=1;k<=NR_;k++)
                                {
                                        odePool->Set(0, reactors_[k]);
					reactors_[k].SolveCSTR_CorrectorDiscrete_Smart(data_.SingleReactor_IntegrationTime(), tmpMatrix, odePool->o(0));
                                }
			}
			else
			{
				cout << "* Corrector (discrete, serial, ode)..." << std::endl;	
				for(int k=1;k<=NR_;k++)
                                {
                                    odePool->Set(0, reactors_[k]);
				    reactors_[k].SolveCSTR_CorrectorDiscrete(a*deltat, odePool->o(0));
                                }
			}
		}
	}

	// 2. CORRECTOR (OpenMP)
	else
	{
		int tid;
                
                BzzMatrix* tmpMatrix = new BzzMatrix[data_.nThreads()];
                for (int kk=0;kk<data_.nThreads();kk++)
                    ChangeDimensions(NumberOfSpecies(), NumberOfSpecies(), &tmpMatrix[kk]);

		if (data_.SequenceCSTR_UpdatingNetwork() == true)       // Continous
		{
			if (data_.iSaveKineticConstants() == true)      // Saved kinetic constants
			{
				if (data_.PredictorCorrector_AlgebraicConstraint() == true)
				{                                        
                                        std::cout << "* Corrector (continous, openmp, dae, saved)..." << std::endl;		
            				#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num(); 
						for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
						{
							int k = indicesReactorsThreads_[tid][kk];
                                                        odePool->Set(tid, reactors_[k]);
                                                        reactors_[k].SolveCSTR_CorrectorContinous_Smart(data_.SingleReactor_IntegrationTime(), *this, tmpMatrix[tid], odePool->o(tid));
                                                }	
					}                            
				}
				else
				{
					cout << "* Corrector (continous, openmp, ode, saved)..." << std::endl;		
            				#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num(); 
						for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
						{
							int k = indicesReactorsThreads_[tid][kk];
                                                        odePool->Set(tid, reactors_[k]);
                                                        reactors_[k].SolveCSTR_CorrectorContinous(a*deltat, *this, odePool->o(tid));
                                                }	
					} 	
				}
			}
                        else      // Recalculated kinetic constants
			{
				if (data_.PredictorCorrector_AlgebraicConstraint() == true)
				{
					cout << "* Corrector (continous, openmp, dae, calculated)..." << std::endl;		
					#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num(); 
						for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
						{
							int k = indicesReactorsThreads_[tid][kk];
							reactors_[k].SetKinetics(&kinetics_[tid]);
                                                        odePool->Set(tid, reactors_[k]);
							reactors_[k].SolveCSTR_CorrectorContinous_Smart(data_.SingleReactor_IntegrationTime(), *this, tmpMatrix[tid], odePool->o(tid));
						}	
					}	
				}
				else
				{
					cout << "* Corrector (continous, openmp, ode, calculated)..." << std::endl;
					#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num(); 
						for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
						{
							int k = indicesReactorsThreads_[tid][kk];
							reactors_[k].SetKinetics(&kinetics_[tid]);
                                                        odePool->Set(tid, reactors_[k]);
							reactors_[k].SolveCSTR_CorrectorContinous(a*deltat, *this, odePool->o(tid));
						}	
					}	
				}
			}
		}
                else    // Discrete
		{
			if (data_.iSaveKineticConstants() == true)      // Saved kinetic constants
			{
				if (data_.PredictorCorrector_AlgebraicConstraint() == true)
				{
                                        std::cout << "* Corrector (discrete, openmp, dae, saved)..." << std::endl;		
            				#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num(); 
						for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
						{
							int k = indicesReactorsThreads_[tid][kk];
                                                        odePool->Set(tid, reactors_[k]);
                                                        reactors_[k].SolveCSTR_CorrectorDiscrete_Smart(data_.SingleReactor_IntegrationTime(), tmpMatrix[tid], odePool->o(tid));
                                                }	
					} 
                                
                                }
				else    if (data_.iSaveKineticConstants() == true && data_.nThreads() == 1)
				{
					cout << "* Corrector (discrete, openmp, ode, saved)..." << std::endl;	
					#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num(); 
						for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
						{
							int k = indicesReactorsThreads_[tid][kk];
                                                        odePool->Set(tid, reactors_[k]);
                                                        reactors_[k].SolveCSTR_CorrectorDiscrete(a*deltat, odePool->o(tid));
                                                }	
					} 
				}
			}
                        else    // Recalculated kinetic constants
			{
				if (data_.PredictorCorrector_AlgebraicConstraint() == true)
				{
					cout << "* Corrector (discrete, openmp, dae, calculated)..." << std::endl;		
					#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num(); 
						for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
						{
							int k = indicesReactorsThreads_[tid][kk];
							reactors_[k].SetKinetics(&kinetics_[tid]);
                                                        odePool->Set(tid, reactors_[k]);
							reactors_[k].SolveCSTR_CorrectorDiscrete_Smart(data_.SingleReactor_IntegrationTime(), tmpMatrix[tid], odePool->o(tid));
						}	
					}
				}
				else
				{
					cout << "* Corrector (discrete, openmp, ode, calculated)..." << std::endl;	
					#pragma omp parallel private(tid)
					{
						tid = omp_get_thread_num(); 
						for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
						{
							int k = indicesReactorsThreads_[tid][kk];
							reactors_[k].SetKinetics(&kinetics_[tid]);
                                                        odePool->Set(tid, reactors_[k]);
							reactors_[k].SolveCSTR_CorrectorDiscrete(a*deltat, odePool->o(tid));
						}
					}
				}		// No algebraic constraints
			}		// Calculated kinetic constants
		}		// Discrete updating
	}		// Corrector (OpenMP)*/
}

void OpenSMOKE_KPP_ReactorNetwork::TimeStepPolicy(double &deltat)
{
	// If no corrections usually occur
	if (statisticsConvective_->iterationsWithoutCorrections_ >= 5)
	{
		deltat *= data_.PredictorCorrector_TimeStepIncrementFactor();
		deltat  = std::min(deltat, data_.PredictorCorrector_MaxTimeStep());
		statisticsConvective_->iterationsWithoutCorrections_ = 0;
		statisticsConvective_->iterationsWithCorrections_    = 0;

		// New factorization
		SetTimeStep(deltat);
	}
	
	if (statisticsConvective_->iterationsWithCorrections_ >= 5)
	{
		deltat *= data_.PredictorCorrector_TimeStepReductionFactor();
		deltat  = std::min(deltat, data_.PredictorCorrector_MaxTimeStep());
		statisticsConvective_->iterationsWithoutCorrections_ = 0;
		statisticsConvective_->iterationsWithCorrections_    = 0;

		// New factorization
		SetTimeStep(deltat);
	}
}

void OpenSMOKE_KPP_ReactorNetwork::SolveSequenceCSTR()
{
	if(nLocalCSTRSequence_ == 0)
	{
	    SendAndBroadcastOmega();
	}

	if(nGlobalCSTRSequence_ == 0 && procrank_ == 0)
	    cpu_->SetStartTimeAllGlobal();

	// Counters
	nLocalCSTRSequence_++;
	nGlobalCSTRSequence_++;

	// Statistics
	if (data_.SequenceCSTR_VerboseStatistics() == true)
		statistics_->Reset();
	
	// Discrete Updating
	{
	    if (data_.SequenceCSTR_UpdatingNetwork() == false)      // Discrete
	    {
		for(int k=1;k<=NR_P[procrank_];k++)
		    reactors_[k].DistributeMassFlowRates(*this);;
//		if (data_.iOpenMP() == true)
		{
                    int tid;
                    BzzMatrix* tmpMatrix = new BzzMatrix[data_.nThreads()];
                    for (int kk=0;kk<data_.nThreads();kk++)
                        ChangeDimensions(NumberOfSpecies(), NumberOfSpecies(), &tmpMatrix[kk]);                                    
		    if (data_.iSaveKineticConstants() == true)
		    {
			if(procrank_ == 0)
			    std::cout << "* Sequence (openmp, discrete, saved)..." << std::endl;

                        #pragma omp parallel private(tid)
			{
			    tid = omp_get_thread_num(); 
			    for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
			    {
				int k = indicesReactorsThreads_[tid][kk];
                                odePool->Set(tid, reactors_[k]);
				reactors_[k].SolveCSTR_CorrectorDiscrete_Smart(data_.SingleReactor_IntegrationTime(), tmpMatrix[tid], odePool->o(tid));
			    }
			}
		    }

		    else
		    {
			if(procrank_ == 0)
			    std::cout << "* Sequence (openmp, discrete, calculated)..." << std::endl;
				
			#pragma omp parallel private(tid)
			{
			    tid = omp_get_thread_num(); 
			    for(int kk=1; kk<=indicesReactorsThreads_[tid].Size(); kk++)
			    {
				int k = indicesReactorsThreads_[tid][kk];
				reactors_[k].SetKinetics(&kinetics_[tid]);
                                odePool->Set(tid, reactors_[k]);
				reactors_[k].SolveCSTR_CorrectorDiscrete_Smart(data_.SingleReactor_IntegrationTime(), tmpMatrix[tid], odePool->o(tid));
			    }
			}
		    }
		    delete [] tmpMatrix;
		}
                
	    }	
            else if (data_.SequenceCSTR_UpdatingNetwork() == true)  // Continous 
	        ErrorMessage("Continuous update not possible in the parallel program!");
	}

}

int OpenSMOKE_KPP_ReactorNetwork::AnalysisSequenceCSTR()
{
	iDirectConvergence = 0;
	iNewtonConvergence = 0;
	iOdeFirstConvergence = 0;
	iOdeSecondConvergence = 0;
	iOdeThirdConvergence = 0;
	iOdeFourthConvergence = 0;
	nJacobianEvaluations = 0;
	nNewtonIterations = 0;
	nFailures = 0;
	F1MaxNewton = 0;
	F1MaxODE = 0;
	F1MaxUnconverged = 0;
	normInf = 0;
	F1Mean = 0;

	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	    {			
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_DIRECT)
			iDirectConvergence++;
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_NEWTON)
		{
			iNewtonConvergence++;
			nJacobianEvaluations+=reactors_[k].status.nJacobianEvaluations;
			nNewtonIterations+=reactors_[k].status.nNewtonIterations;
			F1MaxNewton = std::max(F1MaxNewton, reactors_[k].status.norm1_over_nc);
		}
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_ODE_FIRST)
		{
			iOdeFirstConvergence++;
			F1MaxODE = std::max(F1MaxODE, reactors_[k].status.norm1_over_nc);
		}
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_ODE_SECOND)
		{
			iOdeSecondConvergence++;
			F1MaxODE = std::max(F1MaxODE, reactors_[k].status.norm1_over_nc);
		}
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_ODE_THIRD)
		{
			iOdeThirdConvergence++;
			F1MaxODE = std::max(F1MaxODE, reactors_[k].status.norm1_over_nc);
		}
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_ODE_FOURTH)
		{
			iOdeFourthConvergence++;
			F1MaxODE = std::max(F1MaxODE, reactors_[k].status.norm1_over_nc);
		}

		normInf  = std::max(normInf, reactors_[k].status.normInf);
		F1Mean  += reactors_[k].status.norm1_over_nc;
		nFailures += reactors_[k].status.failure;
	    }
	}

	SendCSTRInfoToMaster();		

	if(procrank_ == 0)
	{
	    std::cout.precision(2);
	    std::cout << "Reactor distribution" << std::endl;
	    std::cout << " * Already OK: " << std::setw(9) << std::fixed << std::right << double(iDirectConvergence)/double(NR_)*100.		<< "%" << std::setw(8) << std::right << iDirectConvergence << "/" << NR_ << std::endl;
	    std::cout << " * Newton:     " << std::setw(9) << std::fixed << std::right << double(iNewtonConvergence)/double(NR_)*100.		<< "%" << std::setw(8) << std::right << iNewtonConvergence << "/" << NR_ << std::endl;
	    std::cout << " * ODE(I):     " << std::setw(9) << std::fixed << std::right << double(iOdeFirstConvergence)/double(NR_)*100.	<< "%" << std::setw(8) << std::right << iOdeFirstConvergence << "/" << NR_ << std::endl;
	    std::cout << " * ODE(II):    " << std::setw(9) << std::fixed << std::right << double(iOdeSecondConvergence)/double(NR_)*100.		<< "%" << std::setw(8) << std::right << iOdeSecondConvergence << "/" << NR_ << std::endl;
	    std::cout << " * ODE(III):   " << std::setw(9) << std::fixed << std::right << double(iOdeThirdConvergence)/double(NR_)*100.	<< "%" << std::setw(8) << std::right << iOdeThirdConvergence << "/" << NR_ << std::endl;
	    std::cout << " * ODE(IV):    " << std::setw(9) << std::fixed << std::right << double(iOdeFourthConvergence)/double(NR_)*100.	<< "%" << std::setw(8) << std::right << iOdeFourthConvergence << "/" << NR_ << std::endl;
	    std::cout << std::endl;

	    std::cout << "Additional data" << std::endl;
	    std::cout << " * Jacobians per reactor:   " << std::setw(9) << std::fixed << std::right << double(nJacobianEvaluations)/double(std::max(1,iNewtonConvergence)) << std::endl;
	    std::cout << " * Newtons per reactor:     " << std::setw(9) << std::fixed << std::right << double(nNewtonIterations)/double(std::max(1,iNewtonConvergence)) << std::endl;
	    std::cout << " * Mass fractions failures: " << std::setw(9) << std::fixed << std::right << nFailures << std::endl;
	    std::cout << std::endl;

	    std::cout << "Residuals" << std::endl;
	    std::setprecision(8);
	    std::cout << " * Max residual Newton:    " << std::setw(9) << std::scientific << F1MaxNewton		<< std::endl;
	    std::cout << " * Max residual ODE:       " << std::setw(9) << std::scientific << F1MaxODE			<< std::endl;
	    std::cout << " * Max residual Failed:    " << std::setw(9) << std::scientific << F1MaxUnconverged	<< std::endl;
	    std::cout << " * Mean residual:          " << std::setw(9) << std::scientific << F1Mean			<< std::endl;
	    std::cout << " * Norm Inf:               " << std::setw(9) << std::scientific << normInf			<< std::endl;
	    std::cout << std::endl;

	    fSequence_ << std::setw(12)  << std::left << iteration_;
	    fSequence_ << std::setw(12)  << std::left << nGlobalCSTRSequence_;
	    fSequence_ << std::setw(12)  << std::left << nLocalCSTRSequence_;
	    fSequence_ << std::setw(16) << std::left << std::scientific << residuals_->normInf();
	    fSequence_ << std::setw(16) << std::left << std::scientific << residuals_->norm1()/NR_;
	    fSequence_ << std::setw(16) << std::left << std::scientific << residuals_->norm2()/NR_;        
	    fSequence_ << std::setw(16) << std::left << std::scientific << residuals_->norm1();
	    fSequence_ << std::setw(16) << std::left << std::scientific << residuals_->norm2();
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << double(iDirectConvergence)/double(NR_)*100.;
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << double(iNewtonConvergence)/double(NR_)*100.;
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << double(iOdeFirstConvergence)/double(NR_)*100.;
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << double(iOdeSecondConvergence)/double(NR_)*100.;
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << double(iOdeThirdConvergence)/double(NR_)*100.;
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << double(iOdeFourthConvergence)/double(NR_)*100.;
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << double(nJacobianEvaluations)/double(std::max(1,iNewtonConvergence));
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << double(nNewtonIterations)/double(std::max(1,iNewtonConvergence));
	    fSequence_ << std::setw(12) << std::left << std::fixed << std::setprecision(2) << nFailures;

            fSequence_ << std::endl;
	}

	MPI::COMM_WORLD.Bcast(&F1Mean, 1, MPI::DOUBLE, MASTER);
	MPI::COMM_WORLD.Bcast(&F1MaxODE, 1, MPI::DOUBLE, MASTER);
	MPI::COMM_WORLD.Bcast(&F1MaxNewton, 1, MPI::DOUBLE, MASTER);
	MPI::COMM_WORLD.Bcast(&nFailures, 1, MPI::DOUBLE, MASTER);

	if ( F1Mean > 1.e-4 && nLocalCSTRSequence_> 150)
	{
		if(procrank_ == 0)
		{
		    fLog_ << "Sequence: Mean residual too large (" << F1Mean << ")" << std::endl;
		    std::cout  << "Convergence Failure in CSTR Sequence: Mean residual too large" << std::endl;
		}
		CSTRConvergengeAlert = true;
		return -1;
	}

	if ( ((F1MaxODE>1.e-1) || (F1MaxNewton>1.e-1)) && nLocalCSTRSequence_> 150 && CSTRConvergengeAlert == false)
	{
		if(procrank_ == 0)
		{
		    fLog_ << "Sequence: Max residual too large (" << F1MaxODE << ", " << F1MaxNewton << ")" << std::endl;
		    fLog_ << "Iteration no. " << nLocalCSTRSequence_ << std::endl;
		    std::cout  << "Convergence Failure in CSTR Sequence: Max residual too large" << std::endl;
		}

		CSTRConvergengeAlert = true;
		return -1;
	}

	if ( double(nFailures)/double(NR_)*100 > 0.5  && nLocalCSTRSequence_>30)
	{
		if(procrank_ == 0)
		{
		    fLog_ << "Sequence: Too many sum of mass fractions failures (" << double(nFailures)/double(NR_)*100 << "%)" << std::endl;
		    std::cout << "Convergence Failure in CSTR Sequence: Too many sum of mass fractions failures..." << std::endl;
		}
		return -1;
	}

	double ratioNormInf, ratioNorm1;
	if(procrank_ == 0)
	{
	    ratioNormInf = Residuals().RatioNormInf();
	    ratioNorm1 = Residuals().RatioNorm1();
	}

	MPI::COMM_WORLD.Bcast(&ratioNormInf, 1, MPI::DOUBLE, MASTER);
	MPI::COMM_WORLD.Bcast(&ratioNorm1, 1, MPI::DOUBLE, MASTER);

	// Checking convergence index
	if ( (ratioNormInf > NormInfConvCriterium) && (ratioNorm1 > Norm1ConvCriterium) && nLocalCSTRSequence_> ReactIterCriterium && ratioNormInf < 1. && ratioNorm1 < 1. &&  nLocalCSTRSequence_ > 15)
	{
		if(procrank_ == 0)
		    fLog_ << "Sequence: Convergence indices satisfied..." << std::endl; 

		return 0;
	}
	// Checking convergence index
	if ( double(iOdeSecondConvergence+iOdeFirstConvergence)/double(NR_)*100.< 0.1 && ratioNormInf >0.96 && Residuals().norm1() > 1e-15)
	{
		if(procrank_ == 0)
		    fLog_ << "Sequence: Small enough number of ODE integrations..." << std::endl; 

		return 0;
	}

	return 1;
}
void OpenSMOKE_KPP_ReactorNetwork::InitializeGlobal(const KPP_SparseLinearSolver kind, const double absoluteTolerance, const double relativeTolerance)
{
	MPI::Status status;
	if(procrank_ == 0)
	{
	    std::cout << "// *********************************************************** // " << std::endl;
	    std::cout << "//                       Global method                         // " << std::endl;
	    std::cout << "// *********************************************************** // " << std::endl;
	}

	if(procrank_ == 0)	
	    ChangeDimensions (NumberOfEquations(),  &globalOmega_);

	// Number of non-zero elements
	int numberOfNonZeroElementsLoc=0;
	long long int numberOfNonZeroElements=0;

	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	        numberOfNonZeroElementsLoc += reactors_[k].nGlobalSparsityPattern();
	}

	for(int p = 0; p <= numworkers_; p++)
	{
	    MPI::COMM_WORLD.Barrier();
	    mtype = FROM_WORKER;

	    if(p == procrank_)
	    {
	        MPI::COMM_WORLD.Send(&numberOfNonZeroElementsLoc, 1, MPI::INT, MASTER, mtype);
	    }

	    if(procrank_ == 0)
	    {
		source = p;
		MPI::COMM_WORLD.Recv(&numberOfNonZeroElementsLoc, 1, MPI::INT, source, mtype, status);
		numberOfNonZeroElements += numberOfNonZeroElementsLoc;
	    }
	}

	MPI::COMM_WORLD.Bcast(&numberOfNonZeroElements, 1, MPI::LONG_LONG, MASTER);


	// Set PARDISO Matrix
	if (kind == KPP_SPARSESOLVER_PARDISO && iAlreadyGlobalPardiso_ == false)
	{
	    ErrorMessage("Pardiso not available in the parallel version");
	}

	// Set MUMPS Matrix
	else if (kind == KPP_SPARSESOLVER_MUMPS && iAlreadyGlobalMumps_ == false)
	{
	    ErrorMessage("MUMPS not available in the parallel version");
	}

	// Set LIS Matrix
	else if (kind == KPP_SPARSESOLVER_LIS && iAlreadyGlobalLis_ == false)
	{		
	    GetSingleReactorSparsityPattern();

	    lisMatrixGlobal->SetLocalRows(*this);
	    lisMatrixGlobal->OpenMatrix(NR_*mix_[0].NumberOfSpecies(), numberOfNonZeroElements);

	    lisMatrixGlobal->ResetAllCounters();

	    for(int k=1;k<=NR_;k++)
	    {
		lisMatrixGlobal->SetLocalSparsityPattern(k, NumberOfSpecies(), mConvDiff_, iConvDiff_, SparsityPattern, *this);
	    }


	    {
    	        for(int k = 1; k <= NR_P[procrank_]; k++)
		{
		    int index = k + offset[procrank_];
	            globalIndicesSparsityValuesSize[index] = reactors_[k].globalIndicesSparsityValues().Size();
		}
	    }

	    communicator_->GatherArray(globalIndicesSparsityValuesSize, NR_P, offset);

	    lisMatrixGlobal->MemoryAllocation(*this);


	    for(int k=1;k<=NR_P[procrank_];k++)
	    {
		int index = k + offset[procrank_];
		lisMatrixGlobal->SetSparsityPattern(index, NumberOfSpecies(), mConvDiff_, iConvDiff_, SparsityPattern, *this);
	    }

	    Petsc_values = new double[lisMatrixGlobal->countLocal()[procrank_]];

	    InitializePlaceArray();
/*
	    InitializePetscMatrix();

	    InitializePetscRHS();*/

	    lisMatrixGlobal->CloseMatrix();

	    iAlreadyGlobalLis_ = true;
	}

	else if (kind == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL && iAlreadyGlobalGaussSiedel_ == false)
	{
	    ErrorMessage("OpenSMOKE_GaussSeidel not available in the parallel version");
	}
}

int OpenSMOKE_KPP_ReactorNetwork::GlobalNLS()
{
        // CPU time
	if(procrank_ == 0)
	{
            cpu_->SetNLSGlobal();
            cpu_->SetStartTimeNLSGlobal();
	}
    
	indexGlobalNLS_++;
        localCounterNLS_=0;

	// First Guess solution
	if(procrank_ == 0)
	{
	    for(int k=1;k<=NR_;k++)
	    {
	        BzzVector omegaLoc(mix_[0].NumberOfSpecies());
	        for(int i = 1; i <= mix_[0].NumberOfSpecies(); i++)
		{
		    omegaLoc[i] = OmegaGlob_[k][i];
		}

 	        globalOmega_.SetBzzVector( (k-1)*mix_[0].NumberOfSpecies()+1, omegaLoc );
	    }
	}

	// Solve Non Linear System
	OpenSMOKE_KPP_NewtonMethod_Manager newton(globalOmega_, externalResiduals, externalSolutionNLS, externalResidualsAnalysis, *this);

	
	// Non linear system option
	newton.SetMethod(data_.GlobalNLS_Method());
	newton.SetRelativeTolerance(data_.GlobalNLS_RelativeTolerance());
	newton.SetAbsoluteTolerance(data_.GlobalNLS_AbsoluteTolerance());
	newton.SetMaximumNumberOfIterations(data_.GlobalNLS_MaxIterations());
	newton.SetMaximumNumberOfArmijioIterations(data_.GlobalNLS_MaxArmijioIterations());

	if(procrank_ == 0)
	{
	    fLog_ << "Index:             " << indexGlobalNLS_ << std::endl;
	    fLog_ << "Rel/Abs Tol:       " << data_.GlobalNLS_RelativeTolerance() << " " << data_.GlobalNLS_AbsoluteTolerance() << std::endl;
	    fLog_ << "Start Iteration:   " << iteration_ << std::endl;
	}
        
	int info = newton.Solve(globalOmega_);

	fLog_ << "End Iteration:     " << iteration_ << std::endl;
	
	if (info == 3 && procrank_ == 0)          fLog_ << "Status:            super convergence..." << std::endl;
	if (info == 2 && procrank_ == 0)          fLog_ << "Status:            slow convergence..." << std::endl;
	else if (info == 1 && procrank_ == 0)     fLog_ << "Status:            std::maximum number of iterations.." << std::endl;
        else if (info == 0 && procrank_ == 0)     fLog_ << "Status:            convergence indices satisfied..." << std::endl;
	else if (info == -1 && procrank_ == 0)    fLog_ << "Status:            linear system failure..." << std::endl;
	else if (info == -2 && procrank_ == 0)    fLog_ << "Status:            Armijio failure..." << std::endl;

	newtonconvergence_ = newton.convergence();
	
	return info;
}

int OpenSMOKE_KPP_ReactorNetwork::GlobalODE()
{
        // CPU time
	if(procrank_ == 0)
	{
            cpu_->SetODEGlobal();
            cpu_->SetStartTimeODEGlobal();
	}
        
	indexGlobalODE_++;
        localCounterODE_=0;

	// First Guess solution
	if(procrank_ == 0)
	{
	    for(int k=1;k<=NR_;k++)
	    {
	        BzzVector omegaLoc(mix_[0].NumberOfSpecies());
	        for(int i = 1; i <= mix_[0].NumberOfSpecies(); i++)
		{
		    omegaLoc[i] = OmegaGlob_[k][i];
		}

 	        globalOmega_.SetBzzVector( (k-1)*mix_[0].NumberOfSpecies()+1, omegaLoc );
	    }
	}

	int maxSubIterations = 1;
	double increasingFactor = 2.;
	double guessTimeStep;

	if(indexGlobalODE_ == 1)
	    guessTimeStep = data_.GlobalODE_InitialTimeStep();
	else if(indexGlobalODE_ > 1)
	    guessTimeStep = lastODEStep_;
	
	double guessIncreasingFactor = data_.GlobalODE_TimeStepIncrementFactor();

	// Solve ODE System
	for(int i=1;i<=maxSubIterations;i++)
	{	
		OpenSMOKE_KPP_ODE_Manager ode(globalOmega_, externalResiduals, externalSolutionODE, externalResidualsAnalysis, *this);
		
		ode.SetRelativeTolerance(data_.GlobalODE_RelativeTolerance());
		ode.SetAbsoluteTolerance(data_.GlobalODE_AbsoluteTolerance());
		ode.SetMaximumNumberOfIterations(data_.GlobalODE_MaxIterations());
		ode.SetStartingDeltaTime(guessTimeStep);
		ode.SetMaximumDeltaTime(guessIncreasingFactor);
		ode.SetIncreasingFactorDeltaTime(data_.GlobalODE_TimeStepIncrementFactor());
		ode.SetUpdatingFrequencyTimeStep(data_.GlobalODE_UpdatingFrequencyTimeStep());

		if(procrank_ == 0)
		{
		    fLog_ << "Index:                " << indexGlobalODE_ << "-" << i << std::endl;
		    fLog_ << "Rel/Abs Tol:          " << data_.GlobalODE_RelativeTolerance() << " " << data_.GlobalODE_AbsoluteTolerance() << std::endl;
		    fLog_ << "Initial time step:    " << guessTimeStep << std::endl;
		    fLog_ << "Increasing factor:    " << guessIncreasingFactor << std::endl;
		    fLog_ << "Start iteration:      " << iteration_ << std::endl;
		}
                
		int info = ode.Solve(globalOmega_);

		lastODEStep_ = ode.LastRequestedDeltaTime();

		if(procrank_ == 0)
		{
		    fLog_ << "End iteration:        " << iteration_ << std::endl;
		    fLog_ << "Last time step (req): " << ode.LastRequestedDeltaTime() << std::endl;

		    if (info == 4)     fLog_ << "Status:               super Convergence..." << std::endl;
                    else if (info == 3)     fLog_ << "Status:               very slow convergence..." << std::endl;
                    else if (info == 2)     fLog_ << "Status:               slow convergence..." << std::endl;
                    else if (info == 1)     fLog_ << "Status:               std::maximum number of iterations..." << std::endl;
		    else if (info == 0)     fLog_ << "Status:               convergence satisfied..." << std::endl;
		    else if (info == -1)    fLog_ << "Status:               linear system failure..." << std::endl;

		    fLog_ << std::endl;
		}

		if (info == 2)
		{
			guessTimeStep *= increasingFactor ;
			guessIncreasingFactor *= 1.25;
			continue;
		}

		if (info != 2)
		{
		    return info;
		}
	}

	return 2;
}

void OpenSMOKE_KPP_ReactorNetwork::ApplyStatistics(BzzOdeStiffObject &o)
{
	statistics_->Analysis(o);
}

void OpenSMOKE_KPP_ReactorNetwork::SolveWithoutReactions()
{
	BzzVector x(NR_);
	BzzMatrix omega(NR_, mix_[0].NumberOfSpecies());
	
	OpenSMOKE_PARDISO_Unsymmetric pardisoSystem(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
	pardisoSystem.SetFillInReducingOrdering(data_.FillInReducingOrdering());
	pardisoSystem.SetSparsityPattern(C_);
	pardisoSystem.UpdateMatrix(C_);
	pardisoSystem.NumericalFactorization();
	pardisoSystem.UnsetMessage();

	for(int j=1;j<=mix_[0].NumberOfSpecies();j++)
	{
		pardisoSystem.Solve(RHS_[j],x);
		omega.SetColumn(j, x);
	}

	UpdateMassFractions(omega);
	ResidualsAnalysis();

	// Write Maps on File
	WriteMassFractionMap();
}

void OpenSMOKE_KPP_ReactorNetwork::WriteMassFractionMap()
{
	if(procrank_ == 0)
	{
	    BzzVector* omega_mass = new BzzVector[NR_ + 1];
	    BzzVector* omega_mole = new BzzVector[NR_ + 1];
	    BzzVector MWmix(NR_);

	    for(int i = 1; i <= NR_; i++)
	    {
	        ChangeDimensions(mix_[0].NumberOfSpecies(), &omega_mass[i]);
	        ChangeDimensions(mix_[0].NumberOfSpecies(), &omega_mole[i]);
		for(int j = 1; j <= mix_[0].NumberOfSpecies(); j++)
		{
		    omega_mass[i][j] = OmegaGlob_[i][j];
		}

		MWmix[i] = mix()[0].GetMWFromMassFractions(omega_mass[i]);

		mix()[0].GetMoleFractionsFromMassFractionsAndMW(omega_mole[i], omega_mass[i], MWmix[i]);
	    }

	    for (int j=1;j<=data_.mapsSpeciesNames().size();j++)
	    {
		int index = mix_[0].recognize_species(data_.mapsSpeciesNames()[j-1]);

                std::ofstream fMap;
		openOutputFileAndControl(fMap, data_.nameFolderOutput() + "/Maps/Mass/Plot.plt." + mix_[0].names[index]);
		fMap.setf(std::ios::scientific);
		
		for (int k=1; k<=NR_; k++)
			fMap << OmegaGlob_[k][index] << std::endl;

		fMap.close();

		openOutputFileAndControl(fMap, data_.nameFolderOutput() + "/Maps/Mole/Plot.plt." + mix_[0].names[index]);
		fMap.setf(std::ios::scientific);
		
		for (int k=1; k<=NR_; k++)
			fMap << omega_mole[k][index] << std::endl;

		fMap.close();
	    }

   	    if (data_.iTraditionalLink()==true)
    	    {
        	BzzVectorInt fromOriginalToCluster;
        	BzzVector volumes;
        	BzzLoad fClusteringTopology(data_.nameFolderInput() + "/ClusteringTopology.out");
        	fClusteringTopology >> fromOriginalToCluster;
        	fClusteringTopology >> volumes;
        	fClusteringTopology.End();

        	for (int j=1;j<=data_.mapsSpeciesNames().size();j++)
        	{
            	    int index = mix_[0].recognize_species(data_.mapsSpeciesNames()[j-1]);

                    std::ofstream fMap;
            	    openOutputFileAndControl(fMap, data_.nameFolderOutput() + "/Maps/Mass/Plot.plt." + mix_[0].names[index] + ".original");
            	    fMap.setf(std::ios::scientific);

            	    for(int k=1;k<=fromOriginalToCluster.Size();k++)
            	    {
                	int i = fromOriginalToCluster[k];
                	fMap << OmegaGlob_[i][index] << std::endl;
            	    }

            	    fMap.close();

            	    openOutputFileAndControl(fMap, data_.nameFolderOutput() + "/Maps/Mole/Plot.plt." + mix_[0].names[index] + ".original");
            	    fMap.setf(std::ios::scientific);

            	    for(int k=1;k<=fromOriginalToCluster.Size();k++)
            	    {
                	int i = fromOriginalToCluster[k];
                	fMap << omega_mole[i][index] << std::endl;
            	    }

            	    fMap.close();
        	}
    	    }

	    delete [] omega_mass;
	    delete [] omega_mole;

	}
        
}

void OpenSMOKE_KPP_ReactorNetwork::ExternalFeeds(double &massFlow, BzzVector &omegaFeeds, BzzVector &omegaElementalFeeds)
{
	communicator_->InitializeVector(React_fIn, mix_[0].NumberOfSpecies());

	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
   	    {
		int index = k + offset[procrank_];
  	        for(int j = 1; j <= mix_[0].NumberOfSpecies(); j++)
   	        {
   		    React_fIn[index][j] = reactors_[k].fIn()[j];
   	        }
   	    }
	}

	communicator_->GatherVector(React_fIn, mix_[0].NumberOfSpecies());

	if(procrank_ == 0)
	{
	    // Species composition
	    omegaFeeds = 0.;
	    massFlow   = 0.;
	    for (int k=1; k <= jExternalFeedReactors.Size(); k++)
	    {
	        massFlow += React_fInTot[jExternalFeedReactors[k]];
		for (int i=1;i<=omegaFeeds.Size();i++)
                    omegaFeeds[i] += React_fIn[jExternalFeedReactors[k]][i];
	    }

	    double sum = omegaFeeds.GetSumElements();
	    omegaFeeds /= sum;

	    // Elemental composition
	    ChangeDimensions(mix_[0].NumberOfElements(), &omegaElementalFeeds);
	    mix_[0].GetElementalMassFractionsFromSpeciesMassFractions(omegaElementalFeeds, omegaFeeds);
	}

	communicator_->DeleteVector(React_fIn, NR_);
}

void OpenSMOKE_KPP_ReactorNetwork::ExternalOutput(double &massFlow, BzzVector &omegaOutput, BzzVector &omegaElementalOutput)
{
	if(procrank_ == 0)
	{
	    // Species composition
	    omegaOutput = 0.;
	    massFlow    = 0.;
	    for (int k=1;k<=jExternalOutputReactors.Size();k++)
	    {
	        double sum = React_fOutTot[jExternalOutputReactors[k]];
		massFlow += sum;
		for (int i=1;i<=omegaOutput.Size();i++)
                    omegaOutput[i] += sum*OmegaGlob_[jExternalOutputReactors[k]][i];
	    }

	    double sum = omegaOutput.GetSumElements();
	    omegaOutput /= sum;

	    // Elemental composition
	    ChangeDimensions(mix_[0].NumberOfElements(), &omegaElementalOutput);	
	    mix_[0].GetElementalMassFractionsFromSpeciesMassFractions(omegaElementalOutput, omegaOutput);
	}
}

void OpenSMOKE_KPP_ReactorNetwork::MinMaxMassFractions(BzzVector &omegaMin, BzzVector &omegaMax)
{
	if(procrank_ == 0)
	{
	    // Species composition
	    omegaMin =  1e16;
	    omegaMax = -1e-16;
	    for (int k=1;k<=NR_;k++)
	    {
        	for (int j=1;j<=NumberOfSpecies();j++)
                {
                    if (OmegaGlob_[k][j] < omegaMin[j]) omegaMin[j] = OmegaGlob_[k][j];  
                    if (OmegaGlob_[k][j] > omegaMax[j]) omegaMax[j] = OmegaGlob_[k][j];
                }
	    }
	}
}

void OpenSMOKE_KPP_ReactorNetwork::UpdateMassFractions(const int jSpecies, BzzVector &omega)
{
	for (int k=1;k<=NR_;k++)
		reactors_[k].setOmega(jSpecies, omega[k]);
}

void OpenSMOKE_KPP_ReactorNetwork::UpdateMassFractions(BzzVector &omega)
{
	ErrorMessage("UpdateMassFractions(BzzVector &omega) not yet implemented");
}

void OpenSMOKE_KPP_ReactorNetwork::UpdateMassFractions(BzzMatrix &omega)
{
//	for (int k=1;k<=NR_;k++)
//		reactors_[k].setOmega(omega.GetRow(k));

	BzzVector aux(omega.Columns());
	for (int k=1;k<=NR_;k++)
	{
		aux = omega.GetRow(k);
		reactors_[k].setOmega(aux);
	}
}

void OpenSMOKE_KPP_ReactorNetwork::ExtractMassFractions(const int jSpecies, BzzVector &omega) const
{
	for (int k=1;k<=NR_;k++)
		omega[k] = reactors_[k].omega()[jSpecies];
}

void OpenSMOKE_KPP_ReactorNetwork::ExtractMassFractions(BzzVector &omega) const 
{
	int count = 1;
	for (int k=1;k<=NR_;k++)
		for (int j=1;j<=mix_[0].NumberOfSpecies();j++)
			omega[count++] = reactors_[k].omega()[j];
}

void OpenSMOKE_KPP_ReactorNetwork::ExtractMassFractions(BzzMatrix &omega) const
{
	for (int k=1;k<=NR_;k++)
		for (int j=1;j<=mix_[0].NumberOfSpecies();j++)
			omega[k][j] = reactors_[k].omega()[j];
}

void OpenSMOKE_KPP_ReactorNetwork::Residuals(BzzVector& x, BzzVector& f)
{
/*	if(procrank_ == 0)
	    std::cout << "Distribution of mass fractions..." << std::endl;
	
	if(procrank_ > 0)
	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	        reactors_[k].setOmegaFromNetwork(k, x, *this); 
	}

	SendAndBroadcastOmega();

	if(procrank_ == 0)
	    std::cout << "Residual evaluation..." << std::endl;

	if(procrank_ > 0)
	{
	    for (int k=1;k<=NR_P[procrank_];k++)
	    {
		reactors_[k].Residuals(*this);

		for(int i = 1; i <= mix_[0].NumberOfSpecies(); i++)
		    EqResidualsLoc[k][i] = reactors_[k].residuals()[i];
	    }
	}

	communicator_->SendVectorToMaster(EqResidualsLoc, EqResiduals, mix_[0].NumberOfSpecies());

	if(procrank_ == 0)
	{
	    for(int k = 1; k <= NR_; k++)
	    {
	        int position = (k-1)*mix_[0].NumberOfSpecies();
		for(int j = 1; j <= mix_[0].NumberOfSpecies(); j++)
		    f[j + position] = EqResiduals[k][j];
	    }
	}*/
}

void OpenSMOKE_KPP_ReactorNetwork::Residuals(BzzVector& x, Vec& f)
{
	if(procrank_ == 0)
	    std::cout << "Distribution of mass fractions..." << std::endl;
	
	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	    {
		int index = k + offset[procrank_];
	        reactors_[k].setOmegaFromNetwork(index, x, *this); 
	    }
	}


	SendAndBroadcastOmega();

	if(procrank_ == 0)
	    std::cout << "Residual evaluation..." << std::endl;

	{
	    int counter = 0;
	    VecGetArray(f, &value);

	    for (int k=1;k<=NR_P[procrank_];k++)
	    {
		reactors_[k].Residuals(*this);
		for(int i = 1; i <= mix_[0].NumberOfSpecies(); i++)
		{
		    value[counter] = reactors_[k].residuals()[i];
		    counter++;
		}
	    }

	    VecRestoreArray(f, &value);
	}
/*
	{
	    int counter = 0;
	    nResiduals_ = NR_P[procrank_] * NumberOfSpecies();
	    double *value = new double[NR_P[procrank_] * mix_[0].NumberOfSpecies()];
	    for (int k=1;k<=NR_P[procrank_];k++)
	    {
		reactors_[k].Residuals(*this);
		for(int i = 1; i <= mix_[0].NumberOfSpecies(); i++)
		{
		    value[counter] = reactors_[k].residuals()[i];
		    counter++;
		}
	    }
	    

	    VecSetValues(f, nResiduals_, place_array, value, INSERT_VALUES);
	    delete [] value;
	}	

	VecAssemblyBegin(f);
	VecAssemblyEnd(f);*/

/*	double pippo;
	VecNorm(f, NORM_INFINITY, &pippo);
	if(procrank_ == 0)
	    std::cout << pippo << std::endl;


	getchar();*/
}

void OpenSMOKE_KPP_ReactorNetwork::Residuals(Vec& x, Vec& f)
{
	if(procrank_ == 0)
	    std::cout << "Distribution of mass fractions..." << std::endl;




	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	    {
	        reactors_[k].setOmegaFromNetwork(k, x, *this); 
	    }
	}


	SendAndBroadcastOmega();

	if(procrank_ == 0)
	    std::cout << "Residual evaluation..." << std::endl;


	{
	    int counter = 0;
	    nResiduals_ = NR_P[procrank_] * NumberOfSpecies();
	    double *value = new double[NR_P[procrank_] * mix_[0].NumberOfSpecies()];
	    for (int k=1;k<=NR_P[procrank_];k++)
	    {
		reactors_[k].Residuals(*this);
		for(int i = 1; i <= mix_[0].NumberOfSpecies(); i++)
		{
		    value[counter] = reactors_[k].residuals()[i];
		    counter++;
		}
	    }
	    

	    VecSetValues(f, nResiduals_, place_array, value, INSERT_VALUES);
	    delete [] value;
	}	

	VecAssemblyBegin(f);
	VecAssemblyEnd(f);
}

void OpenSMOKE_KPP_ReactorNetwork::ResidualsAnalysis()
{
	iteration_++;
	globalTime_ += deltat;

	if(data_.networkStatus() == KPP_NETWORK_STATUS_GLOBALODE || data_.networkStatus() == KPP_NETWORK_STATUS_GLOBALNLS)
	    SendAndBroadcastOmega();

	residuals_->Analysis();

	if(procrank_ == 0)
	    residuals_->WriteResidualsOnVideo();
}

void OpenSMOKE_KPP_ReactorNetwork::WriteReactorNetworkData()
{

	WriteReactorNetworkDataArray();	

	if(procrank_ == 0)
	{
	    std::string fileName = data_.nameFolderOutput() + "/ReactorNetwork.out";
	    std::ofstream fOutput(fileName.c_str(), std::ios::out);
	    fOutput.setf(std::ios::scientific);

/*	fOutput << "*********************************************************" << std::endl;
	fOutput << "*       Diagonal of Convection-Diffusion Matrix C       *" << std::endl;
	fOutput << "*********************************************************" << std::endl;
	for (int k=1;k<=NR_;k++)
	{
		fOutput << std::setw(8) <<  std::left << k;
		fOutput << std::setw(16) << std::left << reactors_[k].M();
		fOutput << std::endl;
	}
*/
/*	fOutput << std::endl;
	fOutput << "*********************************************************" << std::endl;
	fOutput << "*             Convection-Diffusion Matrix C             *" << std::endl;
	fOutput << "*********************************************************" << std::endl;
*/	    {
		double* ptrVal;
		int i, j;
		double val;

		C_.BeginScanning();
		while(ptrVal = C_.Scanning(&i,&j,&val))
		{
			fOutput << std::setw(8)  << std::left << i;
			fOutput << std::setw(8)  << std::left << j;
			fOutput << std::setw(16) << std::left << *ptrVal;
			fOutput << std::endl;
		}
	    }
/*
	fOutput << std::endl;
	fOutput << "*********************************************************" << std::endl;
	fOutput << "*                     External feed                     *" << std::endl;
	fOutput << "*********************************************************" << std::endl;
	for (int k=1;k<=NR_;k++)
	{
		fOutput << std::setw(8)  << std::left << k;
		fOutput << std::setw(16) << std::left << reactors_[k].fInTot();
		fOutput << std::endl;
	}
	*/
/*	fOutput << std::endl;
	fOutput << "*********************************************************" << std::endl;
	fOutput << "*                    External outflow                   *" << std::endl;
	fOutput << "*********************************************************" << std::endl;
	for (int k=1;k<=NR_;k++)
	{
		fOutput << std::setw(8)  << std::left << k;
		fOutput << std::setw(16) << std::left << reactors_[k].fOutTot();
		fOutput << std::endl;
	}

	fOutput << std::endl;
	fOutput << "*********************************************************" << std::endl;
	fOutput << "*          Mass     MassFlow     Mass Umbalances         *" << std::endl;
	fOutput << "*********************************************************" << std::endl;
	for (int k=1;k<=NR_;k++)
	{
		fOutput << std::setw(8)  << std::left << k;
		fOutput << std::setw(16) << std::left << reactors_[k].mass();
		fOutput << std::setw(16) << std::left << reactors_[k].MassFlowIn();
		fOutput << std::setw(16) << std::left << reactors_[k].MassUmbalance();
		fOutput << std::endl;
	}
*/
	    fOutput.close();
	}
}

void OpenSMOKE_KPP_ReactorNetwork::MassFlowErrors()
{
	double *meanError, *maxError;
	MPI::Status status;

	meanError = new double[numworkers_ + 1];
	maxError = new double[numworkers_ + 1];

	for(int i = 0; i <= numworkers_; i++)
	{
	    meanError[i] = 0;
	    maxError[i] = 0;
	}

        for(int k=1; k <= NR_P[procrank_]; k++)
        {
            double mIn  = reactors_[k].MassFlowIn();
            double mOut = reactors_[k].MassFlowOut();
        if (mIn+mOut<1.e-16)
            continue;

            double relative_error = fabs((mOut-mIn)/(mIn+mOut));
            meanError[procrank_]+=relative_error;

            if (relative_error>maxError[procrank_])
                maxError[procrank_] = relative_error;
        }

	if(procrank_ > 0)
	{
	    mtype = FROM_WORKER;
	    MPI::COMM_WORLD.Send(&meanError[procrank_], 1, MPI::DOUBLE, MASTER, mtype);
	    MPI::COMM_WORLD.Send(&maxError[procrank_], 1, MPI::DOUBLE, MASTER, mtype);
	}

        if(procrank_ == 0)
        {
            mtype = FROM_WORKER;
            for(int i = 1; i <= numworkers_; i++)
            {
                source = i;

                MPI::COMM_WORLD.Recv(&meanError[source], 1, MPI::DOUBLE, source, mtype, status);
                MPI::COMM_WORLD.Recv(&maxError[source], 1, MPI::DOUBLE, source, mtype, status);
            }

            for(int i = 1; i <= numworkers_; i++)
            {
                meanError[0] += meanError[i];
                if(maxError[i] > maxError[0])
                    maxError[0] = maxError[i];
            }

            meanError[0] /= double(NR_);

            std::cout << "        Mean relative error in mass balances: " << meanError[0] << std::endl;
            std::cout << "        Max  relative error in mass balances: " <<  maxError[0] << std::endl;
        }

        delete [] meanError;
        delete [] maxError;

}

KPP_Network_Status OpenSMOKE_KPP_ReactorNetwork::status() const
{
	return data_.networkStatus(); 
}

void OpenSMOKE_KPP_ReactorNetwork::CleanMassFractions(BzzVector& omega, const int nr, const int nc, const double threshold)
{

	for (int i=1;i<=nr;i++)
	{
	    int k=(i-1)*nc;
		
	    double sum = 0.;
	    for (int j=k+1;j<=k+nc;j++)
	    {
		// In case of strong negative value
		if (omega[j]<threshold && procrank_ == 0)
		{
		    std::cout << "Negative mass fraction: " << omega[j] << std::endl;
		    std::cout << "Press enter to exit..." << std::endl;
		    getchar();
		    exit(-1);
		}

		// In case of sligthly negative value
		if (omega[j]<0.)	omega[j] = 0.;

		// Sum of mass fractions
		sum+=omega[j];
	    }

	    // Normalization
	    for (int j=k+1;j<=k+nc;j++)
		omega[j] /= sum;
	}
}

int OpenSMOKE_KPP_ReactorNetwork::SolutionNLS(BzzVector& x, BzzVector& direction, const bool jacobianFlag)
{
	MPI::Status status;
        // CPU time
	if(procrank_ == 0)
            cpu_->SetStartTimeNLSLocal();
        
        // Counters
        globalCounterNLS_++;
        localCounterNLS_++;        
        int linearSystemIterations      = 0;
        double linearSystemNorm         = 0.;        
        
	// Cleaning
	CleanMassFractions(x, NR_, NumberOfSpecies(), -1);

	// Statistics
	if(procrank_ == 0)
	statisticsConvective_->Reset();


	// Communication: Update transport
	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	    {
		int index = k + offset[procrank_];
		reactors_[k].setOmegaFromNetwork(index, x, *this); 
	    }
	}

	// Assembling matrices
        BzzMatrix tmpMatrix(NumberOfSpecies(), NumberOfSpecies());
        BzzMatrix diagonalBlockMatrix(NumberOfSpecies(), NumberOfSpecies());
	{
            for(int k=1;k<=NR_P[procrank_];k++)
	    {
		double dummy_deltat=0.;
		reactors_[k].AssemblingLocalContribution(dummy_deltat, jacobianFlag, tmpMatrix, diagonalBlockMatrix);
		reactors_[k].AssemblingNonLocalRHS(dummy_deltat, x);
	    }
	}

	// Update global rhs
	UpdateRHS();

	// Update global matrix
	if (jacobianFlag == true)
	{
                globalCounterNLSJacobians_++;
                        
		if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
		{
		    if(procrank_ == 0)
			ErrorMessage("Pardiso not implemented in the parallel form");
		}

		else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
		{
		    if(procrank_ == 0)
			ErrorMessage("MUMPS not implemented in the parallel form");
		}

		else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
		{
			lisMatrixGlobal->ResetCounters();

			UpdateLISMatrix();

			lisMatrixGlobal->NumericalFactorization();

			lisMatrixGlobal->EnableInitialSolution();

		}

		else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL)
		{
		    if(procrank_ == 0)
			ErrorMessage("Gauss-Seidel not implemented in the parallel form. Use LIS instead");
		}
	}

	// Linear system solution
	if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
	{
	    if(procrank_ == 0)
		direction = lastGlobalNLSSolution_;

	    lisMatrixGlobal->SetInitialSolution(direction);

	    lisMatrixGlobal->Solve(localRHS_, direction);

	    if(procrank_ == 0)
		lastGlobalNLSSolution_ = direction;
	}

	// New values
	if(procrank_ == 0)
	    Sum(x, direction, &auxVector_NRxNC_1);
	
	// Cleaning machine errors
	if(procrank_ == 0)
	{
	    CleanVectorOnTheBottom(0., -1.e-16, auxVector_NRxNC_1);
	    CleanVectorOnTheTop(1., 1.+1.e-16, auxVector_NRxNC_1);
	}

	// Correcting time step
	if(procrank_ == 0)
	    statisticsConvective_->UnsetRobust();

	if(procrank_ == 0)
	{
	    double a = statisticsConvective_->ReductionOfTimeStep(auxVector_NRxNC_1, x);
	    if (statisticsConvective_->iTimeCorrectedMin_ > 0 ||  statisticsConvective_->iTimeCorrectedMax_ > 0)
	    {
		a *= data_.GlobalNLS_SafetyReductionCoefficient();
		Product(a,&direction);
	    }
        
            // CPU time
            cpu_->SetEndTimeNLSLocal();
            cpu_->SetEndTimeNLSGlobal();
            cpu_->SetEndTimeAllGlobal();        
            cpu_->WriteOnFile();
        
            // Write on file
            fGlobalNLS_ << std::setw(12) << std::left << iteration_;
            fGlobalNLS_ << std::setw(12) << std::left << globalCounterNLS_;
            fGlobalNLS_ << std::setw(12) << std::left << localCounterNLS_;
            fGlobalNLS_ << std::setw(16) << std::left << 0;
            fGlobalNLS_ << std::setw(16) << std::left << a;
            fGlobalNLS_ << std::setw(12) << std::left << globalCounterNLSJacobians_;         
            fGlobalNLS_ << std::setw(12) << std::left << linearSystemIterations;
            fGlobalNLS_ << std::setw(16) << std::left << linearSystemNorm;
            fGlobalNLS_ << std::endl;           
	}

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= NumberOfEquations(); i++)
	    {
		x_ptr[i] = x[i];
//		direction_ptr[i] = direction[i];
	    }
	}

	MPI::COMM_WORLD.Bcast(&x_ptr[1], NumberOfEquations(), MPI::DOUBLE, MASTER);
//	MPI::COMM_WORLD.Bcast(&direction_ptr[1], NumberOfEquations(), MPI::DOUBLE, MASTER);


	if(procrank_ > 0)
	{
	    for(int i = 1; i <= NumberOfEquations(); i++)
	    {
		x[i] = x_ptr[i];
//		direction[i] = direction_ptr[i];
	    }
	}


	return 0;
}

int OpenSMOKE_KPP_ReactorNetwork::SolutionODE(BzzVector& omegaOld, BzzVector& step, double& dt, const bool jacobianFlag)
{
	MPI::Status status;

        // CPU time
	if(procrank_ == 0)
	    cpu_->SetStartTimeODELocal();
        
        // Counters
        globalCounterODE_++;
        localCounterODE_++;        
        int linearSystemIterations      = 0;
        double linearSystemNorm         = 0.;

	// Cleaning
	CleanMassFractions(omegaOld, NR_, NumberOfSpecies(), -1.e-6);

	// Statistics
	if(procrank_ == 0)
	    statisticsConvective_->Reset();

	// Communication: Update transport
	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	    {
		int index = k + offset[procrank_];
		reactors_[k].setOmegaFromNetwork(index, omegaOld, *this); 
	    }
	}

	// Assembling matrix and rhs
        BzzMatrix tmpMatrix(NumberOfSpecies(), NumberOfSpecies());
        BzzMatrix diagonalBlockMatrix(NumberOfSpecies(), NumberOfSpecies());

	{
	    for(int k=1;k<=NR_P[procrank_];k++)
	    {
		reactors_[k].AssemblingLocalContribution(dt, jacobianFlag, tmpMatrix, diagonalBlockMatrix);
		reactors_[k].AssemblingNonLocalRHS(dt, omegaOld);
	    }
	}

	if(procrank_ == 0) std::cout << "   Updating RHS...";

	UpdateRHS();

	if(procrank_ == 0) std::cout << "Done!" << std::endl;

	// Update global matrix
	if (jacobianFlag == true)
	{
                globalCounterODEJacobians_++;
                
		if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
		{
		    if(procrank_ == 0)
			ErrorMessage("MUMPS not implemented in the parallel form");
		}

		else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
		{
		    if(procrank_ == 0)
			ErrorMessage("MUMPS not implemented in the parallel form");
		}

		else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
		{
			lisMatrixGlobal->ResetCounters();

			if(procrank_ == 0) std::cout << "   Updating Matrix...";

			UpdateLISMatrix();

			if(procrank_ == 0) std::cout << "Done!" << std::endl;


			lisMatrixGlobal->NumericalFactorization();

			//	lisMatrixGlobal->UnsetMessage();
			lisMatrixGlobal->EnableInitialSolution();
                        
                        //Getting Boolean Sparsity Pattern
                        if (iBooleanMatrix_ == false)
                        {
//                            lisMatrixGlobal->GetBooleanSparsityPattern(NR_ * mix_->NumberOfSpecies(), NR_);
                            iBooleanMatrix_ = true;
                        }
		}

		else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL)
		{
		    if(procrank_ == 0)
			ErrorMessage("OpenSMOKE_GaussSeidel not available in parallel version!");
		}
	}
        
	// Linear system solution

	if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
	{
		if(procrank_ == 0)
		    step = lastGlobalODESolution_;

		lisMatrixGlobal->SetInitialSolution(step);

		lisMatrixGlobal->Solve(localRHS_, step);
		
		if(procrank_ == 0)
		    lastGlobalODESolution_ = step;

		if(lisMatrixGlobal->convergence_index() == false)
	    	    return -1;
	}
     
        // New solution
	if(procrank_ == 0)
	    Sum(omegaOld, step, &auxVector_NRxNC_1);


	// Cleaning machine errors
	if(procrank_ == 0)
	{
	    CleanVectorOnTheBottom(0., -1.e-16, auxVector_NRxNC_1);
	    CleanVectorOnTheTop(1., 1.+1.e-16, auxVector_NRxNC_1);
	}
	
	// Correction of time step
	if(procrank_ == 0)
	{
	    statisticsConvective_->SetRobust();
	}
	if(procrank_ == 0)
	{
	    double a = statisticsConvective_->ReductionOfTimeStep(auxVector_NRxNC_1, omegaOld);
	    if (statisticsConvective_->iTimeCorrectedMin_ > 0 ||  statisticsConvective_->iTimeCorrectedMax_ > 0)
	    {
		a  *= data_.GlobalODE_SafetyReductionCoefficient();
		dt *= a;
		Product(a, &step);
		Sum(&omegaOld, step);
	    }
	    else
		omegaOld = auxVector_NRxNC_1;


	    // Additional Analysis
	    statisticsConvective_->Analysis(dt, a, omegaOld);
       
            fGlobalODE_ << std::setw(12) << std::left << iteration_;
            fGlobalODE_ << std::setw(12) << std::left << globalCounterODE_;
            fGlobalODE_ << std::setw(12) << std::left << localCounterODE_;
            fGlobalODE_ << std::setw(16) << std::left << dt;
            fGlobalODE_ << std::setw(16) << std::left << a;
            fGlobalODE_ << std::setw(12) << std::left << globalCounterODEJacobians_;         
            fGlobalODE_ << std::setw(12) << std::left << linearSystemIterations;
            fGlobalODE_ << std::setw(16) << std::left << linearSystemNorm;
            fGlobalODE_ << std::endl;        
        
            // CPU time
            cpu_->SetEndTimeODELocal();
            cpu_->SetEndTimeODEGlobal();
            cpu_->SetEndTimeAllGlobal();          
            cpu_->WriteOnFile();
        }

	if(procrank_ == 0)
	{
	    for(int i = 1; i <= NumberOfEquations(); i++)
		omega_ptr[i] = omegaOld[i];
	}

	MPI::COMM_WORLD.Bcast(&omega_ptr[1], NumberOfEquations(), MPI::DOUBLE, MASTER);

	if(procrank_ > 0)
	{
	    for(int i = 1; i <= NumberOfEquations(); i++)
		omegaOld[i] = omega_ptr[i];
	}

	MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, MASTER);
        
	return 0;
}

int OpenSMOKE_KPP_ReactorNetwork::SequenceCSTR()
{
        // CPU time
	if(procrank_ == 0)
	{
	    cpu_->SetSequence();
	    cpu_->SetStartTimeSequenceGlobal();
	}
	
	LoopIndex_++;
                
	data_.SetNetworkStatus(KPP_NETWORK_STATUS_SEQUENTIAL_CSTR);
	nLocalCSTRSequence_ = 0;

	for (int k=1;k<=NR_P[procrank_];k++)
		reactors_[k].ResetStatus();

	NormInfConvCriterium = pow(data_.NormInfConvergence(), (1./double(LoopIndex_)));
	Norm1ConvCriterium = pow(data_.Norm1Convergence(), (1./double(LoopIndex_)));
	ReactIterCriterium = data_.MinimumCSTRIterations()/LoopIndex_;

	double timeStartGlobal;

	if(procrank_ == 0)
	    timeStartGlobal = MPI::Wtime();

	if(residuals_->norm2() <= 1.e-11 && LoopIndex_ != 1)
	    return 0;

	for (int i=1; i<=data_.SequenceCSTR_MaxIterations(); i++)		
	{
		if(procrank_ == 0)
		    std::cout << "Iteration no. " << i << std::endl;

		// CPU Time
		if(procrank_ == 0)
		    cpu_->SetStartTimeSequenceLocal();
                
                // Solve the network
		SolveSequenceCSTR();
                
		SendAndBroadcastOmega();
		
		// Residual analysis
		ResidualsAnalysis();
		
		// Analysis of reactors
		int info = AnalysisSequenceCSTR();

                // CPU time
		if(procrank_ == 0)
		{
                    cpu_->SetEndTimeSequenceLocal();
                    cpu_->SetEndTimeSequenceGlobal();
                    cpu_->SetEndTimeAllGlobal();                  
                    cpu_->WriteOnFile();
		}

		if (info == 0)	return info;	// Means that convergence is reached
		if (info  < 0)	return info;	// Means that problems were met
	}

	return 1;	// Means that the std::maximum number of iterations was reached
}

void OpenSMOKE_KPP_ReactorNetwork::ReadBackupFile(const std::string fileName)
{
	std::vector<double*> massFractions;
	massFractions.resize(NR_ + 1);
	for(int i = 1; i <= NR_; i++)
	    massFractions[i] = new double[NumberOfSpecies() + 1];
	
	if (data_.iTraditionalLink() == false)
	{
	    if(procrank_ == 0)
            {
	        char name[Constants::NAME_SIZE];

	        int numberReactorsBackup;
	        int numberSpeciesBackup;
	        BzzVector massFractionsBzz;

	        BzzLoad fBackup('*', fileName.c_str());
	        fBackup >> numberReactorsBackup;
	        fBackup >> numberSpeciesBackup;


	        std::cout << "Number of reactors (current): " << NR_ << std::endl;
	        std::cout << "Number of reactors (backup):  " << numberReactorsBackup << std::endl;
	        std::cout << "Number of species (current):  " << NumberOfSpecies() << std::endl;
	        std::cout << "Number of species (backup):   " << numberSpeciesBackup << std::endl;

	        if (NR_ != numberReactorsBackup)
		    ErrorMessage("The number of reactors in backup file does not fit!");

	        BzzVectorInt indices(numberSpeciesBackup);
	        BzzVector omegaLocal(NumberOfSpecies());
	
	        int countRecognized=0;
	        for(int j=1;j<=numberSpeciesBackup;j++)
	        {
		    fBackup.fileLoad.read((char*) name, sizeof(name));
		    indices[j] = mix_[0].recognize_species_without_exit(name);
		    if (indices[j] > 0)	countRecognized++;
	        }

	        fBackup >> massFractionsBzz;
	        fBackup.End();

	        std::cout << "Number of species recognized:   " << countRecognized << "/" << NumberOfSpecies() << std::endl;

	        double sum_max = 0.;
	        double sum_min = 1.;

	        int i=1;
	        for(int k=1;k<=NR_;k++)
	        {
		    omegaLocal=0.;
		    for(int j=1;j<=numberSpeciesBackup;j++)
		    {
			if (indices[j]>0)	omegaLocal[indices[j]] = massFractionsBzz[i];
			i++;
		    }
		
		    double sum = omegaLocal.GetSumElements();
		    if (sum > sum_max) sum_max = sum;
		    else if (sum < sum_min) sum_min = sum;
		    Product(1./sum, &omegaLocal);
		
		    for(int j = 1; j <= NumberOfSpecies(); j++)
		        massFractions[k][j] = omegaLocal[j];
//		    reactors_[k].setOmega(omegaLocal);
	        }

	        std::cout << "Min/max sums:   " << sum_min << "/" << sum_max << std::endl;
	    }

	}

	if(data_.iTraditionalLink() == true)
        {
	    if(procrank_ == 0)
	    {
        	BzzVectorInt fromOriginalToCluster;
        	BzzVector volumes;
		std::string TopologyFile = data_.nameFolderInput() + "/ClusteringTopology.out";
        	BzzLoad fClusteringTopology(TopologyFile.c_str());
        	fClusteringTopology >> fromOriginalToCluster;
        	fClusteringTopology >> volumes;
        	fClusteringTopology.End();

        	BzzMatrix mass_fractions;
        	int numberOfSpecies, originalNumberOfReactors;
		std::string BackupPath = data_.nameInputBackupFile() + ".start";
        	BzzLoad  fTraditionalBackup('*', BackupPath.c_str());

         	fTraditionalBackup >> originalNumberOfReactors;
        	fTraditionalBackup >> numberOfSpecies;

        	fTraditionalBackup >> mass_fractions;
        	fTraditionalBackup.End();

        	BzzMatrix mass_fractions_star(NR_, NumberOfSpecies());
        	BzzVector volume_star(NR_);

         	for(int k=1;k<=mass_fractions.Rows();k++)
        	{
                    int i = fromOriginalToCluster[k];
                    volume_star[i] += volumes[k];
                    for(int j=1;j<=NumberOfSpecies();j++)
                        mass_fractions_star[i][j] += mass_fractions[k][j]*volumes[k];
                }


        	for(int k=1;k<=NR_;k++)
                    for(int j=1;j<=NumberOfSpecies();j++)
                        mass_fractions_star[k][j] /= volume_star[k];

                double sum_mean = 0.;
                double sum_max = -1e16;
                double sum_min =  1e16;

         	for(int k=1;k<=NR_;k++)
                {
                    BzzVector aux = mass_fractions_star.GetRow(k);
                    for(int j=1;j<=NumberOfSpecies();j++)
                    {
                        double epsilon = 1.e-8;
                        if (aux[j] < 0.)
                        {
                            if (aux[j] < 0.-epsilon)
                            {
                                std::cout << "Negative mass fraction: Reactor " << k << " Species " << j << " Omega " << aux[j] << std::endl;
                                std::cout << "Press enter to exit..." << std::endl;
                                getchar();
                                exit(-1);
                            }
                            else
                                aux[j] = 0.;
                        }
                    }
                    double epsilon = 1.e-5;


                    double sum = aux.GetSumElements();

                    if ( (sum < 1.-epsilon) || (sum > 1.+epsilon) )
                    {
                        std::cout << "Wrong sum of mass fractions: Reactor " << k << " Sum " << sum << std::endl;
                        std::cout << "Press enter to exit..." << std::endl;
                        getchar();
                        exit(-1);
                    }

                    sum_mean += sum;
                    if (sum>sum_max)    sum_max=sum;
                    if (sum<sum_min)    sum_min=sum;

                    aux /= sum;
		    for(int j = 1; j <= NumberOfSpecies(); j++)
		        massFractions[k][j] = aux[j];

//                    reactors_[k].setOmega(aux);
                }

                std::cout << "Statistics from Backup" << std::endl;
                std::cout << std::setprecision(8) << std::fixed << "Mean sums: " << sum_mean/double(NR_) << std::endl;
                std::cout << std::fixed << std::setprecision(8) << "Max sums:  " << sum_max << std::endl;
                std::cout << std::fixed << std::setprecision(8) << "Min sums:  " << sum_min << std::endl;
	    }
        }


	for(int k = 1; k <= NR_; k++)
	    MPI::COMM_WORLD.Bcast(&massFractions[k][1], NumberOfSpecies(), MPI::DOUBLE, MASTER);
	
	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    BzzVector omegaLocal(NumberOfSpecies());
	    int kglob = k + offset[procrank_];
	    for(int j = 1; j <= NumberOfSpecies(); j++)
		omegaLocal(j) = massFractions[kglob][j];
	    reactors_[k].setOmega(omegaLocal);
	}
	

	for(int k = 0; k <= NR_; k++)
	    delete [] massFractions[k];

	//getchar();
			
}

void OpenSMOKE_KPP_ReactorNetwork::WriteBackupFile(const std::string fileName)
{
	SendAndBroadcastOmega();



	if(procrank_ == 0)
	{
	    char name[Constants::NAME_SIZE];

	    int i=1;
	    for(int k=1;k<=NR_;k++)
		for(int j=1;j<=NumberOfSpecies();j++)
		    auxVector_NRxNC_1[i++] = OmegaGlob_[k][j];

	    BzzSave fBackup('*', fileName.c_str());
	    fBackup << NR_;
	    fBackup << NumberOfSpecies();
	    for(int j=1;j<=NumberOfSpecies();j++)
	    {
		strcpy(name, mix_[0].names[j].c_str());
		fBackup.fileSave.write((char*) name, sizeof(name));
	    }
	    fBackup << auxVector_NRxNC_1;
	    fBackup.End();

	    if (data_.iTraditionalLink() == true)
            {
        	BzzVectorInt fromOriginalToCluster;
        	BzzVector volumes;
		std::string TopologyFile = data_.nameFolderInput() + "/ClusteringTopology.out";
        	BzzLoad fClusteringTopology(TopologyFile.c_str());
        	fClusteringTopology >> fromOriginalToCluster;
        	fClusteringTopology >> volumes;
        	fClusteringTopology.End();

        	BzzMatrix mass_fractions(fromOriginalToCluster.Size(), NumberOfSpecies());
        	for(int k=1;k<=fromOriginalToCluster.Size();k++)
		{
		    BzzVector omegaLoc(NumberOfSpecies());
		    for(int j = 1; j <= NumberOfSpecies(); j++)
			omegaLoc[j] = OmegaGlob_[fromOriginalToCluster[k]][j];

                    mass_fractions.SetRow(k, omegaLoc);
		}

		std::string BackupPath = data_.nameOutputBackupFile() + ".end";
        	BzzSave  fTraditionalBackup('*', BackupPath.c_str());
        	fTraditionalBackup << NumberOfSpecies();
        	fTraditionalBackup << fromOriginalToCluster.Size();
         	fTraditionalBackup << mass_fractions;
        	fTraditionalBackup.End();
	    }
        }
}

void OpenSMOKE_KPP_ReactorNetwork::DistributeCells()
{
	 int avecells = NR_/nprocs_;
    	 int extra = NR_%nprocs_;

	 offset[0] = 0;
   	 offset[nprocs_] = NR_;
	
  	  for(int dest = 0; dest <= numworkers_; dest++)
  	  {
   	     NR_P[dest] = (dest < extra) ? avecells + 1 : avecells;
   	     if(dest != numworkers_)
   	         offset[dest+1] = offset[dest] + NR_P[dest];
  	  }

	 ReactorAlternateRank();
	
}

void OpenSMOKE_KPP_ReactorNetwork::TopologyFakeReading()
{
	std::ifstream fInput(data_.nameTopologyFile().c_str(), std::ios::in);

	unsigned int dummy;
	fInput >> dummy;
	if (dummy != NR_)
		ErrorMessage("The number of reactors does not fit!");
	
	Counter_Glob = new int[NR_ + 1];
	counter = new int[numworkers_ + 1];

	if(procrank_ == 0)
	{
	    CountElem = 0;
	    reactors_fake = new OpenSMOKE_KPP_SingleReactor[NR_ + 1];
            for(int k = 1; k <= NR_; k++)
            {
            	reactors_fake[k].Setup(k, k + offset[procrank_], &mix_[0], jReduced);
            	reactors_fake[k].SetDataManager(&data_);
            }
            for (int k = 1; k <= NR_; k++)
                {
                    reactors_fake[k].ReadReactorProperties(fInput);
                    reactors_fake[k].ReadTopology(fInput);
                    Counter_Glob[k] = CountElem;
                    CountElem += reactors_fake[k].count();
                    for(int p = 0; p <= numworkers_; p++)
                    {
                        if(k == offset[p])  counter[p] = CountElem;
                    }
                }

            counter[0] = 0;

	    delete [] reactors_fake;
	}

	MPI::COMM_WORLD.Bcast(counter, numworkers_ + 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(Counter_Glob, NR_+1, MPI::INT, MASTER);

	fInput.close();

}

void OpenSMOKE_KPP_ReactorNetwork::ReactorAlternateRank()
{
	ReactID = new int[NR_ + 1];

    	ReactID = new int[NR_ + 1];
    	for(int i = 1; i <= NR_; i++)
    	{
            for(int j = 0; j <= numworkers_; j++)
            {
            	if(i > offset[j] && i <= offset[j+1])
                ReactID[i] = j;
            }
        }

}

void OpenSMOKE_KPP_ReactorNetwork::DeleteTopologyParametersSize()
{
	communicator_->DeleteArray(ReactOutSize);
	communicator_->DeleteArray(ReactInSize);
	communicator_->DeleteArray(ReactNeighSize);

	communicator_->DeleteVector(React_out, NR_);
	communicator_->DeleteVector(React_in, NR_);
	communicator_->DeleteVector(React_cOut, NR_);
	communicator_->DeleteVector(React_cIn, NR_);
	communicator_->DeleteVector(React_diffusion, NR_);
	communicator_->DeleteVector(React_neighbours, NR_);
}

void OpenSMOKE_KPP_ReactorNetwork::TopologyParametersSize()
{
	communicator_->InitializeArray(ReactOutSize);
	communicator_->InitializeArray(ReactInSize);
	communicator_->InitializeArray(ReactNeighSize);

	for(int k = 1; k <= NR_P[procrank_]; k++)
        {
	    int index = k + offset[procrank_];
            ReactOutSize[index] = reactors_[k].out().Size();
            ReactInSize[index] = reactors_[k].in().Size();
            ReactNeighSize[index] = reactors_[k].neighbours().Size();
        }

	communicator_->GatherArray(ReactOutSize, NR_P, offset);
	communicator_->GatherArray(ReactInSize, NR_P, offset);
	communicator_->GatherArray(ReactNeighSize, NR_P, offset);
	
	communicator_->InitializeVector(React_out, ReactOutSize);
	communicator_->InitializeVector(React_in, ReactInSize);
	communicator_->InitializeVector(React_cOut, ReactOutSize);
	communicator_->InitializeVector(React_cIn, ReactInSize);
	communicator_->InitializeVector(React_diffusion, ReactNeighSize);
	communicator_->InitializeVector(React_neighbours, ReactNeighSize);
}

void OpenSMOKE_KPP_ReactorNetwork::TopologyParameters()
{
	
	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = k + offset[procrank_];
	    for(int j = 1; j <= ReactOutSize[index]; j++)
	    {
		React_out[index][j] = reactors_[k].out()[j];
		React_cOut[index][j] = reactors_[k].cOut()[j];
	    }
	    for(int j = 1; j <= ReactNeighSize[index]; j++)
	    {
		React_diffusion[index][j] = reactors_[k].diffusion()[j];
		React_neighbours[index][j] = reactors_[k].neighbours()[j];
	    }
	    for(int j = 1; j <= ReactInSize[index]; j++)
	    {
		React_in[index][j] = reactors_[k].in()[j];
		React_cIn[index][j] = reactors_[k].cIn()[j];
	    }
	}

	communicator_->GatherVector(React_out, ReactOutSize);
	communicator_->GatherVector(React_cOut, ReactOutSize);
	communicator_->GatherVector(React_diffusion, ReactNeighSize);
	communicator_->GatherVector(React_neighbours, ReactNeighSize);
	communicator_->GatherVector(React_in, ReactInSize);
	communicator_->GatherVector(React_cIn, ReactInSize);
}


void OpenSMOKE_KPP_ReactorNetwork::FeedsAndOutputs()
{	
	
	MPI::Status status;

	{
	    for (int k=1;k<=NR_P[procrank_];k++)
            {
		int index = k + offset[procrank_];
                if (reactors_[k].IsExternalFeedReactor() == true)	jExternalFeedReactors.Append(index);
                if (reactors_[k].IsExternalOutputReactor() == true)	jExternalOutputReactors.Append(index);
            }

            jEFR_size = jExternalFeedReactors.Size();
            jEOR_size = jExternalOutputReactors.Size();

	    mtype = FROM_WORKER;
            MPI::COMM_WORLD.Send(&jEFR_size, 1, MPI::INT, MASTER, mtype);
            MPI::COMM_WORLD.Send(&jEOR_size, 1, MPI::INT, MASTER, mtype);

	    communicator_->InitializeArray(jExternalFeedReactors_Loc, jEFR_size);
	    communicator_->InitializeArray(jExternalOutputReactors_Loc, jEOR_size);

	    for(int i = 1; i <= jEFR_size; i++)
		jExternalFeedReactors_Loc[i] = jExternalFeedReactors[i];

	    for(int i = 1; i <= jEOR_size; i++)
		jExternalOutputReactors_Loc[i] = jExternalOutputReactors[i];

            MPI::COMM_WORLD.Send(&jExternalFeedReactors_Loc[1], jEFR_size, MPI::INT, MASTER, mtype);
            MPI::COMM_WORLD.Send(&jExternalOutputReactors_Loc[1], jEOR_size, MPI::INT, MASTER, mtype);
	}

	if(procrank_ == 0)
	{
	    communicator_->InitializeArray(jEFR_sizeAll, numworkers_);
	    communicator_->InitializeArray(jEOR_sizeAll, numworkers_);

	    mtype = FROM_WORKER;
	    for (int p = 0; p <= numworkers_; p++)
	    {
                source = p;
                MPI::COMM_WORLD.Recv(&jEFR_sizeAll[source], 1, MPI::INT, source, mtype, status);
                MPI::COMM_WORLD.Recv(&jEOR_sizeAll[source], 1, MPI::INT, source, mtype, status);
	    }

	    jEFRsum = 0;
	    jEORsum = 0;
	    communicator_->InitializeArray(jEFRoffset, numworkers_);
	    communicator_->InitializeArray(jEORoffset, numworkers_);

	    for(int i = 0; i <= numworkers_; i++)
            {
                jEFRoffset[i] = jEFRsum;
                jEORoffset[i] = jEORsum;
                jEFRsum += jEFR_sizeAll[i];
                jEORsum += jEOR_sizeAll[i];
            }

	    communicator_->InitializeArray(jExternalFeedReactors_aux, jEFRsum);
	    communicator_->InitializeArray(jExternalOutputReactors_aux, jEORsum);

	    for(int p = 0; p <= numworkers_; p++)
            {
                source = p;
                MPI::COMM_WORLD.Recv(&jExternalFeedReactors_aux[1 + jEFRoffset[source]], jEFR_sizeAll[source], MPI::INT,
                                 source, mtype, status);
                MPI::COMM_WORLD.Recv(&jExternalOutputReactors_aux[1 + jEORoffset[source]], jEOR_sizeAll[source], MPI::INT,
                                 source, mtype, status);
            }

	    ChangeDimensions(jEFRsum, &jExternalFeedReactors);
            ChangeDimensions(jEORsum, &jExternalOutputReactors);

            for(int i = 1; i <= jEFRsum; i++)
            {
                jExternalFeedReactors[i] = jExternalFeedReactors_aux[i];
            }
            for(int i = 1; i <= jEORsum; i++)
            {
                jExternalOutputReactors[i] = jExternalOutputReactors_aux[i];
            }

            Sort(&jExternalFeedReactors);
            Sort(&jExternalOutputReactors);
	}
}

void OpenSMOKE_KPP_ReactorNetwork::VariableTopologyParameters()
{

	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
	    {
		int index = k + offset[procrank_];
		for(int j = 1; j <= ReactOutSize[index]; j++)
		{
		    React_cOut[index][j] = reactors_[k].cOut()[j];
		}
		for(int j = 1; j <= ReactNeighSize[index]; j++)
		{
		    React_diffusion[index][j] = reactors_[k].diffusion()[j];
		}
		for(int j = 1; j <= ReactInSize[index]; j++)
		{
		    React_cIn[index][j] = reactors_[k].cIn()[j];
		}
	    }	
	}	

	communicator_->GatherVector(React_cOut, ReactOutSize);
	communicator_->GatherVector(React_diffusion, ReactNeighSize);
	communicator_->GatherVector(React_cIn, ReactInSize);
}

void OpenSMOKE_KPP_ReactorNetwork::PrepareCDMatrixCommunication()
{
	int p;
	double M_;
	MPI::Status status;
	
	ElementBzzMatrixSparse *elem;

        if(procrank_ == 0)
        {
            elem = C().GetStartingElementInRow(network_index_);
            p = 0;
            while(elem)
            {
                p++;
                elem = elem->next;
            }
            ChangeDimensions(p-1, &iConvDiffMaster);
            ChangeDimensions(p-1, &mConvDiffMaster);

            mtype = FROM_MASTER;
            MPI::COMM_WORLD.Send(&p, 1, MPI::INT, network_ID_, mtype);

	    elem = C().GetStartingElementInRow(network_index_);
            p = 0;

            while(elem)
            {
                int j = elem->column;
                if (j == network_index_)
                {
                    M_ = elem->value;
                }
                else
                {
                    p++;
                    iConvDiffMaster[p] = j;
                    mConvDiffMaster[p] = -elem->value;
                }
                elem = elem->next;
            }

            iConvDiffMaster_aux = new int[p+1];
            mConvDiffMaster_aux = new double[p+1];

            for(int j = 1; j <= p; j++)
            {
                iConvDiffMaster_aux[j] = iConvDiffMaster[j];
                mConvDiffMaster_aux[j] = mConvDiffMaster[j];
            }

             mtype = FROM_MASTER;
             MPI::COMM_WORLD.Send(&iConvDiffMaster_aux[1], p, MPI::INT, network_ID_, mtype);
             MPI::COMM_WORLD.Send(&mConvDiffMaster_aux[1], p, MPI::DOUBLE, network_ID_, mtype);
             MPI::COMM_WORLD.Send(&M_, 1, MPI::DOUBLE, network_ID_, mtype);


             delete [] iConvDiffMaster_aux;
             delete [] mConvDiffMaster_aux;
        }

	if(procrank_ == network_ID_)
        {
	    mtype = FROM_MASTER;
	    MPI::COMM_WORLD.Recv(&p, 1, MPI::INT, MASTER, mtype, status);

	    reactors_[network_index_ - offset[procrank_]].ConvDiffLength = p;

            iConvDiffMaster_aux = new int[p];
            mConvDiffMaster_aux = new double[p];

            MPI::COMM_WORLD.Recv(&iConvDiffMaster_aux[1], p-1, MPI::INT, MASTER, mtype, status);
            MPI::COMM_WORLD.Recv(&mConvDiffMaster_aux[1], p-1, MPI::DOUBLE, MASTER, mtype, status);
            MPI::COMM_WORLD.Recv(&M_, 1, MPI::DOUBLE, MASTER, mtype, status);

            reactors_[network_index_ - offset[procrank_]].M_aux = M_;

            ChangeDimensions(p-1, &reactors_[network_index_ - offset[procrank_]].iConvDiff_aux);
            ChangeDimensions(p-1, &reactors_[network_index_ - offset[procrank_]].mConvDiff_aux);

            for(int j = 1; j <= p-1; j++)
            {
                reactors_[network_index_ - offset[procrank_]].iConvDiff_aux[j] = iConvDiffMaster_aux[j];
                reactors_[network_index_ - offset[procrank_]].mConvDiff_aux[j] = mConvDiffMaster_aux[j];
            }

	    delete [] iConvDiffMaster_aux;
	    delete [] mConvDiffMaster_aux;
	}
}

void OpenSMOKE_KPP_ReactorNetwork::ConvDiffParameters()
{
	communicator_->InitializeArray(ConvDiffSize);
	
	{
 	    for(int k = 1; k <= NR_P[procrank_]; k++)
            {
		int index = k + offset[procrank_];
                ConvDiffSize[index] = reactors_[k].iConvectionDiffusion().Size();
            }
        }

	communicator_->GatherArray(ConvDiffSize, NR_P, offset);

	communicator_->InitializeVector(iConvDiff_, ConvDiffSize);
	communicator_->InitializeVector(mConvDiff_, ConvDiffSize);

        {
            for(int k = 1; k <= NR_P[procrank_]; k++)
            {
		int index = k + offset[procrank_];
                for(int j = 1; j <= ConvDiffSize[index]; j++)
                {
                    iConvDiff_[index][j] = reactors_[k].iConvectionDiffusion()[j];
                    mConvDiff_[index][j] = reactors_[k].mConvectionDiffusion()[j];
                }
            }
        }

	communicator_->GatherVector(iConvDiff_, ConvDiffSize);
	communicator_->GatherVector(mConvDiff_, ConvDiffSize);
}

void OpenSMOKE_KPP_ReactorNetwork::SendAndBroadcastOmega()
{	

	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
            {
		int index = k + offset[procrank_];
                for(int j = 1; j <= mix_[0].NumberOfSpecies(); j++)
                {
                    OmegaGlob_[index][j] = reactors_[k].omega()[j];
                }
            }
	}

	communicator_->GatherVector(OmegaGlob_, mix_[0].NumberOfSpecies());
}

void OpenSMOKE_KPP_ReactorNetwork::WriteReactorNetworkDataArray()
{
	communicator_->InitializeArray(React_M);
	communicator_->InitializeArray(React_fInTot);
	communicator_->InitializeArray(React_fOutTot);
	communicator_->InitializeArray(React_mass);
	communicator_->InitializeArray(React_MassFlowIn);
	communicator_->InitializeArray(React_MassUmbalance);	

	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
            {
		int index = k + offset[procrank_];
                React_M[index] = reactors_[k].M();
                React_fInTot[index] = reactors_[k].fInTot();
                React_fOutTot[index] = reactors_[k].fOutTot();
                React_mass[index] = reactors_[k].mass();
                React_MassFlowIn[index] = reactors_[k].MassFlowIn();
                React_MassUmbalance[index] = reactors_[k].MassUmbalance();
            }
	}

	communicator_->GatherArray(React_M, NR_P, offset);
	communicator_->GatherArray(React_fInTot, NR_P, offset);
	communicator_->GatherArray(React_fOutTot, NR_P, offset);
	communicator_->GatherArray(React_mass, NR_P, offset);
	communicator_->GatherArray(React_MassFlowIn, NR_P, offset);
	communicator_->GatherArray(React_MassUmbalance, NR_P, offset);
}

void OpenSMOKE_KPP_ReactorNetwork::InitializeCSTRVectors()
{
    if(procrank_ == 0)
    {
        ChangeDimensions(nprocs_, &Process_iDirectConvergence);
        ChangeDimensions(nprocs_, &Process_iNewtonConvergence);
        ChangeDimensions(nprocs_, &Process_iOdeFirstConvergence);
        ChangeDimensions(nprocs_, &Process_iOdeSecondConvergence);
        ChangeDimensions(nprocs_, &Process_iOdeThirdConvergence);
        ChangeDimensions(nprocs_, &Process_iOdeFourthConvergence);
        ChangeDimensions(nprocs_, &Process_nJacobianEvaluations);
        ChangeDimensions(nprocs_, &Process_nNewtonIterations);
        ChangeDimensions(nprocs_, &Process_nFailures);
        ChangeDimensions(nprocs_, &Process_F1MaxNewton);
        ChangeDimensions(nprocs_, &Process_F1MaxODE);
        ChangeDimensions(nprocs_, &Process_F1MaxUnconverged);
        ChangeDimensions(nprocs_, &Process_normInf);
        ChangeDimensions(nprocs_, &Process_F1Mean);
    }

}

void OpenSMOKE_KPP_ReactorNetwork::SendCSTRInfoToMaster()
{
	MPI::Status status;

    {
        mtype = FROM_WORKER;
        MPI::COMM_WORLD.Send(&iDirectConvergence, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&iNewtonConvergence, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&iOdeFirstConvergence, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&iOdeSecondConvergence, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&iOdeThirdConvergence, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&iOdeFourthConvergence, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&nJacobianEvaluations, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&nNewtonIterations, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&nFailures, 1, MPI::INT, MASTER, mtype);
        MPI::COMM_WORLD.Send(&F1MaxNewton, 1, MPI::DOUBLE, MASTER, mtype);
        MPI::COMM_WORLD.Send(&F1MaxODE, 1, MPI::DOUBLE, MASTER, mtype);
        MPI::COMM_WORLD.Send(&F1MaxUnconverged, 1, MPI::DOUBLE, MASTER, mtype);
        MPI::COMM_WORLD.Send(&normInf, 1, MPI::DOUBLE, MASTER, mtype);
        MPI::COMM_WORLD.Send(&F1Mean, 1, MPI::DOUBLE, MASTER, mtype);
    }

    if(procrank_ == 0)
    {
        mtype = FROM_WORKER;
        for(int i = 0; i <= numworkers_; i++)
        {
            source = i;
            MPI::COMM_WORLD.Recv(&iDirectConvergence, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&iNewtonConvergence, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&iOdeFirstConvergence, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&iOdeSecondConvergence, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&iOdeThirdConvergence, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&iOdeFourthConvergence, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&nJacobianEvaluations, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&nNewtonIterations, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&nFailures, 1, MPI::INT, source, mtype, status);
            MPI::COMM_WORLD.Recv(&F1MaxNewton, 1, MPI::DOUBLE, source, mtype, status);
            MPI::COMM_WORLD.Recv(&F1MaxODE, 1, MPI::DOUBLE, source, mtype, status);
            MPI::COMM_WORLD.Recv(&F1MaxUnconverged, 1, MPI::DOUBLE, source, mtype, status);
            MPI::COMM_WORLD.Recv(&normInf, 1, MPI::DOUBLE, source, mtype, status);
            MPI::COMM_WORLD.Recv(&F1Mean, 1, MPI::DOUBLE, source, mtype, status);

            Process_iDirectConvergence[i+1] = iDirectConvergence;
            Process_iNewtonConvergence[i+1] = iNewtonConvergence;
            Process_iOdeFirstConvergence[i+1] = iOdeFirstConvergence;
            Process_iOdeSecondConvergence[i+1] = iOdeSecondConvergence;
            Process_iOdeThirdConvergence[i+1] = iOdeThirdConvergence;
            Process_iOdeFourthConvergence[i+1] = iOdeFourthConvergence;
            Process_nJacobianEvaluations[i+1] = nJacobianEvaluations;
            Process_nNewtonIterations[i+1] = nNewtonIterations;
            Process_nFailures[i+1] = nFailures;
            Process_F1MaxNewton[i+1] = F1MaxNewton;
            Process_F1MaxODE[i+1] = F1MaxODE;
            Process_F1MaxUnconverged[i+1] = F1MaxUnconverged;
            Process_normInf[i+1] = normInf;
            Process_F1Mean[i+1] = F1Mean;
        }

        iDirectConvergence = Process_iDirectConvergence.GetSumElements();
        iNewtonConvergence = Process_iNewtonConvergence.GetSumElements();
        iOdeFirstConvergence = Process_iOdeFirstConvergence.GetSumElements();
        iOdeSecondConvergence = Process_iOdeSecondConvergence.GetSumElements();
        iOdeThirdConvergence = Process_iOdeThirdConvergence.GetSumElements();
        iOdeFourthConvergence = Process_iOdeFourthConvergence.GetSumElements();
        nJacobianEvaluations = Process_nJacobianEvaluations.GetSumElements();
        nNewtonIterations = Process_nNewtonIterations.GetSumElements();
        nFailures = Process_nFailures.GetSumElements();
        F1MaxNewton = Process_F1MaxNewton.Max();
        F1MaxODE = Process_F1MaxODE.Max();
        F1MaxUnconverged = Process_F1MaxUnconverged.Max();
        normInf = Process_normInf.Max();
        F1Mean = Process_F1Mean.GetSumElements() / double(NR_);
    }

        MPI::COMM_WORLD.Bcast(&iDirectConvergence, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&iNewtonConvergence, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&iOdeFirstConvergence, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&iOdeSecondConvergence, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&iOdeThirdConvergence, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&iOdeFourthConvergence, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&nJacobianEvaluations, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&nNewtonIterations, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&nFailures, 1, MPI::INT, MASTER);
        MPI::COMM_WORLD.Bcast(&F1MaxNewton, 1, MPI::DOUBLE, MASTER);
        MPI::COMM_WORLD.Bcast(&F1MaxODE, 1, MPI::DOUBLE, MASTER);
        MPI::COMM_WORLD.Bcast(&F1MaxUnconverged, 1, MPI::DOUBLE, MASTER);
        MPI::COMM_WORLD.Bcast(&normInf, 1, MPI::DOUBLE, MASTER);
        MPI::COMM_WORLD.Bcast(&F1Mean, 1, MPI::DOUBLE, MASTER);

}

void OpenSMOKE_KPP_ReactorNetwork::GetSingleReactorSparsityPattern()
{
	communicator_->InitializeArray(SparsityPattern);
	
	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
	    {
		int index = k + offset[procrank_];
	        SparsityPattern[index] = reactors_[k].nGlobalSingleRowSparsityPattern();
	    }
	}
	
	communicator_->GatherArray(SparsityPattern, NR_P, offset);
}

void OpenSMOKE_KPP_ReactorNetwork::InitializeNetworkData()
{

	communicator_->InitializeArray(globalIndicesSparsityValuesSize);
	communicator_->InitializeArray(timeStepsArray);
	communicator_->InitializeArray(gISVSizeCumulated, NR_);
        communicator_->InitializeArray(Hmix_);
        communicator_->InitializeArray(Hout_);
        communicator_->InitializeArray(Hin_);
        communicator_->InitializeArray(delta_T_);
        communicator_->InitializeArray(Delta_H_start);
        communicator_->InitializeArray(Delta_H_end);

	communicator_->InitializeVector(OmegaGlob_, mix_[0].NumberOfSpecies());
	nOmega_ = NR_P[procrank_] * mix_[0].NumberOfSpecies();
	communicator_->InitializePetscVector(OmegaFractions_, nOmega_, highindex_omega, lowindex_omega);
	
	omega_ptr = new double[NumberOfEquations() + 1];
	x_ptr = new double[NumberOfEquations() + 1];
//	direction_ptr = new double[NumberOfEquations() + 1];
	localRHS_ = new double[NR_P[procrank_] * mix_[0].NumberOfSpecies()];
}

void OpenSMOKE_KPP_ReactorNetwork::UpdateRHS()
{

	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
	    {
		int position = (k-1) * NumberOfSpecies();
		for(int i = 0; i < NumberOfSpecies(); i++)
		    localRHS_[position + i] = reactors_[k].LocalRHS()[i+1]+reactors_[k].NonLocalRHS()[i+1];	    
	    }
	}

/*	if(procrank_ > 0)
	{
	    double *value = new double[NR_P[procrank_] * NumberOfSpecies()];
	    for(int k = 1; k <= NR_P[procrank_]; k++)
	    {
		int index = procrank_ + numworkers_ * (k-1);
		int position = (k-1) * NumberOfSpecies();
		for(int i = 0; i < NumberOfSpecies(); i++)
		    value[position + i] = reactors_[k].LocalRHS()[i+1]+reactors_[k].NonLocalRHS()[i+1];
	    
	    }
	    RHS_place = RHS_place_get_array;

	    VecSetValues(RHSget_, NR_P[procrank_] * NumberOfSpecies(), RHS_place_get_array, value, INSERT_VALUES);
	    delete [] value;
	}

	VecAssemblyBegin(RHSget_);
	VecAssemblyEnd(RHSget_);
	VecGetArray(RHSget_, &RHS_value);


	{
	    if(procrank_ > 0)	VecSetValues(RHSset_, nRHSget_, RHS_place_set, RHS_value, INSERT_VALUES);
	}

	VecAssemblyBegin(RHSset_);
	VecRestoreArray(RHSget_, &RHS_value);
	VecAssemblyEnd(RHSset_);*/	
}

void OpenSMOKE_KPP_ReactorNetwork::UpdateLISMatrix()
{	

	int counter = 0;
	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    for(int i = 0; i < reactors_[k].globalIndicesSparsityValues().Size(); i++)
	    {
		Petsc_values[counter] = reactors_[k].globalIndicesSparsityValues()[i+1];
		counter++;
	    }
	}

	lisMatrixGlobal->UpdateMatrix(Petsc_values);


/*	VecSetValues(PetscMatrix_, matrix_place_size, matrix_place_array, Petsc_values, INSERT_VALUES);
	VecAssemblyBegin(PetscMatrix_);
	VecAssemblyEnd(PetscMatrix_);

	lisMatrixGlobal->UpdateMatrix(PetscMatrix_);*/

/*
	delete [] value;

//	double starttime;

//	if(procrank_ == 0) starttime = MPI::Wtime();

	if(procrank_ > 0)
	{
	    double *value = new double[countLocal_Old[procrank_]];
	    int counter_RHS = 0;
	    for(int k = 1; k <= NR_P[procrank_]; k++)
	    {
		for(int i = 0; i < IndicesSparsityValuesSize[k]; i++)
		{
		    value[counter_RHS] = reactors_[k].globalIndicesSparsityValues()[i+1];
		    counter_RHS++;
		}
	    }

	    place = matrix_place_array;
		
	    VecSetValues(MatrixGet_, countLocal_Old[procrank_], place, value, INSERT_VALUES);
	    delete [] value;
	}

//	if(procrank_ == 0)	cout << "Fin qui arrivo...0bis" << std::endl;

	VecAssemblyBegin(MatrixGet_);
	VecGetArray(MatrixGet_, &value);
	VecAssemblyEnd(MatrixGet_);

//	MPI::COMM_WORLD.Barrier();
//	if(procrank_ == 0)	cout << "Time for building 1st vector = " << MPI::Wtime() - starttime << std::endl;
//	if(procrank_ == 0)	starttime = MPI::Wtime();
//	if(procrank_ == 0)	cout << "Fin qui arrivo...A" << std::endl;


	{
	    if(procrank_ > 0)    VecSetValues(MatrixSet_, nGet_, place_set, value, INSERT_VALUES);
	}

//	if(procrank_ == 0)	cout << "Fin qui arrivo...B" << std::endl;

	VecAssemblyBegin(MatrixSet_);
	VecRestoreArray(MatrixGet_, &value);
	VecAssemblyEnd(MatrixSet_);

//	MPI::COMM_WORLD.Barrier();
//	if(procrank_ == 0)	cout << "Assembled" << std::endl;

//	MPI::COMM_WORLD.Barrier();
//	if(procrank_ == 0)	cout << "Time for building 2nd vector = " << MPI::Wtime() - starttime << std::endl;
///	if(procrank_ == 0)	starttime = MPI::Wtime();

	lisMatrixGlobal->UpdateMatrix(MatrixSet_);

//	if(procrank_ == 0)	cout << "Fin qui arrivo...C" << std::endl;

//	MPI::COMM_WORLD.Barrier();
//	if(procrank_ == 0)	cout << "Time for building matrix = " << MPI::Wtime() - starttime << std::endl;

//	getchar();*/
}

void OpenSMOKE_KPP_ReactorNetwork::CommunicationPattern()
{
	MPI::Status status;

	for(int p = 1; p <= numworkers_; p++)
	{
	    mtype = FROM_MASTER;
	    int* process_comm = new int[NR_ + 1];
	    if(procrank_ == 0)
	    {
		for(int k = 1; k <= NR_; k++)
		{
		    if(React_comm_master[p][k] == true)
			process_comm[k] = 1;
		    else
			process_comm[k] = 0;
		}

		MPI::COMM_WORLD.Send(&process_comm[1], NR_, MPI::INT, p, mtype);
	    }
	
	    if(procrank_ == p)
	    {
		MPI::COMM_WORLD.Recv(&process_comm[1], NR_, MPI::INT, MASTER, mtype, status);

		for(int k = 1; k <= NR_; k++)
		{
		    if(process_comm[k] == 1)
			React_communication[k] = true;
		    else
			React_communication[k] = false;
		}
	    }

	    delete [] process_comm;
	}
}

void OpenSMOKE_KPP_ReactorNetwork::InitializePetscMatrix()
{

	nMatrix_ = lisMatrixGlobal->countLocal()[procrank_];

	communicator_->InitializePetscVector(PetscMatrix_, nMatrix_, upper_matrix, lower_matrix);

	matrix_place.resize(NR_P[procrank_] + 1);

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int size = reactors_[k].globalIndicesSparsityValues().Size();
	    matrix_place[k] = new long long int[size];
	}

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = procrank_ + 1 + nprocs_ * (k-1);
	    for(int i = 0; i < reactors_[k].globalIndicesSparsityValues().Size(); i++)
	        matrix_place[k][i] = gISVSizeCumulated[index - 1] + i;
	}

	matrix_place_size = 0;

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    matrix_place_size += reactors_[k].globalIndicesSparsityValues().Size();
	}

	matrix_place_array = new long long int[matrix_place_size];

	long long int counter_matrix = 0;
	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
	    {
	        for(int i = 0; i < reactors_[k].globalIndicesSparsityValues().Size(); i++)
	        {
		    matrix_place_array[counter_matrix] = matrix_place[k][i];
		    if(matrix_place_array[counter_matrix] < 0) std::cout << matrix_place_array[counter_matrix] << std::endl;
		    counter_matrix++;
	        }
	    }
	}

	Petsc_values = new double[matrix_place_size];

	for(int k = 1; k <= NR_P[procrank_]; k++)
	    delete [] matrix_place[k];

/*	VecCreate(PETSC_COMM_WORLD, &MatrixGet_);
	VecCreate(PETSC_COMM_WORLD, &MatrixSet_);
		
	nSet_ = lisMatrixGlobal->countLocal()[procrank_];
	VecSetSizes(MatrixSet_, nSet_, PETSC_DECIDE);
    	VecSetFromOptions(MatrixSet_);
    	VecGetSize(MatrixSet_,&nSetAll_);
	
//	lisMatrixGlobal->MatrixDistribution();

	nGet_ = countLocal_Old[procrank_];
	VecSetSizes(MatrixGet_, nGet_, PETSC_DECIDE);
	VecSetFromOptions(MatrixGet_);
	VecGetSize(MatrixGet_,&nGetAll_);

	VecGetOwnershipRange(MatrixGet_, &lowindex_get, &highindex_get);
	VecGetOwnershipRange(MatrixSet_, &lowindex_set, &highindex_set);

    	place_set = new long long int[nGet_];
    	for(int i = 0; i < nGet_; i++)
    	{
	    place_set[i] = i + lowindex_get;
    	}

	place_get.resize(NR_P[procrank_] + 1);

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int size = reactors_[k].globalIndicesSparsityValues().Size();
	    place_get[k] = new long long int[size];
	}

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = procrank_ + numworkers_ * (k-1);
	    for(int i = 0; i < reactors_[k].globalIndicesSparsityValues().Size(); i++)
	    place_get[k][i] = gISVSizeCumulated[index - 1] + i;
	}

	if(procrank_ > 0) matrix_place_array = new long long int[countLocal_Old[procrank_]];

	long long int counter_matrix = 0;
	if(procrank_ > 0)
	{
	    for(int k = 1; k <= NR_P[procrank_]; k++)
	    {
	        for(int i = 0; i < reactors_[k].globalIndicesSparsityValues().Size(); i++)
	        {
		    matrix_place_array[counter_matrix] = place_get[k][i];
		    if(matrix_place_array[counter_matrix] < 0) std::cout << matrix_place_array[counter_matrix] << std::endl;
		    counter_matrix++;
	        }
	    }
	}

	for(int k = 1; k <= NR_P[procrank_]; k++)
	    delete [] place_get[k];
*/
}

void OpenSMOKE_KPP_ReactorNetwork::InitializePetscRHS()
{
	nRHS_ = lisMatrixGlobal->local_rows()[procrank_];
	communicator_->InitializePetscVector(PetscRHS_, nRHS_, RHS_lowindex, RHS_highindex);


	RHS_place.resize(NR_P[procrank_] + 1);

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int size = NumberOfSpecies();
	    RHS_place[k] = new long long int[size];
	}


	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = procrank_ + 1 + nprocs_ * (k-1);
	    for(int i = 0; i < NumberOfSpecies(); i++)
	        RHS_place[k][i] = (index - 1) * NumberOfSpecies() + i;
	}

	RHS_place_array = new long long int[NR_P[procrank_] * NumberOfSpecies()];
	int counter_RHS = 0;
	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    for(int i = 0; i < NumberOfSpecies(); i++)
	    {
		RHS_place_array[counter_RHS] = RHS_place[k][i];
		counter_RHS++;
	    }
	}

	for(int k = 1; k <= NR_P[procrank_]; k++)
	    delete [] RHS_place[k];
/*
	RHS_place_array = new long long int[NR_P[procrank_] * NumberOfSpecies()];
	int counter_RHS = 0;
	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    for(int i = 0; i < NumberOfSpecies(); i++)
	    {
		    RHS_place_array[counter_RHS] = RHS_place[k][i];
		    counter_RHS++;
	    }
	}


	for(int k = 1; k <= NR_P[procrank_]; k++)
	    delete [] RHS_place[k];*/

/*	VecCreate(PETSC_COMM_WORLD, &RHSget_);
	nRHSget_ = NR_P[procrank_] * NumberOfSpecies();
	VecSetSizes(RHSget_, nRHSget_, PETSC_DECIDE);
	VecSetFromOptions(RHSget_);


	VecCreate(PETSC_COMM_WORLD, &RHSset_);
	nRHSset_ = lisMatrixGlobal->local_rows()[procrank_];
	VecSetSizes(RHSset_, nRHSset_, PETSC_DECIDE);	
	VecSetFromOptions(RHSset_);


	VecGetOwnershipRange(RHSget_, &RHS_lowindex_get, &RHS_highindex_get);
	VecGetOwnershipRange(RHSset_, &RHS_lowindex_set, &RHS_highindex_set);

	RHS_place_get.resize(NR_P[procrank_] + 1);

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int size = NumberOfSpecies();
	    RHS_place_get[k] = new long long int[size];
	}


	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = procrank_ + numworkers_ * (k-1);
	    for(int i = 0; i < NumberOfSpecies(); i++)
	        RHS_place_get[k][i] = (index - 1) * NumberOfSpecies() + i;
	}

	RHS_place_get_array = new long long int[NR_P[procrank_] * NumberOfSpecies()];
	int counter_RHS = 0;
	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    for(int i = 0; i < NumberOfSpecies(); i++)
	    {
		RHS_place_get_array[counter_RHS] = RHS_place_get[k][i];
		counter_RHS++;
	    }
	}

	for(int k = 1; k <= NR_P[procrank_]; k++)
	    delete [] RHS_place_get[k];

    	RHS_place_set = new long long int[nRHSget_];
    	for(int i = 0; i < nRHSget_; i++)
    	{
	    RHS_place_set[i] = i + RHS_lowindex_get;
    	}*/
}

void OpenSMOKE_KPP_ReactorNetwork::CleanMemory()
{
	DeleteTopologyParametersSize();
}

void OpenSMOKE_KPP_ReactorNetwork::InitializePlaceArray()
{
	place_vector.resize(NR_P[procrank_] + 1);

	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int size = NumberOfSpecies();
	    place_vector[k] = new long long int[size];
	}


	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    int index = k + offset[procrank_];
	    for(int i = 0; i < NumberOfSpecies(); i++)
	        place_vector[k][i] = (index - 1) * NumberOfSpecies() + i;
	}

	place_array = new Petsc64bitInt[NR_P[procrank_] * NumberOfSpecies()];
	int counter_place = 0;
	for(int k = 1; k <= NR_P[procrank_]; k++)
	{
	    for(int i = 0; i < NumberOfSpecies(); i++)
	    {
		place_array[counter_place] = place_vector[k][i];
		counter_place++;
	    }
	}

	for(int k = 1; k <= NR_P[procrank_]; k++)
	    delete [] place_vector[k];
}

void OpenSMOKE_KPP_ReactorNetwork::CheckEnergyBalances()
{
    

    for(int k = 1; k <= NR_P[procrank_]; k++)
    {
        reactors_[k].EnergyOutflow(&mix_[0]);
        
        Hmix_[k + offset[procrank_]] = reactors_[k].H_mix();  
        Hout_[k + offset[procrank_]] = reactors_[k].H_out();
    }
    
    communicator_->GatherArray(Hmix_, NR_P, offset);
    communicator_->GatherArray(Hout_, NR_P, offset);

    MPI::COMM_WORLD.Barrier();
    for(int k = 1; k <= NR_P[procrank_]; k++)
    {
        reactors_[k].EnergyInflow(&mix_[0], *this);
        Hin_[k + offset[procrank_]] = reactors_[k].H_in();
    }

    communicator_->GatherArray(Hin_, NR_P, offset);
    
/*    for (int k = 1; k <= NR_P[procrank_]; k++)
        {
            if (k + offset[procrank_] == 210)
            {
                std::ofstream enthalpy("Reattore210.out", ios::out);
                enthalpy << reactors_[k].M() << std::setw(15) << Hmix_[k + offset[procrank_]] << std::endl;
                for(int i = 1; i <= mix_[0].NumberOfSpecies(); i++)
                    enthalpy << mix_[0].names[i] << "   " << reactors_[k].omega()[i] << std::endl;
            
                enthalpy << std::endl;
                BzzVector cell_h(mix_[0].NumberOfSpecies());
                mix_[0].SpeciesEnthalpy(reactors_[k].temperature());
                Product((Constants::R_J_kmol * reactors_[k].temperature()), mix_[0].h, &cell_h); // [J/kmol]
                ElementByElementProduct(cell_h, mix_[0].uM, &cell_h); // [J/kg]
            
                for(int i = 1; i <= mix_[0].NumberOfSpecies(); i++)
                    enthalpy << cell_h[i] << std::endl;
            }
        }*/

    if (procrank_ == 0)
    {
        std::ofstream fEnergyMap;
        if (data_.networkStatus() == KPP_NETWORK_STATUS_START)
        {
            openOutputFileAndControl(fEnergyMap, data_.nameFolderOutput() + "/Additional/Energy/Enthalpy_start.out");
            fEnergyMap.setf(std::ios::scientific);

            fEnergyMap << std::setw(15) << "H_in [J/s]" << std::setw(15) << "H_out [J/s]" << std::setw(15)
                    << "Delta_H [J/s]" << std::setw(15) << "% error" << std::endl;
            
            for(int i = 1; i <= NR_; i++)
            {
                Delta_H_start[i] = Hout_[i] - Hin_[i];
                fEnergyMap << std::setw(15) << Hin_[i] << std::setw(15) <<  Hout_[i] << std::setw(15)
                        << Hout_[i] - Hin_[i] << std::setw(15) << fabs((Hout_[i] - Hin_[i])/Hout_[i]) * 100 << std::endl;
            }
            fEnergyMap.close();
            
            BzzVector EnergyAbs(NR_);
            BzzVector EnergyRel(NR_);

            for (int i = 1; i <= NR_; i++)
            {
                EnergyAbs[i] = Hout_[i] - Hin_[i];
                EnergyRel[i] = fabs((Hout_[i] - Hin_[i])/Hout_[i]) * 100;
            }

            std::cout << std::endl;
            std::cout << "Max Absolute Energy Umbalance = " << std::fixed << EnergyAbs.Max() << " J/s" << std::endl;
//            std::cout << "Max Relative Energy Umbalance = " << std::fixed << EnergyRel.Max() << " %" << std::endl;
            
            if(data_.iTraditionalLink() == false)
            {
                std::string deltaH_path = "Output/Additional/Energy/deltaHstart.bkp";
                
                BzzVector dH_start(NR_);
                
                for(int i = 1; i <= NR_; i++)
                    dH_start[i] = Delta_H_start[i];
                
                BzzSave fDeltaH_Backup('*', deltaH_path.c_str());
                
                fDeltaH_Backup << dH_start;
                
                fDeltaH_Backup.End();
            }
        }

        else
        {
            openOutputFileAndControl(fEnergyMap, data_.nameFolderOutput() + "/Additional/Energy/Enthalpy_end.out");
            fEnergyMap.setf(std::ios::scientific);

            fEnergyMap << std::setw(15) << "H_in [J/s]" << std::setw(15) << "H_out [J/s]" << std::setw(15)
                    << "Delta_H [J/s]" << std::setw(15) << "% error" << std::endl;

            for (int i = 1; i <= NR_; i++)
            {
                Delta_H_end[i] = Hout_[i] - Hin_[i];
                fEnergyMap << std::setw(15) << Hin_[i] << std::setw(15) << Hout_[i] << std::setw(15)
                        << Hout_[i] - Hin_[i] << std::setw(15) << fabs((Hout_[i] - Hin_[i]) / Hout_[i]) << std::endl;
            }
            fEnergyMap.close();
            
            BzzVector EnergyAbs(NR_);
            BzzVector EnergyRel(NR_);

            for (int i = 1; i <= NR_; i++)
            {
                EnergyAbs[i] = Hout_[i] - Hin_[i];
                
                if(fabs(EnergyAbs[i]) > 1e-6)
                    EnergyRel[i] = fabs((Hout_[i] - Hin_[i]) / Hout_[i]) * 100;
                
                else
                    EnergyRel[i] = 0;
            }

            std::cout << std::endl;
            std::cout << "Max Absolute Energy Umbalance = " << std::fixed <<  EnergyAbs.Max() << " J/s" << std::endl;
//            std::cout << "Max Relative Energy Umbalance = " << std::fixed <<  EnergyRel.Max() << " %" << std::endl;
            
            
            
            BzzVector EnergyDeltaAbs(NR_);
            BzzVector EnergyDeltaRel(NR_);
            
            if(data_.iEnergyClustering() == true)
            {
                std::string deltaH_path = "Output/Additional/Energy/deltaHstart.bkp";
                BzzVector dH_start(NR_);
                
                BzzLoad fDeltaH_Backup('*', deltaH_path.c_str());
                fDeltaH_Backup >> dH_start;
                
                fDeltaH_Backup.End();
                
                for(int i = 1; i <= NR_; i++)
                    Delta_H_start[i] = dH_start[i];
            }
            
            for(int i = 1; i <= NR_; i++)
            {
                EnergyDeltaAbs[i] = fabs(Delta_H_end[i] - Delta_H_start[i]);

                if (EnergyDeltaAbs[i] > 1e-6)
                    EnergyDeltaRel[i] = fabs((Delta_H_end[i] - Delta_H_start[i]) / Delta_H_start[i]) * 100;

                else
                    EnergyDeltaRel[i] = 0;
            }
            
            std::cout << "Max Absolute Delta (Start/End) Energy Umbalance = " << EnergyDeltaAbs.MaxAbs() << " J/s" << std::endl;
//            std::cout << "Max Relative Delta (Start/End) Energy Umbalance = " << EnergyDeltaRel.Max() << " %" << std::endl;
            
        }
        
    }
    
    communicator_->BroadcastArray(Delta_H_start, NR_);
    communicator_->BroadcastArray(Delta_H_end, NR_);

    for(int k = 1; k <= NR_P[procrank_]; k++)
    {
        int index = k + offset[procrank_];
        delta_T_[index] = (Delta_H_end[index] - Delta_H_start[index]) / (reactors_[k].M() * reactors_[k].cp_mix());
    }
    
    communicator_->GatherArray(delta_T_, NR_P, offset);
    
    if(procrank_ == 0)
    {
        std::ofstream fEnergyUmbalancesMap;
        openOutputFileAndControl(fEnergyUmbalancesMap, data_.nameFolderOutput() + "/Additional/Energy/Temperature_Umbalance.out");
        fEnergyUmbalancesMap.setf(std::ios::scientific);

        for (int i = 1; i <= NR_; i++) //Aggiustare questa funzione
        {
            fEnergyUmbalancesMap << delta_T_[i] << std::endl;
        }

        fEnergyUmbalancesMap.close();
        
        BzzVector deltaT(NR_);
        for(int i = 1; i <= NR_; i++)
        {
            deltaT[i] = delta_T_[i];
        }
        
        std::cout << "Max Absolute Delta T Umbalance = "  << deltaT.MaxAbs() << " K" << std::endl;
    }
    
    
/*    if(procrank_ == 0)
    {
        std::ofstream borders("Boundaries.out", ios::out);
        
        borders << "Feed" << std::endl;
        for(int i = 1; i <= jExternalFeedReactors.Size(); i++)
            borders << jExternalFeedReactors[i] << std::endl;
        borders << std::endl << "Outputs" << std::endl;
        
        for(int i = 1; i <= jExternalOutputReactors.Size(); i++)
            borders << jExternalOutputReactors[i] << std::endl;
    }
    
    jExternalFeedReactors.BzzPrint();
    jExternalOutputReactors.BzzPrint();*/
    
    
    /*
     	// External feeds
	if (tagExternalFeed_ == true)	mInTot_ = fIn_;
	else							mInTot_ = 0.;

	// Inflow
	for(int j=1;j<=iConvectionDiffusion_.Size();j++)
	{
	    for(int i = 1; i <= numberOfSpecies; i++)
	    	mInTot_[i] += mConvectionDiffusion_[j] * network.OmegaGlob()[iConvectionDiffusion_[j]][i];
	}
     
     */
}
