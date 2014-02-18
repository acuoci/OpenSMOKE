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
#include <sstream>
#include "linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.h"
#include "linear_solvers/OpenSMOKE_MUMPS_Unsymmetric.h"
#include "linear_solvers/OpenSMOKE_LIS_Unsymmetric.h"
#include "OpenSMOKE_KPP_BlockMatrixNetwork.h"
#include "OpenSMOKE_KPP_ReactorNetwork.h"
#include "OpenSMOKE_KPP_DataManager.h"
#include "OpenSMOKE_KPP_ConvectiveNetworkStatistics.h"
#include "OpenSMOKE_KPP_SingleReactorStatistics.h"
#include "OpenSMOKE_KPP_ReactorNetwork_Residuals.h"
#include "OpenSMOKE_KPP_NewtonMethod_Manager.h"
#include "OpenSMOKE_KPP_ODE_Manager.h"
#include "OpenSMOKE_KPP_SingleReactor_KineticsManager.h"

OpenSMOKE_KPP_ReactorNetwork *ptNetwork;

void externalResiduals(BzzVector &x, BzzVector &f)
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

void OpenSMOKE_KPP_ReactorNetwork::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_ReactorNetwork"	<< endl;
    cout << "Error: " << message_				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_ReactorNetwork::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_ReactorNetwork"		<< endl;
    cout << "Warning: " << message_					<< endl;
	cout << endl;
}

OpenSMOKE_KPP_ReactorNetwork::OpenSMOKE_KPP_ReactorNetwork(OpenSMOKE_ReactingGas& mix, OpenSMOKE_KPP_DataManager& data, ofstream& fLog, ofstream& fWarning) :
mix_(mix), data_(data), fLog_(fLog), fWarning_(fWarning)
{
	ptNetwork	= this;

	iteration_   = 0;
	globalTime_  = 0.;
	deltat		 = 0.;

	iAlreadyGlobalPardiso_ = false;
	iAlreadyGlobalMumps_ = false;
	iAlreadyGlobalLis_ = false;
	iAlreadyGlobalGaussSiedel_ = false;

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
        
	currentDeltaTimeGlobalODE_ = data_.GlobalODE_InitialTimeStep()*5.;
}

OpenSMOKE_KPP_ReactorNetwork::~OpenSMOKE_KPP_ReactorNetwork(void)
{
}

void OpenSMOKE_KPP_ReactorNetwork::PreparingIndicesReactorsThreads()
{
	// Preparing
	indicesReactorsThreads_ = new BzzVectorInt[data_.nThreads()];
/*	
	// Policy 1
	for (int k=0;k<data_.nThreads();k++)
	{
		int nBlock = NR_/data_.nThreads();
		int kStart = nBlock*k+1;
		int kEnd   = nBlock*(k+1);
		if (k == data_.nThreads()-1)
			kEnd = NR_;
		
		ChangeDimensions(kEnd-kStart+1, &indicesReactorsThreads_[k]);
		int i=1;
		for (int j=kStart;j<=kEnd;j++)
			indicesReactorsThreads_[k][i++] = j;
	}
   */     
	// Policy 1
        int nBlock = NR_/data_.nThreads();
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

        while(reactor<=NR_)
        {
            for (int k=0;k<data_.nThreads();k++)
            {
                indicesReactorsThreads_[k].Append(reactor++);
                if (reactor>NR_) break;
            }
        }
            
        // Check
        int sum = 0;
        for (int k=0;k<data_.nThreads();k++)
            sum += indicesReactorsThreads_[k].Size();
        if (sum != NR_)
            ErrorMessage("Wrong divide and conquer policy!");
}

void OpenSMOKE_KPP_ReactorNetwork::MemoryAllocation()
{
	ChangeDimensions(NR_, NR_, &C_);
	ChangeDimensions(NR_, NR_, &A_);
	RHS_ = new BzzVector[mix_.NumberOfSpecies()+1];
	for (int j=1;j<=mix_.NumberOfSpecies();j++)
		ChangeDimensions(NR_, &RHS_[j]);
	bStar_ = new BzzVector[mix_.NumberOfSpecies()+1];
	for (int j=1;j<=mix_.NumberOfSpecies();j++)
		ChangeDimensions(NR_, &bStar_[j]);

	ChangeDimensions(NR_*mix_.NumberOfSpecies(), &lastGlobalODESolution_);
	ChangeDimensions(NR_*mix_.NumberOfSpecies(), &lastGlobalNLSSolution_);
	ChangeDimensions(NR_*mix_.NumberOfSpecies(), &auxVector_NRxNC_1);

	residuals_ = new OpenSMOKE_KPP_ReactorNetwork_Residuals(*this);
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
	ifstream fInput(data_.nameFirstGuessFile().c_str(), ios::in);

	unsigned int NCDetailed;
	unsigned int NCReduced;

	fInput >> NCDetailed;
	fInput >> NCReduced;

	ChangeDimensions(NCReduced, &jReduced);
	vector<string> namesReduced;
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
                jReduced[j] = mix_.recognize_species(namesReduced[j]);
        }
        
        for (int j=1;j<=jReduced.Size();j++)
            cout << jReduced[j] << " " << mix_.names[jReduced[j]] << endl;
        
        
	fInput >> NR_;

	// Allocating kinetics
	if (data_.iSaveKineticConstants() == true)
		kinetics_ = new OpenSMOKE_KPP_SingleReactor_KineticsManager[NR_+1];
	else if (data_.iSaveKineticConstants() == false)
	{
		if (data_.iOpenMP() == false)
			kinetics_ = new OpenSMOKE_KPP_SingleReactor_KineticsManager[1];
		else
			kinetics_ = new OpenSMOKE_KPP_SingleReactor_KineticsManager[data_.nThreads()];
	}
        
        // Indices threads
        PreparingIndicesReactorsThreads();



	// Allocating reactors
	reactors_ = new OpenSMOKE_KPP_SingleReactor[NR_+1];

	for (int k=1;k<=NR_;k++)
	{
		reactors_[k].Setup(k, &mix_, jReduced);
		reactors_[k].SetDataManager(&data_);
	}

	ChangeDimensions(NR_, NCDetailed, &residualMatrix_);

	unsigned int dummy;
	fInput >> dummy;
	if (dummy != NCReduced)
		ErrorMessage("The number of reduced species in the FirstGuess file does not fit!");
	if (mix_.NumberOfSpecies() != NCDetailed)
		ErrorMessage("The number of detailed species in the FirstGuess file does not fit with the number of species in the detailed kinetic file!");

	for (int k=1;k<=NR_;k++)
	{
		if ( k%(NR_/10) == 1)	cout << "Reading mass fractions: " << k << "/" << NR_ << endl;
		reactors_[k].ReadMassFractions(fInput);
	}

	if (data_.iBackup() == true)
		ReadBackupFile(data_.nameInputBackupFile());

	fInput.close();
}

void OpenSMOKE_KPP_ReactorNetwork::ReadTopology()
{
	ifstream fInput(data_.nameTopologyFile().c_str(), ios::in);

	unsigned int dummy;
	fInput >> dummy;
	if (dummy != NR_)
		ErrorMessage("The number of reactors does not fit!");

	for (int k=1;k<=NR_;k++)
	{
		if ( k%(NR_/10) == 1)	cout << "Reading reactor topology: " << k << "/" << NR_ << endl;
		reactors_[k].ReadReactorProperties(fInput);
		reactors_[k].ReadTopology(fInput);	
	}

	// Min and Max
	double Tmin = 10000;
	double Tmax = 0;
	for (int k=1;k<=NR_;k++)	
	{	
		if (reactors_[k].temperature() > Tmax) Tmax = reactors_[k].temperature();
		else if (reactors_[k].temperature() < Tmin) Tmin = reactors_[k].temperature();
	}

	if (data_.TemperatureMin() == 0.)	data_.SetMinimumTemperature(Tmin-20., "K");
	if (data_.TemperatureMax() == 0.)	data_.SetMaximumTemperature(Tmax*1.07, "K");

	std::cout << "Minimum temperature:             " << Tmin << " K " << std::endl;
	std::cout << "Maximum temperature:             " << Tmax << " K " << std::endl;
	std::cout << "Minimum temperature (effective): " << data_.TemperatureMin() << " K " << std::endl;
	std::cout << "Minimum temperature (effective): " << data_.TemperatureMax() << " K " << std::endl;

	if (data_.correction() != KPP_CORRECTION_NONE)
	{
		if (data_.TemperatureMin() >= Tmin-10)
			ErrorMessage("The imposed minimum temperature is too large! At least 10 K below than minimum!");

		if ( data_.TemperatureMax() <= Tmax*1.02)
			ErrorMessage("The imposed minimum temperature is too small! At least 2 percent higher than maximum!");
	}
	
	if (data_.iSaveKineticConstants() == true)
	{
		if (data_.correction() == KPP_CORRECTION_NONE)
		{
			for (int k=1;k<=NR_;k++)
			{
				kinetics_[k].SetMinMax(data_.TemperatureMin(), data_.TemperatureMax());
				kinetics_[k].Setup(data_.correction(), &mix_, reactors_[k].temperature(), reactors_[k].pressure(), -1. );
				reactors_[k].SetKinetics(&kinetics_[k]);
			}
		}
		else 
		{
			for (int k=1;k<=NR_;k++)
			{
				kinetics_[k].SetMinMax(data_.TemperatureMin(), data_.TemperatureMax());
				kinetics_[k].Setup(data_.correction(), &mix_, reactors_[k].temperature(), reactors_[k].pressure(), reactors_[k].variance() );
				reactors_[k].SetKinetics(&kinetics_[k]);
			}
		}
	}
	else if (data_.iSaveKineticConstants() == false)
	{
		ErrorMessage("TO UPDATE SaveKineticConstants()==false (especially Corrections)");
		if (data_.iOpenMP() == false)
		{
			kinetics_[0].Setup(&mix_);
			for (int k=1;k<=NR_;k++)
				reactors_[k].SetKinetics(&kinetics_[0]);
		}
		else
		{
			for (int k=0;k<data_.nThreads();k++)
				kinetics_[k].Setup((&(mix_)+k));
			for (int k=1;k<=NR_;k++)
				reactors_[k].SetKinetics(&kinetics_[0]);
		}
	}

	fInput.close();
}

void OpenSMOKE_KPP_ReactorNetwork::BuildNetwork()
{
	// Memory Allocation
	cout << " * Memory allocation..." << endl;
	MemoryAllocation();
    cout << "Done: Memory Allocation" << endl; //getchar();

	// Build network
	cout << " * Assembling mass flow rate topology..." << endl;
	for (int k=1;k<=NR_;k++)
	{
		for(int j=1;j<=reactors_[k].out().Size();j++)
			reactors_[reactors_[k].out()[j]].UpdateInflow(k, reactors_[k].cOut()[j], reactors_[k].diffusion()[j]);
	}
    cout << "Done: Assembling mass flow rate topology" << endl; // getchar();
        
	// Assemblig local reactors
	cout << " * Assembling local topology..." << endl;
	for (int k=1;k<=NR_;k++)
		reactors_[k].Assembling();
    cout << "Done: Assembling local topology" << endl; //getchar();

	// Mass flow rate Errors
	cout << " * Mass flow errors (before corrections)" << endl;
	MassFlowErrors();
    cout << "Done: Mass flow errors" << endl; //getchar();

	// Adjust mass fluxes
	cout << " * Correcting mass flow rates..." << endl;
	CorrectingMassFlowRates();
    cout << "Done: Correcting mass flow rates" << endl; //getchar();

	// Mass flow rate Errors
	cout << " * Mass flow errors (after corrections)" << endl;
	MassFlowErrors();
    cout << "Done: Mass flow errors" << endl; //getchar();
	
	// External feeds and outputs indices
	for (int k=1;k<=NR_;k++)
	{
		if (reactors_[k].IsExternalFeedReactor() == true)	jExternalFeedReactors.Append(k);
		if (reactors_[k].IsExternalOutputReactor() == true)	jExternalOutputReactors.Append(k);
	}
    cout << "Done: External feeds and outputs indices" << endl; //getchar();

	// Convection-Diffusion Matrix
	cout << " * Assembling Convection-Diffusion Matrix..." << endl;
	AssemblingConvectionDiffusionMatrix();
    cout << "Done: Convection-Diffusion Matrix" << endl; //getchar();

	// External feeds Right Hand Sides
	cout << " * Assembling Right Hand Sides..." << endl;
	AssemblingRightHandSides();
    cout << "Done: Assembling Right Hand Sides" << endl; //getchar();

	// Exchange info with single reactors
	cout << " * Convection-Diffusion matrix communication..." << endl;
	for (int k=1;k<=NR_;k++)
		reactors_[k].ConvectionDiffusionMatrixCommunication(*this);
    cout << "Done: Convection-Diffusion matrix communication" << endl; //getchar();
        
	// Statistics on C matrix
	SummaryConvectionDiffusionMatrix();
    cout << "Done: SummaryConvectionDiffusionMatrix" << endl; //getchar();

	// Write additional info
	WriteReactorNetworkData();
    cout << "Done: WriteReactorNetworkData" << endl; //getchar();

	// Initial residuals
	ResidualsAnalysis();
    cout << "Done: ResidualsAnalysis" << endl; //getchar();

	// Initialize OdeSystem
	cout << " * Initialize reactors..." << endl;
    odePool = new ODE_Pool(data_.nThreads(), NumberOfSpecies());
    for(int k=0;k<data_.nThreads();k++)
		odePool->Set(k, reactors_[k+1]);
	odePool->Initialize();
        
//	for (int k=1;k<=NR_;k++)
//		reactors_[k].SetInitialConditions();
    cout << "Done: Initialize reactors" << endl; //getchar();

	cout << " * Initialize convection sparse matrix..." << endl;
	if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
	{
		// PARDISO Convection Linear Systems 
		pardisoMatrixConvection = new OpenSMOKE_PARDISO_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
		pardisoMatrixConvection->SetFillInReducingOrdering(data_.FillInReducingOrdering());
		pardisoMatrixConvection->SetSparsityPattern(C_);
	}
	else if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
	{
		// MUMPS Convection Linear Systems 
		mumpsMatrixConvection = new OpenSMOKE_MUMPS_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
		mumpsMatrixConvection->SetSparsityPattern(C_);
	}
	else if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
	{
		// LIS Convection Linear System
		lisMatrixConvection = new OpenSMOKE_LIS_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
		lisMatrixConvection->SetSparsityPattern(C_);
	}
    cout << "Done: Initialize convection sparse matrix" << endl; //getchar();

	// Open Files
	cout << " * Open statistics files..." << endl;
	statistics_ = new OpenSMOKE_KPP_SingleReactorStatistics(mix_.NumberOfSpecies(), NR_, "SequenceStatistics.out");
	statisticsConvective_ = new OpenSMOKE_KPP_ConvectiveNetworkStatistics(mix_.NumberOfSpecies(), NR_, "SequenceConvection.out");
    cout << "Done: Open statistics files" << endl; //getchar();
	
	if (data_.GlobalODE_SparseLinearSolver()  == KPP_SPARSESOLVER_PARDISO || data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
	{
		// PARDISO Global Linear Systems 
		pardisoMatrixGlobal = new OpenSMOKE_PARDISO_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW);
		pardisoMatrixGlobal->SetFillInReducingOrdering(data_.FillInReducingOrdering());
	}
	
	if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS || data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
	{
		// MUMPS Global Linear Systems
		mumpsMatrixGlobal = new OpenSMOKE_MUMPS_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW);
	}
	
	if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_LIS || data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
	{
		// LIS Global Linear Systems
		lisMatrixGlobal = new OpenSMOKE_LIS_Unsymmetric(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX_ROW_BY_ROW);
	}

	if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL || data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL)
	{
		// OpenSMOKE Gauss Siedel Global Linear Systems
		openSMOKEMatrixGlobal = new OpenSMOKE_KPP_BlockMatrixNetwork();
	}
    cout << "Done: Global Solver" << endl; //getchar();
	
	// Courant number (the maximum value is 1)
	ChangeDimensions(NR_, &timeSteps);
	for(int k=1;k<=NR_;k++)
		timeSteps[k] = reactors_[k].mass()/reactors_[k].M() * 1.;
	cout << "Time step analysis... " << endl;
	cout << " * Maximum time step: " << timeSteps.Max() << " s" << endl;
	cout << " * Minimum time step: " << timeSteps.Min() << " s" << endl;
	cout << " * Mean time step:    " << Mean(timeSteps) << " s" << endl;
	cout << endl;

	timeSteps *= data_.PredictorCorrector_CourantCorrectionCoefficient();

	// Open Output Files
	string sequence_name  = data_.nameFolderOutput() + "/Residuals/Sequence.out";
	string globalnls_name = data_.nameFolderOutput() + "/Residuals/GlobalNLS.out";
	string globalode_name = data_.nameFolderOutput() + "/Residuals/GlobalODE.out";
	

	
	fGlobalNLS_.open(globalnls_name.c_str());
	fGlobalNLS_.setf(ios::scientific);       
    fGlobalNLS_ << setw(12) << left << "Iter.(1)";
    fGlobalNLS_ << setw(12) << left << "Glob-NLS(2)";
    fGlobalNLS_ << setw(12) << left << "Loc-NLS(3)";
    fGlobalNLS_ << setw(16) << left << "Dummy(4)";
    fGlobalNLS_ << setw(16) << left << "Reduction(5)";
    fGlobalNLS_ << setw(12) << left << "Jacob(6)";        
    fGlobalNLS_ << setw(12) << left << "LS-Iter(7)";
    fGlobalNLS_ << setw(16) << left << "LS-Norm1(8)";
    fGlobalNLS_ << endl;
        
	fSequence_.open(sequence_name.c_str());
	fSequence_.setf(ios::scientific);        
        fSequence_ << setw(12) << left << "Iter.(1)";
        fSequence_ << setw(12) << left << "Glob-Seq.(2)";
        fSequence_ << setw(12) << left << "Loc-Seq.(3)";
        fSequence_ << setw(16) << left << "NormInf(4)";
        fSequence_ << setw(16) << left << "Norm1-Mean(5)";
        fSequence_ << setw(16) << left << "Norm2-Mean(6)";
        fSequence_ << setw(16) << left << "Norm1(7)";
        fSequence_ << setw(16) << left << "Norm2(8)";        
        fSequence_ << setw(12) << left << "Direct(9)";
        fSequence_ << setw(12) << left << "Newton(10)";
        fSequence_ << setw(12) << left << "ODE-I(11)";
        fSequence_ << setw(12) << left << "ODE-II(12)";
        fSequence_ << setw(12) << left << "ODE-III(13)";
        fSequence_ << setw(12) << left << "ODE-IV(14)";     
        fSequence_ << setw(12) << left << "Jacobian(15)";        
        fSequence_ << setw(12) << left << "Newton(16)";        
        fSequence_ << setw(12) << left << "Failures(17)";        
        fSequence_ << endl;
        
	fGlobalODE_.open(globalode_name.c_str());
	fGlobalODE_.setf(ios::scientific);        
        fGlobalODE_ << setw(12) << left << "Iter.(1)";
        fGlobalODE_ << setw(12) << left << "Glob-ODE(2)";
        fGlobalODE_ << setw(12) << left << "Loc-ODE(3)";
        fGlobalODE_ << setw(16) << left << "TimeStep(4)";
        fGlobalODE_ << setw(16) << left << "Reduction(5)";
        fGlobalODE_ << setw(12) << left << "Jacob(6)";        
        fGlobalODE_ << setw(12) << left << "LS-Iter(7)";
        fGlobalODE_ << setw(16) << left << "LS-Norm1(8)";
        fGlobalODE_ << endl;
        
        
        cpu_ = new CPUTimer(data_.nameFolderOutput() + "/Additional/CPUTime.out");
        cout << "Done: All" << endl; //getchar();
}

void OpenSMOKE_KPP_ReactorNetwork::AssemblingConvectionDiffusionMatrix()
{
	for (int k=1;k<=NR_;k++)
	{
		C_(k,k) = reactors_[k].M();
		
		for (int i=1;i<=reactors_[k].in().Size();i++)
			C_(k,reactors_[k].in()[i]) -= reactors_[k].cIn()[i];

		for (int i=1;i<=reactors_[k].neighbours().Size();i++)
			C_(k,reactors_[k].neighbours()[i]) -= reactors_[k].diffusion()[i];
	}
}

void OpenSMOKE_KPP_ReactorNetwork::AssemblingRightHandSides()
{
	for(int j=1;j<=mix_.NumberOfSpecies();j++)
		for (int k=1;k<=NR_;k++)
			RHS_[j][k] = reactors_[k].fIn()[j];
}

void OpenSMOKE_KPP_ReactorNetwork::SummaryConvectionDiffusionMatrix()
{
	int lower, upper;
	BzzVectorInt numEquationsForEachVariable;
	BzzVectorInt numVariablesForEachEquation;
	C_.GetNumEquationsForEachVariable(&numEquationsForEachVariable);
	C_.GetNumVariablesForEachEquation(&numVariablesForEachEquation);
	C_.FindBands(&lower,&upper);

	cout << endl;
	cout << "// *********************************************************** // " << endl;
	cout << "//           Convection/Diffusion Matrix Analysis              // " << endl;
	cout << "// *********************************************************** // " << endl;
	cout << "Number of equations:             " << C_.Rows() << endl;
	cout << "Number of non zero coefficients: " << C_.CountElements() << "/" << C_.Rows()*C_.Rows() << endl;
	cout << "Filling:                         " << C_.CountElements()/double(C_.Rows()*C_.Rows())*100. << "%" << endl;
	cout << "Mean number of non zero coefficients per equation: " << double(numVariablesForEachEquation.GetSumElements())/double(C_.Rows()) << endl;
	cout << "Max number of non zero coefficients per equation:  " << numVariablesForEachEquation.Max() << endl;
	cout << "Max number of equations per variables:             " << double(numEquationsForEachVariable.GetSumElements())/double(C_.Rows())<< endl;
	cout << "Mean number of equations per variables:            " << numEquationsForEachVariable.Max() << endl;
	cout << "Band width: upper(" << upper << ")  lower(" << lower << ")" << endl;
	cout << endl;
}

void OpenSMOKE_KPP_ReactorNetwork::CorrectingMassFlowRates()
{
	BzzVector inputMassFlowRate(NR_);
	BzzVector outputMassFlowRate(NR_);
	BzzMatrixSparse splitting(NR_, NR_);
	splitting.SetDiagonal(0, 1.);
	for(int k=1;k<=NR_;k++)
	{
		double mOut = reactors_[k].MassFlowOut();
		for (int j=1;j<=reactors_[k].out().Size();j++)
			splitting(k, reactors_[k].out()[j]) = -reactors_[k].cOut()[j]/mOut;
		inputMassFlowRate[k] = reactors_[k].fInTot();
	}
	Transpose(&splitting);
	
	if (data_.sparseLinearSolverMassFlowRate() == KPP_SPARSESOLVER_BZZ)
	{
		BzzFactorizedSparseGauss splittingFactorized;		
		splittingFactorized = splitting;
		outputMassFlowRate = inputMassFlowRate;
		Solve(&splittingFactorized, &outputMassFlowRate);
	}

	else if (data_.sparseLinearSolverMassFlowRate() == KPP_SPARSESOLVER_PARDISO)
	{
		OpenSMOKE_PARDISO_Unsymmetric pardisoSplitting(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
		pardisoSplitting.SetFillInReducingOrdering(data_.FillInReducingOrdering());
		pardisoSplitting.UnsetMessage();
		pardisoSplitting.SetSparsityPattern(splitting);
		pardisoSplitting.UpdateMatrix(splitting);
		pardisoSplitting.NumericalFactorization();
		pardisoSplitting.Solve(inputMassFlowRate, outputMassFlowRate);
		pardisoSplitting.Delete();
	}

	else if (data_.sparseLinearSolverMassFlowRate() == KPP_SPARSESOLVER_MUMPS)
	{
		OpenSMOKE_MUMPS_Unsymmetric mumpsSplitting(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
		mumpsSplitting.UnsetMessage();
		mumpsSplitting.SetSparsityPattern(splitting);
		mumpsSplitting.UpdateMatrix(splitting);
		mumpsSplitting.NumericalFactorization();
		mumpsSplitting.Solve(inputMassFlowRate, outputMassFlowRate);
		mumpsSplitting.Delete();
	}

	// Update out flow rates from each reactor
	Transpose(&splitting);
	for(int k=1;k<=NR_;k++)
	{
		for (int j=1;j<=reactors_[k].out().Size();j++)
			reactors_[k].setcOut(j, -splitting(k, reactors_[k].out()[j])*outputMassFlowRate[k] );
	}

	// Update in flow rates to each reactor
	for(int k=1;k<=NR_;k++)
	{
		for (int j=1;j<=reactors_[k].in().Size();j++)
		{
			int fromReactor = reactors_[k].in()[j];
			for (int i=1;i<=reactors_[fromReactor].out().Size();i++)
			{
				int toReactor = reactors_[fromReactor].out()[i];
				if (toReactor == k)
				{
					reactors_[k].setcIn(j, reactors_[fromReactor].cOut()[i]);
					break;
				}
			}
		}
	}
	
	// Update external output flows
	for (int k=1;k<=NR_;k++)
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
	else if (data_.PredictorCorrector_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
	{
		mumpsMatrixConvection->ResetCounters();
		mumpsMatrixConvection->UpdateMatrix(A_);
		mumpsMatrixConvection->NumericalFactorization();
		mumpsMatrixConvection->UnsetMessage();
	}
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
	// Statistics
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
		for(int j=1;j<=mix_.NumberOfSpecies();j++)
			for(int k=1;k<=NR_;k++)
				omegaOld[count++] = reactors_[k].omega()[j];
	}
	
	// 1. PREDICTOR 
	cout << " * Predictor..." << endl;
	for(int j=1;j<=mix_.NumberOfSpecies();j++)
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
		for(int j=1;j<=mix_.NumberOfSpecies();j++)
			sum += omegaHalf[(j-1)*NR_+k];
		for(int j=1;j<=mix_.NumberOfSpecies();j++)
			omegaHalf[(j-1)*NR_+k]/=sum;
	}

	// Communication: Update transport
	{
		cout << " * Communication" << endl;
		// a. Updating of mass fractions for every computational cell (k+1/2)
		for(int j=1;j<=mix_.NumberOfSpecies();j++)
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
		for(int j=1;j<=mix_.NumberOfSpecies();j++)
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
				cout << " * Corrector (continous, serial, dae)..." << endl;
				for(int k=1;k<=NR_;k++)
                                {
                                    odePool->Set(0, reactors_[k]);
                                    reactors_[k].SolveCSTR_CorrectorContinous_Smart(data_.SingleReactor_IntegrationTime(), *this, tmpMatrix, odePool->o(0));
                                }
			}
			else
			{
				cout << " * Corrector (continous, serial, ode)..." << endl;
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
				cout << "* Corrector (discrete, serial, dae)..." << endl;		
				for(int k=1;k<=NR_;k++)
                                {
                                        odePool->Set(0, reactors_[k]);
					reactors_[k].SolveCSTR_CorrectorDiscrete_Smart(data_.SingleReactor_IntegrationTime(), tmpMatrix, odePool->o(0));
                                }
			}
			else
			{
				cout << "* Corrector (discrete, serial, ode)..." << endl;	
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
                                        cout << "* Corrector (continous, openmp, dae, saved)..." << endl;		
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
					cout << "* Corrector (continous, openmp, ode, saved)..." << endl;		
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
					cout << "* Corrector (continous, openmp, dae, calculated)..." << endl;		
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
					cout << "* Corrector (continous, openmp, ode, calculated)..." << endl;
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
                                        cout << "* Corrector (discrete, openmp, dae, saved)..." << endl;		
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
					cout << "* Corrector (discrete, openmp, ode, saved)..." << endl;	
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
					cout << "* Corrector (discrete, openmp, dae, calculated)..." << endl;		
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
					cout << "* Corrector (discrete, openmp, ode, calculated)..." << endl;	
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
	}		// Corrector (OpenMP)
}

void OpenSMOKE_KPP_ReactorNetwork::TimeStepPolicy(double &deltat)
{
	// If no corrections usually occur
	if (statisticsConvective_->iterationsWithoutCorrections_ >= 5)
	{
		deltat *= data_.PredictorCorrector_TimeStepIncrementFactor();
		deltat  = min(deltat, data_.PredictorCorrector_MaxTimeStep());
		statisticsConvective_->iterationsWithoutCorrections_ = 0;
		statisticsConvective_->iterationsWithCorrections_    = 0;

		// New factorization
		SetTimeStep(deltat);
	}
	
	if (statisticsConvective_->iterationsWithCorrections_ >= 5)
	{
		deltat *= data_.PredictorCorrector_TimeStepReductionFactor();
		deltat  = min(deltat, data_.PredictorCorrector_MaxTimeStep());
		statisticsConvective_->iterationsWithoutCorrections_ = 0;
		statisticsConvective_->iterationsWithCorrections_    = 0;

		// New factorization
		SetTimeStep(deltat);
	}
}

void OpenSMOKE_KPP_ReactorNetwork::SolveSequenceCSTR()
{
	// Counters
	nLocalCSTRSequence_++;
	nGlobalCSTRSequence_++;

	// Statistics
	if (data_.SequenceCSTR_VerboseStatistics() == true)
		statistics_->Reset();
	
	// Discrete Updating
	if (data_.SequenceCSTR_UpdatingNetwork() == false)      // Discrete
	{
		

		for(int k=1;k<=NR_;k++)
			reactors_[k].DistributeMassFlowRates(*this);

		if (data_.iOpenMP() == true)
		{
                        int tid;
                        BzzMatrix* tmpMatrix = new BzzMatrix[data_.nThreads()];
                        for (int kk=0;kk<data_.nThreads();kk++)
                                ChangeDimensions(NumberOfSpecies(), NumberOfSpecies(), &tmpMatrix[kk]);
                                    
			if (data_.iSaveKineticConstants() == true)
			{
				cout << "* Sequence (openmp, discrete, saved)..." << endl;
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
				cout << "* Sequence (openmp, discrete, calculated)..." << endl;
				
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
		}
                else    ErrorMessage("Sequence Discrete Serial not implemented");
	}
	
        else if (data_.SequenceCSTR_UpdatingNetwork() == true)  // Continous
	{   
		if (data_.iOpenMP() == true)
		{
                        int tid;
                        BzzMatrix* tmpMatrix = new BzzMatrix[data_.nThreads()];
                        for (int kk=0;kk<data_.nThreads();kk++)
                                ChangeDimensions(NumberOfSpecies(), NumberOfSpecies(), &tmpMatrix[kk]);

			if (data_.iSaveKineticConstants() == true)
			{
				cout << "* Sequence (openmp, continous, saved)..." << endl;
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
				cout << "* Sequence (openmp, continous, calculated)..." << endl;	
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
		}
		else
		{
                        BzzMatrix tmpMatrix(NumberOfSpecies(), NumberOfSpecies());
                        
			cout << "* Sequence (serial, continous)..." << endl;
			
			int count=0;
			for(int k=1;k<=NR_;k++)
                        {
				if (count == 500)	{	cout << "Reactor " << k << endl; count=0;}
				count++;

                            odePool->Set(0,reactors_[k]);
				reactors_[k].SolveCSTR_CorrectorContinous_Smart(data_.SingleReactor_IntegrationTime(), *this, tmpMatrix, odePool->o(0));
                        }
			cout << "* Sequence completed..." << endl;
		}
	}

	//PrintDeFalco();
}

int OpenSMOKE_KPP_ReactorNetwork::AnalysisSequenceCSTR()
{
	int iDirectConvergence=0;
	int iNewtonConvergence=0;
	int iOdeFirstConvergence=0;
	int iOdeSecondConvergence=0;
	int iOdeThirdConvergence=0;
	int iOdeFourthConvergence=0;
	int nJacobianEvaluations=0;
	int nNewtonIterations=0;
	int nFailures = 0;
	double F1MaxNewton = 0;
	double F1MaxODE = 0;
	double F1MaxUnconverged = 0;
	double normInf=0;
	double F1Mean=0;

	for(int k=1;k<=NR_;k++)
	{			
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_DIRECT)
			iDirectConvergence++;
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_NEWTON)
		{
			iNewtonConvergence++;
			nJacobianEvaluations+=reactors_[k].status.nJacobianEvaluations;
			nNewtonIterations+=reactors_[k].status.nNewtonIterations;
			F1MaxNewton = max(F1MaxNewton, reactors_[k].status.norm1_over_nc);
		}
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_ODE_FIRST)
		{
			iOdeFirstConvergence++;
			F1MaxODE = max(F1MaxODE, reactors_[k].status.norm1_over_nc);
		}
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_ODE_SECOND)
		{
			iOdeSecondConvergence++;
			F1MaxODE = max(F1MaxODE, reactors_[k].status.norm1_over_nc);
		}
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_ODE_THIRD)
		{
			iOdeThirdConvergence++;
			F1MaxODE = max(F1MaxODE, reactors_[k].status.norm1_over_nc);
		}
		if (reactors_[k].status.convergence == KPP_SINGLEREACTOR_CONVERGENCE_ODE_FOURTH)
		{
			iOdeFourthConvergence++;
			F1MaxODE = max(F1MaxODE, reactors_[k].status.norm1_over_nc);
		}

		normInf  = max(normInf, reactors_[k].status.normInf);
		F1Mean  += reactors_[k].status.norm1_over_nc;
		nFailures += reactors_[k].status.failure;
	}

	// Mean residuals
	F1Mean/=double(NR_);

		
	cout.precision(2);
	cout << "Reactor distribution" << endl;
	cout << " * Already OK: " << setw(9) << fixed << right << double(iDirectConvergence)/double(NR_)*100.		<< "%" << setw(8) << right << iDirectConvergence << "/" << NR_ << endl;
	cout << " * Newton:     " << setw(9) << fixed << right << double(iNewtonConvergence)/double(NR_)*100.		<< "%" << setw(8) << right << iNewtonConvergence << "/" << NR_ << endl;
	cout << " * ODE(I):     " << setw(9) << fixed << right << double(iOdeFirstConvergence)/double(NR_)*100.	<< "%" << setw(8) << right << iOdeFirstConvergence << "/" << NR_ << endl;
	cout << " * ODE(II):    " << setw(9) << fixed << right << double(iOdeSecondConvergence)/double(NR_)*100.		<< "%" << setw(8) << right << iOdeSecondConvergence << "/" << NR_ << endl;
	cout << " * ODE(III):   " << setw(9) << fixed << right << double(iOdeThirdConvergence)/double(NR_)*100.	<< "%" << setw(8) << right << iOdeThirdConvergence << "/" << NR_ << endl;
	cout << " * ODE(IV):    " << setw(9) << fixed << right << double(iOdeFourthConvergence)/double(NR_)*100.	<< "%" << setw(8) << right << iOdeFourthConvergence << "/" << NR_ << endl;
	cout << endl;

	cout << "Additional data" << endl;
	cout << " * Jacobians per reactor:   " << setw(9) << fixed << right << double(nJacobianEvaluations)/double(max(1,iNewtonConvergence)) << endl;
	cout << " * Newtons per reactor:     " << setw(9) << fixed << right << double(nNewtonIterations)/double(max(1,iNewtonConvergence)) << endl;
	cout << " * Mass fractions failures: " << setw(9) << fixed << right << nFailures << endl;
	cout << endl;

	cout << "Residuals" << endl;
	setprecision(8);
	cout << " * Max residual Newton:    " << setw(9) << scientific << F1MaxNewton		<< endl;
	cout << " * Max residual ODE:       " << setw(9) << scientific << F1MaxODE			<< endl;
	cout << " * Max residual Failed:    " << setw(9) << scientific << F1MaxUnconverged	<< endl;
	cout << " * Mean residual:          " << setw(9) << scientific << F1Mean			<< endl;
	cout << " * Norm Inf:               " << setw(9) << scientific << normInf			<< endl;
	cout << endl;

	fSequence_ << setw(12)  << left << iteration_;
	fSequence_ << setw(12)  << left << nGlobalCSTRSequence_;
	fSequence_ << setw(12)  << left << nLocalCSTRSequence_;
	fSequence_ << setw(16) << left << scientific << residuals_->normInf();
	fSequence_ << setw(16) << left << scientific << residuals_->norm1()/NR_;
	fSequence_ << setw(16) << left << scientific << residuals_->norm2()/NR_;        
	fSequence_ << setw(16) << left << scientific << residuals_->norm1();
	fSequence_ << setw(16) << left << scientific << residuals_->norm2();
	fSequence_ << setw(12) << left << fixed << setprecision(2) << double(iDirectConvergence)/double(NR_)*100.;
	fSequence_ << setw(12) << left << fixed << setprecision(2) << double(iNewtonConvergence)/double(NR_)*100.;
	fSequence_ << setw(12) << left << fixed << setprecision(2) << double(iOdeFirstConvergence)/double(NR_)*100.;
	fSequence_ << setw(12) << left << fixed << setprecision(2) << double(iOdeSecondConvergence)/double(NR_)*100.;
	fSequence_ << setw(12) << left << fixed << setprecision(2) << double(iOdeThirdConvergence)/double(NR_)*100.;
	fSequence_ << setw(12) << left << fixed << setprecision(2) << double(iOdeFourthConvergence)/double(NR_)*100.;
	fSequence_ << setw(12) << left << fixed << setprecision(2) << double(nJacobianEvaluations)/double(max(1,iNewtonConvergence));
	fSequence_ << setw(12) << left << fixed << setprecision(2) << double(nNewtonIterations)/double(max(1,iNewtonConvergence));
	fSequence_ << setw(12) << left << fixed << setprecision(2) << nFailures;

        fSequence_ << endl;

	if ( F1Mean > 1.e-4 && nLocalCSTRSequence_>30)
	{
		fLog_ << "Sequence: Mean residual too large (" << F1Mean << ")" << endl;
		cout  << "Convergence Failure in CSTR Sequence: Mean residual too large" << endl;
		return -1;
	}

	if ( ((F1MaxODE>1.e-1) || (F1MaxNewton>1.e-1)) && nLocalCSTRSequence_>30)
	{
		fLog_ << "Sequence: Max residual too large (" << F1MaxODE << ", " << F1MaxNewton << ")" << endl;
		cout  << "Convergence Failure in CSTR Sequence: Max residual too large" << endl;
		return -1;
	}

	if ( double(nFailures)/double(NR_)*100 > 0.5  && nLocalCSTRSequence_>30)
	{
		fLog_ << "Sequence: Too many sum of mass fractions failures (" << double(nFailures)/double(NR_)*100 << "%)" << endl;
		cout << "Convergence Failure in CSTR Sequence: Too many sum of mass fractions failures..." << endl;
		return -1;
	}

	// Checking convergence index
	if ( (Residuals().RatioNormInf()>0.99) && (Residuals().RatioNorm1()>0.99) && nLocalCSTRSequence_>30)
	{
		fLog_ << "Sequence: Convergence indices satisfied..." << endl; 
		return 0;
	}
	// Checking convergence index
	if ( double(iOdeSecondConvergence+iOdeFirstConvergence)/double(NR_)*100.<1. && Residuals().RatioNormInf()>0.96)
	{
		fLog_ << "Sequence: Small enough number of ODE integrations..." << endl; 
		return 0;
	}

	return 1;
}
void OpenSMOKE_KPP_ReactorNetwork::InitializeGlobal(const KPP_SparseLinearSolver kind, const double absoluteTolerance, const double relativeTolerance)
{
	
	ChangeDimensions (NumberOfEquations(),  &globalOmega_);
	ChangeDimensions (NumberOfEquations(),  &globalRHS_);

	// Number of non-zero elements
	int numberOfNonZeroElements=0;
	for(int k=1;k<=NR_;k++)
		numberOfNonZeroElements += reactors_[k].nGlobalSparsityPattern();

	// Set PARDISO Matrix
	if (kind == KPP_SPARSESOLVER_PARDISO && iAlreadyGlobalPardiso_ == false)
	{
		pardisoMatrixGlobal->OpenMatrix(NR_*mix_.NumberOfSpecies(), numberOfNonZeroElements);
		for(int k=1;k<=NR_;k++)
//			pardisoMatrixGlobal->SetSparsityPattern(		reactors_[k].nGlobalSingleRowSparsityPattern(),
//														reactors_[k].globalIndicesSparsityColumns());
			pardisoMatrixGlobal->SetSparsityPattern(k, NumberOfSpecies(), 
														reactors_[k].mConvectionDiffusion(), 
														reactors_[k].iConvectionDiffusion(), 
														reactors_[k].nGlobalSingleRowSparsityPattern());

		pardisoMatrixGlobal->CloseMatrix();

		iAlreadyGlobalPardiso_ = true;
	}

	// Set MUMPS Matrix
	else if (kind == KPP_SPARSESOLVER_MUMPS && iAlreadyGlobalMumps_ == false)
	{
		mumpsMatrixGlobal->OpenMatrix(NR_*mix_.NumberOfSpecies(), numberOfNonZeroElements);
		for(int k=1;k<=NR_;k++)
//			mumpsMatrixGlobal->SetSparsityPattern(	reactors_[k].nGlobalSingleRowSparsityPattern(),
	//												reactors_[k].globalIndicesSparsityColumns() );
		
				mumpsMatrixGlobal->SetSparsityPattern(k, NumberOfSpecies(), 
														reactors_[k].mConvectionDiffusion(), 
														reactors_[k].iConvectionDiffusion(), 
														reactors_[k].nGlobalSingleRowSparsityPattern());

		mumpsMatrixGlobal->CloseMatrix();

		iAlreadyGlobalMumps_ = true;
	}

	// Set LIS Matrix
	else if (kind == KPP_SPARSESOLVER_LIS && iAlreadyGlobalLis_ == false)
	{
		lisMatrixGlobal->OpenMatrix(NR_*mix_.NumberOfSpecies(), numberOfNonZeroElements);
		
		for(int k=1;k<=NR_;k++)
//			lisMatrixGlobal->SetSparsityPattern(		reactors_[k].nGlobalSingleRowSparsityPattern(),
//													reactors_[k].globalIndicesSparsityColumns());
	
			lisMatrixGlobal->SetSparsityPattern(k, NumberOfSpecies(), 
														reactors_[k].mConvectionDiffusion(), 
														reactors_[k].iConvectionDiffusion(), 
														reactors_[k].nGlobalSingleRowSparsityPattern());

		lisMatrixGlobal->CloseMatrix();

		iAlreadyGlobalLis_ = true;
	}

	else if (kind == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL && iAlreadyGlobalGaussSiedel_ == false)
	{
		openSMOKEMatrixGlobal->SetAbsoluteTolerance(absoluteTolerance);
		openSMOKEMatrixGlobal->SetRelativeTolerance(relativeTolerance);
		openSMOKEMatrixGlobal->SetSparsityPattern(NR_, NumberOfSpecies());
		
		for(int k=1;k<=NR_;k++)
			openSMOKEMatrixGlobal->SetSparsityPattern(k, reactors_[k].iConvectionDiffusion());

		iAlreadyGlobalGaussSiedel_ = true;
	}
}

int OpenSMOKE_KPP_ReactorNetwork::GlobalNLS()
{
        // CPU time
        cpu_->SetNLSGlobal();
        cpu_->SetStartTimeNLSGlobal();
        
	indexGlobalNLS_++;
        localCounterNLS_=0;
        
	// First Guess solution
	for(int k=1;k<=NR_;k++)
		globalOmega_.SetBzzVector( (k-1)*mix_.NumberOfSpecies()+1, reactors_[k].omega() );

	// Solve Non Linear System
	OpenSMOKE_KPP_NewtonMethod_Manager *newton = new OpenSMOKE_KPP_NewtonMethod_Manager(globalOmega_, externalResiduals, externalSolutionNLS, externalResidualsAnalysis);
	
	// Non linear system option
	newton->SetMethod(data_.GlobalNLS_Method());
	newton->SetRelativeTolerance(data_.GlobalNLS_RelativeTolerance());
	newton->SetAbsoluteTolerance(data_.GlobalNLS_AbsoluteTolerance());
	newton->SetMaximumNumberOfIterations(data_.GlobalNLS_MaxIterations());
	newton->SetMaximumNumberOfArmijioIterations(data_.GlobalNLS_MaxArmijioIterations());

	fLog_ << "Index:             " << indexGlobalNLS_ << endl;
	fLog_ << "Rel/Abs Tol:       " << data_.GlobalNLS_RelativeTolerance() << " " << data_.GlobalNLS_AbsoluteTolerance() << endl;
	fLog_ << "Start Iteration:   " << iteration_ << endl;
        
	int info = newton->Solve(globalOmega_);

	fLog_ << "End Iteration:     " << iteration_ << endl;
	
	if (info == 3)          fLog_ << "Status:            super convergence..." << endl;
	if (info == 2)          fLog_ << "Status:            slow convergence..." << endl;
	else if (info == 1)     fLog_ << "Status:            maximum number of iterations.." << endl;
        else if (info == 0)     fLog_ << "Status:            convergence indices satisfied..." << endl;
	else if (info == -1)    fLog_ << "Status:            linear system failure..." << endl;
	else if (info == -2)    fLog_ << "Status:            Armijio failure..." << endl;
	
	//PrintDeFalco();

	return info;
}

int OpenSMOKE_KPP_ReactorNetwork::GlobalODE()
{
        // CPU time
        cpu_->SetODEGlobal();
        cpu_->SetStartTimeODEGlobal();
        
	indexGlobalODE_++;
        localCounterODE_=0;

	// First Guess solution
	for(int k=1;k<=NR_;k++)
		globalOmega_.SetBzzVector( (k-1)*mix_.NumberOfSpecies()+1, reactors_[k].omega() );

	int maxSubIterations = 1;
	double increasingFactor = 2.;
	double guessTimeStep = data_.GlobalODE_InitialTimeStep()*pow(10., double(indexGlobalODE_));
	double guessIncreasingFactor = data_.GlobalODE_TimeStepIncrementFactor();
	for(int i=1;i<=maxSubIterations;i++)
	{
		// Solve ODE System
		OpenSMOKE_KPP_ODE_Manager *ode = new OpenSMOKE_KPP_ODE_Manager(globalOmega_, externalResiduals, externalSolutionODE, externalResidualsAnalysis);
	
		ode->SetRelativeTolerance(data_.GlobalODE_RelativeTolerance());
		ode->SetAbsoluteTolerance(data_.GlobalODE_AbsoluteTolerance());
		ode->SetMaximumNumberOfIterations(data_.GlobalODE_MaxIterations());
		ode->SetStartingDeltaTime(guessTimeStep);
		ode->SetMaximumDeltaTime(guessIncreasingFactor);
		ode->SetIncreasingFactorDeltaTime(data_.GlobalODE_TimeStepIncrementFactor());
		ode->SetUpdatingFrequencyTimeStep(data_.GlobalODE_UpdatingFrequencyTimeStep());

		fLog_ << "Index:                " << indexGlobalODE_ << "-" << i << endl;
		fLog_ << "Rel/Abs Tol:          " << data_.GlobalODE_RelativeTolerance() << " " << data_.GlobalODE_AbsoluteTolerance() << endl;
		fLog_ << "Initial time step:    " << guessTimeStep << endl;
		fLog_ << "Increasing factor:    " << guessIncreasingFactor << endl;
		fLog_ << "Start iteration:      " << iteration_ << endl;
                
		int info = ode->Solve(globalOmega_);

		fLog_ << "End iteration:        " << iteration_ << endl;
		fLog_ << "Last time step (req): " << ode->LastRequestedDeltaTime() << endl;

		     if (info == 4)     fLog_ << "Status:               super Convergence..." << endl;
                else if (info == 3)     fLog_ << "Status:               very slow convergence..." << endl;
                else if (info == 2)     fLog_ << "Status:               slow convergence..." << endl;
                else if (info == 1)     fLog_ << "Status:               maximum number of iterations..." << endl;
		else if (info == 0)     fLog_ << "Status:               convergence satisfied..." << endl;
		else if (info == -1)    fLog_ << "Status:               linear system failure..." << endl;

		fLog_ << endl;

		if (info == 2)
		{
			guessTimeStep *= increasingFactor ;
			guessIncreasingFactor *= 1.25;
			continue;
		}

		if (info != 2)
			return info;
	}

//	PrintDeFalco();

	return 2;
}

void OpenSMOKE_KPP_ReactorNetwork::ApplyStatistics(BzzOdeStiffObject &o)
{
	statistics_->Analysis(o);
}

void OpenSMOKE_KPP_ReactorNetwork::SolveWithoutReactions()
{
	BzzVector x(NR_);
	BzzMatrix omega(NR_, mix_.NumberOfSpecies());
	
	OpenSMOKE_PARDISO_Unsymmetric pardisoSystem(OPENSMOKE_DIRECTSOLVER_SQUAREMATRIX);
	pardisoSystem.SetFillInReducingOrdering(data_.FillInReducingOrdering());
	pardisoSystem.SetSparsityPattern(C_);
	pardisoSystem.UpdateMatrix(C_);
	pardisoSystem.NumericalFactorization();
	pardisoSystem.UnsetMessage();

	for(int j=1;j<=mix_.NumberOfSpecies();j++)
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
	for (int j=1;j<=data_.mapsSpeciesNames().size();j++)
	{
		int index = mix_.recognize_species(data_.mapsSpeciesNames()[j-1]);

                ofstream fMap;
		openOutputFileAndControl(fMap, data_.nameFolderOutput() + "/Maps/Plot.plt." + mix_.names[index]);
		fMap.setf(ios::scientific);
		
		for (int k=1; k<=NR_; k++)
			fMap << reactors_[k].omega()[index] << endl;

		fMap.close();
	}

	if (data_.TraditionalLink()==true)
	{
		BzzVectorInt fromOriginalToCluster;
		BzzVector volumes;
		BzzLoad  fClusteringTopology("Input/ClusteringTopology.out");
		fClusteringTopology >> fromOriginalToCluster;
		fClusteringTopology >> volumes;
		fClusteringTopology.End();	

		for (int j=1;j<=data_.mapsSpeciesNames().size();j++)
		{
			int index = mix_.recognize_species(data_.mapsSpeciesNames()[j-1]);

		        ofstream fMap;
			openOutputFileAndControl(fMap, data_.nameFolderOutput() + "/Maps/Plot.plt." + mix_.names[index] + ".original");
			fMap.setf(ios::scientific);
		
			for(int k=1;k<=fromOriginalToCluster.Size();k++)
			{
				int i = fromOriginalToCluster[k];
				fMap << reactors_[i].omega()[index] << endl;
			}

			fMap.close();
		}
	}
}

void OpenSMOKE_KPP_ReactorNetwork::ExternalFeeds(double &massFlow, BzzVector &omegaFeeds, BzzVector &omegaElementalFeeds)
{
	// Species composition
	omegaFeeds = 0.;
	massFlow   = 0.;
	for (int k=1;k<=jExternalFeedReactors.Size();k++)
	{
		massFlow += reactors_[jExternalFeedReactors[k]].fInTot();
		for (int i=1;i<=omegaFeeds.Size();i++)
			omegaFeeds[i] += reactors_[jExternalFeedReactors[k]].fIn()[i];
	}
	
	double sum = omegaFeeds.GetSumElements();
	omegaFeeds /= sum;

	// Elemental composition
	ChangeDimensions(mix_.NumberOfElements(), &omegaElementalFeeds);	
	mix_.GetElementalMassFractionsFromSpeciesMassFractions(omegaElementalFeeds, omegaFeeds);
}

void OpenSMOKE_KPP_ReactorNetwork::ExternalOutput(double &massFlow, BzzVector &omegaOutput, BzzVector &omegaElementalOutput)
{
	// Species composition
	omegaOutput = 0.;
	massFlow    = 0.;
	for (int k=1;k<=jExternalOutputReactors.Size();k++)
	{
		double sum = reactors_[jExternalOutputReactors[k]].fOutTot();
		massFlow += sum;
		for (int i=1;i<=omegaOutput.Size();i++)
			omegaOutput[i] += sum*reactors_[jExternalOutputReactors[k]].omega()[i];
	}

	double sum = omegaOutput.GetSumElements();
	omegaOutput /= sum;

	// Elemental composition
	ChangeDimensions(mix_.NumberOfElements(), &omegaElementalOutput);	
	mix_.GetElementalMassFractionsFromSpeciesMassFractions(omegaElementalOutput, omegaOutput);
}

void OpenSMOKE_KPP_ReactorNetwork::MinMaxMassFractions(BzzVector &omegaMin, BzzVector &omegaMax)
{
	// Species composition
	omegaMin =  1e16;
	omegaMax = -1e-16;
	for (int k=1;k<=NR_;k++)
	{
        	for (int j=1;j<=NumberOfSpecies();j++)
                {
                    if (reactors_[k].omega()[j] < omegaMin[j]) omegaMin[j] = reactors_[k].omega()[j];  
                    if (reactors_[k].omega()[j] > omegaMax[j]) omegaMax[j] = reactors_[k].omega()[j];
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
		for (int j=1;j<=mix_.NumberOfSpecies();j++)
			omega[count++] = reactors_[k].omega()[j];
}

void OpenSMOKE_KPP_ReactorNetwork::ExtractMassFractions(BzzMatrix &omega) const
{
	for (int k=1;k<=NR_;k++)
		for (int j=1;j<=mix_.NumberOfSpecies();j++)
			omega[k][j] = reactors_[k].omega()[j];
}

void OpenSMOKE_KPP_ReactorNetwork::Residuals(BzzVector& x, BzzVector& f)
{
	cout << "Distribution of mass fractions..." << endl;
	for(int k=1;k<=NR_;k++)
		reactors_[k].setOmegaFromNetwork(k, x); 

	cout << "Residual evaluation..." << endl;
	for (int k=1;k<=NR_;k++)
	{
		int position = (k-1)*NumberOfSpecies() + 1;
		reactors_[k].Residuals(position, f, *this);
	}
}

void OpenSMOKE_KPP_ReactorNetwork::ResidualsAnalysis()
{
	iteration_++;
	globalTime_ += deltat;

	residuals_->Analysis();
	residuals_->WriteResidualsOnVideo();
}

void OpenSMOKE_KPP_ReactorNetwork::WriteReactorNetworkData()
{
	string fileName = data_.nameFolderOutput() + "ReactorNetwork.out";
	ofstream fOutput(fileName.c_str(), ios::out);
	fOutput.setf(ios::scientific);

/*	fOutput << "*********************************************************" << endl;
	fOutput << "*       Diagonal of Convection-Diffusion Matrix C       *" << endl;
	fOutput << "*********************************************************" << endl;
	for (int k=1;k<=NR_;k++)
	{
		fOutput << setw(8) <<  left << k;
		fOutput << setw(16) << left << reactors_[k].M();
		fOutput << endl;
	}
*/
/*	fOutput << endl;
	fOutput << "*********************************************************" << endl;
	fOutput << "*             Convection-Diffusion Matrix C             *" << endl;
	fOutput << "*********************************************************" << endl;
*/	{
		double* ptrVal;
		int i, j;
		double val;

		C_.BeginScanning();
		while(ptrVal = C_.Scanning(&i,&j,&val))
		{
			fOutput << setw(8)  << left << i;
			fOutput << setw(8)  << left << j;
			fOutput << setw(16) << left << *ptrVal;
			fOutput << endl;
		}
	}
/*
	fOutput << endl;
	fOutput << "*********************************************************" << endl;
	fOutput << "*                     External feed                     *" << endl;
	fOutput << "*********************************************************" << endl;
	for (int k=1;k<=NR_;k++)
	{
		fOutput << setw(8)  << left << k;
		fOutput << setw(16) << left << reactors_[k].fInTot();
		fOutput << endl;
	}
	*/
/*	fOutput << endl;
	fOutput << "*********************************************************" << endl;
	fOutput << "*                    External outflow                   *" << endl;
	fOutput << "*********************************************************" << endl;
	for (int k=1;k<=NR_;k++)
	{
		fOutput << setw(8)  << left << k;
		fOutput << setw(16) << left << reactors_[k].fOutTot();
		fOutput << endl;
	}

	fOutput << endl;
	fOutput << "*********************************************************" << endl;
	fOutput << "*          Mass     MassFlow     Mass Umbalances         *" << endl;
	fOutput << "*********************************************************" << endl;
	for (int k=1;k<=NR_;k++)
	{
		fOutput << setw(8)  << left << k;
		fOutput << setw(16) << left << reactors_[k].mass();
		fOutput << setw(16) << left << reactors_[k].MassFlowIn();
		fOutput << setw(16) << left << reactors_[k].MassUmbalance();
		fOutput << endl;
	}
*/
	fOutput.close();
}

void OpenSMOKE_KPP_ReactorNetwork::MassFlowErrors()
{
	double meanError=0.;
	double maxError=0.;
	for(int k=1;k<=NR_;k++)
	{
		double mIn  = reactors_[k].MassFlowIn();
		double mOut = reactors_[k].MassFlowOut();
		if (mIn+mOut<1.e-16)
			continue;

		double relative_error = fabs((mOut-mIn)/(mIn+mOut));
		meanError+=relative_error;
		if (relative_error>maxError)
			maxError = relative_error;
	}

	meanError /= double(NR_);

	cout << "        Mean relative error in mass balances: " << meanError << endl;
	cout << "        Max  relative error in mass balances: " << maxError  << endl;
}

KPP_Network_Status OpenSMOKE_KPP_ReactorNetwork::status() const
{
	return data_.networkStatus(); 
}

void CleanMassFractions(BzzVector& omega, const int nr, const int nc, const double threshold)
{
	for (int i=1;i<=nr;i++)
	{
		int k=(i-1)*nc;
		
		double sum = 0.;
		for (int j=k+1;j<=k+nc;j++)
		{
			// In case of strong negative value
			if (omega[j]<threshold)
			{
				cout << "Negative mass fraction: " << omega[j] << endl;
				cout << "Press enter to exit..." << endl;
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
        // CPU time
        cpu_->SetStartTimeNLSLocal();
        
        // Counters
        globalCounterNLS_++;
        localCounterNLS_++;        
        int linearSystemIterations      = 0;
        double linearSystemNorm         = 0.;        
        
	// Cleaning
	CleanMassFractions(x, NR_, NumberOfSpecies(), -1);

	// Statistics
	statisticsConvective_->Reset();

	// Communication: Update transport
	for(int k=1;k<=NR_;k++)
		reactors_[k].setOmegaFromNetwork(k, x); 

	// Assembling matrices
        BzzMatrix tmpMatrix(NumberOfSpecies(), NumberOfSpecies());
        BzzMatrix diagonalBlockMatrix(NumberOfSpecies(), NumberOfSpecies());
        for(int k=1;k<=NR_;k++)
	{
		double dummy_deltat=0.;
		reactors_[k].AssemblingLocalContribution(dummy_deltat, jacobianFlag, tmpMatrix, diagonalBlockMatrix);
		reactors_[k].AssemblingNonLocalRHS(dummy_deltat, x);
	}

	// Update global rhs
	for(int k=1;k<=NR_;k++)
	{
		int jReactor = (k-1)*mix_.NumberOfSpecies();
		for(int i=1;i<=mix_.NumberOfSpecies();i++)
			globalRHS_[jReactor+i] = reactors_[k].LocalRHS()[i]+reactors_[k].NonLocalRHS()[i];
	}

	// Update global matrix
	if (jacobianFlag == true)
	{
                globalCounterNLSJacobians_++;
                        
		if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
		{
			pardisoMatrixGlobal->ResetCounters();

			for(int k=1;k<=NR_;k++)
				pardisoMatrixGlobal->UpdateMatrix(reactors_[k].globalIndicesSparsityValues());

			pardisoMatrixGlobal->NumericalFactorization();
			pardisoMatrixGlobal->UnsetMessage();
		}

		else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
		{
			mumpsMatrixGlobal->ResetCounters();
		
			for(int k=1;k<=NR_;k++)
				mumpsMatrixGlobal->UpdateMatrix(reactors_[k].globalIndicesSparsityValues());

			mumpsMatrixGlobal->NumericalFactorization();
			mumpsMatrixGlobal->UnsetMessage();
		}

		else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
		{
			lisMatrixGlobal->ResetCounters();
		
			for(int k=1;k<=NR_;k++)
				lisMatrixGlobal->UpdateMatrix(reactors_[k].globalIndicesSparsityValues());

			lisMatrixGlobal->NumericalFactorization();
		//	lisMatrixGlobal->UnsetMessage();

			lisMatrixGlobal->EnableInitialSolution();
			lisMatrixGlobal->SetMaximumIterations(5000);
		}

		else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL)
		{
			for(int k=1;k<=NR_;k++)
			{
				BzzMatrix block(NumberOfSpecies(), NumberOfSpecies());
				BzzVector diagonal;
				reactors_[k].ReconstructBlockAndDiagonals(block, diagonal);
				openSMOKEMatrixGlobal->SetBlockAndDiagonals(k, block, diagonal);
				openSMOKEMatrixGlobal->SetRHS(k, reactors_[k].LocalRHS(), reactors_[k].NonLocalRHS());
			}
		}
	}
 
	// Linear system solution
	if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
		pardisoMatrixGlobal->Solve(globalRHS_, direction);
	else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
		mumpsMatrixGlobal->Solve(globalRHS_, direction);
	else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
	{
		direction = lastGlobalNLSSolution_;
		lisMatrixGlobal->SetInitialSolution(direction);
		lisMatrixGlobal->Solve(globalRHS_, direction);
		lastGlobalNLSSolution_ = direction;
	}
	else if (data_.GlobalNLS_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL)
	{
		direction = lastGlobalNLSSolution_;
		int flag = openSMOKEMatrixGlobal->GaussSiedel(direction);
		if (flag<0)	return -1;
		lastGlobalNLSSolution_ = direction;
	}

	// New values
	Sum(x, direction, &auxVector_NRxNC_1);
	
	// Cleaning machine errors
	CleanVectorOnTheBottom(0., -1.e-16, auxVector_NRxNC_1);
	CleanVectorOnTheTop(1., 1.+1.e-16, auxVector_NRxNC_1);

	// Correcting time step
	statisticsConvective_->UnsetRobust();
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
        fGlobalNLS_ << setw(12) << left << iteration_;
        fGlobalNLS_ << setw(12) << left << globalCounterNLS_;
        fGlobalNLS_ << setw(12) << left << localCounterNLS_;
        fGlobalNLS_ << setw(16) << left << 0;
        fGlobalNLS_ << setw(16) << left << a;
        fGlobalNLS_ << setw(12) << left << globalCounterNLSJacobians_;         
        fGlobalNLS_ << setw(12) << left << linearSystemIterations;
        fGlobalNLS_ << setw(16) << left << linearSystemNorm;
        fGlobalNLS_ << endl;           
        
	return 0;
}

int OpenSMOKE_KPP_ReactorNetwork::SolutionODE(BzzVector& omegaOld, BzzVector& step, double& dt, const bool jacobianFlag)
{
	double tFactorizing 		= 0.;
	

	double tStartSolution = BzzGetCpuTime();

        // CPU time
        cpu_->SetStartTimeODELocal();
        
        // Counters
        globalCounterODE_++;
        localCounterODE_++;        
        int linearSystemIterations      = 0;
        double linearSystemNorm         = 0.;
        
	// Cleaning
	CleanMassFractions(omegaOld, NR_, NumberOfSpecies(), -1.e-6);

	// Statistics
	statisticsConvective_->Reset();

	// Communication: Update transport
	double tStartAssembling = BzzGetCpuTime();  
	for(int k=1;k<=NR_;k++)
		reactors_[k].setOmegaFromNetwork(k, omegaOld); 

	// Assembling matrix and rhs
        BzzMatrix tmpMatrix(NumberOfSpecies(), NumberOfSpecies());
        BzzMatrix diagonalBlockMatrix(NumberOfSpecies(), NumberOfSpecies());
	for(int k=1;k<=NR_;k++)
	{
		reactors_[k].AssemblingLocalContribution(dt, jacobianFlag, tmpMatrix, diagonalBlockMatrix);
		reactors_[k].AssemblingNonLocalRHS(dt, omegaOld);
	}

	// Update global rhs
	for(int k=1;k<=NR_;k++)
	{
		int jReactor = (k-1)*mix_.NumberOfSpecies();
		for(int i=1;i<=mix_.NumberOfSpecies();i++)
			globalRHS_[jReactor+i] = reactors_[k].LocalRHS()[i]+reactors_[k].NonLocalRHS()[i];
	}
	double tEndAssembling = BzzGetCpuTime();  

	// Update global matrix
	if (jacobianFlag == true)
	{
                globalCounterODEJacobians_++;
                
		if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
		{
			pardisoMatrixGlobal->ResetCounters();
			for(int k=1;k<=NR_;k++)
				pardisoMatrixGlobal->UpdateMatrix(reactors_[k].globalIndicesSparsityValues());
			double tStartFactorization = BzzGetCpuTime();  
			pardisoMatrixGlobal->NumericalFactorization();
			double tEndFactorization = BzzGetCpuTime();
			tFactorizing = tEndFactorization-tStartFactorization;
			pardisoMatrixGlobal->UnsetMessage();
		}

		else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
		{
			mumpsMatrixGlobal->ResetCounters();
			for(int k=1;k<=NR_;k++)
				mumpsMatrixGlobal->UpdateMatrix(reactors_[k].globalIndicesSparsityValues());
			double tStartFactorization = BzzGetCpuTime();  
			mumpsMatrixGlobal->NumericalFactorization();
			double tEndFactorization = BzzGetCpuTime();
			tFactorizing = tEndFactorization-tStartFactorization;
			mumpsMatrixGlobal->UnsetMessage();
		}

		else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
		{
			lisMatrixGlobal->ResetCounters();
			for(int k=1;k<=NR_;k++)
				lisMatrixGlobal->UpdateMatrix(reactors_[k].globalIndicesSparsityValues());
			double tStartFactorization = BzzGetCpuTime();  
			lisMatrixGlobal->NumericalFactorization();
			double tEndFactorization = BzzGetCpuTime();
			tFactorizing = tEndFactorization-tStartFactorization;

			//	lisMatrixGlobal->UnsetMessage();
			lisMatrixGlobal->EnableInitialSolution();
		}

		else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL)
		{
			for(int k=1;k<=NR_;k++)
			{
				BzzMatrix block(NumberOfSpecies(), NumberOfSpecies());
				BzzVector diagonal;
				reactors_[k].ReconstructBlockAndDiagonals(block, diagonal);
				openSMOKEMatrixGlobal->SetBlockAndDiagonals(k, block, diagonal);
				openSMOKEMatrixGlobal->SetRHS(k, reactors_[k].LocalRHS(),reactors_[k].NonLocalRHS());
			}
		}
	}
        
	// Linear system solution
	double tStartBackSub = BzzGetCpuTime();  
	if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_PARDISO)
		pardisoMatrixGlobal->Solve(globalRHS_, step);
	else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_MUMPS)
		mumpsMatrixGlobal->Solve(globalRHS_, step);
	else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_LIS)
	{
		step = lastGlobalODESolution_;
		lisMatrixGlobal->SetInitialSolution(step);
		lisMatrixGlobal->Solve(globalRHS_, step);
		lastGlobalODESolution_ = step;
	}
	else if (data_.GlobalODE_SparseLinearSolver() == KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL)
	{
		step = lastGlobalODESolution_;
		int flag = openSMOKEMatrixGlobal->GaussSiedel(step);
		if (flag<0)	return -1;
		lastGlobalODESolution_ = step;
	}
	double tEndBackSub = BzzGetCpuTime();  
        
        // New solution
	Sum(omegaOld, step, &auxVector_NRxNC_1);

	// Cleaning machine errors
	CleanVectorOnTheBottom(0., -1.e-16, auxVector_NRxNC_1);
	CleanVectorOnTheTop(1., 1.+1.e-16, auxVector_NRxNC_1);
	
	// Correction of time step
	statisticsConvective_->SetRobust();
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
        
	double tEndSolution = BzzGetCpuTime();
       
        fGlobalODE_ << setw(12) << left << iteration_;
        fGlobalODE_ << setw(12) << left << globalCounterODE_;
        fGlobalODE_ << setw(12) << left << localCounterODE_;
        fGlobalODE_ << setw(16) << left << dt;
        fGlobalODE_ << setw(16) << left << a;
        fGlobalODE_ << setw(12) << left << globalCounterODEJacobians_;         
        fGlobalODE_ << setw(12) << left << linearSystemIterations;
        fGlobalODE_ << setw(16) << left << linearSystemNorm;
        fGlobalODE_ << endl;        
        
        // CPU time
        cpu_->SetEndTimeODELocal();
        cpu_->SetEndTimeODEGlobal();
        cpu_->SetEndTimeAllGlobal();          
        cpu_->WriteOnFile();

	cout << "CPU Times" << endl;
	cout << " * Solution:    " << tEndSolution - tStartSolution << " s" << endl;
	cout << " * Assembling:  " << tEndAssembling-tStartAssembling << " s (" << (tEndAssembling-tStartAssembling)/(tEndSolution - tStartSolution)*100. << "%)" << endl;
	cout << " * Factorizing: " << tFactorizing << " s (" << (tFactorizing)/(tEndSolution - tStartSolution)*100. <<"%)" << endl;
	cout << " * Back Sub:    " << tEndBackSub-tStartBackSub << " s (" << (tEndBackSub-tStartBackSub)/(tEndSolution - tStartSolution)*100. <<"%)" << endl;
	cout << endl;

        
	return 0;
}

int OpenSMOKE_KPP_ReactorNetwork::SequenceCSTR()
{
        // CPU time
        cpu_->SetSequence();
        cpu_->SetStartTimeSequenceGlobal();
                
	data_.SetNetworkStatus(KPP_NETWORK_STATUS_SEQUENTIAL_CSTR);
	nLocalCSTRSequence_ = 0;

	for (int k=1;k<=NR_;k++)
		reactors_[k].ResetStatus();

	double timeStartGlobal = BzzGetCpuTime();
	for (int i=1;i<=data_.SequenceCSTR_MaxIterations();i++)
	{
		// CPU Time
		cpu_->SetStartTimeSequenceLocal();
                
                // Solve the network
		SolveSequenceCSTR();
                
		// Residual analysis
		ResidualsAnalysis();
		
		// Analysis of reactors
		int info = AnalysisSequenceCSTR();

                // CPU time
                cpu_->SetEndTimeSequenceLocal();
                cpu_->SetEndTimeSequenceGlobal();
                cpu_->SetEndTimeAllGlobal();                  
                cpu_->WriteOnFile();

		if (info == 0)	return info;	// Means that convergence is reached
		if (info  < 0)	return info;	// Means that problems were met
	}

	return 1;	// Means that the maximum number of iterations was reached
}

void OpenSMOKE_KPP_ReactorNetwork::ReadBackupFile(const string fileName)
{
	if (data_.TraditionalLink() == false)
	{
		char name[Constants::NAME_SIZE];

		int numberReactorsBackup;
		int numberSpeciesBackup;
		BzzVector massFractions;

		BzzLoad fBackup('*', fileName.c_str());
		fBackup >> numberReactorsBackup;
		fBackup >> numberSpeciesBackup;

		cout << "Number of reactors (current): " << NR_ << endl;
		cout << "Number of reactors (backup):  " << numberReactorsBackup << endl;
		cout << "Number of species (current):  " << NumberOfSpecies() << endl;
		cout << "Number of species (backup):   " << numberSpeciesBackup << endl;

		if (NR_ != numberReactorsBackup)
			ErrorMessage("The number of reactors in backup file does not fit!");

		BzzVectorInt indices(numberSpeciesBackup);
		BzzVector omegaLocal(NumberOfSpecies());
	
		int countRecognized=0;
		for(int j=1;j<=numberSpeciesBackup;j++)
		{
			fBackup.fileLoad.read((char*) name, sizeof(name));
			indices[j] = mix_.recognize_species_without_exit(name);
			if (indices[j] > 0)	countRecognized++;
		}
		fBackup >> massFractions;
		fBackup.End();

		cout << "Number of species recognized:   " << countRecognized << "/" << NumberOfSpecies() << endl;

		double sum_max = 0.;
		double sum_min = 1.;

		int i=1;
		for(int k=1;k<=NR_;k++)
		{
			omegaLocal=0.;
			for(int j=1;j<=numberSpeciesBackup;j++)
			{
				if (indices[j]>0)	omegaLocal[indices[j]] = massFractions[i];
				i++;
			}
		
			double sum = omegaLocal.GetSumElements();
			if (sum > sum_max) sum_max = sum;
			else if (sum < sum_min) sum_min = sum;
			Product(1./sum, &omegaLocal);
		
			reactors_[k].setOmega(omegaLocal);
		}

		cout << "Min/max sums:   " << sum_min << "/" << sum_max << endl;
	}
	else
	{
		BzzVectorInt fromOriginalToCluster;
		BzzVector volumes;
		BzzLoad  fClusteringTopology("Input/ClusteringTopology.out");
		fClusteringTopology >> fromOriginalToCluster;
		fClusteringTopology >> volumes;
		fClusteringTopology.End();	

		BzzMatrix mass_fractions;
		int numberOfSpecies, originalNumberOfReactors;
		BzzLoad  fTraditionalBackup('*', "Output/Backup/Backup.start");
		fTraditionalBackup >> numberOfSpecies;
		fTraditionalBackup >> originalNumberOfReactors;

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
						cout << "Negative mass fraction: Reactor " << k << " Species " << j << " Omega " << aux[j] << endl;
						cout << "Press enter to exit..." << endl;
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
				cout << "Wrong sum of mass fractions: Reactor " << k << " Sum " << sum << endl;
				cout << "Press enter to exit..." << endl;
				getchar();
				exit(-1);
			}
			
			sum_mean += sum;
			if (sum>sum_max)	sum_max=sum;
			if (sum<sum_min)	sum_min=sum;

			aux /= sum;			

			reactors_[k].setOmega(aux);
		}		
		
		cout << "Statistics from Backup" << endl;
		cout << setprecision(8) << fixed << "Mean sums: " << sum_mean/double(NR_) << endl;
		cout << fixed << setprecision(8) << "Max sums:  " << sum_max << endl;
		cout << fixed << setprecision(8) << "Min sums:  " << sum_min << endl;
		
	}		
}

void OpenSMOKE_KPP_ReactorNetwork::WriteBackupFile(const string fileName)
{
	char name[Constants::NAME_SIZE];

	int i=1;
	for(int k=1;k<=NR_;k++)
		for(int j=1;j<=NumberOfSpecies();j++)
			auxVector_NRxNC_1[i++] = reactors_[k].omega()[j];

	BzzSave fBackup('*', fileName.c_str());
	fBackup << NR_;
	fBackup << NumberOfSpecies();
	for(int j=1;j<=NumberOfSpecies();j++)
	{
		strcpy(name, mix_.names[j].c_str());
		fBackup.fileSave.write((char*) name, sizeof(name));
	}
	fBackup << auxVector_NRxNC_1;
	fBackup.End();

	
	if (data_.TraditionalLink() == true)
	{
		BzzVectorInt fromOriginalToCluster;
		BzzVector volumes;
		BzzLoad  fClusteringTopology("Input/ClusteringTopology.out");
		fClusteringTopology >> fromOriginalToCluster;
		fClusteringTopology >> volumes;
		fClusteringTopology.End();	

		BzzMatrix mass_fractions(fromOriginalToCluster.Size(), NumberOfSpecies());
		for(int k=1;k<=fromOriginalToCluster.Size();k++)
			mass_fractions.SetRow(k, reactors_[fromOriginalToCluster[k]].omega());
		
		BzzSave  fTraditionalBackup('*', "Output/Backup/Backup.end");
		fTraditionalBackup << NumberOfSpecies();
		fTraditionalBackup << fromOriginalToCluster.Size();
		fTraditionalBackup << mass_fractions;
		fTraditionalBackup.End();
	}
}

void OpenSMOKE_KPP_ReactorNetwork::PrintDeFalco()
{
	bool iDeFalco = false;
	if (iDeFalco == true)
	{
		ofstream fDeFalco;
		ofstream fDeFalco2;
		fDeFalco.open("Reactors.info");
		fDeFalco2.open("Residuals.info");
		fDeFalco.setf(ios::scientific);
		fDeFalco2.setf(ios::scientific);
		for(int k=1;k<=NR_;k++)
			if (reactors_[k].temperature() > 1950.)
				reactors_[k].PrintDeFalco(k, C_, *this, fDeFalco, fDeFalco2);
		fDeFalco.close();
		fDeFalco2.close();
		cout << "Press enter to continue..." << endl;
//		getchar();
	}
}
