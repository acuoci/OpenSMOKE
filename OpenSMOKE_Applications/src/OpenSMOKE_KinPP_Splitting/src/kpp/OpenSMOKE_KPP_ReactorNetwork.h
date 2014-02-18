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

#ifndef OpenSMOKE_KPP_ReactorNetwork_H
#define OpenSMOKE_KPP_ReactorNetwork_H

#include "BzzMath.hpp"
#include "OpenSMOKE_KPP_Definitions.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "OpenSMOKE_KPP_SingleReactor.h"

class OpenSMOKE_PARDISO_Unsymmetric;
class OpenSMOKE_MUMPS_Unsymmetric;
class OpenSMOKE_LIS_Unsymmetric;
class OpenSMOKE_KPP_BlockMatrixNetwork;

class OpenSMOKE_KPP_DataManager;
class OpenSMOKE_KPP_ReactorNetwork;
class OpenSMOKE_KPP_SingleReactorStatistics;
class OpenSMOKE_KPP_ConvectiveNetworkStatistics;
class OpenSMOKE_KPP_ReactorNetwork_Residuals;

class ODE_Pool
{
public:
    
    ODE_Pool(const int nThreads, const int n)
    {
        n_ = n;
        nThreads_ = nThreads;
        odeSystem_ = new MyOdeSystem_KPP_ContinousReactor[nThreads_];
        o_ = new BzzOdeStiffObject[nThreads_];
    }

    void Set(const int i, OpenSMOKE_KPP_SingleReactor& reactor)
    {
        odeSystem_[i].assignReactor(&reactor);
    }
    
    void Initialize()
    {
        BzzVector omega_(n_);
        omega_ = 1./double(n_);
        for(int i=0;i<nThreads_;i++)
                o_[i].SetInitialConditions(omega_, 0., &odeSystem_[i]);
    }
    
    inline BzzOdeStiffObject& o(const int i) { return o_[i];}
    
    inline MyOdeSystem_KPP_ContinousReactor& odeSystem(const int i) { return odeSystem_[i];}
    
    
private:
    	
    int n_;
    int nThreads_;
    BzzOdeStiffObject* o_;
    MyOdeSystem_KPP_ContinousReactor* odeSystem_;
};

class OpenSMOKE_KPP_ReactorNetwork
{
public:

	// Constructors and destructors
	OpenSMOKE_KPP_ReactorNetwork(OpenSMOKE_ReactingGas& mixture, OpenSMOKE_KPP_DataManager& dictionary, ofstream& fLog, ofstream& fWarning);
	~OpenSMOKE_KPP_ReactorNetwork(void);

	// Access functions (read-only)
	inline int NumberOfReactors()		const { return NR_; }
	inline int NumberOfSpecies()		const { return mix_.NumberOfSpecies(); }
	inline int NumberOfEquations()		const { return NR_*mix_.NumberOfSpecies(); }
	inline int iteration()				const { return iteration_; }
	inline OpenSMOKE_ReactingGas& mix() const { return mix_; };
	KPP_Network_Status status() const; 

	inline OpenSMOKE_KPP_SingleReactor& reactors(const int k) { return reactors_[k]; }
	inline OpenSMOKE_KPP_DataManager& data() { return data_; }
	inline BzzMatrixSparse& C() { return C_; }

	inline OpenSMOKE_KPP_BlockMatrixNetwork& OpenSMOKEMatrixGlobal() { return *openSMOKEMatrixGlobal; };
	inline OpenSMOKE_KPP_ReactorNetwork_Residuals& Residuals() { return *residuals_; } ;

	// Functions which are called once
	void ReadFirstGuess();
	void ReadTopology();
	void BuildNetwork();

	// 1. Network solution: no reactions
	void SolveWithoutReactions();

	// 2. Network solution: only CSTR sequence
	int  SequenceCSTR();
	void ApplyStatistics(BzzOdeStiffObject &o);

	// 3. Network solution: predictor-corrector (differential+algebraic)
	void SetTimeStep(const double deltat_);
	void PredictorCorrector(const double deltat);
	void TimeStepPolicy(double &deltat);

	// 4. Global Systems
	void InitializeGlobal(const KPP_SparseLinearSolver kind, const double absoluteTolerance, const double relativeTolerance);
	int  GlobalODE();
	int  GlobalNLS();

	// Functions which are called more than once
	void WriteMassFractionMap();
	void ResidualsAnalysis();
	
	// BackUp
	void ReadBackupFile(const string fileName);
	void WriteBackupFile(const string fileName);

	// Utilities
	void ExternalFeeds(double &massFlow, BzzVector &omegaFeeds, BzzVector &omegaElementalFeeds);
	void ExternalOutput(double &massFlow, BzzVector &omegaOutput, BzzVector &omegaElementalOutput);
    void MinMaxMassFractions(BzzVector &omegaMin, BzzVector &omegaMax);

	void PrintDeFalco();

private:

	int NR_;
	BzzVectorInt jReduced;
	OpenSMOKE_ReactingGas& mix_;
	OpenSMOKE_KPP_DataManager& data_;
	OpenSMOKE_KPP_SingleReactor_KineticsManager* kinetics_;

	OpenSMOKE_KPP_SingleReactor* reactors_;
	OpenSMOKE_KPP_ReactorNetwork_Residuals* residuals_;

	BzzVectorInt jExternalFeedReactors;
	BzzVectorInt jExternalOutputReactors;

	BzzVectorInt* indicesReactorsThreads_ ;

	BzzMatrix residualMatrix_;
	BzzMatrixSparse C_;
	BzzMatrixSparse A_;
	BzzVector* RHS_;
	BzzVector* bStar_;
	double deltat;

	OpenSMOKE_PARDISO_Unsymmetric*	pardisoMatrixConvection;
	OpenSMOKE_MUMPS_Unsymmetric*	mumpsMatrixConvection;
	OpenSMOKE_LIS_Unsymmetric*		lisMatrixConvection;

	OpenSMOKE_PARDISO_Unsymmetric*	pardisoMatrixGlobal;
	OpenSMOKE_MUMPS_Unsymmetric*	mumpsMatrixGlobal;
	OpenSMOKE_LIS_Unsymmetric*		lisMatrixGlobal;
	OpenSMOKE_KPP_BlockMatrixNetwork* openSMOKEMatrixGlobal;

	void UpdateMassFractions(const int jSpecies, BzzVector &omega);
	void UpdateMassFractions(BzzVector &omega);
	void UpdateMassFractions(BzzMatrix &omega);
	void ExtractMassFractions(const int jSpecies, BzzVector &omega) const;
	void ExtractMassFractions(BzzVector &omega) const;
	void ExtractMassFractions(BzzMatrix &omega) const; 
	
	BzzVector auxVector_NRxNC_1;
	BzzVector lastGlobalNLSSolution_;
	BzzVector lastGlobalODESolution_;

	bool iAlreadyGlobalPardiso_;
	bool iAlreadyGlobalMumps_;
	bool iAlreadyGlobalLis_;
	bool iAlreadyGlobalGaussSiedel_;

	double globalTime_;
	int iteration_;

	OpenSMOKE_KPP_SingleReactorStatistics *statistics_;
	OpenSMOKE_KPP_ConvectiveNetworkStatistics *statisticsConvective_;

	// Global ode
	BzzVector globalOmega_;
	BzzVector globalRHS_;
	BzzVector timeSteps;

	ofstream& fLog_;
	ofstream& fWarning_;

	ofstream fSequence_;
	ofstream fGlobalODE_;
	ofstream fGlobalNLS_;
        
        int globalCounterODE_;
        int localCounterODE_;
        int globalCounterODEJacobians_;     
        
        int globalCounterNLS_;
        int localCounterNLS_;
        int globalCounterNLSJacobians_;           

	void SolveSequenceCSTR();
	int AnalysisSequenceCSTR();

	int nLocalCSTRSequence_;
	int nGlobalCSTRSequence_;
	int nLocalGlobalODE_;
	int indexGlobalODE_;
	int nLocalGlobalNLS_;
	int indexGlobalNLS_;
	double currentDeltaTimeGlobalODE_;
        
        CPUTimer* cpu_;
        
        ODE_Pool* odePool;

private:

	// Functions which are called once
	void MemoryAllocation();
	void AssemblingConvectionDiffusionMatrix();
	void AssemblingRightHandSides();
	void SummaryConvectionDiffusionMatrix();
	void CorrectingMassFlowRates();
	void MassFlowErrors();
	void WriteReactorNetworkData();
	void PreparingIndicesReactorsThreads();

public:

	void Residuals(BzzVector& x, BzzVector& f);
	int  SolutionNLS(BzzVector& x, BzzVector& dir, const bool jacobianFlag);
	int  SolutionODE(BzzVector& xOld, BzzVector& step, double& dt, const bool jacobianFlag);

private:

	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
};


#endif // OpenSMOKE_ReactorNetwork_H