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
#include <petscvec.h>

class OpenSMOKE_PARDISO_Unsymmetric;
class OpenSMOKE_MUMPS_Unsymmetric;
class OpenSMOKE_LIS_Unsymmetric;
class OpenSMOKE_KPP_BlockMatrixNetwork;

class OpenSMOKE_KPP_DataManager;
class OpenSMOKE_KPP_Communicator;
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
	OpenSMOKE_KPP_ReactorNetwork(OpenSMOKE_ReactingGas* mixture, OpenSMOKE_KPP_DataManager& dictionary, std::ofstream& fLog, std::ofstream& fWarning);
	~OpenSMOKE_KPP_ReactorNetwork(void);

	// Access functions (read-only)
	inline int NumberOfReactors()		const { return NR_; }
	inline int* LocalNumberOfReactors()	const { return NR_P; }
	inline int NumberOfSpecies()		const { return mix_[0].NumberOfSpecies(); }
	inline int NumberOfEquations()		const { return NR_*mix_[0].NumberOfSpecies(); }
	inline int iteration()			const { return iteration_; }
	inline int procrank()			const { return procrank_; }
	inline int nprocs()			const { return nprocs_; }
	inline int numworkers()			const { return numworkers_; }
	inline int LoopIndex()			const { return LoopIndex_; }
	inline int* offs()			const { return offset; }
	inline int* CD_Size()			const { return ConvDiffSize; }
	inline int* Loc_Glob()			const { return LocToGlob; }
	inline int* countLocalOld()		const { return countLocal_Old; }
	inline const std:: vector<double*>& OmegaGlob()	const { return OmegaGlob_; }
	inline const std::vector<int*>& iConvDiff()	const { return iConvDiff_; }
	inline const std::vector<double*>& mConvDiff()	const { return mConvDiff_; }
	inline double* Mass_Umbalance()		const { return React_MassUmbalance; }
	inline double* Mass_FlowIn()		const { return React_MassFlowIn; }
        inline double* Hmix()                   const { return Hmix_;}
        inline double* delta_T()                const { return delta_T_;}
	inline long long int*    RHSplace()		const { return RHS_place_get_array; }
	inline Petsc64bitInt*	 placearray()		const { return place_array; }
	inline OpenSMOKE_ReactingGas* mix() const { return mix_; };
	inline OpenSMOKE_KPP_Communicator& communicator() const { return *communicator_; }
	KPP_Network_Status status() const; 
	inline bool newtonconvergence() const {return newtonconvergence_;}

	inline OpenSMOKE_KPP_SingleReactor& reactors(const int k) { return reactors_[k]; }
	inline OpenSMOKE_KPP_DataManager& data() { return data_; }
	inline BzzMatrixSparse& C() { return C_; }

	inline OpenSMOKE_KPP_BlockMatrixNetwork& OpenSMOKEMatrixGlobal() { return *openSMOKEMatrixGlobal; };
	inline OpenSMOKE_KPP_ReactorNetwork_Residuals& Residuals() { return *residuals_; } ;

	// Functions which are called once
	void ReadFirstGuess();
	void ReadTopology();
	void BuildNetwork();
	void SetCommunicator(OpenSMOKE_KPP_Communicator* communicator);

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
        void CheckEnergyBalances();
	
	// BackUp
	void ReadBackupFile(const std::string fileName);
	void WriteBackupFile(const std::string fileName);

	// Utilities
	void ExternalFeeds(double &massFlow, BzzVector &omegaFeeds, BzzVector &omegaElementalFeeds);
	void ExternalOutput(double &massFlow, BzzVector &omegaOutput, BzzVector &omegaElementalOutput);
    void MinMaxMassFractions(BzzVector &omegaMin, BzzVector &omegaMax);

	void PrintDeFalco();

private:

	int NR_;
	int* offset;
	int* NR_P;
	int MASTER;
	int* counter;
	int* Counter_Glob;
	int *ReactID, *LocToGlob;
	int CountElem;
	int* countLocal_Old;
	int mtype, source, FROM_MASTER, FROM_WORKER;
	double *mOut, *inputMassFlowRate_aux;
	double *outputMassFlowRate_aux;
	double *React_transport;
	double *Tmin, *Tmax;
	double Tmin_all, Tmax_all;

	BzzVectorInt jReduced;
	OpenSMOKE_ReactingGas* mix_;
	OpenSMOKE_KPP_DataManager& data_;
	OpenSMOKE_KPP_SingleReactor_KineticsManager* kinetics_;

	OpenSMOKE_KPP_SingleReactor* reactors_;
	OpenSMOKE_KPP_SingleReactor* reactors_fake;
	OpenSMOKE_KPP_ReactorNetwork_Residuals* residuals_;

	BzzVectorInt jExternalFeedReactors;
	BzzVectorInt jExternalOutputReactors;
	int *jExternalFeedReactors_Loc, *jExternalOutputReactors_Loc;
	int *jExternalFeedReactors_aux;
	int *jExternalOutputReactors_aux;
	int jEFR_size, jEOR_size, jEFRsum, jEORsum;
	int *jEFR_sizeAll, *jEOR_sizeAll, *jEFRoffset, *jEORoffset;	

	BzzVectorInt* indicesReactorsThreads_ ;

	BzzMatrix residualMatrix_;
	BzzMatrixSparse C_;
	BzzMatrixSparse A_;
	BzzVector* RHS_;
	BzzVector* bStar_;
	double deltat;
        
        double *Hmix_;
        double *Hin_, *Hout_;
        double *delta_T_;
        double *Delta_H_start, *Delta_H_end;

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

	//Parallel utilities
	void DistributeCells();
	void TopologyFakeReading();
	void DeleteTopologyParametersSize();
	void TopologyParametersSize();
	void TopologyParameters();
	void ReactorAlternateRank();
	void FeedsAndOutputs();
	void VariableTopologyParameters();
	void ConvDiffParameters();
	void SendAndBroadcastOmega();
	void WriteReactorNetworkDataArray();
	void SendCSTRInfoToMaster();
	void InitializeCSTRVectors();
	void GetSingleReactorSparsityPattern();
	void CleanMassFractions(BzzVector& omega, const int nr, const int nc, const double threshold);
	void InitializeNetworkData();
	void UpdateRHS();
	void UpdateLISMatrix();
	void CommunicationPattern();
	void InitializePetscMatrix();
	void InitializePetscRHS();
	void InitializePetscResiduals();
	void CleanMemory();
	void InitializePlaceArray();

	bool* React_communication;
	std::vector<bool*> React_comm_master;
	int* count_comm;

	int* ReactOutSize;
	int* ReactNeighSize;
	int* ReactInSize;
	int ReactfInSize;
	int* ConvDiffSize;
	int* SparsityPattern;
	int* globalIndicesSparsityValuesSize;
	long long int* gISVSizeCumulated;
	int* gISVSize;

	int* IndicesSparsityValuesSize;

	double *React_M;
	double *React_fInTot;
	double *React_fOutTot;
	double *React_mass;
	double *React_MassFlowIn;
	double *React_MassUmbalance;
	double *omega_ptr;
	double *x_ptr;
	double *direction_ptr;
	double *ReactorRHS_;
	double *localRHS_;

	std::vector<int*> React_out;
	std::vector<int*> React_in;
	std::vector<int*> React_neighbours;
   	std::vector<double*> React_cOut;
    	std::vector<double*> React_cIn;
    	std::vector<double*> React_diffusion;
    	std::vector<double*> React_fIn;
	std::vector<double*> Split_Out;
	std::vector<double*> Global_RHS;
	std::vector<double*> globalIndicesSparsityValues;

	std::vector<double*> IndicesSparsityValues;
 
	std::vector<int*> iConvDiff_;
	std::vector<double*> mConvDiff_;

	std::vector<double*> OmegaGlob_;

	Vec PetscMatrix_;
	Vec PetscRHS_;
	Vec OmegaFractions_;
	Vec EqResiduals_;
	PetscInt nGet_, nSet_, nGetAll_, nSetAll_, nRHSset_, nRHSget_, nResiduals_, nRHS_, nMatrix_, nOmega_;
	PetscInt lowindex_get, highindex_get, lowindex_set, highindex_set, localsize_get, localsize_set, upper_matrix, lower_matrix, lowindex_omega, highindex_omega;
	PetscInt RHS_lowindex_get, RHS_highindex_get, RHS_lowindex_set, RHS_highindex_set, RHS_lowindex, RHS_highindex;
	PetscScalar *value, *RHS_value;
	std::vector<long long int*> place_get;
	std::vector<long long int*> matrix_place;
	std::vector<long long int*> place_vector;
	Petsc64bitInt* place_array;
	long long int *place_set;
	long long int *place;
	double* value_get;
	std::vector<long long int*> RHS_place_get;
	std::vector<long long int*> RHS_place;
	long long int* RHS_place_get_array, *matrix_place_array, *RHS_place_array;
	long long int *RHS_place_set;
	long long int matrix_place_size;
	double *Petsc_values;

	//Convergence utilities

	int iDirectConvergence;
	int iNewtonConvergence;
	int iOdeFirstConvergence;
	int iOdeSecondConvergence;
	int iOdeThirdConvergence;
	int iOdeFourthConvergence;
	int nJacobianEvaluations;
	int nNewtonIterations;
	int nFailures;
	bool CSTRConvergengeAlert;
	double F1MaxNewton;
	double F1MaxODE;
	double F1MaxUnconverged;
	double normInf;
	double F1Mean;

        BzzVectorInt Process_iDirectConvergence;
	BzzVectorInt Process_iNewtonConvergence;
	BzzVectorInt Process_iOdeFirstConvergence;
	BzzVectorInt Process_iOdeSecondConvergence;
	BzzVectorInt Process_iOdeThirdConvergence;
	BzzVectorInt Process_iOdeFourthConvergence;
	BzzVectorInt Process_nJacobianEvaluations;
	BzzVectorInt Process_nNewtonIterations;
	BzzVectorInt Process_nFailures;
	BzzVector Process_F1MaxNewton;
	BzzVector Process_F1MaxODE;
	BzzVector Process_F1MaxUnconverged;
	BzzVector Process_normInf;
	BzzVector Process_F1Mean;
	
	BzzVector auxVector_NRxNC_1;
	BzzVector lastGlobalNLSSolution_;
	BzzVector lastGlobalODESolution_;

	bool iAlreadyGlobalPardiso_;
	bool iAlreadyGlobalMumps_;
	bool iAlreadyGlobalLis_;
	bool iAlreadyGlobalGaussSiedel_;
	bool newtonconvergence_;
        bool iBooleanMatrix_;

	double globalTime_;
	int iteration_;
	int LoopIndex_;
	double NormInfConvCriterium, Norm1ConvCriterium, ReactIterCriterium;

	int network_index_, network_ID_;
	BzzVectorInt iConvDiffMaster;
	BzzVector mConvDiffMaster;
	int *iConvDiffMaster_aux;
	double *mConvDiffMaster_aux;

	OpenSMOKE_KPP_SingleReactorStatistics *statistics_;
	OpenSMOKE_KPP_ConvectiveNetworkStatistics *statisticsConvective_;

	// Global ode
	BzzVector globalOmega_;
	BzzVector globalRHS_;
	BzzVector timeSteps;
	double *timeStepsArray;
	double lastODEStep_;

	std::ofstream& fLog_;
	std::ofstream& fWarning_;

	std::ofstream fSequence_;
	std::ofstream fGlobalODE_;
	std::ofstream fGlobalNLS_;
        
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

	int nprocs_, procrank_, numworkers_;
        
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
	void PrepareCDMatrixCommunication();

public:

	void Residuals(BzzVector& x, BzzVector& f);
	void Residuals(BzzVector& x, Vec& f);
	void Residuals(Vec& x, Vec& f);
	int  SolutionNLS(BzzVector& x, BzzVector& dir, const bool jacobianFlag);
	int  SolutionODE(BzzVector& xOld, BzzVector& step, double& dt, const bool jacobianFlag);

private:

	void ErrorMessage(const std::string message_);
	void WarningMessage(const std::string message_);
	OpenSMOKE_KPP_Communicator* communicator_;
};


#endif // OpenSMOKE_ReactorNetwork_H
