/***************************************************************************
 *   Copyright (C) 2003-2008 by                                            *
 *   Guido Buzzi-Ferraris, Alessio Frassoldati and Alberto Cuoci		   *						   *
 *                                                                         *
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

#ifndef OPENSMOKE_CSTRNETWORK_H
#define OPENSMOKE_CSTRNETWORK_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_ReactingGas.h"


#define SYMBOLIC_KINETICS 0
#if SYMBOLIC_KINETICS==1
	#include "symbolickinetics/OpenSMOKE_SymbolicKinetics.h"
#endif

class OpenSMOKE_CSTRNetwork : public BzzBaseClass
{
	// Declaration of friend classes
	// ----------------------------------------------------------------------------
	friend class BzzNonLinearSystemSparse;
	friend class OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR;
	friend class OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR;

private:

	std::string firstGuessFile;
	std::string networkFile;
	std::string name_object;

	// CONSTANTS
	// ----------------------------------------------------------------------------
	static const double			R_CTOT;
	static const double			UR_CTOT;

	// Static Variables
	// ----------------------------------------------------------------------------
	static int count;
	static int countInScope;
	int	whoAmI;
	int printTasks;
	int printSubTasks;


	// Memory Management options
	// ----------------------------------------------------------------------------
	int memoTemperature;	// 1 salva su file 0 memorizza

	// Tolerances
	// ----------------------------------------------------------------------------
	double tolRel;
	double tolAbs;
	double MaxTolRel;
	double MaxTolAbs;

	double TminGlobal;
	double TmaxGlobal;
	double Tmin;
	
	BzzVector Tmax;
	BzzVector DeltaT;

	double MaxCoeffCorr;

	double Fluctuations_DeltaMax;
	double Fluctuations_TMax;
	double Fluctuations_CcMax;

	double MaxCountNewtonIterations;
	

	// REACTOR NETWORK DIMENSION
	// ----------------------------------------------------------------------------
	int	numComponents;
	int numReactions;
	int numCSTRReactors;
	int numComponentsReduced;
	int originalNumCSTRReactors;


	// REACTOR NETWORK VARIABLES
	// ----------------------------------------------------------------------------
	BzzVector	temperature;		// Reactor temperature [K]
	BzzVector	pressure;			// Reactor Pressure [atm]
	BzzVector	volume;				// Reactor Volume [m3]
	BzzVector	qanm;				// Reactor Squared Normalized temperature variance [-]
	BzzVector	massRate;			// Reactor convective mass flow rate [kg/s]
	BzzVector	massInput;			// Reactor External mass flow rate Input [kg/s]
	BzzVector	massOutput;			// Reactor External mass flow rate Output [kg/s]
	BzzVector	cTotm;				// Reactor Total concentration [kmol/m3]

	BzzVector	logTm;				// Reactor Auxiliary variable	log(T)
	BzzVector	loguRTm;			// Reactor Auxiliary variable	log(1/RT)
	BzzVector	feedInkRreactor;

	// Matrix Version
	BzzMatrix	massFractionsInReactorsSolution,		// Mass fractions
				massFractionsInReactors,				// Mass Fractions
				volumeMolecularWeightT,					// V * PM
				feed;									// Input mass flow rate

	// Vector Version
	BzzVector	massFractionsInReactorsSolution_Vector,	// Mass fractions
				massFractionsInReactors_Vector,			// Mass fractions
				volumeMolecularWeightT_Vector,			// V * PM
				feed_Vector;							// Input mass flow rate


	// COEFFICIENTS FOR NON LINEAR SYSTEM OF EQUATIONS
	// ----------------------------------------------------------------------------
	BzzMatrixSparse	Mg;
	BzzMatrixSparse	Ms;
	BzzVector					sM;
	BzzMatrixSparseLockedByRows	Ls;
	BzzMatrixSparseLockedByRows	Ld;
	//BzzFactorizedSparseLockedByRowsGauss Fg;


	// CONVERGENCE VARIABLES
	// ----------------------------------------------------------------------------
	double	F1Stop;							// Stop criterium on the residual
	int		maxNumberOfIterations_Global;	// Maximum number of iterations in Get First


	// DUMMY VARIABLES
	// ----------------------------------------------------------------------------
	BzzMatrix dummyMatrix;			// Auxiliar Matrix
	BzzVector dummyVector;			// Auxiliar Vector


	// SINGLE REACTOR VARIABLES
	// ----------------------------------------------------------------------------
	int		kReactor;	// Reactor Index
	double	T,			// Reactor temperature			(Fixed)
			P,			// Reactor pressure				(Fixed)
			logT,		// Reactor log(T)				(Fixed)
			loguRT,		// Reactor log(1/RT)			(Fixed)
			cTot,		// Reactor Total Concentration	(Fixed)
			wM;			// Reactor Molecular weight		(Variable)


	// SINGLE REACTOR VARIABLES
	// ----------------------------------------------------------------------------
	BzzVector	omegaReactor,		// mass fractions
				cReactor,			// concentrations
				xReactor,			// mole fractions
				RReactor;			// reaction rates


	// CLUSTERING ALGORITHM
	// ----------------------------------------------------------------------------
	BzzVectorInt		cluster;
	BzzVectorInt		giveClusterIndexFromCellIndex;
	BzzVectorIntArray	cstrOut;
	BzzVectorIntArray	cstrConnect;
	BzzVectorInt		iSpec;			// Reduced species indices in the complete scheme
	BzzVectorInt		cstrSequence;	// Reactor Sequence


	// GAS PHASE MIXTURE
	// ----------------------------------------------------------------------------
	OpenSMOKE_ReactingGas *Reactions;


	// INITIALIZATION ALGORITHM
	// ----------------------------------------------------------------------------
	void read_first_guess(std::string first);
	void reading_the_network_file(std::string fileNetwork, int &countInput);
	void build_cluster(double cicloCluster, int numCSTRReactors, int countInput, int &numCluster, BzzVectorInt &cstrClusterSize);
	void calculate_initial_massfractions_in_clusters(int numCluster, BzzMatrix &omegaReduced, BzzMatrix &aux);
	void calculate_clustering(std::string fileNetwork, int cicloDiffusion, int numCluster, BzzVectorInt &cstrClusterSize, BzzMatrix &auxMass, BzzMatrix &aux);
	void read_tolerances(const std::string fileTolerances, int numCSTRReactors, double cicloCluster);
	void complete_clustering(int countInput, int cicloDiffusion, int relaxation, int iaia, BzzMatrix &auxMass);

	void MemoTemperatureFunctions(const int kind);
	void CorrectionCoefficient_None(	BzzVector &u_temperature, BzzVector &log_temperature,
										BzzMatrix &matrix_correction_k1,
										BzzMatrix &matrix_correction_uKeq,
										BzzMatrix &matrix_correction_k2);
	void CorrectionCoefficient_SinExpansion(BzzVector &u_temperature, BzzVector &log_temperature,
											BzzMatrix &matrix_correction_k1, BzzMatrix &matrix_correction_uKeq,
											BzzMatrix &matrix_correction_k2);
	void CorrectionCoefficient_DeltaDirac(BzzVector &u_temperature, BzzVector &log_temperature,
												BzzMatrix &matrix_correction_k1, BzzMatrix &matrix_correction_uKeq,
												BzzMatrix &matrix_correction_k2);
	void CorrectionCoefficient_BetaPDF(	BzzVector &u_temperature, BzzVector &log_temperature,
										BzzMatrix &matrix_correction_k1,
										BzzMatrix &matrix_correction_uKeq, BzzMatrix &matrix_correction_k2);
	void CorrectionCoefficient_ClippedGaussianPDF(	BzzVector &u_temperature, BzzVector &log_temperature,
													BzzMatrix &matrix_correction_k1,
													BzzMatrix &matrix_correction_uKeq, BzzMatrix &matrix_correction_k2);

	void WriteCorrectionMap(const int kReaction, BzzMatrix &matrix_correction_coefficient, const std::string fileName);


	// TOLERANCES FOR CLUSTERING ALGORITHM
	// ----------------------------------------------------------------------------
	int		numTolT;
	BzzVector tInf;
	BzzVector tSup;
	BzzVector dt;
	BzzVector dtBase;
	double	tolRelC;
	double	tolAbsC;
	double	tolRelCBase;
	double	tolAbsCBase;


	// SOLVE NETWORK
	// ----------------------------------------------------------------------------
	void	start_solving_network(const int from_cfd_results);
	void	globalAlgorithm(void);
	void	odeOnlyAlgorithm();
	int		GetFirst(void);
	void	GetSecond(void);
	int		GetThird(void);
	void	GetFourth(const double tEnd);


	// GET RESIDUALS
	// ----------------------------------------------------------------------------
	void GetResiduals(BzzMatrix &massFractionsInReactors,BzzMatrix &residuals);
	void GetAllResiduals(BzzVector &m,BzzVector &residuals);
	void GetResiduals(BzzVector &omega,BzzVector &residuals);
	void GetJacobian(BzzVector &x,BzzMatrix &JJ);
	void GetJacobianForSingleReactor(int iReactor, BzzVector &cRes, BzzVector &R, BzzMatrix &dRC);


	// GET REACTION RATES
	// ----------------------------------------------------------------------------
	void	GetReactionsRateInAllReactorsFromMassFractions
				(BzzMatrix &massFractionsInReactors,
				BzzMatrix &reactionsRateInReactors);
	void	GetReactionsRateInAllReactorsFromMassFractions
				(BzzVector &massFractionsInReactorsV,
				BzzVector &reactionsRateInReactors);
	void	GetReactionRatesInSingleReactor(int iReactor,
				BzzVector &omega, BzzVector &molefractions,
				BzzVector &concentrations, BzzVector &ReactionRates);


	// OTHER FUNCTIONS
	// ----------------------------------------------------------------------------
	void GetDiagonalFactoredMatricesForLinearizedSistem(BzzMatrix &massFractionsInReactors);
	void GetDiagonalMatricesForLinearizedSistem(std::string file, BzzVector &mfV);
	void SwapMassFractionsInReactors(BzzMatrix *massFractionsInReactors);


	void calculate_massflowrate_in_reactors(BzzVector &inputMassFlowRate, BzzVector &outputMassFlowRate);

	int					iAnalyticalJacobian;
	SymbolicKinetics	analyticalJacobian;
	int sinExpansion;
	int iKindOfCorrection;
	
	double Hin_tot, Hout_tot;
	BzzVector  Hinput;
	BzzVectorInt	*FromClusterToCluster_Index;
	BzzVector *FromClusterToCluster_MassFlowRate;
	BzzVector *FromClusterToCluster_DiffusionFlowRate;
	void EnthalpyAnalysis(	const std::string fileNameEnthalpy,				const std::string fileNameDeltaEnthalpy,
							const std::string fileNameDifferenceEnthalpy,	const std::string fileNameInletEnthalpy,
							const std::string fileNameErrorEnthalpy);

	void ConnectionMatrix();

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

	void CheckInputfile(ifstream &file, std::string fileName);

public:

	// CONSTRUCTORS
	// ----------------------------------------------------------------------------
	OpenSMOKE_CSTRNetwork(void);

	void AssignKineticScheme(OpenSMOKE_ReactingGas &_Reactions);

	// OPERATORS
	// ----------------------------------------------------------------------------
	void operator()(const std::string fileNetwork,	const std::string first, const std::string fileTolerances,
					const double cicloCluster,		const int cicloDiffusion,	const int iaia,				const int relaxation, 
					const int _iAnalyticalJacobian, const SymbolicKinetics _analyticalJacobian, const int _iKindOfCorrection);

    void Maps();

	// SAVE NETWORK
	// ----------------------------------------------------------------------------
	void Save(std::string file);					// formatted
	void Save(char, std::string file);		    // binary
	void SaveTemperature(std::string file);		// formatted
	void SaveTemperature(char, std::string file);	// binary

	// UTILITIES
	// ----------------------------------------------------------------------------
	int		WhoAmI(void) const;
	void	SetTolRel(double tolr);
	void	SetTolAbs(double tola);
	void	SetMaxCorrectionCoefficient(const double _MaxCoeffCorr);
	void	SetDeltaTFluctuationsMaxDelta(const double _Fluctuations_DeltaMax);
	void	SetDeltaTFluctuationsMaxUserDefined(const double _Fluctuations_TMax);
	void	SetDeltaTFluctuationsMaxLocal(const double _Fluctuations_CcMax);

	void	SetMaxTolRel(double maxtolr);
	void	SetMaxTolAbs(double maxtola);
	void	SetMaxCountNewtonIterations(int maxCountNewtonIterations);
	void	SetF1Stop(double f1Stop);
	void	SetOdeOnly(const bool value);
	void	SetClusteringOnly(const bool value);

	void	SetMemoTemperature(void);
	int		GetNumComponents(void);
	int		GetNumReactions(void);
	int		GetNumCSTRReactors(void);
	void	SetTasksPrint(void);
	void	SetSubTasksPrint(void);
	void	SetSolution(BzzMatrix &omega);
	void	CleanTemperatureVariance(const int kind);

	void SetFluctuationsList(const std::string speciesFluctuations);
	vector<string> list_of_fluctuating_species;
	BzzVectorInt switchReactions;


	// STATIC FUNCTIONS
	// ----------------------------------------------------------------------------
	static int ObjectCount(void)		{	return count;	}
	static int ObjectCountInScope(void)	{	return countInScope;	}

	// PRINT FUNCTIONS
	// ----------------------------------------------------------------------------
	virtual void ObjectBzzPrint(void);
	void OutputPrint(int ciclo);
	void OutputFinalFile(std::string fileName);
	void GiveMeOutput(ofstream &fOutput, const int index);
	void GiveMeOutputLabel(ofstream &fOutput);

	// ASSIGNMENT OPERATORS
	// ----------------------------------------------------------------------------
	OpenSMOKE_CSTRNetwork &operator = (const OpenSMOKE_CSTRNetwork &rval);

	// MODIFYING FUNCTIONS
	// ----------------------------------------------------------------------------
	friend void Delete(OpenSMOKE_CSTRNetwork *result);

	// DESTRUCTORS
	// ----------------------------------------------------------------------------
	~OpenSMOKE_CSTRNetwork(void);

private:

	ofstream fWarning;
	
	BzzVector Tk_20000;
	BzzVector Tk_40000;
	BzzVector Tk_60000;

	void WarningLargeCorrectionCoefficient( const int k, const double T, const double g, const int iReaction, 
											const std::string stringReaction, const double EsuR, const double n, double &CoeffCorr);

	void WarningSmallCorrectionCoefficient( const int k, const double T, const double g, const int iReaction, 
											const std::string stringReaction, const double EsuR, const double n, double &CoeffCorr);

public:

	ofstream fHistory;
	ofstream fResiduals_1;
	ofstream fResiduals_2;
	ofstream fResiduals_3;
	ofstream fResiduals_4;
	int countSecondTotal;
	int countFourthTotal;
	int countThirdTotal;
	int countTuttoTotal;

	bool OnlyODE;
	bool OnlyClustering;
};

#endif // OPENSMOKE_CSTRNETWORK_H
