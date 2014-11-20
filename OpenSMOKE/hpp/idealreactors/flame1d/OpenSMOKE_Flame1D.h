/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public Licenvse as published by  *
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

#if !defined(FLAME_OPPOSED_PREMIXED)
#define FLAME_OPPOSED_PREMIXED

#include "OpenSMOKE.hpp"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_NLS_DAE_ODE_Classes.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_DataManager.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_OscillatingBoundary.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_ScheduleClass.h"

class OpenSMOKE_Flame1D_Solution;
class OpenSMOKE_Flame1D_OpposedFlameManager;

enum flame1d_ignition { FLAME1D_IGNITION_NONE, FLAME1D_IGNITION_SPARK, FLAME1D_IGNITION_NO_SPARK, FLAME1D_IGNITION_START };

class OpenSMOKE_Flame1D 
{
	friend class OpenSMOKE_Flame1D_Solution;
	friend class OpenSMOKE_Flame1D_OpposedFlameManager;

public:

	// Opposed Flames
	// --------------------------------------------------------------------------
	OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ALL							DAE_Opposed_ALL;
	OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOMOMENTUM					DAE_Opposed_NOMOMENTUM;
	OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOENERGY						DAE_Opposed_NOENERGY;
	OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ONLYMOMENTUM						DAE_Opposed_ONLYMOMENTUM;
	OpenSMOKE_Flame1D_MyDaeSystem_Opposed_SOOT_ALL						DAE_Opposed_SOOT_ALL;
	OpenSMOKE_Flame1D_MyDaeSystem_Opposed_QMOM_ALL						DAE_Opposed_QMOM_ALL;

	OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ALL						NLS_Opposed_ALL;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Opposed_ALL				NLS_Reduced_Opposed_ALL;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYUGH					NLS_Opposed_ONLYUGH;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYMASSFRACTIONS		NLS_Opposed_ONLYMASSFRACTIONS;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_COLDREDUCED				NLS_Opposed_COLDREDUCED;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYT					NLS_Opposed_ONLYT;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_NOENERGY				NLS_Opposed_NOENERGY;

	// Twin Flames
	// --------------------------------------------------------------------------
	OpenSMOKE_Flame1D_MyDaeSystem_Twin_ALL							DAE_Twin_ALL;
	OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOMOMENTUM					DAE_Twin_NOMOMENTUM;
	OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOENERGY						DAE_Twin_NOENERGY;
	OpenSMOKE_Flame1D_MyDaeSystem_Twin_SOOT_ALL						DAE_Twin_SOOT_ALL;
	OpenSMOKE_Flame1D_MyDaeSystem_Twin_QMOM_ALL						DAE_Twin_QMOM_ALL;

	OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ALL					NLS_Twin_ALL;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Twin_ALL			NLS_Reduced_Twin_ALL;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYUGH				NLS_Twin_ONLYUGH;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYMASSFRACTIONS		NLS_Twin_ONLYMASSFRACTIONS;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_COLDREDUCED			NLS_Twin_COLDREDUCED;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYT					NLS_Twin_ONLYT;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_NOENERGY				NLS_Twin_NOENERGY;
	
	// Premixed Flames
	// --------------------------------------------------------------------------
	OpenSMOKE_Flame1D_MyDaeSystem_Premixed_ALL					DAE_Premixed_ALL;
	OpenSMOKE_Flame1D_MyDaeSystem_Premixed_NOENERGY				DAE_Premixed_NOENERGY;
	OpenSMOKE_Flame1D_MyDaeSystem_Premixed_FLAMESPEED			DAE_Premixed_FLAMESPEED;
	OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_ALL				DAE_Premixed_QMOM_ALL;
	OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_NOENERGY		DAE_Premixed_QMOM_NOENERGY;
	OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_ALL				DAE_Premixed_SOOT_ALL;
	OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_NOENERGY		DAE_Premixed_SOOT_NOENERGY;

	OpenSMOKE_Flame1D_MyOdeSystem_Premixed_NOENERGY				ODE_Premixed_NOENERGY;
	OpenSMOKE_Flame1D_MyOdeSystem_Premixed_ALL					ODE_Premixed_ALL;

	OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_ALL			NLS_Premixed_ALL;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_NOENERGY		NLS_Premixed_NOENERGY;
	OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_FLAMESPEED		NLS_Premixed_FLAMESPEED;


	OpenSMOKE_Flame1D_MyOdeSystem_SingleReactor_Isothermal		ODE_SingleReactor_Isothermal;

public:

	OpenSMOKE_Flame1D();
	
	void SetName(const std::string name);
	void Assign(OpenSMOKE_Flame1D_DataManager *_data);
	void Assign(OpenSMOKE_Flame1D_ScheduleClass *_operations);
	void Assign(OpenSMOKE_ReactingGas  *_mix);
	void Assign(OpenSMOKE_GlobalKinetics *_global);

	void Run();
	void FoldersAndFilesManager();
	void FoldersAndFilesManager(const std::string backupFolder);
	
	void ReRun();
	void PasteFromExternalSolution(OpenSMOKE_Flame1D_Solution &solution);

	
	// Premixed Flames
	// --------------------------------------------------------------------------
	int Np, Ni;
	int nBlock, NC, NR;
	

	OpenSMOKE_Grid1D grid;
	OpenSMOKE_Flame1D_DataManager	*data;
	OpenSMOKE_Flame1D_ScheduleClass	*operations;
	OpenSMOKE_ReactingGas			*mix;


	// Physical variables
	// -------------------------------------------------------------------------------
	BzzVector U, G, H, T;
	BzzMatrix W;
	BzzMatrix X;
	BzzVector dU, dG, dH, dT;
	BzzMatrix dW;
	BzzMatrix x_elemental;
	BzzMatrix omega_elemental;
	BzzVector Z;
	BzzVector rho;


	// Initialize
	// -------------------------------------------------------------------------------
	void setup();

	void initializeVariables(const flame1d_ignition kind);

	// Non Linear Systems - Interfaces
	// ---------------------------------------------------------------------------
	int  solveNLS_Opposed(int Hot_or_Cold, const flame1d_model kind);
	int  solveNLS_Twin(int Hot_or_Cold, const flame1d_model kind);
	void solveNLS_Premixed(const flame1d_model kind);

	int solveNLS_Reduced_Twin(int Hot_or_Cold, const flame1d_model NLS_KIND);
	int solveNLS_Reduced_Opposed(int Hot_or_Cold, const flame1d_model NLS_KIND);

	int solveSensitivity(const flame1d_model NLS_KIND, BzzSave &fBinary);

	void nonLinearSystem_Opposed_All(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Reduced_Opposed_All(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Opposed_OnlyUGH(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Opposed_OnlyMassFractions(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Opposed_ColdReduced(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Opposed_OnlyT(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Opposed_NoEnergy(BzzVector &x, BzzVector &f);

	void nonLinearSystem_Twin_All(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Reduced_Twin_All(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Twin_OnlyUGH(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Twin_OnlyMassFractions(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Twin_ColdReduced(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Twin_OnlyT(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Twin_NoEnergy(BzzVector &x, BzzVector &f);

	void nonLinearSystem_Premixed_All(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Premixed_NoEnergy(BzzVector &x, BzzVector &f);
	void nonLinearSystem_Premixed_FlameSpeed(BzzVector &x, BzzVector &f);

	// DAE Systems - Interfaces
	// ---------------------------------------------------------------------------
	void solveDAE_Opposed(int Hot_Or_Cold, const flame1d_model kind, double tEnd);
	void solveDAE_Twin(int Hot_Or_Cold, const flame1d_model kind, double tEnd);
	void solveDAE_Premixed(const flame1d_model kind, double tEnd);

	void DAESystem_Opposed_All(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Opposed_NoMomentum(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Opposed_NoEnergy(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Opposed_OnlyMomentum(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Opposed_SOOT_ALL(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Opposed_QMOM_ALL(BzzVector &x, double t, BzzVector &f);

	void DAESystem_Twin_All(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Twin_NoMomentum(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Twin_NoEnergy(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Twin_SOOT_ALL(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Twin_QMOM_ALL(BzzVector &x, double t, BzzVector &f);

	void DAESystem_Premixed_All(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Premixed_NoEnergy(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Premixed_FlameSpeed(BzzVector &x, double t, BzzVector &f);

	void DAESystem_Premixed_QMOM_All(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Premixed_QMOM_NoEnergy(BzzVector &x, double t, BzzVector &f);
	
	void DAESystem_Premixed_SOOT_ALL(BzzVector &x, double t, BzzVector &f);
	void DAESystem_Premixed_SOOT_NoEnergy(BzzVector &x, double t, BzzVector &f);
	
	// ODE Systems - Interfaces
	// ---------------------------------------------------------------------------
	void solveODE_Premixed(const flame1d_model kind, double tEnd);

	void ODESystem_Premixed_NoEnergy(BzzVector &x, double t, BzzVector &f);
	void ODESystem_Premixed_All(BzzVector &x, double t, BzzVector &f);
	void ODESystem_SingleReactor_Isothermal(BzzVector &x, double t, BzzVector &f);

	// Print Functions
	// ---------------------------------------------------------------------------
	void printOnFile(const std::string fileNameOutput);
	void PrintXMLFile(const std::string file_name);
	void printBackUpOnlyInputData(const std::string fileName);
	void printBackUpOnlyData(const std::string fileName);
	void DAE_ODE_myPrint(BzzVector &v, double time);


	// BackUp Functions
	// ---------------------------------------------------------------------------
	void recoverFromBackUp(const std::string fileName);


	// Refine Grid
	// ---------------------------------------------------------------------------
	void doubleTheGrid();
	void refineGridPeak(double fraction);
	void refineFlameBase();
	void refineFlameBase(const double xA, const double xB);
	bool newPoints(const std::string string_kind, char index);
	void adaptGrid(const int iOption);
	void refineBase();
	void refineAttachPoint();
	void refineStagnationPlane(const double fraction);

	std::string nameOutputFolder; 
	std::string nameFolderBackupData;
	std::string nameFileBackupInputData;
	std::string nameFolderUnsteadyData;
	std::string nameFolderSteadyData;
	std::string nameFolderAdditionalData;
	
	void unsteady_boundary_conditions(double &time);
	OpenSMOKE_Flame1D_OscillatingBoundary unsteady;
	int iUnsteadyFromBackUp;

	OpenSMOKE_GlobalKinetics *global;

	OpenSMOKE_2EModel soot2EModel;
	
	BzzVector phiN;
	BzzVector phiM;
	BzzVector dphiN;
	BzzVector dphiM;
	BzzVector source_phiN;
	BzzVector source_phiM;
	BzzVector diff_phiN;
	BzzVector diff_phiM;
	double phiNC, phiMC;
	double phiNO, phiMO;
	BzzMatrix SootGasCorrection;
	BzzVector DiffusionSoot;
	
	void give_Premixed_SOOT();
	void give_Opposed_SOOT();
	void give_Opposed_QMOM();
	void give_Twin_SOOT();
	void give_Twin_QMOM();
	void allocate_SOOT_Np();
	void initial_conditions_soot_module();
	
	BzzVector BCW_C, BCW_O;

private:

	// Memory allocation
	// ---------------------------------------------------------------------
	void allocate_only_Np();
	void allocate_main_variables_only_Np();
	void allocate_all();

	// Set ODE-DAE-NLS parameters
	// ---------------------------------------------------------------------
	void setMinimumAndMaximumValues(const  flame1d_model string_kind);
	void setMinimumAndMaximumValuesReduced(const  flame1d_model string_kind);
	void setDifferentialAndAlgebraic(const  flame1d_model string_kind);
	
	// Set Initial Conditions
	// ---------------------------------------------------------------------
	void setTemperatureAndMassFractionProfiles(const flame1d_ignition string_kind);

	// Properties of gas mixture
	// ---------------------------------------------------------------------
	void properties(int ReactionON, int jacobianIndex, BzzVectorInt &jacobianVariables, int dimBlock, int indexT);	
	void properties(int ReactionON);
	int indexReactor;
	void propertiesSingleReactors(const int i);	
	void molarFractionsAndPMtot();
	void massFractionsAndPMtot();
	void ElementalAnalysis();
	
	// Diffusional and correction velocities
	// ---------------------------------------------------------------------
	void compute_vStarOpposed();
	void compute_vStarPremixed();

	// Derivatives
	// -------------------------------------------------------------------------------
	BzzVector diffM, diffT, diffTcentral, diffWscalar;
	BzzMatrix diffW;
	double nGeometry;


	// NLS and DAE variables
	// -------------------------------------------------------------------------------
	BzzVectorInt inDerAlg;
	BzzVector xMin, xMax;		
	BzzVector xFirstGuess, fNLS;

	// Print Variables
	// -------------------------------------------------------------------------------	
	ofstream fOutput, fUnsteady, fUnsteadyMax, fUnsteadyQMOMMax;
	ifstream fInput;

	// Service Variables
	// -------------------------------------------------------------------------------	
	BzzVector xVector, wVector, cVector;
	BzzVector DmixVector, RVector, TetaMixVector;
	BzzVector auxNp;

	double Tmax;

	// Residuals
	// -------------------------------------------------------------------------------	
	void give_Opposed_DU_DG_DH();
	void give_Opposed_DW(const std::string string_kind);
	void give_Opposed_DW(const int i);
	void give_Opposed_DT(double t);

	void give_Twin_DU_DG_DH();
	void give_Twin_DW(const std::string string_kind);
	void give_Twin_DT(double t);

	void give_Premixed_DW(const std::string string_kind);
	void give_Premixed_DT();

	
	// BackUp Functions
	// ---------------------------------------------------------------------------
	void recoverVariables(const std::string fileName);

	// Refining Grid
	// ---------------------------------------------------------------------------
	void AddPoints(BzzVectorInt &listPoints);

	// Boundary Conditions
	// -------------------------------------------------------------------------------
	double UO, UC, GO, GC, xST, K;
	BzzVector WC, WO;

	// Auxiliary variables
	// -------------------------------------------------------------------------------
	BzzVector M, sumW, vc;
	BzzMatrix vStar;
	BzzMatrix vThermophoretic;
	BzzMatrix vFick;

	BzzMatrix vStar_e, vStar_w, coeff_e;
	BzzVector vc_e;
	BzzMatrix W_c, X_c;
	BzzVector MW_c;

	double soot_deposition;

	// Gas Properties
	// -------------------------------------------------------------------------------
	BzzVector urho, mu, lambda, Cp, QReaction, PMtot, uPMtot, sumCpDiffusive;
	BzzMatrix Dm, R, Teta, RR;
	BzzVector rhow, rhoe;
	BzzVector A_x_lambdaw, A_x_lambdae;
	BzzVector A_x_rhow, A_x_rhoe;
	BzzMatrix Cpk, CpMap, lambdaMap, muMap, *DjkMap, *TetakjMap;
	BzzMatrix k1Map, k2Map, uKeqMap;
	BzzMatrix logFcentMap, reactionDHMap, reactionDSMap;

	// Initialize
	// -------------------------------------------------------------------------------
	void setupGasMixture();
	void setupBoundaryConditions();

	void updatingProfiles();

	void recoverPhysicalVariables(const flame1d_model string_kind, BzzVector &x);
	void recoverResiduals(const flame1d_model string_kind, BzzVector &f);
	void recoverPhysicalVariables_Reduced(const flame1d_model string_kind, BzzVector &x);
	void recoverResiduals_Reduced(const flame1d_model string_kind, BzzVector &f);

	int  nonLinearSystemSolution(BzzNonLinearSystemSparseObject &o, int dimensionBlock);
	int  nonLinearSystemSolution(BzzNonLinearSystemObject &o);
//	void DAESystemSolution(BzzDaeDoubleSparseObject &o, double tEnd);
	void DAESystemSolution(BzzDaeSparseObject *o, double tEnd);
	void ODESystemSolution(BzzOdeSparseStiffObject &o, double tEnd);

//	void nonLinearSystemSolutionWithControl(const std::string string_kind, BzzNonLinearSystemSparseObject &o, int dimensionBlock);
	void nonLinearSystemSolutionWithControl_FlameSpeed(const double t);
	void nonLinearSystemSolutionWithControl_Opposed(const double time_first);
	void nonLinearSystemSolutionWithControl_Twin(const double time_first);

	int tagConstantT;

	void GnuPlotInterface(ofstream &fOutput);

	void GnuPlotInterfaceUnsteady();
//	void GnuPlotPAHInterface(ofstream &fPAH);
	void GnuPlotSootInterface(ofstream &fSoot);
	void GnuPlotSootDistributionInterface(ofstream &fSoot);

	void GnuPlotLewisInterface(ofstream &fLewis);
	void GnuPlotSoretInterface(ofstream &fSoot);
	void GnuPlotFormationRatesInterface(ofstream &fFormationRates);
	void GnuPlotReactionRatesInterface(ofstream &fReactionRates);
	void GnuPlotSingleContributionsInterface(ofstream &fSingleContributions);
	void GnuPlotMixturePropertiesInterface(ofstream &fProperties);
	void printSummaryForOpposed();
	void printSummaryForPremixed();
	void printSummaryForTwin();

	void nonLinearSystemLabel(BzzNonLinearSystemSparseObject &o, char &control);
	void nonLinearSystemLabel(BzzNonLinearSystemObject &o, char &control);

	void give_Premixed_DH();

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
	
	// QMOM
	OpenSMOKE_QMOM_Module qmom;
	ofstream fSootQMOM;

	BzzMatrix moments;
	BzzMatrix dmoments;
	BzzVector qmom_sources;
	BzzMatrix moments_source;
	BzzMatrix momentsInitial;

	BzzMatrix single_contributions;
	
	BzzMatrix diffMoments;
	BzzMatrix DiffusionMoments;
	BzzMatrix A_x_rho_x_Diffw;
	BzzMatrix A_x_rho_x_Diffe;
	BzzVector momentsC;

	void assign_qmom_module_from_file(OpenSMOKE_ReactingGas &mix, const std::string fileName);
	void allocate_QMOM_Np();
	void allocate_QMOM_N();
	void initial_conditions_qmom_module();
	void give_Premixed_QMOM();

	
	// Radiation
	// ---------------------------------------------------------------------------
	void prepare_radiation();
	void calculate_radiation();
	BzzVector Qrad;
	double sigmaSB , Tenv4;
	int iH2O, iCO2, iCO, iCH4;

	// Sensitivity
	kindOfSensitivityParameter kindOfSensitivity;
	void final_sensitivity_analysis(std::string kindOfProblem,  BzzNonLinearSystemSparseObject &nls, BzzSave &fBinary);
	void final_sensitivity_analysis_transport_properties(std::string kindOfProblem,  BzzNonLinearSystemSparseObject &nls,
									BzzVector &xSolution,  BzzVector &fSolution, BzzSave &fBinary);

	// Reaction Flow Analysis
	OpenSMOKE_RateOfProductionAnalysis ropa;

	// Elemet Flux Analysis
	BzzMatrix rForward;
	BzzMatrix rBackward;
	BzzVector rForwardVector;
	BzzVector rBackwardVector;

	// Save on binary file
	void SaveOnBinaryFile(BzzSave &fOutput);
	void SaveOnBinaryFile(const std::string filename, const bool iSensitivity, const flame1d_model NLS_KIND);

	void SolveODE_SingleReactors(const double tEnd);
};

class OpenSMOKE_Flame1D_Solution
{
public:
	
	BzzVector x;
	BzzVector U;
	BzzVector G;
	BzzVector H;
	BzzVector T;
	BzzMatrix W;
	BzzVector phiN;
	BzzVector phiM;
	BzzMatrix moments;
	BzzVector coordinates;

	void PasteFromExternalSolution(OpenSMOKE_Flame1D &flame, OpenSMOKE_Flame1D_DataManager &data_manager);
	

	double VC, VO, TC, TO;
	BzzVector xO, xC;
	BzzVector XC, XO;
	int Np;
	double MassFlowRate;
	double CrossSectionalArea;
	int iFixedTemperature;
	double fixedTemperature;

};

#endif // !defined(FLAME_OPPOSED_PREMIXED)



