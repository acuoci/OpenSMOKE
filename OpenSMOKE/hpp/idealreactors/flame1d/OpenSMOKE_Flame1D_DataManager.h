/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
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

#if !defined(OpenSMOKE_Flame1D_DataManager_H)
#define OpenSMOKE_Flame1D_DataManager_H

#include "OpenSMOKE.hpp"
#include "preprocessing/OpenSMOKE_PreProcessorReactingGas.h"

class OpenSMOKE_ElementFluxAnalysis;
class OpenSMOKE_ElementFluxAnalysisManager;
class OpenSMOKE_Flame1D;
class OpenSMOKE_LiquidSpecies;

class OpenSMOKE_Dictionary_Flame1D : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_Flame1D();
	void PrepareForOpposedFlames();
	void PrepareForPremixedFlames();
};

enum kind_of_pool_fire {POOL_FIRE_NONE, POOL_FIRE_EQUILIBRIUM, POOL_FIRE_LIQUIDPOOL, POOL_FIRE_TASSIGNED};
enum radiative_soot_model {RADIATIVE_SOOT_MODEL_NONE, RADIATIVE_SOOT_MODEL_FLUENT, RADIATIVE_SOOT_MODEL_BILGER, RADIATIVE_SOOT_MODEL_WIDMANN};

class OpenSMOKE_Flame1D_DataManager
{
public:

	OpenSMOKE_Flame1D_DataManager();
	void SetName(const string name);
	void Assign(OpenSMOKE_ReactingGas *_mix);
	void Assign(OpenSMOKE_Flame1D *_flame);
	void Setup(const string kind);

	flame1d_physics kind_of_flame;
	string flameSpeedAnalysisFileName;	
	string opposedFlameAnalysisFileName;	

	bool iGlobalKinetics;
	bool i2E;
	bool iQMOM;
	bool iUnsteady;
	bool iBackUp;
	bool iFlameSpeedAnalysis;			
	bool iOpposedFlameAnalysis;
	bool iSoretEffect;
	bool iLewisNumbers;
	bool iUnityLewisNumbers;
	bool iTurbulentDiffusivity;	
	bool iCorrectionReactionRates;
	bool iSingleContributions;
	bool iCorrectionDiffusivity;
	bool iCorrectionFormationEnthalpy;
	bool iVerboseMixtureProperties;
	bool iUserDefinedLewisNumbers;
	bool iDepositionWall;

	void readFromFileForOpposed(const string fileName);
	void readFromFileForPremixed(const string fileName);
	void printFileForOpposed(const string fileName,  double vFuel, double vAir);
	void printFileForPremixed(const string fileName, double MassFlowRate_Updated);

	void SetOutputFolder(const string _outputFolderName);
	void DefineFromFileOpposed(const string inputFile);
	void DefineFromFilePremixed(const string inputFile);

	void SetDefaultValues();

	// Pressure
	// -------------------------------------------------------------------------------	
	double P_atm, P_bar, P_Pascal;

	// Boundary Conditions
	// -------------------------------------------------------------------------------	
	double VC, VO, L, TC, TO, xO2;

	// Information
	// -------------------------------------------------------------------------------	
	BzzVector xO, xC;
	BzzVector XC, XO, xPeaks;
	BzzVectorInt iX, iXC, iXO, iOut, iPeaks;
	vector<string> nameO;
	vector<string> nameC;
	vector<string> nameOutput;
	vector<string> namePeaks;

	// Derivates
	// -------------------------------------------------------------------------------	
	char iDerG;
	char iDerT;
	char iDerW;
	char iHOT;

	// Flame
	// -------------------------------------------------------------------------------	
	OpenSMOKE_Flame1D *flame;
	OpenSMOKE_ReactingGas *mix;
	int Np;

	// Print Information
	// -------------------------------------------------------------------------------	
	int iteration, iterationVideoCounter, iterationFileCounter, iterationBackUpCounter, 
		           nStepsVideo, nStepsFile, nStepsBackUp;

	// Information (2)
	// -------------------------------------------------------------------------------	
	int jFUEL, jO2, jINERT;
	double nuINERT;
	string nameFuel;
	string nameOxidizer;
	string nameInert;

	int nDiff, nGrad;
	double deltaDiff, deltaGrad;
	double xcen , wmix;
	string geometry;
	string gridKind;
	double alfa;

	string userDefinedFolderName;
	bool iUserDefinedFolderName;

	// User Defined Profiles
	// -------------------------------------------------------------------------------	
	double Tpeak;
	BzzVector peak;

	flame1d_subphysics  kind_of_subphysics;
	adaptive_grid_model kind_of_adaptive_grid;

	double radialGradientC, radialGradientO;

	double MassFlowRate, CrossSectionalArea;
	int iTemperatureProfile, iAreaProfile;

	OpenSMOKE_UD_Profile	ud_temperature_profile;
	OpenSMOKE_UD_Profile	ud_cross_section_profile;
	string ud_temperature_profile_file_name;
	string ud_cross_section_profile_file_name;

	int nCold;

	int iFixedTemperature;
	double fixedTemperature;

	bool iGasRadiation;
	double environmentTemperature;

	double rel_nls_Tolerances;
	double abs_nls_Tolerances;

	double rel_dae_Tolerances;
	double abs_dae_Tolerances;

	// Sensitivity Analisys
	BzzVectorInt	index_sensitivity;
	vector<string>	names_sensitivity;

	// Reaction Rates Analisys
	BzzVectorInt	index_reaction_rates;
	vector<string>	names_reaction_rates;

	// Formation Rates Analisys
	BzzVectorInt	index_formation_rates;
	vector<string>	names_formation_rates;

	// ROPA
	BzzVectorInt	index_ROPA;
	vector<string>	names_ROPA;

	// Experiment
	BzzVectorInt	index_Experiment;
	vector<string>	names_Experiment;

	// Single Contributions
	BzzVectorInt	index_SingleContributions;
	vector<string>	names_SingleContributions;

	bool iAssignedSensitivity;
	bool iAssignedFormationRates;
	bool iAssignedROPA;
	bool iVerboseAssignedROPA;
	bool iAssignedReactionRates;
	bool iFocusReactionRates;
	bool iFocusElemetFluxAnalysis;
	bool iAssignedAdaptiveGridCoefficients;
	bool iAssignedExperiment;
	bool iAssignedGlobalElementFluxAnalysis;
	bool iFakeTemperatureThermalConductivity;
	
	fittingCoefficientSensitivityMode sensitivityLennardJonesMode;

	double adaptive_grid_alfa;
	double adaptive_grid_beta;
	double adaptive_grid_gamma;
	double fakeTemperatureThermalConductivityIncrement;

	string	unsteady_flame_file_name;
	string	qmom_file_name;
	string  twoEquation_file_name;

	void AssignFuelMassFractions(BzzVector &omega_fuel);


	double correctionReactionRates;
	BzzVectorInt	index_correction_diffusivity;
	BzzVector		correction_diffusivity;
	BzzVectorInt	index_correction_formation_enthalpy;
	BzzVector		correction_formation_enthalpy;
	double diffusivityEnhancingFactor;

	OpenSMOKE_ElementFluxAnalysisManager *element_flux_analysis_manager;
	OpenSMOKE_ElementFluxAnalysis *element_flux_analysis;

	bool iAssignedFixedTemperatureProfileProvisional;
	BzzVector user_defined_lewis_numbers;

	double initial_time_step;
	int max_integration_order;

	unsigned int bin_index_zero;
	double bin_density_A;
	unsigned int bin_index_final;
	double bin_density_B;

private:

	ifstream fInput;
	ofstream fOutput;

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);

	void CheckDictionary(OpenSMOKE_Dictionary_Flame1D &dictionary);
	void CheckDictionaryForOpposedFlames(OpenSMOKE_Dictionary_Flame1D &dictionary);
	void CheckDictionaryForPremixedFlames(OpenSMOKE_Dictionary_Flame1D &dictionary);
	void CopyInputFile(const string inputFile);
	void PasteInputFile(const string fileName);

	void AssignPressure(const string string_value, const double double_value);
    void AssignDistance(const string string_value, const double double_value);
    void AssignFuelVelocity(const string string_value, const double double_value);
    void AssignOxidizerVelocity(const string string_value, const double double_value);
    void AssignFuelTemperature(const string string_value, const double double_value);
    void AssignOxidizerTemperature(const string string_value, const double double_value);
    void AssignTemperaturePeak(const string string_value, const double double_value);
    
    void AssignGeometry(const string string_value);
    void AssignFuel(const string string_value);
    void AssignOxidizer(const string string_value);
    void AssignInert(const string string_value);
	void AssignOutputSpecies(const vector<string> string_vector);
	void AssignGridPoints(const int int_value);
    void AssignFuelMassFractions(const vector<string> string_vector, const vector<double> double_vector);
    void AssignFuelMoleFractions(const vector<string> string_vector, const vector<double> double_vector);
    void AssignOxidizerMassFractions(const vector<string> string_vector, const vector<double> double_vector);
    void AssignOxidizerMoleFractions(const vector<string> string_vector, const vector<double> double_vector);

	void SetGrid(const string string_value);
	void SetFlameThickness(const string string_value, const double double_value);
    void SetFlamePosition(const string string_value, const double double_value);
	void SetPeaks(const vector<string> string_vector, const vector<double> double_vector);
    void SetInitialTemperatureProfile(const string string_value);
    void SetFixedTemperatureProfile(const string string_value);
    void SetFixedTemperatureProfileProvisional(const string string_value);
	void SetVideoSteps(const int int_value);
	void SetFileSteps(const int int_value);
	void SetBackupSteps(const int int_value);

	void SetDaeRelativeTolerance(const double double_value);
	void SetDaeAbsoluteTolerance(const double double_value);
	void SetNlsRelativeTolerance(const double double_value);
	void SetNlsAbsoluteTolerance(const double double_value);

	void SetFuelRadialGradient(const string units, const double value);
	void SetOxidizerRadialGradient(const string units, const double value);
	void SetEnvironmentTemperature(const string units, const double value);
	void SetGasRadiation();
	void SetSootRadiation(const string value);
	void SetSoretEffect();
	void SetThermophoreticEffect();
	void SetLewisNumbers();
	void SetUnityLewisNumbers();
	void SetUserDefinedLewisNumbers(const vector<string> names, const vector<double> values);
	void SetCorrectionReactionRates(const double double_value);
	void SetGridRefineGradient(const int int_value);
	void SetGridRefineCurvature(const int int_value);
	void SetGridRefineGradientStep(const double double_value);
	void SetGridRefineCurvatureStep(const double double_value);
	void SetDerivativeG(const char char_value);
	void SetDerivativeT(const char char_value);
	void SetDerivativeW(const char char_value);
	void SetSensitivityOnFile(const vector<string> _names);
	void SetFormationRatesOnFile(const vector<string> _names);
	void SetVerboseROPAOnFile(const vector<string> _names);
	void SetROPAOnFile();
	void SetReactionRatesOnFile(const vector<string> _names);
	void SetBINDensities(const vector<string> _names);
	void SetSingleContributions(const vector<string> _names);
	void SetAdaptiveGrid(const string name);
	void SetAdaptiveGridCoefficients(const vector<string> values);
	void SetTurbulentDiffusivity(const double value);
	void SetExperimentOnFile(const vector<string> _names);
	void SetChangeFrequencyFactor(const vector<string> _names);
	void SetChangeActivationEnergy(const vector<string> _names);
	void SetChangeDiffusivity(const vector<string> _names);
	void SetChangeFormationEnthalpy(const vector<string> _names);
	void SetLennardJonesMode(const string _option);
	void SetElementFluxAnalysisOnFile(const vector<string> _names);
	void SetGlobalElementFluxAnalysisOnFile(const vector<string> _names);
	void SetChangeFakeTemperatureThermalConductivity(const double value, const string units);

	void SetRobustTemperature();
	void SetVerboseMixtureProperties();

	void SetMaximumIntegrationOrder(const int int_value);
	void SetInitialTimeStep(const string units, const double value);

	void SetDepositionWall();

	// Specific for premixed flames
	void AssignMassFlowRate(const string units, const double value);
	void AssignInletVelocity(const string units, const double value);
	void AssignInletTemperature(const string units, const double value);
	void AssignCrossSection(const string units, const double value);
	void AssignInletMassFractions(const vector<string> names, const vector<double> values);
	void AssignInletMoleFractions(const vector<string> names, const vector<double> values);
	
	void SetOutletTemperature(const string units, const double value);
	void SetOutletMassFractions(const vector<string> names, const vector<double> values);
	void SetOutletMoleFractions(const vector<string> names, const vector<double> values);
	void SetFlameSpeedIndex(const int int_value);
	void SetFlameSpeedTemperature(const string units, const double value);
	void SetStretchingFactor(const double double_value);
	void SetCrossSectionProfile(const string string_value);
	void SetCrossSectionMIT();

	void AssignEquivalenceRatioForPremixedFlames(const double value);
	void AssignFuelMassFractionsForPremixedFlames(const vector<string> names, const vector<double> values);
	void AssignFuelMoleFractionsForPremixedFlames(const vector<string> names, const vector<double> values);
	void AssignOxidizerMassFractionsForPremixedFlames(const vector<string> names, const vector<double> values);
	void AssignOxidizerMoleFractionsForPremixedFlames(const vector<string> names, const vector<double> values);
	void AssignInletStreamCompositionForPremixedFlames();

	double	equivalence_ratio;

	bool iAssignedFlameThickness;
	bool iAssignedFlamePosition;
	bool iAssignedGrid;
	bool iAssignedPeaks;
	bool iAssignedInitialTemperatureProfile;
	bool iAssignedFixedTemperatureProfile;
	bool iAssignedVideoSteps;
	bool iAssignedFileSteps;
	bool iAssignedBackupSteps;

	bool iAssignedDaeRelativeTolerance;
	bool iAssignedDaeAbsoluteTolerance;
	bool iAssignedNlsRelativeTolerance;
	bool iAssignedNlsAbsoluteTolerance;

	bool iAssignedFuelRadialGradient;
	bool iAssignedOxidizerRadialGradient;
	bool iAssignedEnvironmentTemperature;
	bool iAssignedGasRadiation;
	bool iAssignedSootRadiation;
	bool iAssignedGridRefineGradient;
	bool iAssignedGridRefineCurvature;
	bool iAssignedGridRefineGradientStep;
	bool iAssignedGridRefineCurvatureStep;
	bool iAssignedDerivativeG;
	bool iAssignedDerivativeT;
	bool iAssignedDerivativeW;

	bool iAssignedOutletTemperature;
	bool iAssignedOutletMassFractions;
	bool iAssignedOutletMoleFractions;
	bool iAssignedFlameSpeedIndex;
	bool iAssignedFlameSpeedTemperature;
	bool iAssignedStretchingFactor;
	bool iAssignedCrossSectionProfile;
	bool iAssignedCrossSectionMIT;

	void LockForOpposedFlames();
	void LockForPremixedFlames();

	vector<string> input_file_lines;
	int indexLinePoints;
	int indexLineFuelVelocity;
	int indexLineOxidizerVelocity;
	int indexLineEquivalenceRatio;
	int indexLineFlameSpeedIndex;



public:

	bool	iEquivalenceRatioForPremixedFlames;
	vector<string> fuel_names;
	BzzVector moles_fuel;
	vector<string> oxidizer_names;
	BzzVector moles_oxidizer;

public:

	OpenSMOKE_LiquidSpecies *pool_fire_liquid_species;
	void SetPoolFire(const vector<string> values);
	void SetPoolFireGridOptions(const vector<string> values);
	void SetPoolFireCorrectionFactorVaporizationHeat(const double value);
	void SetPoolFireCorrectionFactorVaporPressure(const double value);
	void SetPoolFireCorrectionFactorThermalConductivity(const double value);
	void SetPoolFireCorrectionFactorSpecificHeat(const double value);

	kind_of_pool_fire iPoolFire;
	double pool_fire_temperature;
	double pool_fire_feed_temperature;
	double pool_fire_depth;

	radiative_soot_model iRadiativeSootModel;
	bool iThermophoreticEffect;

	double poolfire_grid_alfa_fuel; 
	double poolfire_grid_alfa_oxidizer; 
	double poolfire_grid_point_fraction; 
	double poolfire_grid_distance_fraction; 
	double correctionFactorVaporizationHeat;
	double correctionFactorVaporPressure;
	double correctionFactorSpecificHeat;
	double correctionFactorThermalConductivity;

public:

	int counterUnchanged;
	double TMaxOld;
	bool iRobustTemperature;
};

#endif // !defined(OpenSMOKE_Flame1D_DataManager_H)
