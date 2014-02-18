/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_DROPLETMICROGRAVITY_DATAMANAGER)
#define OPENSMOKE_DROPLETMICROGRAVITY_DATAMANAGER

#include "OpenSMOKE.hpp"
#include "liquid/OpenSMOKE_LiquidProperties_Database.h"
#include "liquid/OpenSMOKE_LiquidSpecies.h"

class OpenSMOKE_DropletMicrogravity;

enum 	DropletBoundaryConditions { BCS_NEUMANN, BCS_DIRICHLET };
enum    DropletInterfaceTemperature { INTERFACE_USER_DEFINED, INTERFACE_EQUILIBRIUM };
enum	kind_of_system { DROPLET_NONE, EIGENVALUE, UNSTEADY_BATCH };
enum    DropletLiquidPhaseModel {LIQUID_DROPLET_PERFECTLY_STIRRED, LIQUID_DROPLET_PDE_ONLYT, LIQUID_DROPLET_PDE_NOMOMENTUM, LIQUID_DROPLET_PDE_ALL};

class OpenSMOKE_Dictionary_DropletMicrogravity : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_DropletMicrogravity();
};

class OpenSMOKE_DropletMicrogravity_DataManager  
{
	friend class OpenSMOKE_DropletMicrogravity;

public:

	OpenSMOKE_LiquidProperties_Database *liquid_database;
	OpenSMOKE_LiquidSpecies *liquid_species;

	OpenSMOKE_DropletMicrogravity_DataManager();
	void SetName(const string name);

	void Assign(OpenSMOKE_DropletMicrogravity *_droplet);
	void Assign(OpenSMOKE_ReactingGas *_mix);
	void Assign(OpenSMOKE_LiquidProperties_Database *_liquid_database);
	void ReadFromFile(const string fileName);
	void SetOutputFolderName(const string _name);

	void DefineFromFile(const string inputFile);
	void Lock();

	// Pressure
	double P_atm;
	double P_bar;
	double P_Pascal;
	
	// Boundary Conditions
	double diameterDroplet0;
	double massDroplet0;
	double ratioRadii;
	double stretchingFactor;
	double enhancingFactor;
	double TDroplet0;
	double TEnvironment;
	double rEnvironment;
	double tEnd;

	// Information
	BzzVector OmegaDroplet;
	BzzVector XDroplet;
	BzzVector OmegaEnvironment;
	vector <string> nameOutput;
	string twoEquation_file_name;

	// Flame
	OpenSMOKE_DropletMicrogravity	*droplet;
	OpenSMOKE_ReactingGas			*mix;
	int N;

	// Print Information
	int iteration;
	int iterationVideoCounter;
	int iterationFileCounter;
	int iterationUnsteadyFileCounter;
	int iterationBackUpCounter ;
	int nStepsVideo;
	int nStepsFile;
	int nStepsBackUp;

	// Derivatives
	char iDerT;
	char iDerW;
	char iDerU;

	// Grid refining
	int		nDiff;
	int		nGrad;
	double	deltaDiff;
	double	deltaGrad;

	// User Defined Profiles
	double	Tpeak;

	// Miscellanea
	bool	iBackUp;
	bool	iGlobalKinetics;
	bool	iVelocity;

	BzzVector x_elemental_C;
	BzzVector omega_elemental_C;
	BzzVector x_elemental_O;
	BzzVector omega_elemental_O;

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

	// Additional variables
	bool i2E;

	int iRadiation;
	double absTolerances;
	double relTolerances;

	vector<string> nameFuels;
	BzzVectorInt jFuel;
	int jOxidizer;
	int jInert;
	int jH2O;
	int jCO2;
	int jCO;
	int jCH4;

	bool iSpark;
	vector<double> SparkRatio;
	bool iBackupFromBinaryFile;
	string nameBackupFile;
	DropletBoundaryConditions boundaryConditions;
	DropletInterfaceTemperature dropletInterfaceTemperature;
	double interfaceTemperature;
	vector<string> radiationOptions;
	string nameOutputFolder;

	int NDroplet;
	double dropletStretchingFactor;
	string dropletGridMode;
	DropletLiquidPhaseModel iLiquidPhase;

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

	string name_object;

	void AssignDropletDiameter(const string units, const double value);
	void AssignEnvironmentExtension(const double double_value);
	void AssignStretchingFactor(const double double_value);
	void SetEnhancingFactor(const double double_value);

	void AssignMode(const string _name);
	void AssignVelocity(const string _name);
	void AssignFuel(const vector<string> _names);
	void AssignPressure(const string units, const double value);
	void AssignFlameTemperature(const string units, const double value);
	void AssignDropletTemperature(const string units, const double value);
	void AssignEnvironmentTemperature(const string units, const double value);
	void AssignDropletMassFractions(const vector<string> _names, const vector<double> _values);
	void AssignDropletMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignEnvironmentMassFractions(const vector<string> _names, const vector<double> _values);
	void AssignEnvironmentMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignGridPoints(const int int_value);
	void AssignInterfaceTemperature(const vector<string> _names);
	void AssignLiquidPhase(const string _name);

	void SetVideoSteps(const int int_value);
	void SetFileSteps(const int int_value);
	void SetBackupSteps(const int int_value);
	void SetRelativeTolerance(const double double_value);
	void SetAbsoluteTolerance(const double double_value);
	void SetGasRadiation(const vector<string> _options);
	void SetSootRadiation();
	void SetGridRefineGradient(const int int_value);
	void SetGridRefineCurvature(const int int_value);
	void SetGridRefineGradientStep(const double double_value);
	void SetGridRefineCurvatureStep(const double double_value);
	void SetDerivativeT(const char char_value);
	void SetDerivativeW(const char char_value);
	void SetDerivativeU(const char char_value);
	void SetSensitivityOnFile(const vector<string> _names);
	void SetFormationRatesOnFile(const vector<string> _names);
	void SetROPAOnFile(const vector<string> _names);
	void SetReactionRatesOnFile(const vector<string> _names);
	void UnsetSoretEffect();
	void UnsetReactions();
	void SetIntegrationTime(const string units, const double value);
	void SetSpark(const vector<string> _names);
	void SetBackup(const string _name);
	void SetDropletGrid(const vector<string> string_vector);


	bool iAssignedVideoSteps;
	bool iAssignedFileSteps;
	bool iAssignedBackupSteps;
	bool iAssignedRelativeTolerance;
	bool iAssignedAbsoluteTolerance;
	bool iAssignedGasRadiation;
	bool iAssignedSootRadiation;
	bool iAssignedGridRefineGradient;
	bool iAssignedGridRefineCurvature;
	bool iAssignedGridRefineGradientStep;
	bool iAssignedGridRefineCurvatureStep;
	bool iAssignedSensitivity;
	bool iAssignedFormationRates;
	bool iAssignedROPA;
	bool iAssignedReactionRates;
	bool iSoretEffect;
	bool iReactions;
	kind_of_system iMode;
};

#endif // DATAMANAGER

