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

#include <sstream>
#include "droplet/OpenSMOKE_DropletMicrogravity_DataManager.h"
#include "droplet/OpenSMOKE_DropletMicrogravity.h"

void OpenSMOKE_DropletMicrogravity_DataManager::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_DropletMicrogravity_DataManager"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_DropletMicrogravity_DataManager::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_DropletMicrogravity_DataManager"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_DropletMicrogravity_DataManager::OpenSMOKE_DropletMicrogravity_DataManager()
{
	name_object			= "[Name not assigned]";
	iBackUp				= false;
	iGlobalKinetics		= false;
	i2E					= false;

	iAssignedVideoSteps = false;
	iAssignedFileSteps = false;
	iAssignedBackupSteps = false;
	iAssignedRelativeTolerance = false;
	iAssignedAbsoluteTolerance = false;
	iAssignedGasRadiation = false;
	iAssignedSootRadiation = false;
	iAssignedGridRefineGradient = false;
	iAssignedGridRefineCurvature = false;
	iAssignedGridRefineGradientStep = false;
	iAssignedGridRefineCurvatureStep = false;
	iAssignedSensitivity = false;
	iAssignedFormationRates = false;
	iAssignedROPA = false;
	iAssignedReactionRates = false;
	iSoretEffect = true;
	iReactions = true;
	iSpark = false;
	iBackupFromBinaryFile = false;
	iMode = DROPLET_NONE;
	iVelocity = false;
	
	// Grid refinement
	nDiff = 4;
	nGrad = 2;	
	deltaDiff = 0.05;
	deltaGrad = 0.05;

	// Output
	nStepsVideo  = 10;		
	nStepsFile	 = 1000;
	nStepsBackUp = 2000;
	tEnd		 = 1.e8;

	// Radiation
	iRadiation = 0;
	
	// Tolerances
	relTolerances = 100.*MachEpsFloat();			
	absTolerances = 1.e-10;	

	// Derivatives
	iDerT = 'C';
	iDerW = 'C';
	iDerU = 'C';

	// Vaporization rate
	SparkRatio.resize(4);
	SparkRatio[0] = 1.2;
	SparkRatio[1] = 2.;
	SparkRatio[2] = 10.;
	SparkRatio[3] = 15.;
	boundaryConditions = BCS_DIRICHLET;
	dropletInterfaceTemperature = INTERFACE_USER_DEFINED;

	nameOutputFolder = "Output";

	// Liquid Droplet
	dropletGridMode = "fixedRatio";
	dropletStretchingFactor = 1.;
	NDroplet = 20;
	iLiquidPhase = LIQUID_DROPLET_PERFECTLY_STIRRED;

	enhancingFactor = 0;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetName(const string name)
{
	name_object = name;
}

void OpenSMOKE_DropletMicrogravity_DataManager::Assign(OpenSMOKE_ReactingGas *_mix)
{
	mix = _mix;
}

void OpenSMOKE_DropletMicrogravity_DataManager::Assign(OpenSMOKE_DropletMicrogravity *_droplet)
{
	droplet = _droplet;
}

void OpenSMOKE_DropletMicrogravity_DataManager::Assign(OpenSMOKE_LiquidProperties_Database *_liquid_database)
{
	liquid_database = _liquid_database;
}

void OpenSMOKE_DropletMicrogravity_DataManager::ReadFromFile(const string fileName)
{
	DefineFromFile(fileName);

	iteration				= 1;
	iterationVideoCounter	= nStepsVideo;
	iterationFileCounter	= nStepsFile;
	iterationBackUpCounter	= nStepsBackUp;
	iterationUnsteadyFileCounter	= 50;
/*
	StoichiometricMixtureFraction();

	if (iAssignedFlamePosition == false)	xcen = 1.-zStoichiometric;

	// Enthalpy
	{
		double pmc = mix->GetMWFromMoleFractions(XC);
		enthalpyC  = mix->GetMixEnthalpy_Mole(TC, XC)/pmc;		// [J/kg]

		double pmo = mix->GetMWFromMoleFractions(XO);
		enthalpyO  = mix->GetMixEnthalpy_Mole(TO, XO)/pmo;		// [J/kg]

		double tc = mix->GetTemperatureFromMassEnthalpyAndMoleFractions(300., enthalpyC, XC);
		double to = mix->GetTemperatureFromMassEnthalpyAndMoleFractions(300., enthalpyO, XO);

		cout << endl;
		cout << "Fuel Enthalpy        " << enthalpyC << " J/kg" << endl;
		cout << "Oxidizer Enthalpy    " << enthalpyO << " J/kg" << endl;
		cout << "Fuel Temperature     " << tc << " K" << endl;
		cout << "Oxidizer Temperature " << to << " K" << endl;
		cout << endl;
	}

	// Element composition
	{
		int j;

		ChangeDimensions(mix->NumberOfElements(), &x_elemental_C);
		ChangeDimensions(mix->NumberOfElements(), &x_elemental_O);
		ChangeDimensions(mix->NumberOfElements(), &omega_elemental_C);
		ChangeDimensions(mix->NumberOfElements(), &omega_elemental_O);

		// Elemental analysis
		mix->GetElementalMoleFractionsFromSpeciesMoleFractions(x_elemental_C, XC);
		mix->GetElementalMassFractionsFromSpeciesMoleFractions(omega_elemental_C, XC);
		mix->GetElementalMoleFractionsFromSpeciesMoleFractions(x_elemental_O, XO);
		mix->GetElementalMassFractionsFromSpeciesMoleFractions(omega_elemental_O, XO);


		cout << "Fuel elemental composition: x and w" << endl;
		for(j=1;j<=mix->NumberOfElements();j++)
			cout << mix->list_of_elements[j-1] << " " << x_elemental_C[j] << " " << omega_elemental_C[j] << endl;
		cout << endl;

		cout << "Oxidizer elemental composition: x and w" << endl;
		for(j=1;j<=mix->NumberOfElements();j++)
			cout << mix->list_of_elements[j-1] << " " << x_elemental_O[j] << " " << omega_elemental_O[j] << endl;
		cout << endl;
	}	
	*/
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//									DICTIONARY													   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_DropletMicrogravity::OpenSMOKE_Dictionary_DropletMicrogravity()
{
    SetupBase();
	SetName("OpenSMOKE_DropletMicrogravity Dictionary");

	Add("#Mode",					    'C', 'S', "Mode: Steady || Unsteady");
	Add("#Velocity",					'O', 'S', "Velocity: on || off");
	Add("#Pressure",					'C', 'M', "Pressure");
	Add("#Fuel",					    'C', 'V', "Fuel names");
	Add("#LiquidPhase",				    'C', 'S', "LiquidPhase: stirred || pde-only-T || pde-no-momentum || pde-all");
	

	Add("#DropletGrid",					'O', 'V', "Droplet grid: Npoints Grid Stretching");

	Add("#DropletDiameter",				'C', 'M', "Droplet diameter");
	Add("#EnvironmentExtension",		'C', 'D', "Environment number of diameters");
	Add("#DropletTemperature",			'C', 'M', "Droplet temperature");
	Add("#EnvironmentTemperature",		'C', 'M', "Environment temperature");

	Add("#DropletMassFractions",		'O', 'L', "Droplet mass fractions");
    Add("#DropletMoleFractions",		'O', 'L', "Droplet mole fractions");
	Add("#EnvironmentMassFractions",    'O', 'L', "Environment mass fractions");
    Add("#EnvironmentMoleFractions",	'O', 'L', "Environment mole fractions");

	Add("#Points",						'C', 'I', "Number of grid points");
	Add("#StretchingFactor",			'C', 'D', "Stretching factor");
	Add("#EnhancingFactor",				'O', 'D', "Enhancing factor (from 1 to 10)");
	
	Add("#FlameTemperature",		'C', 'M', "Flame temperature");
	Add("#Spark",					'O', 'V', "Spark");
	Add("#Backup",					'O', 'S', "File name");
	Add("#InterfaceTemperature",	'C', 'V', "Droplet interface temperature: user defined or equilibrium");

	Add("#GridRefineGradient",      'O', 'I', "Maximum number of new grid points in grid refinement (gradient)");
	Add("#GridRefineCurvature",     'O', 'I', "Maximum number of new grid points in grid refinement (curvature)");
	Add("#GridRefineGradientStep",  'O', 'D', "Grid refinement control parameter (gradient)");
	Add("#GridRefineCurvatureStep", 'O', 'D', "Grid refinement control parameter (curvature)");
 
	Add("#RelativeTolerance",		'O', 'D', "Relative tolerance");
	Add("#AbsoluteTolerance",		'O', 'D', "Absolute tolerance");

	Add("#GasRadiation",			'O', 'V', "Radiative heat transfer from gas phase");
	Add("#SootRadiation",			'O', 'N', "Radiative heat transfer from gas phase and soot");

	Add("#ReactionRates",			'O', 'V', "Reaction rates on file (list of reaction indices)");
	Add("#FormationRates",			'O', 'V', "Formation rates on file (list of species names)");
	Add("#ROPA",					'O', 'V', "Rate of Production Analysis (list of species names)");
	Add("#Sensitivity",				'O', 'V', "Sensitivity Analysis (list of species names)");

	Add("#NoReactions",				'O', 'N', "Unset gas phase reactions");
	Add("#NoSoretEffect",			'O', 'N', "Unset Soret effect");
	Add("#IntegrationTime",			'O', 'M', "Integration time");

	Add("#DerivativeT",             'O', 'H', "Temperature derivative: U || C || B || F");
	Add("#DerivativeW",             'O', 'H', "Mass fractions derivative: U || C || B || F");
	Add("#DerivativeU",             'O', 'H', "Velocity derivative: U || C || B || F");

	Compulsory("#DropletMassFractions",     "#DropletMoleFractions");
	Compulsory("#EnvironmentMassFractions", "#EnvironmentMoleFractions");

    Conflict("#DropletMassFractions",		"#DropletMoleFractions");
    Conflict("#EnvironmentMassFractions",	"#EnvironmentMoleFractions");
	Conflict("#GasRadiation",				"#SootRadiation");

    Lock();
}

void OpenSMOKE_DropletMicrogravity_DataManager::DefineFromFile(const string inputFile)
{
    int     int_value;
	double  double_value;
    string  string_value;
	char    char_value;
	vector<string> string_vector;
	vector<double> double_vector;

	OpenSMOKE_Dictionary_DropletMicrogravity dictionary;
    dictionary.ParseFile(inputFile);

	if (dictionary.Return("#Velocity", string_value))
        AssignVelocity(string_value);

	if (dictionary.Return("#Mode", string_value))
        AssignMode(string_value);

	if (dictionary.Return("#Fuel", string_vector))
        AssignFuel(string_vector);

	if (dictionary.Return("#LiquidPhase", string_value))
        AssignLiquidPhase(string_value);

	if (dictionary.Return("#DropletGrid", string_vector))
        SetDropletGrid(string_vector);

	if (dictionary.Return("#Pressure", double_value, string_value))
        AssignPressure(string_value, double_value);

	if (dictionary.Return("#DropletDiameter", double_value, string_value))
        AssignDropletDiameter(string_value, double_value);

	if (dictionary.Return("#EnvironmentExtension", double_value))
        AssignEnvironmentExtension(double_value);

	if (dictionary.Return("#StretchingFactor", double_value))
        AssignStretchingFactor(double_value);

	if (dictionary.Return("#EnhancingFactor", double_value))
        SetEnhancingFactor(double_value);

	if (dictionary.Return("#IntegrationTime", double_value, string_value))
        SetIntegrationTime(string_value, double_value);

	if (dictionary.Return("#Spark", string_vector))	
		SetSpark(string_vector);
	
	if (dictionary.Return("#Backup", string_value))	
		SetBackup(string_value);

	if (dictionary.Return("#InterfaceTemperature", string_vector))	
		AssignInterfaceTemperature(string_vector);

	if (dictionary.Return("#DropletTemperature", double_value, string_value))
        AssignDropletTemperature(string_value, double_value);

	if (dictionary.Return("#EnvironmentTemperature", double_value, string_value))
        AssignEnvironmentTemperature(string_value, double_value);

    if (dictionary.Return("#DropletMassFractions", double_vector, string_vector))
		AssignDropletMassFractions(string_vector, double_vector);
    if (dictionary.Return("#DropletMoleFractions", double_vector, string_vector))
		AssignDropletMoleFractions(string_vector, double_vector);

    if (dictionary.Return("#EnvironmentMassFractions", double_vector, string_vector))
		AssignEnvironmentMassFractions(string_vector, double_vector);
    if (dictionary.Return("#EnvironmentMoleFractions", double_vector, string_vector))
		AssignEnvironmentMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#Points", int_value))
        AssignGridPoints(int_value);

	if (dictionary.Return("#FlameTemperature", double_value, string_value))
        AssignFlameTemperature(string_value, double_value);

	if (dictionary.Return("#GridRefineGradient", int_value))
		SetGridRefineGradient(int_value);

	if (dictionary.Return("#GridRefineCurvature", int_value))
		SetGridRefineCurvature(int_value);

	if (dictionary.Return("#GridRefineGradientStep", double_value))
		SetGridRefineGradientStep(double_value);

	if (dictionary.Return("#GridRefineCurvatureStep", double_value))
		SetGridRefineCurvatureStep(double_value);

	if (dictionary.Return("#RelativeTolerance", double_value))
		SetRelativeTolerance(double_value);

	if (dictionary.Return("#AbsoluteTolerance", double_value))
		SetAbsoluteTolerance(double_value);

   if (dictionary.Return("#GasRadiation", string_vector))
		SetGasRadiation(string_vector);

	if (dictionary.Return("#SootRadiation"))	
		SetSootRadiation();

    if (dictionary.Return("#ReactionRates", string_vector))
		SetReactionRatesOnFile(string_vector);

    if (dictionary.Return("#FormationRates", string_vector))
		SetFormationRatesOnFile(string_vector);

    if (dictionary.Return("#ROPA", string_vector))
		SetROPAOnFile(string_vector);

    if (dictionary.Return("#Sensitivity", string_vector))
		SetSensitivityOnFile(string_vector); 

	if (dictionary.Return("#Nvideo", int_value))
		SetVideoSteps(int_value);

	if (dictionary.Return("#Nfile", int_value))
		SetFileSteps(int_value);

	if (dictionary.Return("#Nbackup", int_value))
		SetBackupSteps(int_value);

	if (dictionary.Return("#NoReactions"))
		UnsetReactions();

	if (dictionary.Return("#NoSoretEffect"))
		UnsetSoretEffect();

	if (dictionary.Return("#DerivativeT", char_value))
		SetDerivativeT(char_value);

	if (dictionary.Return("#DerivativeW", char_value))
		SetDerivativeW( char_value);

	if (dictionary.Return("#DerivativeU", char_value))
		SetDerivativeU( char_value);

	Lock();
}

void OpenSMOKE_DropletMicrogravity_DataManager::Lock()
{
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetOutputFolderName(const string _name)
{
	nameOutputFolder = _name;
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignMode(const string _name)
{
	if (_name == "Steady")			iMode = EIGENVALUE;	
	else if (_name == "Unsteady")	iMode = UNSTEADY_BATCH;
	else ErrorMessage ("Wrong mode!");
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignVelocity(const string _name)
{
	if (_name == "on")			iVelocity = true;	
	else if (_name == "off")	iVelocity = false;
	else ErrorMessage ("Wrong velocity!");
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignLiquidPhase(const string _name)
{
	if (_name == "stirred")				iLiquidPhase = LIQUID_DROPLET_PERFECTLY_STIRRED;	
	else if (_name == "pde-only-T")		iLiquidPhase = LIQUID_DROPLET_PDE_ONLYT;	
	else if (_name == "pde-nomomentum")	iLiquidPhase = LIQUID_DROPLET_PDE_NOMOMENTUM;	
	else if (_name == "pde-all")		iLiquidPhase = LIQUID_DROPLET_PDE_ALL;	
	else ErrorMessage ("Wrong #LiquidPhase options!");
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetDropletGrid(const vector<string> string_vector)
{
	NDroplet = atoi(string_vector[1].c_str());
	dropletGridMode = string_vector[0];
	dropletStretchingFactor = atof(string_vector[2].c_str());
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignFuel(const vector<string> string_vector)
{
	liquid_species = new OpenSMOKE_LiquidSpecies[string_vector.size()+1];

	nameFuels = string_vector;
	for (int j=1;j<=string_vector.size();j++)
	{
		liquid_species[j].SetName(nameFuels[j-1]);
		liquid_species[j].SetProperties(*liquid_database);
		liquid_species[j].Summary();
	}
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignPressure(const string units, const double value)
{
    P_Pascal = OpenSMOKE_Conversions::conversion_pressure(value, units);
	P_atm	 = P_Pascal / OpenSMOKE_Conversions::Pa_from_atm;
	P_bar	 = P_Pascal / OpenSMOKE_Conversions::Pa_from_bar;
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignFlameTemperature(const string units, const double value)
{
    Tpeak = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignDropletTemperature(const string units, const double value)
{
    TDroplet0 = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignEnvironmentTemperature(const string units, const double value)
{
    TEnvironment = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignDropletMassFractions(const vector<string> _names, const vector<double> _values)
{
	ChangeDimensions(mix->NumberOfSpecies(), &OmegaDroplet);
	for(int i=1;i<=int(_names.size());i++)
		OmegaDroplet[mix->recognize_species(_names[i-1])] = _values[i-1];

	if (OmegaDroplet.GetSumElements() != 1.)
		ErrorMessage("Wrong droplet composition...");

	double MWmix;
	mix->GetMWAndMoleFractionsFromMassFractions(MWmix, XDroplet, OmegaDroplet);
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignDropletMoleFractions(const vector<string> _names, const vector<double> _values)
{
	BzzVector xDroplet(mix->NumberOfSpecies());
	for(int i=1;i<=int(_names.size());i++)
		xDroplet[mix->recognize_species(_names[i-1])] = _values[i-1];

	if (xDroplet.GetSumElements() != 1.)
		ErrorMessage("Wrong droplet composition...");

	double sum=0.;
	for(int i=1;i<=mix->NumberOfSpecies();i++)
		sum += mix->M[i]*xDroplet[i];

	ChangeDimensions(mix->NumberOfSpecies(), &OmegaDroplet);
	for(int i=1;i<=mix->NumberOfSpecies();i++)
		OmegaDroplet[i] = xDroplet[i]*mix->M[i]/sum;
	XDroplet = xDroplet;
}
	
void OpenSMOKE_DropletMicrogravity_DataManager::AssignEnvironmentMassFractions(const vector<string> _names, const vector<double> _values)
{
	ChangeDimensions(mix->NumberOfSpecies(), &OmegaEnvironment);
	for(int i=1;i<=int(_names.size());i++)
		OmegaEnvironment[mix->recognize_species(_names[i-1])] = _values[i-1];

	if (OmegaEnvironment.GetSumElements() != 1.)
		ErrorMessage("Wrong environment composition...");
}
    
void OpenSMOKE_DropletMicrogravity_DataManager::AssignEnvironmentMoleFractions(const vector<string> _names, const vector<double> _values)
{
	BzzVector XEnvironment(mix->NumberOfSpecies());
	for(int i=1;i<=int(_names.size());i++)
		XEnvironment[mix->recognize_species(_names[i-1])] = _values[i-1];

	if (XEnvironment.GetSumElements() != 1.)
		ErrorMessage("Wrong environment composition...");

	double sum=0.;
	for(int i=1;i<=mix->NumberOfSpecies();i++)
		sum += mix->M[i]*XEnvironment[i];

	ChangeDimensions(mix->NumberOfSpecies(), &OmegaEnvironment);
	for(int i=1;i<=mix->NumberOfSpecies();i++)
		OmegaEnvironment[i] = XEnvironment[i]*mix->M[i]/sum;	
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignGridPoints(const int int_value)
{
	N = int_value;
}

//void OpenSMOKE_DropletMicrogravity_DataManager::SetFlamePosition(const double value)
//{
//    flamePosition = value;
//}

void OpenSMOKE_DropletMicrogravity_DataManager::SetVideoSteps(const int int_value)
{
	nStepsVideo = int_value;
	iAssignedVideoSteps = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetFileSteps(const int int_value)
{
	nStepsFile = int_value;
	iAssignedFileSteps = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetBackupSteps(const int int_value)
{
	nStepsBackUp = int_value;
	iAssignedBackupSteps = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetRelativeTolerance(const double double_value)
{
	relTolerances = double_value;
	iAssignedRelativeTolerance = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetAbsoluteTolerance(const double double_value)
{
	absTolerances = double_value;
	iAssignedAbsoluteTolerance = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetGasRadiation(const vector<string> options)
{
	iRadiation = 1;
	iAssignedGasRadiation = true;
	radiationOptions = options;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetSootRadiation()
{
	iRadiation = 2;
	iAssignedSootRadiation = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetGridRefineGradient(const int int_value)
{
	nDiff = int_value;
	iAssignedGridRefineGradient = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetGridRefineCurvature(const int int_value)
{
	nGrad = int_value;
	iAssignedGridRefineCurvature = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetGridRefineGradientStep(const double double_value)
{
	deltaDiff = double_value;
	iAssignedGridRefineGradientStep = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetGridRefineCurvatureStep(const double double_value)
{
	deltaGrad = double_value;
	iAssignedGridRefineCurvatureStep = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetSensitivityOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_sensitivity, names_sensitivity);
	iAssignedSensitivity = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetFormationRatesOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_formation_rates, names_formation_rates);
	iAssignedFormationRates = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignInterfaceTemperature(const vector<string> _names)
{
	if (_names.size() == 2)
	{
		dropletInterfaceTemperature = INTERFACE_USER_DEFINED;
		interfaceTemperature = OpenSMOKE_Conversions::conversion_temperature(atof(_names[0].c_str()), _names[1]);
	}
	else
	{
		dropletInterfaceTemperature = INTERFACE_EQUILIBRIUM;
		interfaceTemperature = 0.;
	}
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetBackup(const string _name)
{
	iBackupFromBinaryFile = true;
	nameBackupFile = _name;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetROPAOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_ROPA, names_ROPA);
	iAssignedROPA = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetReactionRatesOnFile(const vector<string> _names)
{
	if (_names[0] == "ALL")
	{
		ChangeDimensions(mix->NumberOfReactions(), &index_reaction_rates);
		for(int k=1;k<=mix->NumberOfReactions();k++)
		{
			index_reaction_rates[k] = k;
			names_reaction_rates.push_back(mix->strReaction[index_reaction_rates[k]]);
		}
	}
	else
	{
		ChangeDimensions(_names.size(), &index_reaction_rates);
		for(int k=0;k<int(_names.size());k++)
		{
			stringstream index_string;
			index_string << _names[k];
			index_string >> index_reaction_rates[k+1];
			if (index_reaction_rates[k+1] > mix->NumberOfReactions())
				ErrorMessage("The requested reaction for the reaction rates analysis is not available: " + _names[k]);
			names_reaction_rates.push_back(mix->strReaction[index_reaction_rates[k+1]]);
		}
	}
	
	iAssignedReactionRates = true;
}

void OpenSMOKE_DropletMicrogravity_DataManager::AssignDropletDiameter(const string units, const double value)
{
    diameterDroplet0 = OpenSMOKE_Conversions::conversion_length(value, units);
}

void OpenSMOKE_DropletMicrogravity_DataManager:: AssignEnvironmentExtension(const double double_value)
{
	ratioRadii = double_value;
}

void OpenSMOKE_DropletMicrogravity_DataManager:: AssignStretchingFactor(const double double_value)
{
	stretchingFactor = double_value;
}

void OpenSMOKE_DropletMicrogravity_DataManager:: SetEnhancingFactor(const double double_value)
{
	enhancingFactor = double_value;
}

void OpenSMOKE_DropletMicrogravity_DataManager:: SetIntegrationTime(const string units, const double value)
{
	tEnd = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_DropletMicrogravity_DataManager::UnsetSoretEffect()
{
	iSoretEffect = false;
}

void OpenSMOKE_DropletMicrogravity_DataManager::UnsetReactions()
{
	iReactions = false;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetDerivativeT(const char char_value)
{
	if (char_value != 'U' && char_value != 'C' && 
		char_value != 'B' && char_value != 'F')
			ErrorMessage("#DerivativeT options: U || C || B || F");

	iDerT = char_value;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetDerivativeW(const char char_value)
{
	if (char_value != 'U' && char_value != 'C' && 
		char_value != 'B' && char_value != 'F')
			ErrorMessage("#DerivativeW options: U || C || B || F");

	iDerW = char_value;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetDerivativeU(const char char_value)
{
	if (char_value != 'U' && char_value != 'C' && 
		char_value != 'B' && char_value != 'F')
			ErrorMessage("#DerivativeU options: U || C || B || F");

	iDerU = char_value;
}

void OpenSMOKE_DropletMicrogravity_DataManager::SetSpark(const vector<string> _names)
{
	iSpark = true;

	if (_names[0] == "default")
		return;
	else
	{
		SparkRatio[0] = atof(_names[0].c_str());
		SparkRatio[1] = atof(_names[1].c_str());
		SparkRatio[2] = atof(_names[2].c_str());
		SparkRatio[3] = atof(_names[3].c_str());

		//int count = 4;
		//for (int j=1;j<=(_names.size()-4)/2;j++)
		//{
		//	SparkSpeciesNames.push_back(_names[count++]);
		//	SparkSpeciesMassFractions.push_back(atof(_names[count++].c_str()));
		//}
	}
}
