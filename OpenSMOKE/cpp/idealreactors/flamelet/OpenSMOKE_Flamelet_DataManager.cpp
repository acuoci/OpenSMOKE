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
#include "idealreactors/flamelet/OpenSMOKE_Flamelet_DataManager.h"
#include "idealreactors/flamelet/OpenSMOKE_Flamelet.h"

void OpenSMOKE_Flamelet_DataManager::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flamelet_DataManager"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flamelet_DataManager::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flamelet_DataManager"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_Flamelet_DataManager::OpenSMOKE_Flamelet_DataManager()
{
	name_object			= "[Name not assigned]";
	iBackUp				= false;
	iGlobalKinetics		= false;
	i2E					= false;
	iDensityCorrection	= false;
	iEnthalpy			= true;
	iEnthalpyDefect		= false;

	iAssignedFlamePosition = false;
	iAssignedVideoSteps = false;
	iAssignedFileSteps = false;
	iAssignedBackupSteps = false;
	iAssignedRelativeTolerance = false;
	iAssignedAbsoluteTolerance = false;
	iAssignedEnvironmentTemperature = false;
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

	// Grid refinement
	nDiff = 4;
	nGrad = 2;	
	deltaDiff = 0.05;
	deltaGrad = 0.05;

	// Output
	nStepsVideo  = 10;		
	nStepsFile	 = 1000;
	nStepsBackUp = 2000;

	// Radiation
	iRadiation = 0;
	environmentTemperature = 298.15;	
	
	// Tolerances
	relTolerances = 1.e2;				
	absTolerances = 1.e-10;	

	// Enthalpy defect
	enthalpyDefect = 0.;
	minimumTemperature = Constants::T_Zero;
}

void OpenSMOKE_Flamelet_DataManager::SetName(const string name)
{
	name_object = name;
}

void OpenSMOKE_Flamelet_DataManager::Assign(OpenSMOKE_ReactingGas *_mix)
{
	mix = _mix;
}

void OpenSMOKE_Flamelet_DataManager::Assign(OpenSMOKE_Flamelet *_flamelet)
{
	flamelet = _flamelet;
}

void OpenSMOKE_Flamelet_DataManager::ReadFromFile(const string fileName)
{
	DefineFromFile(fileName);

	iteration				= 1;
	iterationVideoCounter	= nStepsVideo;
	iterationFileCounter	= nStepsFile;

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
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//									DICTIONARY													   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_Flamelet::OpenSMOKE_Dictionary_Flamelet()
{
    SetupBase();
	SetName("OpenSMOKE_Flamelet Dictionary");

	Add("#Fuel",          'C', 'S', "Fuel name");
	Add("#Oxidizer",      'C', 'S', "Oxidizer name");
	Add("#Inert",         'C', 'S', "Inert name");

	Add("#Pressure",				'C', 'M', "Pressure");

	Add("#FuelTemperature",			'C', 'M', "Fuel stream temperature");
	Add("#OxidizerTemperature",		'C', 'M', "Oxidizer stream temperature");

	Add("#FuelMassFractions",		'O', 'L', "Fuel stream mass fractions");
    Add("#FuelMoleFractions",		'O', 'L', "Fuel stream mole fractions");
	Add("#OxidizerMassFractions",   'O', 'L', "Oxidizer stream mass fractions");
    Add("#OxidizerMoleFractions",	'O', 'L', "Oxidizer stream mole fractions");

	Add("#Points",					'C', 'I', "Number of grid points");
	
	Add("#OutputSpecies",			'C', 'V', "Output species");

	Add("#FlameTemperature",		'C', 'M', "Flame temperature");
   	Add("#FlamePosition",			'O', 'D', "Flame position from fuel side");

	Add("#ScalarDissipationRates",	'C', 'V', "List of scalar dissipation rates");

	Add("#GridRefineGradient",      'O', 'I', "Maximum number of new grid points in grid refinement (gradient)");
	Add("#GridRefineCurvature",     'O', 'I', "Maximum number of new grid points in grid refinement (curvature)");
	Add("#GridRefineGradientStep",  'O', 'D', "Grid refinement control parameter (gradient)");
	Add("#GridRefineCurvatureStep", 'O', 'D', "Grid refinement control parameter (curvature)");
 
	Add("#RelativeTolerance",		'O', 'D', "Relative tolerance");
	Add("#AbsoluteTolerance",		'O', 'D', "Absolute tolerance");

	Add("#GasRadiation",			'O', 'N', "Radiative heat transfer from gas phase");
	Add("#SootRadiation",			'O', 'N', "Radiative heat transfer from gas phase and soot");
	Add("#EnvironmentTemperature",	'O', 'M', "Environment temperature");

	Add("#EnthalpyDefect",			'O', 'M', "Enthalpy defect");
	Add("#MinimumTemperature",		'O', 'M', "Minimum temperature");

	Add("#ReactionRates",			'O', 'V', "Reaction rates on file (list of reaction indices)");
	Add("#FormationRates",			'O', 'V', "Formation rates on file (list of species names)");
	Add("#ROPA",					'O', 'V', "Rate of Production Analysis (list of species names)");
	Add("#Sensitivity",				'O', 'V', "Sensitivity Analysis (list of species names)");

	Add("#DensityCorrection",		'O', 'N', "Density correction on scalar dissipation rate");

	Compulsory("#FuelMassFractions",     "#FuelMoleFractions");
	Compulsory("#OxidizerMassFractions", "#OxidizerMoleFractions");

    Conflict("#FuelMassFractions",		"#FuelMoleFractions");
    Conflict("#OxidizerMassFractions",	"#OxidizerMoleFractions");
	Conflict("#GasRadiation",			"#SootRadiation");

    Lock();
}

void OpenSMOKE_Flamelet_DataManager::DefineFromFile(const string inputFile)
{
    int     int_value;
	double  double_value;
    string  string_value;
	vector<string> string_vector;
	vector<double> double_vector;


	OpenSMOKE_Dictionary_Flamelet dictionary;
    dictionary.ParseFile(inputFile);

	if (dictionary.Return("#Fuel", string_value))
        AssignFuel(string_value);

	if (dictionary.Return("#Oxidizer", string_value))
        AssignOxidizer(string_value);

	if (dictionary.Return("#Inert", string_value))
        AssignInert(string_value);

	if (dictionary.Return("#Pressure", double_value, string_value))
        AssignPressure(string_value, double_value);

	if (dictionary.Return("#FuelTemperature", double_value, string_value))
        AssignFuelTemperature(string_value, double_value);

	if (dictionary.Return("#OxidizerTemperature", double_value, string_value))
        AssignOxidizerTemperature(string_value, double_value);

    if (dictionary.Return("#FuelMassFractions", double_vector, string_vector))
		AssignFuelMassFractions(string_vector, double_vector);
    if (dictionary.Return("#FuelMoleFractions", double_vector, string_vector))
		AssignFuelMoleFractions(string_vector, double_vector);

    if (dictionary.Return("#OxidizerMassFractions", double_vector, string_vector))
		AssignOxidizerMassFractions(string_vector, double_vector);
    if (dictionary.Return("#OxidizerMoleFractions", double_vector, string_vector))
		AssignOxidizerMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#Points", int_value))
        AssignGridPoints(int_value);
	
	if (dictionary.Return("#OutputSpecies", string_vector))
		AssignOutputSpecies(string_vector);

	if (dictionary.Return("#FlameTemperature", double_value, string_value))
        AssignFlameTemperature(string_value, double_value);

	if (dictionary.Return("#FlamePosition", double_value))
        SetFlamePosition(double_value);

	if (dictionary.Return("#ScalarDissipationRates", string_vector))
		AssignScalarDissipationRates(string_vector);

	
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

	if (dictionary.Return("#EnvironmentTemperature", double_value, string_value))
		SetEnvironmentTemperature(string_value, double_value);

	if (dictionary.Return("#EnthalpyDefect", double_value, string_value))
		SetEnthalpyDefect(string_value, double_value);

	if (dictionary.Return("#MinimumTemperature", double_value, string_value))
		SetMinimumTemperature(string_value, double_value);

	if (dictionary.Return("#DensityCorrection"))	
		SetDensityCorrection();

	if (dictionary.Return("#GasRadiation"))	
		SetGasRadiation();

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

	Lock();
}

void OpenSMOKE_Flamelet_DataManager::Lock()
{
}


void OpenSMOKE_Flamelet_DataManager::AssignFuel(const string string_value)
{
	nameFuel = string_value;
	jFUEL = mix->recognize_species(nameFuel);
}

void OpenSMOKE_Flamelet_DataManager::AssignOxidizer(const string string_value)
{
	nameOxidizer = string_value;
	jO2 = mix->recognize_species(nameOxidizer);
}

void OpenSMOKE_Flamelet_DataManager::AssignInert(const string string_value)
{
	nameInert = string_value;
	jINERT = mix->recognize_species(nameInert);
}

void OpenSMOKE_Flamelet_DataManager::AssignPressure(const string units, const double value)
{
    P_Pascal = OpenSMOKE_Conversions::conversion_pressure(value, units);
	P_atm	 = P_Pascal / OpenSMOKE_Conversions::Pa_from_atm;
	P_bar	 = P_Pascal / OpenSMOKE_Conversions::Pa_from_bar;
}

void OpenSMOKE_Flamelet_DataManager::AssignFlameTemperature(const string units, const double value)
{
    Tpeak = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Flamelet_DataManager::AssignFuelTemperature(const string units, const double value)
{
    TC = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Flamelet_DataManager::AssignOxidizerTemperature(const string units, const double value)
{
    TO = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Flamelet_DataManager::AssignFuelMassFractions(const vector<string> _names, const vector<double> _values)
{
	int i;
	double MWmix;
	BzzVector y(mix->NumberOfSpecies());
	ChangeDimensions(mix->NumberOfSpecies(), &XC);
	
	for(i=1;i<=int(_names.size());i++)
		y[mix->recognize_species(_names[i-1])] = _values[i-1];
	mix->GetMWAndMoleFractionsFromMassFractions(MWmix, XC, y);

	ChangeDimensions(0, &iXC);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XC[i]!=0.)	iXC.Append(i);
}

void OpenSMOKE_Flamelet_DataManager::AssignFuelMoleFractions(const vector<string> _names, const vector<double> _values)
{
	int i;
	ChangeDimensions(mix->NumberOfSpecies(), &XC);

	for(i=1;i<=int(_names.size());i++)
	{
		cout << _names[i-1] << "  " << _values[i-1] << endl;
		XC[mix->recognize_species(_names[i-1])] = _values[i-1];
	}

	ChangeDimensions(0, &iXC);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XC[i]!=0.)	iXC.Append(i); 
}
	
void OpenSMOKE_Flamelet_DataManager::AssignOxidizerMassFractions(const vector<string> _names, const vector<double> _values)
{
	int i;
	double MWmix;
	BzzVector y(mix->NumberOfSpecies());
	ChangeDimensions(mix->NumberOfSpecies(), &XO);
	
	for(i=1;i<=int(_names.size());i++)
		y[mix->recognize_species(_names[i-1])] = _values[i-1];
	mix->GetMWAndMoleFractionsFromMassFractions(MWmix, XO, y);

	ChangeDimensions(0, &iXO);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XO[i]!=0.)	iXO.Append(i); 	
}
    
void OpenSMOKE_Flamelet_DataManager::AssignOxidizerMoleFractions(const vector<string> _names, const vector<double> _values)
{
	int i;
	ChangeDimensions(mix->NumberOfSpecies(), &XO);

	for(i=1;i<=int(_names.size());i++)
		XO[mix->recognize_species(_names[i-1])] = _values[i-1];

	ChangeDimensions(0, &iXO);
	for(i=1;i<=mix->NumberOfSpecies();i++)
		if (XO[i]!=0.)	iXO.Append(i); 	
}

void OpenSMOKE_Flamelet_DataManager::AssignGridPoints(const int int_value)
{
	Np = int_value;
}

void OpenSMOKE_Flamelet_DataManager::AssignOutputSpecies(const vector<string> string_vector)
{
	if (string_vector.size()==1 && string_vector[0] == "ALL")
	{
		int iOutput = mix->NumberOfSpecies();
		ChangeDimensions(iOutput, &iOut);
		nameOutput.resize(iOutput+1); 
		for(int i=1;i<=iOutput;i++)
		{
			nameOutput[i]	= mix->names[i];
			iOut[i]			= i;
		}
	}
	else
	{
		int iOutput = string_vector.size();
		ChangeDimensions(iOutput, &iOut);
		nameOutput.resize(iOutput+1);
		for(int i=1;i<=iOutput;i++)
		{
			nameOutput[i] = string_vector[i-1];
			iOut[i] = mix->recognize_species(string_vector[i-1]);
		}
	}
}

void OpenSMOKE_Flamelet_DataManager::SetFlamePosition(const double value)
{
    xcen = value;
	iAssignedFlamePosition = true;
}

void OpenSMOKE_Flamelet_DataManager::SetVideoSteps(const int int_value)
{
	nStepsVideo = int_value;
	iAssignedVideoSteps = true;
}

void OpenSMOKE_Flamelet_DataManager::SetFileSteps(const int int_value)
{
	nStepsFile = int_value;
	iAssignedFileSteps = true;
}

void OpenSMOKE_Flamelet_DataManager::SetBackupSteps(const int int_value)
{
	nStepsBackUp = int_value;
	iAssignedBackupSteps = true;
}

void OpenSMOKE_Flamelet_DataManager::SetRelativeTolerance(const double double_value)
{
	relTolerances = double_value;
	iAssignedRelativeTolerance = true;
}

void OpenSMOKE_Flamelet_DataManager::SetAbsoluteTolerance(const double double_value)
{
	absTolerances = double_value;
	iAssignedAbsoluteTolerance = true;
}

void OpenSMOKE_Flamelet_DataManager::SetEnvironmentTemperature(const string units, const double value)
{
	environmentTemperature = OpenSMOKE_Conversions::conversion_temperature(value, units);
	iAssignedEnvironmentTemperature = true;
}

void OpenSMOKE_Flamelet_DataManager::SetEnthalpyDefect(const string units, const double value)
{
	enthalpyDefect = OpenSMOKE_Conversions::conversion_specificEnergy(value, units);
	iEnthalpyDefect = true;
}

void OpenSMOKE_Flamelet_DataManager::SetMinimumTemperature(const string units, const double value)
{
	minimumTemperature = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Flamelet_DataManager::SetDensityCorrection()
{
	iDensityCorrection = true;
}

void OpenSMOKE_Flamelet_DataManager::SetGasRadiation()
{
	iRadiation = 1;
	iAssignedGasRadiation = true;
}

void OpenSMOKE_Flamelet_DataManager::SetSootRadiation()
{
	iRadiation = 2;
	iAssignedSootRadiation = true;
}

void OpenSMOKE_Flamelet_DataManager::SetGridRefineGradient(const int int_value)
{
	nDiff = int_value;
	iAssignedGridRefineGradient = true;
}

void OpenSMOKE_Flamelet_DataManager::SetGridRefineCurvature(const int int_value)
{
	nGrad = int_value;
	iAssignedGridRefineCurvature = true;
}

void OpenSMOKE_Flamelet_DataManager::SetGridRefineGradientStep(const double double_value)
{
	deltaDiff = double_value;
	iAssignedGridRefineGradientStep = true;
}

void OpenSMOKE_Flamelet_DataManager::SetGridRefineCurvatureStep(const double double_value)
{
	deltaGrad = double_value;
	iAssignedGridRefineCurvatureStep = true;
}

void OpenSMOKE_Flamelet_DataManager::SetSensitivityOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_sensitivity, names_sensitivity);
	iAssignedSensitivity = true;
}

void OpenSMOKE_Flamelet_DataManager::SetFormationRatesOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_formation_rates, names_formation_rates);
	iAssignedFormationRates = true;
}

void OpenSMOKE_Flamelet_DataManager::SetROPAOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_ROPA, names_ROPA);
	iAssignedROPA = true;
}

void OpenSMOKE_Flamelet_DataManager::SetReactionRatesOnFile(const vector<string> _names)
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

void OpenSMOKE_Flamelet_DataManager::AssignScalarDissipationRates(const vector<string> string_vector)
{
	int nScalarDissipationRates = string_vector.size()-1;
	string units = string_vector[nScalarDissipationRates];
	ChangeDimensions(nScalarDissipationRates, &listScalarDissipationRates);

	for(int i=1;i<=nScalarDissipationRates;i++)
		listScalarDissipationRates[i] = OpenSMOKE_Conversions::conversion_frequency(atof(string_vector[i-1].c_str()), units);
	chiSt = listScalarDissipationRates[1];	
}
	
void OpenSMOKE_Flamelet_DataManager::StoichiometricMixtureFraction()
{
//	int j;

/*
	double nC = 0.;
	double nH = 0.;

	int jC = mix->recognize_element_without_error("c");
	int jH = mix->recognize_element_without_error("h");
	
	if (jC>0)
		for(j=1;j<=mix->NumberOfSpecies();j++)
			nC += mix->elements[jC][j]*XC[j];

	if (jH>0)
		for(j=1;j<=mix->NumberOfSpecies();j++)
			nH += mix->elements[jH][j]*XC[j];

	double MWO = mix->M[jO2];
	double MWF = mix->GetMWFromMoleFractions(XC);
	double nuF = 1.;
	double nuO = nC*1. + nH*0.25;
	double s = nuO*MWO/nuF/MWF;
	double omegaFuel	 = 1.;
	double omegaOxidizer = XO[jO2]*MWO/mix->GetMWFromMoleFractions(XO);
	double phi = s*omegaFuel/omegaOxidizer;
	zStoichiometric = 1./(1.+phi);

	cout << "Stoichiometric Ratio:            " << s << endl;
	cout << "Equivalence Ratio:               " << phi << endl;
	cout << "Stoichiometric Mixture Fraction: " << zStoichiometric << endl;
*/

	// Correction coefficient evaluation
	OpenSMOKE_GasStream stream;
	vector<string> fuel_names;
	vector<string> oxidizer_names;
	vector<double> fuel_values;
	vector<double> oxidizer_values;
	BzzVector omega_elemental_vector(mix->NumberOfElements());
	BzzVector omega_elemental_fuel_vector(mix->NumberOfElements());
	BzzVector omega_elemental_oxidizer_vector(mix->NumberOfElements());
	BzzVector fuel_mole_fractions(mix->NumberOfSpecies());
	BzzVector oxidizer_mole_fractions(mix->NumberOfSpecies());

	for(int j=1;j<=mix->NumberOfSpecies();j++)
	{
		if (XC[j] > 0.)	{ fuel_names.push_back(mix->names[j]); fuel_values.push_back(XC[j]); fuel_mole_fractions[j]=XC[j];}
		if (XO[j] > 0.)	{ oxidizer_names.push_back(mix->names[j]); oxidizer_values.push_back(XO[j]); oxidizer_mole_fractions[j]=XO[j];}
	}

	stream.AssignKineticScheme(*mix);
	stream.AssignTemperature(Constants::T_Reference, "K");
	stream.AssignPressure(Constants::P_Reference, "Pa");
	stream.ChangeMassFlowRate(1.0, "kg/s");
	stream.AssignFuelMoleFractions(fuel_names, fuel_values);
	stream.AssignOxidizerMoleFractions(oxidizer_names, oxidizer_values);
	stream.AssignEquivalenceRatio(1.0);
	stream.lock();
	mix->GetElementalMassFractionsFromSpeciesMassFractions(omega_elemental_vector, stream.omega);
	mix->GetElementalMassFractionsFromSpeciesMoleFractions(omega_elemental_fuel_vector, fuel_mole_fractions);
	mix->GetElementalMassFractionsFromSpeciesMoleFractions(omega_elemental_oxidizer_vector, oxidizer_mole_fractions);
	zStoichiometric = mix->GetMixtureFraction(omega_elemental_vector, omega_elemental_fuel_vector, omega_elemental_oxidizer_vector);
	double Zst = 1.-zStoichiometric;

	BzzInverseErrorFunction erf;
	correctionChi = exp(-2.*BzzPow2(erf.at(1.-2.*zStoichiometric)));

	cout << "Stoichiometric mixture fraction (fluent): " << zStoichiometric << endl;
	cout << "Stoichiometric mixture fraction (os):     " << Zst << endl;
	cout << "Reduction coefficient:                    " << correctionChi << endl;
	cout << "Scalar dissipation rate (max.):           " << chiSt/correctionChi << " 1/s" << endl;
	cout << "Scalar dissipation rate (st.):            " << chiSt				<< " 1/s" << endl;
	cout << "Strain rate:                              " << chiSt*Constants::pi << " 1/s" << endl;
}	