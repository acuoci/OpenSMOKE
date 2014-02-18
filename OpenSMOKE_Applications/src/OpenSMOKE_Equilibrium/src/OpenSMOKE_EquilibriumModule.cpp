/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
 *   alberto.cuoci@polimi.it   						                       *
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

#include <iomanip>
#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_IdealGas.h"
#include "engine/OpenSMOKE_EquilibriumStanjan.h"
#include "OpenSMOKE_EquilibriumModule.h"

void OpenSMOKE_EquilibriumModule::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_EquilibriumModule"	<< endl;
    cout << "Object: " << name_object				<< endl;
    cout << "Error:  " << message					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_EquilibriumModule::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_EquilibriumModule"	<< endl;
    cout << "Object:  " << name_object				<< endl;
    cout << "Warning: " << message					<< endl;
    cout << "Press enter to continue..."			<< endl;
	getchar();
	cout << endl;
}

OpenSMOKE_EquilibriumModule::OpenSMOKE_EquilibriumModule()
{
	name_object	= "[Name not assigned]";
	Reset();
}

void OpenSMOKE_EquilibriumModule::Reset()
{
	iVerbose				= false;
	iSetTemperature			= false;
	iSetPressure			= false;
	iSetEnthalpy			= false;
	iSetEntropy				= false;
	iSetInternalEnergy		= false;
	iSetDensity				= false;
	iSetVolume				= false;
	iSetSpeciesComposition	= false;
	iSetEquivalenceRatio	= false;
	iSetFuelComposition		= false;
	iSetOxidizerComposition	= false;


	T_0			= 0.;	// [K]
	P_0_Pascal	= 0.;	// [Pa]
	MWmix_0		= 0.;	// [kg/kmol]
	rho_0		= 0.;	// [kg/m3]
	V_0			= 0.;	// [m3/kg]
	N_0			= 1.;	// [kmol]
	Cp_0		= 0.;	// [J/kg/K]
	H_0			= 0.;	// [J/kg]
	U_0			= 0.;	// [J/kg]
	S_0			= 0.;	// [J/kg/K]
	G_0			= 0.;	// [J/kg]
	A_0			= 0.;	// [J/kg]

	ChangeDimensions(0,	&x_0);
	ChangeDimensions(0,	&omega_0);
}

void OpenSMOKE_EquilibriumModule::SetName(const string name)
{
	name_object = name;
}

void OpenSMOKE_EquilibriumModule::SetInputFileName(const string name)
{
	name_input_file = name;
}

void OpenSMOKE_EquilibriumModule::SetVerbose()
{
	iVerbose = true;
}

void OpenSMOKE_EquilibriumModule::UnsetVerbose()
{
	iVerbose = false;
}

void OpenSMOKE_EquilibriumModule::Assign(OpenSMOKE_IdealGas *ideal_gas)
{
	gas = ideal_gas;
}

OpenSMOKE_Dictionary_EquilibriumModule::OpenSMOKE_Dictionary_EquilibriumModule()
{
    SetupBase();
	SetName("OpenSMOKE_EquilibriumModule Dictionary");

	Add("#Constant",				'C', 'V', "Assigned thermodynamic state functions");

	Add("#ElementalMassFractions",	'O', 'L', "Elemental mass fractions");
	Add("#ElementalMoleFractions",	'O', 'L', "Elemental mole fractions");
	
	Add("#MassFractions",			'O', 'L', "Species mass fractions");
	Add("#MoleFractions",			'O', 'L', "Species mole fractions");

	Add("#EquivalenceRatio",		'O', 'D', "Equivalence Ratio");
	Add("#FuelMoleFractions",		'O', 'L', "Fuel mole fractions");
	Add("#OxidizerMoleFractions",	'O', 'L', "Oxidizer mole fractions");
	Add("#FuelMassFractions",		'O', 'L', "Fuel mass fractions");
	Add("#OxidizerMassFractions",	'O', 'L', "Oxidizer mass fractions");

	Add("#Temperature",      'O', 'M', "Temperature");
	Add("#Pressure",         'O', 'M', "Pressure");
	Add("#Enthalpy",         'O', 'M', "Specific enthalpy");
	Add("#Entropy",          'O', 'M', "Specific entropy");
	Add("#InternalEnergy",   'O', 'M', "Internal energy");
	Add("#Volume",			 'O', 'M', "Specific volume");
	Add("#Density",			 'O', 'M', "Density");

	Add("#Verbose",			 'O', 'N', "Verbose");

	Conflict("#MassFractions",			"#MoleFractions");
	Conflict("#ElementalMassFractions",	"#ElementalMoleFractions");
	Conflict("#MassFractions",			"#ElementalMoleFractions");
	Conflict("#MoleFractions",			"#ElementalMassFractions");
	Conflict("#MassFractions",			"#ElementalMassFractions");
	Conflict("#MoleFractions",			"#ElementalMoleFractions");
	Conflict("#FuelMoleFractions",		"#FuelMassFractions");
	Conflict("#OxidizerMoleFractions",	"#OxidizerMassFractions");
	Conflict("#FuelMoleFractions",		"#MassFractions");
	Conflict("#FuelMoleFractions",		"#MoleFractions");
	Conflict("#FuelMoleFractions",		"#ElementalMoleFractions");
	Conflict("#FuelMoleFractions",		"#ElementalMassFractions");
	Conflict("#EquivalenceRatio",		"#MassFractions");
	Conflict("#EquivalenceRatio",		"#MoleFractions");
	Conflict("#EquivalenceRatio",		"#ElementalMoleFractions");
	Conflict("#EquivalenceRatio",		"#ElementalMassFractions");

    Lock();
}

void OpenSMOKE_EquilibriumModule::Lock()
{
	int counter = 0;
	if (iSetEnthalpy		== true)	counter++;
	if (iSetTemperature		== true)	counter++;
	if (iSetPressure		== true)	counter++;
	if (iSetDensity			== true)	counter++;
	if (iSetEntropy			== true)	counter++;
	if (iSetInternalEnergy	== true)	counter++;
	if (iSetVolume			== true)	counter++;

	if (counter != 2)
		ErrorMessage("Only 2 thermodynamic properties between (T, P, RHO,V, H, U, S) must be assigned...");

	
	if (iSetSpeciesComposition	== true)	
	{
		GetInitialConditions();
		AllProperties_0();
	}
	
	else
	{
		if(iEquilibriumProblem == EQUILIBRIUM_TP)
		{
			if (iSetPressure	== false)	ErrorMessage("Pressure is required for a (T,P) equilibrium problem if composition is not assigned...");
			if (iSetTemperature == false)	ErrorMessage("Temperature is required for a (T,P) equilibrium problem if composition is not assigned...");
		}

		else if(iEquilibriumProblem == EQUILIBRIUM_HP)
		{
			if (iSetPressure == false)	ErrorMessage("Pressure is required for a (H,P) equilibrium problem if composition is not assigned...");
			if (iSetEnthalpy == false)	ErrorMessage("Enthalpy is required for a (H,P) equilibrium problem if composition is not assigned...");
		}

		else if(iEquilibriumProblem == EQUILIBRIUM_SP)
		{
			if (iSetPressure == false)	ErrorMessage("Pressure is required for a (S,P) equilibrium problem if composition is not assigned...");
			if (iSetEntropy== false)	ErrorMessage("Entropy is required for a (S,P) equilibrium problem if composition is not assigned...");
		}

		else if(iEquilibriumProblem == EQUILIBRIUM_TV || iEquilibriumProblem == EQUILIBRIUM_TRHO)
		{
			if (iSetTemperature == false)						ErrorMessage("Temperature is required for a (T,V) or (T,RHO) equilibrium problem if composition is not assigned...");
			if (iSetVolume  == false &&	iSetDensity == false)	ErrorMessage("Density (or Specific Volume) is required for a (T,V) or (T,RHO) equilibrium problem if composition is not assigned...");
		}

		else if(iEquilibriumProblem == EQUILIBRIUM_UV || iEquilibriumProblem == EQUILIBRIUM_URHO)
		{
			if (iSetInternalEnergy == false)					ErrorMessage("Internal Energy is required for a (U,V) or (U,RHO) equilibrium problem if composition is not assigned...");
			if (iSetVolume  == false && iSetDensity == false)	ErrorMessage("Density (or Specific Volume) is required for a (U,V) or (U,RHO) equilibrium problem if composition is not assigned...");
		}

		else if(iEquilibriumProblem == EQUILIBRIUM_SV || iEquilibriumProblem == EQUILIBRIUM_SRHO)
		{
			if (iSetEntropy == false)							ErrorMessage("Entropy is required for a (S,V) or (S,RHO) equilibrium problem if composition is not assigned...");
			if (iSetVolume  == false &&	iSetDensity == false)	ErrorMessage("Density (or Specific volume) is required for a (S,V) or (S,RHO) equilibrium problem if composition is not assigned...");
		}
	}
}

void OpenSMOKE_EquilibriumModule::AllProperties_0()
{
	rho_0	 = P_0_Pascal*MWmix_0/Constants::R_J_kmol/T_0;					// [kg/m3]
	V_0		 = 1./rho_0;													// [m3/kg]
	N_0		 = 1.;															// [kmol]

	gas->SpeciesCp(T_0);
	Cp_0  = gas->MixCp_FromMoleFractions(x_0);									// [J/kg/K]
	H_0   = gas->GetMixEnthalpy_Mole(T_0, x_0) / MWmix_0;						// [J/kg]
	U_0   = gas->GetMixInternalEnergy_Mole(T_0, x_0) / MWmix_0;					// [J/kg]
	S_0   = gas->GetMixEntropy_Mole(P_0_Pascal, T_0, x_0) / MWmix_0;			// [J/kg/K]
	G_0   = gas->GetMixGibbsFreeEnergy_Mole(P_0_Pascal, T_0, x_0) / MWmix_0;	// [J/kg]
	A_0   = gas->GetMixHelmotzFreeEnergy_Mole(P_0_Pascal, T_0, x_0) / MWmix_0;	// [J/kg]
}

void OpenSMOKE_EquilibriumModule::AllProperties_E()
{
	gas->GetMWAndMassFractionsFromMoleFractions(MWmix_E, omega_E, x_E);

	rho_E = P_E_Pascal*MWmix_E/Constants::R_J_kmol/T_E;					// [kg/m3]
	V_E	  = 1./rho_E;													// [m3/kg]

	gas->SpeciesCp(T_E);
	Cp_E  = gas->MixCp_FromMoleFractions(x_E);									// [J/kg/K]

	H_E   = gas->GetMixEnthalpy_Mole(T_E, x_E) / MWmix_E;						// [J/kg]
	U_E   = gas->GetMixInternalEnergy_Mole(T_E, x_E) / MWmix_E;					// [J/kg]
	S_E   = gas->GetMixEntropy_Mole(P_E_Pascal, T_E, x_E) / MWmix_E;			// [J/kg/K]
	G_E   = gas->GetMixGibbsFreeEnergy_Mole(P_E_Pascal, T_E, x_0) / MWmix_E;	// [J/kg]
	A_E   = gas->GetMixHelmotzFreeEnergy_Mole(P_E_Pascal, T_E, x_0) / MWmix_E;	// [J/kg]
}

void OpenSMOKE_EquilibriumModule::DefineFromFile()
{
    OpenSMOKE_Dictionary_EquilibriumModule dictionary;

    dictionary.ParseFile(name_input_file);
	CheckDictionary(dictionary);
	Lock();
}

void OpenSMOKE_EquilibriumModule::CheckDictionary(OpenSMOKE_Dictionary_EquilibriumModule &dictionary)
{
	double  double_value;
    string  string_value;
	vector<string> string_vector;
	vector<double> double_vector;

	if (dictionary.Return("#Constant", string_vector))
        AssignConstant(string_vector);

    if (dictionary.Return("#ElementalMassFractions", double_vector, string_vector))
		AssignElementalMassFractions(string_vector, double_vector);
    if (dictionary.Return("#ElementalMoleFractions", double_vector, string_vector))
		AssignElementalMoleFractions(string_vector, double_vector);

    if (dictionary.Return("#MassFractions", double_vector, string_vector))
		AssignMassFractions(string_vector, double_vector);
    if (dictionary.Return("#MoleFractions", double_vector, string_vector))
		AssignMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#FuelMoleFractions", double_vector, string_vector))
		AssignFuelMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#FuelMassFractions", double_vector, string_vector))
		AssignFuelMassFractions(string_vector, double_vector);

	if (dictionary.Return("#OxidizerMoleFractions", double_vector, string_vector))
		AssignOxidizerMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#OxidizerMassFractions", double_vector, string_vector))
		AssignOxidizerMassFractions(string_vector, double_vector);

	if (dictionary.Return("#EquivalenceRatio", double_value))
		AssignEquivalenceRatio(double_value);

	if (dictionary.Return("#Temperature", double_value, string_value))
        SetTemperature(string_value, double_value);

	if (dictionary.Return("#Pressure", double_value, string_value))
        SetPressure(string_value, double_value);

	if (dictionary.Return("#Enthalpy", double_value, string_value))
        SetEnthalpy(string_value, double_value);

	if (dictionary.Return("#Entropy", double_value, string_value))
        SetEntropy(string_value, double_value);

	if (dictionary.Return("#InternalEnergy", double_value, string_value))
        SetInternalEnergy(string_value, double_value);

	if (dictionary.Return("#Volume", double_value, string_value))
        SetVolume(string_value, double_value);

	if (dictionary.Return("#Density", double_value, string_value))
        SetDensity(string_value, double_value);

	if (dictionary.Return("#Verbose"))
        SetVerbose();
}

void OpenSMOKE_EquilibriumModule::AssignConstant(const vector<string> string_vector)
{
	if (string_vector.size() > 2)
		ErrorMessage("Only 2 thermodynamic state functions must be assigned: (T,P) || (H,P) || (S,P) || (T,V) || (T,RHO) || (U,V) || (U,RHO) || (S,V) || (S,RHO)");
	if (string_vector[0] == string_vector[1])
		ErrorMessage("The assigned thermodynamic state functions must be different: (T,P) || (H,P) || (S,P) || (T,V) || (T,RHO) || (U,V) || (U,RHO) || (S,V) || (S,RHO)");

		 if ( (string_vector[0] == "T" && string_vector[1] == "P") || (string_vector[0] == "P" && string_vector[1] == "T"))
		iEquilibriumProblem = EQUILIBRIUM_TP;
	else if ( (string_vector[0] == "H" && string_vector[1] == "P") || (string_vector[0] == "P" && string_vector[1] == "H"))
		iEquilibriumProblem = EQUILIBRIUM_HP;
	else if ( (string_vector[0] == "S" && string_vector[1] == "P") || (string_vector[0] == "P" && string_vector[1] == "S"))
		iEquilibriumProblem = EQUILIBRIUM_SP;
	else if ( (string_vector[0] == "T" && string_vector[1] == "V") || (string_vector[0] == "V" && string_vector[1] == "T"))
		iEquilibriumProblem = EQUILIBRIUM_TV;
	else if ( (string_vector[0] == "T" && string_vector[1] == "RHO") || (string_vector[0] == "RHO" && string_vector[1] == "T"))
		iEquilibriumProblem = EQUILIBRIUM_TRHO;
	else if ( (string_vector[0] == "U" && string_vector[1] == "V") || (string_vector[0] == "V" && string_vector[1] == "U"))
		iEquilibriumProblem = EQUILIBRIUM_UV;
	else if ( (string_vector[0] == "U" && string_vector[1] == "RHO") || (string_vector[0] == "RHO" && string_vector[1] == "U"))
		iEquilibriumProblem = EQUILIBRIUM_URHO;
	else if ( (string_vector[0] == "S" && string_vector[1] == "V") || (string_vector[0] == "V" && string_vector[1] == "S"))
		iEquilibriumProblem = EQUILIBRIUM_SV;
	else if ( (string_vector[0] == "S" && string_vector[1] == "RHO") || (string_vector[0] == "RHO" && string_vector[1] == "S"))
		iEquilibriumProblem = EQUILIBRIUM_SRHO;
	else
		ErrorMessage("The assigned thermodynamic state functions must be chosen between the following couples: (T,P) || (H,P) || (S,P) || (T,V) || (T,RHO) || (U,V) || (U,RHO) || (S,V) || (S,RHO)");
}

void OpenSMOKE_EquilibriumModule::AssignElementalMassFractions(const vector<string> _names, const vector<double> _values)
{
	for(int i=1;i<=int(_names.size());i++)
		omega_elements[gas->recognize_element(_names[i-1])] = _values[i-1];
	gas->GetMWAndElementalMoleFractionsFromElementalMassFractions(MWmix_0, x_elements, omega_elements);
}

void OpenSMOKE_EquilibriumModule::AssignElementalMoleFractions(const vector<string> _names, const vector<double> _values)
{
	for(int i=1;i<=int(_names.size());i++)
		x_elements[gas->recognize_element(_names[i-1])] = _values[i-1];
	gas->GetMWAndElementalMassFractionsFromElementalMoleFractions(MWmix_0, omega_elements, x_elements);
}

void OpenSMOKE_EquilibriumModule::AssignMassFractions(const vector<string> _names, const vector<double> _values)
{
	for(int i=1;i<=int(_names.size());i++)
		omega_0[gas->recognize_species(_names[i-1])] = _values[i-1];
	gas->GetMWAndMoleFractionsFromMassFractions(MWmix_0, x_0, omega_0);

	gas->GetElementalMoleFractionsFromSpeciesMoleFractions(x_elements, x_0);
	gas->GetElementalMassFractionsFromSpeciesMassFractions(omega_elements, omega_0);

	iSetSpeciesComposition = true;
}

void OpenSMOKE_EquilibriumModule::AssignMoleFractions(const vector<string> _names, const vector<double> _values)
{
	for(int i=1;i<=int(_names.size());i++)
		x_0[gas->recognize_species(_names[i-1])] = _values[i-1];
	gas->GetMWAndMassFractionsFromMoleFractions(MWmix_0, omega_0, x_0);

	gas->GetElementalMoleFractionsFromSpeciesMoleFractions(x_elements, x_0);
	gas->GetElementalMassFractionsFromSpeciesMassFractions(omega_elements, omega_0);

	iSetSpeciesComposition = true;
}

void OpenSMOKE_EquilibriumModule::AssignEquivalenceRatio(const double equivalence_ratio)
{
	vector<string>	fuel_names;
	BzzVector		fuel_moles;
	vector<string>  oxidizer_names;
	BzzVector		oxidizer_moles;

	if (iSetFuelComposition == false)
		ErrorMessage("EquivalenceRatio option requires fuel composition...");
	else
	{
		fuel_names.push_back("fuel names");
		for(int j=1;j<=gas->NumberOfSpecies();j++)
			if (x_Fuel[j] > 0.)
			{
				fuel_names.push_back(gas->names[j]);
				fuel_moles.Append(x_Fuel[j]);
			}
	}

	if (iSetOxidizerComposition == false)
		ErrorMessage("If #EquivalenceRatio option is used, oxidizer compisition must be specified!");
	else
	{
		oxidizer_names.push_back("oxidizer names");
		for(int j=1;j<=gas->NumberOfSpecies();j++)
			if (x_Oxidizer[j] > 0.)
			{
				oxidizer_names.push_back(gas->names[j]);
				oxidizer_moles.Append(x_Oxidizer[j]);
			}
	}

	BzzVector moles = gas->GetMoleFractionsFromEquivalenceRatio(equivalence_ratio, fuel_names, fuel_moles, oxidizer_names, oxidizer_moles);
	
	vector<string> _names;
	vector<double> _moles;
	for(int j=1;j<=gas->NumberOfSpecies();j++)
		if (moles[j] > 0.)
		{
			_names.push_back(gas->names[j]);
			_moles.push_back(moles[j]);
		}

	AssignMoleFractions(_names, _moles);

	iSetEquivalenceRatio = true;
}

void OpenSMOKE_EquilibriumModule::AssignFuelMoleFractions(const vector<string> _names, const vector<double> _values)
{
	ChangeDimensions(gas->NumberOfSpecies(), &x_Fuel);
	for(int i=1;i<=int(_names.size());i++)
		x_Fuel[gas->recognize_species(_names[i-1])] = _values[i-1];
	iSetFuelComposition = true;
}

void OpenSMOKE_EquilibriumModule::AssignOxidizerMoleFractions(const vector<string> _names, const vector<double> _values)
{
	ChangeDimensions(gas->NumberOfSpecies(), &x_Oxidizer);
	for(int i=1;i<=int(_names.size());i++)
		x_Oxidizer[gas->recognize_species(_names[i-1])] = _values[i-1];
	iSetOxidizerComposition = true;
}

void OpenSMOKE_EquilibriumModule::AssignFuelMassFractions(const vector<string> _names, const vector<double> _values)
{
	double MW_Fuel;
	BzzVector omega_Fuel(gas->NumberOfSpecies());
	ChangeDimensions(gas->NumberOfSpecies(), &x_Fuel);
	for(int i=1;i<=int(_names.size());i++)
		omega_Fuel[gas->recognize_species(_names[i-1])] = _values[i-1];
	gas->GetMWAndMoleFractionsFromMassFractions(MW_Fuel, x_Fuel, omega_Fuel);
	iSetFuelComposition = true;
}

void OpenSMOKE_EquilibriumModule::AssignOxidizerMassFractions(const vector<string> _names, const vector<double> _values)
{
	double MW_Oxidizer;
	BzzVector omega_Oxidizer(gas->NumberOfSpecies());
	ChangeDimensions(gas->NumberOfSpecies(), &x_Oxidizer);
	for(int i=1;i<=int(_names.size());i++)
		omega_Oxidizer[gas->recognize_species(_names[i-1])] = _values[i-1];
	gas->GetMWAndMoleFractionsFromMassFractions(MW_Oxidizer, x_Oxidizer, omega_Oxidizer);
	iSetOxidizerComposition = true;
}

void OpenSMOKE_EquilibriumModule::SetTemperature(const string units, const double value)
{
	T_0 = OpenSMOKE_Conversions::conversion_temperature(value, units);
	iSetTemperature = true;
}

void OpenSMOKE_EquilibriumModule::SetPressure(const string units, const double value)
{
	P_0_Pascal = OpenSMOKE_Conversions::conversion_pressure(value, units);
	iSetPressure = true;
}

void OpenSMOKE_EquilibriumModule::SetEnthalpy(const string units, const double value)
{
	H_0 = OpenSMOKE_Conversions::conversion_specificEnergy(value, units);
	iSetEnthalpy = true;
}

void OpenSMOKE_EquilibriumModule::SetEntropy(const string units, const double value)
{
	S_0 = OpenSMOKE_Conversions::conversion_specificEntropy(value, units);
	iSetEntropy = true;
}

void OpenSMOKE_EquilibriumModule::SetInternalEnergy(const string units, const double value)
{
	U_0 = OpenSMOKE_Conversions::conversion_specificEnergy(value, units);
	iSetInternalEnergy = true;
}

void OpenSMOKE_EquilibriumModule::SetDensity(const string units, const double value)
{
	rho_0	= OpenSMOKE_Conversions::conversion_density(value, units);
	V_0		= 1./rho_0;
	iSetDensity = true;
}

void OpenSMOKE_EquilibriumModule::SetVolume(const string units, const double value)
{
	V_0		= OpenSMOKE_Conversions::conversion_specificVolume(value, units);
	rho_0	= 1./V_0;
	iSetVolume = true;
}

void OpenSMOKE_EquilibriumModule::Run()
{
	ChangeDimensions(gas->NumberOfElements(),	&x_elements);
	ChangeDimensions(gas->NumberOfElements(),	&omega_elements);
	ChangeDimensions(gas->NumberOfSpecies(),	&x_0);
	ChangeDimensions(gas->NumberOfSpecies(),	&omega_0);
	ChangeDimensions(gas->NumberOfSpecies(),	&x_E);
	ChangeDimensions(gas->NumberOfSpecies(),	&omega_E);

	DefineFromFile();
	Solve();
	AllProperties_E();
	PrintSummary();
}

void OpenSMOKE_EquilibriumModule::Solve()
{
		 if (iEquilibriumProblem == EQUILIBRIUM_TP)		Solve_TP();
	else if (iEquilibriumProblem == EQUILIBRIUM_HP)		Solve_HP();
	else if (iEquilibriumProblem == EQUILIBRIUM_SP)		Solve_SP();
	else if (iEquilibriumProblem == EQUILIBRIUM_TV)		Solve_TRHO();
	else if (iEquilibriumProblem == EQUILIBRIUM_TRHO)	Solve_TRHO();
	else if (iEquilibriumProblem == EQUILIBRIUM_UV)		Solve_URHO();
	else if (iEquilibriumProblem == EQUILIBRIUM_URHO)	Solve_URHO();
	else if (iEquilibriumProblem == EQUILIBRIUM_SV)		Solve_SRHO();
	else if (iEquilibriumProblem == EQUILIBRIUM_SRHO)	Solve_SRHO();
}

void OpenSMOKE_EquilibriumModule::Solve_TP()
{
	double tStart = BzzGetCpuTime();
	gas->Equilibrium_TP(x_E, N_E, T_0, P_0_Pascal, x_elements, iVerbose);
	double tEnd = BzzGetCpuTime();
	cout << "Equilibrium CPU time: " <<  tEnd-tStart << " s" << endl;

	T_E = T_0;
	P_E_Pascal = P_0_Pascal;
}

void OpenSMOKE_EquilibriumModule::Solve_HP()
{
	// First guess
	T_E = T_0;
	N_E = 1.;

	int flag;
	double tStart = BzzGetCpuTime();
	gas->Equilibrium_HP(T_E, x_E, N_E, H_0, P_0_Pascal, x_elements, iVerbose, flag);
	double tEnd = BzzGetCpuTime();
	cout << "Equilibrium CPU time: " <<  tEnd-tStart << " s" << endl;

	P_E_Pascal = P_0_Pascal;
}

void OpenSMOKE_EquilibriumModule::Solve_SP()
{
	// First guess
	T_E = T_0;
	N_E = 1.;

	double tStart = BzzGetCpuTime();
	gas->Equilibrium_SP(T_E, x_E, N_E, S_0, P_0_Pascal, x_elements, iVerbose);
	double tEnd = BzzGetCpuTime();
	cout << "Equilibrium CPU time: " <<  tEnd-tStart << " s" << endl;

	P_E_Pascal = P_0_Pascal;
}

void OpenSMOKE_EquilibriumModule::Solve_TRHO()
{
	// First guess
	N_E = 1.;

	double tStart = BzzGetCpuTime();
	gas->Equilibrium_TRHO(P_E_Pascal, x_E, N_E, T_0, rho_0, x_elements, iVerbose);
	double tEnd = BzzGetCpuTime();
	cout << "Equilibrium CPU time: " <<  tEnd-tStart << " s" << endl;

	T_E = T_0;
}

void OpenSMOKE_EquilibriumModule::Solve_URHO()
{
	double tStart = BzzGetCpuTime();
	gas->Equilibrium_URHO(T_E, P_E_Pascal, x_E, N_E, U_0, rho_0, x_elements, iVerbose);
	double tEnd = BzzGetCpuTime();
	cout << "Equilibrium CPU time: " <<  tEnd-tStart << " s" << endl;
}

void OpenSMOKE_EquilibriumModule::Solve_SRHO()
{
	double tStart = BzzGetCpuTime();
	gas->Equilibrium_SRHO(T_E, P_E_Pascal, x_E, N_E, S_0, rho_0, x_elements, iVerbose);
	double tEnd = BzzGetCpuTime();
	cout << "Equilibrium CPU time: " <<  tEnd-tStart << " s" << endl;
}

void OpenSMOKE_EquilibriumModule::PrintSummary()
{
	int j;
	ofstream fSummary;
	openOutputFileAndControl(fSummary, "Summary.out");

	OpenSMOKE_logo(fSummary, "OpenSMOKE_Equilibrium");


	fSummary << setprecision(9);
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << "                         Equilibrium         Initial                 " << endl;
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << " Temperature        " << setw(16) << right << T_E			<< setw(16) << right << T_0			<< "   " << left << "[K]" << endl;
	fSummary << " Pressure           " << setw(16) << right << P_E_Pascal	<< setw(16) << right << P_0_Pascal	<< "   " << left << "[Pa]" << endl;
	fSummary << " Molecular Weight   " << setw(16) << right << MWmix_E		<< setw(16) << right << MWmix_0		<< "   " << left << "[kg/kmol]" << endl;
	fSummary << " Moles              " << setw(16) << right << N_E			<< setw(16) << right << N_0			<< "   " << left << "[kmol]" << endl;
	fSummary << " Density            " << setw(16) << right << rho_E		<< setw(16) << right << rho_0		<< "   " << left << "[kg/m3]" << endl;
	fSummary << " Volume             " << setw(16) << right << V_E			<< setw(16) << right << V_0			<< "   " << left << "[m3/kg]" << endl;
	fSummary << " Specific Heat      " << setw(16) << right << Cp_E/1e3		<< setw(16) << right << Cp_0/1e3	<< "   " << left << "[kJ/kg/K]" << endl;
	fSummary << " Enthalpy           " << setw(16) << right << H_E/1e3		<< setw(16) << right << H_0/1e3		<< "   " << left << "[kJ/kg]" << endl;
	fSummary << " Entropy            " << setw(16) << right << S_E/1e3		<< setw(16) << right << S_0/1e3		<< "   " << left << "[kJ/kg/K]" << endl;
	fSummary << " Internal Energy    " << setw(16) << right << U_E/1e3		<< setw(16) << right << U_0/1e3		<< "   " << left << "[kJ/kg]" << endl;
	fSummary << " Helmotz Energy     " << setw(16) << right << A_E/1e3		<< setw(16) << right << A_0/1e3		<< "   " << left << "[kJ/kg]" << endl;
	fSummary << " Gibbs Energy       " << setw(16) << right << G_E/1e3		<< setw(16) << right << G_0/1e3		<< "   " << left << "[kJ/kg]" << endl;
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << endl;

	fSummary << setprecision(6);
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << "                       Mole Fraction   Mass Fraction                 " << endl;
	fSummary << "---------------------------------------------------------------------" << endl;
	for(j=1;j<=x_elements.Size();j++)
		fSummary << " " << setw(19) << left  << gas->list_of_elements[j-1] 
				        << setw(16) << right << scientific << x_elements[j] 
						<< setw(16) << right << scientific << omega_elements[j] 
						<< endl;
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << endl;

	fSummary << setprecision(6);
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << " Mole fractions          Equilibrium         Initial                 " << endl;
	fSummary << "---------------------------------------------------------------------" << endl;
	for(j=1;j<=gas->NumberOfSpecies();j++)
		fSummary << " " << setw(19) << left  << gas->names[j] 
				        << setw(16) << right << scientific << x_E[j] 
						<< setw(16) << right << scientific << x_0[j] 
						<< endl;
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << endl;

	fSummary << setprecision(6);
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << " Mass fractions          Equilibrium         Initial                 " << endl;
	fSummary << "---------------------------------------------------------------------" << endl;
	for(j=1;j<=gas->NumberOfSpecies();j++)
		fSummary << " " << setw(19) << left  << gas->names[j] 
				        << setw(16) << right << scientific << omega_E[j] 
						<< setw(16) << right << scientific << omega_0[j] 
						<< endl;
	fSummary << "---------------------------------------------------------------------" << endl;
	fSummary << endl;
	
	fSummary.close();
}

void OpenSMOKE_EquilibriumModule::GetInitialConditions()
{
	// 1. Temperature
		 if (iSetTemperature == true && iSetPressure == true)
		 {	/* Nothing to do */	 }

	else if (iSetTemperature == true && iSetDensity == true)
	{	P_0_Pascal = rho_0*Constants::R_J_kmol*T_0/MWmix_0;	}

	else if (iSetTemperature == true && iSetVolume == true)
	{	P_0_Pascal = Constants::R_J_kmol*T_0/MWmix_0/V_0;	}

	else if (iSetTemperature == true && iSetInternalEnergy == true)
	{	ErrorMessage("Temperature and Internal Energy are mutually exclusive...");	}

	else if (iSetTemperature == true && iSetEnthalpy == true)
	{	ErrorMessage("Temperature and Enthalpy are mutually exclusive...");	}

	else if (iSetTemperature == true && iSetEntropy == true)
	{	P_0_Pascal = gas->GetPressureFromMassEntropyAndMoleFractions(T_0, S_0, x_0);	}


	// 2. Pressure
	else if (iSetPressure == true && iSetDensity == true)
	{	T_0 = P_0_Pascal*MWmix_0/Constants::R_J_kmol/rho_0;}

	else if (iSetPressure == true && iSetVolume == true)
	{	T_0 = P_0_Pascal*MWmix_0/Constants::R_J_kmol*V_0;}

	else if (iSetPressure == true && iSetInternalEnergy == true)
	{	T_0 = gas->GetTemperatureFromMassInternalEnergyAndMoleFractions(Constants::T_Reference, U_0, x_0);}

	else if (iSetPressure == true && iSetEnthalpy == true)
	{	T_0 = gas->GetTemperatureFromMassEnthalpyAndMoleFractions(Constants::T_Reference, H_0, x_0);}

	else if (iSetPressure == true && iSetEntropy == true)
	{	T_0 = gas->GetTemperatureFromMassEntropyAndMoleFractions(Constants::T_Reference, P_0_Pascal, S_0, x_0);	}


	// 3. Density
	else if (iSetDensity == true && iSetVolume == true)
	{	ErrorMessage("Density and Specific Volume are mutually exclusive...");	}

	else if (iSetDensity == true && iSetInternalEnergy == true)
	{
		T_0 = gas->GetTemperatureFromMassInternalEnergyAndMoleFractions(Constants::T_Reference, U_0, x_0);
		P_0_Pascal = rho_0*Constants::R_J_kmol*T_0/MWmix_0;
	}

	else if (iSetDensity == true && iSetEnthalpy == true)
	{
		T_0 = gas->GetTemperatureFromMassEnthalpyAndMoleFractions(Constants::T_Reference, H_0, x_0);
		P_0_Pascal = rho_0*Constants::R_J_kmol*T_0/MWmix_0;
	}

	else if (iSetDensity == true && iSetEntropy == true)
	{	ErrorMessage("Density and Entropy are mutually exclusive...");	}


	// 4. Volume
	else if (iSetVolume == true && iSetInternalEnergy == true)
	{
		T_0 = gas->GetTemperatureFromMassInternalEnergyAndMoleFractions(Constants::T_Reference, U_0, x_0);
		P_0_Pascal = Constants::R_J_kmol*T_0/MWmix_0/V_0;
	}

	else if (iSetVolume == true && iSetEnthalpy == true)
	{
		T_0 = gas->GetTemperatureFromMassEnthalpyAndMoleFractions(Constants::T_Reference, H_0, x_0);
		P_0_Pascal = Constants::R_J_kmol*T_0/MWmix_0/V_0;
	}

	else if (iSetVolume == true && iSetEntropy == true)
	{	ErrorMessage("Specific Volume and Entropy are mutually exclusive...");	}


	// 5. Internal energy
	else if (iSetInternalEnergy == true && iSetEnthalpy == true)
	{	ErrorMessage("Internal Energy and Enthalpy are mutually exclusive...");	}

	else if (iSetInternalEnergy == true && iSetEntropy == true)
	{
		T_0 = gas->GetTemperatureFromMassInternalEnergyAndMoleFractions(Constants::T_Reference, U_0, x_0);
		P_0_Pascal = gas->GetPressureFromMassEntropyAndMoleFractions(T_0, S_0, x_0);
	}		


	// 6. Enthalpy
	else if (iSetEnthalpy == true && iSetEntropy == true)
	{
		T_0 = gas->GetTemperatureFromMassEnthalpyAndMoleFractions(Constants::T_Reference, H_0, x_0);
		P_0_Pascal = gas->GetPressureFromMassEntropyAndMoleFractions(T_0, S_0, x_0);
	}		

	else 
		ErrorMessage("The specified couple of thermodynamic functions cannot be used for this equilibrium problem... ");
}