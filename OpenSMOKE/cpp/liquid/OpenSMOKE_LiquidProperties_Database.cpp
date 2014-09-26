/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci								   *
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
#include "liquid/OpenSMOKE_LiquidProperties_Database.h"

void OpenSMOKE_LiquidProperties_Database::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidProperties_Database"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LiquidProperties_Database::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidProperties_Database"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_LiquidProperties_Database::OpenSMOKE_LiquidProperties_Database()
{
	name_object = "[not assigned]";

	dictionary_critical_constants	= new OpenSMOKE_LiquidCriticalConstants_Dictionary();
	dictionary_density				= new OpenSMOKE_LiquidDensity_Dictionary();
	dictionary_vaporization_heat	= new OpenSMOKE_LiquidVaporizationHeat_Dictionary();
	dictionary_vapor_pressure		= new OpenSMOKE_LiquidVaporPressure_Dictionary();
	dictionary_specific_heat		= new OpenSMOKE_LiquidSpecificHeat_Dictionary();
	dictionary_thermal_conductivity	= new OpenSMOKE_LiquidThermalConductivity_Dictionary();
}

void OpenSMOKE_LiquidProperties_Database::ReadFromFolder(const std::string folder_name)
{
	dictionary_critical_constants->ReadFromFile(folder_name + "/CriticalConstants.liq");
	dictionary_density->ReadFromFile(folder_name + "/Density.liq");
	dictionary_vaporization_heat->ReadFromFile(folder_name + "/VaporizationHeat.liq");
	dictionary_vapor_pressure->ReadFromFile(folder_name + "/VaporPressure.liq");
	dictionary_specific_heat->ReadFromFile(folder_name + "/SpecificHeat.liq");
	dictionary_thermal_conductivity->ReadFromFile(folder_name + "/ThermalConductivity.liq");

	CheckConsistency();

	dictionary_density->WriteToFile(folder_name + "/Density.liq.new");
	dictionary_critical_constants->WriteToFile(folder_name + "/CriticalConstants.liq.new");
	dictionary_vaporization_heat->WriteToFile(folder_name + "/VaporizationHeat.liq.new");
	dictionary_vapor_pressure->WriteToFile(folder_name + "/VaporPressure.liq.new");
	dictionary_specific_heat->WriteToFile(folder_name + "/SpecificHeat.liq.new");
	dictionary_thermal_conductivity->WriteToFile(folder_name + "/ThermalConductivity.liq.new");

	SaveOnFile(folder_name + "/Database.liq");
	LoadFromFile(folder_name + "/Database.liq");
}

void OpenSMOKE_LiquidProperties_Database::CheckConsistency()
{
	if (dictionary_critical_constants->CAS.Size() != dictionary_density->CAS.Size() ||
		dictionary_critical_constants->CAS.Size() != dictionary_vaporization_heat->CAS.Size() ||
		dictionary_critical_constants->CAS.Size() != dictionary_vapor_pressure->CAS.Size() ||
		dictionary_critical_constants->CAS.Size() != dictionary_thermal_conductivity->CAS.Size() ||
		dictionary_critical_constants->CAS.Size() != dictionary_specific_heat->CAS.Size() )
		ErrorMessage("LiquidProperties Dictionaries have inconsistent number of species...");

	for(int i=1;i<=dictionary_critical_constants->CAS.Size();i++)
	{
		stringstream index;
		index << i;

		if (dictionary_critical_constants->CAS[i] != dictionary_density->CAS[i] ||
			dictionary_critical_constants->CAS[i] != dictionary_vaporization_heat->CAS[i] ||
			dictionary_critical_constants->CAS[i] != dictionary_vapor_pressure->CAS[i] ||
			dictionary_critical_constants->CAS[i] != dictionary_thermal_conductivity->CAS[i] ||
			dictionary_critical_constants->CAS[i] != dictionary_specific_heat->CAS[i] )
			ErrorMessage("Inconsistent CAS for species: " + index.str());

		if (dictionary_critical_constants->name_first[i] != dictionary_density->name_first[i] ||
			dictionary_critical_constants->name_first[i] != dictionary_vaporization_heat->name_first[i] ||
			dictionary_critical_constants->name_first[i] != dictionary_vapor_pressure->name_first[i] ||
			dictionary_critical_constants->name_first[i] != dictionary_thermal_conductivity->name_first[i] ||
			dictionary_critical_constants->name_first[i] != dictionary_specific_heat->name_first[i] )
			ErrorMessage("Inconsistent name_first for species: " + index.str());

		if (dictionary_critical_constants->name_second[i] != dictionary_density->name_second[i] ||
			dictionary_critical_constants->name_second[i] != dictionary_vaporization_heat->name_second[i] ||
			dictionary_critical_constants->name_second[i] != dictionary_vapor_pressure->name_second[i] ||
			dictionary_critical_constants->name_second[i] != dictionary_thermal_conductivity->name_second[i] ||
			dictionary_critical_constants->name_second[i] != dictionary_specific_heat->name_second[i] )
			ErrorMessage("Inconsistent name_second for species: " + index.str());

		if (dictionary_critical_constants->name_extended[i] != dictionary_density->name_extended[i] ||
			dictionary_critical_constants->name_extended[i] != dictionary_vaporization_heat->name_extended[i] ||
			dictionary_critical_constants->name_extended[i] != dictionary_vapor_pressure->name_extended[i] ||
			dictionary_critical_constants->name_extended[i] != dictionary_thermal_conductivity->name_extended[i] ||
			dictionary_critical_constants->name_extended[i] != dictionary_specific_heat->name_extended[i] )
			ErrorMessage("Inconsistent name_extended for species: " + index.str());
	}
}

void OpenSMOKE_LiquidProperties_Database::SaveOnFile(const std::string file_name)
{
	BzzSave fOutput('*', file_name);
	dictionary_critical_constants->SaveToFile(fOutput);
	dictionary_density->SaveToFile(fOutput);
	dictionary_vapor_pressure->SaveToFile(fOutput);
	dictionary_vaporization_heat->SaveToFile(fOutput);
	dictionary_specific_heat->SaveToFile(fOutput);
	dictionary_thermal_conductivity->SaveToFile(fOutput);
	fOutput.End();
}

void OpenSMOKE_LiquidProperties_Database::LoadFromFile(const std::string file_name)
{
	BzzVectorInt equation_int;
	char dummy[Constants::NAME_SIZE];

	BzzLoad fInput('*', file_name);

	fInput.fileLoad.read((char*) dummy, sizeof(dummy));
	if (strcmp(dummy, "CONSTANTS"))	WarningMessage("Expected: CONSTANTS - Found: " + std::string(dummy));

	// Number of species
	fInput >> N;

	// Memory allocation
	name_extended.resize(1+N);
	name_first.resize(1+N);
	name_second.resize(1+N);
	Rho_equation.resize(1+N);
	Pv_equation.resize(1+N);
	Hv_equation.resize(1+N);
	Cp_equation.resize(1+N);
	Lambda_equation.resize(1+N);

	// Name extended
	for (int i=1;i<=N;i++)
	{
		fInput.fileLoad.read((char*) dummy, sizeof(dummy));
		name_extended[i] = dummy;
	}

	// Formula(I)
	for (int i=1;i<=N;i++)
	{
		fInput.fileLoad.read((char*) dummy, sizeof(dummy));
		name_first[i] = dummy;
	}

	// Formula(II)
	for (int i=1;i<=N;i++)
	{
		fInput.fileLoad.read((char*) dummy, sizeof(dummy));
		name_second[i] = dummy;
	}

	fInput >> CAS;
	fInput >> MW;
	fInput >> Tc;
	fInput >> Pc;
	fInput >> Vc;
	fInput >> Zc;
	fInput >> omega;

	fInput.fileLoad.read((char*) dummy, sizeof(dummy));
	if (strcmp(dummy, "DENSITY"))	WarningMessage("Expected: DENSITY - Found: " + std::string(dummy));
	fInput >> Rho_C1;
	fInput >> Rho_C2;
	fInput >> Rho_C3;
	fInput >> Rho_C4;
	fInput >> Rho_Tmin;
	fInput >> Rho_Tmax;
	fInput >> equation_int;
	for(int i=1;i<=N;i++)
	{
		     if (equation_int[i] == 0)	Rho_equation[i] = liquid_density_equation::EQ0;
		else if (equation_int[i] == 1)	Rho_equation[i] = liquid_density_equation::EQ1;
		else if (equation_int[i] == 7)	Rho_equation[i] = liquid_density_equation::EQ7;
	}

	fInput.fileLoad.read((char*) dummy, sizeof(dummy));
	if (strcmp(dummy, "VAPORPRESSURE"))	WarningMessage("Expected: VAPORPRESSURE - Found: " + std::string(dummy));
	fInput >> Pv_C1;
	fInput >> Pv_C2;
	fInput >> Pv_C3;
	fInput >> Pv_C4;
	fInput >> Pv_C5;
	fInput >> Pv_Tmin;
	fInput >> Pv_Tmax;
	fInput >> equation_int;
	for(int i=1;i<=N;i++)
	     if (equation_int[i] == 1)	Pv_equation[i] = liquid_vaporpressure_equation::EQ1;

	fInput.fileLoad.read((char*) dummy, sizeof(dummy));
	if (strcmp(dummy, "VAPORIZATIONHEAT"))	WarningMessage("Expected: VAPORIZATIONHEAT - Found: " + std::string(dummy));
	fInput >> Hv_C1;
	fInput >> Hv_C2;
	fInput >> Hv_C3;
	fInput >> Hv_C4;
	fInput >> Hv_Tmin;
	fInput >> Hv_Tmax;
	fInput >> equation_int;
	for(int i=1;i<=N;i++)
	     if (equation_int[i] == 1)	Hv_equation[i] = liquid_vaporizationheat_equation::EQ1;

	fInput.fileLoad.read((char*) dummy, sizeof(dummy));
	if (strcmp(dummy, "SPECIFICHEAT"))	WarningMessage("Expected: SPECIFICHEAT - Found: " + std::string(dummy));
	fInput >> Cp_C1;
	fInput >> Cp_C2;
	fInput >> Cp_C3;
	fInput >> Cp_C4;
	fInput >> Cp_C5;
	fInput >> Cp_Tmin;
	fInput >> Cp_Tmax;
	fInput >> equation_int;
	for(int i=1;i<=N;i++)
	{
			 if (equation_int[i] == 1)	Cp_equation[i] = liquid_specificheat_equation::EQ1;
		else if (equation_int[i] == 2)	Cp_equation[i] = liquid_specificheat_equation::EQ2;
	}

	fInput.fileLoad.read((char*) dummy, sizeof(dummy));
	if (strcmp(dummy, "THERMALCOND"))	WarningMessage("Expected: THERMALCOND - Found: " + std::string(dummy));
	fInput >> Lambda_C1;
	fInput >> Lambda_C2;
	fInput >> Lambda_C3;
	fInput >> Lambda_C4;
	fInput >> Lambda_C5;
	fInput >> Lambda_Tmin;
	fInput >> Lambda_Tmax;
	fInput >> equation_int;
	for(int i=1;i<=N;i++)
	{
			 if (equation_int[i] == 0)	Lambda_equation[i] = liquid_thermalconductivity_equation::EQ0;
		else if (equation_int[i] == 1)	Lambda_equation[i] = liquid_thermalconductivity_equation::EQ1;
		else if (equation_int[i] == 2)	Lambda_equation[i] = liquid_thermalconductivity_equation::EQ2;
	}


	fInput.End();
}

int OpenSMOKE_LiquidProperties_Database::RecognizeSpecies(const std::string name)
{
	for(int i=1;i<=N;i++)
		if (name == name_first[i] || name == name_second[i])
			return i;

	ErrorMessage("The following species is not included: " + name);

	return 0;
}