/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci								   *
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

#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "basic/OpenSMOKE_Conversions.h"
#include "basic/OpenSMOKE_WarningFile.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_UnitsData.h"

void OpenSMOKE_CHEMKINInterpreter_UnitsData::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_UnitsData"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_UnitsData::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_UnitsData"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_CHEMKINInterpreter_UnitsData::Setup(ofstream *_fLog, OpenSMOKE_WarningFile *_fWarning)
{
	fLog		= _fLog;
	fWarning	= _fWarning; 
}

OpenSMOKE_CHEMKINInterpreter_UnitsData::OpenSMOKE_CHEMKINInterpreter_UnitsData()
{
	name_object = "Units parser";

	energy_units.push_back("Allowed energy units");
	frequency_factor_units.push_back("Allowed frequency factor units");
	energy_units_original.push_back("Allowed energy units - Original");
	frequency_factor_units_original.push_back("Allowed frequency factor units - Original");

	energy_units.push_back("CAL#MOLE");		energy_units_original.push_back("CAL-MOLE");
	energy_units.push_back("KCAL#MOLE");	energy_units_original.push_back("KCAL-MOLE");
	energy_units.push_back("JOULES#MOLE");	energy_units_original.push_back("JOULES-MOLE");
	energy_units.push_back("KJOULES#MOLE");	energy_units_original.push_back("KJOULES-MOLE");
	energy_units.push_back("KELVINS");		energy_units_original.push_back("KELVINS");
	energy_units.push_back("EVOLTS");		energy_units_original.push_back("EVOLTS");

	energy_units.push_back("CAL#");			energy_units_original.push_back("CAL-");
	energy_units.push_back("KCAL");			energy_units_original.push_back("KCAL");
	energy_units.push_back("JOUL");			energy_units_original.push_back("JOUL");
	energy_units.push_back("KJOU");			energy_units_original.push_back("KJOU");
	energy_units.push_back("KELV");			energy_units_original.push_back("KELV");
	energy_units.push_back("EVOL");			energy_units_original.push_back("EVOL");


	frequency_factor_units.push_back("MOLES");		frequency_factor_units_original.push_back("MOLES");
	frequency_factor_units.push_back("MOLECULES");	frequency_factor_units_original.push_back("MOLECULES");

	frequency_factor_units.push_back("MOLE");		frequency_factor_units_original.push_back("MOLE");
	frequency_factor_units.push_back("MOLEC");		frequency_factor_units_original.push_back("MOLEC");

	frequency_factor_units.push_back("MOLE#CM3");		frequency_factor_units_original.push_back("MOLE-CM3");
	frequency_factor_units.push_back("MOLECULE#CM3");	frequency_factor_units_original.push_back("MOLECULE-CM3");

	frequency_factor_units.push_back("OPENSMOKE");	frequency_factor_units_original.push_back("OPENSMOKE");

	current_energy				= "CAL#MOLE";
	current_frequency_factor	= "MOLES";

	iCurrentEnergy				= false;
	iCurrentFrequencyFactor		= false;
}

OpenSMOKE_CHEMKINInterpreter_UnitsData::~OpenSMOKE_CHEMKINInterpreter_UnitsData()
{

}


void OpenSMOKE_CHEMKINInterpreter_UnitsData::Parse_Units(const string units_string)
{
	int i;
	bool success = false;
	for(i=1;i<=int(energy_units.size())-1;i++)
		if (caseInsCompare(units_string, energy_units[i]) == true)
		{
			if (iCurrentEnergy == true)	ErrorMessage("Activation energy UNITS for all the reactions are assigned more than once");
			iCurrentEnergy = true;
			current_energy = energy_units[i];
			success = true;
			break;
		}
	
	for(i=1;i<=int(frequency_factor_units.size())-1;i++)
		if (caseInsCompare(units_string, frequency_factor_units[i]) == true)
		{
			if (iCurrentFrequencyFactor == true)	ErrorMessage("Frequency factor UNITS for all the reactions are assigned more than once");
			iCurrentFrequencyFactor  = true;
			current_frequency_factor = frequency_factor_units[i];
			success = true;
			break;
		}

	if (success == false)
	{
		string units_modified = units_string;
		StringSubstitutionAll(units_modified, "#",  "-");
		ErrorMessage("The " + units_modified + " units are not allowed");
	}
}

void OpenSMOKE_CHEMKINInterpreter_UnitsData::GiveMeConversionFactor(const string units_string, double &conversion, bool &iEnergy)
{
	if		(units_string == "CAL#MOLE")		{ conversion = 1.e0; iEnergy=true;}
	else if (units_string == "KCAL#MOLE")		{ conversion = 1.e3; iEnergy=true;}
	else if (units_string == "JOULES#MOLE")		{ conversion = 1.e0/OpenSMOKE_Conversions::J_from_cal; iEnergy=true;}
	else if (units_string == "KJOULES#MOLE")	{ conversion = 1.e3/OpenSMOKE_Conversions::J_from_cal; iEnergy=true;}
	else if (units_string == "KELVINS")			{ conversion = Constants::R_cal_mol; iEnergy=true;}
	else if (units_string == "EVOLTS")			{ conversion = Constants::Nav_mol*OpenSMOKE_Conversions::J_from_eV/OpenSMOKE_Conversions::J_from_cal; iEnergy=true;}
	else if (units_string == "CAL#")			{ conversion = 1.e0; iEnergy=true;}
	else if (units_string == "KCAL")			{ conversion = 1.e3; iEnergy=true;}
	else if (units_string == "JOUL")			{ conversion = 1.e0/OpenSMOKE_Conversions::J_from_cal; iEnergy=true;}
	else if (units_string == "KJOU")			{ conversion = 1.e3/OpenSMOKE_Conversions::J_from_cal; iEnergy=true;}
	else if (units_string == "KELV")			{ conversion = Constants::R_cal_mol; iEnergy=true;}
	else if (units_string == "EVOL")			{ conversion = Constants::Nav_mol*OpenSMOKE_Conversions::J_from_eV/OpenSMOKE_Conversions::J_from_cal; iEnergy=true; iEnergy=true;}
	
	else if	(units_string == "MOLE")			{ conversion = 1.e0; iEnergy=false;}
	else if	(units_string == "MOLES")			{ conversion = 1.e0; iEnergy=false;}
	else if	(units_string == "MOLE#CM3")		{ conversion = 1.e0; iEnergy=false;}
	else if	(units_string == "MOLECULES")		{ conversion = Constants::Nav_mol; iEnergy=false;}
	else if	(units_string == "MOLEC")			{ conversion = Constants::Nav_mol; iEnergy=false;}
	else if	(units_string == "MOLECULE#CM3")	{ conversion = Constants::Nav_mol; iEnergy=false;}

	else if	(units_string == "OPENSMOKE")		{ conversion = 1.e3; iEnergy=false;}

	else ErrorMessage("The " + units_string + " units are not allowed");
}

void OpenSMOKE_CHEMKINInterpreter_UnitsData::Summary()
{
	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << "                           Units                                 " << endl;
	*fLog << " ----------------------------------------------------------------" << endl;
	*fLog << " Energy:           " << current_energy				<< endl;
	*fLog << " Frequency factor: " << current_frequency_factor	<< endl;
	*fLog << endl;
}