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
#include <iomanip>
#include "liquid/OpenSMOKE_LiquidThermalConductivity_Dictionary.h"

void OpenSMOKE_LiquidThermalConductivity_Dictionary::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidThermalConductivity_Dictionary"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LiquidThermalConductivity_Dictionary::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidThermalConductivity_Dictionary"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

void OpenSMOKE_LiquidThermalConductivity_Dictionary::ErrorMessage(const string message, const int iLine)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidThermalConductivity_Dictionary"	<< endl;
    cout << "Line:   " << iLine					<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

OpenSMOKE_LiquidThermalConductivity_Dictionary::OpenSMOKE_LiquidThermalConductivity_Dictionary()
{
	name_object = "[not assigned]";
}

void OpenSMOKE_LiquidThermalConductivity_Dictionary::ReadFromFile(const string file_name)
{
	vector<string> lines;
	BzzVectorInt indexLines;

	PrepareLiquidPropertiesDictionary("Liquid Thermal Conductivity Database", file_name, lines, indexLines);

	// Memory allocation
	int total_number_of_species = indexLines.Size();
	ChangeDimensions(total_number_of_species, &CAS);
	ChangeDimensions(total_number_of_species, &MW);
	ChangeDimensions(total_number_of_species, &C1);
	ChangeDimensions(total_number_of_species, &C2);
	ChangeDimensions(total_number_of_species, &C3);
	ChangeDimensions(total_number_of_species, &C4);
	ChangeDimensions(total_number_of_species, &C5);
	ChangeDimensions(total_number_of_species, &Tmin);
	ChangeDimensions(total_number_of_species, &Tmax);
	ChangeDimensions(total_number_of_species, &Lambdamin);
	ChangeDimensions(total_number_of_species, &Lambdamax);
	name_extended.resize(total_number_of_species+1);
	name_first.resize(total_number_of_species+1);
	name_second.resize(total_number_of_species+1);
	equation.resize(total_number_of_species+1);

	for(int i=1;i<=total_number_of_species;i++)
		ProcessData(lines[indexLines[i]], indexLines[i], i);
}

void OpenSMOKE_LiquidThermalConductivity_Dictionary::ProcessData(const string line, const int iLine, const int index)
{
	vector<string> instructions;
	instructions.push_back("instructions");
		
	string dummy;
	stringstream parsed_string(line);

	for(;;)
	{
		parsed_string >> dummy;
		if (parsed_string.fail())
			break;
		instructions.push_back(dummy);
	}

	if (instructions.size()-1 != 16)
		ErrorMessage("Syntax error in thermal conductivity coefficient definition", iLine);

	name_extended[index]	= instructions[2];
	name_first[index]		= instructions[3];
	name_second[index]		= instructions[4];
	CAS[index]				= atoi(instructions[5].c_str());
	MW[index]				= atof(instructions[6].c_str());
	C1[index]				= atof(instructions[7].c_str());
	C2[index]				= atof(instructions[8].c_str());
	C3[index]				= atof(instructions[9].c_str());
	C4[index]				= atof(instructions[10].c_str());
	C5[index]				= atof(instructions[11].c_str());
	Tmin[index]				= atof(instructions[12].c_str());
	Lambdamin[index]		= atof(instructions[13].c_str());
	Tmax[index]				= atof(instructions[14].c_str());
	Lambdamax[index]		= atof(instructions[15].c_str());

	int dummy_int			= atoi(instructions[16].c_str());
	
	if (dummy_int == 0)	equation[index] = liquid_thermalconductivity_equation::EQ0;
	else if (dummy_int == 1)	equation[index] = liquid_thermalconductivity_equation::EQ1;
	else if (dummy_int == 2)	equation[index] = liquid_thermalconductivity_equation::EQ2;
	else ErrorMessage("Unknown equation...", iLine);
}

void OpenSMOKE_LiquidThermalConductivity_Dictionary::WriteToFile(const string file_name)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);

	fOutput << setw(5)  << left << "!No.";
	fOutput << setw(40) << left << "Name";
	fOutput << setw(16) << left << "Formula";
	fOutput << setw(16) << left << "Formula(II)";
	fOutput << setw(16) << left << "CAS";
	fOutput << setw(16) << left << "MW";
	fOutput << setw(16) << left << "C1";
	fOutput << setw(16) << left << "C2";
	fOutput << setw(16) << left << "C3";
	fOutput << setw(16) << left << "C4";
	fOutput << setw(16) << left << "C5";
	fOutput << setw(16) << left << "Tmin";
	fOutput << setw(16) << left << "Lambdamin";
	fOutput << setw(16) << left << "Tmax";
	fOutput << setw(16) << left << "Lambdamax";
	fOutput << setw(4)  << left << "Eq";
	fOutput << endl;

	for(int i=1;i<=CAS.Size();i++)
	{
		fOutput << setw(5)  << left << i;
		fOutput << setw(40) << left << name_extended[i];
		fOutput << setw(16) << left << name_first[i];
		fOutput << setw(16) << left << name_second[i];
		fOutput << setw(16) << left << CAS[i];
		fOutput << setw(16) << left << scientific << MW[i];
		fOutput << setw(16) << left << scientific << C1[i];
		fOutput << setw(16) << left << scientific << C2[i];
		fOutput << setw(16) << left << scientific << C3[i];
		fOutput << setw(16) << left << scientific << C4[i];
		fOutput << setw(16) << left << scientific << C5[i];
		fOutput << setw(16) << left << scientific << Tmin[i];
		fOutput << setw(16) << left << scientific << Lambdamin[i];
		fOutput << setw(16) << left << scientific << Tmax[i];
		fOutput << setw(16) << left << scientific << Lambdamax[i];

			 if (equation[i] == liquid_thermalconductivity_equation::EQ0)	fOutput << setw(4) << left << fixed << 0;
		else if (equation[i] == liquid_thermalconductivity_equation::EQ1)	fOutput << setw(4) << left << fixed << 1;
		else if (equation[i] == liquid_thermalconductivity_equation::EQ2)	fOutput << setw(4) << left << fixed << 2;

		fOutput << endl;
	}
}

void OpenSMOKE_LiquidThermalConductivity_Dictionary::SaveToFile(BzzSave &fSave)
{
	string dummy= "THERMALCOND";
	char name[Constants::NAME_SIZE];

	BzzVectorInt equation_int(C1.Size());
	for(int i=1;i<=C1.Size();i++)
	{
			 if (equation[i] == liquid_thermalconductivity_equation::EQ0)	equation_int[i] = 0;
		else if (equation[i] == liquid_thermalconductivity_equation::EQ1)	equation_int[i] = 1;
		else if (equation[i] == liquid_thermalconductivity_equation::EQ2)	equation_int[i] = 2;
	}		

	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));

	fSave << C1;
	fSave << C2;
	fSave << C3;
	fSave << C4;
	fSave << C5;
	fSave << Tmin;
	fSave << Tmax;
	fSave << equation_int;
}