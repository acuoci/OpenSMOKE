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
#include "liquid/OpenSMOKE_LiquidCriticalConstants_Dictionary.h"

void OpenSMOKE_LiquidCriticalConstants_Dictionary::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidCriticalConstants_Dictionary"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LiquidCriticalConstants_Dictionary::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidCriticalConstants_Dictionary"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

void OpenSMOKE_LiquidCriticalConstants_Dictionary::ErrorMessage(const string message, const int iLine)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LiquidCriticalConstants_Dictionary"	<< endl;
    cout << "Line:   " << iLine					<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

OpenSMOKE_LiquidCriticalConstants_Dictionary::OpenSMOKE_LiquidCriticalConstants_Dictionary()
{
	name_object = "[not assigned]";
}

void OpenSMOKE_LiquidCriticalConstants_Dictionary::ReadFromFile(const string file_name)
{
	vector<string> lines;
	BzzVectorInt indexLines;

	PrepareLiquidPropertiesDictionary("Liquid Critical Constants Database", file_name, lines, indexLines);

	// Memory allocation
	int total_number_of_species = indexLines.Size();
	ChangeDimensions(total_number_of_species, &CAS);
	ChangeDimensions(total_number_of_species, &MW);
	ChangeDimensions(total_number_of_species, &Tc);
	ChangeDimensions(total_number_of_species, &Pc);
	ChangeDimensions(total_number_of_species, &Vc);
	ChangeDimensions(total_number_of_species, &Zc);
	ChangeDimensions(total_number_of_species, &omega);

	name_extended.resize(total_number_of_species+1);
	name_first.resize(total_number_of_species+1);
	name_second.resize(total_number_of_species+1);

	for(int i=1;i<=total_number_of_species;i++)
		ProcessData(lines[indexLines[i]], indexLines[i], i);
}

void OpenSMOKE_LiquidCriticalConstants_Dictionary::ProcessData(const string line, const int iLine, const int index)
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

	if (instructions.size()-1 != 11)
		ErrorMessage("Syntax error in critical constants definition", iLine);

	name_extended[index]	= instructions[2];
	name_first[index]		= instructions[3];
	name_second[index]		= instructions[4];
	CAS[index]				= atoi(instructions[5].c_str());
	MW[index]				= atof(instructions[6].c_str());			// [kg/kmol]
	Tc[index]				= atof(instructions[7].c_str());			// [K]
	Pc[index]				= atof(instructions[8].c_str())*1.e6;		// [Pa]
	Vc[index]				= atof(instructions[9].c_str());			// [m3/kmol]
	Zc[index]				= atof(instructions[10].c_str());			// [-]
	omega[index]			= atof(instructions[11].c_str());			// [-]
}


void OpenSMOKE_LiquidCriticalConstants_Dictionary::WriteToFile(const string file_name)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);

	fOutput << setw(5)  << left << "!No.";
	fOutput << setw(40) << left << "Name";
	fOutput << setw(16) << left << "Formula";
	fOutput << setw(16) << left << "Formula(II)";
	fOutput << setw(16) << left << "CAS";
	fOutput << setw(16) << left << "MW";
	fOutput << setw(16) << left << "Tc[K]";
	fOutput << setw(16) << left << "Pc*1e-6[Pa]";
	fOutput << setw(16) << left << "Vc[m3/kmol]";
	fOutput << setw(16) << left << "Zc";
	fOutput << setw(16) << left << "omega";
	fOutput << endl;

	for(int i=1;i<=CAS.Size();i++)
	{
		fOutput << setw(5)  << left << i;
		fOutput << setw(40) << left << name_extended[i];
		fOutput << setw(16) << left << name_first[i];
		fOutput << setw(16) << left << name_second[i];
		fOutput << setw(16) << left << CAS[i];
		fOutput << setw(16) << left << scientific << MW[i];
		fOutput << setw(16) << left << scientific << Tc[i];
		fOutput << setw(16) << left << scientific << Pc[i]*1e-6;
		fOutput << setw(16) << left << scientific << Vc[i];
		fOutput << setw(16) << left << scientific << Zc[i];
		fOutput << setw(16) << left << scientific << omega[i];
		fOutput << endl;
	}
}

void OpenSMOKE_LiquidCriticalConstants_Dictionary::SaveToFile(BzzSave &fSave)
{
	string dummy= "CONSTANTS";
	char name[Constants::NAME_SIZE];

	strcpy(name, dummy.c_str());
	fSave.fileSave.write((char*) name, sizeof(name));

	// Number of species
	fSave << CAS.Size();

	// Name extended
	for(int i=1;i<=CAS.Size();i++)
	{
		strcpy(name, name_extended[i].c_str());
		fSave.fileSave.write((char*) name, sizeof(name));
	}	

	// Formula
	for(int i=1;i<=CAS.Size();i++)
	{
		strcpy(name, name_first[i].c_str());
		fSave.fileSave.write((char*) name, sizeof(name));
	}	

	// Formula(II)
	for(int i=1;i<=CAS.Size();i++)
	{
		strcpy(name, name_second[i].c_str());
		fSave.fileSave.write((char*) name, sizeof(name));
	}	

	fSave << CAS;
	fSave << MW;
	fSave << Tc;
	fSave << Pc;
	fSave << Vc;
	fSave << Zc;
	fSave << omega;
}