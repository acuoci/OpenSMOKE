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
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportSpecies.h"
#include <sstream>

void OpenSMOKE_CHEMKINInterpreter_TransportSpecies::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_TransportSpecies"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Line:   " << index_line        << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_TransportSpecies::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_TransportSpecies"	<< endl;
    cout << "Object:   "	<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
	cout << "Line:     "	<< index_line   << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}
OpenSMOKE_CHEMKINInterpreter_TransportSpecies::OpenSMOKE_CHEMKINInterpreter_TransportSpecies()
{
}

OpenSMOKE_CHEMKINInterpreter_TransportSpecies::~OpenSMOKE_CHEMKINInterpreter_TransportSpecies()
{
}

void OpenSMOKE_CHEMKINInterpreter_TransportSpecies::ReadMainData(const std::string line, const int iLine)
{
	vector<string> instructions;
	instructions.push_back("instructions");
		
	std::string dummy;
	stringstream parsed_string(line);

	for(;;)
	{
		parsed_string >> dummy;
		if (parsed_string.fail())
			break;
		instructions.push_back(dummy);
	}

	index_line = iLine;
	if (instructions.size()-1 != 7)
		ErrorMessage("Syntax error in transport properties definition");

	name			= instructions[1].c_str();
	shape_factor	= atoi(instructions[2].c_str());
	epsylon_over_kb	= atof(instructions[3].c_str());
	sigma			= atof(instructions[4].c_str());
	mu				= atof(instructions[5].c_str())*1.e-6;
	alfa			= atof(instructions[6].c_str());
	zRot298			= atof(instructions[7].c_str());			
}


void OpenSMOKE_CHEMKINInterpreter_TransportSpecies::Analyze()
{
/*	vector<string> instructions;
	SeparateInstructions(elements, instructions);
	element_names.push_back("names");
	element_numbers.push_back(0.);
	for(int i=1;i<=instructions.size()-1;i++)
	{
		std::string name;
		double number;
		int flag = SeparateNumberFromStringForElements(instructions[i], name, number);
		if(flag == 1)	element_numbers.push_back(number);
		if(flag == 2)	element_names.push_back(name);
		if(flag == 3)	{ element_numbers.push_back(number);element_names.push_back(name);}
	}

	if (element_numbers.size() != element_names.size())
		ErrorMessage("Syntax error in elemental composition definition");*/
}