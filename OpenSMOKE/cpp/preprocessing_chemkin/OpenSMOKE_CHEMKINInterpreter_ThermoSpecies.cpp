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
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.h"
#include <sstream>

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ThermoSpecies"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Line:   " << index_line        << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ThermoSpecies"	<< endl;
    cout << "Object:   "	<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
	cout << "Line:     "	<< index_line   << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}
OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::OpenSMOKE_CHEMKINInterpreter_ThermoSpecies()
{
	iContinuation = false;
	ChangeDimensions(7, &upper);
	ChangeDimensions(7, &lower);
}

OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::~OpenSMOKE_CHEMKINInterpreter_ThermoSpecies()
{
}

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::AssignTemperatures(const double _t_min, const double _t_mean, const double _t_max)
{
	t_min  = _t_min;
	t_mean = _t_mean;
	t_max  = _t_max;
}

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::ReadMainData(const std::string line, const int iLine)
{
	index_line = iLine;

	int l = line.length();
	if (l<80)	ErrorMessage("The length must be at least 80 chars");
	if (l>120)	ErrorMessage("The length exceeds 120 chars");
	if (line.substr(79, 1) != "1")	ErrorMessage("Integer 1 must appear in column 80");
	
	name = line.substr(0, 16); 
	CleanFromBlanksForThermo(name);				// name (1)
	date = line.substr(18, 6);					// date (2)
	elements = line.substr(24, 20);				// elements (3)
	phase = line.substr(44, 1);					// phase (4)
	low = line.substr(45, 10);					// low temperature (5)
	high = line.substr(55, 10);					// high temperature (6)
	mean = line.substr(65, 8);					// mean temperature (7)
	elements += line.substr(73, 5);				// elements (8)
	index = line.substr(79, 1);					// index (9)
	if (l>80)
	{
		if (line.substr(80, 1) == "&")	iContinuation = true;
		else 
		{
			bool empty = true;
			for(int i=1;i<=l-80;i++)
			{
				if (line[79+i] != ' ' && line[79+i] != '\t')
				{
					empty = false;
					break;
				}
			}

			if (empty == false) elements += line.substr(80, l-80);
		}
	}

	if (CheckForBlankLine(low) == true || CheckForBlankLine(high) == true)
			ErrorMessage("Syntax Error in Defining Thermodynamic Temperature Intervals");

	if (CheckForBlankLine(mean) == true)
	{
		t_min  = atof(low.c_str());
		t_max  = atof(high.c_str());
	}
	else
	{
		t_min  = atof(low.c_str());
		t_mean = atof(mean.c_str());
		t_max  = atof(high.c_str());
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::ReadFirstLine(const std::string line)
{		
	stringstream parsed_string(line);

	int last_char;
	for(int i=1;i<=5;i++)
		parsed_string >> upper[i];
	parsed_string >> last_char;
	if (last_char != 2)	
	{
		cout << "Second line: ***" << line << "***" << endl;
		ErrorMessage("Integer 2 must appear in column 80");
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::ReadSecondLine(const std::string line)
{		
	stringstream parsed_string(line);

	int i;
	int last_char;
	for(i=1;i<=2;i++)
		parsed_string >> upper[5+i];
	for(i=1;i<=3;i++)
		parsed_string >> lower[i];

	parsed_string >> last_char;
	if (last_char != 3)	ErrorMessage("Integer 3 must appear in column 80");
}

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::ReadThirdLine(const std::string line)
{		
	stringstream parsed_string(line);

	int last_char;
	std::string additional;
	for(int i=1;i<=4;i++)
		parsed_string >> lower[3+i];
	parsed_string >> additional; // TODO
	
	if (additional != "4")	
	{	
		parsed_string >> last_char;
		if (last_char != 4)
			ErrorMessage("Integer 4 must appear in column 80");
	}
}

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::ReadAdditionalLine(const std::string line)
{
	elements += " " + line;
}


void stripSpace(std::string &str) 
{
	for (int i=0;i<str.length();i++)
	if (str[i]==' ') 
	{
		str.erase(i,1);
		i--;
	}
}

void ParsingAdditionalElements(vector<string> &subelements, const std::string elements_additional)
{
	vector<int> pattern(elements_additional.length());
	for(int i=0;i<elements_additional.length();i++)
	{
		if (isdigit(elements_additional[i]) == 0)		
		{
			if (elements_additional[i] == ' ')	pattern[i] = -1;	// is a white space
			else								pattern[i] =  0;	// is a letter
		}
		else
			pattern[i]	= 1;							// is a number
	}

	vector<int> change_pattern;
	int k=0;
	do
	{
		if(pattern[k] == 0)
		{
			change_pattern.push_back(k);
			for(int i=k+1;i<pattern.size();i++)
				if (pattern[i] != 0)
				{
					k=i-1;
					break;
				}
		}
		k++;
	}
	while(k<pattern.size());
	change_pattern.push_back(elements_additional.length());

	for(int i=0;i<change_pattern.size()-1;i++)
		subelements.push_back(elements_additional.substr(change_pattern[i], change_pattern[i+1]-change_pattern[i]));
}	

void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::Analyze()
{
	vector<string> subelements;

	subelements.push_back(elements.substr(0, 5));	
	subelements.push_back(elements.substr(5, 5));	
	subelements.push_back(elements.substr(10, 5));	
	subelements.push_back(elements.substr(15, 5));

	if (elements.length() > 20)
	{
		subelements.push_back(elements.substr(20, 5));
		if (elements.length() > 25)
		{
			std::string elements_additional = elements.substr(25, elements.length()-25);
			ParsingAdditionalElements(subelements, elements_additional);
		}
	}

	element_names.push_back("names");
	element_numbers.push_back(0.);
	for(int i=0;i<subelements.size();i++)
	{
		std::string nameelement;
		double number;
		int flag = SeparateNumberFromStringForElements(subelements[i], nameelement, number);
		
		/*if (name == "C2H5O2")
		{
			cout << "Ex: " << subelements[i] << "####" <<endl;
			cout << "Na: " << nameelement << "####" <<endl;
			cout << "Nu: " << number << endl;
		}*/

		if (flag == 1 || flag == 2)	
			ErrorMessage("Syntax error in elemental composition definition");
		else if(flag == 3)	
		{ 
			if (number != 0)
			{
				// Recognize if element already available
				bool iAlreadyAvailable = false;
				for(int j=1;j<=element_names.size()-1;j++)
					if (element_names[j]==nameelement)
					{
						element_numbers[j] += number;
						iAlreadyAvailable = true;
						break;
					}

				if (iAlreadyAvailable == false)
				{
					element_numbers.push_back(number);
					element_names.push_back(nameelement);
				}
			}
		}
	}

	if (element_numbers.size() != element_names.size())
		ErrorMessage("Syntax error in elemental composition definition");
}




/*
void OpenSMOKE_CHEMKINInterpreter_ThermoSpecies::Analyze()
{
	vector<int> pattern(elements.length());
	for(int i=0;i<elements.length();i++)
	{
		if (isdigit(elements[i]) == 0)		
		{
			if (elements[i] == ' ')	pattern[i] = -1;	// is a white space
			else					pattern[i] =  0;	// is a letter
		}
		else
			pattern[i]	= 1;							// is a number
	}

	vector<int> change_pattern;
	int k=0;
	do
	{
		if(pattern[k] == 0)
		{
			change_pattern.push_back(k);
			for(int i=k+1;i<pattern.size();i++)
				if (pattern[i] != 0)
				{
					k=i-1;
					break;
				}
		}
		k++;
	}
	while(k<pattern.size());
	change_pattern.push_back(elements.length());

	vector<string> subelements;
	for(int i=0;i<change_pattern.size()-1;i++)
		subelements.push_back(elements.substr(change_pattern[i], change_pattern[i+1]-change_pattern[i]));
	
	element_names.push_back("names");
	element_numbers.push_back(0.);
	for(int i=0;i<subelements.size();i++)
	{
		std::string name;
		double number;
		int flag = SeparateNumberFromStringForElements(subelements[i], name, number);
		
		//cout << "Ex: " << subelements[i] << "####" <<endl;
		//cout << "Na: " << name << "####" <<endl;
		//cout << "Nu: " << number << endl;
		
		if (flag == 1 || flag == 2)	
			ErrorMessage("Syntax error in elemental composition definition");
		else if(flag == 3)	
		{ 
			if (number != 0)
			{
				// Recognize if element already available
				bool iAlreadyAvailable = false;
				for(int j=1;j<=element_names.size()-1;j++)
					if (element_names[j]==name)
					{
						element_numbers[j] += number;
						iAlreadyAvailable = true;
						break;
					}

				if (iAlreadyAvailable == false)
				{
					element_numbers.push_back(number);
					element_names.push_back(name);
				}
			}
		}
	}

	if (element_numbers.size() != element_names.size())
		ErrorMessage("Syntax error in elemental composition definition");
}

*/