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

#if !defined(OPENSMOKE_LIQUIDSPECIFICHEAT_DICTIONARY_H)
#define OPENSMOKE_LIQUIDSPECIFICHEAT_DICTIONARY_H

#include "OpenSMOKE.hpp"

class OpenSMOKE_LiquidSpecificHeat_Dictionary
{

friend class OpenSMOKE_LiquidProperties_Database;

public:

	OpenSMOKE_LiquidSpecificHeat_Dictionary();
	void ReadFromFile(const std::string file_name);
	void WriteToFile(const std::string file_name);
	void SaveToFile(BzzSave &fSave);

private:

	BzzVectorInt CAS;

	BzzVector MW;
	BzzVector C1;
	BzzVector C2;
	BzzVector C3;
	BzzVector C4;
	BzzVector C5;
	BzzVector Tmin;
	BzzVector Tmax;
	BzzVector Cpmin;
	BzzVector Cpmax;

	vector<string> name_extended;
	vector<string> name_first;
	vector<string> name_second;
	vector<liquid_specificheat_equation::equation> equation;

	void ProcessData(const std::string line, const int iLine, const int index);

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void ErrorMessage(const std::string message, const int iLine);
	void WarningMessage(const std::string message);
};

#endif // OPENSMOKE_LIQUIDSPECIFICHEAT_DICTIONARY_H