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

#if !defined(OPENSMOKE_LIQUIDDENSITY_DICTIONARY_H)
#define OPENSMOKE_LIQUIDDENSITY_DICTIONARY_H

#include "OpenSMOKE.hpp"

class OpenSMOKE_LiquidDensity_Dictionary
{

friend class OpenSMOKE_LiquidProperties_Database;

public:

	OpenSMOKE_LiquidDensity_Dictionary();
	void ReadFromFile(const string file_name);
	void WriteToFile(const string file_name);
	void SaveToFile(BzzSave &fSave);

private:

	BzzVectorInt CAS;

	BzzVector MW;
	BzzVector C1;
	BzzVector C2;
	BzzVector C3;
	BzzVector C4;
	BzzVector Tmin;
	BzzVector Tmax;
	BzzVector Rhomin;
	BzzVector Rhomax;

	vector<string> name_extended;
	vector<string> name_first;
	vector<string> name_second;
	vector<liquid_density_equation::equation> equation;

	void ProcessData(const string line, const int iLine, const int index);

private:

	string name_object;
	void ErrorMessage(const string message);
	void ErrorMessage(const string message, const int iLine);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_LIQUIDDENSITY_DICTIONARY_H

/*	
231	Water	H2O	H2O	7732185	18.015	4.9684	0.2779	647.1300	0.1874	333.15	54.703	403.15	51.934	7
231	Water	H2O	H2O	7732185	18.015	4.3580	0.2487	647.1300	0.2534	403.15	51.934	647.15	20.576	7
*/