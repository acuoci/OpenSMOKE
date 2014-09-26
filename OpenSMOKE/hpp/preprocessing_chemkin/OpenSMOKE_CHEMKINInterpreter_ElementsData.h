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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_ELEMENTSDATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_ELEMENTSDATA_H

#include <vector>
#include "BzzMath.hpp"

class OpenSMOKE_WarningFile;

class OpenSMOKE_CHEMKINInterpreter_ElementsData
{
public:

	OpenSMOKE_CHEMKINInterpreter_ElementsData();

	int elements_in_database;
	vector<string> elements_name;
	vector<string> elements_description;
	vector<double> elements_mw;

	void Setup(ofstream *_fLog, OpenSMOKE_WarningFile *_fWarning);
	bool Parse_Element_Name(const std::string element_name);
	void Parse_Isotope_Name(const std::string isotope_name, const double mw);
	void Summary();
	void SummaryOnFile(ofstream &fOutput);


	vector<string>	elements_in_list;
	vector<string>	isotope_in_list;
	vector<double>	isotope_mw_in_list;
	vector<int>		elements_in_list_indices;

	ofstream *fLog;
	OpenSMOKE_WarningFile *fWarning;

private:
	
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
	std::string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_ELEMENTSDATA_H)
