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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_SPECIESDATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_SPECIESDATA_H

#include <vector>
#include "BzzMath.hpp"

class OpenSMOKE_WarningFile;

class OpenSMOKE_CHEMKINInterpreter_SpeciesData
{
public:

	OpenSMOKE_CHEMKINInterpreter_SpeciesData();

	bool Parse_Species_Name(const string species_name);
	void Setup(ofstream *_fLog, OpenSMOKE_WarningFile *_fWarning);
	void Summary();

	vector<string>	species_in_list;
	ofstream *fLog;
	OpenSMOKE_WarningFile *fWarning;

private:
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_SPECIESDATA_H)
