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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_UNITSDATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_UNITSDATA_H

#include <vector>
#include "BzzMath.hpp"

class OpenSMOKE_WarningFile;

class OpenSMOKE_CHEMKINInterpreter_UnitsData  
{
public:
	OpenSMOKE_CHEMKINInterpreter_UnitsData();
	virtual ~OpenSMOKE_CHEMKINInterpreter_UnitsData();

	vector<string> energy_units;
	vector<string> frequency_factor_units;
	vector<string> energy_units_original;
	vector<string> frequency_factor_units_original;

	string current_energy;
	string current_frequency_factor;

	bool iCurrentEnergy;
	bool iCurrentFrequencyFactor;

	void Parse_Units(const string units_string);
	void GiveMeConversionFactor(const string units_string, double &conversion, bool &iEnergy);
	void Setup(ofstream *_fLog);
	void Setup(ofstream *_fLog, OpenSMOKE_WarningFile *_fWarning);
	void Summary();

	ofstream *fLog;
	OpenSMOKE_WarningFile *fWarning;

private:
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_UNITSDATA_H)
