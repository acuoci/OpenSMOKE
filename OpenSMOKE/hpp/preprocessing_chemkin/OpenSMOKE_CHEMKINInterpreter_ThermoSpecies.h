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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_THERMOSPECIES_H)
#define OPENSMOKE_CHEMKININTERPRETER_THERMOSPECIES_H

#include "BzzMath.hpp"
#include <vector>

class OpenSMOKE_CHEMKINInterpreter_ThermoSpecies  
{
public:
	OpenSMOKE_CHEMKINInterpreter_ThermoSpecies();
	virtual ~OpenSMOKE_CHEMKINInterpreter_ThermoSpecies();

	bool iContinuation;
	int index_line;
	
	string name;
	string date;
	string elements;
	string phase;
	string low;
	string high;
	string mean;
	string index;


	BzzVector upper;
	BzzVector lower;
	BzzVector element_indices;
	BzzVector isotope_indices;

	double mw;
	double t_min;
	double t_mean;
	double t_max;

	vector<string> element_names;
	vector<double> element_numbers;

	void AssignTemperatures(const double _t_min, const double _t_mean, const double _t_max);
	void ReadMainData(const string line, const int  iLine);
	void ReadFirstLine(const string line);
	void ReadSecondLine(const string line);
	void ReadThirdLine(const string line);
	void ReadAdditionalLine(const string line);
	void Analyze();

private:
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

class OpenSMOKE_CHEMKINInterpreter_SurfaceThermoSpecies : OpenSMOKE_CHEMKINInterpreter_ThermoSpecies
{
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_THERMOSPECIES_H)
