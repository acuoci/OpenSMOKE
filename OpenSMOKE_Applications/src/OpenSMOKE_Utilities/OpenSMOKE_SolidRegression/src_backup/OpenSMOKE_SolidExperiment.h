/***************************************************************************
 *   Copyright (C) 2009 by Tiziano Maffei and Alberto Cuoci  	           *
 *   tiziano.maffei@chem.polimi.it   				                       *
 *   alberto.cuoci@polimi.it   						                       *
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
#ifndef OPENSMOKE_SOLIDEXPERIMENT
#define OPENSMOKE_SOLIDEXPERIMENT

#include "BzzMath.hpp"


class OpenSMOKE_SolidExperiment
{
public:

	string name_of_file;

	double temperature;
	double Sg0;
	double Sigma;
	double pressure;
	vector<string>	gas_names;
	vector<double>	gas_mole_fraction;
	vector<int>		gas_index;

	BzzVector		gas_c;


	int nPoints;
	BzzVector x;
	BzzVector y;

	

public:

	void PrepareConcentrations(const vector<string> list_of_gas_names);
	void ReadFromFile(const string filename);
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_SOLIDEXPERIMENT