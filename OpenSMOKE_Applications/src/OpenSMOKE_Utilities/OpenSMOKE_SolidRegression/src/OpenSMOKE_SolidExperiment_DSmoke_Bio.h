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
#ifndef OPENSMOKE_SOLIDEXPERIMENT_DSMOKE_BIO
#define OPENSMOKE_SOLIDEXPERIMENT_DSMOKE_BIO

#include "BzzMath.hpp"
class OpenSMOKE_ReactingGas;

enum OpenSMOKE_SolidExperiment_DSmoke_Bio_Experimental_Data {DSMOKE_BIO_MASS_RATIO, DSMOKE_BIO_D_MASS_RATIO};
class OpenSMOKE_SolidExperiment_DSmoke_Bio
{
public:

	string name_of_file;

	BzzVector times;
	BzzVector temperatures;
	BzzVector dtemperatures;

	double temperature_mean;
	double pressure_atm;
	vector<string>	solid_names;
	vector<double>	solid_mass_fraction;
	vector<int>		solid_index;

	int nPoints;
	BzzVector x;
	BzzVector y;
	BzzVector initial_masses;

	BzzVectorInt	mass_ratio_indices;
	vector<string>	mass_ratio_names;

	void PrepareMasses(OpenSMOKE_ReactingGas &mix);
	double	GetTemperature(const double t);
	void ReadFromFile(const string filename);
	
	OpenSMOKE_SolidExperiment_DSmoke_Bio_Experimental_Data	kind_of_experimental_data;

public:

	void CheckUserInputs();
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_SOLIDEXPERIMENT_DSMOKE_BIO