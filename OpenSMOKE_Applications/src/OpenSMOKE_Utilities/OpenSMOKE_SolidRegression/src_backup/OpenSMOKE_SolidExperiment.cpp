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

#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "OpenSMOKE_SolidExperiment.h"

void OpenSMOKE_SolidExperiment::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SolidExperiment"		<< endl;
    cout << "File:   " << name_of_file			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SolidExperiment::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SolidExperiment"	<< endl;
    cout << "File:    " << name_of_file			<< endl;
    cout << "Warning: " << message				<< endl;
	cout << endl;
}

void OpenSMOKE_SolidExperiment::ReadFromFile(const string filename)
{
	string dummy;

	name_of_file = filename;
	ifstream fInput;
	openInputFileAndControl(fInput, filename);

	gas_names.push_back("gas names");
	gas_mole_fraction.push_back(0.);
	gas_index.push_back(-1);

	fInput >> dummy;
	if (dummy!="TEMPERATURE")	ErrorMessage("Expected: TEMPERATURE - Found: " + dummy);
	fInput >> temperature;
	fInput >> dummy;
	if (dummy!="K" && dummy!="C")	ErrorMessage("Expected: K || C - Found: " + dummy);
	if (dummy=="C")	temperature += 273.15;

	fInput >> dummy;
	if (dummy!="PRESSURE")	ErrorMessage("Expected: PRESSURE - Found: " + dummy);
	fInput >> pressure;
	fInput >> dummy;
	if (dummy!="atm" && dummy!="bar" && dummy!="Pa")	ErrorMessage("Expected: atm || bar || Pa - Found: " + dummy);
	if (dummy=="bar")	pressure /= 1.01325;
	if (dummy=="Pa")	pressure /= 101325.;

	for(;;)
	{
		fInput >> dummy;
		if (dummy == "SG0") break;
		if (dummy != "GAS_MOLE_FRACTION")	ErrorMessage("Expected: GAS_MOLE_FRACTION - Found: " + dummy);
		fInput >> dummy;
		gas_names.push_back(dummy);
		fInput >> dummy;
		gas_mole_fraction.push_back(atof(dummy.c_str()));
	}

	if (dummy!="SG0")	ErrorMessage("Expected: SG0 - Found: " + dummy);
	fInput >> Sg0;
	fInput >> dummy;
	if (dummy!="m2/kg" && dummy!="m2/g")	ErrorMessage("Expected: m2/kg || m2/g - Found: " + dummy);
	if (dummy=="m2/g")	Sg0 *= 1.e3;

	fInput >> dummy;
	if (dummy!="SIGMA")	ErrorMessage("Expected: SIGMA - Found: " + dummy);
	fInput >> Sigma;
	fInput >> dummy;
	if (dummy!="1/m2")	ErrorMessage("Expected: 1/m2 - Found: " + dummy);

	for(;;)
	{
		fInput >> dummy;
		if (dummy == "//")	break;
		x.Append(atof(dummy.c_str()));

		fInput >> dummy;
		y.Append(atof(dummy.c_str()));
	}
	fInput.close();

	nPoints = x.Size();

	cout << "-------------------------------------------"		<< endl;
	cout << "Experiment:      " << filename					<< endl;
	cout << "-------------------------------------------"		<< endl;
	cout << "  * Temperature:  " << temperature	<< " K"		<< endl;
	cout << "  * Pressure:     " << pressure	<< " atm"		<< endl;
	for(int j=1;j<=int(gas_names.size())-1;j++)
		cout << "  * x_mole " << gas_names[j] << "\t" << gas_mole_fraction[j] << endl;
	cout << "  * Sg0:          " << Sg0	<< " m2/kg"			<< endl;
	cout << "  * Sigma:        " << Sigma	<< " 1/m2"			<< endl;
	cout << "  * #points:      " << nPoints					<< endl;
	cout << "  * x0:           " << x[1]						<< endl;
	cout << "  * xF:           " << x[nPoints]				<< endl;
	cout << "  * yMean:        " << Mean(y)					<< endl;
	cout << "-------------------------------------------"		<< endl;
	cout << endl;
}


void OpenSMOKE_SolidExperiment::PrepareConcentrations(const vector<string> list_of_gas_names)
{
	ChangeDimensions(list_of_gas_names.size()-1, &gas_c);

	for(int j=1;j<=int(gas_names.size())-1;j++)
	{
		bool iFound = false;
		for(int i=1;i<=int(list_of_gas_names.size())-1;i++)
		{
			if (gas_names[j] == list_of_gas_names[i])
			{
				gas_c[i] = gas_mole_fraction[j]*(101325.*pressure)/Constants::R_J_kmol/temperature;	// [kmol/m3]
				iFound = true;
				break;
			}
		}

		if (iFound == false)
			ErrorMessage("This species was not declared in the kinetic scheme: " + gas_names[j]);
	}
}