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
#include "basic/OpenSMOKE_Conversions.h"
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
	double dummy_double;

	string x_type;
	string y_type;
	string x_units;
	string y_units;

	name_of_file = filename;
	ifstream fInput;
	openInputFileAndControl(fInput, filename);

	solid_names.push_back("solid names");
	solid_mass_fraction.push_back(0.);
	gas_names.push_back("gas names");
	gas_mole_fraction.push_back(0.);
	gas_index.push_back(-1);

	// TODO
	for(;;)
	{
		fInput >> dummy;
		if (dummy == "PRESSURE") break;

		if (dummy!="TEMPERATURE")	ErrorMessage("Expected: TEMPERATURE - Found: " + dummy);

		fInput >> dummy_double;
		fInput >> dummy;
		times.Append(OpenSMOKE_Conversions::conversion_time(dummy_double, dummy));
		fInput >> dummy_double;
		fInput >> dummy;
		temperatures.Append(OpenSMOKE_Conversions::conversion_temperature(dummy_double, dummy));
	}
	

	if (dummy!="PRESSURE")	ErrorMessage("Expected: PRESSURE - Found: " + dummy);
	fInput >> pressure_atm;
	fInput >> dummy;
	pressure_atm = OpenSMOKE_Conversions::conversion_pressure(pressure_atm, dummy)/101325.;


	for(;;)
	{
		fInput >> dummy;
		if (dummy == "SOLID_MASS_FRACTION") break;
		if (dummy != "GAS_MOLE_FRACTION")	ErrorMessage("Expected: GAS_MOLE_FRACTION - Found: " + dummy);
		fInput >> dummy;
		gas_names.push_back(dummy);
		fInput >> dummy;
		gas_mole_fraction.push_back(atof(dummy.c_str()));
	}

	for(;;)
	{
		if (dummy == "SG0") break;
		if (dummy != "SOLID_MASS_FRACTION")	ErrorMessage("Expected: SOLID_MASS_FRACTION - Found: " + dummy);
		fInput >> dummy;
		solid_names.push_back(dummy);
		fInput >> dummy;
		solid_mass_fraction.push_back(atof(dummy.c_str()));
		fInput >> dummy;
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

	fInput >> dummy;
	if (dummy!="PSI")	ErrorMessage("Expected: PSI - Found: " + dummy);
	fInput >> Psi;

	for(;;)
	{
		fInput >> dummy;
		if (dummy == "X") break;
		if (dummy != "SOLID_MASS_FRACTION")	ErrorMessage("Expected: SOLID_MASS_FRACTION - Found: " + dummy);
		fInput >> dummy;
		solid_names.push_back(dummy);
		fInput >> dummy;
		solid_mass_fraction.push_back(atof(dummy.c_str()));
	}

	if (dummy != "X") ErrorMessage("Expected: X - Found: " + dummy);
	fInput >> x_type;	if (x_type != "TIME")	ErrorMessage("Expected: TIME - Found: " + x_type);
	fInput >> dummy;	if (dummy  != "UNITS")	ErrorMessage("Expected: UNITS - Found: " + dummy);
	fInput >> x_units;

	// Read x and y units
	fInput >> dummy;
	if (dummy != "Y")			ErrorMessage("Expected: Y - Found: " + dummy);
	fInput >> y_type;
	if (y_type == "CONVERSION_CHAR")
	{
		for(;;)
		{
			fInput >> dummy;
			break;
			// ... //
		}
	}
	else	ErrorMessage("Expected: CONVERSION_CHAR - Found: " + y_type);
	
	if (dummy != "UNITS")		ErrorMessage("Expected: UNITS - Found: " + dummy);
	fInput >> y_units;


	for(;;)
	{
		fInput >> dummy;
		if (dummy == "//")	break;
		if (x_type == "TIME")	x.Append( OpenSMOKE_Conversions::conversion_time(atof(dummy.c_str()), x_units) );

		fInput >> dummy;
		if (y_type == "CONVERSION_CHAR")	y.Append(atof(dummy.c_str()));
	}
	fInput.close();

	nPoints = x.Size();

	CheckUserInputs();

	cout << "-------------------------------------------"		<< endl;
	cout << "Experiment:      " << filename					<< endl;
	cout << "-------------------------------------------"		<< endl;
	cout << "  * Temperature (initial): " << temperatures[1]					<< " K"			<< endl;
	cout << "  * Temperature (final):   " << temperatures[temperatures.Size()]	<< " K"			<< endl;
	cout << "  * Pressure:              " << pressure_atm						<< " atm"		<< endl;
	for(int j=1;j<=int(gas_names.size())-1;j++)
		cout << "  * x_mole " << gas_names[j] << "\t" << gas_mole_fraction[j] << endl;
	cout << "  * Sg0:          " << Sg0		<< " m2/kg"			<< endl;
	cout << "  * Sigma:        " << Sigma	<< " 1/m2"			<< endl;
	cout << "  * Psi:          " << Psi  						<< endl;
	cout << "  * #points:      " << nPoints						<< endl;
	cout << "  * x0:           " << x[1]		<< " s"		<< endl;
	cout << "  * xF:           " << x[nPoints]	<< " s"		<< endl;
	cout << "  * yMean:        " << Mean(y)						<< endl;
	cout << "-------------------------------------------"		<< endl;
	cout << endl;
}

void OpenSMOKE_SolidExperiment::CheckUserInputs()
{
	int j;

	if (times.Size()<2)	ErrorMessage("At least a couple of temperatures must be specified...");
	
	for(j=2;j<=times.Size();j++)
	{
		if (times[j]<=times[j-1])								ErrorMessage("Error in the temperature profile definition...");
		if (times[j]<0.)										ErrorMessage("Error in the temperature profile definition...");
		if (temperatures[j]<=200. || temperatures[j]>=5000.)	ErrorMessage("Error in the temperature profile definition...");
	}
	
	if (pressure_atm <= 1e-10 || pressure_atm >= 1000.)			ErrorMessage("Error in the pressure definition...");

	if (x[x.Size()] >= times[times.Size()])						ErrorMessage("The last time in the temperature profile is not enough to cover the whole experiment...");

	for(j=2;j<=x.Size();j++)
		if (x[j] <= x[j-1])										ErrorMessage("Error in the x profile definition...");
}

void OpenSMOKE_SolidExperiment::PrepareMasses(OpenSMOKE_CharKineticScheme &kinetics)
{
	int j;

	double sum = 0.;
	for(j=1;j<=int(solid_names.size())-1;j++)
		sum += solid_mass_fraction[j];

	if ( sum > 1.0001 || sum < 0.9999)
		ErrorMessage("The sum of solid mass fractions must be equal to 1...");
	
	for(j=1;j<=int(solid_names.size())-1;j++)
		solid_mass_fraction[j] /= sum;

	ChangeDimensions(kinetics.nCarbon, &initial_carbon_mass_fractions);
	for(j=1;j<=int(solid_names.size())-1;j++)
		initial_carbon_mass_fractions[kinetics.RecognizeCarbonSpecies(solid_names[j])] = solid_mass_fraction[j];

//	for(j=1;j<=int(mass_ratio_names.size())-1;j++)
//		mass_ratio_indices.Append(mix.recognize_species(mass_ratio_names[j]));

	temperature_mean = temperatures.GetSumElements()/temperatures.Size();

	ChangeDimensions(times.Size()-1, &dtemperatures);
	for(j=1;j<=times.Size()-1;j++)
		dtemperatures[j] = (temperatures[j+1]-temperatures[j])/(times[j+1]-times[j]);
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
				index_gas_species.Append(i);
				iFound = true;
				break;
			}
		}

		if (iFound == false)
			ErrorMessage("This species was not declared in the kinetic scheme: " + gas_names[j]);
	}
}

void OpenSMOKE_SolidExperiment::UpdateGasConcentrations(const double temperature)
{
	for(int j=1;j<=index_gas_species.Size();j++)
	{
		int i = index_gas_species[j];
		gas_c[i] = gas_mole_fraction[j]*(101325.*pressure_atm)/Constants::R_J_kmol/temperature;	// [kmol/m3]
	}
}

double OpenSMOKE_SolidExperiment::GetTemperature(const double t)
{
	for(int j=1;j<=times.Size()-1;j++)
		if (t<=times[j+1])
			return temperatures[j] + dtemperatures[j]*(t-times[j]);
	if (t>times[times.Size()])	return temperatures[times.Size()];
	ErrorMessage("Time over the minimum value in the defined temperature profile...");
	return -1;
}