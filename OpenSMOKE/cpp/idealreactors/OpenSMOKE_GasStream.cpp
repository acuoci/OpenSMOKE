/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
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

#include <iomanip>
#include "basic/OpenSMOKE_Conversions.h"
#include "basic/OpenSMOKE_Constants.h"
#include "idealreactors/OpenSMOKE_GasStream.h"

const double    OpenSMOKE_GasStream::EPSILON			=  1.e-6;
const double    OpenSMOKE_GasStream::_1_PLUS_EPSILON	=  1.+EPSILON;
const double    OpenSMOKE_GasStream::_1_MINUS_EPSILON   =  1.-EPSILON;

void OpenSMOKE_GasStream::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_GasStream"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_GasStream::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_GasStream"		<< endl;
    cout << "Object:  " << name_object			<< endl;
    cout << "Warning: " << message				<< endl;
	cout << endl;
}


OpenSMOKE_GasStream::OpenSMOKE_GasStream()
{
	name_object					= "[Name not assigned]";

	assignedKineticScheme		= false;
	assignedPressure			= false;
	assignedDensity				= false;
	assignedMassFractions		= false;
	assignedMoleFractions		= false;
	assignedTemperature			= false;
	assignedMassFlowRate		= false;
	assignedMoleFlowRate		= false;
	assignedVolumetricFlowRate	= false;
	iSetEquivalenceRatio		= false;
	iSetFuelComposition			= false;
	iSetOxidizerComposition		= false;
	iUndefinedFlowRate			= false;
	locked						= false;
}

void OpenSMOKE_GasStream::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_GasStream::AssignKineticScheme(OpenSMOKE_ReactingGas &_mix)
{
	mix = &_mix;

	int NC = mix->NumberOfSpecies();
	
	ChangeDimensions(NC, &omega);
	ChangeDimensions(NC, &x);
	ChangeDimensions(NC, &c);			// [kmol/m3]
	ChangeDimensions(NC, &h_mass);		// [J/kg]
	ChangeDimensions(NC, &h_mole);		// [J/kmol]

	assignedKineticScheme = true;
}

void OpenSMOKE_GasStream::AssignMassFlowRate(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The mass flow rate must be greater than zero!");

	massFlowRate = OpenSMOKE_Conversions::conversion_massFlowRate(_value, _units);

	assignedVolumetricFlowRate	= false;
	assignedMassFlowRate		= true;
	assignedMoleFlowRate		= false;
}

void OpenSMOKE_GasStream::AssignMoleFlowRate(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The mole flow rate must be greater than zero!");

	moleFlowRate = OpenSMOKE_Conversions::conversion_moleFlowRate(_value, _units);
	
	assignedVolumetricFlowRate	= false;
	assignedMassFlowRate		= false;
	assignedMoleFlowRate		= true;
}

void OpenSMOKE_GasStream::AssignVolumetricFlowRate(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The volumetric flow rate must be greater than zero!");

	volumetricFlowRate = OpenSMOKE_Conversions::conversion_volumetricFlowRate(_value, _units);
	
	assignedVolumetricFlowRate	= true;
	assignedMassFlowRate		= false;
	assignedMoleFlowRate		= false;
}

void OpenSMOKE_GasStream::AssignTemperature(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The temperature must be greater than zero!");

	T = OpenSMOKE_Conversions::conversion_temperature(_value, _units);
	assignedTemperature = true;
}

void OpenSMOKE_GasStream::AssignDensity(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The density must be greater than zero!");

	rho = OpenSMOKE_Conversions::conversion_density(_value, _units);
	assignedDensity = true;
}

void OpenSMOKE_GasStream::AssignPressure(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The pressure must be greater than zero!");

	P = OpenSMOKE_Conversions::conversion_pressure(_value, _units);
	assignedPressure = true;
}

void OpenSMOKE_GasStream::AssignMassFractions(BzzVector &_values)
{
	double sum = _values.GetSumElements();

	if ( sum < _1_MINUS_EPSILON || sum > _1_PLUS_EPSILON)
		ErrorMessage("The sum of mass fractions must be 1.0!");

	if (_values.Min() < 0.0)
		ErrorMessage("The mass fraction of each species must be greater than zero!");

	if (_values.Size() != mix->NumberOfSpecies())
		ErrorMessage("The mass fractions must be defined for all the species!");

	omega = _values;
	omega /= sum;

	assignedMassFractions = true;
	assignedMoleFractions = false;
}

void OpenSMOKE_GasStream::AssignMoleFractions(BzzVector &_values)
{
	double sum = _values.GetSumElements();

	if ( sum < _1_MINUS_EPSILON || sum > _1_PLUS_EPSILON)
		ErrorMessage("The sum of mole fractions must be 1.0!");

	if (_values.Min() < 0.0)
		ErrorMessage("The mole fraction of each species must be greater than zero!");

	if (_values.Size() != mix->NumberOfSpecies())
		ErrorMessage("The mole fractions must be defined for all the species in the kinetic scheme!");

	x = _values;
	x /= sum;

	assignedMoleFractions = true;
	assignedMassFractions = false;
}

void OpenSMOKE_GasStream::AssignMoleFractions(const vector<string> _names, const vector<double> _values)
{
	BzzVector mole(mix->NumberOfSpecies());

	for(int k=0; k<int(_names.size()); k++)
		mole[mix->recognize_species(_names[k])]	= _values[k];

	AssignMoleFractions(mole);
}

void OpenSMOKE_GasStream::AssignMassFractions(const vector<string> _names, const vector<double> _values)
{
	BzzVector mass(mix->NumberOfSpecies());

	for(int k=0; k<int(_names.size()); k++)
		mass[mix->recognize_species(_names[k])]	= _values[k];

	AssignMassFractions(mass);
}

/*
// Equivalence Ratio in air
void OpenSMOKE_GasStream::AssignEquivalenceRatio(const double equivalence_ratio, const std::string fuel_name)
{
	double nC	= mix->elements[mix->recognize_element("c")][mix->recognize_species(fuel_name)];
	double nH	= mix->elements[mix->recognize_element("h")][mix->recognize_species(fuel_name)];

	double nO2		= (2*nC + 0.50*nH)/2.;
	double nN2		= 0.79/0.21*nO2;
	double nFuel	= equivalence_ratio;
	double n		= nO2 + nN2 + nFuel;


	BzzVector mole(mix->NumberOfSpecies());	
	mole[mix->recognize_species(fuel_name)]	= nFuel/n;
	mole[mix->recognize_species("O2")]		= nO2/n;
	mole[mix->recognize_species("N2")]		= nN2/n;

	AssignMoleFractions(mole);
}*/

void OpenSMOKE_GasStream::AssignEquivalenceRatio(const double equivalence_ratio)
{
	vector<string>	fuel_names;
	BzzVector		fuel_moles;
	vector<string>  oxidizer_names;
	BzzVector		oxidizer_moles;

	if (iSetFuelComposition == false)
		ErrorMessage("EquivalenceRatio option requires fuel composition...");
	else
	{
		fuel_names.push_back("fuel names");
		for(int j=1;j<=mix->NumberOfSpecies();j++)
			if (x_Fuel[j] > 0.)
			{
				fuel_names.push_back(mix->names[j]);
				fuel_moles.Append(x_Fuel[j]);
			}
	}

	if (iSetOxidizerComposition == false)
	{
		oxidizer_names.push_back("oxidizer names");
		oxidizer_names.push_back("O2");
		oxidizer_names.push_back("N2");
		oxidizer_moles.Append(0.21);
		oxidizer_moles.Append(0.79);
	}
	else
	{
		oxidizer_names.push_back("oxidizer names");
		for(int j=1;j<=mix->NumberOfSpecies();j++)
			if (x_Oxidizer[j] > 0.)
			{
				oxidizer_names.push_back(mix->names[j]);
				oxidizer_moles.Append(x_Oxidizer[j]);
			}
	}

	BzzVector moles = mix->GetMoleFractionsFromEquivalenceRatio(equivalence_ratio, fuel_names, fuel_moles, oxidizer_names, oxidizer_moles);
	
	vector<string> _names;
	vector<double> _moles;
	for(int j=1;j<=mix->NumberOfSpecies();j++)
		if (moles[j] > 0.)
		{
			_names.push_back(mix->names[j]);
			_moles.push_back(moles[j]);
		}

	AssignMoleFractions(_names, _moles);

	iSetEquivalenceRatio = true;
}

void OpenSMOKE_GasStream::AssignFuelMoleFractions(const vector<string> _names, const vector<double> _values)
{
	ChangeDimensions(mix->NumberOfSpecies(), &x_Fuel);
	for(int i=1;i<=int(_names.size());i++)
		x_Fuel[mix->recognize_species(_names[i-1])] = _values[i-1];
	iSetFuelComposition = true;
}

void OpenSMOKE_GasStream::AssignOxidizerMoleFractions(const vector<string> _names, const vector<double> _values)
{
	ChangeDimensions(mix->NumberOfSpecies(), &x_Oxidizer);
	for(int i=1;i<=int(_names.size());i++)
		x_Oxidizer[mix->recognize_species(_names[i-1])] = _values[i-1];
	iSetOxidizerComposition = true;
}

void OpenSMOKE_GasStream::AssignFuelMoleFractions(BzzVector &_x_Fuel)
{
	ChangeDimensions(mix->NumberOfSpecies(), &x_Fuel);
	x_Fuel = _x_Fuel;
	iSetFuelComposition = true;
}

void OpenSMOKE_GasStream::AssignOxidizerMoleFractions(BzzVector &_x_Oxidizer)
{
	ChangeDimensions(mix->NumberOfSpecies(), &x_Oxidizer);
	x_Oxidizer = _x_Oxidizer;
	iSetOxidizerComposition = true;
}

void OpenSMOKE_GasStream::AssignFuelMassFractions(const vector<string> _names, const vector<double> _values)
{
	double MW_Fuel;
	BzzVector omega_Fuel(mix->NumberOfSpecies());
	ChangeDimensions(mix->NumberOfSpecies(), &x_Fuel);
	for(int i=1;i<=int(_names.size());i++)
		omega_Fuel[mix->recognize_species(_names[i-1])] = _values[i-1];
	mix->GetMWAndMoleFractionsFromMassFractions(MW_Fuel, x_Fuel, omega_Fuel);
	iSetFuelComposition = true;
}

void OpenSMOKE_GasStream::AssignOxidizerMassFractions(const vector<string> _names, const vector<double> _values)
{
	double MW_Oxidizer;
	BzzVector omega_Oxidizer(mix->NumberOfSpecies());
	ChangeDimensions(mix->NumberOfSpecies(), &x_Oxidizer);
	for(int i=1;i<=int(_names.size());i++)
		omega_Oxidizer[mix->recognize_species(_names[i-1])] = _values[i-1];
	mix->GetMWAndMoleFractionsFromMassFractions(MW_Oxidizer, x_Oxidizer, omega_Oxidizer);
	iSetOxidizerComposition = true;
}

void OpenSMOKE_GasStream::AssignFuelMassFractions(BzzVector &omega_Fuel)
{
	double MW_Fuel;
	ChangeDimensions(mix->NumberOfSpecies(), &x_Fuel);
	mix->GetMWAndMoleFractionsFromMassFractions(MW_Fuel, x_Fuel, omega_Fuel);
	iSetFuelComposition = true;
}

void OpenSMOKE_GasStream::AssignOxidizerMassFractions(BzzVector &omega_Oxidizer)
{
	double MW_Oxidizer;
	ChangeDimensions(mix->NumberOfSpecies(), &x_Oxidizer);
	mix->GetMWAndMoleFractionsFromMassFractions(MW_Oxidizer, x_Oxidizer, omega_Oxidizer);
	iSetOxidizerComposition = true;
}

// Equivalence Ratio (most general)
void OpenSMOKE_GasStream::AssignEquivalenceRatio(const double equivalence_ratio,	const vector<string> fuel_names, BzzVector &moles_fuel,
																					const vector<string> oxidizer_names, BzzVector &moles_oxidizer)
{
	BzzVector mole = EquivalenceRatio(equivalence_ratio,	fuel_names, moles_fuel, 
															oxidizer_names, moles_oxidizer);

	AssignMoleFractions(mole);
}

void OpenSMOKE_GasStream::ChangeEquivalenceRatio(const double equivalence_ratio,	const vector<string> fuel_names, BzzVector &moles_fuel,
																					const vector<string> oxidizer_names, BzzVector &moles_oxidizer)
{
	BzzVector mole = EquivalenceRatio(equivalence_ratio,	fuel_names, moles_fuel, 
															oxidizer_names, moles_oxidizer);
	ChangeMoleFractions(mole);
}

// Equivalence Ratio (most general)
BzzVector OpenSMOKE_GasStream::EquivalenceRatio(const double equivalence_ratio,	const vector<string> fuel_names, BzzVector &moles_fuel,
																			const vector<string> oxidizer_names, BzzVector &moles_oxidizer)
{
	int j;
	int number_of_fuels		= moles_fuel.Size();
	int number_of_oxidizers	= moles_oxidizer.Size();

	BzzVector nC(number_of_fuels);
	BzzVector nH(number_of_fuels);
	BzzVector nO(number_of_fuels);
	BzzVector nOxidizers(number_of_oxidizers);

	int jC = mix->recognize_element_without_error("c");
	int jH = mix->recognize_element_without_error("h");
	int jO = mix->recognize_element_without_error("o");
	
	if (jC>0)
		for(j=1;j<=number_of_fuels;j++)
			nC[j]	= mix->elements[jC][mix->recognize_species(fuel_names[j])];

	if (jH>0)
		for(j=1;j<=number_of_fuels;j++)
			nH[j]	= mix->elements[jH][mix->recognize_species(fuel_names[j])];

	if (jO>0)
		for(j=1;j<=number_of_fuels;j++)
			nO[j]	= mix->elements[jO][mix->recognize_species(fuel_names[j])];

	double nFuel	= equivalence_ratio*moles_fuel.GetSumElements();

	int jO2 = 0;
	for(j=1;j<=number_of_oxidizers;j++)
		if (oxidizer_names[j] == "O2")	{ jO2 = j; break;}

	double nO2		= (	2.*Dot(nC,moles_fuel) + 
						0.50*Dot(nH,moles_fuel) -
						Dot(nO,moles_fuel))/2.;

	for(j=1;j<=number_of_oxidizers;j++)
		nOxidizers[j] = moles_oxidizer[j]/moles_oxidizer[jO2]*nO2;
	
	double n		= nFuel + nOxidizers.GetSumElements();

	BzzVector mole(mix->NumberOfSpecies());	
	for(j=1;j<=number_of_fuels;j++)
		mole[mix->recognize_species(fuel_names[j])] += equivalence_ratio*moles_fuel[j]/n;
	for(j=1;j<=number_of_oxidizers;j++)
		mole[mix->recognize_species(oxidizer_names[j])] += nOxidizers[j]/n;

	//for(j=1;j<=number_of_fuels;j++)
	//	cout << fuel_names[j] << " " << mole[mix->recognize_species(fuel_names[j])] << endl;
	//for(j=1;j<=number_of_oxidizers;j++)
	//	cout << oxidizer_names[j] << " " << mole[mix->recognize_species(oxidizer_names[j])] << endl;

	return mole;
}

void OpenSMOKE_GasStream::GiveMeStoichiometricProducts(BzzVector &x_reactants, BzzVector &fuel,
													   BzzVector &mole_reactants, BzzVector &mole_products)
{
	double nC =0;
	double nH =0;
	double nO =0;
	double nN =0;
	double nHE=0;
	double nAR=0;

	int jC = mix->recognize_element_without_error("c");
	int jH = mix->recognize_element_without_error("h");
	int jO = mix->recognize_element_without_error("o");
	int jN = mix->recognize_element_without_error("n");
	int jAR = mix->recognize_element_without_error("ar");
	int jHE = mix->recognize_element_without_error("he");
	
	if (jC>0)
		for(int j=1;j<=mix->NumberOfSpecies();j++)
			nC += mix->elements[jC][j]*x_reactants[j];
	if (jH>0)
		for(int j=1;j<=mix->NumberOfSpecies();j++)
			nH += mix->elements[jH][j]*x_reactants[j];
	if (jO>0)
		for(int j=1;j<=mix->NumberOfSpecies();j++)
			nO += mix->elements[jO][j]*x_reactants[j];
	if (jN>0)
		for(int j=1;j<=mix->NumberOfSpecies();j++)
			nN += mix->elements[jN][j]*x_reactants[j];
	if (jHE>0)
		for(int j=1;j<=mix->NumberOfSpecies();j++)
			nHE += mix->elements[jHE][j]*x_reactants[j];
	if (jAR>0)
		for(int j=1;j<=mix->NumberOfSpecies();j++)
			nAR += mix->elements[jAR][j]*x_reactants[j];

	BzzVector x_products(mix->NumberOfSpecies());

	int jCO2 = mix->recognize_species_without_exit("CO2");
	if (jCO2>0)	x_products[jCO2]=nC;
	int jH2O = mix->recognize_species_without_exit("H2O");
	if (jH2O>0)	x_products[jH2O]=nH/2;
	int jN2 = mix->recognize_species_without_exit("N2");
	if (jN2>0)	x_products[jN2]=nN/2;
	int jjHE = mix->recognize_species_without_exit("HE");
	if (jjHE>0)	x_products[jjHE]=nHE/2;
	int jjAR = mix->recognize_species_without_exit("AR");
	if (jjAR>0)	x_products[jjAR]=nAR/2;

	double sum =0.;
	for(int j=1;j<=mix->NumberOfSpecies();j++)
		if (fuel[j]>0.)	sum+=x_reactants[j];

	mole_reactants = x_reactants;
	mole_reactants /= sum;
	mole_products = x_products;
	mole_products /= sum;

	cout << "Reactants:" << endl;
	for(int i=1;i<=mix->NumberOfSpecies();i++)
			if (mole_reactants[i] > 0.)	
			{
				cout << setw(20) << left << mix->names[i];
				cout << setw(10) << left << mole_reactants[i] << endl;
			}

	cout << "Products:" << endl;
	for(int i=1;i<=mix->NumberOfSpecies();i++)
			if (mole_products[i] > 0.)	
			{
				cout << setw(20) << left << mix->names[i];
				cout << setw(10) << left << mole_products[i] << endl;
			}
}

double OpenSMOKE_GasStream::GiveMeReactionEnthalpy(BzzVector &mole_reactants, BzzVector &mole_products)
{
	BzzVector H(mix->NumberOfSpecies());
	mix->GetStandardEnthalpy_Mole(H, T);			// [J/kmol]

	double Hreactants = 0.;
	for(int j=1;j<=mix->NumberOfSpecies();j++)
		Hreactants += H[j]*mole_reactants[j];
	double Hproducts = 0.;
	for(int j=1;j<=mix->NumberOfSpecies();j++)
		Hproducts += H[j]*mole_products[j];

	return (Hproducts - Hreactants);
}


void OpenSMOKE_GasStream::lock()
{
	if (assignedKineticScheme == false)
		ErrorMessage("The kinetic scheme was not defined!!");
	if (assignedMassFlowRate == false && assignedMoleFlowRate == false && assignedVolumetricFlowRate == false)
	{
		//	ErrorMessage("The flow rate was not defined!!");
		AssignMassFlowRate(1.e-10, "kg/s");
		iUndefinedFlowRate = true;
	}
	if (assignedMoleFractions == false && assignedMassFractions == false)
		ErrorMessage("The composition was not defined!!");

	Composition();

	if (assignedTemperature == true && assignedPressure == true && assignedDensity == true)
		ErrorMessage("Only 2 between: T || P || rho");
	if (assignedTemperature == true && assignedPressure == true)
		{ /* Nothing to do */						}
	else if (assignedDensity == true && assignedPressure == true)
		{	T = P*MW/Constants::R_J_kmol/rho;	}
	else if (assignedDensity == true && assignedTemperature == true)
		{	P = rho*Constants::R_J_kmol*T/MW;	}
	else ErrorMessage("2 between must be assigned: T || P || rho");
		
	Concentrations();
	Density();
	FlowRates();
	SpecificEnthalpies();
	Enthalpy();
	SpecificEntropies();
	Entropy();
}

void OpenSMOKE_GasStream::SetVelocity(const double _value, const std::string _units)
{
	velocity = OpenSMOKE_Conversions::conversion_velocity(_value, _units);
}

void OpenSMOKE_GasStream::VideoSummary()
{
	cout.setf(ios::scientific);
	cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
	cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
	cout << endl;
	
	cout << "---------------------------------------------------------------------"		<< endl;
	cout << "Gas Stream Summary                                                   "		<< endl;
	cout << "---------------------------------------------------------------------"		<< endl;
	cout << "Pressure:\t\t"				<< P					<< " [Pa]"				<< endl;
	cout << "Temperature:\t\t"			<< T					<< " [K]"				<< endl;
	cout << "Density:\t\t"				<< rho					<< " [kg/m3]"			<< endl;
	cout << "Molecular weight:\t"		<< MW					<< " [kg/kmol]"			<< endl;
	cout << "Total concentration:\t"	<< cTot					<< " [kmol/m3]"			<< endl;
	cout << "Mass flow rate:\t\t"		<< massFlowRate			<< " [kg/s]"			<< endl;
	cout << "Mole flow rate:\t\t"		<< moleFlowRate			<< " [kmol/s]"			<< endl;
	cout << "Volumetric flow rate:\t"	<< volumetricFlowRate	<< " [m3/s]"			<< endl;
	cout << "Specific enthalpy:\t"		<< massSpecificEnthalpy	<< " [J/kg]"			<< endl;
	cout << "Specific enthalpy:\t"		<< moleSpecificEnthalpy	<< " [J/kmol]"			<< endl;
	cout << "Enthalpy:\t\t"				<< enthalpy				<< " [J/s]"				<< endl;
	cout << endl;

	cout << "#\tName\t\tx\t\tomega\t\tc[kmol/m3]" << endl;
	cout << "---------------------------------------------------------------------" << endl;
	for(int i=1;i<=mix->NumberOfSpecies();i++)
		if (x[i]!=0.) 
			cout	<< i				<< "\t" 
					<< mix->names[i]	<< "\t\t" 
					<< x[i]				<< "\t" 
					<< omega[i]			<< "\t" 
					<< c[i]				<< "\t" 
					<< endl;
	cout << endl;	
	cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
	cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
	cout << endl;

}

void OpenSMOKE_GasStream::ChangeTemperature(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The temperature must be greater than zero!");

	T = OpenSMOKE_Conversions::conversion_temperature(_value, _units);

	Concentrations();
	Density();
	FlowRates();
	SpecificEnthalpies();
	Enthalpy();
	SpecificEntropies();
	Entropy();
}

void OpenSMOKE_GasStream::ChangePressure(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The pressure must be greater than zero!");

	P = OpenSMOKE_Conversions::conversion_pressure(_value, _units);

	Concentrations();
	Density();
	FlowRates();
	Enthalpy();
}

void OpenSMOKE_GasStream::ChangeMassFlowRate(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The mass flow rate must be greater than zero!");

	massFlowRate = OpenSMOKE_Conversions::conversion_massFlowRate(_value, _units);
	
	assignedVolumetricFlowRate	= false;
	assignedMassFlowRate		= true;
	assignedMoleFlowRate		= false;

	FlowRates();
	Enthalpy();
}

void OpenSMOKE_GasStream::ChangeMoleFlowRate(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The mole flow rate must be greater than zero!");

	moleFlowRate = OpenSMOKE_Conversions::conversion_moleFlowRate(_value, _units);
	
	assignedVolumetricFlowRate	= false;
	assignedMassFlowRate		= false;
	assignedMoleFlowRate		= true;

	FlowRates();
	Enthalpy();
}

void OpenSMOKE_GasStream::ChangeVolumetricFlowRate(const double _value, const std::string _units)
{
	if (_value<=0.)
		ErrorMessage("The volumetric flow rate must be greater than zero!");

	volumetricFlowRate = OpenSMOKE_Conversions::conversion_volumetricFlowRate(_value, _units);
	
	assignedVolumetricFlowRate	= true;
	assignedMassFlowRate		= false;
	assignedMoleFlowRate		= false;

	FlowRates();
	Enthalpy();
}

void OpenSMOKE_GasStream::ChangeMassFractions(BzzVector &_values)
{
	AssignMassFractions(_values);

	Composition();
	Concentrations();
	Density();
	FlowRates();
	SpecificEnthalpies();
	Enthalpy();
	SpecificEntropies();
	Entropy();
}

void OpenSMOKE_GasStream::ChangeMoleFractions(BzzVector &_values)
{
	AssignMoleFractions(_values);

	Composition();
	Concentrations();
	Density();
	FlowRates();
	SpecificEnthalpies();
	Enthalpy();
	SpecificEntropies();
	Entropy();
}

void OpenSMOKE_GasStream::Composition()
{	
	if(assignedMoleFractions == true)
		mix->GetMWAndMassFractionsFromMoleFractions(MW, omega, x);
	
	else if(assignedMassFractions == true)
		mix->GetMWAndMoleFractionsFromMassFractions(MW, x, omega);
}

void OpenSMOKE_GasStream::Concentrations()
{
	cTot	= P/(Constants::R_J_kmol*T);		// [kmol/m3]
	c		= cTot*x;							// [kmol/m3]
}

void OpenSMOKE_GasStream::Density()
{
	rho		= cTot*MW;							// [kg/m3]
}

void OpenSMOKE_GasStream::FlowRates()
{	
	if (assignedMassFlowRate == true)
	{
		moleFlowRate		= massFlowRate/MW;		// [kmol/s]
		volumetricFlowRate	= massFlowRate/rho;		// [m3/s]
	}
	else if (assignedMoleFlowRate == true)
	{
		massFlowRate		= moleFlowRate*MW;		// [kg/s]
		volumetricFlowRate	= massFlowRate/rho;		// [m3/s]
	}
	else if (assignedVolumetricFlowRate == true)
	{
		massFlowRate	= volumetricFlowRate*rho;	// [kg/s]
		moleFlowRate	= massFlowRate/MW;			// [kmol/s]
	}
}
	
void OpenSMOKE_GasStream::SpecificEnthalpies()
{	
	mix->GetMixAveragedEnthalpy_Mole(h_mole, T);			// [J/kmol]
	ElementByElementProduct(h_mole, mix->uM, &h_mass);		// [J/kg]
	
	massSpecificEnthalpy = Dot(omega, h_mass);				// [J/kg]
	moleSpecificEnthalpy = Dot(x, h_mole);					// [J/kmol]

	massSpecificPressureEnergy = P/rho;
	massSpecificInternalEnergy = massSpecificEnthalpy - massSpecificPressureEnergy;	// [J/kg]
}

void OpenSMOKE_GasStream::SpecificEntropies()
{	
	massSpecificEntropy = mix->GetMixEntropy_Mass(P, T, omega);	// [J/kg/K]
}

void OpenSMOKE_GasStream::Enthalpy()
{	
	enthalpy = massSpecificEnthalpy * massFlowRate;			// [J/s]
}

void OpenSMOKE_GasStream::Entropy()
{	
	entropy = massSpecificEntropy * massFlowRate;			// [J/s/K]
}
void OpenSMOKE_GasStream::DefineFromFile(const std::string inputFile)
{
    double			double_value;
    std::string			string_value;
    vector<double>  double_vector;
    vector<string>  string_vector;
    

    OpenSMOKE_Dictionary_GasStream dictionary;
    dictionary.ParseFile(inputFile);

    // COMPULSORY: Temperature
    if (dictionary.Return("#Temperature", double_value, string_value))
		AssignTemperature(double_value, string_value);

    // COMPULSORY: Pressure
    if (dictionary.Return("#Pressure", double_value, string_value))
		AssignPressure(double_value, string_value);

    // COMPULSORY: Pressure
    if (dictionary.Return("#Density", double_value, string_value))
		AssignDensity(double_value, string_value);

    // SEMI-COMPULSORY: Mass Flow Rate
    if (dictionary.Return("#MassFlowRate", double_value, string_value))
		AssignMassFlowRate(double_value, string_value);

    // SEMI-COMPULSORY: Mole Flow Rate
    else if (dictionary.Return("#MoleFlowRate", double_value, string_value))
		AssignMoleFlowRate(double_value, string_value);

    // SEMI-COMPULSORY: Volumetric Flow Rate
    else if (dictionary.Return("#VolumetricFlowRate", double_value, string_value))
		AssignVolumetricFlowRate(double_value, string_value);

    // SEMI-COMPULSORY: Mass fractions
    if (dictionary.Return("#MassFractions", double_vector, string_vector))
		AssignMassFractions(string_vector, double_vector);

    // SEMI-COMPULSORY: Mole fractions
    else if (dictionary.Return("#MoleFractions", double_vector, string_vector))
		AssignMoleFractions(string_vector, double_vector);

    // SEMI-COMPULSORY: EquivalenceRatio
	if (dictionary.Return("#FuelMoleFractions", double_vector, string_vector))
		AssignFuelMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#FuelMassFractions", double_vector, string_vector))
		AssignFuelMassFractions(string_vector, double_vector);

	if (dictionary.Return("#OxidizerMoleFractions", double_vector, string_vector))
		AssignOxidizerMoleFractions(string_vector, double_vector);

	if (dictionary.Return("#OxidizerMassFractions", double_vector, string_vector))
		AssignOxidizerMassFractions(string_vector, double_vector);

	if (dictionary.Return("#EquivalenceRatio", double_value))
		AssignEquivalenceRatio(double_value);

	lock();
}

OpenSMOKE_Dictionary_GasStream::OpenSMOKE_Dictionary_GasStream()
{
    SetupBase();
	SetName("OpenSMOKE_GasStream Dictionary");

    Add("#Temperature",			'O', 'M', "Gas stream temperature");
    Add("#Pressure",			'O', 'M', "Gas stream pressure");
    Add("#Density",				'O', 'M', "Gas stream density");
    
	Add("#MassFlowRate",		'O', 'M', "Gas stream mass flow rate");
    Add("#MoleFlowRate",		'O', 'M', "Gas stream mole flow rate");
    Add("#VolumetricFlowRate",  'O', 'M', "Gas stream volumetric flow rate");

    Add("#MassFractions",		'O', 'L', "Gas stream mass fractions");
    Add("#MoleFractions",		'O', 'L', "Gas stream mole fractions");

	Add("#EquivalenceRatio",	'O', 'D', "Equivalence Ratio");
	Add("#FuelMoleFractions",	'O', 'L', "Fuel composition (mole fractions)");
	Add("#OxidizerMoleFractions",	'O', 'L', "Oxidizer composition (mole fractions)");
	Add("#FuelMassFractions",	'O', 'L', "Fuel composition (mass fractions)");
	Add("#OxidizerMassFractions",	'O', 'L', "Oxidizer composition (mass fractions)");

    Conflict("#MassFlowRate",   "#MoleFlowRate");
    Conflict("#MassFlowRate",   "#VolumetricFlowRate");
    Conflict("#MoleFlowRate",   "#VolumetricFlowRate");

	Conflict("#MassFractions",			"#MoleFractions");
	Conflict("#FuelMoleFractions",		"#FuelMassFractions");
	Conflict("#OxidizerMoleFractions",		"#OxidizerMassFractions");
	Conflict("#FuelMoleFractions",		"#MassFractions");
	Conflict("#FuelMoleFractions",		"#MoleFractions");
	Conflict("#EquivalenceRatio",		"#MassFractions");
	Conflict("#EquivalenceRatio",		"#MoleFractions");


	Compulsory("#MassFractions", "#MoleFractions", "#EquivalenceRatio");
	// Compulsory("#MassFlowRate",  "#MoleFlowRate", "#VolumetricFlowRate");

    Lock();
}

void OpenSMOKE_GasStream::Sum(OpenSMOKE_GasStream& gas1, OpenSMOKE_GasStream& gas2)
{
	if (fabs(gas1.P - gas2.P)/gas1.P > 1.e-5)
		ErrorMessage("The mixing gas stream do not have the same pressure!");

	// Presure   
	double PSum = gas1.P;
	
	// MassFlowRate
	double massFlowRateSum = gas1.massFlowRate + gas2.massFlowRate;

	// Mass Fractions
	BzzVector omegaSum(mix->NumberOfSpecies());
	for(int i=1;i<=mix->NumberOfSpecies();i++)
		omegaSum[i] = gas1.massFlowRate*gas1.omega[i] + gas2.massFlowRate*gas2.omega[i];
	Product(1./massFlowRateSum, &omegaSum);

	// Temperature
	gas1.Enthalpy();
	gas2.Enthalpy();
	double enthalpySum = (gas1.enthalpy + gas2.enthalpy) / massFlowRateSum;										// [J/kg]
	double TFirstGuess = (gas1.T*gas1.massFlowRate+gas2.T*gas2.massFlowRate)/massFlowRateSum;
	double TSum = mix->GetTemperatureFromMassEnthalpyAndMassFractions(TFirstGuess, enthalpySum, omegaSum, 1.e-8);	// temperature [K]

	
    // COMPULSORY: Temperature
	AssignTemperature(TSum, "K");

	// COMPULSORY: Pressure
	AssignPressure(PSum, "Pa");

    // SEMI-COMPULSORY: Mass Flow Rate
	AssignMassFlowRate(massFlowRateSum, "kg/s");

    // SEMI-COMPULSORY: Mass fractions
	AssignMassFractions(omegaSum);

	lock();
}