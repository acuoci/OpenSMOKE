/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
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

#ifndef OPENSMOKE_GASSTREAM
#define OPENSMOKE_GASSTREAM

#include <vector>
#include "engine/OpenSMOKE_ReactingGas.h"
#include "basic/OpenSMOKE_Dictionary.h"

class OpenSMOKE_GasStream
{
private:

	// Name
	string name_object;

	// Control variables
	bool assignedKineticScheme;
	bool assignedPressure;
	bool assignedDensity;
	bool assignedMassFractions;
	bool assignedMoleFractions;
	bool assignedTemperature;
	bool assignedMassFlowRate;
	bool assignedMoleFlowRate;
	bool assignedVolumetricFlowRate;
	bool locked;

	// Gas stream properties
	void Composition();
	void Concentrations();
	void Density();
	void FlowRates();
	void SpecificEnthalpies();
	void Enthalpy();
	void SpecificEntropies();
	void Entropy();
	
	BzzVector EquivalenceRatio( const double equivalence_ratio,	const vector<string> fuel_names, BzzVector &moles_fuel,
								const vector<string> oxidizer_names, BzzVector &moles_oxidizer);

	// Error and warning messages
	void ErrorMessage(string message);
	void WarningMessage(string message);

	static const double EPSILON;
	static const double _1_PLUS_EPSILON;
	static const double _1_MINUS_EPSILON;

	bool iSetEquivalenceRatio;
	bool iSetFuelComposition;
	bool iSetOxidizerComposition;

	BzzVector x_Fuel;
	BzzVector x_Oxidizer;

public:

	OpenSMOKE_GasStream();

	BzzVector omega;
	BzzVector x;
	BzzVector c;
	BzzVector h_mass;
	BzzVector h_mole;

	double P;
	double T;
	double rho;
	double MW;
	double cTot;
	double massFlowRate;
	double moleFlowRate;
	double volumetricFlowRate;
	double massSpecificEnthalpy;
	double moleSpecificEnthalpy;
	double massSpecificEntropy;
	double moleSpecificEntropy;
	double enthalpy;
	double entropy;
	double massSpecificPressureEnergy;
	double massSpecificInternalEnergy;
	double velocity;

	OpenSMOKE_ReactingGas *mix;

	// Assignements (compulsory)
	void AssignKineticScheme(OpenSMOKE_ReactingGas &_mix);
	void AssignTemperature(const double _value, const string _units);
	void AssignDensity(const double _value, const string _units);
	void AssignPressure(const double _value, const string _units);
	void AssignMassFlowRate(const double _value, const string _units);
	void AssignMoleFlowRate(const double _value, const string _units);
	void AssignVolumetricFlowRate(const double _value, const string _units);
	void AssignMoleFractions(BzzVector &_values);
	void AssignMassFractions(BzzVector &_values);
	void AssignMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignMassFractions(const vector<string> _names, const vector<double> _values);
	void AssignEquivalenceRatio(const double equivalence_ratio);
	void AssignFuelMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignOxidizerMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignFuelMassFractions(const vector<string> _names, const vector<double> _values);
	void AssignOxidizerMassFractions(const vector<string> _names, const vector<double> _values);
	void AssignFuelMoleFractions(BzzVector &_values);
	void AssignOxidizerMoleFractions(BzzVector &_values);
	void AssignFuelMassFractions(BzzVector &_values);
	void AssignOxidizerMassFractions(BzzVector &_values);
	
	void AssignEquivalenceRatio(const double equivalence_ratio,	const vector<string> fuel_names, BzzVector &moles_fuel,
																const vector<string> oxidizer_names, BzzVector &moles_oxidizer);
	void ChangeEquivalenceRatio(const double equivalence_ratio,	const vector<string> fuel_names, BzzVector &moles_fuel,
																const vector<string> oxidizer_names, BzzVector &moles_oxidizer);
	void GiveMeStoichiometricProducts(BzzVector &x_reactants, BzzVector &fuel,
									  BzzVector &mole_reactants, BzzVector &mole_products);
	double GiveMeReactionEnthalpy(BzzVector &mole_reactants, BzzVector &mole_products);

	void lock();

	// Assignements (non-compulsory)
	void SetVelocity(const double _value, const string _units);
	void SetName(const string name);

	// Define from file
	void DefineFromFile(const string inputFile);
	
	// Change stream parameters
	void ChangeTemperature(const double _value, const string _units);
	void ChangePressure(const double _value, const string _units);
	void ChangeMassFlowRate(const double _value, const string _units);
	void ChangeMoleFlowRate(const double _value, const string _units);
	void ChangeVolumetricFlowRate(const double _value, const string _units);
	void ChangeMassFractions(BzzVector &_values);
	void ChangeMoleFractions(BzzVector &_values);	

	void Sum(OpenSMOKE_GasStream& gas1, OpenSMOKE_GasStream& gas2);

	void VideoSummary();

	bool iUndefinedFlowRate;
};

class OpenSMOKE_Dictionary_GasStream : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_GasStream();
};

#endif // OPENSMOKE_GASSTREAM


