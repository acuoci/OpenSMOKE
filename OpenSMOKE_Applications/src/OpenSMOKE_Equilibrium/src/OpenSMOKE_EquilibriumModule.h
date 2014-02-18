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

#ifndef OPENSMOKE_EQUILIBRIUMMODULE_H
#define OPENSMOKE_EQUILIBRIUMmodule_h

#include "BzzMath.hpp"
#include "OpenSMOKE.hpp"

class OpenSMOKE_IdealGas;
class OpenSMOKE_EquilibriumStanjan;

class OpenSMOKE_Dictionary_EquilibriumModule : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_EquilibriumModule();
};

class OpenSMOKE_EquilibriumModule
{
public:

	// Default constructor
	OpenSMOKE_EquilibriumModule();

	void SetName(const string name);
	void SetVerbose();
	void UnsetVerbose();

	void SetInputFileName(const string name);
	void Assign(OpenSMOKE_IdealGas *ideal_gas);
	void Reset();

	void DefineFromFile();

	void AssignConstant(const vector<string> string_vector);
	void AssignElementalMassFractions(const vector<string> string_vector, const vector<double> double_vector);
    void AssignElementalMoleFractions(const vector<string> string_vector, const vector<double> double_vector);
	void AssignMassFractions(const vector<string> string_vector, const vector<double> double_vector);
    void AssignMoleFractions(const vector<string> string_vector, const vector<double> double_vector);
	void AssignEquivalenceRatio(const double _value);
	void AssignFuelMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignFuelMassFractions(const vector<string> _names, const vector<double> _values);
	void AssignOxidizerMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignOxidizerMassFractions(const vector<string> _names, const vector<double> _values);

	void SetTemperature(const string string_value, const double double_value);
	void SetPressure(const string string_value, const double double_value);
	void SetEnthalpy(const string string_value, const double double_value);
	void SetEntropy(const string string_value, const double double_value);
	void SetInternalEnergy(const string string_value, const double double_value);
	void SetDensity(const string string_value, const double double_value);
	void SetVolume(const string string_value, const double double_value);

	void Run();
	void Solve();

private:
	
	string name_object;
	string name_input_file;
	OpenSMOKE_IdealGas *gas;

	void CheckDictionary(OpenSMOKE_Dictionary_EquilibriumModule &dictionary);
	void GetInitialConditions();
	void Lock();
	void PrintSummary();

	equilibriumProblem iEquilibriumProblem;

	bool iVerbose;
	bool iSetTemperature;
	bool iSetPressure;
	bool iSetEnthalpy;
	bool iSetEntropy;
	bool iSetInternalEnergy;
	bool iSetDensity;
	bool iSetVolume;
	bool iSetSpeciesComposition;
	bool iSetEquivalenceRatio;
	bool iSetFuelComposition;
	bool iSetOxidizerComposition;

	double T_0;
	double P_0_Pascal;
	double H_0;
	double U_0;
	double S_0;
	double V_0;
	double rho_0;
	double A_0;
	double G_0;
	double Cp_0;
	double N_0;

	double T_E;
	double P_E_Pascal;
	double H_E;
	double U_E;
	double S_E;
	double V_E;
	double rho_E;
	double A_E;
	double G_E;
	double Cp_E;
	double N_E;

	BzzVector x_elements;
	BzzVector omega_elements;

	BzzVector x_0;
	BzzVector omega_0;
	double MWmix_0;

	BzzVector x_E;
	BzzVector omega_E;
	double MWmix_E;

	BzzVector x_Fuel;
	BzzVector x_Oxidizer;

	void Solve_TP();
	void Solve_HP();
	void Solve_SP();
	void Solve_TRHO();
	void Solve_URHO();
	void Solve_SRHO();

private:

	void AllProperties_0();
	void AllProperties_E();

    void ErrorMessage(const string message);
    void WarningMessage(const string message);
};

#endif // OPENSMOKE_EQUILIBRIUMMODULE_H