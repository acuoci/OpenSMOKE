/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_FLAMELET_DATAMANAGER)
#define OPENSMOKE_FLAMELET_DATAMANAGER

#include "OpenSMOKE.hpp"

class OpenSMOKE_Flamelet;

class OpenSMOKE_Dictionary_Flamelet : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_Flamelet();
};

class OpenSMOKE_Flamelet_DataManager  
{
public:

	OpenSMOKE_Flamelet_DataManager();
	void SetName(const string name);

	void Assign(OpenSMOKE_Flamelet *_flamelet);
	void Assign(OpenSMOKE_ReactingGas *_mix);
	void ReadFromFile(const string fileName);

	void DefineFromFile(const string inputFile);
	void Lock();

	// Pressure
	double P_atm;
	double P_bar;
	double P_Pascal;
	
	// Strain rate	
	double chiSt;
	double correctionChi;

	// Boundary Conditions
	double TC;
	double TO;

	// Information
	BzzVector XC;
	BzzVector XO;
	BzzVectorInt iXC;
	BzzVectorInt iXO;
	BzzVectorInt iOut;
	vector <string> nameOutput;
	string twoEquation_file_name;

	// Flame
	OpenSMOKE_Flamelet		*flamelet;
	OpenSMOKE_ReactingGas	*mix;
	int Np;


	// Print Information
	int iteration;
	int iterationVideoCounter;
	int iterationFileCounter;
	int nStepsVideo;
	int nStepsFile;
	int nStepsBackUp;


	// Information
	int		jFUEL;
	int		jO2;
	int		jINERT;
	string	nameFuel;
	string	nameOxidizer;
	string	nameInert;

	// Grid refining
	int		nDiff;
	int		nGrad;
	double	deltaDiff;
	double	deltaGrad;


	// User Defined Profiles
	double	Tpeak;
	double	xcen;
	BzzVector peak;


	// Miscellanea
	bool	iBackUp;
	int		iGlobalKinetics;

	// File Name
	BzzVector	listScalarDissipationRates;

	double enthalpyC;
	double enthalpyO;

	BzzVector x_elemental_C;
	BzzVector omega_elemental_C;
	BzzVector x_elemental_O;
	BzzVector omega_elemental_O;

	// Sensitivity Analisys
	BzzVectorInt	index_sensitivity;
	vector<string>	names_sensitivity;

	// Reaction Rates Analisys
	BzzVectorInt	index_reaction_rates;
	vector<string>	names_reaction_rates;

	// Formation Rates Analisys
	BzzVectorInt	index_formation_rates;
	vector<string>	names_formation_rates;

	// ROPA
	BzzVectorInt	index_ROPA;
	vector<string>	names_ROPA;

	// Additional variables
	bool iEnthalpy;
	bool iEnthalpyDefect;
	bool iDensityCorrection;
	bool i2E;

	int iRadiation;
	double environmentTemperature;
	double absTolerances;
	double relTolerances;
	double zStoichiometric;
	double enthalpyDefect;
	double minimumTemperature;

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

	string name_object;


	void AssignFuel(const string string_value);
	void AssignOxidizer(const string string_value);
	void AssignInert(const string string_value);
	void AssignPressure(const string units, const double value);
	void AssignFlameTemperature(const string units, const double value);
	void AssignFuelTemperature(const string units, const double value);
	void AssignOxidizerTemperature(const string units, const double value);
	void AssignFuelMassFractions(const vector<string> _names, const vector<double> _values);
	void AssignFuelMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignOxidizerMassFractions(const vector<string> _names, const vector<double> _values);
	void AssignOxidizerMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignGridPoints(const int int_value);
	void AssignOutputSpecies(const vector<string> string_vector);
	void AssignScalarDissipationRates(const vector<string> string_vector);

	void SetFlamePosition(const double double_value);
	void SetVideoSteps(const int int_value);
	void SetFileSteps(const int int_value);
	void SetBackupSteps(const int int_value);
	void SetRelativeTolerance(const double double_value);
	void SetAbsoluteTolerance(const double double_value);
	void SetEnvironmentTemperature(const string units, const double value);
	void SetEnthalpyDefect(const string units, const double value);
	void SetMinimumTemperature(const string units, const double value);
	void SetGasRadiation();
	void SetSootRadiation();
	void SetGridRefineGradient(const int int_value);
	void SetGridRefineCurvature(const int int_value);
	void SetGridRefineGradientStep(const double double_value);
	void SetGridRefineCurvatureStep(const double double_value);
	void SetDerivativeT(const char char_value);
	void SetDerivativeW(const char char_value);
	void SetSensitivityOnFile(const vector<string> _names);
	void SetFormationRatesOnFile(const vector<string> _names);
	void SetROPAOnFile(const vector<string> _names);
	void SetReactionRatesOnFile(const vector<string> _names);
	void SetDensityCorrection();

	bool iAssignedFlamePosition;
	bool iAssignedVideoSteps;
	bool iAssignedFileSteps;
	bool iAssignedBackupSteps;
	bool iAssignedRelativeTolerance;
	bool iAssignedAbsoluteTolerance;
	bool iAssignedEnvironmentTemperature;
	bool iAssignedGasRadiation;
	bool iAssignedSootRadiation;
	bool iAssignedGridRefineGradient;
	bool iAssignedGridRefineCurvature;
	bool iAssignedGridRefineGradientStep;
	bool iAssignedGridRefineCurvatureStep;
	bool iAssignedDerivativeT;
	bool iAssignedDerivativeW;
	bool iAssignedSensitivity;
	bool iAssignedFormationRates;
	bool iAssignedROPA;
	bool iAssignedReactionRates;

	void StoichiometricMixtureFraction();
};

#endif // DATAMANAGER

