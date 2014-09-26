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

#include <sstream>
#include <iomanip>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_WarningFile.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ReactionData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_KineticsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ThermoData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_UnitsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportData.h"

void CompactStoichiometries(BzzVectorInt &index, BzzVector &nu, int &n, vector<string> &names);
void CompleteGlobalStoichiometries(BzzVectorInt &index, BzzVectorInt &indexGlobal, BzzVector &nu, BzzVector &lambdaGlobal);

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_CHEMKINInterpreter_ReactionData"	<< endl;
    cout << "Error:    " << message         << endl;
    cout << "Reaction: " << iReaction       << endl;
    cout << "Line:     " << iLine           << endl;    
	cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_ReactionData"	<< endl;
    cout << "Warning:  " << message         << endl;
    cout << "Reaction: " << iReaction       << endl;
    cout << "Line:     " << iLine           << endl;    
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::Setup(const int _iReaction, const int _iLine, OpenSMOKE_CHEMKINInterpreter_KineticsData *_ptKinetics, OpenSMOKE_CHEMKINInterpreter_UnitsData  *_ptUnits)
{
	iReaction	= _iReaction;
	iLine		= _iLine;
	ptKinetics	= _ptKinetics;
	ptUnits		= _ptUnits;
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::Setup(const double _energy_factor, const double _frequency_factor)
{
	energy_factor = _energy_factor;
	frequency_factor = _frequency_factor;
}

OpenSMOKE_CHEMKINInterpreter_ReactionData::OpenSMOKE_CHEMKINInterpreter_ReactionData()
{
	A		= 0.;
	Beta	= 0.;
	E		= 0.;

	nReactants = 0;
	nameDirect.push_back("reactant names");

	nProducts = 0;
	nameInverse.push_back("product names");

	nGlobalReactants = 0;
	nameGlobalDirect.push_back("reactant names (global)");
	
	nGlobalProducts = 0;
	nameGlobalInverse.push_back("product names (global)");
	
	sumNu			= 0.;
	lambdaTotal		= 0.;

	iReversible		= false;
	iGlobal			= false;

	kindPressureDependence		= PRESSURE_LINDEMANN;

	iThirdBody					= false;
	
	iHighPressure				= false;	// cabr
	iLowPressure				= false;	// fall-off
	iLandauTeller				= false;
	iJanevLanger				= false;
	iPowerSeries				= false;
	iPressureLogarithmic		= false;
	iChebishevPolynomials		= false;
	iReverseRate				= false;
	iReverseRateLandauTeller	= false;
	iDuplicate					= false;
	iTAR						= false;

	iConversionFrequencyFactor	= false;
	iConversionEnergy			= false;
	iThermodynamicConsistency	= false;
	iCollisionEfficiency		= false;

	iReactionBis				= 0;
}

OpenSMOKE_CHEMKINInterpreter_ReactionData::~OpenSMOKE_CHEMKINInterpreter_ReactionData()
{
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ReverseRateFromReaction(const int _iReactionBis, OpenSMOKE_CHEMKINInterpreter_ReactionData &r)
{
	if (r.iHighPressure == true)			ErrorMessage("REV and HIGH options are ambigous. Please write the reaction as a couple of irreversible reactions...");
	if (r.iLowPressure  == true)			ErrorMessage("REV and LOW options are ambigous. Please write the reaction as a couple of irreversible reactions...");
	if (r.iChebishevPolynomials == true)	ErrorMessage("REV and CHEB options are ambigous. Please write the reaction as a couple of irreversible reactions...");
	if (r.iPressureLogarithmic == true)		ErrorMessage("REV and PLOG options are ambigous. Please write the reaction as a couple of irreversible reactions...");

	iReverseRateLandauTeller = r.iReverseRateLandauTeller;

	if (iReverseRateLandauTeller == true)
	{
		iLandauTeller = true;
		landau_teller = r.reverse_rate_landau_teller;
	}

	A		= r.reverse_rate[1];
	Beta	= r.reverse_rate[2];
	E		= r.reverse_rate[3];

	energy_factor		= r.energy_factor;
	frequency_factor	= r.frequency_factor;

	nReactants			= r.nProducts;
	indexDirect			= r.indexInverse;
	nuDirect			= r.nuInverse;
	lambdaGlobalDirect	= r.lambdaGlobalInverse;
	
	nProducts			= r.nReactants;
	indexInverse		= r.indexDirect;
	nuInverse			= r.nuDirect;
	lambdaGlobalInverse = r.lambdaGlobalDirect;

	int i;
	for(i=1;i<=nReactants;i++)
		nameDirect.push_back(r.nameInverse[i]);
	for(i=1;i<=nProducts;i++)
		nameInverse.push_back(r.nameDirect[i]);

	// Third-body
	iThirdBody = r.iThirdBody;
	efficiencies_coefficients = r.efficiencies_coefficients;
	efficiencies_indices = r.efficiencies_indices;

	// Power series and Janev-Langer
	iPowerSeries			= false;
	iJanevLanger			= false;
	iCollisionEfficiency	= false;
	iTAR					= false;

	iLine				= r.iLine;
	iReaction			= r.iReaction;
	iReactionBis		= _iReactionBis;
	r.iReverseReaction	= _iReactionBis;

	iGlobal					= r.iGlobal;
	iHighPressure			= r.iHighPressure;
	iLowPressure			= r.iLowPressure;
	iPressureLogarithmic	= r.iPressureLogarithmic;
	iChebishevPolynomials	= r.iChebishevPolynomials;

	iThermodynamicConsistency	= false;
	iReversible					= false;
	iReverseRate				= false;
	iReverseRateLandauTeller	= false;

	iConversionEnergy			= r.iConversionEnergy;
	iConversionFrequencyFactor	= r.iConversionFrequencyFactor;

	ptKinetics	= r.ptKinetics;
	ptUnits		= r.ptUnits;

	Complete();
	ReactionString();
	r.ReactionString();
}


void OpenSMOKE_CHEMKINInterpreter_ReactionData::Check()
{
	if (iHighPressure==true && iLowPressure==true)
		ErrorMessage("Only LOW or HIGH pressure limits have to be specified");
	
	if (kindPressureDependence == PRESSURE_TROE)
		if (iHighPressure==false && iLowPressure==false)
			ErrorMessage("LOW or HIGH Pressure limits must be specified for TROE pressure dependence");

	if (kindPressureDependence == PRESSURE_SRI)
		if (iHighPressure==false && iLowPressure==false)
			ErrorMessage("LOW or HIGH Pressure limits must be specified for SRI pressure dependence");


	if (iPressureLogarithmic == true)
	{
		BzzVector values = pressure_logarithmic.GetColumn(1);
		BzzVectorInt indices(values.Size());
		Sort(&values, &indices);
		ReorderByRows(&pressure_logarithmic, indices);
	}

	if (iChebishevPolynomials == true)
	{
		if (int(chebishev_vector[1])*int(chebishev_vector[2]) != chebishev_vector.Size()-2)
			ErrorMessage("Wrong number of arguments in CHEB option!");

		if (chebishev_pressure.Size() == 0)
		{
			ChangeDimensions(2, &chebishev_pressure);
			chebishev_pressure[1] = 0.001;		// atm
			chebishev_pressure[2] = 100.0;		// atm
		}

		if (chebishev_temperature.Size() == 0)
		{
			ChangeDimensions(2, &chebishev_temperature);
			chebishev_temperature[1] = 300.;	// K
			chebishev_temperature[2] = 2500.;	// K
		}

		if (iConversionEnergy == true || iConversionFrequencyFactor == true)
			ErrorMessage("UNITS option has no effect on CHEB reactions");
	}

	if (iReverseRateLandauTeller == true && iReverseRate == false)
		ErrorMessage("RLT option requires REV option");

	if (iLandauTeller == true && iReverseRate==true && iReverseRateLandauTeller==false)
		ErrorMessage("LT and REV options require RLT option");

	if (iGlobal == true)
	{
		int i;

		for(i=1;i<=indexGlobalDirect.Size();i++)
			for(int j=1;j<=indexGlobalDirect.Size();j++)
				if (i != j)
				{
					if (indexGlobalDirect[i] == indexGlobalDirect[j])
						ErrorMessage("FORD option is used more than once for the following species: " + nameGlobalDirect[i]);
				}
	}

	if (iGlobal == true && iReversible == true)
	{
		int i;
		BzzVector aux_lambda(ptKinetics->species_list.size()-1);
		
		for(i=1;i<=nProducts;i++)
			aux_lambda[indexInverse[i]] += nuInverse[i];

		if (iThermodynamicConsistency == true)
		{
			for(i=1;i<=nReactants;i++)
				aux_lambda[indexDirect[i]]  += -nuDirect[i];

			for(i=1;i<=indexGlobalDirect.Size();i++)
				aux_lambda[indexGlobalDirect[i]]  += lambdaGlobalDirect[i];
		}

		for(i=1;i<=aux_lambda.Size();i++)
			if (aux_lambda[i] != 0.)
			{
				indexGlobalInverse.Append(i);
				lambdaGlobalInverse.Append(aux_lambda[i]);
			}
	}

	if (iJanevLanger == true)
		WriteWarningMessageOnFile("Please note that Janev-Langer kinetics is implemented in a different way\nwith respect the CHEMKIN formulation.\nMore details on the OpenSMOKE User's Guide");

	if (iChebishevPolynomials == true)
		WriteWarningMessageOnFile("Chebishev Polynomials kinetics uses always the CHEMKIN\nstandard units: [mol, cm3,s] and [cal/mol]. Different units are not allowed.\nMore details on the OpenSMOKE User's Guide");

	if (iGlobal == true && iReversible == true && iThermodynamicConsistency == false)
		WriteWarningMessageOnFile("The reaction is reversible, but the thermodynamic consistency is not forced.\nThis could lead to unexpected and unphysical results.\nIt is recommended to use the CONSISTENCY option to force thermodynamic consistency.\nMore details on the OpenSMOKE User's Guide");

	if (iCollisionEfficiency == true)
		WriteWarningMessageOnFile("Bimolecular Collision Efficiency option (COLLEFF) needs to be checked.");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::CheckStoichiometry(OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, BzzVectorInt &indices, vector<string> &element_names, vector<string> &isotope_names)
{
	int nElements = thermo.species[indices[1]].element_indices.Size();
	int nIsotopes = thermo.species[indices[1]].isotope_indices.Size();

	// Elements
	{
		int i;	
		BzzVector sumReactants(nElements);
		BzzVector sumProducts(nElements);

		for(i=1;i<=nReactants;i++)
		{
			int h = indices[indexDirect[i]];
			for(int j=1;j<=nElements;j++)
				sumReactants[j] += thermo.species[h].element_indices[j] * nuDirect[i];
		}

		for(i=1;i<=nProducts;i++)
		{
			int h = indices[indexInverse[i]];
			for(int j=1;j<=nElements;j++)
				sumProducts[j] += thermo.species[h].element_indices[j] * nuInverse[i];
		}

		for(int j=1;j<=nElements;j++)
			if (sumReactants[j] > sumProducts[j]*1.0001 || sumReactants[j] < sumProducts[j]*0.9999)
			{
				stringstream sumR; sumR << sumReactants[j];
				stringstream sumP; sumP << sumProducts[j];
				std::string message = "Umbalanced stoichiometry on " + element_names[j] + " element\n";
				message += "  Reactants: " + sumR.str() + " - Products: " + sumP.str();
//				ErrorMessage(message);
				WarningMessage(message);
			}
	}

	// Isotopes
	if (nIsotopes !=0)
	{
		int i;
		BzzVector sumReactants(nIsotopes);
		BzzVector sumProducts(nIsotopes);

		for(i=1;i<=nReactants;i++)
		{
			int h = indices[indexDirect[i]];
			for(int j=1;j<=nIsotopes;j++)
				sumReactants[j] += thermo.species[h].isotope_indices[j] * nuDirect[i];
		}

		for(i=1;i<=nProducts;i++)
		{
			int h = indices[indexInverse[i]];
			for(int j=1;j<=nIsotopes;j++)
				sumProducts[j] += thermo.species[h].isotope_indices[j] * nuInverse[i];
		}

		for(int j=1;j<=nIsotopes;j++)
			if (sumReactants[j] != sumProducts[j])
//				ErrorMessage("Umbalanced stoichiometry on " + isotope_names[j] + " element");	
				WarningMessage("Umbalanced stoichiometry on " + isotope_names[j] + " element");	
	}
}


void OpenSMOKE_CHEMKINInterpreter_ReactionData::Chebishev(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Chebishev polynomials (CHEB) option";

	// Checking for syntax errors
	if (instructions[2] != "/")					ErrorMessage(error_message);
	if (instructions[nInstructions] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iHighPressure			== true)		ErrorMessage("CHEB and other pressure-depency (HIGH) are mutually exclusive...");
	if (iLowPressure			== true)		ErrorMessage("CHEB and other pressure-depency (LOW) are mutually exclusive...");
	if (iLandauTeller			== true)		ErrorMessage("CHEB and other temperature-depency (LT) are mutually exclusive...");
	if (iCollisionEfficiency	== true)		ErrorMessage("CHEB and other temperature-depency (COLLEFF) are mutually exclusive...");
	if (iJanevLanger			== true)		ErrorMessage("CHEB and other temperature-depency (JAN) are mutually exclusive...");
	if (iPowerSeries			== true)		ErrorMessage("CHEB and other temperature-depency (FIT1) are mutually exclusive...");
	if (iPressureLogarithmic	== true)		ErrorMessage("CHEB and other pressure-depency (PLOG) are mutually exclusive...");
	if (iTAR					== true)		ErrorMessage("CHEB and other temperature-depency (REACTAR) are mutually exclusive...");

	// Assign properties
	iChebishevPolynomials	= true;
	iThirdBody				= false;

	for(int i=1;i<=nInstructions-3;i++)
		chebishev_vector.Append(atof(instructions[2+i].c_str()));
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::CollisionalEfficiency(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Efficiency of Collision Frequency (COLLEFF) option";

	if (nInstructions != 1)			ErrorMessage(error_message);
	
	if (iHighPressure			== true)		ErrorMessage("COLLEFF and other pressure-depency (HIGH) are mutually exclusive...");
	if (iLowPressure			== true)		ErrorMessage("COLLEFF and other pressure-depency (LOW) are mutually exclusive...");
	if (iLandauTeller			== true)		ErrorMessage("COLLEFF and other temperature-depency (LT) are mutually exclusive...");
	if (iChebishevPolynomials	== true)		ErrorMessage("COLLEFF and other temperature-depency (CHEB) are mutually exclusive...");
	if (iJanevLanger			== true)		ErrorMessage("COLLEFF and other temperature-depency (JAN) are mutually exclusive...");
	if (iPowerSeries			== true)		ErrorMessage("COLLEFF and other temperature-depency (FIT1) are mutually exclusive...");
	if (iPressureLogarithmic	== true)		ErrorMessage("COLLEFF and other pressure-depency (PLOG) are mutually exclusive...");
	if (iGlobal					== true)		ErrorMessage("COLLEFF and non stoichiometric reaction orders (FORD) are mutually exclusive...");
	if (iTAR					== true)		ErrorMessage("COLLEFF and other temperature-depency (REACTAR) are mutually exclusive...");

	iCollisionEfficiency = true;
	
//	ErrorMessage("Collisional efficiencies (COLLEFF) not yet implemented in OpenSMOKE!");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::Duplicate(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Duplicate (DUP || DUPLICATE) option";

	// Checking for syntax errors
	if (nInstructions != 1)			ErrorMessage(error_message);
	if (iDuplicate == true)			ErrorMessage(error_message);

	iDuplicate = true;
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::EnergyLossParameter(vector<string> instructions)
{
	ErrorMessage("Energy Loss Parameter (EXCI) not yet implemented in OpenSMOKE!");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::Fit1(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Power Series (FIT1) option";

	// Checking for syntax errors
	if (nInstructions != 7)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[7] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iPowerSeries			== true)		ErrorMessage("FIT1 option is used more than once");
	if (iJanevLanger			== true)		ErrorMessage("FIT1 and other temperature-depency (JAN) are mutually exclusive...");
	if (iLandauTeller			== true)		ErrorMessage("FIT1 and other temperature-depency (LT) are mutually exclusive...");
	if (iPressureLogarithmic	== true)		ErrorMessage("FIT1 and other pressure-depency (PLOG) are mutually exclusive...");
	if (iChebishevPolynomials	== true)		ErrorMessage("FIT1 and other pressure-depency (CHEB) are mutually exclusive...");
	if (iHighPressure			== true)		ErrorMessage("FIT1 and other pressure-depency (HIGH) are mutually exclusive...");
	if (iLowPressure			== true)		ErrorMessage("FIT1 and other pressure-depency (LOW) are mutually exclusive...");
	if (iCollisionEfficiency	== true)		ErrorMessage("FIT1 and other temperature-depency (COLLEFF) are mutually exclusive...");
	if (iTAR					== true)		ErrorMessage("FIT1 and other temperature-depency (REACTAR) are mutually exclusive...");

	// Additional checks
	if (E != 0.)	ErrorMessage("The reaction energy must be equal to zero for Power-Series reactions");	
	
	// Assign properties
	iPowerSeries = true;
	ChangeDimensions(4, &power_series);
	power_series[1] = atof(instructions[3].c_str());
	power_series[2] = atof(instructions[4].c_str());
	power_series[3] = atof(instructions[5].c_str());
	power_series[4] = atof(instructions[6].c_str());
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::TAR(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in TAR Reaction (TAR) option";

	// Checking for syntax errors
	if (nInstructions != 9)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[9] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iPowerSeries			== true)		ErrorMessage("REACTAR option is used more than once");
	if (iJanevLanger			== true)		ErrorMessage("REACTAR and other temperature-depency (JAN) are mutually exclusive...");
	if (iLandauTeller			== true)		ErrorMessage("REACTAR and other temperature-depency (LT) are mutually exclusive...");
	if (iPressureLogarithmic	== true)		ErrorMessage("REACTAR and other pressure-depency (PLOG) are mutually exclusive...");
	if (iChebishevPolynomials	== true)		ErrorMessage("REACTAR and other pressure-depency (CHEB) are mutually exclusive...");
	if (iHighPressure			== true)		ErrorMessage("REACTAR and other pressure-depency (HIGH) are mutually exclusive...");
	if (iLowPressure			== true)		ErrorMessage("REACTAR and other pressure-depency (LOW) are mutually exclusive...");
	if (iCollisionEfficiency	== true)		ErrorMessage("REACTAR and other temperature-depency (COLLEFF) are mutually exclusive...");

	// Additional checks
	if (E != 1.)	ErrorMessage("The reaction energy must be equal to 1 for TAR reactions");	
	if (A != 1.)	ErrorMessage("The reaction frequency factor must be equal to 1 for TAR reactions");	
	if (Beta != 1.)	ErrorMessage("The reaction temperature exponent must be equal to 1 for TAR reactions");	
	
	// Assign properties
	iTAR = true;
	ChangeDimensions(6, &TAR_series);
	TAR_series[1] = atof(instructions[3].c_str());
	TAR_series[2] = atof(instructions[4].c_str());
	TAR_series[3] = atof(instructions[5].c_str());
	TAR_series[4] = atof(instructions[6].c_str());
	TAR_series[5] = atof(instructions[7].c_str());
	TAR_series[6] = atof(instructions[8].c_str());
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ForwardReactionKineticParameter(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in forward reaction orders (FORD)";

	// Checking for syntax errors
	if (nInstructions != 5)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[5] != "/")		ErrorMessage(error_message);

	iGlobal = true;
	lambdaGlobalDirect.Append(atof(instructions[4].c_str()));
	indexGlobalDirect.Append(ptKinetics->RecognizeSpecies(instructions[3]));
	nameGlobalDirect.push_back(instructions[3]);
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::HighPressureLimit(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in High-pressure limits (HIGH) option";

	// Checking for syntax errors
	if (nInstructions != 6)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[6] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iHighPressure			== true)		ErrorMessage("HIGH option is used more than once");
	if (iLowPressure			== true)		ErrorMessage("HIGH and other pressure-depency (LOW) are mutually exclusive...");
	if (iPressureLogarithmic	== true)		ErrorMessage("HIGH and other pressure-depency (PLOG) are mutually exclusive...");
	if (iChebishevPolynomials	== true)		ErrorMessage("HIGH and other pressure-depency (CHEB) are mutually exclusive...");
	if (iLandauTeller			== true)		ErrorMessage("HIGH and other temperature-depency (LT) are mutually exclusive...");
	if (iJanevLanger			== true)		ErrorMessage("HIGH and other temperature-depency (JAN) are mutually exclusive...");
	if (iPowerSeries			== true)		ErrorMessage("HIGH and other temperature-depency (FIT1) are mutually exclusive...");
	if (iCollisionEfficiency	== true)		ErrorMessage("HIGH and other temperature-depency (COLLEFF) are mutually exclusive...");
	if (iTAR					== true)		ErrorMessage("HIGH and other temperature-depency (REACTAR) are mutually exclusive...");

	// Assign properties
	iHighPressure	= true;
	ChangeDimensions(3, &high_pressure);
	ChangeDimensions(3, &low_pressure);
	high_pressure[1] = atof(instructions[3].c_str());
	high_pressure[2] = atof(instructions[4].c_str());
	high_pressure[3] = atof(instructions[5].c_str());
	low_pressure[1]  = A;
	low_pressure[2]  = Beta;
	low_pressure[3]  = E;

	// Single third body
	if (thirdbodysingle != "null")	ThirdBodySingle();
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::JanevLangerReactionRate(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Janev-Langer Reactions (JAN) option";

	// Checking for syntax errors
	if (nInstructions != 12)		ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[12] != "/")	ErrorMessage(error_message);

	// Checking for duplicates
	if (iJanevLanger			== true)		ErrorMessage("JAN option is used more than once");
	if (iLandauTeller			== true)		ErrorMessage("JAN and other pressure-depency (LT) are mutually exclusive...");
	if (iPowerSeries			== true)		ErrorMessage("JAN and other pressure-depency (FIT1) are mutually exclusive...");
	if (iPressureLogarithmic	== true)		ErrorMessage("JAN and other pressure-depency (PLOG) are mutually exclusive...");
	if (iChebishevPolynomials	== true)		ErrorMessage("JAN and other pressure-depency (CHEB) are mutually exclusive...");
	if (iHighPressure			== true)		ErrorMessage("JAN and other pressure-depency (HIGH) are mutually exclusive...");
	if (iLowPressure			== true)		ErrorMessage("JAN and other pressure-depency (LOW) are mutually exclusive...");
	if (iCollisionEfficiency	== true)		ErrorMessage("JAN and other temperature-depency (COLLEFF) are mutually exclusive...");
	if (iTAR					== true)		ErrorMessage("JAN and other temperature-depency (REACTAR) are mutually exclusive...");

	// Assign properties
	iJanevLanger = true;
	ChangeDimensions(9, &janev_langer);
	for(int i=1;i<=9;i++)
		janev_langer[i] = atof(instructions[2+i].c_str());

	// Single third body
	if (thirdbodysingle != "null")	ErrorMessage("Single third body species can be used only with LOW || HIGH options");
	
	//ErrorMessage("Janev-Langer reaction orders (JAN) not yet implemented in OpenSMOKE!");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::LowPressureLimit(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Low-pressure limits (LOW) option";

	// Checking for syntax errors
	if (nInstructions != 6)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[6] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iLowPressure			== true)		ErrorMessage("LOW option is used more than once");
	if (iHighPressure			== true)		ErrorMessage("LOW and other pressure-depency (HIGH) are mutually exclusive...");
	if (iPressureLogarithmic	== true)		ErrorMessage("LOW and other pressure-depency (PLOG) are mutually exclusive...");
	if (iChebishevPolynomials	== true)		ErrorMessage("LOW and other pressure-depency (CHEB) are mutually exclusive...");
	if (iLandauTeller			== true)		ErrorMessage("LOW and other temperature-depency (LT) are mutually exclusive...");
	if (iJanevLanger			== true)		ErrorMessage("LOW and other temperature-depency (JAN) are mutually exclusive...");
	if (iPowerSeries			== true)		ErrorMessage("LOW and other temperature-depency (FIT1) are mutually exclusive...");
	if (iCollisionEfficiency	== true)		ErrorMessage("LOW and other temperature-depency (COLLEFF) are mutually exclusive...");
	if (iTAR					== true)		ErrorMessage("LOW and other temperature-depency (REACTAR) are mutually exclusive...");

	// Assign properties
	iLowPressure = true;
	ChangeDimensions(3, &low_pressure);
	ChangeDimensions(3, &high_pressure);
	low_pressure[1] = atof(instructions[3].c_str());
	low_pressure[2] = atof(instructions[4].c_str());
	low_pressure[3] = atof(instructions[5].c_str());
	high_pressure[1] = A;
	high_pressure[2] = Beta;
	high_pressure[3] = E;

	// Single third body
	if (thirdbodysingle != "null")	ThirdBodySingle();
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::LandauTellerReaction(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Landau-Teller (LT) option";

	// Checking for syntax errors
	if (nInstructions != 5)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[5] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iLandauTeller			== true)		ErrorMessage("LT option is used more than once");
	if (iHighPressure			== true)		ErrorMessage("LT and other pressure-depency (HIGH) are mutually exclusive...");
	if (iLowPressure			== true)		ErrorMessage("LT and other pressure-depency (LOW) are mutually exclusive...");
	if (iJanevLanger			== true)		ErrorMessage("LT and other pressure-depency (JAN) are mutually exclusive...");
	if (iPowerSeries			== true)		ErrorMessage("LT and other pressure-depency (FIT1) are mutually exclusive...");
	if (iPressureLogarithmic	== true)		ErrorMessage("LT and other pressure-depency (PLOG) are mutually exclusive...");
	if (iChebishevPolynomials	== true)		ErrorMessage("LT and other pressure-depency (CHEB) are mutually exclusive...");
	if (iCollisionEfficiency	== true)		ErrorMessage("LT and other temperature-depency (COLLEFF) are mutually exclusive...");
	if (iTAR					== true)		ErrorMessage("LT and other temperature-depency (REACTAR) are mutually exclusive...");

	// Assign properties
	iLandauTeller = true;
	ChangeDimensions(2, &landau_teller);
	landau_teller[1] = atof(instructions[3].c_str());
	landau_teller[2] = atof(instructions[4].c_str());
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::MomentumTransferCollisionFrequency(vector<string> instructions)
{
	ErrorMessage("Momentum Transfer Collision Frequency (MOME) not yet implemented in OpenSMOKE!");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ChebishevPressureLimits(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Chebishev Pressure Limits (PCHEB) option";

	if (nInstructions == 5)
	{		
		if (instructions[2] != "/")				ErrorMessage(error_message);
		if (instructions[5] != "/")				ErrorMessage(error_message);

		// Assign properties
		ChangeDimensions(2, &chebishev_pressure);
		chebishev_pressure[1] = atof(instructions[3].c_str());
		chebishev_pressure[2] = atof(instructions[4].c_str());

		if (chebishev_pressure[1] >= chebishev_pressure[2])		ErrorMessage(error_message);
		if (chebishev_pressure[1] <=0. || chebishev_pressure[2] <=0.)	ErrorMessage(error_message);
	}

	else if (nInstructions == 10)
	{		
		if (instructions[6] != "TCHEB")				ErrorMessage(error_message);
		
		vector<string> instructions_additional;
		instructions_additional.push_back("instructions_additional");
		for(int i=1;i<=5;i++)
			instructions_additional.push_back(instructions[5+i]);
		ChebishevTemperatureLimits(instructions_additional);
	}

	else	ErrorMessage(error_message);



}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ChebishevTemperatureLimits(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Chebishev Temperature Limits (TCHEB) option";	

	if (nInstructions == 5)
	{		
		if (instructions[2] != "/")				ErrorMessage(error_message);
		if (instructions[5] != "/")				ErrorMessage(error_message);

		// Assign properties
		ChangeDimensions(2, &chebishev_temperature);
		chebishev_temperature[1] = atof(instructions[3].c_str());
		chebishev_temperature[2] = atof(instructions[4].c_str());

		if (chebishev_temperature[1] >= chebishev_temperature[2])			ErrorMessage(error_message);
		if (chebishev_temperature[1] <=0. || chebishev_temperature[2] <=0.)		ErrorMessage(error_message);
	}

	else if (nInstructions == 10)
	{		
		if (instructions[6] != "PCHEB")				ErrorMessage(error_message);
		
		vector<string> instructions_additional;
		instructions_additional.push_back("instructions_additional");
		for(int i=1;i<=5;i++)
			instructions_additional.push_back(instructions[5+i]);
		ChebishevPressureLimits(instructions_additional);
	}

	else	ErrorMessage(error_message);

}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::PressureLogarithmicInterpolation(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Pressure Logarithmic Interpolation (PLOG) option";

	// Checking for syntax errors
	if (nInstructions != 7)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[7] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iThirdBody				== true)		ErrorMessage("PLOG and other pressure-depency (THREE-BODY) are mutually exclusive...");
	if (iHighPressure			== true)		ErrorMessage("PLOG and other pressure-depency (HIGH) are mutually exclusive...");
	if (iLowPressure			== true)		ErrorMessage("PLOG and other pressure-depency (LOW) are mutually exclusive...");
	if (iLandauTeller			== true)		ErrorMessage("PLOG and other temperature-depency (LT) are mutually exclusive...");
	if (iJanevLanger			== true)		ErrorMessage("PLOG and other temperature-depency (JAN) are mutually exclusive...");
	if (iPowerSeries			== true)		ErrorMessage("PLOG and other temperature-depency (FIT1) are mutually exclusive...");
	if (iChebishevPolynomials	== true)		ErrorMessage("PLOG and other pressure-depency (CHEB) are mutually exclusive...");
	if (iCollisionEfficiency	== true)		ErrorMessage("PLOG and other temperature-depency (COLLEFF) are mutually exclusive...");
	if (iTAR					== true)		ErrorMessage("PLOG and other temperature-depency (REACTAR) are mutually exclusive...");

	// Assign properties
	iPressureLogarithmic = true;

	BzzVector plog(4);
	plog[1] = atof(instructions[3].c_str());
	plog[2] = atof(instructions[4].c_str());
	plog[3] = atof(instructions[5].c_str());
	plog[4] = atof(instructions[6].c_str());
	pressure_logarithmic.AppendRow(plog);
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ReverseRate(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Reverse Rate (REV) option";

	// Checking for syntax errors
	if (nInstructions != 6)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[6] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iReversible == false)		ErrorMessage("REV and irreversible reactions are mutually exclusive...");

	// Assign properties
	iReverseRate = true;

	ChangeDimensions(3, &reverse_rate);
	reverse_rate[1] = atof(instructions[3].c_str());
	reverse_rate[2] = atof(instructions[4].c_str());
	reverse_rate[3] = atof(instructions[5].c_str());
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::LandauTellerReverseReaction(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Reverse Landau Teller reactions (RLT) option";

	// Checking for syntax errors
	if (nInstructions != 5)			ErrorMessage(error_message);
	if (instructions[2] != "/")		ErrorMessage(error_message);
	if (instructions[5] != "/")		ErrorMessage(error_message);

	// Checking for duplicates
	if (iReversible == false)		ErrorMessage("The RLT option cannot be used for irreversible reactions");
	if (iReverseRate == false)		ErrorMessage("The RLT option cannot be used without the REV option");

	// Assign properties
	iReverseRateLandauTeller = true;

	ChangeDimensions(3, &reverse_rate_landau_teller);
	reverse_rate_landau_teller[1] = atof(instructions[3].c_str());
	reverse_rate_landau_teller[2] = atof(instructions[4].c_str());
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::BackwardReactionKineticParameter(vector<string> instructions)
{
	ErrorMessage("Backward reaction orders (RORD) not yet implemented in OpenSMOKE!");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::SRIFallOffReaction(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in SRI parameters (SRI) option";

	// Checking for syntax errors
	if (nInstructions != 6 && nInstructions != 8)	ErrorMessage(error_message);
	if (instructions[2] != "/")						ErrorMessage(error_message);
	if (instructions[nInstructions] != "/")			ErrorMessage(error_message);

	// Checking for duplicates
	if (kindPressureDependence == PRESSURE_SRI)		ErrorMessage("The SRI option is used more than once");

	// Assign properties
	kindPressureDependence  = PRESSURE_SRI;
	ChangeDimensions(nInstructions-3, &sri);
	sri[1] = atof(instructions[3].c_str());
	sri[2] = atof(instructions[4].c_str());
	sri[3] = atof(instructions[5].c_str());
	if (nInstructions == 8)	
	{
		sri[4] = atof(instructions[6].c_str());
		sri[5] = atof(instructions[7].c_str());
	}
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::SpeciesTemperatureDependence(vector<string> instructions)
{
	ErrorMessage("Species temperature dependence (TDEP) not yet implemented in OpenSMOKE!");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ThermodynamicConsistency(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in Thermodynamic Consistency (CONSISTENCY) option";

	// Checking for syntax errors
	if (nInstructions != 1)	ErrorMessage(error_message);

	// Checking for duplicates
	if (iGlobal	== false)		ErrorMessage("CONSISTENCY option can be used only for global reactions...");

	// Assign properties
	iThermodynamicConsistency = true;
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::TROEFallOffReaction(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in TROE parameters (TROE) option";

	// Checking for syntax errors
	if (nInstructions != 6 && nInstructions != 7)	ErrorMessage(error_message);
	if (instructions[2] != "/")						ErrorMessage(error_message);
	if (instructions[nInstructions] != "/")			ErrorMessage(error_message);

	// Checking for duplicates
	if (kindPressureDependence == PRESSURE_TROE)		ErrorMessage("The TROE option is used more than once");

	// Assign properties
	kindPressureDependence  = PRESSURE_TROE;
	ChangeDimensions(nInstructions-3, &troe);
	troe[1] = atof(instructions[3].c_str());
	troe[2] = atof(instructions[4].c_str());
	troe[3] = atof(instructions[5].c_str());
	if (nInstructions==7)	troe[4] = atof(instructions[6].c_str());
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::Units(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in UNITS option";

	// Checking for syntax errors
	if (nInstructions != 4)	ErrorMessage(error_message);

	double conversion_factor;
	bool iEnergy;
	ptUnits->GiveMeConversionFactor(instructions[3], conversion_factor, iEnergy);

	if (iEnergy == true) 
	{
		if (iConversionEnergy == true)	ErrorMessage("Activation energy UNITS are assigned more than once");
		energy_factor		= conversion_factor;
		iConversionEnergy   = true;
	}
	else				
	{
		if (iConversionFrequencyFactor == true)	ErrorMessage("Frequency factor UNITS are assigned more than once");
		frequency_factor				= conversion_factor;
		iConversionFrequencyFactor		= true;
	}
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::UserRateSubRoutine(vector<string> instructions)
{
	ErrorMessage("User Subroutines (USRPROG) not supported in OpenSMOKE!");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::CollisionCrossSection(vector<string> instructions)
{
	ErrorMessage("Collision cross section (XSMI) not yet implemented in OpenSMOKE!");
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ThirdBodyEfficiencies(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	std::string error_message = "Syntax error in third-body efficiencies";

	// Checking for syntax errors
	if (nInstructions%4 != 0)	ErrorMessage(error_message);

	// Checking for precedences
	if (iThirdBody == false)			ErrorMessage("The reaction is not a THREE-BODY reaction");	
	if (iChebishevPolynomials == true)	ErrorMessage("THREE-BODY and other pressure-depency (CHEB) are mutually exclusive...");	
	if (iPressureLogarithmic  == true)	ErrorMessage("THREE-BODY and other pressure-depency (PLOG) are mutually exclusive...");	
	if (iCollisionEfficiency  == true)	ErrorMessage("THREE-BODY and other temperature-depency (COLLEFF) are mutually exclusive...");

/*	ChangeDimensions(nInstructions/4, &efficiencies_coefficients);
	ChangeDimensions(nInstructions/4, &efficiencies_indices);
	for(int j=1;j<=nInstructions/4;j++)
	{
		if (instructions[(j-1)*4+2] != "/")		ErrorMessage(error_message);
		if (instructions[(j-1)*4+4] != "/")		ErrorMessage(error_message);

		efficiencies_coefficients[j] = atof(instructions[(j-1)*4+3].c_str());
		efficiencies_indices[j] = ptKinetics->RecognizeSpecies(instructions[(j-1)*4+1]);
	}
*/
	for(int j=1;j<=nInstructions/4;j++)
	{
		if (instructions[(j-1)*4+2] != "/")		ErrorMessage(error_message);
		if (instructions[(j-1)*4+4] != "/")		ErrorMessage(error_message);

		efficiencies_coefficients.Append(atof(instructions[(j-1)*4+3].c_str()));
		efficiencies_indices.Append(ptKinetics->RecognizeSpecies(instructions[(j-1)*4+1]));
	}

}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ThirdBodySingle()
{
	int j = ptKinetics->RecognizeSpecies(thirdbodysingle);
	iThirdBody = true;
	ChangeDimensions(1, &efficiencies_coefficients);
	ChangeDimensions(1, &efficiencies_indices);
	efficiencies_coefficients[1] = 1.;
	efficiencies_indices[1] = j;
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::ReactionString()
{
	int i;
	reaction_string_clean = "";
	
	for(i=1;i<=nReactants;i++)
	{
		if (nuDirect[i]!=1.)	
		{
			stringstream number;
			number << nuDirect[i];
			reaction_string_clean += number.str();
		}
		reaction_string_clean += nameDirect[i];
		if (i<nReactants)		reaction_string_clean += "+";
	}
	if (iThirdBody==true)		reaction_string_clean += "+M";

	if (iReversible==true)		reaction_string_clean += "=";
	else						reaction_string_clean += "=>";

	for(i=1;i<=nProducts;i++)
	{
		if (nuInverse[i]!=1.)
		{
			stringstream number;
			number << nuInverse[i];
			reaction_string_clean += number.str();
		}
		reaction_string_clean += nameInverse[i];
		if (i<nProducts)			reaction_string_clean += "+";
	}
	if (iThirdBody==true)		reaction_string_clean += "+M";

	sumNu = 0.;
	for(i=1;i<=nProducts;i++)	sumNu += nuInverse[i];
	for(i=1;i<=nReactants;i++)	sumNu -= nuDirect[i];

	reaction_string_complete = reaction_string_clean;
	
	if(iReactionBis !=0)	
	{
		stringstream index;
		index << iReaction;
		reaction_string_complete  += "#";
		reaction_string_complete += index.str();
		reaction_string_complete += "bis";
		reaction_string_complete += "#";
	}
	{
		std::string dummy;
		stringstream index;
		if(iReactionBis != 0)	index << iReactionBis;
		else					index << iReaction;

		dummy  = "#";
		dummy += index.str();
		dummy += "#";

		reaction_string_complete.insert(0,dummy);
	}

}


void OpenSMOKE_CHEMKINInterpreter_ReactionData::Complete()
{
	int i;

	if (iGlobal == true)
		CompleteGlobalStoichiometries(indexDirect, indexGlobalDirect, nuDirect, lambdaGlobalDirect);

	if (iGlobal == false)	for(i=1;i<=nReactants;i++)					lambdaTotal += nuDirect[i];
	else					for(i=1;i<=indexGlobalDirect.Size();i++)	lambdaTotal += lambdaGlobalDirect[i];
	
	
	if (iThirdBody==true)	lambdaTotal += 1.0;

		 if (lambdaTotal == 1.)		conversion_frequency_factor = 1.e0;
	else if (lambdaTotal == 2.)		conversion_frequency_factor = 1.e3;
	else if (lambdaTotal == 3.)		conversion_frequency_factor = 1.e6;
	else							conversion_frequency_factor = pow(1.e3, lambdaTotal-1.);
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::PrintOnFile(ofstream &fOutput)
{
	// Provisional TODO
	//fOutput << reaction_string_complete << "\t" << sumNu << "\t" << lambdaTotal << endl;;

	fOutput << reaction_string_complete << endl;

	if(iReversible == true)	fOutput << "1\t\t";
	else					fOutput << "0\t\t";

		 if (iThirdBody   == false && iLandauTeller == true)		fOutput << "100\t\t";
	else if (iThirdBody   == true  && iLandauTeller == true)		fOutput << "10\t\t";

	else if (iThirdBody == false && iJanevLanger == true)			fOutput << "110\t\t";
	else if (iThirdBody == true  && iJanevLanger == true)			fOutput << "11\t\t";
	else if (iThirdBody == false && iPowerSeries == true)			fOutput << "120\t\t";
	else if (iThirdBody == true  && iPowerSeries == true)			fOutput << "12\t\t";
	else if (iThirdBody == false && iChebishevPolynomials == true)	fOutput << "130\t\t";
	else if (iThirdBody == true  && iChebishevPolynomials == true)	ErrorMessage("CHEBISHEV and THREE-BODY option are mutually exclusive...");
	else if (iThirdBody == false && iPressureLogarithmic == true)	fOutput << "140\t\t";
	else if (iThirdBody == true  && iPressureLogarithmic == true)	ErrorMessage("PLOG and THREE-BODY option are mutually exclusive...");
	else if (iThirdBody == false && iCollisionEfficiency == true)	fOutput << "150\t\t";
	else if (iThirdBody == true  && iCollisionEfficiency == true)	ErrorMessage("COLLEFF and THREE-BODY option are mutually exclusive...");
	else if (iThirdBody == false && iTAR == true)					fOutput << "160\t\t";
	else if (iThirdBody == true && iTAR == true)					ErrorMessage("REACTAR and THREE-BODY option are mutually exclusive...");
	
	else if(iLowPressure==true)			// Fall-Off Reactions (2,3,4)
	{
		if (kindPressureDependence == PRESSURE_LINDEMANN)	fOutput << "2\t\t";
		if (kindPressureDependence == PRESSURE_TROE)		fOutput << "3\t\t";
		if (kindPressureDependence == PRESSURE_SRI)			fOutput << "4\t\t";
	}
	else if(iHighPressure==true)	// Chemically-Activated Bimolecular Reactions (5,6,7)
	{
		if (kindPressureDependence == PRESSURE_LINDEMANN)	fOutput << "5\t\t";
		if (kindPressureDependence == PRESSURE_TROE)		fOutput << "6\t\t";
		if (kindPressureDependence == PRESSURE_SRI)			fOutput << "7\t\t";
	}
	
	else
	{
		if (iThirdBody == true)					fOutput << "1\t\t";
		if (iThirdBody == false)				fOutput << "0\t\t";
	}

	if (efficiencies_indices.Size()==0)			fOutput << "0\t\t";
	else										fOutput << efficiencies_indices.Size();
	
	fOutput << endl;

	if (efficiencies_indices.Size()!=0)	
	{
		for(int j=1;j<=efficiencies_indices.Size();j++)
			fOutput << efficiencies_indices[j]		<< "\t" 
					<< efficiencies_coefficients[j]	<< "\t";
		fOutput << endl;
	}

	if (iLowPressure == true)	
	{
		fOutput << low_pressure[1]/pow(1.e3, lambdaTotal-1.)*pow(frequency_factor, lambdaTotal-1.)	<< "\t\t" 
				<< low_pressure[2]												<< "\t\t" 
				<< low_pressure[3]*energy_factor								<< endl;

		fOutput << high_pressure[1]/pow(1.e3, lambdaTotal-2.)*pow(frequency_factor, lambdaTotal-2.)	<< "\t\t" 
				<< high_pressure[2]												<< "\t\t" 
				<< high_pressure[3]*energy_factor								<< endl;
	}
	else if (iHighPressure == true)	
	{
		fOutput << low_pressure[1]/pow(1.e3, lambdaTotal-2.)*pow(frequency_factor, lambdaTotal-2.)	<< "\t\t" 
				<< low_pressure[2]												<< "\t\t" 
				<< low_pressure[3]*energy_factor								<< endl;

		fOutput << high_pressure[1]/pow(1.e3, lambdaTotal-3.)*pow(frequency_factor, lambdaTotal-3.)	<< "\t\t" 
				<< high_pressure[2]												<< "\t\t" 
				<< high_pressure[3]*energy_factor								<< endl;
	}
	else	// Conventional reactions
	{
		fOutput << A/conversion_frequency_factor*pow(frequency_factor, lambdaTotal-1.)	<< "\t\t" 
				<< Beta																	<< "\t\t" 
				<< E*energy_factor														<< endl;
	}

	// Fall-Off Reactions - CABR
	if(iLowPressure==true || iHighPressure==true)
	{
		if (kindPressureDependence == PRESSURE_TROE)
		{
			fOutput << troe[1] << "\t\t" << troe[2] << "\t\t" << troe[3] << "\t\t"; 
			if (troe.Size()==3)	fOutput << 0		<< "\t\t" << 0 << endl;
			if (troe.Size()==4)	fOutput << troe[4]	<< "\t\t" << 0 << endl;
		}
		
		if (kindPressureDependence == PRESSURE_SRI)
		{
			fOutput << sri[1] << "\t\t" << sri[2] << "\t\t" << sri[3] << "\t\t"; 
			if (sri.Size()==3)	fOutput << 1.0		<< "\t\t" << 0		<< endl;
			if (sri.Size()==5)	fOutput << sri[4]	<< "\t\t" << sri[5] << endl;
		}

		if (kindPressureDependence == PRESSURE_LINDEMANN)
		{
			fOutput << 0 << "\t\t" << 0 << "\t\t" << 0 << "\t\t" << 0 << "\t\t" << 0 << endl;
		}
	}

	// Landau-Teller
	if (iLandauTeller == true)
		fOutput << landau_teller[1] << "\t\t" << landau_teller[2] << endl;

	// Janev-Langer (TODO units)
	if (iJanevLanger == true)
	{	
		for (int j=1;j<=9;j++)
			fOutput << janev_langer[j] << "\t";
		fOutput << endl;
	}

	// Power-Series
	if (iPowerSeries == true)
	{	
		for (int j=1;j<=4;j++)
			fOutput << power_series[j] << "\t";
		fOutput << endl;
	}

	// Tar-Series
	if (iTAR == true)
	{	
		for (int j=1;j<=6;j++)
			fOutput << TAR_series[j] << "\t";
		fOutput << endl;
	}

	// Chebishev
	if (iChebishevPolynomials == true)
	{	
		fOutput << 	int(chebishev_vector[1]) << "\t";
		fOutput << 	int(chebishev_vector[2]) << endl;
		int m=3;
		for (int j=1;j<=int(chebishev_vector[1]);j++)
		{
			for (int k=1;k<=int(chebishev_vector[2]);k++)
				fOutput << 	chebishev_vector[m++] << "\t";
			fOutput << endl;
		}
		// The units are always [cm, mol] (so, no additional conversion factor is required)
		fOutput << conversion_frequency_factor	<< endl;	// TODO (units)
		fOutput << chebishev_temperature[1]		<< "\t" << chebishev_temperature[2] << endl;
		fOutput << chebishev_pressure[1]		<< "\t" << chebishev_pressure[2]	<< endl;
	}

	// Logarithmic-Pressure Dependence
	if (iPressureLogarithmic == true)
	{	
		int j;
		fOutput << 	pressure_logarithmic.Rows() << endl;
		for (j=1;j<=pressure_logarithmic.Rows();j++)
			fOutput << pressure_logarithmic[j][1] << "\t";
		fOutput << endl;
		for (j=1;j<=pressure_logarithmic.Rows();j++)
			fOutput << pressure_logarithmic[j][2]/conversion_frequency_factor*pow(frequency_factor, lambdaTotal-1.) << "\t";
		fOutput << endl;
		for (j=1;j<=pressure_logarithmic.Rows();j++)
			fOutput << pressure_logarithmic[j][3] << "\t";
		fOutput << endl;
		for (j=1;j<=pressure_logarithmic.Rows();j++)
			fOutput << pressure_logarithmic[j][4]*energy_factor << "\t";
		fOutput << endl;
	}

	// Collision Efficiency
	if (iCollisionEfficiency == true)
		fOutput << kStar_collision_frequency;
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::SummaryOnFile(ofstream &fOutput)
{
	int i;
	if (iReactionBis != 0)
		fOutput << setw(7) << right << iReactionBis;
	else
		fOutput << setw(7) << right << iReaction;
	fOutput << ". ";
	fOutput << reaction_string_clean << endl;

	if (iLowPressure == true)	
	{
		fOutput << setw(9) << " ";
		fOutput << "Fall-Off Reaction" << endl;
		
		fOutput << setw(9) << " ";
		fOutput << setw(9) << left << "k0:";
		fOutput << scientific << setprecision(6) << right << low_pressure[1]/pow(1.e3, lambdaTotal-1.)*pow(frequency_factor, lambdaTotal-1.)	<< "\t";
		fOutput << setw(8)    << setprecision(2) << fixed << right << low_pressure[2];
		fOutput << setw(14)	  << setprecision(2) << fixed << right << low_pressure[3]*energy_factor   << endl;

		fOutput << setw(9) << " ";
		fOutput << setw(9) << left << "kInf:";
		fOutput << scientific << setprecision(6) << right << high_pressure[1]/pow(1.e3, lambdaTotal-2.)*pow(frequency_factor, lambdaTotal-2.)	<< "\t"; 
		fOutput << setw(8)    << setprecision(2) << fixed << right << high_pressure[2];
		fOutput << setw(14)	  << setprecision(2) << fixed << right << high_pressure[3]*energy_factor << endl;
	}
	else if (iHighPressure == true)	
	{
		fOutput << setw(9) << " ";
		fOutput << "Chemically Activated Bimolecular Reaction" << endl;

		fOutput << setw(9) << " ";
		fOutput << setw(9) << left << "k0:";
		fOutput << scientific << setprecision(6) << right << low_pressure[1]/pow(1.e3, lambdaTotal-2.)*pow(frequency_factor, lambdaTotal-2.)	<< "\t";
		fOutput << setw(8) << setprecision(2) << fixed << right << low_pressure[2];
		fOutput << setw(14) << setprecision(2) << fixed << right << low_pressure[3]*energy_factor  << endl;

		fOutput << setw(9) << " ";
		fOutput << setw(9) << left << "kInf:";
		fOutput << scientific	<< setprecision(6) << right << high_pressure[1]/pow(1.e3, lambdaTotal-3.)*pow(frequency_factor, lambdaTotal-3.)	<< "\t"; 
		fOutput << setw(8) << setprecision(2) << fixed << right << high_pressure[2];
		fOutput << setw(14) << setprecision(2) << fixed << right << high_pressure[3]*energy_factor << endl;
	}
	else	// Conventional reactions
	{
		fOutput << setw(9)  << " ";
		fOutput << setw(9)  << left << "k:";
		fOutput << scientific	<< setprecision(6) << right << A/conversion_frequency_factor*pow(frequency_factor, lambdaTotal-1.) << "\t";
		fOutput << setw(8) << setprecision(2) << fixed << right << Beta;
		fOutput << setw(14) << setprecision(2) << fixed << right << E*energy_factor << endl;
	}

	if (iGlobal == true)
	{
		fOutput << setw(9) << " "; fOutput << "Global reaction" << endl;
		for (int j=1;j<=indexGlobalDirect.Size();j++)
			fOutput << setw(9) << " " << setw(20) << left << ptKinetics->species_list[indexGlobalDirect[j]] << "\t" << fixed << setw(8) << right << setprecision(3) << lambdaGlobalDirect[j] << endl;
		if (iReversible == true)
		{
			if (iThermodynamicConsistency == false)
			{ fOutput << setw(9) << " "; fOutput << "Inverse reaction without thermodynamic consistency" << endl; }
			if (iThermodynamicConsistency == true)
			{ fOutput << setw(9) << " "; fOutput << "Inverse reaction with thermodynamic consistency" << endl; }
			for (int j=1;j<=indexGlobalInverse.Size();j++)
				fOutput << setw(9) << " " << setw(20) << left << ptKinetics->species_list[indexGlobalInverse[j]] << "\t" << fixed << setw(8) << right << setprecision(3) << lambdaGlobalInverse[j] << endl;	
		}
	}

	// Fall-Off Reactions - CABR
	if(iLowPressure==true || iHighPressure==true)
	{
		if (kindPressureDependence == PRESSURE_TROE)
		{

			fOutput << setw(9) << " "; fOutput << "Troe Parameters" << endl;
			fOutput << setw(9) << " "; fOutput << scientific << "a\t" << troe[1] << endl;
			fOutput << setw(9) << " "; fOutput << scientific << "b\t" << troe[2] << endl;
			fOutput << setw(9) << " "; fOutput << scientific << "c\t" << troe[3] << endl;
			
			if (troe.Size()==3)	
			{	fOutput << setw(9) << " "; fOutput << scientific << "d\t" << 0. << endl;}

			if (troe.Size()==4)	
			{	fOutput << setw(9) << " "; fOutput << scientific << "d\t" << troe[4] << endl;}
		}
		
		if (kindPressureDependence == PRESSURE_SRI)
		{
			fOutput << setw(9) << " "; fOutput << "SRI Parameters" << endl;
			fOutput << setw(9) << " "; fOutput << scientific << "a\t" << sri[1] << endl;
			fOutput << setw(9) << " "; fOutput << scientific << "b\t" << sri[2] << endl;
			fOutput << setw(9) << " "; fOutput << scientific << "c\t" << sri[3] << endl;
			
			if (sri.Size()==3)	
			{
				fOutput << setw(9) << " "; fOutput << scientific << "d\t" << 1. << endl;
				fOutput << setw(9) << " "; fOutput << scientific << "e\t" << 0. << endl;
			}

			if (sri.Size()==5)	
			{
				fOutput << setw(9) << " "; fOutput << scientific << "d\t" << sri[4] << endl;
				fOutput << setw(9) << " "; fOutput << scientific << "e\t" << sri[5] << endl;
			}
		}
	}

	// Landau-Teller
	if (iLandauTeller == true)
	{
		fOutput << setw(9) << " "; fOutput << "Landau-Teller Parameters" << endl;
		fOutput << setw(9) << " "; fOutput << "B\t" << scientific << landau_teller[1] << endl;
		fOutput << setw(9) << " "; fOutput << "C\t" << scientific << landau_teller[2] << endl;
	}

	// Janev-Langer (TODO units)
	if (iJanevLanger == true)
	{	
		fOutput << setw(9) << " "; 
		fOutput << "Janev-Langer Parameters" << endl;
		for (int j=1;j<=9;j++)
		{	fOutput << setw(9) << " "; fOutput << "b" << j << "\t" << scientific << janev_langer[j] << endl;}
	}

	// Power-Series
	if (iPowerSeries == true)
	{	
		fOutput << setw(9) << " "; 
		fOutput << "Power-Series Parameters" << endl;
		for (int j=1;j<=4;j++)
		{	fOutput << setw(9) << " "; fOutput << "b" << j << "\t" << scientific << power_series[j] << endl;}
	}

	// Power-Series
	if (iTAR == true)
	{	
		fOutput << setw(9) << " "; 
		fOutput << "TAR Reaction Parameters" << endl;
		for (int j=1;j<=6;j++)
		{	fOutput << setw(9) << " "; fOutput << "b" << j << "\t" << scientific << TAR_series[j] << endl;}
	}

	// Chebishev
	if (iChebishevPolynomials == true)
	{	
		fOutput << setw(9) << " "; 
		fOutput << "Chebichev Polynomial Expansion" << endl;
		int m=3;
		for (int j=1;j<=int(chebishev_vector[1]);j++)
		{
			for (int k=1;k<=int(chebishev_vector[2]);k++)
			{	fOutput << setw(9) << " "; fOutput << "a(" << j << "," << k <<")\t" << chebishev_vector[m++] << endl;}
		}

		// The units are always [cm, mol] (so, no additional conversion factor is required)
		fOutput << setw(9) << " "; fOutput << "Conv.\t" << conversion_frequency_factor	<< endl;	// TODO (units)
		fOutput << setw(9) << " "; fOutput << "Tmin:\t" << chebishev_temperature[1]	<< endl;
		fOutput << setw(9) << " "; fOutput << "Tmax:\t" << chebishev_temperature[2]	<< endl;
		fOutput << setw(9) << " "; fOutput << "Pmin:\t" << chebishev_pressure[1]		<< endl;
		fOutput << setw(9) << " "; fOutput << "Pmax:\t" << chebishev_pressure[2]		<< endl;		
	}

	// Logarithmic-Pressure Dependence
	if (iPressureLogarithmic == true)
	{	
		fOutput << setw(9) << " "; 
		fOutput << "Pressure Logarithmic Interpolation" << endl;

		for (int j=1;j<=pressure_logarithmic.Rows();j++)
		{
			fOutput << setw(9) << " "; 
			fOutput << j << "\t";
			fOutput << pressure_logarithmic[j][1] << "\t";
			fOutput << pressure_logarithmic[j][2]/conversion_frequency_factor*pow(frequency_factor, lambdaTotal-1.) << "\t";
			fOutput << pressure_logarithmic[j][3] << "\t";
			fOutput << pressure_logarithmic[j][4]*energy_factor << "\t";
			fOutput << endl;
		}
	}

	if (iCollisionEfficiency == true)
	{
		fOutput << setw(9) << " "; 
		fOutput << "Bimolecular Collision Efficiency Reaction" << endl;
		fOutput << setw(11) << " "; fOutput << "  Collision diameter:       " << setw(7) << dAB*1e10	<< " [A]" << endl;
		fOutput << setw(11) << " "; fOutput << "  Reduced Molecular weight: " << setw(7) << WAB			<< " [kg/kmol]" << endl;
		fOutput << setw(11) << " "; fOutput << "  Kinetic constant:         " << scientific << kStar_collision_frequency << " * T^0.50   [m3/kmol/s]" << endl;
	}

	// Third Body Efficiencies
	if (efficiencies_indices.Size()!=0)
		for(i=1;i<=efficiencies_indices.Size();i++)
		{
			fOutput << setw(9) << " ";
			fOutput << setw(20) << left << ptKinetics->species_list[efficiencies_indices[i]] << "\tenhanced by\t" << fixed << setw(8) << right << setprecision(3) << efficiencies_coefficients[i] <<endl;
		}

	// Duplicate
	if (iDuplicate == true)
	{
		int j;
		fOutput << setw(9) << " ";
		fOutput << "Duplicate: ";
		BzzVectorInt aux;
		for(j=1;j<=ptKinetics->aDuplicate.Size();j++)
			if(ptKinetics->aDuplicate[j] == iReaction) aux.Append(ptKinetics->bDuplicate[j]);
		for(j=1;j<=ptKinetics->bDuplicate.Size();j++)
			if(ptKinetics->bDuplicate[j] == iReaction) aux.Append(ptKinetics->aDuplicate[j]);
		Sort(&aux);
		for(j=1;j<=aux.Size();j++)
			fOutput << setw(6) << right << aux[j];
		fOutput << endl;
	}

	// Reverse Rate
	if (iReverseRate != 0)
	{
		fOutput << setw(9) << " ";
		fOutput << "This reaction has the explicit reverse reaction " << iReverseReaction << endl;
	}

	// Reverse Rate
	if (iReactionBis != 0)
	{
		fOutput << setw(9) << " ";
		fOutput << "This is the reverse reaction of reaction " << iReaction << endl;
	}

	fOutput << endl << endl;
}

bool OpenSMOKE_CHEMKINInterpreter_ReactionData::IsReversibleReaction()
{
	return iReversible; 
}

bool OpenSMOKE_CHEMKINInterpreter_ReactionData::IsThirdBodyReaction()
{
	return iThirdBody; 
}

bool OpenSMOKE_CHEMKINInterpreter_ReactionData::IsLandauTellerReaction()
{
	return iLandauTeller; 
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::PrintOnBinaryFile(BzzSave &outputFile, OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_CHEMKINInterpreter_TransportData &transport)
{
	if(iReversible == true)	outputFile << 1;
	else					outputFile << 0;

		 if (iThirdBody	== false && iLandauTeller == true)		    outputFile << 100;
	else if (iThirdBody	== true  && iLandauTeller == true)		    outputFile << 10;
	else if (iThirdBody == false && iJanevLanger == true)			outputFile << 110;
	else if (iThirdBody == true  && iJanevLanger == true)			outputFile << 11;
	else if (iThirdBody == false && iPowerSeries == true)			outputFile << 120;
	else if (iThirdBody == true  && iPowerSeries == true)			outputFile << 12;
	else if (iThirdBody == false && iChebishevPolynomials == true)	outputFile << 130;
	else if (iThirdBody == true  && iChebishevPolynomials == true)	ErrorMessage("CHEBISHEV and THREE-BODY option are mutually exclusive...");
	else if (iThirdBody == false && iPressureLogarithmic == true)	outputFile << 140;
	else if (iThirdBody == true  && iPressureLogarithmic == true)	ErrorMessage("PLOG and THREE-BODY option are mutually exclusive...");
	else if (iThirdBody == false && iCollisionEfficiency == true)	outputFile << 150;
	else if (iThirdBody == true  && iCollisionEfficiency == true)	ErrorMessage("COLLEFF and THREE-BODY option are mutually exclusive...");
	else if (iThirdBody == false && iTAR == true)					outputFile << 160;
	else if (iThirdBody == true  && iTAR == true)					ErrorMessage("REACTAR and THREE-BODY option are mutually exclusive...");

	else if(iLowPressure==true)			// Fall-Off Reactions (2,3,4)
	{
		if (kindPressureDependence == PRESSURE_LINDEMANN)	outputFile << 2;
		if (kindPressureDependence == PRESSURE_TROE)		outputFile << 3;
		if (kindPressureDependence == PRESSURE_SRI)			outputFile << 4;
	}
	else if(iHighPressure==true)	// Chemically-Activated Bimolecular Reactions (5,6,7)
	{
		if (kindPressureDependence == PRESSURE_LINDEMANN)	outputFile << 5;
		if (kindPressureDependence == PRESSURE_TROE)		outputFile << 6;
		if (kindPressureDependence == PRESSURE_SRI)			outputFile << 7;
	}
	
	else
	{
		if (iThirdBody == true)					outputFile << 1;
		if (iThirdBody == false)				outputFile << 0;
	}

	if (efficiencies_indices.Size()==0)			outputFile << 0;
	else										outputFile << efficiencies_indices.Size();
	
	if (efficiencies_indices.Size()!=0)	
	{
		for(int j=1;j<=efficiencies_indices.Size();j++)
			outputFile << efficiencies_indices[j] << efficiencies_coefficients[j];
	}

	if (iLowPressure == true)	
	{
		outputFile	<< low_pressure[1]/pow(1.e3, lambdaTotal-1.)*pow(frequency_factor, lambdaTotal-1.) 
					<< low_pressure[2] 
					<< low_pressure[3]*energy_factor;

		outputFile	<< high_pressure[1]/pow(1.e3, lambdaTotal-2.)*pow(frequency_factor, lambdaTotal-2.)
					<< high_pressure[2]	 
					<< high_pressure[3]*energy_factor;
	}
	else if (iHighPressure == true)	
	{
		outputFile	<< low_pressure[1]/pow(1.e3, lambdaTotal-2.)*pow(frequency_factor, lambdaTotal-2.)
					<< low_pressure[2]
					<< low_pressure[3]*energy_factor;

		outputFile	<< high_pressure[1]/pow(1.e3, lambdaTotal-3.)*pow(frequency_factor, lambdaTotal-3.)
					<< high_pressure[2]
					<< high_pressure[3]*energy_factor;
	}
	else	// Conventional reactions
	{
		outputFile	<< A/conversion_frequency_factor*pow(frequency_factor, lambdaTotal-1.) 
					<< Beta
					<< E*energy_factor;
	}

	// Fall-Off Reactions - CABR
	if(iLowPressure==true || iHighPressure==true)
	{
		if (kindPressureDependence == PRESSURE_TROE)
		{
			outputFile << troe[1] << troe[2] << troe[3]; 
			if (troe.Size()==3)	outputFile << 0. << 0.;
			if (troe.Size()==4)	outputFile << troe[4] << 0.;
		}
		
		if (kindPressureDependence == PRESSURE_SRI)
		{
			outputFile << sri[1] << sri[2] << sri[3]; 
			if (sri.Size()==3)	outputFile << 1.0 << 0.;
			if (sri.Size()==5)	outputFile << sri[4] << sri[5];
		}

		if (kindPressureDependence == PRESSURE_LINDEMANN)
		{
			outputFile << 0. << 0. << 0. << 0. << 0.;
		}
	}

	// Landau-Teller
	if (iLandauTeller == true)
		outputFile << landau_teller[1] << landau_teller[2];

	// Janev-Langer (TODO units)
	if (iJanevLanger == true)
		for (int j=1;j<=9;j++)
			outputFile << janev_langer[j];

	// Power-Series
	if (iPowerSeries == true)
		for (int j=1;j<=4;j++)
			outputFile << power_series[j];

	// Power-Series
	if (iTAR == true)
		for (int j=1;j<=6;j++)
			outputFile << TAR_series[j];

	// Chebishev
	if (iChebishevPolynomials == true)
	{	
		outputFile << int(chebishev_vector[1]);
		outputFile << int(chebishev_vector[2]);
		int m=3;
		for (int j=1;j<=int(chebishev_vector[1]);j++)
		{
			for (int k=1;k<=int(chebishev_vector[2]);k++)
				outputFile << 	chebishev_vector[m++];
			outputFile;
		}
		// The units are always [cm, mol] (so, no additional conversion factor is required)
		outputFile << conversion_frequency_factor;	// TODO (units)
		outputFile << chebishev_temperature[1] << chebishev_temperature[2];
		outputFile << chebishev_pressure[1] << chebishev_pressure[2];
	}

	// Logarithmic-Pressure Dependence
	if (iPressureLogarithmic == true)
	{	
		int j;
		outputFile << 	pressure_logarithmic.Rows();
		for (j=1;j<=pressure_logarithmic.Rows();j++)
			outputFile << pressure_logarithmic[j][1];
		for (j=1;j<=pressure_logarithmic.Rows();j++)
			outputFile << pressure_logarithmic[j][2]/conversion_frequency_factor*pow(frequency_factor, lambdaTotal-1.);
		for (j=1;j<=pressure_logarithmic.Rows();j++)
			outputFile << pressure_logarithmic[j][3];
		for (j=1;j<=pressure_logarithmic.Rows();j++)
			outputFile << pressure_logarithmic[j][4]*energy_factor;
	}

	// Efficiency of Collision Frequency
	if (iCollisionEfficiency == true)
	{
		if (transport.IsActivated() == true)
		{
			if (indexDirect.Size()>2)	ErrorMessage("COLLEFF option can be applied only to bimolecular reactions!");
			if (indexDirect.Size()==1)
				if (nuDirect[1] != 2.)	ErrorMessage("COLLEFF option can be applied only to bimolecular reactions!");
			if (indexDirect.Size()==2)
				if (nuDirect[1] != 1. || nuDirect[2] != 1.)	ErrorMessage("COLLEFF option can be applied only to bimolecular reactions!");

			if (indexDirect.Size()==1)
			{
				dAB = transport.ReducedDiameter(indexDirect[1], indexDirect[1])*1e-10;
				WAB = thermo.ReducedMolecularWeight(nameDirect[1], nameDirect[1]);
			}
			if (indexDirect.Size()==2)
			{
				dAB = transport.ReducedDiameter(indexDirect[1], indexDirect[2])*1.e-10;
				WAB = thermo.ReducedMolecularWeight(nameDirect[1], nameDirect[2]);
			}

			kStar_collision_frequency = Constants::Nav_kmol*(dAB*dAB)*sqrt(8.*Constants::pi*Constants::R_J_kmol/WAB);
			outputFile << kStar_collision_frequency;
		}
		else
			ErrorMessage("COLLEFF option can be used only if transport properties are available...");
	}
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::CompactStoichiometry()
{
	CompactStoichiometries(indexDirect, nuDirect, nReactants, nameDirect);
	CompactStoichiometries(indexInverse, nuInverse, nProducts, nameInverse);
}

void OpenSMOKE_CHEMKINInterpreter_ReactionData::WriteWarningMessageOnFile(const std::string message)
{
	stringstream string_iReaction;	string_iReaction << iReaction;
	stringstream string_iLine;		string_iLine << iLine;

	std::string warning_message;
	warning_message += "Reaction: " + string_iReaction.str() + "\n";
	warning_message += "Line:     " + string_iLine.str() + "\n";
	warning_message += message;

	ptKinetics->fWarning->WriteWarningMessage(warning_message);
}

void CompactStoichiometries(BzzVectorInt &index, BzzVector &nu, int &n, vector<string> &names)
{
	int i;
	
	for(i=1;i<=index.Size();i++)
	{
		if (index[i] != 0)
		for(int j=1;j<=index.Size();j++)
			if (index[i] == index[j])
				if (i != j)
				{
					index[j]   = 0;
					nu[i]	  += nu[j];
					nu[j]	   = 0;
				}	
	}

	BzzVectorInt indices;
	for(i=1;i<=index.Size();i++)
		if (index[i] == 0)	indices.Append(i);

	if (indices.Size() != 0)
	{
		vector<string> aux_names; aux_names.push_back("names");
		for(i=1;i<=index.Size();i++)
			if (index[i] != 0)	aux_names.push_back(names[i]);

		names.clear();
		names = aux_names;

		index.DeleteElements(indices);
		nu.DeleteElements(indices);

		n = index.Size();
	}
}

void CompleteGlobalStoichiometries(BzzVectorInt &index, BzzVectorInt &indexGlobal, BzzVector &nu, BzzVector &lambdaGlobal)
{
		for(int i=1;i<=index.Size();i++)
		{
			bool iFound = false;
			for(int j=1;j<=indexGlobal.Size();j++)
				if (index[i] == indexGlobal[j])
				{
					iFound = true;
					break;
				}
			
			if (iFound == false)	
			{
				indexGlobal.Append(index[i]);
				lambdaGlobal.Append(nu[i]);
			}
		}
}
