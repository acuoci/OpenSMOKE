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
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ThermoData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ThermoSpecies.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_UnitsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter.h"

void CompactStoichiometries(BzzVectorInt &index, vector<char> &phase, BzzVector &nu, int &n, vector<string> &names);
void SurfaceCompleteGlobalStoichiometries(BzzVectorInt &index, BzzVectorInt &indexGlobal, BzzVector &nu, BzzVector &lambdaGlobal);

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:    OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData"	<< endl;
    cout << "Error:    " << message         << endl;
    cout << "Reaction: " << iReaction       << endl;
    cout << "Line:     " << iLine           << endl;    
	cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData"	<< endl;
    cout << "Warning:  " << message         << endl;
    cout << "Reaction: " << iReaction       << endl;
    cout << "Line:     " << iLine           << endl;    
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData()
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
	sumNuGas		= 0.;
	sumNuSite		= 0.;
	sumNuBulk		= 0.;
	lambdaTotal		= 0.;

	iConventionalReaction			= true;
	iStickyReaction					= false;
	iLangmuirHinshelwoodReaction	= false;
	iCoverageDependentReaction		= false;

	iReversible		= false;
	iGlobal			= false;
	iDuplicate		= false;

	iReverseRate	= false;
	iReactionBis	= 0;

	nCoverageDependentReaction		= 0;

	nLangmuirHinshelwoodReaction	= 0;
	langmuirHinshelwoodDenominatorExponentParameter = 2.;
	langmuirHinshelwoodEquilibriumPressure = "none";

	stickyGasSpeciesMW = 0;
	stickyPowerOccupancySites = 0;
}

OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::~OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData()
{
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::Setup(const int _iReaction, const int _iLine, OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData *_ptKinetics, OpenSMOKE_CHEMKINInterpreter_UnitsData  *_ptUnits)
{
	iReaction	= _iReaction;
	iLine		= _iLine;
	ptKinetics	= _ptKinetics;
	ptUnits		= _ptUnits;
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::Setup(const double _energy_factor, const double _frequency_factor)
{
	energy_factor = _energy_factor;
	frequency_factor = _frequency_factor;
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::CompactStoichiometry()
{
	CompactStoichiometries(indexDirect, phaseDirect, nuDirect, nReactants, nameDirect);
	CompactStoichiometries(indexInverse, phaseInverse, nuInverse, nProducts, nameInverse);
}

void CompactStoichiometries(BzzVectorInt &index, vector<char> &phase, BzzVector &nu, int &n, vector<string> &names)
{
	int i;
	
	for(i=1;i<=index.Size();i++)
	{
		if (index[i] != 0)
		for(int j=1;j<=index.Size();j++)
			if (index[i] == index[j] && phase[i]==phase[j])
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

		vector<char> aux_phases; aux_phases.push_back('0');
		for(i=1;i<=index.Size();i++)
			if (index[i] != 0)	aux_phases.push_back(phase[i]);

		names.clear();
		names = aux_names;

		phase.clear();
		phase = aux_phases;

		index.DeleteElements(indices);
		nu.DeleteElements(indices);

		n = index.Size();
	}
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::Complete()
{
	int i;

	if (iGlobal == true)
		SurfaceCompleteGlobalStoichiometries(indexDirect, indexGlobalDirect, nuDirect, lambdaGlobalDirect);

	// 
	if (iGlobal == false)	for(i=1;i<=nReactants;i++)					lambdaTotal += nuDirect[i];
	else					for(i=1;i<=indexGlobalDirect.Size();i++)	lambdaTotal += lambdaGlobalDirect[i];
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::Check()
{
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
		ErrorMessage("iGlobal == true && iReversible == true");
		/*
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
		*/
	}
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::SummaryOnFile(ofstream &fOutput)
{
	if (iReactionBis != 0)
		fOutput << setw(7) << right << iReactionBis;
	else
		fOutput << setw(7) << right << iReaction;
	fOutput << ". ";
	fOutput << reaction_string_clean << endl;

	// Conventional reactions
	{
		fOutput << setw(9)  << " ";
		fOutput << setw(9)  << left << "k:";
		fOutput << scientific	<< setprecision(6) << right << A << "\t";
		fOutput << setw(8)      << setprecision(2) << fixed << right << Beta;
		fOutput << setw(14)     << setprecision(2) << fixed << right << E*energy_factor << endl;
	}

	fOutput << endl << endl;
}

void SurfaceCompleteGlobalStoichiometries(BzzVectorInt &index, BzzVectorInt &indexGlobal, BzzVector &nu, BzzVector &lambdaGlobal)
{
	{
		cout << "SurfaceCompleteGlobalStoichiometries not yet available" << endl;
		cout << "Press enter to exit..." << endl;
		getchar();
		exit(-1);
	}
	/*
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
	*/
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::ReactionString()
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

	sumNu = 0.;
	for(i=1;i<=nProducts;i++)	sumNu += nuInverse[i];
	for(i=1;i<=nReactants;i++)	sumNu -= nuDirect[i];

	reaction_string_complete = reaction_string_clean;
	
	if(iReactionBis != 0)	
	{
		stringstream index;
		index << iReaction;
		reaction_string_complete  += "#";
		reaction_string_complete += index.str();
		reaction_string_complete += "bis";
		reaction_string_complete += "#";
	}
	
	{
		string dummy;
		stringstream index;
		if(iReactionBis != 0)	index << iReactionBis;
		else					index << iReaction;

		dummy  = "#";
		dummy += index.str();
		dummy += "#";

		reaction_string_complete.insert(0,dummy);
	}

}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::ReverseRateFromReaction(const int _iReactionBis, OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData &r)
{
	ErrorMessage("ReverseRateFromReaction not yet");
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::CheckStoichiometry(OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, BzzVectorInt &gas_indices, BzzVectorInt &site_indices, BzzVectorInt &bulk_indices, vector<string> &element_names, vector<string> &isotope_names)
{
	int nElements = thermo.species[gas_indices[1]].element_indices.Size();
	int nIsotopes = thermo.species[gas_indices[1]].isotope_indices.Size();

	// Elements
	{
		int i;	
		BzzVector sumReactants(nElements);
		BzzVector sumProducts(nElements);

		for(i=1;i<=nReactants;i++)
		{
			char phase = phaseDirect[i];
			
			int h=0;
			if (phase == 'G')		h = gas_indices[indexDirect[i]];
			else if (phase == 'S')	h = site_indices[indexDirect[i]];
			else if (phase == 'B')	h = bulk_indices[indexDirect[i]];
			
			for(int j=1;j<=nElements;j++)
				sumReactants[j] += thermo.species[h].element_indices[j] * nuDirect[i];
		}

		for(i=1;i<=nProducts;i++)
		{	
			char phase = phaseInverse[i];

			int h =0;
			if (phase == 'G')      h= gas_indices[indexInverse[i]];
			else if (phase == 'S') h= site_indices[indexInverse[i]];
			else if (phase == 'B') h= bulk_indices[indexInverse[i]];
			
			for(int j=1;j<=nElements;j++)
				sumProducts[j] += thermo.species[h].element_indices[j] * nuInverse[i];
		}

		for(int j=1;j<=nElements;j++)
			if (sumReactants[j] > sumProducts[j]*1.0001 || sumReactants[j] < sumProducts[j]*0.9999)
			{
				stringstream sumR; sumR << sumReactants[j];
				stringstream sumP; sumP << sumProducts[j];
				string message = "Umbalanced stoichiometry on " + element_names[j] + " element\n";
				message += "  Reactants: " + sumR.str() + " - Products: " + sumP.str();
				ErrorMessage(message);
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
			char phase = phaseDirect[i];

			int h=0;
			if (phase == 'G')      h = gas_indices[indexDirect[i]];
			else if (phase == 'S') h = site_indices[indexDirect[i]];
			else if (phase == 'B') h = bulk_indices[indexDirect[i]];

			for(int j=1;j<=nIsotopes;j++)
				sumReactants[j] += thermo.species[h].isotope_indices[j] * nuDirect[i];
		}

		for(i=1;i<=nProducts;i++)
		{
			char phase = phaseInverse[i];

			int h=0;
			if (phase == 'G')      h = gas_indices[indexInverse[i]];
			else if (phase == 'S') h = site_indices[indexInverse[i]];
			else if (phase == 'B') h = bulk_indices[indexInverse[i]];

			for(int j=1;j<=nIsotopes;j++)
				sumProducts[j] += thermo.species[h].isotope_indices[j] * nuInverse[i];
		}

		for(int j=1;j<=nIsotopes;j++)
			if (sumReactants[j] != sumProducts[j])
				ErrorMessage("Umbalanced stoichiometry on " + isotope_names[j] + " element");	
	}
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::PrintOnBinaryFile(BzzSave &outputFile, BzzSave &asciiFile, OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_PreProcessorSurfaceMaterial &material)
{
	if(iReversible == true)	outputFile << 1;
	else					outputFile << 0;

		 if (iStickyReaction				== true)	outputFile << 1;
	else if (iLangmuirHinshelwoodReaction	== true)	outputFile << 2;
	else if (iCoverageDependentReaction		== true)	outputFile << 3;
	else												outputFile << 0;

	if(iReversible == true)	asciiFile << 1 ;
	else					asciiFile << 0 ;

		 if (iStickyReaction				== true)	asciiFile << 1 ;
	else if (iLangmuirHinshelwoodReaction	== true)	asciiFile << 2 ;
	else if (iCoverageDependentReaction		== true)	asciiFile << 3 ;
	else												asciiFile << 0 ;
		
	// Stoichiometry
	int nGasReactants  = 0;
	int nSiteReactants = 0;
	int nBulkReactants = 0;

	for(int j=1;j<=nReactants;j++)
	{
		if (phaseDirect[j] == 'G')	{ nGasReactants++;	sumNuGas  -= nuDirect[j]; }
		if (phaseDirect[j] == 'S')	{ nSiteReactants++;	sumNuSite -= nuDirect[j]; }
		if (phaseDirect[j] == 'B')	{ nBulkReactants++;	sumNuBulk -= nuDirect[j]; }
	}

	int nGasProducts  = 0;
	int nSiteProducts = 0;
	int nBulkProducts = 0;

	for(int j=1;j<=nProducts;j++)
	{
		if (phaseInverse[j] == 'G')	{ nGasProducts++;	sumNuGas  += nuInverse[j]; }
		if (phaseInverse[j] == 'S')	{ nSiteProducts++;	sumNuSite += nuInverse[j]; }
		if (phaseInverse[j] == 'B')	{ nBulkProducts++;	sumNuBulk += nuInverse[j]; }
	}

	outputFile << nReactants;
	outputFile << nProducts;
	outputFile << nGasReactants;
	outputFile << nSiteReactants;
	outputFile << nBulkReactants;
	outputFile << nGasProducts;
	outputFile << nSiteProducts;
	outputFile << nBulkProducts;
	outputFile << sumNuGas;
	outputFile << sumNuSite;
	outputFile << sumNuBulk;

	asciiFile << nReactants;
	asciiFile << nProducts;
	asciiFile << nGasReactants;
	asciiFile << nSiteReactants;
	asciiFile << nBulkReactants;
	asciiFile << nGasProducts;
	asciiFile << nSiteProducts;
	asciiFile << nBulkProducts;
	asciiFile << sumNuGas;
	asciiFile << sumNuSite;
	asciiFile << sumNuBulk;

	BzzVector lambdaDirect  = nuDirect;
	BzzVector lambdaInverse = nuInverse;

	for(int j=1;j<=nReactants;j++)
		if (phaseDirect[j] == 'G')	outputFile << indexDirect[j] << nuDirect[j] << lambdaDirect[j];

	for(int j=1;j<=nReactants;j++)
		if (phaseDirect[j] == 'S')	outputFile << indexDirect[j] << nuDirect[j] << lambdaDirect[j];

	for(int j=1;j<=nReactants;j++)
		if (phaseDirect[j] == 'B')	outputFile << indexDirect[j] << nuDirect[j] << lambdaDirect[j];


	for(int j=1;j<=nProducts;j++)
		if (phaseInverse[j] == 'G')	outputFile << indexInverse[j] << nuInverse[j] << lambdaInverse[j];

	for(int j=1;j<=nProducts;j++)
		if (phaseInverse[j] == 'S')	outputFile << indexInverse[j] << nuInverse[j] << lambdaInverse[j];

	for(int j=1;j<=nProducts;j++)
		if (phaseInverse[j] == 'B')	outputFile << indexInverse[j] << nuInverse[j] << lambdaInverse[j];


	for(int j=1;j<=nReactants;j++)
		if (phaseDirect[j] == 'G')	asciiFile << indexDirect[j] << nuDirect[j] <<   lambdaDirect[j];

	for(int j=1;j<=nReactants;j++)
		if (phaseDirect[j] == 'S')	asciiFile << indexDirect[j] << nuDirect[j] <<  lambdaDirect[j];

	for(int j=1;j<=nReactants;j++)
		if (phaseDirect[j] == 'B')	asciiFile << indexDirect[j] << nuDirect[j] << lambdaDirect[j];


	for(int j=1;j<=nProducts;j++)
		if (phaseInverse[j] == 'G')	asciiFile << indexInverse[j] << nuInverse[j] << lambdaInverse[j];

	for(int j=1;j<=nProducts;j++)
		if (phaseInverse[j] == 'S')	asciiFile << indexInverse[j] << nuInverse[j] <<  lambdaInverse[j];

	for(int j=1;j<=nProducts;j++)
		if (phaseInverse[j] == 'B')	asciiFile << indexInverse[j] << nuInverse[j] <<  lambdaInverse[j];

	// Writing reaction names on file
	char name[Constants::REACTION_NAME_SIZE];
	strcpy(name, reaction_string_clean.c_str());
	outputFile.fileSave.write((char*) name, sizeof(name));
	asciiFile << reaction_string_clean;

	// Conversion factors
	{
		double alfa = 0.;
		double beta = 0.;
		
		for(int j=1;j<=nReactants;j++)
		{
			if (phaseDirect[j] == 'G')		beta += lambdaDirect[j];
			else if (phaseDirect[j] == 'S')	alfa += lambdaDirect[j];
		}

		conversion_frequency_factor = pow(10., 1.-alfa-3.*beta);

		if (iStickyReaction == true)	conversion_frequency_factor = 1.;
	}

	// Kinetic data
	outputFile	<< A*conversion_frequency_factor
				<< Beta
				<< E*energy_factor;

	asciiFile	<< A*conversion_frequency_factor 
				<< Beta 
				<< E*energy_factor;

	if (iReversible == true)
	{
		cout << "Rev" << endl;
		powerSiteDensity = pow(10. * material.site_density[1], sumNuSite);
		cout << sumNuSite << " " << powerSiteDensity << endl;
		powerOccupancySites = 1.;
		for(int j=1;j<=nReactants;j++)
		{
			if (phaseDirect[j] == 'S')	
				powerOccupancySites *= pow(material.site_occupancy[1][indexDirect[j]], nuDirect[j]);
		}
		for(int j=1;j<=nProducts;j++)
		{
			if (phaseInverse[j] == 'S')	
				powerOccupancySites *= pow(material.site_occupancy[1][indexInverse[j]], -nuInverse[j]);
		}
		cout << powerOccupancySites << endl;

		outputFile << powerSiteDensity << powerOccupancySites;
		asciiFile << powerSiteDensity << " " << powerOccupancySites;
		ErrorMessage("No equilibrium reactions");
		getchar();
	}

	if (iStickyReaction == true)
	{
		outputFile << stickyGasSpeciesMW << stickyPowerOccupancySites;
		asciiFile  << stickyGasSpeciesMW << " " << stickyPowerOccupancySites;
	}

	if (iCoverageDependentReaction == true)
	{
		outputFile << coverageDependentReactionSpecies << coverageDependentReactionParameters;
		asciiFile  << coverageDependentReactionSpecies << coverageDependentReactionParameters;
	}

	if (iLangmuirHinshelwoodReaction == true)
	{
		outputFile	<< langmuirHinshelwoodReactionSpecies;
		outputFile	<< langmuirHinshelwoodReactionParameters;
		outputFile	<< langmuirHinshelwoodReactantSpecies;
		outputFile	<< langmuirHinshelwoodReactantExponent;
		outputFile	<< langmuirHinshelwoodDenominatorExponentParameter;
		outputFile	<< langmuirHinshelwoodEquilibriumPressure ;

		asciiFile	<< langmuirHinshelwoodReactionSpecies;
		asciiFile	<< langmuirHinshelwoodReactionParameters;
		asciiFile	<< langmuirHinshelwoodReactantSpecies;
		asciiFile	<< langmuirHinshelwoodReactantExponent;
		asciiFile	<< langmuirHinshelwoodDenominatorExponentParameter;
		asciiFile	<< langmuirHinshelwoodEquilibriumPressure;
	}

}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::StickyReaction(vector<string> instructions, OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_PreProcessorSurfaceMaterial &material, BzzVectorInt &gas_indices, BzzVectorInt &site_indices)
{
	int nInstructions = instructions.size()-1;
	string error_message = "Syntax error in Sticky Reaction (STICK) option";
	
	// Checking for syntax errors
	if (nInstructions != 1)	ErrorMessage(error_message);

	// Checking 
	bool foundGasPhase = false;
	for(int j=1;j<=nReactants;j++)
	{
		if (phaseDirect[j] == 'G')	
		{
			if (foundGasPhase == true)	ErrorMessage(error_message + " Only 1 gas species must be specified");
			if (nuDirect[j] != 1.)	ErrorMessage(error_message + " The gas specie must have a stoichiometric coefficient equal to 1");
			int h = gas_indices[indexDirect[j]];
			stickyGasSpeciesMW = thermo.species[h].mw;
			foundGasPhase = true;
		}
	}
	if (foundGasPhase == false)	
		ErrorMessage(error_message + " A gas species must be specified");

	// Sticky data
	stickyPowerOccupancySites = 1.;
	double m = 0;
	for(int j=1;j<=nReactants;j++)
	{
		if (phaseDirect[j] == 'S')	
		{			
			// TODO: lambda--nu
			double lambdaDirect = nuDirect[j];
			stickyPowerOccupancySites *= pow(material.site_occupancy[1][indexDirect[j]], lambdaDirect);
			m += nuDirect[j];
		}
	}
	// Conversion from mol/cm2 to kmol/m2
	stickyPowerOccupancySites /= pow(10.*material.site_density[1], m);

	// Assign properties
	iStickyReaction = true;
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::CoverageDependentReactions(vector<string> instructions)
{
	
	int nInstructions = instructions.size()-1;
	string error_message = "Syntax error in CoverageDependent Reaction (COV) option";
	
	// Checking for syntax errors
	if (nInstructions != 7)							ErrorMessage(error_message);
	if (instructions[2] != "/")						ErrorMessage(error_message);
	if (instructions[nInstructions] != "/")			ErrorMessage(error_message);
	
	// Checking for duplicates
	iCoverageDependentReaction = true;
	nCoverageDependentReaction++;
	
	// Read parameters
	{
		int index;
		char phase;
		ptKinetics->RecognizeSpecies(instructions[3], index, phase);
		if (phase != 'S')	ErrorMessage(error_message);
		coverageDependentReactionSpecies.Append(index);

		// Assign properties
		coverageDependentReactionParameters.Append(atof(instructions[4].c_str()));
		coverageDependentReactionParameters.Append(atof(instructions[5].c_str()));
		coverageDependentReactionParameters.Append(atof(instructions[6].c_str()));
	}
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::LangmuirHinshelwoodReactions(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	string error_message = "Syntax error in Langmuir-Hinshelwood Reaction (LANG) option";
	
	// Checking for syntax errors
	if (nInstructions != 8)							ErrorMessage(error_message);
	if (instructions[2] != "/")						ErrorMessage(error_message);
	if (instructions[nInstructions] != "/")			ErrorMessage(error_message);
	
	// Checkings
	if (iReversible == true)
		ErrorMessage("The LANG option can be used only for irreversible reactions");

	// Reactants
	ChangeDimensions(nReactants, &langmuirHinshelwoodReactantSpecies);
	ChangeDimensions(nReactants, &langmuirHinshelwoodReactantExponent);
	for(int i=1;i<=nReactants;i++)
	{
		int index;
		char phase;
		ptKinetics->RecognizeSpecies(instructions[3], index, phase);
		if (phase != 'G')	ErrorMessage(error_message);
		langmuirHinshelwoodReactantSpecies[i] = index;
		langmuirHinshelwoodReactantExponent[i] = nuDirect[i];
	}

	// Counters
	iLangmuirHinshelwoodReaction = true;
	nLangmuirHinshelwoodReaction++;
	
	// Read parameters
	{
		int index;
		char phase;
		ptKinetics->RecognizeSpecies(instructions[3], index, phase);
		if (phase != 'G')	ErrorMessage(error_message);
		langmuirHinshelwoodReactionSpecies.Append(index);

		// Assign properties
		double conversion_factor = pow(1000., -atof(instructions[7].c_str()));
		langmuirHinshelwoodReactionParameters.Append(atof(instructions[4].c_str())*conversion_factor);
		langmuirHinshelwoodReactionParameters.Append(atof(instructions[5].c_str()));
		langmuirHinshelwoodReactionParameters.Append(atof(instructions[6].c_str()));
		langmuirHinshelwoodReactionParameters.Append(atof(instructions[7].c_str()));
	}
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::LangmuirHinshelwoodDenominatorExponentParameter(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	string error_message = "Syntax error in Langmuir-Hinshelwood Denominator Exponent Parameter (LHDE) option";
	
	// Checking for syntax errors
	if (nInstructions != 4)							ErrorMessage(error_message);
	if (instructions[2] != "/")						ErrorMessage(error_message);
	if (instructions[nInstructions] != "/")			ErrorMessage(error_message);
	
	// Checkings
	if (iLangmuirHinshelwoodReaction == false)
		ErrorMessage("The LHDE option can be used only for LANG reactions");
	
	// Assign properties
	langmuirHinshelwoodDenominatorExponentParameter = atof(instructions[3].c_str());
}

void OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData::LangmuirHinshelwoodEquilibriumPressure(vector<string> instructions)
{
	int nInstructions = instructions.size()-1;
	string error_message = "Syntax error in Langmuir-Hinshelwood Equilibrium Pressure (LHPR) option";
	
	// Checking for syntax errors
	if (nInstructions != 4)							ErrorMessage(error_message);
	if (instructions[2] != "/")						ErrorMessage(error_message);
	if (instructions[nInstructions] != "/")			ErrorMessage(error_message);
	
	// Checkings
	if (iLangmuirHinshelwoodReaction == false)
		ErrorMessage("The LHPR option can be used only for LANG reactions");
	if (instructions[3] != "ATM"  && instructions[3] != "atm" &&
	    instructions[3] != "BAR"  && instructions[3] != "bar" &&
	    instructions[3] != "TORR" && instructions[3] != "torr" &&
	    instructions[3] != "PASC" && instructions[3] != "pasc" &&
	    instructions[3] != "DYNE" && instructions[3] != "dyne" )
		ErrorMessage("The pressure units in the LHPR options are not recognizable");
	
	// Assign properties
	langmuirHinshelwoodEquilibriumPressure = instructions[3];
}