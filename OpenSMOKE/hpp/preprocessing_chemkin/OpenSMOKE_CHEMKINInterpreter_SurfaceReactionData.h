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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_SURFACEREACTIONDATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_SURFACEREACTIONDATA_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"

class OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData;
class OpenSMOKE_CHEMKINInterpreter_ThermoData;
class OpenSMOKE_CHEMKINInterpreter_UnitsData;
class OpenSMOKE_CHEMKINInterpreter_TransportData;
class OpenSMOKE_PreProcessorSurfaceMaterial;

class OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData
{
public:

	OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData();
	void Setup(const int _iReaction, const int _iLine, OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData *_ptKinetics, OpenSMOKE_CHEMKINInterpreter_UnitsData  *_ptUnits);
	void Setup(const double _energy_factor, const double _frequency_factor);
	void ReverseRateFromReaction(const int _iReactionBis, OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData &r);

	virtual ~OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData();

	string reaction_string_clean;
	string reaction_string_complete;

	double A;
	double Beta;
	double E;

	double energy_factor;
	double frequency_factor;

	int nReactants;
	vector<string>	nameDirect;
	BzzVectorInt	indexDirect;
	BzzVector		nuDirect;
	vector<char>	phaseDirect;

	int nProducts;
	vector<string>	nameInverse;
	BzzVectorInt	indexInverse;
	BzzVector		nuInverse;
	vector<char>	phaseInverse;

	int nGlobalReactants;
	vector<string>	nameGlobalDirect;
	BzzVectorInt	indexGlobalDirect;
	BzzVector		lambdaGlobalDirect;

	int nGlobalProducts;
	vector<string>	nameGlobalInverse;
	BzzVectorInt	indexGlobalInverse;
	BzzVector		lambdaGlobalInverse;

	bool iReversible;
	bool iGlobal;

	BzzVectorInt	indexDirectOrdered;
	BzzVectorInt	indexInverseOrdered;
	BzzVectorInt	indexDirectOrderedSequence;
	BzzVectorInt	indexInverseOrderedSequence;

	void Check();
	void CheckStoichiometry(OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, BzzVectorInt &gas_indices, BzzVectorInt &site_indices, BzzVectorInt &bulk_indices, vector<string> &element_names, vector<string> &isotope_names);
	void CompactStoichiometry();
	void ReactionString();
	void Complete();
	void PrintOnFile(ofstream &fOutput);
	void PrintOnBinaryFile(BzzSave &outputFile, BzzSave &asciiFile, OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_PreProcessorSurfaceMaterial &material);
	void SummaryOnFile(ofstream &fOutput);
	void WriteWarningMessageOnFile(const string message);

	OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData	*ptKinetics;
	OpenSMOKE_CHEMKINInterpreter_UnitsData				*ptUnits;

	int iLine;
	int iReaction;
	int iReactionBis;
	int iReverseReaction;
	bool iReverseRate;
	bool iDuplicate;

	bool iConventionalReaction;
	bool iStickyReaction;
	bool iCoverageDependentReaction;
	bool iLangmuirHinshelwoodReaction;

	double sumNu;
	double sumNuGas;
	double sumNuSite;
	double sumNuBulk;
	double lambdaTotal;
	double conversion_frequency_factor;

	void StickyReaction(vector<string> instructions, OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_PreProcessorSurfaceMaterial &material, BzzVectorInt &gas_indices, BzzVectorInt &site_indices);
	void CoverageDependentReactions(vector<string> instructions);
	void LangmuirHinshelwoodReactions(vector<string> instructions);	
	void LangmuirHinshelwoodDenominatorExponentParameter(vector<string> instructions);
	void LangmuirHinshelwoodEquilibriumPressure(vector<string> instructions);
	
	int nLangmuirHinshelwoodReaction;
	double langmuirHinshelwoodDenominatorExponentParameter;
	string langmuirHinshelwoodEquilibriumPressure;
	BzzVectorInt langmuirHinshelwoodReactionSpecies;
	BzzVector langmuirHinshelwoodReactionParameters;
	BzzVectorInt langmuirHinshelwoodReactantSpecies;
	BzzVector    langmuirHinshelwoodReactantExponent;

	int nCoverageDependentReaction;
	BzzVectorInt coverageDependentReactionSpecies;
	BzzVector coverageDependentReactionParameters;

	double powerSiteDensity;
	double powerOccupancySites;

	double stickyGasSpeciesMW;
	double stickyPowerOccupancySites;

private:
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_REACTIONDATA_H)
