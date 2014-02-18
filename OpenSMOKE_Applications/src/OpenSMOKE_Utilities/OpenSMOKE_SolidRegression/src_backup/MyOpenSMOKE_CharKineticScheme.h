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

#ifndef OPENSMOKE_CHARKINETICSCHEME
#define OPENSMOKE_CHARKINETICSCHEME

#include "BzzMath.hpp"

class OpenSMOKE_CharKineticScheme;
class MyOpenSMOKE_SolidRegression;

class OpenSMOKE_CharKineticScheme_ReactionData
{
public:

	OpenSMOKE_CharKineticScheme_ReactionData();
	void Setup(const int _iReaction, const int _iLine, OpenSMOKE_CharKineticScheme *_ptKinetics);

	void CompactStoichiometry();
	void ReactionString();
	void SummaryOnFile(ofstream &fOutput);

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
	
	int nProducts;
	vector<string>	nameInverse;
	BzzVectorInt	indexInverse;
	BzzVector		nuInverse;

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
	bool iThirdBody;
	int  iReactionBis;

	OpenSMOKE_CharKineticScheme	*ptKinetics;
	int iLine;
	int iReaction;

	double sumNu;
	double lambdaTotal;

private:
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

class OpenSMOKE_CharKineticScheme
{

friend class MyOpenSMOKE_SolidRegression;

public:

	OpenSMOKE_CharKineticScheme();
	void SetName(const string name);

	void ReadDatabase(const string file_name);
	void ReadKineticScheme(const string file_name);
	
	void UpdateKineticConstants(const double T);
	void UpdateKineticConstants(BzzVector &Ccarbon, BzzVector &Cbulk, BzzVector &Cgas, BzzVector &Cadsorbed);

protected:

	BzzVector	kappa;			// Kinetic constnat [kmol, m, s]
	BzzVector	rr;				// Reaction rates [kmol/m2/s]
	BzzVector	R;				// Formation rates [kmol/m2/s]
	double		Rchar;			// Formation rates [kmol/m2/s]

	BzzVector A;				// Frequency factor [kmol, m, s]
	BzzVector Beta;				// Temperature exponent [K]
	BzzVector E;				// Activation energy [J/kmol]

private:

	string name_kinetics_file;
	vector<string> lines;
	int number_of_lines;
	int number_of_reactions;

	BzzVectorInt indexBlankLines;
	BzzVectorInt indexCommentLines;
	BzzVectorInt indexLines;
	BzzVectorInt lineReactionIndex;
	BzzVectorInt additionalLineReactionIndex;
	
	ofstream fLog;

	int startGas;
	int startCarbon;
	int startAdsorbed;
	int startBulk;
	int startReactions;
	int startEnd;

	int endGas;
	int endCarbon;
	int endAdsorbed;
	int endBulk;
	int endReactions;
	int endEnd;

	vector<string> speciesBulk;
	vector<string> speciesGas;
	vector<string> speciesAdsorbed;
	vector<string> speciesCarbon;
	vector<string> list_of_species;

	int nGas;
	int nBulk;
	int nAdsorbed;
	int nCarbon;

	OpenSMOKE_CharKineticScheme_ReactionData *reactions;

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	void ErrorMessage(const string message, const int iLine);
	void WarningMessage(const string message, const int iLine);

private:

	BzzVector eta_adsorbed;		// Number of char surface sites occupied by an adsorbed molecule

	BzzMatrix lambdaAdsorbed;	// Reaction orders: Adsorbed Species 
	BzzMatrix lambdaCarbon;		// Reaction orders: Free sites
	BzzMatrix lambdaGas;		// Reaction orders: Gas Species
	BzzMatrix lambdaBulk;		// Reaction orders: Bulk Species

	BzzMatrix nuCarbon;			// Stoichiometric Coefficients: Free sites
	BzzMatrix nuAdsorbed;		// Stoichiometric Coefficients: Adsorbed Species
	BzzMatrix nuBulk;			// Stoichiometric Coefficients: Bulk Species
	BzzMatrix nuGas;			// Stoichiometric Coefficients: Gas Species
	BzzVector nuChar;			// Stoichiometric Coefficients: ????????????

private:
	
	vector<string> database_species;
	vector<double> database_species_eta;
	double GetEta(const string species_name);

	BzzVector moleDimension;
	BzzVector lengthDimension;
	BzzVector conversion_from_kmol_m_s__to__mol_m_s;
	BzzVector conversion_from_kmol_m_s__to__mol_cm_s;

private:

	void ParsingSections();
	void SyntaxCorrections();
	void CleaningLines();
	void ParsingSpecies();
	void ParsingReactions();
	void AddSpaces(string &line, const char symbol);
	void ParsingSpecies(const int start, const int end, vector<string> &species);
	void ParsingSpecies(vector<string> instructions, vector<string> &species);
	void ParsingReactions(const int iReaction, const int iLine, vector<string> instructions);
	int  RecognizeSpecies(const string name_species);
	void ParsingReactionsAdditionalData();
	void MemoryAllocation();
	void Summary(ofstream &fOutput);
	string DimensionString(const int i);

};


#endif // OPENSMOKE_CHARKINETICSCHEME