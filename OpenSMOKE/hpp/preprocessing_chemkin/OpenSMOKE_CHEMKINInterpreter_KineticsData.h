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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_KINETICSDATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_KINETICSDATA_H

#include "BzzMath.hpp"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ReactionData.h"

class OpenSMOKE_WarningFile;
class OpenSMOKE_CHEMKINInterpreter;
class OpenSMOKE_CHEMKINInterpreter_TransportData;

class OpenSMOKE_CHEMKINInterpreter_KineticsData
{
public:
	OpenSMOKE_CHEMKINInterpreter_KineticsData();
	virtual ~OpenSMOKE_CHEMKINInterpreter_KineticsData();

	BzzVectorInt lineReactionIndex;
	BzzVectorInt additionalLineReactionIndex;

	int number_of_reactions;
	vector<string> species_list;

	OpenSMOKE_CHEMKINInterpreter_ReactionData *reaction;
	OpenSMOKE_CHEMKINInterpreter_ReactionData *reaction_reverserate;
	void Setup(ofstream *_fLog, OpenSMOKE_WarningFile *fWarning);
	int  RecognizeSpecies(const string name_species);
	void CheckDuplicateReactions();
	void CheckReactionRatesReactions();

	bool CheckDuplicateReactions(const int i, const int j);
	bool CheckDuplicateReactionsGivenOrder(const int i, const int j);
	bool CheckDuplicateReactionsInverseOrder(const int i, const int j);

	void PrintNamesFile(const string file_name);
	void PrintStoichiometryFile(const string file_name);
	void PrintReactionsFile(const string file_name);
	void PrintBinaryFile(const string file_name, OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_CHEMKINInterpreter_TransportData &transport, OpenSMOKE_CHEMKINInterpreter &interp);

	void SummaryOnFile(ofstream &fOutput);
	void Statistics();
	void PrepareStoichiometry();
	void PrepareReactionOrders();
	void FinalChecks();

	ofstream *fLog;
	OpenSMOKE_WarningFile *fWarning;

	BzzVectorInt aDuplicate;
	BzzVectorInt bDuplicate;

private:
	
	BzzVectorInt iCABR;
	BzzVectorInt iFallOff;
	BzzVectorInt iCABR_Lindemann;
	BzzVectorInt iCABR_Troe;
	BzzVectorInt iCABR_SRI;
	BzzVectorInt iFallOff_Lindemann;
	BzzVectorInt iFallOff_Troe;
	BzzVectorInt iFallOff_SRI;
	BzzVectorInt iLandauTeller;
	BzzVectorInt iPowerSeries;
	BzzVectorInt iJanevLanger;
	BzzVectorInt iChebishevPolynomials;
	BzzVectorInt iPressureLogarithmic;
	BzzVectorInt iCollisionEfficiency;

	BzzVectorInt iReversible;
	BzzVectorInt iReverseRate;
	BzzVectorInt iDuplicate;
	BzzVectorInt iGlobal;

	BzzVectorInt	*species_reaction;
	BzzVector		*species_nu;

	BzzVectorInt	numDir1;
	BzzVectorInt	numDir2;
	BzzVectorInt	numDir3;
	BzzVectorInt	numDir4;
	BzzVectorInt	numDir5;
	BzzVectorInt	numInvTot1;
	BzzVectorInt	numInvTot2;
	BzzVectorInt	numInvTot3;
	BzzVectorInt	numInvTot4;
	BzzVectorInt	numInvTot5;
	BzzVectorInt	numInvEq1;
	BzzVectorInt	numInvEq2;
	BzzVectorInt	numInvEq3;
	BzzVectorInt	numInvEq4;
	BzzVectorInt	numInvEq5;

	BzzVectorInt	jDir1;
	BzzVectorInt	jDir2;
	BzzVectorInt	jDir3;
	BzzVectorInt	jDir4;
	BzzVectorInt	jDir5;
	BzzVectorInt	jInvTot1;
	BzzVectorInt	jInvTot2;
	BzzVectorInt	jInvTot3;
	BzzVectorInt	jInvTot4;
	BzzVectorInt	jInvTot5;
	BzzVectorInt	jInvEq1;
	BzzVectorInt	jInvEq2;
	BzzVectorInt	jInvEq3;
	BzzVectorInt	jInvEq4;
	BzzVectorInt	jInvEq5;

	BzzVector		valDir5;
	BzzVector		valInvTot5;
	BzzVector		valInvEq5;
	BzzVector		sumNuij;
	BzzVector		forwardOrders;
	BzzVector		backwardOrders;

	// Reaction orders

	BzzVectorInt	*lambda_species_reaction;
	BzzVector		*lambda_species_nu;

	BzzVectorInt	lambda_numDir1;
	BzzVectorInt	lambda_numDir2;
	BzzVectorInt	lambda_numDir3;
	BzzVectorInt	lambda_numDir4;
	BzzVectorInt	lambda_numDir5;
	BzzVectorInt	lambda_numInvEq1;
	BzzVectorInt	lambda_numInvEq2;
	BzzVectorInt	lambda_numInvEq3;
	BzzVectorInt	lambda_numInvEq4;
	BzzVectorInt	lambda_numInvEq5;

	BzzVectorInt	lambda_jDir1;
	BzzVectorInt	lambda_jDir2;
	BzzVectorInt	lambda_jDir3;
	BzzVectorInt	lambda_jDir4;
	BzzVectorInt	lambda_jDir5;
	BzzVectorInt	lambda_jInvEq1;
	BzzVectorInt	lambda_jInvEq2;
	BzzVectorInt	lambda_jInvEq3;
	BzzVectorInt	lambda_jInvEq4;
	BzzVectorInt	lambda_jInvEq5;

	BzzVector		lambda_valDir5;
	BzzVector		lambda_valInvEq5;

private:
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_KINETICSDATA_H)
