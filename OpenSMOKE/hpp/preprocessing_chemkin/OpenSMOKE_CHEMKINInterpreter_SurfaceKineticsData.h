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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_SURFACEKINETICSDATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_SURFACEKINETICSDATA_H

#include "BzzMath.hpp"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData.h"

class OpenSMOKE_WarningFile;
class OpenSMOKE_CHEMKINInterpreter;
class OpenSMOKE_CHEMKINInterpreter_TransportData;

class OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData
{
public:
	
	OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData();
	virtual ~OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData();

	BzzVectorInt lineReactionIndex;
	BzzVectorInt additionalLineReactionIndex;

	int number_of_reactions;
	vector<string> gas_species_list;
	vector<string> site_species_list;
	vector<string> bulk_species_list;

	OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData *reaction;
	OpenSMOKE_CHEMKINInterpreter_SurfaceReactionData *reaction_reverserate;
	void Setup(ofstream *_fLog, OpenSMOKE_WarningFile *fWarning);
	void RecognizeSpecies(const string name_species, int &index, char &phase);
	void FinalChecks();
	void Statistics();
	void SummaryOnFile(ofstream &fOutput);
	void PrintBinaryFile(BzzSave &fOutput, BzzSave &fASCII,OpenSMOKE_CHEMKINInterpreter_ThermoData &thermo, OpenSMOKE_CHEMKINInterpreter &interp, OpenSMOKE_PreProcessorSurfaceMaterial &material);

private:

	ofstream *fLog;
	OpenSMOKE_WarningFile *fWarning;

	void CheckDuplicateReactions();
	void CheckReactionRatesReactions();

	BzzVectorInt iReversible;
	BzzVectorInt iDuplicate;
	BzzVectorInt iGlobal;

private:
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData_H)
