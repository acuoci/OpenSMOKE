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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_THERMODATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_THERMODATA_H

#include "BzzMath.hpp"
#include <vector>

class OpenSMOKE_CHEMKINInterpreter_ThermoSpecies;
class OpenSMOKE_CHEMKINInterpreter_ElementsData;

class OpenSMOKE_CHEMKINInterpreter_ThermoData  
{
public:
	OpenSMOKE_CHEMKINInterpreter_ThermoData();
	virtual ~OpenSMOKE_CHEMKINInterpreter_ThermoData();

	void ReadThermoData(const string file_name, ofstream *_fLog);
	BzzVectorInt GiveMeSpeciesIndices(vector<string> &list);
	void PrintThermoFile(const string file_name, BzzVectorInt &indices);
	void GiveMeMolecularWeights(BzzVectorInt &indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements);
	void SummaryOnFile(ofstream &fOutput, BzzVectorInt &indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements);
	void SummaryOnFile(ofstream &fOutput, BzzVectorInt &site_indices, BzzVectorInt &bulk_indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements);

	void ProcessThermoData(BzzVectorInt &indices);
	void ProcessThermoData(BzzVectorInt &site_indices, BzzVectorInt &bulk_indices);
	void ProcessElementData(BzzVectorInt &indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements);
	void ProcessElementData(BzzVectorInt &site_indices, BzzVectorInt &bulk_indices, OpenSMOKE_CHEMKINInterpreter_ElementsData &elements);

	double ReducedMolecularWeight(const string name1, const string name2);

	void WriteSurfaceThermodynamicData(BzzSave &binaryFile, BzzSave &asciiFile, BzzVectorInt &site_indices, BzzVectorInt &bulk_indices);
	void WriteSurfaceElementsData(BzzSave &binaryFile, BzzSave &asciiFile);

	// Surface
	BzzVector site_elements_mw;
	BzzMatrix site_elements_matrix;
	BzzVector bulk_elements_mw;
	BzzMatrix bulk_elements_matrix;

	vector<string> lines;
	int number_of_lines;
	int total_number_of_species;
	BzzVectorInt indexBlankLines;
	BzzVectorInt indexCommentLines;
	BzzVectorInt indexLines;

	double tmin;
	double tmean;
	double tmax;

	OpenSMOKE_CHEMKINInterpreter_ThermoSpecies *species;
	OpenSMOKE_CHEMKINInterpreter_ThermoSpecies *surfaceSpecies;

	ofstream *fLog;

public:
	
	BzzVector *CpHT;
	BzzVector *CpLT;
	BzzVector *aDH;
	BzzVector *bDH;
	BzzVector *aDS;
	BzzVector *bDS;
	BzzVector T1;
	BzzVector T2;
	BzzVector T3;
	BzzVector M;

	BzzVector elements_mw;
	BzzMatrix elements_matrix;
	vector<string>	list_of_elements;


	BzzVector *site_CpHT;
	BzzVector *site_CpLT;
	BzzVector *site_aDH;
	BzzVector *site_bDH;
	BzzVector *site_aDS;
	BzzVector *site_bDS;
	BzzVector site_T1;
	BzzVector site_T2;
	BzzVector site_T3;
	BzzVector site_M;

	BzzVector *bulk_CpHT;
	BzzVector *bulk_CpLT;
	BzzVector *bulk_aDH;
	BzzVector *bulk_bDH;
	BzzVector *bulk_aDS;
	BzzVector *bulk_bDS;
	BzzVector bulk_T1;
	BzzVector bulk_T2;
	BzzVector bulk_T3;
	BzzVector bulk_M;

private:
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_THERMODATA_H)
