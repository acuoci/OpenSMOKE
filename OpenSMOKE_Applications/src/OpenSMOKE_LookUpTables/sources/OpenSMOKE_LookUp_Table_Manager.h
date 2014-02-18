/***************************************************************************
 *   Copyright (C) 2006-2009 by Alberto Cuoci						       *
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

#if !defined(OPENSMOKE_DICTIONARY_LOOKUP_TABLE)
#define OPENSMOKE_DICTIONARY_LOOKUP_TABLE

#include "engine/OpenSMOKE_ReactingGas.h"
#include "basic/OpenSMOKE_Dictionary.h"

class OpenSMOKE_Dictionary_LookUp_Table: public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_LookUp_Table();
};



class OpenSMOKE_LookUp_Table_Manager
{
public:

	OpenSMOKE_LookUp_Table_Manager();

private:

	bool iImportingFlameletFiles;		// from OpenSMOKE flamelets to FLUENT flamelets
	bool iFlameletLibraryConstruction;	// applying the PDF
	
	bool iWriteLookUpTableTest;			// write PDF LookUp Table (test) 
	bool iWriteLookUpTableFLUENT;		// write PDF LookUp Table for FLUENT (optional)

	bool iSootSourceTermsLibraryConstruction;
	bool iSootProfilesLibraryConstruction;

	bool iImportingSootFlameletFiles;

private:

	vector<string> flamelet_names;
	string flamelet_library_name;
	string lookup_table_library_name;
	string lookup_table_folder_name;
	string fluent_lookup_table_library_name;
	vector<string> soot_flamelet_names;
	
	int number_of_mixture_fraction_variance_points;
	double stretching_factor;
	string kind_of_pdf;
	bool iPerfectlyCorrelated;
	bool verboseFormationRates;

	OpenSMOKE_ReactingGas *mix;

	int nProgressVariable;
	BzzVectorInt*  progress_variable_index;
	BzzVector*     progress_variable_value;
	vector<string> progress_variable_name;

private:

	void SetFlameletList(const vector<string> _flamelet_names);

	void SetSootFlameletList(const vector<string> _soot_flamelet_names);
	void SetFlameletLibraryName(const string  _flamelet_library_name);
	void SetLookUpTableName(const string _lookup_table_library_name);
	void SetLookUpTableFolderName(const string  _lookup_table_folder_name);
	void SetFLUENTLookUpTableName(const string _fluent_lookup_table_library_name);
	void SetMixtureFractionVariancePoints(const int _number_of_mixture_fraction_variance_points);
	void SetStretchingFactor(const double _stretching_factor);
	void SetKindOfPDF(const string _kind_of_pdf);
	void SetSootClosureModel(const string _iSootMode);
	void SetSootSourceTermsLibrary();
	void SetSootProfilesLibrary();
	void SetProgressVariable(const vector<string> _names);

public:

	void DefineFromFile(const string file_name);
	void SetKineticScheme(OpenSMOKE_ReactingGas &_mix);
	void SetVerboseFormationRates(const bool value) { verboseFormationRates = value; }
	void Run();
	
private:
	void ImportFlamelets();
	void BuildFlameletLibrary();
	void WriteLookUpTableTest();
	void WriteFLUENTLookUpTable();

	void SootSourceTermsLibraryConstruction();
	void SootProfilesLibraryConstruction();

private:

	nucleation_models	nucleation_model;
	growth_models		growth_model;
	aggregation_models	aggregation_model;
	oxidation_models	oxidation_model;

	void SetNucleationModel(const string nucleation);
	void SetGrowthModel(const string growth);
	void SetAggregationModel(const string aggregation);
	void SetOxidationModel(const string oxidation);

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_DICTIONARY_LOOKUP_TABLE)

