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

#include "OpenSMOKE_LookUp_Table_Manager.h"
#include "flamelet_group.h"
#include "Utilities.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//									OpenSMOKE_LookUp_Table_Manager							              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_LookUp_Table_Manager::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_Manager"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LookUp_Table_Manager::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_Manager"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_LookUp_Table_Manager::OpenSMOKE_LookUp_Table_Manager()
{
	name_object						= "[not assigned]";

	iImportingFlameletFiles				= false;	// Step 01: from OpenSMOKE flamelets to FLUENT flamelets
	iFlameletLibraryConstruction		= false;	// Step 02: applying the Beta PDF
	iWriteLookUpTableTest				= false;	// Step 03: write PDF Look Up Table 
	iWriteLookUpTableFLUENT				= false;	// Step 04: (optional) write PDF Look Up Table for FLUENT
	
	iSootSourceTermsLibraryConstruction = false;	// Uncorrelated approach

	iSootProfilesLibraryConstruction	= false;	// Correlated approach
	iImportingSootFlameletFiles			= false;

	verboseFormationRates = false;

	number_of_mixture_fraction_variance_points	= 32;
	stretching_factor							= 1.18;

	nucleation_model	= NUCLEATION_NONE;
	growth_model		= GROWTH_NONE;
	aggregation_model   = AGGREGATION_NONE;
	oxidation_model     = OXIDATION_NONE;

	nProgressVariable   = 0;
}

void OpenSMOKE_LookUp_Table_Manager::SetKineticScheme(OpenSMOKE_ReactingGas &_mix)
{
    mix = &_mix;
}

void OpenSMOKE_LookUp_Table_Manager::SetFlameletList(const vector<string> _flamelet_names)
{
	flamelet_names = _flamelet_names;
	iImportingFlameletFiles = true;
}

void OpenSMOKE_LookUp_Table_Manager::SetSootFlameletList(const vector<string> _soot_flamelet_names)
{
	soot_flamelet_names = _soot_flamelet_names;
	iImportingSootFlameletFiles = true;
}

void OpenSMOKE_LookUp_Table_Manager::SetFlameletLibraryName(const string  _flamelet_library_name)
{
	flamelet_library_name = _flamelet_library_name;
	iFlameletLibraryConstruction = true;
}

void OpenSMOKE_LookUp_Table_Manager::SetLookUpTableName(const string  _lookup_table_library_name)
{
	lookup_table_library_name = _lookup_table_library_name;
}

void OpenSMOKE_LookUp_Table_Manager::SetSootSourceTermsLibrary()
{
	iSootSourceTermsLibraryConstruction = true;
}

void OpenSMOKE_LookUp_Table_Manager::SetSootProfilesLibrary()
{
	iSootProfilesLibraryConstruction = true;
}

void OpenSMOKE_LookUp_Table_Manager::SetFLUENTLookUpTableName(const string  _fluent_lookup_table_library_name)
{
	fluent_lookup_table_library_name = _fluent_lookup_table_library_name;
	iWriteLookUpTableFLUENT = true;
}

void OpenSMOKE_LookUp_Table_Manager::SetLookUpTableFolderName(const string  _lookup_table_folder_name)
{
	lookup_table_folder_name = _lookup_table_folder_name;
	iWriteLookUpTableTest = true;
}

void OpenSMOKE_LookUp_Table_Manager::SetMixtureFractionVariancePoints(const int _number_of_mixture_fraction_variance_points)
{
	number_of_mixture_fraction_variance_points = _number_of_mixture_fraction_variance_points;
}

void OpenSMOKE_LookUp_Table_Manager::SetStretchingFactor(const double _stretching_factor)
{
	stretching_factor = _stretching_factor;
}

void OpenSMOKE_LookUp_Table_Manager::SetKindOfPDF(const string _kind_of_pdf)
{
	kind_of_pdf = _kind_of_pdf;
}

void OpenSMOKE_LookUp_Table_Manager::SetSootClosureModel(const string _iSootMode)
{
		 if (_iSootMode == "CORRELATED")	iPerfectlyCorrelated = true;
	else if (_iSootMode == "UNCORRELATED")	iPerfectlyCorrelated = false;
	else ErrorMessage("Wrong soot closure model...");
}

void OpenSMOKE_LookUp_Table_Manager::SetNucleationModel(const string nucleation)
{
		 if ( nucleation == "None" )						nucleation_model = NUCLEATION_NONE;
	else if ( nucleation == "Liu_2001" )					nucleation_model = NUCLEATION_LIU_2001;
	else if ( nucleation == "Liu_2002" )					nucleation_model = NUCLEATION_LIU_2002;		
	else if ( nucleation == "Liu_2003" )					nucleation_model = NUCLEATION_LIU_2003;		
	else if ( nucleation == "Moss_1999" )					nucleation_model = NUCLEATION_MOSS_1999;		
	else if ( nucleation == "Wen_2003" )					nucleation_model = NUCLEATION_WEN_2003;		
	else if ( nucleation == "Lindstedt_1994" )				nucleation_model = NUCLEATION_LINDSTEDT_1994;		
	else if ( nucleation == "Leung_1991" )					nucleation_model = NUCLEATION_LEUNG_1991;		
	else if ( nucleation == "Hall_1997" )					nucleation_model = NUCLEATION_HALL_1997;		
	else ErrorMessage("Wrong nucleation model...");
}

void OpenSMOKE_LookUp_Table_Manager::SetGrowthModel(const string growth)
{
		 if ( growth == "None" )							growth_model = GROWTH_NONE;
	else if ( growth == "Liu_2001" )						growth_model = GROWTH_LIU_2001;
	else if ( growth == "Liu_2002" )						growth_model = GROWTH_LIU_2002;		
	else if ( growth == "Liu_2003" )						growth_model = GROWTH_LIU_2003;		
	else if ( growth == "Moss_1999" )						growth_model = GROWTH_MOSS_1999;		
	else if ( growth == "Wen_2003" )						growth_model = GROWTH_WEN_2003;		
	else if ( growth == "Lindstedt_1994" )					growth_model = GROWTH_LINDSTEDT_1994;		
	else if ( growth == "Leung_1991" )						growth_model = GROWTH_LEUNG_1991;		
	else ErrorMessage("Wrong growth model...");
}


void OpenSMOKE_LookUp_Table_Manager::SetAggregationModel(const string aggregation)
{
		 if ( aggregation == "None" )						aggregation_model = AGGREGATION_NONE;
	else if ( aggregation =="Smoluchowski" )				aggregation_model = AGGREGATION_SMOLUCHOWSKI;
	else if ( aggregation == "Moss" )						aggregation_model = AGGREGATION_MOSS;		
	else ErrorMessage("Wrong aggregation model...");
}

void OpenSMOKE_LookUp_Table_Manager::SetOxidationModel(const string oxidation)
{
	     if ( oxidation == "None" )							oxidation_model = OXIDATION_NONE;
	else if ( oxidation == "Lee" )							oxidation_model = OXIDATION_LEE;
	else if ( oxidation == "NSC" )							oxidation_model = OXIDATION_NSC;		
	else if ( oxidation == "Neoh" )							oxidation_model = OXIDATION_NEOH;		
	else ErrorMessage("Wrong oxidation model...");
}

void OpenSMOKE_LookUp_Table_Manager::SetProgressVariable(const vector<string> _names)
{
	nProgressVariable = atoi(_names[0].c_str());
	progress_variable_index = new BzzVectorInt[nProgressVariable+1];
	progress_variable_value = new BzzVector[nProgressVariable+1];
	progress_variable_name.resize(nProgressVariable+1);

	int j=1;
	for(int k=1;k<=nProgressVariable;k++)
	{
		progress_variable_name[k] = _names[j++];
		int n = atoi(_names[j++].c_str());
		ChangeDimensions(n, &progress_variable_index[k]);
		ChangeDimensions(n, &progress_variable_value[k]);
		for(int i=1;i<=n;i++)
		{
			progress_variable_index[k][i] = mix->recognize_species(_names[j++]);
			progress_variable_value[k][i] = atof(_names[j++].c_str());
		}

		cout << "Progress Variable: " << progress_variable_name[k] << endl;
		for(int i=1;i<=n;i++)
			cout << mix->names[progress_variable_index[k][i]] << " " << progress_variable_value[k][i] << endl;
		cout << endl;
	}
}

void OpenSMOKE_LookUp_Table_Manager::DefineFromFile(const string inputFile)
{
    double			double_value;
    string			string_value;
    int				int_value;
	vector<string>  string_vector;

    OpenSMOKE_Dictionary_LookUp_Table dictionary;
    dictionary.ParseFile(inputFile);

	if (dictionary.Return("#ListOfFlamelets", string_vector))
		SetFlameletList(string_vector);

	if (dictionary.Return("#ListOfSootFlamelets", string_vector))
		SetSootFlameletList(string_vector);

	if (dictionary.Return("#FlameletLibrary", string_value))
		SetFlameletLibraryName(string_value);

	if (dictionary.Return("#LookUpTable", string_value))
		SetLookUpTableName(string_value);

	if (dictionary.Return("#LookUpTableOpenFOAM", string_value))
		SetLookUpTableFolderName(string_value);

	if (dictionary.Return("#LookUpTableForFLUENT", string_value))
		SetFLUENTLookUpTableName(string_value);

	if (dictionary.Return("#Pdf", string_value))
		SetKindOfPDF(string_value);

	if (dictionary.Return("#NumberOfVariances", int_value))
		SetMixtureFractionVariancePoints(int_value);

	if (dictionary.Return("#StretchingFactor", double_value))
		SetStretchingFactor(double_value);

	if (dictionary.Return("#SootClosure", string_value))
		SetSootClosureModel(string_value);

	if (dictionary.Return("#SootSources"))
		SetSootSourceTermsLibrary();

	if (dictionary.Return("#SootProfiles"))
		SetSootProfilesLibrary();

	if (dictionary.Return("#NucleationModel", string_value))
		SetNucleationModel(string_value);

	if (dictionary.Return("#GrowthModel", string_value))
		SetGrowthModel(string_value);

	if (dictionary.Return("#AggregationModel", string_value))
		SetAggregationModel(string_value);

	if (dictionary.Return("#OxidationModel", string_value))
		SetOxidationModel(string_value);

	if (dictionary.Return("#ProgressVariable", string_vector))
		SetProgressVariable(string_vector);
}

void OpenSMOKE_LookUp_Table_Manager::Run()
{
	if( iImportingFlameletFiles == true )					ImportFlamelets();
	if( iFlameletLibraryConstruction == true )				BuildFlameletLibrary();	
	if( iWriteLookUpTableTest == true )						WriteLookUpTableTest();
	if( iWriteLookUpTableFLUENT == true )					WriteFLUENTLookUpTable();
	if( iSootSourceTermsLibraryConstruction == true )		SootSourceTermsLibraryConstruction();
	if( iSootProfilesLibraryConstruction == true )			SootProfilesLibraryConstruction();
}

void OpenSMOKE_LookUp_Table_Manager::ImportFlamelets()
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, flamelet_library_name);
	fOutput.setf(ios::scientific);

	for(int i=1;i<=int(flamelet_names.size());i++)
	{
		cout << "Flamelet n." << i << " - " << flamelet_names[i-1] << endl;
		conversion_from_FlameCRECK_to_FLUENT(*mix, flamelet_names[i-1], fOutput, verboseFormationRates);
		cout << "DONE!" << endl;
	}
	fOutput << endl;
	fOutput << endl;
	fOutput.close();
}

void OpenSMOKE_LookUp_Table_Manager::BuildFlameletLibrary()
{
	flamelet_group flameletGroup;
	import_flamelet_file_from_FLUENT(*mix, flamelet_library_name, flameletGroup);
	flameletGroup.write_progress_variables(*mix, nProgressVariable, progress_variable_name, progress_variable_index, progress_variable_value);
	flameletGroup.build_library(*mix, kind_of_pdf, lookup_table_library_name, number_of_mixture_fraction_variance_points, stretching_factor);
}

void OpenSMOKE_LookUp_Table_Manager::WriteLookUpTableTest()
{
	flamelet_group flameletGroup;
	flameletGroup.read_flamelet_library(*mix, lookup_table_library_name);
	flameletGroup.print_on_file_look_up_table(lookup_table_folder_name, *mix);	
}

void OpenSMOKE_LookUp_Table_Manager::WriteFLUENTLookUpTable()
{
	flamelet_group flameletGroup;
	flameletGroup.read_flamelet_library(*mix, lookup_table_library_name);
	flameletGroup.print_on_file_FLUENT(fluent_lookup_table_library_name, *mix);
}

void OpenSMOKE_LookUp_Table_Manager::SootSourceTermsLibraryConstruction()
{
	BzzVector *m0_N, *fv_N;
	flamelet_group flameletGroup;
	flameletGroup.read_flamelet_library(*mix, lookup_table_library_name);
	flameletGroup.apply_source_betaPDF(*mix, iPerfectlyCorrelated, 4, m0_N, fv_N,
										nucleation_model, growth_model, aggregation_model, oxidation_model);

	//flameletGroup.soot_post_processing();
	//flameletGroup.apply_soot();
}	

void OpenSMOKE_LookUp_Table_Manager::SootProfilesLibraryConstruction()
{
	flamelet_group flameletGroup;
	BzzVector *m0;
	BzzVector *fv;
		
	// Memory Allocation for Normalized Profiles
	// ------------------------------------------------------------------------------
	cout << "Total number of flamelets: " << soot_flamelet_names.size() << endl;
	m0 = new BzzVector[soot_flamelet_names.size()+1]; 
	fv = new BzzVector[soot_flamelet_names.size()+1]; 

	// Read Soot Properties From CRECK Files: m0 and fv
	// ------------------------------------------------------------------------------		
	for(int i=1;i<=int(soot_flamelet_names.size());i++)
	{
		cout << "Reading strain rate " << i << "... ";
		read_SootProperties_from_FlameCRECK(soot_flamelet_names[i-1].c_str(), m0[i], fv[i]);
		cout << "DONE!! " << endl;
	}

	// Normalization of soot properties profiles
	// ----------------------------------------------------------------------------
	if (iPerfectlyCorrelated == true)
	{
		cout << "Normalization of soot properties profiles..." << endl;
		for (int i=1;i<=int(soot_flamelet_names.size());i++)
		{
			double m0_MaximumValue = m0[i].Max();
			double fv_MaximumValue = fv[i].Max();
			
			m0[i] /= m0_MaximumValue;
			fv[i] /= fv_MaximumValue;

			cout << "Flamelet " << i << " - Max Soot Volume Fraction " << fv_MaximumValue << endl; 
		}
		cout << "Successfully DONE!" << endl;
	}

	// Importing Flamelet Library
	// ------------------------------------------------------------------------------
	cout << "Importing Flamelet Library... ";
	flameletGroup.read_flamelet_library(*mix, lookup_table_library_name);
	cout << "Successfully DONE!!" << endl;

	// Generation of soot (normalized) profiles - Library
	// ------------------------------------------------------------------------------
	cout << "Generation of soot (normalized) profiles - Library" << endl;
	flameletGroup.apply_source_betaPDF_SootProperties(*mix, 2, m0, fv);	
	cout << "Successfully DONE!!" << endl;

	// Generation of source terms - Library
	// ------------------------------------------------------------------------------
	if (iPerfectlyCorrelated == true)
	{
		cout << "Generation of source terms - Library" << endl;
		flameletGroup.apply_source_betaPDF(*mix, iPerfectlyCorrelated, 4, m0, fv,
											nucleation_model, growth_model, aggregation_model, oxidation_model);	
		cout << "Successfully DONE!!" << endl;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											DICTIONARY										              //
////////////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_LookUp_Table::OpenSMOKE_Dictionary_LookUp_Table()
{
    SetupBase();

    Add("#ListOfFlamelets",				'O', 'V', "List of flamelets");
	Add("#Pdf",							'O', 'S', "PDF: BETA, GAUSSIAN, DIRAC");
	Add("#LookUpTable",					'O', 'S', "LookUp Table name");
	Add("#FlameletLibrary",				'O', 'S', "Flamelet library name (FLUENT format)");

	Add("#LookUpTableOpenFOAM",			'O', 'S', "LookUp Table Test folder name");
	Add("#LookUpTableForFLUENT",		'O', 'S', "LookUp Table name (FLUENT format)");

	Add("#NumberOfVariances",			'O', 'I', "Number of points in the mixture fraction variance space");
	Add("#StretchingFactor",			'O', 'D', "Stretching factor");

	Add("#ListOfSootFlamelets",			'O', 'V', "List of soot flamelets");
	Add("#SootClosure",					'O', 'S', "Soot closure model: CORRELATED || UNCORRELATED");
	Add("#SootSources",					'O', 'N', "Soot source library");
	Add("#SootProfiles",				'O', 'N', "Soot profile library");

	Add("#NucleationModel",				'O', 'S', "Nucleation model: None Liu_2001 Liu_2002 Liu_2003 Moss_1999 Wen_2003 Lindstedt_1994 Leung_1991 Hall_1997");
	Add("#GrowthModel",					'O', 'S', "Growth model: None Liu_2001 Liu_2002 Liu_2003 Moss_1999 Wen_2003 Lindstedt_1994 Leung_1991");
	Add("#AggregationModel",			'O', 'S', "Growth model: None Smoluchowski Moss");
	Add("#OxidationModel",				'O', 'S', "Growth model: None Lee NSC Neoh");

	Add("#ProgressVariable",			'O', 'V', "Progress variables: nPV, [name, N, species, mass fraction]");

    Lock();
}