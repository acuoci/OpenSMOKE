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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_H)
#define OPENSMOKE_CHEMKININTERPRETER_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_WarningFile.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ElementsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_SpeciesData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_UnitsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_KineticsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_ThermoData.h"
#include "preprocessing_chemkin/OpenSMOKE_CHEMKINInterpreter_TransportData.h"

class OpenSMOKE_PreProcessorSurfaceMaterial
{
public:
	
	string name;

	int nSurfaceSites;
	int nSurfaceBulks;
	int nSurfaceReactions;

	vector<string> site_name;
	BzzVector  site_density;
	BzzVector *site_occupancy;	
	vector<string> *site_species_names;

	vector<string> bulk_name;
	BzzVector *bulk_density;	
	vector<string> *bulk_species_names;

	int start_material_section;
	int end_material_section;
	BzzVectorInt start_sites;
	BzzVectorInt start_bulks;
	BzzVectorInt start_reactions;

	BzzMatrix site_occupancy_matrix;
	BzzMatrix bulk_density_matrix;

	OpenSMOKE_CHEMKINInterpreter_SurfaceKineticsData	surfaceKinetics;
};

class OpenSMOKE_CHEMKINInterpreter  
{
friend class OpenSMOKE_PreProcessorIdealGas;

public:

	OpenSMOKE_CHEMKINInterpreter();
	virtual ~OpenSMOKE_CHEMKINInterpreter();

	bool iSootMode;
	int number_of_lines;
	vector<string> lines;
	int surface_number_of_lines;
	vector<string> surface_lines;

	void Run();
	void PrintIdealGasBinaryFile(const string file_name);
	void PrintSurfaceBinaryFile(const string file_name);

	OpenSMOKE_CHEMKINInterpreter_ElementsData	elements;
	OpenSMOKE_CHEMKINInterpreter_SpeciesData	species;
	OpenSMOKE_CHEMKINInterpreter_KineticsData	kinetics;
	OpenSMOKE_CHEMKINInterpreter_UnitsData		units;
	OpenSMOKE_CHEMKINInterpreter_ThermoData		thermo;
	OpenSMOKE_CHEMKINInterpreter_TransportData	transport;

	void SetMinimumTemperature(const double tmin);
	void SetMaximumTemperature(const double tmax);
	void SetFittingPoints(const int fittingpoints);
	void SetName(const string name);
	void SetAuthorName(const string name);
	void SetKineticFileName(const string name);
	void SetSurfaceKineticFileName(const string name);
	void SetThermodynamicFileName(const string name);
	void SetTransportFileName(const string name);
	void SetReducedListFileName(const string name);
	void SetOutputFolderName(const string name);
	void SetBuzziMode();
	void SetVerboseMode();
	void SetNoVerboseMode();
	void SetSootMode();

	string nameKineticsFile();
	string nameSurfaceKineticsFile();
	string nameThermoFile();
	string nameTransportFile();
	string nameOutputFolder();

private:

	void ParsingSpecies(const int iLine, vector<string> instructions);
	void ParsingElements(const int iLine, vector<string> instructions);
	void ParsingReactions(const int iReaction, const int iLine, vector<string> instructions);
	void ParsingReactionsAdditional(const int i, const int iLine, vector<string> instructions);
	void AddSpaces(string &line, const char symbol);

	// Surface specific parsing functions
	int nSurfaceMaterials;
	OpenSMOKE_PreProcessorSurfaceMaterial *material;

	void SurfaceParsingSections();
//	void SurfaceParsingUnits();
	void SurfaceParsingMaterials();
	void SurfaceParsingSites();
	void SurfaceParsingBulk();
	void SurfaceParsingReactions();
	void SurfaceParsingReactionsAdditionalData();
	void SurfaceParsingReactions(const int k, const int iReaction, const int iLine, vector<string> instructions);
	void SurfaceParsingReactionsAdditional(const int k, const int i, const int iLine, vector<string> instructions);

	void SurfaceSyntaxCorrections();
	void SurfaceCleaningLines();

	vector<string> bulk_species_global_list;
	vector<string> site_species_global_list;
	BzzVectorInt indices_site_thermo_database;
	BzzVectorInt indices_bulk_thermo_database;

	void WriteReducedScheme(string fileName);

	//

	void PrintSummaryFile(const string file_name);
	void PrintBuzziFiles();
	void WriteAdditionalFiles();

	void ReadKineticScheme();
	void ReadSurfaceKineticScheme();
	void ParsingSections();
	void SyntaxCorrections();
	void ParsingUnits();
	void CleaningLines();
	void ParsingElements();
	void ParsingSpecies();
	void ParsingReactions();
	void ParsingReactionsAdditionalData();
	void SurfaceSummaryOnFile();

	BzzVectorInt indexLines;
	BzzVectorInt indexBlankLines;
	BzzVectorInt indexCommentLines;

	BzzVectorInt surfaceIndexLines;
	BzzVectorInt surfaceIndexBlankLines;
	BzzVectorInt surfaceIndexCommentLines;

	BzzVectorInt indices_species_thermo_database;
	BzzVectorInt indices_species_transport_database;
	
	ofstream				fLog;
	OpenSMOKE_WarningFile	fWarning;

	int startElements;
	int startSpecies;
	int startReactions;
	BzzVectorInt startUnits;
	BzzVectorInt finalUnits;

private:

	double TMIN;
	double TMAX;
	int FITTINGPOINTS;

	bool iSurfaceKinetics;
	bool iBuzziMode;
	bool iVerboseMode;

	string name_author;
	string building_date;
	string name_object;

	string name_kinetics_file;
	string name_surface_kinetics_file;
	string name_thermo_file;
	string name_transport_file;
	string name_output_folder;
	string name_ascii_folder;
	string name_reduced_list_file;

	void ErrorMessage(const int iLine, const string message);
	void WarningMessage(const int iLine, const string message);
};




#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_H)
