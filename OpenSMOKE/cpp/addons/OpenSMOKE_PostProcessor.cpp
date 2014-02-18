/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci         	                       *
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

#include <vector>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "addons/OpenSMOKE_PostProcessor.h"
#include "addons/OpenSMOKE_PostProcessor_RateOfProductionAnalysis.h"
#include "addons/OpenSMOKE_PostProcessor_SensitivityAnalysis_General.h"
#include "addons/OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D.h"
#include "addons/OpenSMOKE_PostProcessor_Flame1D.h"
#include "addons/OpenSMOKE_PostProcessor_PFR.h"
#include "addons/OpenSMOKE_PostProcessor_CSTR.h"
#include "addons/OpenSMOKE_PostProcessor_Batch.h"
#include "addons/OpenSMOKE_PostProcessor_ShockTube.h"
#include "addons/OpenSMOKE_PostProcessor_ICEM_MultiZone.h"

void OpenSMOKE_PostProcessor::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PostProcessor"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_PostProcessor::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PostProcessor"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

void OpenSMOKE_PostProcessor::SetFocusSensitivityProfiles(const bool additional, const int index)
{
	focus_sensitivity_index = index;
	focus_sensitivity_additional = additional;
}

double	OpenSMOKE_PostProcessor::get_T(const int j)		
{
	return general->T[j];
}

double	OpenSMOKE_PostProcessor::get_P_Pa(const int j)	
{
	return general->P_Pa[j];
}

double	OpenSMOKE_PostProcessor::get_Ctot(const int j)	
{
	return general->Ctot[j];
}

const BzzVector OpenSMOKE_PostProcessor::get_omega_profile(const int j_species)
{ 
	BzzVector aux(general->omega[j_species].Size());
	for(int i=1;i<=general->omega[j_species].Size();i++)
		aux[i] = general->omega[j_species][i];
	return aux;
}

const BzzVector OpenSMOKE_PostProcessor::get_temperature_profile()
{ 
	BzzVector aux(general->T.Size());
	for(int i=1;i<=general->T.Size();i++)
		aux[i] = general->T[i];
	return aux;
}

const BzzVector OpenSMOKE_PostProcessor::get_C(const int j)	
{
	BzzVector aux(NC);
	for(int i=1;i<=NC;i++)
		aux[i] = general->C[i][j];
	return aux;
}

BzzVector OpenSMOKE_PostProcessor::get_x()
{ 
	BzzVector aux;
	aux = get_general().x;
	return aux;
}


void OpenSMOKE_PostProcessor::SetName(const string name)
{
	name_object = name;
}

OpenSMOKE_PostProcessor::OpenSMOKE_PostProcessor()
{
	name_object					= "name not assigned";
	kind						= POST_PROCESSOR_NONE;
	iROPA_						= false;
	iSensitivity_				= false;
	iSensitivity_Diffusivity_	= false;
	iChart_						= false;
}

void OpenSMOKE_PostProcessor::Check(const string flag)
{
	char dummy[Constants::NAME_SIZE];
	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	if (strcmp(dummy, flag.c_str()))	
		WarningMessage("Expected: " + flag + " - Found: " + string(dummy));
}

string OpenSMOKE_PostProcessor::Next()
{
	char dummy[Constants::NAME_SIZE];
	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	return dummy;
}

void OpenSMOKE_PostProcessor::ReadFromBinaryFile(const string fileName)
{
	// Kinetic Scheme name
	size_t found;
	found = fileName.find_last_of("/\\");
	string folder_path	= fileName.substr(0,found);
	mix = new OpenSMOKE_ReactingGas();
	mix->SetupBinary(folder_path);

	nameFileLoad = fileName;

	sensitivity = new OpenSMOKE_PostProcessor_SensitivityAnalysis_General *[2]; // 2 panels

	int j;
	char dummy[Constants::NAME_SIZE];
	char name[Constants::NAME_SIZE];
	char name_reaction[Constants::REACTION_NAME_SIZE];

	fLoad('*', nameFileLoad);

	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	string version = dummy;
	if (version != "V20100417")
		ErrorMessage("This version of post processing file is not supported: " + version);
	cout << "Version: " << version << endl;

	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	cout << "Date: "  << dummy << endl;

	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	cout << "Kinetic scheme: "  << dummy << endl;

	// Number of species
	Check("NC");
	fLoad >> NC;
	cout << "Number of species:   "  << NC << endl;

	// Number of reactions
	Check("NR");
	fLoad >> NR;
	cout << "Number of reactions: "  << NR << endl;

	// Molecular weights
	Check("M");
	fLoad >> M;

	// Species Names
	Check("SPECIES");
	names = new string[NC + 1];
	for(j=1;j<=NC;j++)
	{
		fLoad.fileLoad.read((char*) name, sizeof(name));
		names[j] = name;
	}

	// Reaction names
	Check("REACTIONS");
	reactions = new string[NR + 1];
	for(j=1;j<=NR;j++)
	{
		fLoad.fileLoad.read((char*) name_reaction, sizeof(name_reaction));
		reactions[j] = name_reaction;
	}

	{
		string tag = Next();
		cout << tag << endl;
		
		if (tag == "FLAME1D")
		{
			kind	= POST_PROCESSOR_FLAME1D;
			iChart_	= true;
			flame1d = new OpenSMOKE_PostProcessor_Flame1D(this);
			general = flame1d;
		}
		else if (tag == "PFR")
		{
			kind	= POST_PROCESSOR_PFR;
			iChart_	= true;
			pfr = new OpenSMOKE_PostProcessor_PFR(this);
			general = pfr;
		}
		else if (tag == "CSTR")
		{
			kind	= POST_PROCESSOR_CSTR;
			iChart_	= true;
			cstr = new OpenSMOKE_PostProcessor_CSTR(this);
			general = cstr;
		}
		else if (tag == "BATCH")
		{
			kind	= POST_PROCESSOR_BATCH;
			iChart_	= true;
			batch = new OpenSMOKE_PostProcessor_Batch(this);
			general = batch;
		}
		else if (tag == "SHOCKTUBE")
		{
			kind	= POST_PROCESSOR_SHOCKTUBE;
			iChart_	= true;
			shocktube = new OpenSMOKE_PostProcessor_ShockTube(this);
			general = shocktube;
		}
		else if (tag == "ICEM-MULTI")
		{
			kind	= POST_PROCESSOR_ICEM_MULTIZONE;
			iChart_	= true;
			icem_multizone = new OpenSMOKE_PostProcessor_ICEM_MultiZone(this);
			general = icem_multizone;
		}
		else
			ErrorMessage("No post processing available for the following system: " + tag);

		general->ReadFromBinaryFile(fLoad);
	}
	
	for(;;)
	{
		string tag = Next();
		cout << tag << endl;
		if (tag == "END") 
			break;
		else if (tag == "ROPA")
		{
			iROPA_ = true;
			ropa = new OpenSMOKE_PostProcessor_RateOfProductionAnalysis(this);
			ropa->ReadFromBinaryFile(fLoad);
		}
		else if (tag == "SENSITIVITY")
		{
			iSensitivity_ = true;
			sensitivity_frequency_factor = new OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D(this);
			sensitivity[0] = sensitivity_frequency_factor;
			sensitivity[0]->ReadFromBinaryFile(fLoad, 0 );
		}
		else if (tag == "SENSITIVITY-DIFF")
		{
			iSensitivity_Diffusivity_ = true;
			sensitivity_diffusivity = new OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D(this);
			sensitivity[1] = sensitivity_diffusivity;
			sensitivity[1]->ReadFromBinaryFile(fLoad, 1 );
		}
		else
			ErrorMessage("No post processing available for the following option: " + tag);
	}

	fLoad.End();
}

void OpenSMOKE_PostProcessor::ExportAvailableXAxis(vector<string> &x_axis)
{
	general->ExportAvailableXAxis(x_axis);
}

void OpenSMOKE_PostProcessor::ExportAvailableYAxis(vector<string> &y_axis)
{
	general->ExportAvailableYAxis(y_axis);
}

void OpenSMOKE_PostProcessor::ExportAvailableYAxisReactionRates(vector<string> &y_axis)
{
	ropa->ExportAvailableYAxisReactionRates(y_axis);
}

void OpenSMOKE_PostProcessor::ExportAvailableYAxisFormationRates(vector<string> &y_axis)
{
	ropa->ExportAvailableYAxisFormationRates(y_axis);
}

void OpenSMOKE_PostProcessor::ExportAvailableYAxisSensitivityCoefficients(vector<string> &y_axis)
{
	// TODO
	sensitivity[0]->ExportAvailableYAxisSensitivityCoefficients(y_axis);
}

string OpenSMOKE_PostProcessor::ExportMainXAxis()
{
	return general->ExportMainXAxis();
}

void OpenSMOKE_PostProcessor::CallIntegralAnalysis()
{
	ropa->CallIntegralAnalysis();
}

void OpenSMOKE_PostProcessor::CallLocalAnalysis(const double local_coordinate)
{
	ropa->CallLocalAnalysis(local_coordinate);
}
	
void OpenSMOKE_PostProcessor::CallRegionAnalysis(const double xA, const double xB)
{
	ropa->CallRegionAnalysis(xA, xB);
}

void OpenSMOKE_PostProcessor::ExportAllSpeciesNames(vector<string> &species_names)
{
	species_names.resize(NC);
	for(int j=1;j<=NC;j++)
		species_names[j-1] = names[j];;
}

void OpenSMOKE_PostProcessor::ExportAdditionalNames(vector<string> &additional_names, const int index)
{	
	sensitivity[index]->ExportAdditionalNames(additional_names);
}

void OpenSMOKE_PostProcessor::ExportSensitivitySpeciesNames(vector<string> &species_names, const int index)
{
	sensitivity[index]->ExportSensitivitySpeciesNames(species_names);
}

void OpenSMOKE_PostProcessor::ImportSelectedAxis(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, string &name_x, string &name_y, vector<string> &names_lines)
{
	general->ImportSelectedAxis(x_axis, y_axis, x, y, name_x, name_y, names_lines);
}

void OpenSMOKE_PostProcessor::ImportSelectedAxisReactionRates(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, string &name_x, string &name_y, vector<string> &names_lines)
{
	ropa->ImportSelectedAxisReactionRates(x_axis, y_axis, x, y, name_x, name_y, names_lines);
}

void OpenSMOKE_PostProcessor::ImportSelectedAxisFormationRates(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, string &name_x, string &name_y, vector<string> &names_lines)
{
	ropa->ImportSelectedAxisFormationRates(x_axis, y_axis, x, y, name_x, name_y, names_lines);
}

void OpenSMOKE_PostProcessor::ImportSelectedAxisSensitivityCoefficients(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, string &name_x, string &name_y, vector<string> &names_lines)
{
	// TODO
	sensitivity[0]->ImportSelectedAxisSensitivityCoefficients(x_axis, y_axis, x, y, name_x, name_y, names_lines);
}

void OpenSMOKE_PostProcessor::ImportIntegralROPA(const int index, vector<double> &p, vector<double> &d, vector<int> &ip, vector<int> &id, vector<string> &names_p, vector<string> &names_d)
{
	ropa->MostImportantReactions(index, p, d, ip, id, names_p, names_d);
}

void OpenSMOKE_PostProcessor::ImportIntegralROPA(const int index, vector<double> &t, vector<int> &it, vector<string> &names_t)
{
	ropa->MostImportantReactions(index, t, it, names_t);
}

void OpenSMOKE_PostProcessor::ImportUnimportantReactions(const double eps_threshold_percentage, stringstream &string_out)
{
	ropa->UnimportantReactions(eps_threshold_percentage, string_out);
}

void OpenSMOKE_PostProcessor::ImportAdditionalSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, vector<double> &t, vector<int> &it, vector<string> &names_t, const int index, const int index_selection)
{
	sensitivity[index]->GetSensitivityBars(iTotal, iLocal, coordinate, index_selection+startSensitivityAdditional, t, it, names_t);
}

void OpenSMOKE_PostProcessor::ImportMassFractionSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, const string name, vector<double> &t, vector<int> &it, vector<string> &names_t, const int index)
{
	sensitivity[index]->GetMassFractionSensitivityBars(iTotal, iLocal, coordinate, name, t, it, names_t);
}

void OpenSMOKE_PostProcessor::ImportAdditionalSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t, const int index, const int index_selection)
{
	sensitivity[index]->GetSensitivityProfiles(iTotal, iLocal, coordinate, index_selection+startSensitivityAdditional, SMatrix, indices, names_t);
}

void OpenSMOKE_PostProcessor::ImportMassFractionSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, const string name,BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t, const int index)
{
	sensitivity[index]->GetMassFractionSensitivityProfiles(iTotal, iLocal,  coordinate, name, SMatrix, indices, names_t);
}

string OpenSMOKE_PostProcessor::sensitivity_list(const int index) 
{
	return sensitivity[index]->sensitivity_list();
}
