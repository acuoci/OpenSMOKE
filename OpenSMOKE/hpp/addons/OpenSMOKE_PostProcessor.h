/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
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

#ifndef OPENSMOKE_POSTPROCESSOR
#define OPENSMOKE_POSTPROCESSOR

#include <sstream>
#include "BzzMath.hpp"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_PostProcessor_RateOfProductionAnalysis;

class OpenSMOKE_PostProcessor_General;
class OpenSMOKE_PostProcessor_Flame1D;
class OpenSMOKE_PostProcessor_PFR;
class OpenSMOKE_PostProcessor_CSTR;
class OpenSMOKE_PostProcessor_Batch;
class OpenSMOKE_PostProcessor_ShockTube;
class OpenSMOKE_PostProcessor_ICEM_MultiZone;

class OpenSMOKE_PostProcessor_SensitivityAnalysis_General;
class OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D;

enum kind_of_reacting_system {POST_PROCESSOR_NONE, POST_PROCESSOR_FLAME1D, POST_PROCESSOR_PFR, POST_PROCESSOR_CSTR, POST_PROCESSOR_BATCH, POST_PROCESSOR_FLAMELET, POST_PROCESSOR_SHOCKTUBE, POST_PROCESSOR_ICEM_MULTIZONE};

class OpenSMOKE_PostProcessor
{
public:

	OpenSMOKE_PostProcessor();
	void SetName(const std::string name);

	void ReadFromBinaryFile(const std::string fileName);

	int NC;							// number of species
	int NR;							// number of reactions
	std::string *names;					// names of species
	std::string *reactions;				// names of reactions
	BzzVector M;
	kind_of_reacting_system kind;
	OpenSMOKE_ReactingGas	*mix;	// Pointer to gas mixture

	std::string ExportMainXAxis();
	void ExportAvailableXAxis(vector<string> &x_axis);
	void ExportAvailableYAxis(vector<string> &y_axis);
	void ExportAvailableYAxisReactionRates(vector<string> &y_axis);
	void ExportAvailableYAxisFormationRates(vector<string> &y_axis);
	void ExportAvailableYAxisSensitivityCoefficients(vector<string> &y_axis);


	void ExportAllSpeciesNames(vector<string> &names_species);
	void ExportAdditionalNames(vector<string> &additional_names, const int index);
	void ExportSensitivitySpeciesNames(vector<string> &species_names, const int index);
	
	// ROPA
	void ImportSelectedAxis(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, std::string &name_x, std::string &name_y, vector<string> &names_lines);
	void ImportSelectedAxisReactionRates(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, std::string &name_x, std::string &name_y, vector<string> &names_lines);
	void ImportSelectedAxisFormationRates(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, std::string &name_x, std::string &name_y, vector<string> &names_lines);
	void ImportSelectedAxisSensitivityCoefficients(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, std::string &name_x, std::string &name_y, vector<string> &names_lines);
	
	void ImportIntegralROPA(const int index, vector<double> &p, vector<double> &d, vector<int> &ip, vector<int> &id, vector<string> &names_p, vector<string> &names_d);
	void ImportIntegralROPA(const int index, vector<double> &t, vector<int> &it, vector<string> &names_t);
	void ImportUnimportantReactions(const double eps_threshold_percentage, stringstream &string_out);
	void CallIntegralAnalysis();
	void CallLocalAnalysis(const double local_coordinate);
	void CallRegionAnalysis(const double xA, const double xB);

	// Sensitivity - Bars
	void ImportAdditionalSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, vector<double> &t, vector<int> &it, vector<string> &names_t, const int index, const int index_selection);
	void ImportMassFractionSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, const std::string name, vector<double> &t, vector<int> &it, vector<string> &names_t, const int index);

	// Sensitivity - Profiles
	void ImportAdditionalSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t, const int index, const int index_selection);
	void ImportMassFractionSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, const std::string name,BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t, const int index);

	void SetFocusSensitivityProfiles(const bool additional_index, const int index);

	bool	focus_sensitivity_additional;
	int		focus_sensitivity_index;

	inline bool iROPA()						{ return iROPA_;} 
	inline bool iSensitivity()				{ return iSensitivity_;} 
	inline bool iSensitivityDiffusivity()	{ return iSensitivity_Diffusivity_;} 
	inline bool iChart()					{ return iChart_;} 

	inline const OpenSMOKE_PostProcessor_General	 &get_general()				{return *general;}
	inline const OpenSMOKE_PostProcessor_Flame1D	 &get_flame()				{return *flame1d;}
	
	double	get_T(const int j);
	double	get_P_Pa(const int j);
	double	get_Ctot(const int j);
	const BzzVector get_C(const int j);
	const BzzVector get_omega_profile(const int j_species);
	const BzzVector get_temperature_profile();

	BzzVector get_x();

	std::string sensitivity_list(const int index);
	int startSensitivityAdditional;

private:

	OpenSMOKE_PostProcessor_Flame1D			*flame1d;
	OpenSMOKE_PostProcessor_PFR				*pfr;
	OpenSMOKE_PostProcessor_CSTR			*cstr;
	OpenSMOKE_PostProcessor_Batch			*batch;
	OpenSMOKE_PostProcessor_ShockTube		*shocktube;
	OpenSMOKE_PostProcessor_ICEM_MultiZone	*icem_multizone;
	OpenSMOKE_PostProcessor_General			*general;

	OpenSMOKE_PostProcessor_SensitivityAnalysis_General **sensitivity;
	OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D *sensitivity_frequency_factor;
	OpenSMOKE_PostProcessor_SensitivityAnalysis_Flame1D *sensitivity_diffusivity;

	OpenSMOKE_PostProcessor_RateOfProductionAnalysis *ropa;

	std::string nameFileLoad;
	BzzLoad fLoad;

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

	void Check(const std::string flag);
	std::string Next();

	bool iROPA_;
	bool iSensitivity_;
	bool iSensitivity_Diffusivity_;
	bool iChart_;
};

#endif // OPENSMOKE_POSTPROCESSOR

