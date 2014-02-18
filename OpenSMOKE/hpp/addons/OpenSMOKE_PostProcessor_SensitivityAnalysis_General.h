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

#ifndef OPENSMOKE_POSTPROCESSOR_SENSITIVITYANALYSIS_GENERAL
#define OPENSMOKE_POSTPROCESSOR_SENSITIVITYANALYSIS_GENERAL

#include <sstream>
#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_PostProcessor; 

class OpenSMOKE_PostProcessor_SensitivityAnalysis_General
{
public:

	OpenSMOKE_PostProcessor_SensitivityAnalysis_General() {};
	OpenSMOKE_PostProcessor_SensitivityAnalysis_General(OpenSMOKE_PostProcessor *post_processor_);
	void SetName(const string name);
	string sensitivity_list();
	
	virtual void ReadFromBinaryFile(BzzLoad &fLoad, const int index) {};

	// Get Bars
	void GetAdditionalSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, vector<double> &t, vector<int> &it, vector<string> &names_t, const int index_selection);
	void GetMassFractionSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, const string name, vector<double> &t, vector<int> &it, vector<string> &names_t);

	// Get Profiles
	void GetAdditionalSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t, const int index_selection);
	void GetMassFractionSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, const string name, BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t);

	// Export functions
	void ExportAdditionalNames(vector<string> &_additional_names);
	void ExportSensitivitySpeciesNames(vector<string> &species_names);
	void ExportAvailableYAxisSensitivityCoefficients(vector<string> &y_available);

	// Import functions
	void ImportSelectedAxisSensitivityCoefficients(int x_axis, vector<int> y_axis, BzzMatrix &xAxis, BzzMatrix &yAxis, string &name_x, string &name_y, vector<string> &names_lines);

	BzzVector GetProfileNormalized(const int index_local, const int index_parameter);
	BzzVector GetProfile(const int index_local, const int index_parameter);

	void GetSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, const int index_local, BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t);
	void GetSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, const int index_local, vector<double> &t, vector<int> &it, vector<string> &names_t);

protected:

	OpenSMOKE_PostProcessor *post_processor;

	int N;					// Number of points
	int NC;					// Number of species
	int NR;					// Number of reactions
	int NP;					// Number of parameters
	int NV;					// Number of variables
	int S_NC;				// Number of species (reduced)
	int S_NV;				// Number of variables (reduced)
	string *reactions;		// Reaction strings

	BzzVector parameters;			// Parameters
	BzzMatrix S;					// Sensitivity coefficients
	BzzVector vectorMax;			// Max values
	BzzVector vectorMin;			// Min values
	BzzVectorInt indexTransfer;		// Transfer indices

	BzzMatrix Slocal;
	BzzVector SlocalVector;
	BzzMatrix matrixVariables;

	int	index_temperature;
	int	index_massflowrate;
	int	index_species;
	int	index_velocity;
	int	index_pressurecurvature;

	int	nExtracted;

	string sensitivity_list_;

protected:

	vector<string> additional_names;
	vector<string> species_names;
	vector<string> list_of_y_available_reaction_rates;
	vector<string> list_of_y_labels_reaction_rates;

	void Prepare();
	void MinMaxCalculations();
	void PrepareSensitivitySpeciesNames();

	virtual void PrepareAdditionalNames()			{};

protected:

	string name_object;
	virtual void ErrorMessage(const string message);
	virtual void WarningMessage(const string message);
};

#endif // OPENSMOKE_POSTPROCESSOR_SENSITIVITYANALYSIS_GENERAL
