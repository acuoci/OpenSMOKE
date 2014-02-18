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

#if !defined(FLAMELET_GROUP)
#define FLAMELET_GROUP

#include "basic/OpenSMOKE_Utilities.h"
#include "flamelet.h"
#include "Utilities.h"

class OpenSMOKE_ReactingGas;

class flamelet_group
{
public:
	int nFlamelets;
	flamelet *flames;
	void initialize();
	void append(flamelet newflame);
	void reorder();
	void lock();

	void build_library(OpenSMOKE_ReactingGas &mix, const string option, const string file_name, int NVariance, double alfa);
	void write_progress_variables(OpenSMOKE_ReactingGas &mix, const int nProgressVariable, vector<string>& progress_variable_name, BzzVectorInt* progress_variable_index, BzzVector* progress_variable_value);
	void print_on_file(char *fileName);
	void print_on_file_FLUENT(string fileName, OpenSMOKE_ReactingGas &mix);
	void print_on_file_look_up_table(string folderName, OpenSMOKE_ReactingGas &mix);

	void apply_source_betaPDF(int nSources);
	void read_flamelet_library(OpenSMOKE_ReactingGas &mix, const string fileName);

	void soot_post_processing(OpenSMOKE_ReactingGas &mix);
	void apply_soot();

	void apply_source_betaPDF_SootProperties(OpenSMOKE_ReactingGas &mix, int nSources, BzzVector *m0, BzzVector *fv);
	void apply_source_betaPDF(	OpenSMOKE_ReactingGas &mix, bool iPerfectlyCorrelatedApproach, int nSources, 	BzzVector *m0_N, 	BzzVector *fv_N,
								nucleation_models nucleation_model, growth_models growth_model, aggregation_models aggregation_model, oxidation_models oxidation_model);


private:
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(FLAMELET_GROUP)
