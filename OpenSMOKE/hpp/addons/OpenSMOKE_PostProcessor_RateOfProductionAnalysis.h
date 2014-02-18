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

#ifndef OPENSMOKE_POSTPROCESSOR_RATEOFPRODUCTIONANALYSIS
#define OPENSMOKE_POSTPROCESSOR_RATEOFPRODUCTIONANALYSIS

#include <sstream>
#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_PostProcessor; 
class OpenSMOKE_RateOfProductionAnalysis;

class OpenSMOKE_PostProcessor_RateOfProductionAnalysis
{
public:

	OpenSMOKE_PostProcessor_RateOfProductionAnalysis(OpenSMOKE_PostProcessor *post_processor_);
	void SetName(const string name);
	void ReadFromBinaryFile(BzzLoad &fLoad);

	void CallIntegralAnalysis();
	void CallLocalAnalysis(const double local_coordinate);
	void CallRegionAnalysis(const double xA, const double xB);

	void MostImportantReactions(const int index, vector<double> &p, vector<double> &d, vector<int> &ip, vector<int> &id, 
								vector<string> &names_p, vector<string> &names_d);
	void MostImportantReactions(const int index, vector<double> &t, vector<int> &it, vector<string> &names_t);
	void UnimportantReactions(const double eps_threshold_percentage, stringstream &string_out);

	void Prepare();
	void ExportAvailableYAxisReactionRates(vector<string> &y_available);
	void ImportSelectedAxisReactionRates(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, string &name_x, string &name_y, vector<string> &names_lines);
	void ExportAvailableYAxisFormationRates(vector<string> &y_available);
	void ImportSelectedAxisFormationRates(int x_axis, vector<int> y_axis, BzzMatrix &x, BzzMatrix &y, string &name_x, string &name_y, vector<string> &names_lines);

private:

	OpenSMOKE_PostProcessor					*post_processor;
	OpenSMOKE_RateOfProductionAnalysis		*ropa;

	BzzMatrix		reactionRates;
	BzzMatrix		formationRates;
	BzzVectorInt	indices;

	vector<string> list_of_y_available_reaction_rates;
	vector<string> list_of_y_labels_reaction_rates;
	vector<string> list_of_y_available_formation_rates;
	vector<string> list_of_y_labels_formation_rates;

	int NC;					// Number of species
	int NR;					// Number of reactions

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_POSTPROCESSOR_RATEOFPRODUCTIONANALYSIS
