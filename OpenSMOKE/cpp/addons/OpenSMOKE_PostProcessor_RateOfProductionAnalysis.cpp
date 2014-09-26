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
#include <iomanip>
#include <sstream>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "addons/OpenSMOKE_PostProcessor.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "addons/OpenSMOKE_RateOfProductionAnalysis.h"
#include "addons/OpenSMOKE_PostProcessor_RateOfProductionAnalysis.h"

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PostProcessor_RateOfProductionAnalysis"	<< endl;
    cout << "Object: " << name_object									<< endl;
    cout << "Error:  " << message										<< endl;
    cout << "Press a key to continue... "								<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PostProcessor_RateOfProductionAnalysis"	<< endl;
    cout << "Object: "		<< name_object								<< endl;
    cout << "Warning:  "	<< message									<< endl;
    cout << "Press a key to continue... "								<< endl;
    getchar();
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::SetName(const std::string name)
{
	name_object = name;
}

OpenSMOKE_PostProcessor_RateOfProductionAnalysis::OpenSMOKE_PostProcessor_RateOfProductionAnalysis(OpenSMOKE_PostProcessor *post_processor_)
{
	name_object = "name not assigned";
	post_processor = post_processor_;
	NC = post_processor->NC;
	NR = post_processor->NR;
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::ReadFromBinaryFile(BzzLoad &fLoad)
{
	char dummy[Constants::NAME_SIZE];
	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	std::string version = dummy;
	if (version != "V20100417")
		ErrorMessage("This version post processing file is not supported: " + version);
	cout << "Version: " << version << endl;

	// Indices of species
	CheckInBinaryFile(fLoad, "INDICES");
	fLoad >> indices;	

	// Reaction rates
	CheckInBinaryFile(fLoad, "REACTIONRATES");
	fLoad >> reactionRates;

	// Preparing data
	ropa = new OpenSMOKE_RateOfProductionAnalysis();
	BzzVector aux_x = post_processor->get_x();
	ropa->Initialize(post_processor->mix, indices);
	ropa->SetNumberOfPoints(aux_x.Size());
	ropa->Run(reactionRates, aux_x);

	// Formation rates
	{
		ChangeDimensions(post_processor->get_x().Size(), NC, &formationRates);
		BzzVector RVector(NC);
		BzzVector C(NC);
		for(int j=1;j<=post_processor->get_x().Size();j++)
		{
			double T    = post_processor->get_T(j);
			double P_Pa = post_processor->get_P_Pa(j);
			double Ctot = post_processor->get_Ctot(j);
			C = post_processor->get_C(j);
			post_processor->mix->ComputeKineticParameters( T, log(T), 1./T, P_Pa);
			post_processor->mix->ComputeFromConcentrations( T, C, Ctot, &RVector);	// [kmol/m3/s]
			formationRates.SetRow(j, RVector);										// [kmol/m3/s]
		}
	}
	Prepare();
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::Prepare()
{
	list_of_y_available_reaction_rates.resize(0);
	list_of_y_labels_reaction_rates.resize(0);
	for(int j=1;j<=NR;j++)
	{
		stringstream number;
		number << j;
		list_of_y_available_reaction_rates.push_back(number.str() + " - " + post_processor->reactions[j]);
		list_of_y_labels_reaction_rates.push_back(number.str());
	}

	list_of_y_available_formation_rates.resize(0);
	list_of_y_labels_formation_rates.resize(0);
	for(int j=1;j<=NC;j++)
	{
		list_of_y_available_formation_rates.push_back(post_processor->names[j]);
		list_of_y_labels_formation_rates.push_back(post_processor->names[j]);
	}
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::ExportAvailableYAxisReactionRates(vector<string> &y_available)
{
	y_available = list_of_y_available_reaction_rates;
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::ExportAvailableYAxisFormationRates(vector<string> &y_available)
{
	y_available = list_of_y_available_formation_rates;
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::ImportSelectedAxisReactionRates(int x_axis, vector<int> y_axis, BzzMatrix &xAxis, BzzMatrix &yAxis, std::string &name_x, std::string &name_y, vector<string> &names_lines)
{
	ChangeDimensions(y_axis.size(), post_processor->get_x().Size(), &xAxis);
	ChangeDimensions(y_axis.size(), post_processor->get_x().Size(), &yAxis);

	name_x = "x coordinate";
	name_y = "reaction rate [kmol/m3/s]";

	names_lines.resize(0);
	for(int j=1;j<=int(y_axis.size());j++)
	{
		names_lines.push_back(list_of_y_labels_reaction_rates[y_axis[j-1]]);

		BzzVector aux = reactionRates.GetColumn(y_axis[j-1]+1);
		xAxis.SetRow(j,post_processor->get_x());
		yAxis.SetRow(j,aux);
	}
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::ImportSelectedAxisFormationRates(int x_axis, vector<int> y_axis, BzzMatrix &xAxis, BzzMatrix &yAxis, std::string &name_x, std::string &name_y, vector<string> &names_lines)
{
	ChangeDimensions(y_axis.size(), post_processor->get_x().Size(), &xAxis);
	ChangeDimensions(y_axis.size(), post_processor->get_x().Size(), &yAxis);

	name_x = "x coordinate";
	name_y = "formation rate [kmol/m3/s]";

	names_lines.resize(0);
	for(int j=1;j<=int(y_axis.size());j++)
	{
		names_lines.push_back(list_of_y_labels_formation_rates[y_axis[j-1]]);

		BzzVector aux = formationRates.GetColumn(y_axis[j-1]+1);
		xAxis.SetRow(j,post_processor->get_x());
		yAxis.SetRow(j,aux);
	}
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::CallIntegralAnalysis()
{
	BzzVector aux_x = post_processor->get_x();
	ropa->Run(reactionRates, aux_x);
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::CallLocalAnalysis(const double local_coordinate)
{
	BzzVector aux_x = post_processor->get_x();
	ropa->Run(reactionRates, aux_x, local_coordinate);
}
	
void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::CallRegionAnalysis(const double xA, const double xB)
{
	BzzVector aux_x = post_processor->get_x();
	ropa->Run(reactionRates, aux_x, xA, xB);
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::MostImportantReactions(const int index, vector<double> &p, vector<double> &d, vector<int> &ip, vector<int> &id, 
																							   vector<string> &names_p, vector<string> &names_d)
{
	ropa->MostImportantReactions(index, p, d, ip, id, names_p, names_d);
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::MostImportantReactions(const int index, vector<double> &t, vector<int> &it, vector<string> &names_t)
{
	ropa->MostImportantReactions(index, t, it, names_t);
}

void OpenSMOKE_PostProcessor_RateOfProductionAnalysis::UnimportantReactions(const double eps_threshold, stringstream &string_out)
{
	ropa->UnimportantReactions(string_out, eps_threshold);
}
