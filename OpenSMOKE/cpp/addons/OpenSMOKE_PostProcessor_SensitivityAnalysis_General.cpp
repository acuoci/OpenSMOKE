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
#include "addons/OpenSMOKE_PostProcessor_SensitivityAnalysis_General.h"

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::ErrorMessage(const string message)
{
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::WarningMessage(const string message)
{
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::SetName(const string name)
{
	name_object = name;
}

string OpenSMOKE_PostProcessor_SensitivityAnalysis_General::sensitivity_list() 
{
	return sensitivity_list_;
}


OpenSMOKE_PostProcessor_SensitivityAnalysis_General::OpenSMOKE_PostProcessor_SensitivityAnalysis_General(OpenSMOKE_PostProcessor *post_processor_)
{
	name_object = "name not assigned";
	post_processor = post_processor_;
	NC = post_processor->NC;
	NR = post_processor->NR;
	
	N		= 0;
	NP		= 0;
	NV		= 0;
	S_NC	= 0;
	S_NV	= 0;
	
	nExtracted					= 20;
	index_velocity				= 0;
	index_massflowrate			= 0;
	index_pressurecurvature		= 0;
	index_temperature			= 0;
	index_species				= 0;
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::MinMaxCalculations()
{
	cout << " * Min-Max calculations... " << endl;
	for(int j=1;j<=S_NV;j++)
		vectorMax[j] = matrixVariables.GetColumn(j).MaxAbs();
	for(int j=1;j<=S_NV;j++)
		vectorMin[j] = matrixVariables.GetColumn(j).MinAbs();
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::ExportAdditionalNames(vector<string> &_additional_names)
{
	_additional_names = additional_names;
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::ExportSensitivitySpeciesNames(vector<string> &_species_names)
{
	_species_names = species_names;
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::PrepareSensitivitySpeciesNames()
{
	species_names.resize(indexTransfer.Size());
	for(int j=1;j<=indexTransfer.Size();j++)
		species_names[j-1] = post_processor->names[indexTransfer[j]];
}


void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::GetSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, const int index_local, vector<double> &t, vector<int> &it, vector<string> &names_t)
{
	BzzVector SVector;
	BzzVectorInt indices;
	BzzMatrix SMatrix;
	GetSensitivityProfiles(iTotal, iLocal, coordinate, index_local, SMatrix, indices, names_t);

	ChangeDimensions(nExtracted, &SVector);
	if (iTotal == true)
	{
		for(int j=1;j<=nExtracted;j++)
		{
			double a = SMatrix.GetColumn(j).Max();
			double b = SMatrix.GetColumn(j).Min();
			SVector[j] = fabs(a) >= fabs(b) ? a : b;
		}
	}
	else
	{
		BzzVector x;
		x = post_processor->get_x();

		int iPoint;
		for(int i=1;i<=N;i++)
			if ( x[i] >= coordinate)	{	iPoint = i; break;}

		for(int j=1;j<=nExtracted;j++)
			SVector[j] =SMatrix[iPoint][j];
	}

	// From BzzMath to STL
	t.resize(indices.Size());
	it.resize(indices.Size());
	for(int i=0;i<=indices.Size()-1;i++)
	{
		t[i]  = SVector[i+1];
		it[i] = indices[i+1];
	}
}

BzzVector OpenSMOKE_PostProcessor_SensitivityAnalysis_General::GetProfileNormalized(const int index_local, const int index_parameter)
{
	BzzVector y = GetProfile(index_local, index_parameter);
	y *= (parameters[index_parameter] / vectorMax[index_local]);
	return y;
}

BzzVector OpenSMOKE_PostProcessor_SensitivityAnalysis_General::GetProfile(const int index_local, const int index_parameter)
{
	BzzVector y(N);
	for(int i=1;i<=N;i++)
	{
		int index_global = (i-1)*S_NV + index_local;
		y[i]  = S[index_global][index_parameter];
	}
	return y;
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::GetSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, const int index_local, BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t)
{
	int i,j;

	for(i=1;i<=N;i++)
	{
		int index_global = (i-1)*S_NV + index_local;
		Slocal.SetRow(i, S.GetRow(index_global));
	}


	double threshold = 1.e-10;

	if (iLocal == false)
	{
		for(i=1;i<=N;i++)
			for(j=1;j<=NP;j++)
			{
				if (vectorMax[index_local] >= threshold)	Slocal[i][j] *= parameters[j] / vectorMax[index_local];
				else										Slocal[i][j] *= 0;
			}	
	}
	
	if (iLocal == true)
	{
		for(i=1;i<=N;i++)
			for(j=1;j<=NP;j++)
			{
				if (vectorMin[index_local] >= threshold)	Slocal[i][j] *= parameters[j] / matrixVariables[i][index_local];
				else																		Slocal[i][j] *= 0;
			}	
	}

	BzzVector yMax(NP);
	BzzVector yMaxWithSign(NP);
	BzzVectorInt iMax(NP);

	if (iTotal == true)
	{
		for(j=1;j<=NP;j++)
		{	
			int iiMax;
			yMax[j] = Slocal.GetColumn(j).MaxAbs(&iiMax); 
			yMaxWithSign[j] = Slocal.GetColumn(j)[iiMax];
		}
	}
	else
	{
		BzzVector x;
		x = post_processor->get_x();

		int iPoint;
		for(i=1;i<=N;i++)
			if (x[i]>=coordinate)	{	iPoint = i; break;}
		for(j=1;j<=NP;j++)
		{
			yMax[j] = fabs(Slocal[iPoint][j]);
			yMaxWithSign[j] = Slocal[iPoint][j];
		}
	}

	Sort(&yMax, &iMax);
	Reorder(&yMaxWithSign, iMax);
	Reverse(&iMax);
	Reverse(&yMax);
	Reverse(&yMaxWithSign);

	ChangeDimensions(N, nExtracted, &SMatrix);
	ChangeDimensions(nExtracted, &indices);
	for(j=1;j<=nExtracted;j++)
		SMatrix.SetColumn(j, Slocal.GetColumn(iMax[j]));
	
	for(j=1;j<=nExtracted;j++)
		indices[j] = iMax[j];

	// From BzzMath to STL
	names_t.resize(indices.Size());
	for(j=0;j<=indices.Size()-1;j++)
		names_t[j] =reactions[indices[j+1]];

	
	int nVideo = (NP <= 200) ? NP : 200;
//	cout.setf(ios::scientific);
	sensitivity_list_  = "-----------------------------------------------------------------------------------------------\n";
	sensitivity_list_ += " Sensitivity list (first 200 parameters)                                                       \n";
	sensitivity_list_ += "-----------------------------------------------------------------------------------------------\n";
	for(j=1;j<=nVideo;j++)
	{
		stringstream str_j; str_j << j;
		stringstream str_iMax; str_iMax << iMax[j];
		stringstream str_yMax; str_yMax << yMaxWithSign[j];
		sensitivity_list_ += " " + str_j.str() + "\t" + str_iMax.str() + "\t" + str_yMax.str() + "\t" + reactions[iMax[j]] + "\n";
	}
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::GetAdditionalSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, vector<double> &t, vector<int> &it, vector<string> &names_t, const int index_selection)
{
	GetSensitivityBars(iTotal, iLocal, coordinate, index_selection, t, it, names_t);
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::GetAdditionalSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t, const int index_selection)
{
	GetSensitivityProfiles(iTotal, iLocal, coordinate, index_selection, SMatrix, indices, names_t);
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::GetMassFractionSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, const string name, vector<double> &t, vector<int> &it, vector<string> &names_t)
{
	int index = 0;
	for(int j=1;j<=S_NC;j++)
	{
		if (post_processor->names[indexTransfer[j]] == name)
		{
			index = j;
			break;
		}
	}

	int index_local = index_species + index-1;

	GetSensitivityBars(iTotal, iLocal, coordinate, index_local, t, it, names_t);
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::GetMassFractionSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, const string name, BzzMatrix &SMatrix, BzzVectorInt &indices, vector<string> &names_t)
{
	cout << "Transfer" << endl;
	cout << S_NC << endl;
	cout << index_species << endl;
	cout << indexTransfer.Size() << endl;
	for(int j=1;j<=indexTransfer.Size();j++)
		cout << "a" << j << " " << indexTransfer[j] << endl;
	getchar();

	int index = 0;
	for(int j=1;j<=S_NC;j++)
		if (post_processor->names[indexTransfer[j]] == name)
		{
			index = j;
			break;
		}

	int index_local = index_species + index-1;

	GetSensitivityProfiles(iTotal, iLocal, coordinate, index_local, SMatrix, indices, names_t);
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::Prepare()
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
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::ExportAvailableYAxisSensitivityCoefficients(vector<string> &y_available)
{
	y_available = list_of_y_available_reaction_rates;
}

void OpenSMOKE_PostProcessor_SensitivityAnalysis_General::ImportSelectedAxisSensitivityCoefficients(int x_axis, vector<int> y_axis, BzzMatrix &xAxis, BzzMatrix &yAxis, string &name_x, string &name_y, vector<string> &names_lines)
{
	ChangeDimensions(y_axis.size(), post_processor->get_x().Size(), &xAxis);
	ChangeDimensions(y_axis.size(), post_processor->get_x().Size(), &yAxis);

	name_x = "x coordinate";
	name_y = "sensitivity coefficient";

	names_lines.resize(0);
	for(int j=1;j<=int(y_axis.size());j++)
	{
		names_lines.push_back(list_of_y_labels_reaction_rates[y_axis[j-1]]);

		if (post_processor->focus_sensitivity_additional == true)
		{
			BzzVector aux = GetProfileNormalized(post_processor->focus_sensitivity_index+post_processor->startSensitivityAdditional, y_axis[j-1]+1);
			xAxis.SetRow(j,post_processor->get_x());
			yAxis.SetRow(j,aux);
		}
		else
		{
			BzzVector aux = GetProfileNormalized(index_species+post_processor->focus_sensitivity_index-1, y_axis[j-1]+1);
			xAxis.SetRow(j,post_processor->get_x());
			yAxis.SetRow(j,aux);
		}		
	}
}


