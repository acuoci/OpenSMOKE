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
#include "addons/OpenSMOKE_PostProcessor.h"
#include "addons/OpenSMOKE_PostProcessor_Batch.h"

void OpenSMOKE_PostProcessor_Batch::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PostProcessor_Batch"		<< endl;
    cout << "Object: " << name_object					<< endl;
    cout << "Error:  " << message						<< endl;
    cout << "Press a key to continue... "				<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_PostProcessor_Batch::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_PostProcessor_Batch"	<< endl;
    cout << "Object: "		<< name_object				<< endl;
    cout << "Warning:  "	<< message					<< endl;
    cout << "Press a key to continue... "				<< endl;
    getchar();
}

OpenSMOKE_PostProcessor_Batch::OpenSMOKE_PostProcessor_Batch(OpenSMOKE_PostProcessor *post_processor_) 
: OpenSMOKE_PostProcessor_General(post_processor_)
{
}

void OpenSMOKE_PostProcessor_Batch::ReadFromBinaryFile(BzzLoad &fLoad)
{
	char dummy[Constants::NAME_SIZE];
	fLoad.fileLoad.read((char*) dummy, sizeof(dummy));
	string version = dummy;
	if (version != "V20100417")
		ErrorMessage("This version post processing file is not supported: " + version);
	cout << "Version: " << version << endl;

	// Grid
	cout << "Reading TIME..." << endl;
	CheckInBinaryFile(fLoad, "TIME");
	fLoad >> x;
	
	// Temperature
	CheckInBinaryFile(fLoad, "TEMPERATURE");
	fLoad >> T;

	// Pressure
	CheckInBinaryFile(fLoad, "PRESSURE");
	fLoad >> P_Pa;

	// Molecular weight
	CheckInBinaryFile(fLoad, "MW");
	fLoad >> MW;

	// Indices
	CheckInBinaryFile(fLoad, "INDICES");
	fLoad >> indices;
	X		= new BzzVector[NC+1];
	omega	= new BzzVector[NC+1];
	C		= new BzzVector[NC+1];

	// Mole Fractions
	cout << "Reading MOLEFRACTIONS..." << endl;
	CheckInBinaryFile(fLoad, "MOLEFRACTIONS");
	for(int j=1;j<=NC;j++)
		fLoad >> X[j];

	// Mass fractions
	cout << "Calculating MASSFRACTIONS..." << endl;
	for(int j=1;j<=NC;j++)
		ChangeDimensions(x.Size(), &omega[j]);
	for(int k=1;k<=x.Size();k++)
		for(int j=1;j<=NC;j++)
			omega[j][k] = X[j][k] * post_processor->M[j] / MW[k];

	// Concentrations
	cout << "Calculating CONCENTRATIONS..." << endl;
	ChangeDimensions(x.Size(), &Ctot);
	ChangeDimensions(x.Size(), &rho);
	for(int j=1;j<=NC;j++)
		ChangeDimensions(x.Size(), &C[j]);
	for(int k=1;k<=x.Size();k++)
	{
		Ctot[k] = P_Pa[k] / (Constants::R_J_kmol*T[k]);		// [kmol/m3]
		rho[k]  = Ctot[k] * MW[k];							// [kg/m3]
		for(int j=1;j<=NC;j++)
			C[j][k] = Ctot[k] * X[j][k];					// [kmol/m3]
	}

	cout << "BATCH succesfully loaded..." << endl;

	Prepare();
}

void OpenSMOKE_PostProcessor_Batch::Prepare()
{
	list_of_x_available.resize(0);
	list_of_x_available.push_back("time");

	list_of_y_available.resize(0);
	list_of_y_available.push_back("temperature");
	list_of_y_available.push_back("concentration");
	list_of_y_available.push_back("density");
	list_of_y_available.push_back("mol. weight");
	list_of_y_available.push_back("pressure");
	for(int j=1;j<=indices.Size();j++)
		list_of_y_available.push_back(post_processor->names[indices[j]]);


	list_of_x_labels.resize(0);
	list_of_x_labels.push_back("time [s]");

	list_of_y_labels.resize(0);
	list_of_y_labels.push_back("temperature [K]");
	list_of_y_labels.push_back("concent. [kmol/m3]");
	list_of_y_labels.push_back("density [kg/m3]");
	list_of_y_labels.push_back("mol. weight [kg/kmol]");
	list_of_y_labels.push_back("pressure [Pa]");
	for(int j=1;j<=indices.Size();j++)
		list_of_y_labels.push_back(post_processor->names[indices[j]]);

	y_start_species = 5;	// from 0
}

void OpenSMOKE_PostProcessor_Batch::ExportAvailableXAxis(vector<string> &x_available)
{
	x_available = list_of_x_available;
}	

void OpenSMOKE_PostProcessor_Batch::ExportAvailableYAxis(vector<string> &y_available)
{
	y_available = list_of_y_available;
}

void OpenSMOKE_PostProcessor_Batch::ImportSelectedAxis(int x_axis, vector<int> y_axis, BzzMatrix &xAxis, BzzMatrix &yAxis, string &name_x, string &name_y, vector<string> &names_lines)
{
	ChangeDimensions(y_axis.size(), x.Size(), &xAxis);
	ChangeDimensions(y_axis.size(), x.Size(), &yAxis);

	name_x = list_of_x_labels.at(x_axis);
	name_y = list_of_y_labels.at(y_axis[0]);

	if (y_axis.size() > 1)
	{
		name_y = "mole fraction";
		for(int j=1;j<=int(y_axis.size());j++)
			if (y_axis[j-1] < y_start_species)
			{
				name_y = "mixed";
				break;
			}
	}

	int iMassFractions = false;
	names_lines.resize(0);
	for(int j=1;j<=int(y_axis.size());j++)
	{
		names_lines.push_back(list_of_y_available[y_axis[j-1]]);

		if (x_axis == 0)			// x: contact time			
			xAxis.SetRow(j,x);

		if (x_axis == 1)			// x: space		
			xAxis.SetRow(j,csi);

		if (y_axis[j-1] == 0)		// y: temperature
			yAxis.SetRow(j,T);

		else if (y_axis[j-1] == 1)		// y: concentration
			yAxis.SetRow(j,Ctot);	

		else if (y_axis[j-1] == 2)		// y: density
			yAxis.SetRow(j,rho);	

		else if (y_axis[j-1] == 3)		// y: molecular weight
			yAxis.SetRow(j,MW);		
	
		else if (y_axis[j-1] == 4)		// y: pressure
			yAxis.SetRow(j,P_Pa);	

		else						// y: species mole fractions
		{
			if (iMassFractions == false)	yAxis.SetRow(j,X[indices[y_axis[j-1]-(y_start_species-1)]]);
			else							yAxis.SetRow(j,omega[indices[y_axis[j-1]-(y_start_species-1)]]);
		}
	}
}