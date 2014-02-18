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

#include "basic/OpenSMOKE_Utilities.h"
#include "addons/OpenSMOKE_SolidRegression.h"

OpenSMOKE_SolidRegression* ptRegression;
void BzzModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);
void BzzOdeModel_01(BzzVector &y, double t, BzzVector &f);

BzzVector kOde(3);

void OpenSMOKE_SolidRegression::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SolidRegression"		<< endl;
    //cout << "File:   " << name_of_file			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SolidRegression::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SolidRegression"	<< endl;
    //cout << "File:    " << name_of_file			<< endl;
    cout << "Warning: " << message				<< endl;
	cout << endl;
}


void OpenSMOKE_SolidRegression::Setup(const string filename)
{
	ptRegression = this;

	int i;
	string dummy;
	ifstream fInput;
	openInputFileAndControl(fInput, filename);

	while (fInput.eof()==0)
	{
		fInput >> dummy;
		list_of_names_of_files.push_back(dummy);
	}

	fInput.close();

	nCases = list_of_names_of_files.size();
	ChangeDimensions(nCases, &indices);
	ChangeDimensions(nCases, &initialconditions_x);
	ChangeDimensions(nCases, &initialconditions_y);
	ChangeDimensions(nCases, &yMin);

	cout << "List of files: " << endl;
	for(i=1;i<=nCases;i++)
		cout << list_of_names_of_files[i-1] << endl;

	experiments = new OpenSMOKE_SolidExperiment[nCases+1];
	for(i=1;i<=nCases;i++)
		experiments[i].ReadFromFile(list_of_names_of_files[i-1]);

	numTotal = 0;
	for(i=1;i<=nCases;i++)
	{
		indices[i]				= numTotal+1;
		initialconditions_x[i]	= experiments[i].x[1]; 
		initialconditions_y[i]	= experiments[i].y[1]; 
		numTotal  += experiments[i].nPoints;
	}

	ChangeDimensions(numTotal, 1, &YY);

	for(i=1;i<=nCases;i++)
	{
		cout << "-----------------------------------------"		<< endl;
		cout << "Experiment #" << i << endl;
		cout << "-----------------------------------------"		<< endl;
		cout << "  * Index:  " << indices[i] << endl;
		cout << "  * x0:     " << initialconditions_x[i] << endl;
		cout << "  * y0:     " << initialconditions_y[i] << endl;
		cout << endl;
	}
	
}

void OpenSMOKE_SolidRegression::Run()
{
	int numModels	= 1;
	int numX		= 1;
	int numY		= 1;

	BzzMatrix X(numTotal, numX);
	BzzMatrix Y(numTotal, numY);

	int k=1;
	for(int i=1;i<=nCases;i++)
	{
		for(int j=1;j<=experiments[i].nPoints;j++)
		{
			X[k][1]		= experiments[i].x[j]; 
			Y[k++][1]	= experiments[i].y[j]; 
		}
	}
	
	BzzVector b(3, 1., 1., 1.);

	BzzNonLinearRegression nonLinReg(numModels, X, Y, BzzModelOdeRegression);
	nonLinReg.InitializeModel(1,b);
	nonLinReg.LeastSquaresAnalysis();
}

void BzzModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	ptRegression->ModelOdeRegression(model, ex, b, x, y);
}

void OpenSMOKE_SolidRegression::ModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	int i;

	for(i=1;i<=nCases;i++)
		if(ex == indices[i])
		{
			cout << "Initialize " << i << "(" << ex << ")" << endl;
			kOde[1] = b[1];
			kOde[2] = b[2]/experiments[i].temperature;
			kOde[3] = pow(experiments[i].pressureO2, b[3]);
			o.Deinitialize();
			BzzVector y0(1, initialconditions_y[i]);
			o.SetInitialConditions(y0, initialconditions_x[i], BzzOdeModel_01);
			o.SetMinimumConstraints(&yMin);
			//break;
		}

	for(i=1;i<=nCases;i++)
		if(ex == indices[i])
		{
			cout << "Calculating " << i << "(" << ex << ")" << endl;

			for(int j=1;j<=experiments[i].nPoints;j++)
			{
				cout << indices[i]+j-1 << "\t" << experiments[i].x[j] << endl;
				yy = o(experiments[i].x[j]);
				YY.SetRow(indices[i]+j-1,yy); 
			}
			break;
		}
	
	YY.GetRow(ex, &y);
}

void BzzOdeModel_01(BzzVector &y, double t, BzzVector &f)
{
	cout << t << "\t" << y[1] << "\t" << kOde[1] << "\t" << kOde[2] << "\t" << kOde[3] <<endl;
	double r = kOde[1]*exp(-kOde[2])*kOde[3];
	f[1] = -y[1] * r;
}


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

void OpenSMOKE_SolidExperiment::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SolidExperiment"		<< endl;
    cout << "File:   " << name_of_file			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press enter to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SolidExperiment::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SolidExperiment"	<< endl;
    cout << "File:    " << name_of_file			<< endl;
    cout << "Warning: " << message				<< endl;
	cout << endl;
}

void OpenSMOKE_SolidExperiment::ReadFromFile(const string filename)
{
	string dummy;

	name_of_file = filename;
	ifstream fInput;
	openInputFileAndControl(fInput, filename);

	fInput >> dummy;
	if (dummy!="TEMPERATURE")	ErrorMessage("Expected: TEMPERATURE - Found: " + dummy);
	fInput >> temperature;
	fInput >> dummy;
	if (dummy!="K" && dummy!="C")	ErrorMessage("Expected: K || C - Found: " + dummy);
	if (dummy=="C")	temperature += 273.15;

	fInput >> dummy;
	if (dummy!="PRESSURE_O2")	ErrorMessage("Expected: PRESSURE_O2 - Found: " + dummy);
	fInput >> pressureO2;
	fInput >> dummy;
	if (dummy!="atm" && dummy!="bar" && dummy!="Pa")	ErrorMessage("Expected: atm || bar || Pa - Found: " + dummy);
	if (dummy=="bar")	pressureO2 /= 1.01325;
	if (dummy=="Pa")	pressureO2 /= 101325.;

	for(;;)
	{
		fInput >> dummy;
		if (dummy == "//")	break;
		x.Append(atof(dummy.c_str()));

		fInput >> dummy;
		y.Append(atof(dummy.c_str()));
	}
	fInput.close();

	nPoints = x.Size();

	cout << "-----------------------------------------"		<< endl;
	cout << "Experiment:      " << filename					<< endl;
	cout << "-----------------------------------------"		<< endl;
	cout << "  * Temperature: " << temperature	<< " K"		<< endl;
	cout << "  * Pressure O2: " << pressureO2	<< " atm"	<< endl;
	cout << "  * #points:     " << nPoints					<< endl;
	cout << "  * x0:          " << x[1]						<< endl;
	cout << "  * xF:          " << x[nPoints]				<< endl;
	cout << "  * yMean:       " << Mean(y)					<< endl;
	cout << "-----------------------------------------"		<< endl;
	cout << endl;
}

void nonlinearregression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
}