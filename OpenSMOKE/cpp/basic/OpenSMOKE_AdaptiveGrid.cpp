/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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

#include "OpenSMOKE.hpp"
#include "basic/OpenSMOKE_AdaptiveGrid.h"


void OpenSMOKE_AdaptiveGrid::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_AdaptiveGrid"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_AdaptiveGrid::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_AdaptiveGrid"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

void OpenSMOKE_AdaptiveGrid::SetName(const std::string name)
{
	name_object = name;
}


OpenSMOKE_AdaptiveGrid::OpenSMOKE_AdaptiveGrid()
{
	model = ADAPTIVE_CHEMKIN_POW;
	alfa  = 1.00;
	beta  = 0.60;
	gamma = 0.15;
	name_object = "[not assigned]";
}

void OpenSMOKE_AdaptiveGrid::SetModel(adaptive_grid_model _model)
{
	model = _model;
	
	if (model == ADAPTIVE_GRADIENT)
	{
		alfa  = 1.00;
		beta  = 0.;
		gamma = 0.;
	}
	else if (model == ADAPTIVE_CHEMKIN)
	{
		alfa  = 1.00;	// 0.1 0.001 0.
		beta  = 0.60;
		gamma = 0.;
	}
	else if (model == ADAPTIVE_CHEMKIN_SQRT)
	{
		alfa  = 0.60;	// 0.60 0.05 0.
		beta  = 0.05;
		gamma = 0.;
	}
	else if (model == ADAPTIVE_CHEMKIN_POW)
	{
		alfa  = 1.00;	// 1. 0.001 0.25
		beta  = 0.60;
		gamma = 0.15;
	}
	else if (model == ADAPTIVE_EISEMAN)
	{
		alfa  = 0.10;	// 0.1 0.01 0.
		beta  = 0.01;
		gamma = 0.;
	}
}

void OpenSMOKE_AdaptiveGrid::SetConstants(const double _alfa, const double _beta, const double _gamma)
{
	alfa  = _alfa;
	beta  = _beta;
	gamma = _gamma;
}

void OpenSMOKE_AdaptiveGrid::SetAlfa(const double _alfa)
{
	alfa  = _alfa;
}

void OpenSMOKE_AdaptiveGrid::SetBeta(const double _beta)
{
	beta  = _beta;
}

void OpenSMOKE_AdaptiveGrid::SetGamma(const double _gamma)
{
	gamma  = _gamma;
}

void OpenSMOKE_AdaptiveGrid::Setup(BzzVector &_xOld, BzzVector &_yOld)
{
	N = _xOld.Size();
	
	ChangeDimensions(N, &d2x_over_dcsi2);
	ChangeDimensions(N, &dx_over_dcsi);
	ChangeDimensions(N, &d2y_over_dx2);
	ChangeDimensions(N, &dy_over_dx);
	ChangeDimensions(N, &W);
	ChangeDimensions(N, &dW_over_dcsi);
	ChangeDimensions(N, &y);

	xOld	= _xOld;
	yOld	= _yOld;

	grid.Construct(xOld);

	UpdateProfileDerivatives(xOld);
	UpdateWeights();
}

void OpenSMOKE_AdaptiveGrid::UpdateWeights()
{
	int j;

	if (model == ADAPTIVE_GRADIENT)
		for(j=1;j<=N;j++)	W[j] = sqrt(1. + alfa*dy_over_dx[j]*dy_over_dx[j]);

	else if (model == ADAPTIVE_CHEMKIN)
		for(j=1;j<=N;j++)	W[j] = 1. + alfa*fabs(dy_over_dx[j]) + beta*fabs(d2y_over_dx2[j]);
	
	else if (model == ADAPTIVE_CHEMKIN_SQRT)
		for(j=1;j<=N;j++)	W[j] = sqrt(1.+alfa*fabs(dy_over_dx[j])+beta*fabs(d2y_over_dx2[j]));
	
	else if (model == ADAPTIVE_CHEMKIN_POW)	
		for(j=1;j<=N;j++)	W[j] = pow(1.+alfa*fabs(dy_over_dx[j])+beta*fabs(d2y_over_dx2[j]), gamma);
	
	else if (model == ADAPTIVE_EISEMAN)	
		for(j=1;j<=N;j++)	
		{
			double kappa = d2y_over_dx2[j] / pow(1.+dy_over_dx[j]*dy_over_dx[j], 1.5);
			W[j] = sqrt(1.+alfa*alfa*dy_over_dx[j]*dy_over_dx[j]) * (1.+beta*beta*fabs(kappa));
		}

	UpdateWeightDerivatives();
}

void OpenSMOKE_AdaptiveGrid::UpdateWeightDerivatives()
{
	dW_over_dcsi[1] = W[2]-W[1];
	dW_over_dcsi[N] = W[N]-W[N-1];

	for(int j=2;j<=N-1;j++)	dW_over_dcsi[j] = (W[j+1]-W[j-1])/2.;
}

void OpenSMOKE_AdaptiveGrid::UpdateProfileDerivatives(BzzVector &xNew)
{
	bool iFind;
	int  jOld;

	y[1] = yOld[1];
	y[N] = yOld[N];

	jOld = 1;
	for(int j=2;j<=N-1;j++)
	{
		iFind = false;
		for(int k=jOld;k<=N-1;k++)
		{
			if (iFind == true) break;

			if (xNew[j] > xOld[k] && xNew[j] <= xOld[k+1])
			{
				iFind = true;
				jOld = 1;
				double m = (xNew[j]-xOld[k])/(xOld[k+1]-xOld[k]);
				y[j]  = yOld[k]+m*(yOld[k+1]-yOld[k]);
			}
		}
		if (iFind == false)
			ErrorMessage("Point outside the boundaries!");
	}

	grid.FirstDerivative('C', dummy, y, dy_over_dx);
	grid.SecondDerivative(y, d2y_over_dx2);
	d2y_over_dx2[1] = d2y_over_dx2[2];
	d2y_over_dx2[N]  = d2y_over_dx2[N-1];
}

void OpenSMOKE_AdaptiveGrid::Ode(BzzVector &x, double t, BzzVector &f)
{
	int j;

	for(j=2;j<=N-1;j++)
		d2x_over_dcsi2[j] = (x[j+1]-2.*x[j]+x[j-1]);

	for(j=2;j<=N-1;j++)
		dx_over_dcsi[j] = (x[j+1]-x[j-1])/2.;

	UpdateProfileDerivatives(x);
	UpdateWeights();

	f[1] = 0.;
	for(j=2;j<=N-1;j++)
		f[j] = W[j]*d2x_over_dcsi2[j] + dW_over_dcsi[j]*dx_over_dcsi[j];
	f[N] = 0.;
}

void OpenSMOKE_AdaptiveGrid::PrintOnFile(BzzVector &x)
{
	ofstream fOut;
	openOutputFileAndControl(fOut, "GridProfile.out");
	fOut.setf(ios::scientific);
	for(int j=1;j<=N;j++)
		fOut	<< xOld[j]			<< "\t" << yOld[j]				<< "\t"	
				<< x[j]				<< "\t" << y[j]					<< "\t"
				<< dy_over_dx[j]	<< "\t" << d2y_over_dx2[j]		<< endl;
	fOut.close();
}