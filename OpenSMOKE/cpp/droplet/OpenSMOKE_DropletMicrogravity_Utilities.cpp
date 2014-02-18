/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci								   *
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
#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"

void InterpolateBetweenGrids(const BzzVector &xold, BzzVector &xnew, BzzVector &yold, BzzVector &ynew)
{
	// Checking boundaries
	BzzVector xStar;
	xStar = xold;
	if (xold[1] != xnew[1])
	{
		cout << " * Warning: non conformal domain..." << endl;
		cout << " * Backup domain: " << xold[1]*1000. << " " << xold[xold.Size()]*1000. << endl;
		cout << " * New domain:    " << xnew[1]*1000. << " " << xnew[xnew.Size()]*1000. << endl;
		
		double ratio = (xold[xold.Size()]-xnew[1])/(xold[xold.Size()]-xold[1]);
		xStar[1] = xnew[1];
		for(int j=2;j<=xold.Size();j++)
			xStar[j] = xStar[j-1] + (xold[j]-xold[j-1])*ratio;
		if ( fabs(xStar[xStar.Size()] - xold[xold.Size()])/xold[xold.Size()] > 1.e-6 )
		{
			cout << "Domain was not correctly conformed..." << endl;
			cout << "Press enter to exit..." << endl;
			getchar();
			exit(-1);
		}
	}

	LinearInterpolation linear;
	linear(xStar, yold);
	
	double xStar_max = xStar[xStar.Size()];
	for(int i=1;i<=xnew.Size();i++)
	{
		if (xnew[i] > xStar_max)
			ynew[i] = yold[yold.Size()];
		else
			ynew[i] = linear(xnew[i]);
	}
}

void FindMaximumUsingSecondOrderPolynomial(BzzVector &x, BzzVector &y, double &xMax, double &yMax)
{
	int iMax;
	y.Max(&iMax);
	
	if ((iMax == 1) || (iMax == x.Size()))
	{
		xMax = 0.;
		yMax = 0.;
		return;
	}	

	double X1 = x[iMax-1]; double X2 = x[iMax]; double X3 = x[iMax+1];
	double Y1 = y[iMax-1]; double Y2 = y[iMax]; double Y3 = y[iMax+1];
	double ratiox = (X1+X2+X3)/3.;
	double ratioy = (Y1+Y2+Y3)/3.;

	X1 /= ratiox; X2 /= ratiox; X3 /= ratiox;
	Y1 /= ratioy; Y2 /= ratioy; Y3 /= ratioy;

	if ( (fabs(Y2-Y1)/Y2< 1e-6) || (fabs(Y2-Y3)/Y2< 1e-6) )
	{
		xMax = 0.;
		yMax = 0.;
		return;
	}

	BzzMatrix A(3,3);
	A[1][1] = X1*X1; A[1][2] = X1; A[1][3] = 1.; 
	A[2][1] = X2*X2; A[2][2] = X2; A[2][3] = 1.; 
	A[3][1] = X3*X3; A[3][2] = X3; A[3][3] = 1.; 

	BzzVector b(3);
	b[1] = Y1; b[2] = Y2; b[3] = Y3;

	BzzFactorizedGauss Alfa(A);
	Solve(&Alfa, &b);

	xMax = -b[2]/2./b[1] ;
	yMax = b[1]*xMax*xMax + b[2]*xMax + b[3];
	xMax *= ratiox;
	yMax *= ratioy;
}

double FindWidth(BzzVector &x, BzzVector &y)
{
	int iMax;
	double yLeft = y[1];
	double yMax = y.Max(&iMax);
	double yHalf = yLeft + (yMax-yLeft)/2.;

	int iLeft = 0;
	for(int j=1;j<=iMax;j++)
		if (y[j] > yHalf)	
		{
			iLeft = j;
			break;
		}

	int iRight = 0;
	for(int j=iMax;j<=x.Size();j++)
		if (y[j] < yHalf)	
		{
			iRight = j;
			break;
		}

	if (iLeft != 0 && iRight != 0)
	{
		double xLeft, xRight;
		{
			double a = (y[iLeft-1] - y[iLeft])/(x[iLeft-1] - x[iLeft]);
			double b = y[iLeft]-a*x[iLeft];
			xLeft = (yHalf-b)/a;
		}
		{
			double a = (y[iRight-1] - y[iRight])/(x[iRight-1] - x[iRight]);
			double b = y[iRight]-a*x[iRight];
			xRight = (yHalf-b)/a;
		}
		return (xRight - xLeft);
	}
	else
		return 0.;
}