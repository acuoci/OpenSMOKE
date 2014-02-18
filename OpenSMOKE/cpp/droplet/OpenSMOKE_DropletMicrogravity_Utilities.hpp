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

#if !defined(OPENSMOKE_DROPLETMICROGAVITY_UTILITIES_H)
#define OPENSMOKE_DROPLETMICROGAVITY_UTILITIES_H

#include "BzzMath.hpp"
#include "engine/OpenSMOKE_ReactingGas.h"

double AbsorptionCoefficient_H2O(const double T)
{
	double uT = 1000./T;
	return  -0.23093 +uT*(-1.1239+uT*(9.4153 +uT*(-2.9988 +uT*( 0.51382 + uT*(-1.8684e-5)))));
}

double AbsorptionCoefficient_CO2(const double T)
{
	double uT = 1000./T;
	return  18.741  +uT*(-121.31+uT*(273.5  +uT*(-194.05 +uT*( 56.31   + uT*(-5.8169)))));
}

double AbsorptionCoefficient_CO(const double T)
{
	if( T < 750. )	
		return 4.7869+T*(-0.06953 + T*(2.95775e-4 + T*(-4.25732e-7 + T*2.02894e-10)));
	else
		return 10.09+T*(-0.01183 + T*(4.7753e-6 + T*(-5.87209e-10 + T*-2.5334e-14)));
}

double AbsorptionCoefficient_CH4(const double T)
{
	return 6.6334 +T*(- 0.0035686+T*(1.6682e-08+T*(2.5611e-10-2.6558e-14*T)));
}

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

void SparkProfileComplete(BzzVector &x, BzzVector &T, vector<double> &sparkRatio, const double Tinternal, const double Tpeak, const double Tenvironment)
{
	int N = x.Size();
	double radius = x[1]/2.;

	int iA = 0;
	int iB = 0;
	int iC = 0;
	int iD = 0;

	for(int i=1;i<=N;i++)
		if ( x[i]/radius >= sparkRatio[0] )
		{
			iA = i;
			break;
		}
	for(int i=1;i<=N;i++)
		if ( x[i]/radius >= sparkRatio[1] )
		{
			iB = i;
			break;
		}
	for(int i=1;i<=N;i++)
		if ( x[i]/radius >= sparkRatio[2] )
		{
			iC = i;
			break;
		}
	for(int i=1;i<=N;i++)
		if ( x[i]/radius >= sparkRatio[3] )
		{
			iD = i;
			break;
		}

	// Internal zone
	for(int i=2;i<iA;i++)
		T[i] = Tinternal;

	// Ramp (I)
	for(int i=iA;i<iB;i++)
		T[i] = Tinternal + (Tpeak-Tinternal)/(x[iB]-x[iA])*(x[i]-x[iA]);

	// Plateau
	for(int i=iB;i<iC;i++)
		T[i] = Tpeak;

	// Ramp (II)
	for(int i=iC;i<iD;i++)
		T[i] = Tpeak + (Tenvironment-Tpeak)/(x[iD]-x[iC])*(x[i]-x[iC]);

	// External
	for(int i=iD;i<=N;i++)
		T[i] = Tenvironment;

	if ( (iA==iB) || (iB==iC) || (iC==iD))
	{
		cout << "Wrong definition of flame spark: [" << iA << " " << iB << " " << iC << " " << iD << "]" << endl;
		cout << "Press enter to exit..." << endl;
		getchar();
		exit(-1);
	}
	else
	{
		cout << "Spark: [" << iA << " " << iB << " " << iC << " " << iD << "]" << endl;
		cout << "       [" << x[iA] << " " << x[iB] << " " << x[iC] << " " << x[iD] << "]" << endl; 
		cout << "       [" << x[iA]/radius << " " << x[iB]/radius << " " << x[iC]/radius << " " << x[iD]/radius << "]" << endl; 
	}
}

void SparkProfilePartial(BzzVector &x, BzzVector &T, vector<double> &sparkRatio, const double Tpeak)
{
	int N = x.Size();
	double radius = x[1]/2.;

	int iA = 0;
	int iB = 0;
	int iC = 0;
	int iD = 0;

	for(int i=1;i<=N;i++)
		if ( x[i]/radius >= sparkRatio[0] )
		{
			iA = i;
			break;
		}
	for(int i=1;i<=N;i++)
		if ( x[i]/radius >= sparkRatio[1] )
		{
			iB = i;
			break;
		}
	for(int i=1;i<=N;i++)
		if ( x[i]/radius >= sparkRatio[2] )
		{
			iC = i;
			break;
		}
	for(int i=1;i<=N;i++)
		if ( x[i]/radius >= sparkRatio[3] )
		{
			iD = i;
			break;
		}


	// Ramp (I)
	for(int i=iA;i<iB;i++)
		T[i] = T[iA] + (Tpeak-T[iA])/(x[iB]-x[iA])*(x[i]-x[iA]);

	// Plateau
	for(int i=iB;i<iC;i++)
		T[i] = Tpeak;

	// Ramp (II)
	for(int i=iC;i<iD;i++)
		T[i] = Tpeak + (T[iD]-Tpeak)/(x[iD]-x[iC])*(x[i]-x[iC]);
}

void ReadFromBackupFile(OpenSMOKE_ReactingGas &mix, const string nameFile, BzzVector &x, BzzVector &T, BzzVector &m, BzzMatrix &omega)
{
		int N = x.Size();
		int NCGas = omega.Columns();
		int NCGasBackup;
		vector<string> namesBackup;

		BzzVector aux_x;
		BzzVector aux_m;
		BzzVector aux_T;
		BzzMatrix aux_OmegaGas;
		
		{
			cout << "Reading Backup file..." << endl;
			BzzLoad fBackup(nameFile);

			fBackup >> NCGasBackup;
			namesBackup.resize(NCGasBackup+1);
			cout << " * Number of species in backup file: " << NCGasBackup << endl;
			
			for (int j=1;j<=NCGas;j++)
			{
				string dummy;
				fBackup >> dummy;
				namesBackup[j] = dummy;
			}

			fBackup >> aux_x;
			fBackup >> aux_m;
			fBackup >> aux_T;
			fBackup >> aux_OmegaGas;
			fBackup.End();
		}

		cout << " * Interpolating from backup file" << endl;
		InterpolateBetweenGrids(aux_x, x, aux_m, m);
		InterpolateBetweenGrids(aux_x, x, aux_T, T);
			
		int countSpecies = 0;
		for(int j=1;j<=NCGasBackup;j++)
		{
			BzzVector aux(N);
			int k = mix.recognize_species_without_exit(namesBackup[j]);
			if (k>0)
			{
				countSpecies++;
				BzzVector auxvector = aux_OmegaGas.GetColumn(j);
				InterpolateBetweenGrids(aux_x, x, auxvector, aux);
				omega.SetColumn(k, aux);
				cout << " * Recognized species:   " << namesBackup[j] << " " << j << " => " << k << endl;
			}
			else
				cout << " * Unrecognized species: " << namesBackup[j] << endl;
		}
		cout << " * Number of recognized species: " << countSpecies << " over " << NCGasBackup << endl;

		cout << " * Normalization..." << endl;
		for(int j=1;j<=N;j++)
		{
			BzzVector aux = omega.GetRow(j);
			double sum = aux.GetSumElements();

			const double eps = 1.e-4;
			if ( (sum > 1.+eps) || (sum < 1.-eps))
			{
				cout << "Wrong backup file. The sum of mass fractions in point " << j << " is " << sum << endl;
				cout << "Press enter to exit..." << endl;
				getchar();
				exit(-1);
			}

			aux /= sum;
			omega.SetRow(j, aux);
		}
}

#endif  // OPENSMOKE_DROPLETMICROGAVITY_UTILITIES_H