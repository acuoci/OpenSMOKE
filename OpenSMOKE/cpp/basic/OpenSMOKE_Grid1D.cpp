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
#include "basic/OpenSMOKE_Grid1D.h"

void OpenSMOKE_Grid1D::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Grid1D"		<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Grid1D::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Grid1D"		<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
} 

OpenSMOKE_Grid1D::OpenSMOKE_Grid1D()
{
	name_object = "Undefined name";
	iSpherical = false;
}

void OpenSMOKE_Grid1D::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Grid1D::SetSpherical()
{
	iSpherical = true;
}

// 1. Construct equispaced grid
void OpenSMOKE_Grid1D::Construct(const int NP, const double LL, const double xAA)
{
	double deltax;

	Np = NP;
	Ni = Np-1;

	L = LL;
	xA = xAA;
	xB = xA+L;

	Allocate();

	deltax = L/Ni;

	dxw = deltax;
	dxw[1] = 0.;

	Build();
}

// 2. Construct stretched grid (on left side)
void OpenSMOKE_Grid1D::Construct(const int NP, const double LL, const double alfa, const double xAA)
{
	Np = NP;
	Ni = Np-1;

	L = LL;
	xA = xAA;
	xB = xA+L;

	Allocate();

	if (iSpherical == false)
	{
		int i;
		double sum, dx1;
		double fraction_of_points_A		= 0.80;
		double fraction_of_distance_A	= 0.70;
		double alfaB = 2.;
		int NpA = int(fraction_of_points_A*double(Np));
		int NiA = NpA-1;
		int NpB = Np-NiA;
		int NiB = NpB-1;

		sum=1.;
		for(i=2;i<=NiA;i++)
			sum+=pow(alfa, (i-1));

		dx1=L*fraction_of_distance_A/sum;
		dxw[1]=0.;
		for(i=2;i<=NpA;i++)
			dxw[i]=dx1*pow(alfa, i-2);

		sum=1.;
		for(i=2;i<=NiB;i++)
			sum+=pow(alfaB, (i-1));

		dx1=L*(1.-fraction_of_distance_A)/sum;
		dxw[Np]=dx1;
		for(i=Ni;i>=NpA+1;i--)
		{
			int j=Np-i;
			dxw[i]=dx1*pow(alfaB, j);
		}
	}
	else
	{
		double sum=1.;
		for(int i=1;i<=Np-2;i++)
			sum+=pow(alfa, double(i));
		double dx0=L/sum;

		dxw[1]=0.;
		for(int i=2;i<=Np;i++)
			dxw[i]=dx0*pow(alfa, double(i-2));
	}

	Build();
}

// 2. Construct stretched grid (on left side)
void OpenSMOKE_Grid1D::ConstructStretchedStagnation(const int NP, const double LL, const double alfa, const double xAA)
{
	Np = NP;
	Ni = Np-1;

	L = LL;
	xA = xAA;
	xB = xA+L;

	Allocate();

	{
		double sum=1.;
		for(int i=1;i<=Np-2;i++)
			sum+=pow(alfa, double(i));
		double dx0=L/sum;

		dxw[1]=0.;
		for(int i=2;i<=Np;i++)
			dxw[Np+2-i]=dx0*pow(alfa, double(i-2));
	}

	Build();
}


// 2. Construct stretched grid (on both sides)
void OpenSMOKE_Grid1D::Construct(const int NP, const double xStart, const double Length, 
								 const double alfaA, const double alfaB,
								 const double fraction_of_points_A, double fraction_of_distance_A)
{
	Np = NP;
	Ni = Np-1;

	L = Length;
	xA = xStart;
	xB = xA+L;

	Allocate();

	{
		int i;
		double sum, dx1;
		int NpA = int(fraction_of_points_A*double(Np));
		int NiA = NpA-1;
		int NpB = Np-NiA;
		int NiB = NpB-1;

		sum=1.;
		for(i=2;i<=NiA;i++)
			sum+=pow(alfaA, (i-1));

		dx1=L*fraction_of_distance_A/sum;
		dxw[1]=0.;
		for(i=2;i<=NpA;i++)
			dxw[i]=dx1*pow(alfaA, i-2);

		sum=1.;
		for(i=2;i<=NiB;i++)
			sum+=pow(alfaB, (i-1));

		dx1=L*(1.-fraction_of_distance_A)/sum;
		dxw[Np]=dx1;
		for(i=Ni;i>=NpA+1;i--)
		{
			int j=Np-i;
			dxw[i]=dx1*pow(alfaB, j);
		}
	}
	Build();
}


// 2. Construct centered grid
void OpenSMOKE_Grid1D::Construct(const int NP, const double LL, const double xcen, const double wmix, const double xAA)
{
	int i;
	double deltax;

	Np = NP;
	Ni = Np-1;

	L = LL;
	xA = xAA;
	xB = xA+L;

	Allocate();

	// Controlli
	// -----------------------------------------------------
	if (wmix>0.80*L)
	{
		cout << "Error: wmix>0.80L" << endl; exit(1);
	}
	if ((xcen-wmix*0.50)<=0. || (xcen+wmix*0.50)>=L)
	{
		cout << "Error: wmix too large" << endl; exit(1);
	}
	// ------------------------------------------------------

	int n = 3;
	// Lato combustibile
	// ----------------------------
	deltax = (xcen-wmix*0.50)/double(n);
	dxw[1] = 0.;
	for(i=2;i<=n+1;i++)
		dxw[i] = deltax;

	// Lato ossidante
	// ----------------------------
	deltax = (L-(xcen+wmix*0.50))/double(n);
	for(i=Np-n+1;i<=Np;i++)
		dxw[i] = deltax;

	// Mixing zone
	// ----------------------------
	deltax = (wmix)/(Ni-n*2);
	for(i=n+2;i<=Np-n;i++)
		dxw[i] = deltax;

	Build();
}

// 4. Construct grid from file
void OpenSMOKE_Grid1D::Construct(const std::string fileName)
{
	BzzVector xVector;
	std::string dummy, units;

	ifstream inputFile;
	openInputFileAndControl(inputFile, fileName);

	
	inputFile >> dummy;
	if (dummy != "#X")
		ErrorMessage("Expected #X keyword");
	inputFile >> units;
	inputFile >> dummy;
	if (dummy != "#LIST")
		ErrorMessage("Expected #LIST keyword");

	for(;;)
	{
		inputFile >> dummy;
		if (dummy == "#END")
			break;
		xVector.Append(OpenSMOKE_Conversions::conversion_length(atof(dummy.c_str()), units));
	} 

	inputFile.close();

	Np = xVector.Size();
	Ni = Np-1;

	L = xVector[Np];
	xA = xVector[1];
	xB = xA+L;

	Allocate();

	dxw[1]=0.;
	for(int i=2;i<=Np;i++)
		dxw[i]=xVector[i]-xVector[i-1];

	Build();
}

// 5. Construct from a vector of points
void OpenSMOKE_Grid1D::Construct(BzzVector &coordinates)
{
	int i;

	Np = coordinates.Size();
	Ni = Np-1;


	L = coordinates[Np];
	xA = coordinates[1];
	xB = xA+L;

	Allocate();

	dxw[1]=0.;
	for(i=2;i<=Np;i++)
		dxw[i]=coordinates[i]-coordinates[i-1];

	Build();
}

// 1. Construct equispaced grid
void OpenSMOKE_Grid1D::ConstructExponential(const int NP, const double LL, const double Beta)
{
	Np = NP;
	Ni = Np-1;

	L = LL;
	xA = 0.;
	xB = xA+L;

	Allocate();

	BzzVector z(Np);
	for (int j=1;j<=Np;j++)
		z[Np-j+1] = 1.+Beta*(1.-(1.-exp(-1./Beta))*double(Np-j)/double(Np-1));
	z[1] = 0.;
	z *= LL;

	dxw[1] = 0.;
	for(int j=2;j<=Np;j++)
		dxw[j]=z[j]-z[j-1];

	Build();
}

void OpenSMOKE_Grid1D::Build()
{
	int i;

	// Costruzione di tutte le altre informazioni
	x[1]=xA;
	for(i=2;i<=Np;i++)
		x[i]=x[i-1]+dxw[i];

	for(i=1;i<=Np;i++)
		x2[i] = x[i]*x[i];

	if(x[1]!=0.) ux[1] = 1./x[1];
	else ux[1] = 0.;
	for(i=2;i<=Np;i++)
		ux[i] = 1./x[i];

	for(i=2;i<=Np;i++)
		udxw[i] = 1./dxw[i];

	for(i=1;i<=Np-1;i++)
	{
		dxe[i] = x[i+1]-x[i];
		udxe[i] = 1./dxe[i];
	}

	for(i=2;i<=Ni;i++)
	{
		dxc[i] = (x[i+1]-x[i-1]);
		udxc[i] = 1./dxc[i];
		dxc_over_2[i] = 0.50*dxc[i];
		udxc_over_2[i] = 1./dxc_over_2[i];
	}

	for(i=2;i<=Ni;i++)
	{
		cenp[i] = dxw[i]*udxe[i]*udxc[i];
		cenc[i] = (dxe[i]-dxw[i])*udxe[i]*udxw[i];
		cenm[i] = dxe[i]*udxw[i]*udxc[i];
	}


	if (iSpherical == true)
	{
		for(int j=1;j<Np;j++)
			spherical[j] = ( 3./4.*( BzzPow4(x[j+1])-BzzPow4(x[j]) ) /( BzzPow3(x[j+1])-BzzPow3(x[j]) ) - x[j] ) 
															/ (x[j+1]-x[j]);
		for(int j=1;j<Np;j++)
			one_minus_spherical[j] = 1.-spherical[j];

		for(int j=1;j<=Np;j++)
			A[j] = 4.*Constants::pi*x2[j];
	}
}

void OpenSMOKE_Grid1D::Allocate()
{
	ChangeDimensions(Np, &x);
	ChangeDimensions(Np, &dxw);
	ChangeDimensions(Np, &dxe);
	ChangeDimensions(Np, &ux);
	ChangeDimensions(Np, &x2);
	ChangeDimensions(Np, &udxw);
	ChangeDimensions(Np, &udxe);
	ChangeDimensions(Np, &dxc);
	ChangeDimensions(Np, &udxc);
	ChangeDimensions(Np, &dxc_over_2);
	ChangeDimensions(Np, &udxc_over_2);
	ChangeDimensions(Np, &cenp);
	ChangeDimensions(Np, &cenc);
	ChangeDimensions(Np, &cenm);

	if (iSpherical == true)
	{
		ChangeDimensions(Np,	&A);
		ChangeDimensions(Np-1,	&spherical);
		ChangeDimensions(Np-1,	&one_minus_spherical);
	}
}

void OpenSMOKE_Grid1D::RefineDouble()
{
	int i;

	BzzVector dxwOld = dxw;
	int NpOld = Np;

	Np = NpOld+Ni;
	Ni = Np-1;

	Allocate();

	dxw[1] = 0.;
	for(i=2;i<=NpOld;i++)
	{
		dxw[(i-1)*2] = dxwOld[i]*0.50;
		dxw[(i-1)*2+1] = dxwOld[i]*0.50;
	}

	Build();
}

void OpenSMOKE_Grid1D::Refine(const int j)
{
	int i;

	BzzVector dxwOld = dxw;
	int NpOld = Np;

	Np = NpOld+1;
	Ni = Np-1;

	Allocate();

	dxw[1] = 0.;
	for(i=2;i<=j;i++)
		dxw[i]	= dxwOld[i];

	dxw[j+1]	= dxwOld[j+1] * 0.50;
	dxw[j+2]	= dxwOld[j+1] * 0.50;
	
	for(i=j+3;i<=Np;i++)
		dxw[i] = dxwOld[i-1];

	Build();
}

void OpenSMOKE_Grid1D::RecoverFromBackUp(const std::string fileName)
{
	int dummy;

	BzzLoad fInput('*', fileName);

	fInput >> Np;		
	Ni = Np-1;
	fInput >> dummy;
	Allocate();
	fInput >> x;
	fInput.End();

	L  = x[Np];
	xA = x[1];
	xB = xA+L;

	dxw[1]=0.;
	for(int i=2;i<=Np;i++)
		dxw[i] = x[i] - x[i-1];

	Build();
}

void OpenSMOKE_Grid1D::RecoverFromBackUp(BzzVector &_x)
{

	Np = _x.Size();;		
	Ni = Np-1;
	Allocate();
	x  = _x;

	L = x[Np];
	xA = x[1];
	xB = xA+L;

	dxw[1]=0.;
	for(int i=2;i<=Np;i++)
		dxw[i] = x[i] - x[i-1];

	Build();
}

// First Approach
BzzVectorInt OpenSMOKE_Grid1D::QueryNewPointsDifference(const int nDiff, const double delta, BzzVector &phi)
{
	int i;

	BzzVector diffPhi(Ni);
	BzzVectorInt pointList;
	BzzVectorInt iList;
	BzzVector vList;

	int index;
	double maximumValue;
	double maxValue = phi.Max();
	double minValue = phi.Min();
	double coefficient = delta*(maxValue-minValue);

	for(i=1;i<=Ni;i++)
		diffPhi[i] = fabs(phi[i+1]-phi[i]);

	for(i=1;i<=Ni;i++)
		if (diffPhi[i] > coefficient)
		{
			iList.Append(i);
			vList.Append(diffPhi[i]);
		}

	for(i=1;i<=min(iList.Size(), nDiff);i++)
	{
		maximumValue = vList.Max(&index);
		pointList.Append(iList[index]);
		vList[index] = -1.;
	}
	
	return pointList;
}

// Second Approach
BzzVectorInt OpenSMOKE_Grid1D::QueryNewPointsGradient(const int nGrad, const double delta, BzzVector &phi)
{
	int i;

	BzzVector diffPhi(Ni);
	BzzVector dphi(Np);
	BzzVectorInt pointList;
	BzzVectorInt iList;
	BzzVector vList;

	int index;
	double maximumValue;

	FirstDerivative('C', phi, dphi);
	
	double maxValue = dphi.Max();
	double minValue = dphi.Min();
	double coefficient = delta*(maxValue-minValue);


	for(i=1;i<=Ni;i++)
		diffPhi[i] = fabs(dphi[i+1]-dphi[i]);


	for(i=1;i<=Ni;i++)
		if (diffPhi[i] > coefficient)
		{
			iList.Append(i);
			vList.Append(diffPhi[i]);
		}

	for(i=1;i<=min(iList.Size(), nGrad);i++)
	{
		maximumValue = vList.Max(&index);
		pointList.Append(iList[index]);
		vList[index] = -1.;
	}

	return pointList;
}


void OpenSMOKE_Grid1D::FirstDerivative(const char KIND, BzzVector &u, BzzVector &v, BzzVector &dv)
{
	int i;

	if (KIND == 'U')
	{
		dv[1] = (v[2]-v[1])*udxe[1];
		for(i=2;i<=Ni;i++)
			dv[i] = (u[i]<=0.) ?  (-v[i]+v[i+1])*udxe[i] : (-v[i-1]+v[i])*udxw[i];
		dv[Np] = (v[Np]-v[Ni])*udxw[Np];
	}

	else
		FirstDerivative(KIND, v, dv);
}

void OpenSMOKE_Grid1D::FirstDerivative(const char KIND, BzzVector &u, BzzMatrix &v, BzzMatrix &dv)
{
	int i, j;
	int nCol = v.Columns();

	if (KIND == 'U')
	{
		for (j=1;j<=nCol;j++)
			dv[1][j] = (v[2][j]-v[1][j])*udxe[1];
		
		for(i=2;i<=Ni;i++)
			for (j=1;j<=nCol;j++)
			dv[i][j] = (u[i]<=0.) ?  (-v[i][j]+v[i+1][j])*udxe[i] : (-v[i-1][j]+v[i][j])*udxw[i];

		for (j=1;j<=nCol;j++)
			dv[Np][j] = (v[Np][j]-v[Ni][j])*udxw[Np];
	}
	else
		FirstDerivative(KIND, v, dv);
}


void OpenSMOKE_Grid1D::FirstDerivative(const char KIND, BzzVector &v, BzzVector &dv)
{
	int i;

	if (KIND == 'C')
	{
		dv[1] = (v[2]-v[1])*udxe[1];
		for(i=2;i<=Ni;i++)
			dv[i] = - cenm[i]*v[i-1] + cenc[i]*v[i] + cenp[i]*v[i+1];
		dv[Np] = (v[Np]-v[Ni])*udxw[Np];
	}
	else if (KIND == 'B')
	{
		dv[1] = (v[2]-v[1])*udxe[1];
		for(i=2;i<=Np;i++)
			dv[i] =  (-v[i-1]+v[i])*udxw[i];
	}
	else if (KIND == 'F')
	{
		for(i=1;i<=Ni;i++)
			dv[i] =  (-v[i]+v[i+1])*udxe[i];
		dv[Np] = (v[Np]-v[Ni])*udxw[Np];
	}
	else
		ErrorMessage("Wrong FirstDerivative option: " + KIND);
}

void OpenSMOKE_Grid1D::SecondDerivative(BzzVector &v, BzzVector &d2v)
{
	int i;
	for(i=2;i<=Ni;i++)
		d2v[i] = ((v[i+1]-v[i])*udxe[i]-(v[i]-v[i-1])*udxw[i]) *udxc_over_2[i];
}

void OpenSMOKE_Grid1D::SecondDerivative(BzzMatrix &v, BzzMatrix &d2v)
{
	int i, j;
	int nCol = v.Columns();

	for(i=2;i<=Ni;i++)
		for(j=1;j<=nCol;j++)
			d2v[i][j] = ((v[i+1][j]-v[i][j])*udxe[i]-(v[i][j]-v[i-1][j])*udxw[i])*udxc_over_2[i];
}

void OpenSMOKE_Grid1D::FirstDerivative(const char KIND, BzzMatrix &v, BzzMatrix &dv)
{
	int i, j;
	int nCol = v.Columns();

	if (KIND == 'C')
	{
		for (j=1;j<=nCol;j++)
			dv[1][j] = (v[2][j]-v[1][j])*udxe[1];

		for(i=2;i<=Ni;i++)
			for (j=1;j<=nCol;j++)
				dv[i][j] = - cenm[i]*v[i-1][j] + cenc[i]*v[i][j] + cenp[i]*v[i+1][j];

		for (j=1;j<=nCol;j++)
			dv[Np][j] = (v[Np][j]-v[Ni][j])*udxw[Np];
	}
	else if (KIND == 'B')
	{
		for (j=1;j<=nCol;j++)
			dv[1][j] = (v[2][j]-v[1][j])*udxe[1];
		for(i=2;i<=Np;i++)
			for (j=1;j<=nCol;j++)
				dv[i][j] =  (-v[i-1][j]+v[i][j])*udxw[i];
	}
	else if (KIND == 'F')
	{
		for(i=1;i<=Ni;i++)
			for (j=1;j<=nCol;j++)
				dv[i][j] =  (-v[i][j]+v[i+1][j])*udxe[i];
		dv[Np][j] = (v[Np][j]-v[Ni][j])*udxw[Np];
	}
}

void OpenSMOKE_Grid1D::DoubleField(BzzVector &phi)
{
	int NpOld;

	BzzVector phiOld = phi;
	NpOld = phiOld.Size();
	ChangeDimensions(Np, &phi);

	// Linear Interpolation
	phi[1] = phiOld[1];
	for(int i=2;i<=NpOld;i++)
	{
		int index		= (i-1)*2;
		phi[index]		= 0.50*(phiOld[i-1]+phiOld[i]);
		phi[index+1]	= phiOld[i];
	}
}

void OpenSMOKE_Grid1D::AdaptField(BzzVector &phi, BzzVector &xNew)
{
	bool iFind;
	BzzVector phiOld = phi;

	int j, k;
	int jOld = 1;
	
	for(j=2;j<=Np-1;j++)
	{
		iFind = false;
		for(k=jOld;k<=Np-1;k++)
		{
			if (iFind == true) break;

			if (xNew[j] > x[k] && xNew[j] <= x[k+1])
			{
				double m = (xNew[j]-x[k])/(x[k+1]-x[k]);
				phi[j]  = phiOld[k]+m*(phiOld[k+1]-phiOld[k]);
				iFind = true;
				jOld = k;
			}
		}

		if (iFind == false)
		{
			cout << "Error in point " << j << " xNew= " << xNew[j] << endl;
			cout << "xNew from" << xNew.Min() << " to " << xNew.Max() << endl;
			cout << "xOld from" << x.Min() << " to " << x.Max() << endl;
			for(int kk=1;kk<=Np;kk++)
				cout << kk << " " << xNew[kk] << " " << x[kk] << endl;
			ErrorMessage("Point outside the boundaries!");
		}
	}
}

void OpenSMOKE_Grid1D::AdaptField(BzzMatrix &phi, BzzVector &xNew)
{
	bool iFind;
	BzzMatrix phiOld = phi;

	int j, k;
	int jOld = 1;
	
	for(j=2;j<=Np-1;j++)
	{
		iFind = false;
		for(k=jOld;k<=Np-1;k++)
		{
			if (iFind == true) break;

			if (xNew[j] > x[k] && xNew[j] <= x[k+1])
			{
				double m = (xNew[j]-x[k])/(x[k+1]-x[k]);
				for(int i=1;i<=phi.Columns();i++)
					phi[j][i]  = phiOld[k][i]+m*(phiOld[k+1][i]-phiOld[k][i]);
				iFind = true;
				jOld = k;
			}
		}
		if (iFind == false)
			ErrorMessage("Point outside the boundaries!");
	}
}

void OpenSMOKE_Grid1D::DoubleField(BzzMatrix &phi)
{
	int NpOld;
	int NCols;

	BzzMatrix phiOld = phi;
	NpOld = phiOld.Rows();
	NCols = phiOld.Columns();
	ChangeDimensions(Np, NCols, &phi);

	// Linear Interpolation
	for(int j=1;j<=NCols;j++)	phi[1][j] = phiOld[1][j];
	for(int i=2;i<=NpOld;i++)
	{
		int index = (i-1)*2;
		for(int j=1;j<=NCols;j++)
		{
			phi[index][j]	= 0.50*(phiOld[i-1][j]+phiOld[i][j]);
			phi[index+1][j] = phiOld[i][j];
		}
	}
}

void OpenSMOKE_Grid1D::AddPointsField(BzzVector &phi, BzzVectorInt &listPoints)
{
	int i, k;
	int NpOld;

	BzzVector phiOld = phi;
	NpOld = phiOld.Size();
	ChangeDimensions(Np, &phi);

	int count = 0;
	for (i=1;i<=NpOld;i++)
	{
		phi[i+count] = phiOld[i];
		if (count < listPoints.Size())
			if(listPoints[count+1] == i)	count++;
	}

	BzzVectorInt indices(listPoints.Size());
	for (k=1;k<=listPoints.Size();k++)
		indices[k]= listPoints[k]+k;
	
	for (k=1;k<=listPoints.Size();k++)
		phi[indices[k]] = 0.50*(phi[indices[k]-1]+phi[indices[k]+1]);
}

void OpenSMOKE_Grid1D::AddPointsField(BzzMatrix &phi, BzzVectorInt &listPoints)
{
	int i, j, k;
	int NpOld;
	int NCols;

	BzzMatrix phiOld = phi;
	NpOld = phiOld.Rows();
	NCols = phiOld.Columns();
	ChangeDimensions(Np, NCols, &phi);

	int count = 0;
	for (i=1;i<=NpOld;i++)
	{
		for(j=1;j<=NCols;j++)	phi[i+count][j] = phiOld[i][j];
		if (count < listPoints.Size())
			if(listPoints[count+1] == i)	count++;
	}

	BzzVectorInt indices(listPoints.Size());
	for (k=1;k<=listPoints.Size();k++)
		indices[k]= listPoints[k]+k;
	
	for (k=1;k<=listPoints.Size();k++)
		for(j=1;j<=NCols;j++)	phi[indices[k]][j] = 0.50*(phi[indices[k]-1][j]+phi[indices[k]+1][j]);
}