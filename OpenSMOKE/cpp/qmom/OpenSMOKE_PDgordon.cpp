/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci   	                               *
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
#include "qmom/OpenSMOKE_PDgordon.h"


void OpenSMOKE_PDgordon::initialize()
{
	iSetup = 0;
	NLS.assignPDgordon(this);
}

void OpenSMOKE_PDgordon::setup(BzzVector &moments) 
{
	m0 = moments[1];
	mr = moments/m0;;

	if (iSetup!=moments.Size())
	{
		// ----------------------------------------------------------
		// Dimensions of array and matrices
		// ----------------------------------------------------------
		N = moments.Size() / 2;
		Nd = N;
		_2N_plus_1	= 2*N+1;
		_2N_minus_1 = 2*N-1;

		//int kappa;
		//if(N==2)				kappa = 1;
		//else if(N==3 || N==4)	kappa = 2;
		//else if(N>=5 && N<=8)	kappa = 3;
		//LWORK	= 2*N*N + (3+2*kappa)*N+1;
		//LIWORK	= 5*N+3;

		LWORK  = 1+4*N+N*N;
		LIWORK = 3+5*N;

		// ----------------------------------------------------------
		// Memory allocation
		// ----------------------------------------------------------
		allocateMemory();

		iSetup = moments.Size();
	}

	iCalculate = 0;
}

void OpenSMOKE_PDgordon::allocateMemory() 
{
	ChangeDimensions(_2N_plus_1, _2N_plus_1, &P);
	ChangeDimensions(_2N_minus_1, _2N_minus_1, &J);
	ChangeDimensions(2*N, &alfa);
	ChangeDimensions(N, &a);
	ChangeDimensions(N-1, &b);
	ChangeDimensions(N, &csi);
	ChangeDimensions(N, &w);
	ChangeDimensions(N, &chi);

	diagonal		= (double*) malloc(N * sizeof(double));
	subdiagonal		= (double*) malloc((N) * sizeof(double));
	eigenvectors	= (double*) malloc((N*N) * sizeof(double));
	work			= (double*) malloc((LWORK) * sizeof(double));
	iwork			= (int*) malloc((LIWORK) * sizeof(int));
}



void OpenSMOKE_PDgordon::calculateAbscissasAndWeigths() 
{
	int i, j;
	
	if (iCalculate == 0)
	{	
		if (N == 2)
		{
			alfa[1] = 0.;
			alfa[2] = mr[2];
			alfa[3] = (mr[3]-mr[2]*mr[2])/mr[2];
			alfa[4] = (mr[2]*mr[4]-mr[3]*mr[3]) /
				      (mr[2]*mr[3]-mr[2]*mr[2]*mr[2]);
		}
		
		else
		{
			// ----------------------------------------------------------
			// Construction of matrix P
			// ----------------------------------------------------------
			P[1][1] = 1.;
			P[1][2] = 1.;
			double coefficient=1;
			for(i=2;i<=_2N_plus_1;i++)
			{
				coefficient*=-1;
				P[i][2] = coefficient*mr[i];
			}

			for(j=3;j<=_2N_plus_1;j++)
				for(i=1;i<=2*N+2-j;i++)
					P[i][j] = P[1][j-1]*P[i+1][j-2]-P[1][j-2]*P[i+1][j-1];

			// ----------------------------------------------------------
			// Construction of vector alfa
			// ----------------------------------------------------------
			alfa[1] = 0.;
			for(i=2;i<=2*N;i++)
				alfa[i] = P[1][i+1]/(P[1][i]*P[1][i-1]);
		}

		// ----------------------------------------------------------
		// Construction of Jacobi matrix
		// ----------------------------------------------------------
		for(i=1;i<=N;i++)
			a[i] = alfa[2*i-1]+alfa[2*i];
	
		for(i=1;i<=N-1;i++)
			b[i] = sqrt(fabs(alfa[2*i]*alfa[2*i+1]));	// See McGraw


		// ----------------------------------------------------------
		// Weigths and abscissas
		// ----------------------------------------------------------
		if (N == 2)
		{
			double coefficient = a[2]*a[2] - 2.*a[1]*a[2] + a[1]*a[1] + 4.*b[1]*b[1];
	
			csi[1]	= 0.50*a[2]+0.50*a[1]+0.50*sqrt(coefficient); 
			csi[2]	= 0.50*a[2]+0.50*a[1]-0.50*sqrt(coefficient);

			w[1]	= 0.25*m0*BzzPow2(-a[2]+a[1]+sqrt(coefficient))/b[1]/b[1];
			w[2]	= 0.25*m0*BzzPow2(-a[2]+a[1]-sqrt(coefficient))/b[1]/b[1];
		}

		else
		{
			FromBzzVectorToC(N, a, diagonal);
			FromBzzVectorToC(N-1, b, subdiagonal);

			#if MKL==0
				cout << "The MKL Library is not available!" << endl;
				cout << "Press enter to continue...       " << endl;
				getchar();
				exit(-1);
			#else
				char kind = 'I';
				dstedc(&kind, &Nd, diagonal, subdiagonal, eigenvectors, &Nd, work, &LWORK, iwork, &LIWORK, &info);
			#endif

			if (info != 0)
			{
				cout << "Error in computing eigenvalues/eigenvectors!" << endl;
				cout << "Error kind: info = " << info << endl;
				cout << "Press enter to continue... " << endl;
				getchar();
				exit(-1);
			}

			// ----------------------------------------------------------
			// Calculation of abscissas and weights
			// ----------------------------------------------------------
			FromCToBzzVector(N, csi, diagonal);
			for (i=1;i<=N;i++)
				w[i] = m0*eigenvectors[N*(i-1)]*eigenvectors[N*(i-1)];
		}

		// ----------------------------------------------------------
		// Calculation of weigthed abscissas
		// ----------------------------------------------------------
		ElementByElementProduct(w,csi, &chi);
		iCalculate = 1;
	}
}

void OpenSMOKE_PDgordon::printVideo() 
{
	int i;

	cout.setf(ios::scientific);

	cout << endl;
	cout << "n\tm\t\ta\t\tb\t\tcsi\t\tw\t" << endl;
	cout << "-------------------------------------------------------------------------------------" << endl;
		cout << "0" << "\t" << mr[1] << endl;
	for(i=1;i<=N-1;i++)
		cout <<  i  << "\t" << mr[i+1] << "\t" << a[i] << "\t" << b[i]
		            << "\t" << csi[i]  << "\t" << w[i] << endl;
		cout <<  N  << "\t" << mr[N+1] << "\t" << a[N] << "\t" << "\t"
		            << "\t" << csi[N]  << "\t" << w[N] << endl;

	for(i=N+1;i<=2*N-1;i++)
		cout <<  i  << "\t" << mr[i+1] << endl;

	cout << endl;
}

// NON LINEAR SYSTEM
// ---------------------------------------------------------------------------------
void MyNonLinearSystem::ObjectBzzPrint(void)
{
}
void MyNonLinearSystem::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptPD->nonLinearSystem(x, f);
}
void MyNonLinearSystem::assignPDgordon(OpenSMOKE_PDgordon *PD)
{
	ptPD = PD;
}

void OpenSMOKE_PDgordon::nonLinearSystem(BzzVector &x, BzzVector &f)
{
	int i, j, k;
	
	k=1;
	for(i=1;i<=N;i++)
		w[i] = x[k++];
	for(i=1;i<=N;i++)
		csi[i] = x[k++];

	f = -mr;
	for(i=1;i<=2*N;i++)
		for(j=1;j<=N;j++)
			f[i] += w[j]*pow(csi[j], i-1);
}