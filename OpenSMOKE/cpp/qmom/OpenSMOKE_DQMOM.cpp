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

#include "qmom/OpenSMOKE_DQMOM.h"
#include "qmom/OpenSMOKE_PhysicalModels.h"

void OpenSMOKE_DQMOM::setup(int N)
{
	// Dimensions
	_N = N;
	_2N = 2*N;
	_N_plus_1 = N+1;
	_2N_minus_1 = _2N - 1;

	// Memory Allocation
	ChangeDimensions(_2N, _2N,	&A);
	ChangeDimensions(_2N, _N,	&powerOfCsi);
	ChangeDimensions(_2N, &coeffA1);
	ChangeDimensions(_2N, &coeffA2);
	ChangeDimensions(_2N, &Source);
	ChangeDimensions(_2N, &alfa);
	ChangeDimensions(_N, &a);
	ChangeDimensions(_N, &b);

	// Matrix A preparation
	prepareCoefficients();
	prepareA();
}

void OpenSMOKE_DQMOM::update(OpenSMOKE_PhysicalModels &_models)
{
	// -----------------------------------------------------
	// Updating of matrix A and Vector S : ATTENTION
	// -----------------------------------------------------
	ChangeDimensions(_2N, _2N, &A);
	Source = 0.;

	powerOfCsi	= _models.powerOfCsi;
	
	prepareA();	
	assemblingA();
}


void OpenSMOKE_DQMOM::prepareCoefficients()
{
	int i;

	// -----------------------------------------------------
	// Coefficients
	// -----------------------------------------------------
	for(i=2;i<=_2N_minus_1;i++)
	{
		coeffA1[i+1] = (1-i);
		coeffA2[i+1] = i;
	}
}

void OpenSMOKE_DQMOM::prepareA()
{
	int i;

	// -----------------------------------------------------
	// First two rows matrix A
	// -----------------------------------------------------
	for(i=1;i<=_N;i++)				// a. 1st Row - A1 portion
		A[1][i] = 1.;
	for(i=_N_plus_1;i<=_2N;i++)		// b. 1st Row - A2 portion
		A[1][i] = 0.;
	for(i=1;i<=_N;i++)				// c. 2nd Row - A1 portion
		A[2][i] = 0.;
	for(i=_N_plus_1;i<=_2N;i++)		// d. 2nd Row - A2 portion
		A[2][i] = 1.;
}

void OpenSMOKE_DQMOM::assemblingA()
{
	int i, j;
	int row;

	// The first two rows of A are constant and therefore it is not
	// necessary to calculate them

	for(i=3;i<=_2N;i++)
	{
		// Matrix A1
		row = i;		
		for(j=1;j<=_N;j++)
			A[i][j] = coeffA1[i]*powerOfCsi[row][j];

		// Matrix A2
		row = i-1;
		for(j=1;j<=_N;j++)
			A[i][j+_N] = coeffA2[i]*powerOfCsi[row][j];
	}
}

void OpenSMOKE_DQMOM::solveLinearSystem()
{
	int i;
	
	if (_N==20)
	{
		double L1 = A[3][2] - A[3][1];
		double L2 = -A[3][4] + A[3][3];
		double L3 = A[4][2] - A[4][1];
		double L4 = -A[4][4] + A[4][3];
		double D1 = Source[3] - A[3][3]*Source[2] - A[3][1]*Source[1];
		double D2 = Source[4] - A[4][3]*Source[2] - A[4][1]*Source[1];
	
		alfa[4] = (D2*L1 - D1*L3)/(-L1*L4+L2*L3);
		alfa[2] = (D1+L2*alfa[4])/L1;
		alfa[3] = Source[2]-alfa[4];
		alfa[1] = Source[1]-alfa[2];		
	}
	else
		Solve(A, Source, &alfa);


	for(i=1;i<=_N;i++)
		a[i] = alfa[i];

	for(i=1;i<=_N;i++)
		b[i] = alfa[i+_N];

}





