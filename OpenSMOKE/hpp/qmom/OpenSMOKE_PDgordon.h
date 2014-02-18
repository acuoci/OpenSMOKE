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

#ifndef		OpenSMOKE_QMOM_PDGORDON
#define		OpenSMOKE_QMOM_PDGORDON

#if MKL==0
	// Nothing to do
#else
	#include "mkl.h"
#endif
			
#include "BzzMath.hpp"
#include "qmom/OpenSMOKE_Distributions.h"

class OpenSMOKE_PDgordon;

class MyNonLinearSystem : public BzzMyNonLinearSystemObject
{
public:
	void assignPDgordon(OpenSMOKE_PDgordon *PD);

	OpenSMOKE_PDgordon *ptPD;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};


class OpenSMOKE_PDgordon
{
private:

	// Constants
	int N;
	int _2N_plus_1;
	int _2N_minus_1;

	// Work array and thei dimensions
	// (eigenvalues and eigenvectors computation)
	double *work;
	int LWORK;
	int *iwork;
	int LIWORK;

	// Variables to calculate eigenvalues/eigenvectors
	int Nd;
	int info;
	double *diagonal;
	double *subdiagonal;
	double *eigenvectors;

	// Arrays and matrices
	BzzVector mr;
	BzzMatrix P;
	BzzMatrix J;

	BzzVector alfa;
	BzzVector a;
	BzzVector b;

	void allocateMemory();

	// Flags
	int iSetup;
	int iCalculate;
	double m0;

	OpenSMOKE_MultiDiracDistribution mdirac;

	MyNonLinearSystem NLS;

public:

	BzzVector csi;
	BzzVector w;
	BzzVector chi;

	void setup(BzzVector &moments);
	void calculateAbscissasAndWeigths();
	void printVideo();

	void initialize();

	void nonLinearSystem(BzzVector &x, BzzVector &f);
};

#endif // PDGORDON
