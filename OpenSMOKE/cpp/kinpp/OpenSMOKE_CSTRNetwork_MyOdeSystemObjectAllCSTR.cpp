/***************************************************************************
 *   Copyright (C) 2003-2008 by                                            *
 *   Guido Buzzi-Ferraris, Alessio Frassoldati and Alberto Cuoci		   *						   *
 *                                                                         *
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

#include "kinpp/OpenSMOKE_CSTRNetwork.h"
#include "kinpp/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR.h"

OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR::OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR(OpenSMOKE_CSTRNetwork *cstrN)
{
	cstrNewtwork = cstrN;
	numComponents = cstrNewtwork->Reactions->NumberOfSpecies();
}

void OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR::ObjectBzzPrint(void)
{
	::BzzPrint("\nObject Print for Ode Jacobian");
}

void OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR::GetSystemFunctions(BzzVector &x,double t, BzzVector &f)
{
	double FMax, FMean;
	cstrNewtwork->countFourthTotal++;

	cstrNewtwork->GetAllResiduals(x,f);

	FMean = f.GetSumAbsElements() / double(x.Size());
	FMax = f.MaxAbs();

	// Write on video
	cout	<< "# = "		<< cstrNewtwork->countFourthTotal	<< "\t"
			<< "time = "	<< t								<< "\t"
			<< "Mean "		<< FMean							<< "\t"
			<< "Max "		<< FMax								<< "\t"
			<< endl;

	if (cstrNewtwork->OnlyODE == true)
		cstrNewtwork->fResiduals_4 << cstrNewtwork->countFourthTotal++ << "\t" << t << "\t" << FMean << "\t" << FMax << endl;
}

void OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR::GetJacobian(BzzVector &y, double t)
{
	// Viene scritta su file la matrice
	// corrispondente allo jacobiano del sistema non fattorizzato; in realta la matrice viene
	// approssimata come fosse perfettamente diagonale a blocchi, cioe' non sono stati
	// presi in considerazione gli elementi extradiagonale
	// y rappresenta il vettore delle frazioni massive di tuttOpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR::e le specie e di tutti
	// i reattori
	cstrNewtwork->GetDiagonalMatricesForLinearizedSistem(fileDiagonal,y);
}

void OpenSMOKE_CSTRNetwork_MyOdeSystemObjectAllCSTR::GetBuildJacobian(double hr, BzzMatrixSparseLockedByRows *Sh)
{
	int i,j,n;

	// Viene aperto il file contenente lo jacobiano, approssimato, del sistema
	// Sono stati completamente posti pari a zero i termini al di fuori della diagonale
	// (a blocchi) principale e viene letto il numero di blocchi
	BzzLoad load ('*', fileDiagonal);
	load >> n;

	// Viene aperto il file su cui verra' scritto lo jacobiano aperto sopra ma in
	// forma fattorizzata
	BzzSave save ('*', fileDiagonalFactorized);
	save << n;

	for(i = 1;i <= n;i++)
	{
		// Viene caricata la matrice jacobiana corrispondente al blocco i
		// Si tratta della matrice non fattorizzata
		load >> A;

		//
		Product(-hr,&A);
		for(j = 1;j <= A.Rows();j++)
			A[j][j] += 1.;

		// Viene scritta su file la stessa matrice, ma stavolta fattorizzata
		ReplaceBzzMatrixWithBzzFactorized(&A,&G);
		Factorize(&G);
		save << G;
	}
	save.End();
	load.End();

	Product(-hr,SL,Sh);
}

