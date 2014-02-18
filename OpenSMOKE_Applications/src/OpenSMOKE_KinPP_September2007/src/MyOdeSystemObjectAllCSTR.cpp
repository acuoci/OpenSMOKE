#include "MyOdeSystemObjectAllCSTR.h"
#include "myBzzCSTRNetwork.hpp"


void MyOdeSystemObjectAllCSTR::ObjectBzzPrint(void)
{
	::BzzPrint("\nObject Print for Ode Jacobian");
}

void MyOdeSystemObjectAllCSTR::GetSystemFunctions(BzzVector &x,double t, BzzVector &f)
{
	double FMax, FMean;

	cstrNewtwork->GetAllResiduals(x,f);
	
	FMean = f.GetSumAbsElements() / double(x.Size());
	FMax = f.MaxAbs();
	cout	<< "time = "	<< t		<< "\t"
			<< "Mean "		<< FMean	<< "\t"
			<< "Max "		<< FMax		<< "\t" 
			<< endl;
}

void MyOdeSystemObjectAllCSTR:: GetJacobian(BzzVector &y, double t)
{
	// Viene scritta su file la matrice
	// corrispondente allo jacobiano del sistema non fattorizzato; in realta la matrice viene
	// approssimata come fosse perfettamente diagonale a blocchi, cioe' non sono stati
	// presi in considerazione gli elementi extradiagonale
	// y rappresenta il vettore delle frazioni massive di tutte le specie e di tutti
	// i reattori
	cstrNewtwork->GetDiagonalMatricesForLinearizedSistem(fileDiagonal,y);
}

void MyOdeSystemObjectAllCSTR::GetBuildJacobian(double hr, BzzMatrixSparseLockedByRows *Sh)
{
	int i,j,n;

	// Viene aperto il file contenente lo jacobiano, approssimato, del sistema
	// Sono stati completamente posti pari a zero i termini al di fuori della diagonale
	// (a blocchi) principale e viene letto il numero di blocchi
	BzzLoad load ('*',fileDiagonal);
	load >> n;

	// Viene aperto il file su cui verra' scritto lo jacobiano aperto sopra ma in 
	// forma fattorizzata
	BzzSave save ('*',fileDiagonalFactorized);
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

MyOdeSystemObjectAllCSTR::MyOdeSystemObjectAllCSTR(void)
{
	cstrNewtwork = 0;
}

void MyOdeSystemObjectAllCSTR::operator()(myBzzCSTRNetwork *cstrN)
{
	cstrNewtwork = cstrN;
	numComponents = cstrNewtwork->Reactions.NumberOfSpecies();
}