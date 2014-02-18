/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci   	   *
 *   alberto.cuoci@polimi.it   						   *
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

#ifndef OPENSMOKE_PREPROCESSORKINETICS_H
#define OPENSMOKE_PREPROCESSORKINETICS_H

#include "BzzMath.hpp"

class OpenSMOKE_PreProcessorReactingGas;

class OpenSMOKE_PreProcessorKinetics
{
public:

	// Variabili per contare il numero di reazioni di un certo tipo
	int numEquilibrium;
	int	numMortoCheParla;
	int	numThirdBodyOnly;
	int	numFallOff;
	int	numCABR;
	int	numLandauTeller;
	int	numJanevLanger;
	int	numPowerSeries;
	int	numChebishev;
	int	numThirdBody;
	int	numConventional;
	int	numIrreversible;
	int	numLogarithmicPressure;
	int	numCollisionEfficiency;
	int	numTAR;

	// Vettori contenenti le informazioni stechiometriche organizate a seconda
	// dei coeficienti stechiometrici con cui le spcie compaiono all'interno
	// delle reazioni
	BzzVectorInt	numDir1,	numDir2,	numDir3,	numDir4,	numDir5,
					numInvTot1, numInvTot2,	numInvTot3,	numInvTot4, numInvTot5,
					numInvEq1,	numInvEq2,	numInvEq3,	numInvEq4,	numInvEq5,
					jDir1,		jDir2,		jDir3,		jDir4,		jDir5,
					jInvTot1,	jInvTot2,	jInvTot3,	jInvTot4,	jInvTot5,
					jInvEq1,	jInvEq2,	jInvEq3,	jInvEq4,	jInvEq5;

	BzzVector valDir5,
		            valInvTot5,
					valInvEq5;

	BzzVector	sumNuij;
	BzzVector	forwardOrders;
	BzzVector	backwardOrders;


	// Vettori per la gestione delle reazioni con Equilibrio
	BzzVectorInt	jEquil;				// 0 se non c'e' equilibrio
	BzzVectorInt	reactionWithEquil;	// indice j della reazione // dimensionato per jEquil=1

	// Vettori per la gestione delle reazioni con Terzo Corpo
	BzzVector	*efficiency;
	BzzVectorInt jThirdBody;
	BzzVectorInt jNumEfficiency;
	BzzVectorInt *iEfficiency;

	// Informazioni cinetiche per le reazioni dirette e quelle di FallOff
	BzzVector	beta1, E1;
	BzzVector	beta2, E2;
	BzzVector	k01,k02;
	BzzVector	aPressure,bPressure,cPressure,dPressure,ePressure;
	BzzVector	exp_k01;

	// Vettori contenenti gli indici delle reazioni "non convenzionali"
	BzzVectorInt iThirdBodyOnly;
	BzzVectorInt iThirdBody;
	BzzVectorInt iFallOff;
	BzzVectorInt iCABR;
	BzzVectorInt iLandauTeller;
	BzzVectorInt iJanevLanger;
	BzzVectorInt iPowerSeries;
	BzzVectorInt iChebishev;
	BzzVectorInt iLogarithmicPressure;
	BzzVectorInt iCollisionEfficiency;
	BzzVectorInt iTAR;
	BzzVectorInt negativeSigns;

	void readStoichiometricFile(const string fileSt);
	void readStoichiometricFile_2(const string fileSt);

public:

	OpenSMOKE_PreProcessorReactingGas *reactionRates;

	// Private variables
	int		NC;
	int		NR;
	string	name_object;

	// Lettura delle informazioni dello schema cinetico da file
	void readFromFile(const string fileSt, const string fileKin);
	void SparsityStructures(const string fileSt);
	void CheckingStoichiometry();

	void ErrorMessage(const string message);
	void WarningMessage(const string message);

private:

	static const int NUM_MAX;
};

#endif


