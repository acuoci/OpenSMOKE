/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
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

#ifndef OPENSMOKE_EQUILIBRIUMSTANJAN
#define OPENSMOKE_EQUILIBRIUMSTANJAN

#include "BzzMath.hpp"

class OpenSMOKE_IdealGas;
class OpenSMOKE_EquilibriumStanjan;

enum exit_status
{
	CONVERGED, 
	FAILED, 
	GO_ON,
};


class OpenSMOKE_EquilibriumStanjan
{
public:

	// Default constructor
	OpenSMOKE_EquilibriumStanjan();

	void Reset();
	void SetName(string name);
	void SetVerbose();
	void UnsetVerbose();
	void Setup(OpenSMOKE_IdealGas *ideal_gas);
	void SetElementalComposition(BzzVector &_xInitialElemental);
	void SetElementalCompositionFromSpeciesMoleFractions(BzzVector &_xInitial);
	void Equilibrate(const double _T, const double _P_Pa, BzzVector &xFinal, double &Nfinal);

	void Mode1(double &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);
	void Mode2(double &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);
	void Mode3(double &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);
	void Mode2bis(double &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);


private:

	int NE;
	int NC;
	int NER;
	int NCR;
	
	string name_object;
	bool assignedElementalComposition;

	OpenSMOKE_IdealGas *gas;
	
	
	BzzVector xGas;
	BzzVector omegaGas;
	BzzVector cGas;
	BzzVector NGas;
	
	double		TInitial;
	double		PInitial;

	// Original variables
	BzzVector p;
	BzzMatrix E;
	BzzVector gtilde;


	// Reduced variables
	BzzVector pR;
	BzzVector molesR;
	BzzVector gtildeR;
	
	BzzMatrix			ER;
	BzzMatrix			ERT;
	BzzMatrix			ERT_x_ER;
	BzzFactorizedGauss	ERT_x_ER_Factorized;
	BzzMatrix			ER_x_ERT;
	BzzFactorizedGauss	ER_x_ERT_Factorized;
	BzzVector           *ERT_Row;

	BzzVector	BR;
	BzzVector	bR;

	BzzVector H;
	BzzVector f;
	BzzVector deltalambda;
	BzzMatrix Q;
	BzzFactorizedGauss Q_Factorized;
	
	BzzFactorizedGauss QExtended_Factorized;
	BzzVector D;
	BzzVector EE;
	double A;
	BzzVector bExtended;
	BzzVector deltaExtended;
	
	BzzVector lambdaFinalSolutionR;

	BzzVectorInt indexNoSpecies;
	BzzVectorInt indexIncludeSpecies;
	BzzVectorInt indexSpecies_Original;
	BzzVectorInt indexIncludeElement;
	BzzVectorInt indexElement_Original;

	void ReorderingSpecies();
	void GetGtilde(const double T, const double P_Pa);
	void GetFirstGuess(const double T, const double P_Pa);
	void GetLambdaFirstGuessFromLinearProgramming(double &NtotFirstGuess, BzzVector &xFirstGuessR, BzzVector &lambdaFirstGuessR);
	void GetLambdaFirstGuessFromLeastSquares(BzzVector &xFirstGuessR, BzzVector &lambdaFirstGuessR);

    void ErrorMessage(const string message);
    void WarningMessageStrong(const string message);
    void WarningMessageSoft(const string message);

	static const double threshold_element;

	bool	  iAlreadySolved;
	double    Database_NtotFirstGuess;
	BzzVector Database_lambdaFirstGuessR;
	BzzVector Database_xFirstGuessR;

private:

	int countMode1;
	int countMode2;
	int countMode2bis;
	int countMode3;

	bool iVerboseFile;
	ofstream fLog;

private:
	
	void GiveMeQ(double &Ntot, BzzVector &x);
	void GiveMeH(double &Ntot, BzzVector &x);
	void GiveMeD(BzzVector &x);
	void GiveMeA();
	void GiveMeZandW(BzzVector &x, BzzVector &lambda, double &Ntot, double &Z, double &W);
	void GiveMeZandWandV(BzzVector &x, BzzVector &lambda, double &Ntot, double &Z, double &W, double &V);

	void LinearProgramming(BzzVector &Nsolution);
	void LinearProgrammingMinMax(BzzVector &Nsolution);
	void Blending(double &NtotFirstGuess, BzzVector &xFirstGuessR, BzzVector &Nmin, BzzVector &Nmm);

	exit_status ModuleMode1(double &Ntot, BzzVector &lambda, BzzVector &x);
	exit_status ModuleMode2(double &Ntot, BzzVector &lambda, BzzVector &x);
	exit_status ModuleMode2bis(double &Ntot, BzzVector &lambda, BzzVector &x);
	exit_status Solution(double &Ntot, BzzVector &lambda, BzzVector &x);

	void FinalSummary(double &Ntot, BzzVector &xR, const double P_Pa, const double T);

private:
	
	int NmaxIterationsMode1;
	int NmaxIterationsMode2;
	int NmaxIterations;

	double toleranceMode1;
	double toleranceMode2;
	double toleranceSteeping;

	int NCR_Gas;
};

#endif // OPENSMOKE_EQUILIBRIUM


