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

class MyNonLinearSystem_EquilibriumStanjan : public BzzMyNonLinearSystemObject
{
public:
	void AssignEquilibrium(OpenSMOKE_EquilibriumStanjan *equilibrium);

	OpenSMOKE_EquilibriumStanjan *equilibrium;
	virtual void GetResiduals(BzzVector &lambda, BzzVector &f);
	virtual void ObjectBzzPrint(void);
};

enum exit_status
{
	CONVERGED, 
	FAILED, 
	GO_ON,
	NO_SOLID_PHASE
};


class OpenSMOKE_EquilibriumStanjan
{
public:

	// Default constructor
	OpenSMOKE_EquilibriumStanjan();

	void SetName(string name);
	void Setup(OpenSMOKE_IdealGas *ideal_gas);
	void SetElementalCompositionFromSpeciesMoleFractions(BzzVector &_xInitial);
	void Equilibrate(const double _T, const double _P_Pa, BzzVector &_xFinal);

	void Mode1(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);
	void Mode2(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);
	void Mode2bis(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);
	void Mode3(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);
	void Mode2bisSinglePhase(BzzVector &Ntot, BzzVector &lambda, BzzVector &x, exit_status &flag);


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

	BzzVector xSolid;
	BzzVector omegaSolid;
	BzzVector NSolid;

	BzzVector xLiquid;
	BzzVector omegaLiquid;
	BzzVector NLiquid;

	
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
	BzzMatrix D;
	BzzMatrix EE;
	BzzMatrix A;
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

private:

	int countMode1;
	int countMode2;
	int countMode2bis;
	int countMode3;

	bool iVerboseFile;
	ofstream fLog;

private:
	
	void GiveMeQ(BzzVector &Ntot, BzzVector &x);
	void GiveMeH(BzzVector &Ntot, BzzVector &x);
	void GiveMeD(BzzVector &x);
	void GiveMeA();
	void GiveMeZandW(BzzVector &x, BzzVector &lambda, BzzVector &Ntot, BzzVector &Z, double &W);
	void GiveMeZandWandV(BzzVector &x, BzzVector &lambda, BzzVector &Ntot, BzzVector &Z, double &W, BzzVector &V);

	void LinearProgramming(BzzVector &Nsolution);
	void LinearProgrammingMinMax(BzzVector &Nsolution);
	void Blending(BzzVector &NtotFirstGuess, BzzVector &xFirstGuessR, BzzVector &Nmin, BzzVector &Nmm);

	exit_status ModuleMode1(BzzVector &Ntot, BzzVector &lambda, BzzVector &x);
	exit_status ModuleMode2(BzzVector &Ntot, BzzVector &lambda, BzzVector &x);
	exit_status ModuleMode2bis(BzzVector &Ntot, BzzVector &lambda, BzzVector &x);
	exit_status Solution(BzzVector &Ntot, BzzVector &lambda, BzzVector &x);

	void FinalSummary(BzzVector &Ntot, BzzVector &xR, const double P_Pa, const double T);

private:
	
	int NmaxIterationsMode1;
	int NmaxIterationsMode2;
	int NmaxIterations;

	double toleranceMode1;
	double toleranceMode2;
	double toleranceSteeping;

	bool iLiquidPhase;
	bool iSolidPhase;
	bool iMultiPhase;
	int nPhases;

	int NCR_Gas;
	int NCR_Solid;
	int NCR_Liquid;

	BzzVectorInt listOfSolid;
	BzzVectorInt listOfLiquid;

	bool iSolidAbsent;
	bool iLiquidAbsent;
};

#endif // OPENSMOKE_EQUILIBRIUM


