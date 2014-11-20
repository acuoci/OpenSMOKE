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

#ifndef OPENSMOKE_KINETICS_H
#define OPENSMOKE_KINETICS_H

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"
#include <vector>

class OpenSMOKE_ReactingGas;
class OpenSMOKE_NuManager;
class OpenSMOKE_RateOfProductionAnalysis;
class OpenSMOKE_SensitivityAnalysis;
class OpenSMOKE_ChebishevPolynomialsReaction;
class OpenSMOKE_LogarithmicPressureReaction;

class OpenSMOKE_Kinetics
{
friend class  OpenSMOKE_FluxAnalysis;
friend  class OpenSMOKE_Nu;
friend  class OpenSMOKE_SensitivityAnalysis;
friend  class OpenSMOKE_SensitivityAnalysis_Flame1D;
friend  class OpenSMOKE_SensitivityAnalysis_Fast_Flame1D;
friend  class OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D;
friend  class OpenSMOKE_ReactingGas;
friend  class OpenSMOKE_CSTRNetwork;
friend  class myBzzCSTRNetwork;
friend void CorrectionCoefficient_DeltaDirac( OpenSMOKE_ReactingGas *mix, const double T, const double cTot, const double Sigma, const double Tmin, const double Tmax, const double MaxCoeffCorr, BzzVector &correction_k1, BzzVector &correction_uKeq, BzzVector &correction_k2);

public:

	OpenSMOKE_Kinetics();
	void SetName(const std::string name);
	OpenSMOKE_ReactingGas *reactionRates;

	void StoichiometricMatrix(BzzMatrix &gamma);

	int numEquilibrium;
	int numThirdBodyOnly;
	int numFallOff;
	int numThirdBody;
	int numConventional;
	int numIrreversible;
	int numCABR;
	int numLandauTeller;
	int numJanevLanger;
	int numPowerSeries;
	int numChebishev;
	int numLogarithmicPressure;
	int numGlobal;
	int numCollisionEfficiency;
	int numTAR;

	BzzVectorInt	reactionWithEquil;	// indice j della reazione // dimensionato per jEquil=1
	BzzVectorInt	iFallOff;
	BzzVector		k01;
	BzzVector		exp_k01;
	BzzVector		exp_k02;

	double GiveMe_A(const int j);
	double GiveMe_E(const int j);
	double GiveMe_Beta(const int j);
	double GiveMe_ForwardOrder(const int j);

	void UpdateKineticParameters(const int iModel, solid_regression_parameters &parameters, BzzVector &b, BzzMatrixInt &user);
	void UpdateOptimizerParameters(const int iModel, solid_regression_parameters &parameters, BzzVector &b, BzzMatrixInt &user);

	void ChangeFrequencyFactor(const int j, const double variation);
	void ChangeActivationEnergy(const int j, const double new_value_cal_mol);


private:

	int NR;
	int	NC;

	BzzVectorInt	numDir1,	numDir2,	numDir3,	numDir4,	numDir5;
	BzzVectorInt	numInvTot1, numInvTot2,	numInvTot3,	numInvTot4, numInvTot5;
	BzzVectorInt	numInvEq1,	numInvEq2,	numInvEq3,	numInvEq4,	numInvEq5;
	BzzVectorInt	jDir1,		jDir2,		jDir3,		jDir4,		jDir5;
	BzzVectorInt	jInvTot1,	jInvTot2,	jInvTot3,	jInvTot4,	jInvTot5;
	BzzVectorInt	jInvEq1,	jInvEq2,	jInvEq3,	jInvEq4,	jInvEq5;

	BzzVector valDir5;
	BzzVector valInvTot5;
	BzzVector valInvEq5;

	BzzVectorInt	lambda_jDir1, lambda_jDir2, lambda_jDir3, lambda_jDir4, lambda_jDir5;
	BzzVectorInt	lambda_jInvEq1, lambda_jInvEq2, lambda_jInvEq3, lambda_jInvEq4, lambda_jInvEq5;
	BzzVectorInt	lambda_numDir1, lambda_numDir2, lambda_numDir3, lambda_numDir4, lambda_numDir5;
	BzzVectorInt	lambda_numInvEq1, lambda_numInvEq2, lambda_numInvEq3, lambda_numInvEq4, lambda_numInvEq5;
	BzzVector		lambda_valDir5;
	BzzVector		lambda_valInvEq5;

	BzzVector		*efficiency;

	BzzVectorInt	jNumEfficiency;
	BzzVectorInt	*iEfficiency;
	
	BzzVectorInt	jEquil;				// 0 se non c'e' equilibrio

	BzzVector		forwardOrders;
	BzzVector		backwardOrders;

public:
	BzzVectorInt	jThirdBody;
   BzzVector		sumNuij;
	BzzVector	beta1;
	BzzVector	E1;
	BzzVector	beta2;
	BzzVector	E2;
	BzzVector	k02;
	
private:

	BzzVector	aPressure;
	BzzVector	bPressure;
	BzzVector	cPressure;
	BzzVector	dPressure;
	BzzVector	ePressure;

	BzzMatrix	LandauTeller;
	BzzMatrix	JanevLanger;
	BzzMatrix	PowerSeries;
	BzzMatrix	TAR_series;
	OpenSMOKE_ChebishevPolynomialsReaction* ChebishevPolynomials;
	OpenSMOKE_LogarithmicPressureReaction*	LogarithmicPressure;	
	BzzVector	CollisionEfficiency;

	BzzVectorInt	iThirdBodyOnly;
	BzzVectorInt	iThirdBody;
	BzzVectorInt	iCABR;
	BzzVectorInt	iLandauTeller;
	BzzVectorInt	iJanevLanger;
	BzzVectorInt	iPowerSeries;
	BzzVectorInt	iChebishev;
	BzzVectorInt	iLogarithmicPressure;
	BzzVectorInt	iCollisionEfficiency;
	BzzVectorInt	iTAR;
	BzzVectorInt	negativeSigns;

private:

	BzzVectorInt	 iGlobal;
	BzzVectorInt	*iGlobalDirect;
	BzzVectorInt	*iGlobalInverse;
	BzzVector		*lambdaGlobalDirect;
	BzzVector		*lambdaGlobalInverse;

public:
	BzzVector mc;
	BzzVector mR;
	BzzVector mr;
	BzzVector mrDirC;
	BzzVector mrDirT;
	BzzVector mrInvC;

public:

	// Lettura delle informazioni dello schema cinetico da file
	void readFromFileBinary(const std::string fileKin);
	void readStoichiometricFileBinary(BzzLoad &fInput);
	BzzMatrix constructGamma(const std::string fileSt);

	// Calcolo dei contributi alle velocita di reazione
	void ComputeDirectAndInverse(BzzVector& c , BzzVector &rDirC, BzzVector &rInvC);
	void GlobalReactionRates(BzzVector& c, BzzVector &rDirC, BzzVector &rInvC);
	void ThirdBody(BzzVector &c, double cTot, BzzVector &coeffM);
	void FallOff(double T, BzzVector &coeffM, BzzVector &k1, BzzVector &k2, BzzVector &logFcent, BzzVector &coeffFallOff);
	void ChemicallyActivatedBimolecularReactions(double T, BzzVector &coeffM, BzzVector &k1, BzzVector &k2, BzzVector &logFcent, BzzVector &coeffFallOff);

	// Assemblaggio dei contributi per la costruzione delle velocita di reazione
	void reactionsWithEquilibrium(BzzVector &rDirC, BzzVector &rInvC,BzzVector &uKeq,BzzVector &r);
	void reactionsWithThirdBody(BzzVector &coeffM,BzzVector &r);

	void reactionsWithEquilibriumForwardAndBackward(BzzVector &rBackward, BzzVector &uKeq);
	void reactionsWithThirdBodyForwardAndBackward(BzzVector &coeffM, BzzVector &rForward, BzzVector &rBackward);

	// Calcolo delle velocita di reazione per le specie 1...NC
	void compositionReactionRates(BzzVector &r, BzzVector *R);

	// Calcolo dei parametri cinetici delle reazioni
	void ComputeKineticParameters(BzzVector &reactionDH, BzzVector &reactionDS, double T, double P_Pa, double logT, double uT, double loguRT, BzzVector &uKeq, BzzVector &k1, BzzVector &k2, BzzVector &logFcent);

	// Calcolo dei deltaH e deltaS delle reazioni dello schema cinetico
	void ComputeDHandDSreaction(BzzVector &componentDH,BzzVector &componentDS,BzzVector &reactionDH,BzzVector &reactionDS);



	void GetDerivativesC(double T, double cTot, BzzMatrix *dRC, BzzVector &cRes, BzzVector &R);
	
	void GiveMeIndexOfSpeciesInEachReaction(const std::string fileName, BzzVectorIntArray &indices);

	int IsAReversibleReaction(const int index);
	int IsAFallOffReaction(const int index);
	int IsACABReaction(const int index);
	int IsALandauTellerReaction(const int index);
	int IsAJanevLangerReaction(const int index);
	int IsAPowerSeriesReaction(const int index);
	int IsAChebishevPolynomialsReaction(const int index);
	int IsALogarithmicPressureReaction(const int index);
	int IsACollisionEfficiencyReaction(const int index);
	int IsATARReaction(const int index);

	void KineticExpressionString(ofstream &fOutput, const int j);

	BzzVectorInt ReactionIndices(const vector<string> list_of_species);

	BzzVector GiveMeSumNuDirect();

	void UpdateOmegaC(const double omegaC);

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

private:
	static const int NUM_MAX;
	static const int NUM_MAX_CHEBISHEV;
	static const int NUM_MAX_LOGPRESSURE;
	static const double S1;
	static const double S2;
	static const double S3;
	static const double SSQ;

	double TAR_omegaC;
};

#endif


