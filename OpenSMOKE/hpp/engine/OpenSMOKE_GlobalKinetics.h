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

#ifndef OpenSMOKE_GLOBALKINETICS
#define OpenSMOKE_GLOBALKINETICS

#include "BzzMath.hpp"

class OpenSMOKE_ReactingGas;


class OpenSMOKE_GlobalKinetics
{
public:

	int NC;
	int nReactions;

	void assign_mix(OpenSMOKE_ReactingGas *mixture);
	void read_from_file(std::string fileName);
	void GiveMeReactionHeat(const double T, BzzVector &R, double &QReaction);
	void GiveMeFormationRates(const double T, BzzVector &c, BzzVector &R);
	void GiveMeSpecificEnthalpy(const double T, BzzVector &omega, double &specificEnthalpy);

	void GiveMeFormationRatesPASR(const double T, BzzVector &c, BzzVector &R, double TauMix, BzzVector &TauC, BzzVector &TauR);
	void GiveMeFormationRatesED(const double T, const double rho, BzzVector &omega, BzzVector &c, const double TauMix, BzzVector &R);

    void GiveMeKEquilibrium(const double T);
    BzzVector DHReaction;
    BzzVector DSReaction;
    BzzVector sumNuij;
    BzzVector uKeq;
    BzzMatrix lambdaInverse;
    BzzVectorInt    iEquilibrium;
    BzzVector rInverse;

	// TODO
	BzzVector A;
	BzzVector Tatt;
	BzzVector Beta;
	BzzMatrix nu;
	BzzMatrix lambda;
	OpenSMOKE_ReactingGas *mix;

	BzzVector kappa;
	BzzVector r;
	BzzVector h;
	BzzVector s;

	// Sensitivity
	void GiveMedRdk(char kind, const int index, const double T, BzzVector &c, BzzVector &R);

  // Optimizer
	void ChangeKineticParameters(BzzVector &newParameters);
	void SetupOptimization(BzzVectorInt &_iKind, BzzVectorInt &_iReaction, BzzVectorInt &_iSpecies);

    void SensitivityCoefficients(const std::string fileName, BzzVector &grid, BzzVector &T, BzzMatrix &C);
    void SensitivityAnalysis(const std::string fileName, int indexT, int dimBlock,
                             BzzFactorizedTridiagonalBlocksGauss &J,
                             BzzVector &grid,
                             BzzVector &T, BzzMatrix &C,
                             BzzVector &rho, BzzVector &Cp);

	void KineticConstants(const double T);

	int IsTattAKineticParameter(const int index_reaction);
	void IsLambdaAKineticParameter(const int index_reaction, BzzVectorInt &index_parameters, BzzVectorInt &index_species);

	int rotation_index;
	double Cstar;
	double ALFA;
	double delta;

	void default_set();
	void change_set();

private:

	static const int MAXREACTIONS;

	void ErrorMessage(std::string message);

	void checkInputFile(const char *found, const char *expected);
	void checkInputFile(const int found, const int expected);

	void ReactionRates(const double T, BzzVector &c);
	void FormationRates(BzzVector &R);

	// Sensitivity Analysis
	void KineticConstants(const int index, const double T);

private:

	BzzVectorInt 	iKind;
	BzzVectorInt 	iSpecies;
	BzzVectorInt 	iReaction;

};

#endif // GLOBALKINETICS

/*

//	void SetupKineticParametersOptimization(std::string fileName);
public:
	int nParameters;
	BzzVector minValues;
	BzzVector maxValues;
	BzzVector startingPoint;

private:

	BzzVectorInt 	iKind;
	BzzVectorInt 	iSpecies;
	BzzVectorInt 	iReaction;
	void CheckKineticParametersOptimization();
	void 		InitializeStartingPoint();



*/
