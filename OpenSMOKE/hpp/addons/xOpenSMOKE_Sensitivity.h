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

#ifndef	OPENSMOKE_SENSITIVITY
#define OPENSMOKE_SENSITIVITY

#include "BzzMath.hpp"

class OpenSMOKE_GlobalKinetics;
class OpenSMOKE_ReactingGas;


class OpenSMOKE_Sensitivity
{
public:

	OpenSMOKE_Sensitivity();
	void SetName(const string name);

	void GiveMe_Jalfa_Matrix(OpenSMOKE_GlobalKinetics &global, const double T, const double rho, BzzVector &c);
	void GiveMe_Jalfa_Matrix(OpenSMOKE_ReactingGas &mix, const double T, const double rho, const double cTot, BzzVector &c);

	void GiveMe_SensitivityCoefficients(BzzMatrix &J, BzzVector &omega);
	void GiveMe_SensitivityCoefficients(const double deltat, BzzMatrix &J, BzzVector &omega);

	void setup(char _kind, int NSpecies, string *listOfNames, OpenSMOKE_ReactingGas &mix);
	void setup(OpenSMOKE_ReactingGas &mix, OpenSMOKE_GlobalKinetics &global);

	void writeForBars(string fileName, OpenSMOKE_ReactingGas &mix);
	void writeForBars(const double coordinate);


    // ------------------------------------------------------------------------------------------------------------------------------
    // One dimensional problems
    // ------------------------------------------------------------------------------------------------------------------------------
    void GiveMe_Jalfa_Matrix(OpenSMOKE_ReactingGas &mix, BzzVector &T, BzzMatrix &W, BzzMatrix &X,
                             BzzVector &rho, BzzVector &PMtot);

	void GiveMe_SensitivityCoefficients(BzzFactorizedTridiagonalBlocksGauss &J, BzzMatrix &omega, BzzMatrix &X, BzzVector &T);
	void writeForBars(string fileName, OpenSMOKE_ReactingGas &mix, BzzVector &coordinate);

    void preProcessing(	string _kindOfProblem,  int numberOfPoints, int numberOfSpecies, int _iNormalization);
    void fromMassFlowRateToFlameSpeed(OpenSMOKE_ReactingGas &mix, BzzVector &u, BzzVector &T, BzzVector &rho, BzzVector &PMtot);


    void close();

private:

	string name_object;

    int NC;
    int NR;

	int Np;
	int dimBlock;
	int nAdd;
	int indexT;
	int indexH;

	int iNormalization;
	string kindOfProblem;
    char kind;

	int iAllFlag;
	BzzVectorInt listOfIndexSpecies;

    // Zero dimensional problems
    BzzVector JVector;
    BzzVector kappa;

    // One dimensional problems
    BzzMatrix JalfaTridiagonal;
	BzzMatrix kappaTridiagonal;

    BzzMatrix Eold;
    BzzMatrix E;
    BzzMatrix nE;

    BzzMatrix Jalfa;
	BzzMatrixDiagonal I;

	ofstream fSensitivity;

	void Normalize(BzzVector &omega);
	void NormalizeMassFractions_Local(BzzMatrix &omega, BzzVector &T);
	void NormalizeMassFractions_Maximum(BzzMatrix &omega, BzzVector &T);

    void allocate(const int numberOfSpecies, const int numberOfReactions);

	void setupOutputFile(string fileName, OpenSMOKE_ReactingGas &mix);
	void setupOutputFile(string fileName, OpenSMOKE_ReactingGas &mix, OpenSMOKE_GlobalKinetics &global);

	void printJalfaTridiagonalMatrixOnFile(char *fileName);

	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_SENSITIVITY