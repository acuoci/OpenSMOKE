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

#ifndef OPENSMOKE_SENSITIVITYANALYSIS_FLAME1D
#define OPENSMOKE_SENSITIVITYANALYSIS_FLAME1D

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_Kinetics;
class OpenSMOKE_NuManager;

class OpenSMOKE_SensitivityAnalysis_Flame1D
{

friend  class OpenSMOKE_Kinetics;

public:

	OpenSMOKE_SensitivityAnalysis_Flame1D();
	void SetName(const string name);
	void Initialize(const int kind_of_flame, OpenSMOKE_ReactingGas *_mix, BzzVectorInt &_indices, const int _N, const kindOfSensitivityParameter _kindOfParameter);
	void BuildJAlfaMatrix(const int iPoint, const double T, const double rho, const double Cp);
	void BuildJAlfaMatrixBoundary(const int iPoint, const double T, const double rho, const double Cp);
	void GiveMe_SensitivityCoefficients();
	
	void LocalNormalization(const int iParameter, BzzMatrix &omega, BzzVector &MWmix, BzzVector &T, BzzVector &u, BzzVector &rho, BzzVector &H);
	void GlobalNormalization(const int iParameter, BzzMatrix &omega, BzzMatrix &x, BzzVector &T, BzzVector &u, BzzVector &H);
	
	void PrintOnFileSpeciesGlobalSensitivityCoefficients(const string name_of_species, BzzVector &x);
	void PrintOnFileTemperatureGlobalSensitivityCoefficients(BzzVector &x);
	void PrintOnFileFlameSpeedGlobalSensitivityCoefficients(BzzVector &x);
	void PrintOnFileMassFlowRateLocalSensitivityCoefficients(BzzVector &x);


	void PrintOnFileSpeciesLocalSensitivityCoefficients(const string name_of_species, BzzVector &x);
	void PrintOnFileTemperatureLocalSensitivityCoefficients(BzzVector &x);
	void PrintOnFileFlameSpeedLocalSensitivityCoefficients(BzzVector &x);
	void PrintOnFile(BzzVector &x, BzzVector &T, BzzVector &H, BzzMatrix &omega, 
					 BzzVector &U, BzzVector &G, BzzVector &rho, BzzVector &MWtot,
					 vector<string> list_of_names);

	void VideoSummary();

	BzzMatrix S;
	
	BzzMatrix SLx;
	BzzMatrix SLomega;
	BzzMatrix SLT;
	BzzMatrix SLu;
	BzzMatrix SLH;

	BzzMatrix SGx;
	BzzMatrix SGomega;
	BzzMatrix SGT;
	BzzMatrix SGu;
	BzzMatrix SGH;

	void GiveMe_SensitivityCoefficients(BzzFactorizedTridiagonalBlocksGauss &_Jacobian, BzzMatrix &_Jalfa, BzzVector &_parameter, BzzVector &x);

	BzzFactorizedTridiagonalBlocksGauss *AFactorized;
	
private:

	string name_object;
	string reacting_system;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	void BuildNuMatrix(OpenSMOKE_Kinetics *kinetics);

	kindOfSensitivityParameter kindOfParameter;

private:

	OpenSMOKE_ReactingGas	*mix;
	OpenSMOKE_NuManager		*nu;

	int NC;
	int NR;
	int NP;
	int NV;
	int NE;
	int N;

	int	 indexTemperature;
	int	 indexSpecies;
	int	 indexMassFlowRate;
	int	 indexVelocity;							// index of velocity
	int	 indexPressureCurvature;		// index of pressure curvature		

	BzzVectorInt indices_print_species;

	BzzMatrix JAlfa;
	BzzMatrix JAlfaPoint;

	BzzVector M;
	BzzVector uM;
	BzzVector parameters;

	double	threshold_normalization;
};

#endif // OPENSMOKE_SENSITIVITYANALYSIS_FLAME1D

