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

#ifndef OPENSMOKE_SENSITIVITYANALYSIS
#define OPENSMOKE_SENSITIVITYANALYSIS

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_Kinetics;


class OpenSMOKE_SensitivityAnalysis
{

friend  class OpenSMOKE_Kinetics;

public:

	OpenSMOKE_SensitivityAnalysis();
	void Initialize(OpenSMOKE_ReactingGas *_mix, BzzVectorInt &indices, bool _iTemperature);
	void BuildJAlfaMatrix(BzzVector &r, const double T);
	void GiveMe_SensitivityCoefficients();
	void GiveMe_SensitivityCoefficients(const double deltat);
	void Normalize_SensitivityCoefficients(BzzVector &omega, BzzVector &x, const double MWmix, const double T);
	void PrintOnFile_SensitivityCoefficients(ofstream &fOutput, int count);
	void PrintOnFile_SensitivityCoefficients(ofstream &fOutput);

	BzzMatrix Jacobian;
	BzzMatrix Sx;
	BzzMatrix Somega;
	BzzMatrix S;
	BzzVector scaling;

private:

	OpenSMOKE_ReactingGas *mix;

	int NC;
	int NR;

	int NP;
	int NV;

	bool iTemperature;
	int	 indexTemperature;
	int	 indexSpecies;

	OpenSMOKE_NuManager *nu;

	BzzVectorInt indices_print_species;

	BzzMatrix JAlfa;

	BzzMatrix Sold;


	BzzMatrixDiagonal I;

	BzzVector A;
	BzzVector M;
	BzzVector uM;


	double	threshold_normalization;
	bool	iImplicit;

	void BuildNuMatrix(OpenSMOKE_Kinetics *kinetics);
};

#endif // OPENSMOKE_SENSITIVITYANALYSIS

