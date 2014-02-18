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

#ifndef OPENSMOKE_SENSITIVITYANALYSIS_DIFUSIONE_FLAME1D
#define OPENSMOKE_SENSITIVITYANALYSIS_DIFUSIONE_FLAME1D

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_Kinetics;

class OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D
{

public:

	OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D();
	void SetName(const string name);
	void Initialize(const int kind_of_flame, OpenSMOKE_ReactingGas *_mix, BzzVectorInt &_indices, const int _N);

	void PrintOnFile(BzzVector &x, BzzVector &T, BzzVector &H, BzzMatrix &omega, 
					 BzzVector &U, BzzVector &G, BzzVector &rho, BzzVector &MWtot,
					 vector<string> list_of_names, BzzMatrix &S, BzzMatrix &Diffusivity);

	void VideoSummary();

	
private:

	string name_object;
	string reacting_system;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);

private:

	OpenSMOKE_ReactingGas	*mix;

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
};

#endif // OPENSMOKE_SENSITIVITYANALYSIS_DIFFUSION_FLAME1D

