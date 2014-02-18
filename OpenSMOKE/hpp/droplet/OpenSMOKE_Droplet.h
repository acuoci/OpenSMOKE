/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
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

#if !defined(OPENSMOKE_DROPLET_H)
#define OPENSMOKE_DROPLET_H

#include "OpenSMOKE.hpp"
#include "droplet/OpenSMOKE_DropletMicrogravity_DAE_Objects.h"
#include "droplet/OpenSMOKE_DropletMicrogravity_DataManager.h"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_GlobalKinetics;

class OpenSMOKE_Droplet
{

public:

	OpenSMOKE_Droplet();
	void SetName(const string name);
	void Assign(OpenSMOKE_ReactingGas *_mix);
	void Assign(OpenSMOKE_DropletMicrogravity_DataManager *_data);

	void Setup();
	void SetInitialConditions();
	void Run();
	double TSurface()	{return T[N];}
	double conductionFluxInterface() { return lambda[N]*dT_over_dr[N]*grid.A[N]; }

	void UpdateFromGasPhase(const double tStart_, const double tEnd_);
//	void UpdateFromGasPhase(const double mGas_, const double conductionGas_, const double radiationFlux_);
	void UpdateFromGasPhase(const double TGasInterface_);
	void DAESystemNoMomentum(BzzVector &x, double t, BzzVector &f);
	void DAESystemOnlyTemperature(BzzVector &x, double t, BzzVector &f);
	void DAE_myPrint(BzzVector &x, double t);
	BzzVector T;
private:

	int N;
	int NE;
	int dimensionBlock;
	int NC;

//	double mGas;
//	double conductionGas;
//	double radiationFlux;
	double TGasInterface;

	BzzVector rho;
	BzzVector cp;
	BzzVector lambda;


	BzzVector dT_over_dt;
	BzzVector dT_over_dr;
	BzzVector T_conduction;

	BzzMatrix Omega_diffusion;

	BzzMatrix Omega;
	BzzMatrix dOmega_over_dt;
	BzzMatrix dOmega_over_dr;

	BzzVector OmegaVector;
	BzzVector XVector;

	BzzVector xMin;
	BzzVector xMax;
	BzzVector xFirstGuess;
	BzzVectorInt inDerAlg;

	void MemoryAllocation();
	void GiveMeDT(const double t);
	void GiveMeDOmega(const double t);
	
	void setMinimumAndMaximumValues();
	void setDifferentialAndAlgebraic();
	void setFirstGuess();
	void PrepareSystem();
	void recoverUnknowns(BzzVector &x);
	void recoverResiduals(BzzVector &f);

	void UpdateProperties();
	void GiveMeSpatialDerivatives();
	void MoleFractionsAndMolecularWeights();
	void DiffusionVelocities();

	void DAESystemSolution(BzzDaeSparseObject *o);

	double GetDensity(const double t, const BzzVector &omega);
	double GetSpecificHeat(const double t, const BzzVector &omega);
	double GetThermalConductivity(const double t, const BzzVector &omega);
	BzzVector GetDiffusionCoefficients(const double t, const BzzVector &omega);

private:

	BzzVector A_x_lambdae;
	BzzVector A_x_lambdaw;
	BzzVector A_x_rhoe;
	BzzVector A_x_rhow;

	BzzMatrix X;
	BzzMatrix Dm;

	BzzVector MW;
	BzzVector uMW;

	BzzMatrix vStar_e, vStar_w, coeff_e;
	BzzVector vc_e;

	double tStart;
	double tEnd;

	ofstream fUnsteady;

private:

	OpenSMOKE_Grid1D grid;
	OpenSMOKE_ReactingGas *mix;
	OpenSMOKE_DropletMicrogravity_DataManager *data;

	OpenSMOKE_Droplet_MyDaeSystemNoMomentum daeSystemNoMomentum;
	OpenSMOKE_Droplet_MyDaeSystemNoMomentum daeSystemOnlyTemperature;

	enum dropletGridMode {LIQUID_DROPLET_GRID_FIXED_RATIO, LIQUID_DROPLET_GRID_EXPONENTIAL} gridMode;
	double stretchingFactor;

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_DROPLET_H

		
