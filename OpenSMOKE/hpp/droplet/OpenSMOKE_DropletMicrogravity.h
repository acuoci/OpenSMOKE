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

#if !defined(OPENSMOKE_DROPLETMICROGAVITY_H)
#define OPENSMOKE_DROPLETMICROGAVITY_H

#include "OpenSMOKE.hpp"
#include "droplet/OpenSMOKE_DropletMicrogravity_DAE_Objects.h"
#include "droplet/OpenSMOKE_DropletMicrogravity_DataManager.h"
#include "droplet/OpenSMOKE_DropletMicrogravity_GasRadiation.h"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_GlobalKinetics;

class OpenSMOKE_DropletMicrogravity
{

public:

	OpenSMOKE_DropletMicrogravity();
	void SetName(const std::string name);
	void Assign(OpenSMOKE_ReactingGas *_mix);
	void Assign(OpenSMOKE_DropletMicrogravity_DataManager *_data);
	void Assign(OpenSMOKE_GlobalKinetics *_global);

	void Setup();
	void Run();

	void DAESystemEigenValue(BzzVector &x, double t, BzzVector &f);
	void DAESystemUnsteadyBatch(BzzVector &x, double t, BzzVector &f);
	void DAE_myPrint(BzzVector &x, double t);

private:

	int N;
	int NE;
	int dimensionBlock;
	int NCGas;
	int NCLiquid;

	BzzVector rhoDroplet;
	BzzVector dHVaporization;
	BzzVector pvDroplet;
	BzzVector cpDroplet;

	BzzVector lambda;
	BzzVector rho;
	BzzVector Cp;
	BzzVector Kp;
	BzzVector QReaction;
	BzzMatrix RGas;

	BzzVector m;
	BzzVector dm_over_dt;
	BzzVector du_over_dr;

	BzzVector mass;
	BzzVector dmass_over_dt;

	BzzVector T;
	BzzVector dT_over_dt;
	BzzVector dT_over_dr;
	BzzVector dT_over_dr_central;

	BzzVector u;

	BzzVector T_convection;
	BzzVector T_conduction;
	BzzVector T_diffusion_fluxes;
	BzzVector T_reaction;

	BzzMatrix OmegaGas_diffusion;
	BzzMatrix OmegaGas_convection;
	BzzMatrix OmegaGas_reaction;

	BzzMatrix OmegaGas;
	BzzMatrix dOmegaGas_over_dt;
	BzzMatrix dOmegaGas_over_dr;

	BzzVector xMin;
	BzzVector xMax;
	BzzVector xFirstGuess;
	BzzVectorInt inDerAlg;

	BzzVector TInitialEstimation;
	BzzMatrix OmegaInitialEstimation;

	double tOld;

	void MemoryAllocation();
	void MemoryAllocationKinetics();
	void GiveMeDT(const double t);
	void GiveMeDOmegaGas(const double t);
	void GiveMeContinuityEquation(const double t);
	void GiveMeDropletMassEquation(const double t);
	
	void setMinimumAndMaximumValues();
	void setDifferentialAndAlgebraic();
	void setFirstGuess();
	void PrepareSystem();
	void recoverUnknowns(BzzVector &x);
	void recoverResiduals(BzzVector &f);

	void setInitialConditions();
	void UpdateProperties();
	void UpdateLiquidProperties();
	void GiveMeSpatialDerivatives();
	void GiveMeTimeDerivatives(const double t);
	void MoleFractionsAndMolecularWeights();
	void DiffusionVelocities();
	void AbsorptionCoefficients();
	void Velocities();
	void MassFlow();
	void InitialEstimations();

	void DAESystemSolution(BzzDaeSparseObject *o, double tEnd);
	int nonLinearSystemSolution(BzzNonLinearSystemSparseObject &o, int dimensionBlock);
	void PrintFinalSolution();
	void PrintFinalSolutionProperties();
	void PrintBackup();
	void PrintFinalSummary();
	void PrintTagInterfaceFile();
	void PrintTagFlameFile();
	void Summary();

	double SootPropensity(const int j);
	double SootPropensity();
	double SpeciesIntegralMass(const int j, const int indexSpecies);
	double MassIntegral(const std::string names);

	double InterfaceLiquidDensity();
	double InterfaceLiquidVaporizationHeat(const double t);
	double InterfaceLiquidMolecularWeight();
	double InterfaceLiquidSpecificHeat();
	
	double DropletMass();
	double DropletVolume();
	double DropletDiameter();

	OpenSMOKE_DropletMicrogravity_GasRadiation radiation;
	OpenSMOKE_Droplet* droplet;

private:

	BzzVector A_x_lambdae;
	BzzVector A_x_lambdaw;
	BzzVector A_x_rhoe;
	BzzVector A_x_rhow;

	BzzVector OmegaGasVector;
	BzzVector XGasVector;
	BzzVector CGasVector;
	BzzVector RGasVector;
	BzzVector DmixVector;
	BzzVector TetaMixVector;

	BzzMatrix XGas;
	BzzMatrix Dm;
	BzzMatrix Teta;

	BzzVector MW;
	BzzVector uMW;

	BzzMatrix Cpk;

	BzzMatrix vStar_e, vStar_w, coeff_e;
	BzzVector vc_e;

	ofstream fUnsteady;
	ofstream fUnsteadyRadiation;
	ofstream fInterface;
	ofstream fFlame;

private:

	OpenSMOKE_Grid1D grid;
	OpenSMOKE_ReactingGas *mix;
	OpenSMOKE_DropletMicrogravity_DataManager *data;
	OpenSMOKE_GlobalKinetics *global;

	OpenSMOKE_DropletMicrogravity_MyDaeSystemEigenValue daeSystemEigenValue;
	OpenSMOKE_DropletMicrogravity_MyDaeSystemUnsteadyBatch daeSystemUnsteadyBatch;

	static const double MMIN;

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};

#endif // OPENSMOKE_DROPLETMICROGAVITY_H

		
