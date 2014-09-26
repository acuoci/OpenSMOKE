/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
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

#include <iomanip>
#include "sstream"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_FlameSpeedManager.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_OpposedFlameManager.h"
#include "basic/OpenSMOKE_AdaptiveGrid.h"
//#include "addons/OpenSMOKE_SensitivityAnalysis_Flame1D.h"
#include "addons/OpenSMOKE_SensitivityAnalysis_Fast_Flame1D.h"
//#include "addons/OpenSMOKE_SensitivityAnalysis_Flame1D_PostProcessor.h"
#include "addons/OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D.h"
#include "addons/OpenSMOKE_PostProcessor.h"
#include "liquid/OpenSMOKE_LiquidSpecies.h"
#include "addons/OpenSMOKE_PolimiSoot.h"

const int OnConstantT	 = -33333; 
const int OffConstantT   = -44444;
const int ResetConstantT = -55555;
const double epsilonDiffusion = 0.01;

//////////////////////////////////////////////////////////////////////
// Static Variables													//
//////////////////////////////////////////////////////////////////////

OpenSMOKE_Flame1D *ptFlame;

void DAE_ODE_Print(BzzVector &x, double t)
{
	ptFlame->DAE_ODE_myPrint(x, t);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								MEMORY ALLOCATION												   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Caronte"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flame1D::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Caronte"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_Flame1D::OpenSMOKE_Flame1D()
{
	name_object			= "[Name not assigned]";
	iUnsteadyFromBackUp = 0;
}

void OpenSMOKE_Flame1D::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Flame1D::Assign(OpenSMOKE_Flame1D_DataManager *_data)
{
	data = _data;
}

void OpenSMOKE_Flame1D::Assign(OpenSMOKE_Flame1D_ScheduleClass *_operations)
{
	operations = _operations;
}

void OpenSMOKE_Flame1D::Assign(OpenSMOKE_ReactingGas *_mix)
{
	mix = _mix;
}

void OpenSMOKE_Flame1D::Assign(OpenSMOKE_GlobalKinetics *_global)
{
	global = _global;
}

void OpenSMOKE_Flame1D::allocate_main_variables_only_Np()
{
	ChangeDimensions(Np, &U);
	ChangeDimensions(Np, &G);
	ChangeDimensions(Np, &H);
	ChangeDimensions(Np, &T);
	ChangeDimensions(Np, NC, &W);
	
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)
		ChangeDimensions(Np, 2*qmom.N, &moments);
		
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)
	{	
		ChangeDimensions(Np, &phiN);
		ChangeDimensions(Np, &phiM);
	}
}

void OpenSMOKE_Flame1D::allocate_only_Np()
{
	// Physical variables
	// --------------------------------------------
	ChangeDimensions(Np, NC, &X);
	ChangeDimensions(Np, NC, &vStar);
	ChangeDimensions(Np, NC, &vFick);
	ChangeDimensions(Np, NC, &vThermophoretic);

	ChangeDimensions(Np, NC, &vStar_e);
	ChangeDimensions(Np, NC, &vStar_w);
	ChangeDimensions(Np, NC, &coeff_e);

	// Corrections
	ChangeDimensions(Np, NC, &W_c);
	ChangeDimensions(Np, NC, &X_c);
	ChangeDimensions(Np,     &MW_c);

	// Elemental
	// --------------------------------------------
	ChangeDimensions(Np, mix->NumberOfElements(), &x_elemental);
	ChangeDimensions(Np, mix->NumberOfElements(), &omega_elemental);
	ChangeDimensions(Np, &Z);

	// Differentials
	// --------------------------------------------
	ChangeDimensions(Np, &dU);
	ChangeDimensions(Np, &dG);
	ChangeDimensions(Np, &dH);
	ChangeDimensions(Np, &dT);
	ChangeDimensions(Np, NC, &dW);

	// Properties
	// --------------------------------------------
	ChangeDimensions(Np, &mu);
	ChangeDimensions(Np, &rho);
	ChangeDimensions(Np, &urho);
	ChangeDimensions(Np, &Cp);
	ChangeDimensions(Np, &lambda);
	ChangeDimensions(Np, &QReaction);
	ChangeDimensions(Np, &PMtot);
	ChangeDimensions(Np, &uPMtot);
	ChangeDimensions(Np, NC, &Dm);
	ChangeDimensions(Np, NC, &R);
	ChangeDimensions(Np, &rhow);
	ChangeDimensions(Np, &rhoe);
	ChangeDimensions(Np, &sumCpDiffusive);
		
	ChangeDimensions(Np, &A_x_lambdaw);
	ChangeDimensions(Np, &A_x_lambdae);

	ChangeDimensions(Np, &A_x_rhow);
	ChangeDimensions(Np, &A_x_rhoe);

	ChangeDimensions(Np, &Qrad);

	ChangeDimensions(Np, &auxNp);
	
	ChangeDimensions(Np, NC, &CpMap);
	ChangeDimensions(Np, NC, &Cpk);
	ChangeDimensions(Np, NC, &lambdaMap);
	ChangeDimensions(Np, NC, &muMap);
	ChangeDimensions(Np, NR, &k1Map);
	ChangeDimensions(Np, NR, &k2Map);
	ChangeDimensions(Np, NR, &uKeqMap);
	ChangeDimensions(Np, NR, &logFcentMap);
	ChangeDimensions(Np, NR, &reactionDSMap);
	ChangeDimensions(Np, NR, &reactionDHMap);


	DjkMap = new BzzMatrix[Np+1];
	for(int i=0;i<=Np;i++)
		ChangeDimensions(NC, NC, &DjkMap[i]);

	if (data->iSoretEffect == true)
	{
		ChangeDimensions(Np, NC, &Teta);
		TetakjMap = new BzzMatrix[Np+1];
		for(int j=0;j<=Np;j++)
			ChangeDimensions(NC, NC, &TetakjMap[j]);
	}


	// Variabili ausiliarie dipendenti solo da NP
	// --------------------------------------------
	ChangeDimensions(Np, &M);
	ChangeDimensions(Np, &diffM);
	ChangeDimensions(Np, &diffT);
	ChangeDimensions(Np, &diffTcentral);
	ChangeDimensions(Np, NC, &diffW);
	ChangeDimensions(Np, &diffWscalar);
	ChangeDimensions(Np, &sumW);
	ChangeDimensions(Np, &vc);
	ChangeDimensions(Np, &vc_e);


	// Variabili per la risoluzione dei sistemi
	// --------------------------------------------
	ChangeDimensions(Np*nBlock, &inDerAlg);
	ChangeDimensions(Np*nBlock, &xMin);
	ChangeDimensions(Np*nBlock, &xMax);

	// QMOM Variables
	// --------------------------------------------
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)
		allocate_QMOM_Np();
		
	// Soot Variables 
	// --------------------------------------------
	
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)
		allocate_SOOT_Np();

	if (data->iAssignedReactionRates == true || data->iAssignedROPA == true)
		ChangeDimensions(Np, mix->NumberOfReactions(), &RR);

	if (data->iSingleContributions == true)
		ChangeDimensions(Np, data->index_SingleContributions.Size()*3, &single_contributions);

	if (data->iCorrectionDiffusivity == true)
	{
		for(int kk=1;kk<=data->correction_diffusivity.Size();kk++)
			mix->CorrectBinaryDiffusionCoefficients(data->index_correction_diffusivity[2*(kk-1)+1], 
													data->index_correction_diffusivity[2*(kk-1)+2], data->correction_diffusivity[kk]);
	
		data->iCorrectionDiffusivity = false;
	}

	if (data->iCorrectionFormationEnthalpy == true)
	{
		for(int kk=1;kk<=data->correction_formation_enthalpy.Size();kk++)
			mix->CorrectFormationEnthalpy(	data->index_correction_formation_enthalpy[kk], 
											data->correction_formation_enthalpy[kk]);
	
		data->iCorrectionFormationEnthalpy = false;
	}
}


void OpenSMOKE_Flame1D::allocate_all()
{
	std::cout << " Allocating main variables (only Np)" << std::endl;
	allocate_main_variables_only_Np();
	std::cout << " Allocating secondary variables (only Np)" << std::endl;
	allocate_only_Np();

	// Variabili ausiliarie dipendenti solo da NC
	// --------------------------------------------
	std::cout << " Allocating secondary variables (only NC)" << std::endl;
	ChangeDimensions(NC, &xVector);
	ChangeDimensions(NC, &wVector);
	ChangeDimensions(NC, &cVector);
	ChangeDimensions(NC, &RVector);
	ChangeDimensions(NC, &DmixVector);
	ChangeDimensions(NC, &WC);
	ChangeDimensions(NC, &WO);
	ChangeDimensions(NC, &BCW_C);
	ChangeDimensions(NC, &BCW_O);

	if (data->iSoretEffect == true)
		ChangeDimensions(NC, &TetaMixVector);

	std::cout << " Initializing ROPA" << std::endl;
	if (data->iVerboseAssignedROPA == true)
		ropa.Initialize(mix, data->index_ROPA);
}

void OpenSMOKE_Flame1D::allocate_QMOM_Np()
{
	ChangeDimensions(Np, 2*qmom.N, &dmoments);
	ChangeDimensions(Np, 2*qmom.N, &moments_source);

	ChangeDimensions(Np, 2*qmom.N, &diffMoments);
	ChangeDimensions(Np, 2*qmom.N, &DiffusionMoments);
	ChangeDimensions(Np, 2*qmom.N, &A_x_rho_x_Diffe);
	ChangeDimensions(Np, 2*qmom.N, &A_x_rho_x_Diffw);
}

void OpenSMOKE_Flame1D::allocate_QMOM_N()
{
	ChangeDimensions(2*qmom.N, &qmom_sources);
	ChangeDimensions(2*qmom.N, &momentsC);
}

void OpenSMOKE_Flame1D::allocate_SOOT_Np()
{
	ChangeDimensions(Np, &dphiN);
	ChangeDimensions(Np, &dphiM);
	ChangeDimensions(Np, &source_phiN);
	ChangeDimensions(Np, &source_phiM);
	ChangeDimensions(Np, NC, &SootGasCorrection);

	ChangeDimensions(Np, &diff_phiN);
	ChangeDimensions(Np, &diff_phiM);
	ChangeDimensions(Np, &DiffusionSoot);
	ChangeDimensions(Np, 1, &A_x_rho_x_Diffe);
	ChangeDimensions(Np, 1, &A_x_rho_x_Diffw);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								UTILITIES										   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::setMinimumAndMaximumValues(const flame1d_model string_kind)
{
	// 1. OPPOSED_ALL
	// 2. OPPOSED_ONLY_MOMENTUM
	// 3. OPPOSED_ONLY_TEMPERATURE
	// 4. OPPOSED_ONLY_MASS_FRACTIONS
	// 5. OPPOSED_NO_ENERGY
	// 6. OPPOSED_NO_MOMENTUM
	// 7. OPPOSED_COLD_REDUCED
	// ----------------------------------
	// 8.  PREMIXED_ALL
	// 9.  PREMIXED_NOENERGY
	// 10. PREMIXED_FLAMESPEED

	int i, j, k;
	double ZERO =  0.;
	double ONE  =  1.;
	double UMAX =  1e32;
	double GMAX =  1e32;
	double HMAX =  1e32;
	double TMIN  =  200.;
	double TMAX  =  9000.;
	double MINFLOWRATE =  1.e-12;	// kg/s/m2
	double MAXFLOWRATE =  1.e+12;	// kg/s/m2

	// -------------------------------------------------------------------------
	// Case 1: OPPOSED: Complete System
	//         Velocity, Mass Flow Rate, EigenValue, Temperature, Species
	// -------------------------------------------------------------------------
	if ((string_kind == OPPOSED_ALL) || (string_kind == TWIN_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = -UMAX;		// Velocity
			xMin[k++] = -GMAX;		// Mass Flow Rate
			xMin[k++] = -HMAX;		// Eigen Value
			xMin[k++] =  TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass fractions
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = UMAX;		// Velocity
			xMax[k++] = GMAX;		// Mass Flow Rate
			xMax[k++] = HMAX;		// EigenValue
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
		}
	}

	// -------------------------------------------------------------------------
	// Case 2: OPPOSED: Just Temperature and Species
	//         Temperature, Species
	// -------------------------------------------------------------------------
	else if ((string_kind == OPPOSED_NO_MOMENTUM) || (string_kind == TWIN_NO_MOMENTUM))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] =  TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass fractions
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
		}
	}

	// -------------------------------------------------------------------------
	// Case 3: OPPOSED: Complete Without Temperature
	//         Velocity, Mass Flow Rate, EigenValue, Species
	// -------------------------------------------------------------------------
	else if ((string_kind == OPPOSED_NO_ENERGY) || (string_kind == TWIN_NO_ENERGY))
	{	
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = -UMAX;		// Velocity
			xMin[k++] = -GMAX;		// Mass Flow Rate
			xMin[k++] = -HMAX;		// Eigen Value
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass fractions
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = UMAX;		// Velocity
			xMax[k++] = GMAX;		// Mass Flow Rate
			xMax[k++] = HMAX;		// EigenValue
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
		}
	}

	// -------------------------------------------------------------------------
	// Case 4: OPPOSED: Only Temperature
	//         Temperature
	// -------------------------------------------------------------------------
	else if ((string_kind == OPPOSED_ONLY_TEMPERATURE) || (string_kind == TWIN_ONLY_TEMPERATURE))
	{	
		k=1;
		for(i=1;i<=Np;i++)
			xMin[k++] = TMIN;		// Temperature

		k=1;
		for(i=1;i<=Np;i++)
			xMax[k++] = TMAX;		// Temperature
	}

	// -------------------------------------------------------------------------
	// Case 5: OPPOSED: Only Momentum
	//         Velocity, Mass Flow rate, EigenValue
	// -------------------------------------------------------------------------
	else if ((string_kind == OPPOSED_ONLY_MOMENTUM) || (string_kind == TWIN_ONLY_MOMENTUM))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = -UMAX;		// Velocity
			xMin[k++] = -GMAX;		// Mass Flow Rate
			xMin[k++] = -HMAX;		// EigenValue
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = UMAX;		// Velocity
			xMax[k++] = GMAX;		// Mass flow rate
			xMax[k++] = HMAX;		// Eigenvalue
		}
	}

	// -------------------------------------------------------------------------
	// Case 6: PREMIXED: Complete System
	//         Temperature and Species
	// -------------------------------------------------------------------------
	else if (string_kind == PREMIXED_ALL)
	{	
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass Fractions
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
		}
	}

	// -------------------------------------------------------------------------
	// Case 7: PREMIXED: Flame Speed Calculations
	//         Temperature, Species, Mass Flow Rate
	// -------------------------------------------------------------------------
	else if (string_kind == PREMIXED_FLAMESPEED)
	{	
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = TMIN;			// Temperature
			xMin[k++] = MINFLOWRATE;	// Mass flow rate
			for(j=1;j<=NC;j++)
				xMin[k++] = 0.;			// Mass fractions
		}
		

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = TMAX;			// Temperature
			xMax[k++] = MAXFLOWRATE;	// Mass flow rate
			for(j=1;j<=NC;j++)
				xMax[k++] = 1.;			// Mass fractions
		}

	}

	else if ((string_kind == OPPOSED_ONLY_MASS_FRACTIONS) || (string_kind == TWIN_ONLY_MASS_FRACTIONS))
	{	
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;		// Mass Fractions

		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;		// Mass Fractions
	}

	else if ((string_kind == OPPOSED_COLD_REDUCED) || (string_kind == TWIN_COLD_REDUCED))
	{	
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=3;j++)

				xMin[k++] = ZERO;		// Mass Fractions

		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=3;j++)
				xMax[k++] = ONE;		// Mass Fractions
	}

	// -------------------------------------------------------------------------
	// Case 8: PREMIXED Case: No Energy
	//         Species
	// -------------------------------------------------------------------------
	else if (   string_kind == PREMIXED_NOENERGY
		     || string_kind == OPPOSED_ONLY_MASS_FRACTIONS
		     || string_kind == OPPOSED_COLD_REDUCED
			 || string_kind == TWIN_ONLY_MASS_FRACTIONS
			 || string_kind == TWIN_COLD_REDUCED )
  {
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;		// Mass Fractions

		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;		// Mass Fractions
	}

	else if (string_kind == PREMIXED_QMOM_ALL)
	{	
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass Fractions
			for(j=1;j<=2*qmom.N;j++)
				xMin[k++] = ZERO;	// Moments
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
			for(j=1;j<=2*qmom.N;j++)
				xMax[k++] = 1e16;	// Moments
		}
	}
	
	else if (string_kind == PREMIXED_QMOM_NOENERGY)
	{	
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass Fractions
			for(j=1;j<=2*qmom.N;j++)
				xMin[k++] = ZERO;	// Moments
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
			for(j=1;j<=2*qmom.N;j++)
				xMax[k++] = 1e16;	// Moments
		}
	}

	else if (string_kind == PREMIXED_SOOT_ALL)
	{	
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass Fractions
			xMin[k++] = ZERO;		// PhiN
			xMin[k++] = ZERO;		// PhiM
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
			xMax[k++] = ONE;		// PhiN
			xMax[k++] = ONE;		// PhiM
		}
	}
	
	else if (string_kind == PREMIXED_SOOT_NOENERGY)
	{	
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass Fractions
			xMin[k++] = ZERO;		// PhiN
			xMin[k++] = ZERO;		// PhiM
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
			xMax[k++] = ONE;		// PhiN
			xMax[k++] = ONE;		// PhiM
		}
	}

	// -------------------------------------------------------------------------
	// Case 11: OPPOSED: Complete System + Two Equation Model
	//          Velocity, Mass Flow Rate, EigenValue, Temperature, Species
	// -------------------------------------------------------------------------
	else if ((string_kind == OPPOSED_SOOT_ALL) || (string_kind == TWIN_SOOT_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = -UMAX;		// Velocity
			xMin[k++] = -GMAX;		// Mass Flow Rate
			xMin[k++] = -HMAX;		// Eigen Value
			xMin[k++] =  TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass fractions
			xMin[k++] = ZERO;		// PhiN
			xMin[k++] = ZERO;		// PhiM
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = UMAX;		// Velocity
			xMax[k++] = GMAX;		// Mass Flow Rate
			xMax[k++] = HMAX;		// EigenValue
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
			xMax[k++] = ONE;		// PhiN
			xMax[k++] = ONE;		// PhiM
		}
	}

	else if ((string_kind == OPPOSED_QMOM_ALL) || (string_kind == TWIN_QMOM_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = -UMAX;		// Velocity
			xMin[k++] = -GMAX;		// Mass Flow Rate
			xMin[k++] = -HMAX;		// Eigen Value
			xMin[k++] =  TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass fractions
			for(j=1;j<=2*qmom.N;j++)
				xMin[k++] = ZERO;	// Moments
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = UMAX;		// Velocity
			xMax[k++] = GMAX;		// Mass Flow Rate
			xMax[k++] = HMAX;		// EigenValue
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
			for(j=1;j<=2*qmom.N;j++)
				xMax[k++] = 1.e16;	// Moments
		}
	}

	else
		ErrorMessage("Error in Minimum/Maximum setup, Wrong Type: " + string_kind);
}

void OpenSMOKE_Flame1D::setMinimumAndMaximumValuesReduced(const flame1d_model string_kind)
{
	// 1. OPPOSED_ALL
	// 2. OPPOSED_ONLY_MOMENTUM
	// 3. OPPOSED_ONLY_TEMPERATURE
	// 4. OPPOSED_ONLY_MASS_FRACTIONS
	// 5. OPPOSED_NO_ENERGY
	// 6. OPPOSED_NO_MOMENTUM
	// 7. OPPOSED_COLD_REDUCED
	// ----------------------------------
	// 8.  PREMIXED_ALL
	// 9.  PREMIXED_NOENERGY
	// 10. PREMIXED_FLAMESPEED

	int i, j, k;
	double ZERO =  0.;
	double ONE  =  1.;
	double UMAX =  1e32;
	double GMAX =  1e32;
	double HMAX =  1e32;
	double TMIN  =  200.;
	double TMAX  =  9000.;
	double MINFLOWRATE =  1.e-12;	// kg/s/m2
	double MAXFLOWRATE =  1.e+12;	// kg/s/m2

	// -------------------------------------------------------------------------
	// Case 1: OPPOSED: Complete System
	//         Velocity, Mass Flow Rate, EigenValue, Temperature, Species
	// -------------------------------------------------------------------------
	if ((string_kind == OPPOSED_ALL) || (string_kind == TWIN_ALL))
	{
		// Minimum
		k=1;
		{
			xMin[k++] = -UMAX;		// Velocity
			xMin[k++] = -GMAX;		// Mass Flow Rate
			xMin[k++] = -HMAX;		// Eigen Value
			xMin[k++] =  TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass fractions
		}
		for(i=2;i<=Ni;i++)
		{
			xMin[k++] = -UMAX;		// Velocity
			xMin[k++] = -GMAX;		// Mass Flow Rate
			xMin[k++] = -HMAX;		// Eigen Value
		}
		{
			xMin[k++] = -UMAX;		// Velocity
			xMin[k++] = -GMAX;		// Mass Flow Rate
			xMin[k++] = -HMAX;		// Eigen Value
			xMin[k++] =  TMIN;		// Temperature
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;	// Mass fractions
		}

		// Maximum
		k=1;
		{
			xMax[k++] = UMAX;		// Velocity
			xMax[k++] = GMAX;		// Mass Flow Rate
			xMax[k++] = HMAX;		// EigenValue
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
		}
		for(i=2;i<=Ni;i++)
		{
			xMax[k++] = UMAX;		// Velocity
			xMax[k++] = GMAX;		// Mass Flow Rate
			xMax[k++] = HMAX;		// EigenValue
		}
		{
			xMax[k++] = UMAX;		// Velocity
			xMax[k++] = GMAX;		// Mass Flow Rate
			xMax[k++] = HMAX;		// EigenValue
			xMax[k++] = TMAX;		// Temperature
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;	// Mass Fractions
		}
	}

	else
		ErrorMessage("Error in Minimum/Maximum setup, Wrong Type: " + string_kind);
}

void OpenSMOKE_Flame1D::setDifferentialAndAlgebraic(const flame1d_model string_kind)
{
	// 1. OPPOSED_ALL
	// 5. OPPOSED_NO_ENERGY
	// 6. OPPOSED_NO_MOMENTUM
	// ----------------------------------
	// 8.  PREMIXED_ALL
	// 9.  PREMIXED_NOENERGY
	// 10. PREMIXED_FLAMESPEED		

	// This options are used only for non linear systems
	// -------------------------------------------------
	// 2. OPPOSED_ONLY_MOMENTUM			
	// 3. OPPOSED_ONLY_TEMPERATURE
	// 4. OPPOSED_ONLY_MASS_FRACTIONS
	// 7. OPPOSED_COLD REDUCED

	int i, j, k;

	int DIFFERENTIAL = 1;
	int ALGEBRAIC = 0;

	if (string_kind == OPPOSED_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		inDerAlg[k++] = ALGEBRAIC;					// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = ALGEBRAIC;				
			inDerAlg[k++] = ALGEBRAIC;				// Algebraic is really better? // TODO
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
		}

		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
	}

	else if (string_kind == TWIN_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		inDerAlg[k++] = ALGEBRAIC;					// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = ALGEBRAIC;				// U
			inDerAlg[k++] = ALGEBRAIC;				// G
			inDerAlg[k++] = ALGEBRAIC;				// H
			inDerAlg[k++] = DIFFERENTIAL;			// T
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;		// W
		}

		inDerAlg[k++] = ALGEBRAIC;				// U
		inDerAlg[k++] = ALGEBRAIC;				// G
		inDerAlg[k++] = ALGEBRAIC;				// H
		inDerAlg[k++] = ALGEBRAIC;			// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;			// W
	}

	else if (string_kind == OPPOSED_NO_MOMENTUM)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
		}

		inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
	}

	else if (string_kind == TWIN_NO_MOMENTUM)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
		}

		inDerAlg[k++] = DIFFERENTIAL;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
	}


	else if (string_kind == OPPOSED_NO_ENERGY)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = ALGEBRAIC;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
		}

		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
	} 

	else if (string_kind == TWIN_NO_ENERGY)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = ALGEBRAIC;				// U
			inDerAlg[k++] = ALGEBRAIC;				// G
			inDerAlg[k++] = ALGEBRAIC;				// H
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;		// W
		}

		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W

	} 

	else if (string_kind == PREMIXED_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;				     	// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;					// W

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
		}

		inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
	}

	else if (string_kind == PREMIXED_FLAMESPEED)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;				     	// T
		inDerAlg[k++] = ALGEBRAIC;						// H
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;					// W

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = ALGEBRAIC;					// H
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
		}

		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;						// H
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
	}

	else if (string_kind == PREMIXED_NOENERGY)
	{
		k=1;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;					// W

		for(i=2;i<=Ni;i++)
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;

		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
	} 

	else if (string_kind == PREMIXED_QMOM_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;				     	// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;					// W
		for(j=1;j<=2*qmom.N;j++)
			inDerAlg[k++] = ALGEBRAIC;					// Moments
			
		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=2*qmom.N;j++)
				inDerAlg[k++] = DIFFERENTIAL;			// Moments
		}

		inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=2*qmom.N;j++)
			inDerAlg[k++] = DIFFERENTIAL;					// Moments
	}
	
	else if (string_kind == PREMIXED_QMOM_NOENERGY)
	{
		k=1;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;					// W
		for(j=1;j<=2*qmom.N;j++)
			inDerAlg[k++] = ALGEBRAIC;					// Moments
			
		for(i=2;i<=Ni;i++)
		{
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=2*qmom.N;j++)
				inDerAlg[k++] = DIFFERENTIAL;			// Moments
		}

		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=2*qmom.N;j++)
			inDerAlg[k++] = DIFFERENTIAL;					// Moments
	}
	
	else if (string_kind == PREMIXED_SOOT_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;				     	// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;					// W
		inDerAlg[k++] = ALGEBRAIC;						// PhiN
		inDerAlg[k++] = ALGEBRAIC;						// PhiM
			
		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = DIFFERENTIAL;				// PhiN
			inDerAlg[k++] = DIFFERENTIAL;				// PhiM
		}

		inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;						// PhiN
		inDerAlg[k++] = ALGEBRAIC;						// PhiM
	}
	
	else if (string_kind == PREMIXED_SOOT_NOENERGY)
	{
		k=1;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;					// W
		inDerAlg[k++] = ALGEBRAIC;						// PhiN
		inDerAlg[k++] = ALGEBRAIC;						// PhiM
			
		for(i=2;i<=Ni;i++)
		{
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = DIFFERENTIAL;				// PhiN
			inDerAlg[k++] = DIFFERENTIAL;				// PhiM
		}

		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;						// PhiN
		inDerAlg[k++] = ALGEBRAIC;						// PhiM
	}

	else if (string_kind == OPPOSED_SOOT_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		inDerAlg[k++] = ALGEBRAIC;					// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W
		inDerAlg[k++] = ALGEBRAIC;					// PhiN
		inDerAlg[k++] = ALGEBRAIC;					// PhiM

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = DIFFERENTIAL;			// PhiN
			inDerAlg[k++] = DIFFERENTIAL;			// PhiM
		}

		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;					// PhiN
		inDerAlg[k++] = ALGEBRAIC;					// PhiM
	}

	else if (string_kind == TWIN_SOOT_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		inDerAlg[k++] = ALGEBRAIC;					// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W
		inDerAlg[k++] = ALGEBRAIC;					// PhiN
		inDerAlg[k++] = ALGEBRAIC;					// PhiM

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = DIFFERENTIAL;			// PhiN
			inDerAlg[k++] = DIFFERENTIAL;			// PhiM
		}

		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = DIFFERENTIAL;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;					// PhiN
		inDerAlg[k++] = ALGEBRAIC;					// PhiM
	}

	else if (string_kind == OPPOSED_QMOM_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		inDerAlg[k++] = ALGEBRAIC;					// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W
		for(j=1;j<=2*qmom.N;j++)
			inDerAlg[k++] = ALGEBRAIC;				// Moments

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=2*qmom.N;j++)
				inDerAlg[k++] = DIFFERENTIAL;			// Moments
		}

		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=2*qmom.N;j++)
			inDerAlg[k++] = ALGEBRAIC;				// Moments
	}

	else if (string_kind == TWIN_QMOM_ALL)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H
		inDerAlg[k++] = ALGEBRAIC;					// T
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;				// W
		for(j=1;j<=2*qmom.N;j++)
			inDerAlg[k++] = DIFFERENTIAL;				// Moments

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = ALGEBRAIC;
			inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;
			for(j=1;j<=2*qmom.N;j++)
				inDerAlg[k++] = DIFFERENTIAL;			// Moments
		}

		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = DIFFERENTIAL;
		for(j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;
		for(j=1;j<=2*qmom.N;j++)
			inDerAlg[k++] = DIFFERENTIAL;				// Moments
	}

	else if (string_kind == OPPOSED_ONLY_MOMENTUM)
	{
		k=1;
		inDerAlg[k++] = ALGEBRAIC;					// U
		inDerAlg[k++] = ALGEBRAIC;					// G
		inDerAlg[k++] = ALGEBRAIC;					// H

		for(i=2;i<=Ni;i++)
		{
			inDerAlg[k++] = ALGEBRAIC;			// Algebraic is really better? // TODO
			inDerAlg[k++] = DIFFERENTIAL;
			inDerAlg[k++] = ALGEBRAIC;
		}

		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
		inDerAlg[k++] = ALGEBRAIC;
	}

	else
		ErrorMessage("Error in Differential/Algebraic setup, Wrong Type: " + string_kind);
}


void OpenSMOKE_Flame1D::recoverPhysicalVariables(const flame1d_model string_kind, BzzVector &x)
{
	// 1. OPPOSED_ALL
	// 2. OPPOSED_ONLY_MOMENTUM
	// 3. OPPOSED_ONLY_TEMPERATURE
	// 4. OPPOSED_ONLY_MASS_FRACTIONS
	// 5. OPPOSED_NO_ENERGY
	// 6. OPPOSED_NO_MOMENTUM
	// 7. OPPOSED_COLD_REDUCED
	// ----------------------------------
	// 8.  PREMIXED_ALL
	// 9.  PREMIXED_NOENERGY
	// 10. PREMIXED_FLAMESPEED

	int i, k, j;


	if ((string_kind == OPPOSED_ALL) || (string_kind == TWIN_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			U[i] = x[k++];
			G[i] = x[k++];
			H[i] = x[k++];
			T[i] = x[k++];
	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];
		}
	}

	else if (string_kind == PREMIXED_NOENERGY	||
		     string_kind == OPPOSED_ONLY_MASS_FRACTIONS || 
			 string_kind == TWIN_ONLY_MASS_FRACTIONS)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];
		}
	}

	else if ((string_kind == OPPOSED_ONLY_MOMENTUM) || (string_kind == TWIN_ONLY_MOMENTUM))
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			U[i] = x[k++];
			G[i] = x[k++];
			H[i] = x[k++];
		}
	}

	else if ((string_kind == OPPOSED_ONLY_TEMPERATURE) || (string_kind == TWIN_ONLY_TEMPERATURE))
	{		
		k=1;
		for(i=1;i<=Np;i++)
			T[i] = x[k++];
	}

	else if ((string_kind == OPPOSED_COLD_REDUCED) || (string_kind == TWIN_COLD_REDUCED))
	{	
		sumW = 1.; W = 0.;
		
		k=1; 
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=data->nCold-1;j++)
			{
				W[i][data->iX[j]] = x[k++];
				sumW[i] -= W[i][data->iX[j]];
			}
			W[i][data->iX[data->nCold]] = x[k++];
		}
	}

	else if ( string_kind == PREMIXED_ALL		|| 
		      string_kind == OPPOSED_NO_MOMENTUM ||
			  string_kind == TWIN_NO_MOMENTUM)
	{		
		k=1; 
		for(i=1;i<=Np;i++)
		{
			T[i] = x[k++];
	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];	
		}
	}

	else if ((string_kind == OPPOSED_NO_ENERGY) || (string_kind == TWIN_NO_ENERGY))
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			U[i] = x[k++];
			G[i] = x[k++];
			H[i] = x[k++];
	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];
		}
	}

	else if (string_kind == PREMIXED_FLAMESPEED)
	{		
		k=1; 
		for(i=1;i<=Np;i++)
		{
			// Temperature
			T[i] = x[k++];

			// Flow Rate
			H[i] = x[k++];
			
			// Mass Fractions
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];
		}
	}

	else if (string_kind == PREMIXED_QMOM_ALL)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			T[i] = x[k++];
	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];

			for(j=1;j<=2*qmom.N;j++)
				moments[i][j] = x[k++];
		}
	}
	
	else if (string_kind == PREMIXED_QMOM_NOENERGY)
	{		
		k=1; 
		for(i=1;i<=Np;i++)
		{	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];

			for(j=1;j<=2*qmom.N;j++)
				moments[i][j] = x[k++];
		}
	}	
	
	else if (string_kind == PREMIXED_SOOT_ALL)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			T[i] = x[k++];
	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];
	
			phiN[i] = x[k++];	// soot particle number density
			phiM[i] = x[k++];	// soot mass fraction
		}
	}
	
	else if (string_kind == PREMIXED_SOOT_NOENERGY)
	{		
		k=1; 
		for(i=1;i<=Np;i++)
		{	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];	

			phiN[i] = x[k++];		// soot particle number density
			phiM[i] = x[k++];		// soot mass fraction
		}
	}

	else if ((string_kind == OPPOSED_SOOT_ALL) || (string_kind == TWIN_SOOT_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			U[i] = x[k++];
			G[i] = x[k++];
			H[i] = x[k++];
			T[i] = x[k++];
	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];

			phiN[i] = x[k++];		// soot particle number density
			phiM[i] = x[k++];		// soot mass fraction
		}
	}

	else if ((string_kind == OPPOSED_QMOM_ALL) || (string_kind == TWIN_QMOM_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			U[i] = x[k++];
			G[i] = x[k++];
			H[i] = x[k++];
			T[i] = x[k++];
	
			for(j=1;j<=NC;j++)
				W[i][j] = x[k++];

			for(j=1;j<=2*qmom.N;j++)
				moments[i][j] = x[k++];
		}
	}

	else if (string_kind == SINGLEREACTOR_ISOTHERMAL)
	{
		k=1;
		for(j=1;j<=NC;j++)
			W[indexReactor][j] = x[k++];
	}
	
	else 
		ErrorMessage("Error in Initial Values, Wrong Type: " + string_kind);
}

void OpenSMOKE_Flame1D::recoverPhysicalVariables_Reduced(const flame1d_model string_kind, BzzVector &x)
{
	if ((string_kind == OPPOSED_ALL) || (string_kind == TWIN_ALL))
	{
		int k=1;
		{
			U[1] = x[k++];
			G[1] = x[k++];
			H[1] = x[k++];
			T[1] = x[k++];
	
			for(int j=1;j<=NC;j++)
				W[1][j] = x[k++];
		}
		for(int i=2;i<=Ni;i++)
		{
			U[i] = x[k++];
			G[i] = x[k++];
			H[i] = x[k++];
		}
		{
			U[Np] = x[k++];
			G[Np] = x[k++];
			H[Np] = x[k++];
			T[Np] = x[k++];
	
			for(int j=1;j<=NC;j++)
				W[Np][j] = x[k++];
		}
	}
	else 
		ErrorMessage("Error in Initial Values, Wrong Type: " + string_kind);
}


void OpenSMOKE_Flame1D::recoverResiduals(const flame1d_model string_kind, BzzVector &f)
{
	// 1. OPPOSED_ALL
	// 2. OPPOSED_ONLY_MOMENTUM
	// 3. OPPOSED_ONLY_TEMPERATURE
	// 4. OPPOSED_ONLY_MASS_FRACTIONS
	// 5. OPPOSED_NO_ENERGY
	// 6. OPPOSED_NO_MOMENTUM
	// 7. OPPOSED_COLD_REDUCED
	// ----------------------------------
	// 8.  PREMIXED_ALL
	// 9.  PREMIXED_NOENERGY
	// 10. PREMIXED_FLAMESPEED

	int i, k, j;


	if ((string_kind == OPPOSED_ALL) || (string_kind == TWIN_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dU[i];
			f[k++] = dG[i];
			f[k++] = dH[i];
			f[k++] = dT[i];
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
		}
	}

	else if ( string_kind == PREMIXED_NOENERGY	||
		      string_kind == OPPOSED_ONLY_MASS_FRACTIONS ||
			  string_kind == TWIN_ONLY_MASS_FRACTIONS) 
	{		
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
	}

	else if ((string_kind == OPPOSED_ONLY_MOMENTUM) || (string_kind == TWIN_ONLY_MOMENTUM))
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dU[i];
			f[k++] = dG[i];
			f[k++] = dH[i];
		}
	}

	else if ((string_kind == OPPOSED_ONLY_TEMPERATURE) || (string_kind == TWIN_ONLY_TEMPERATURE))
	{		
		k=1;
		for(i=1;i<=Np;i++)
			f[k++] = dT[i];
	}

	else if ((string_kind == OPPOSED_COLD_REDUCED) || (string_kind == TWIN_COLD_REDUCED))
	{	
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=data->nCold;j++)
				f[k++] = dW[i][data->iX[j]];
		}
	}

	else if (	string_kind == PREMIXED_ALL	||  
				string_kind == OPPOSED_NO_MOMENTUM ||
				string_kind == TWIN_NO_MOMENTUM)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dT[i];
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
		}
	}

	else if ((string_kind == OPPOSED_NO_ENERGY) || (string_kind == TWIN_NO_ENERGY))
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dU[i];
			f[k++] = dG[i];
			f[k++] = dH[i];
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
		}
	}

	else if (string_kind == PREMIXED_FLAMESPEED)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dT[i];
			f[k++] = dH[i];
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
		}
	}

	else if (string_kind == PREMIXED_QMOM_ALL)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dT[i];
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
			for(j=1;j<=2*qmom.N;j++)
				f[k++] = dmoments[i][j];
		}
	}
	
	else if (string_kind == PREMIXED_QMOM_NOENERGY)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
			for(j=1;j<=2*qmom.N;j++)
				f[k++] = dmoments[i][j];
		}
	}
	
	else if (string_kind == PREMIXED_SOOT_ALL)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dT[i];
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
			f[k++] = dphiN[i];			// soot particle number density
			f[k++] = dphiM[i];			// soot mass fraction
		}
	}
	
	else if (string_kind == PREMIXED_SOOT_NOENERGY)
	{		
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
			f[k++] = dphiN[i];			// soot particle number density
			f[k++] = dphiM[i];			// soot mass fraction
		}
	}	

	else if ((string_kind == OPPOSED_SOOT_ALL) || (string_kind == TWIN_SOOT_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dU[i];
			f[k++] = dG[i];
			f[k++] = dH[i];
			f[k++] = dT[i];
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
			f[k++] = dphiN[i];			// soot particle number density
			f[k++] = dphiM[i];			// soot mass fraction
		}
	}

	else if ((string_kind == OPPOSED_QMOM_ALL) || (string_kind == TWIN_QMOM_ALL))
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			f[k++] = dU[i];
			f[k++] = dG[i];
			f[k++] = dH[i];
			f[k++] = dT[i];
			for(j=1;j<=NC;j++)
				f[k++] = dW[i][j];
			for(j=1;j<=2*qmom.N;j++)
				f[k++] = dmoments[i][j];
		}
	}

	else if (string_kind == SINGLEREACTOR_ISOTHERMAL)
	{
		k=1;
		for(j=1;j<=NC;j++)
			f[k++] = W[indexReactor][j];
	}
	
	else
		ErrorMessage("Error in Recover Residuals, Wrong Type: " + string_kind);
}

void OpenSMOKE_Flame1D::recoverResiduals_Reduced(const flame1d_model string_kind, BzzVector &f)
{
	if ((string_kind == OPPOSED_ALL) || (string_kind == TWIN_ALL))
	{
		int k=1;
		{
			f[k++] = dU[1];
			f[k++] = dG[1];
			f[k++] = dH[1];
			f[k++] = dT[1];
			for(int j=1;j<=NC;j++)
				f[k++] = dW[1][j];
		}
		for(int i=2;i<=Ni;i++)
		{
			f[k++] = dU[i];
			f[k++] = dG[i];
			f[k++] = dH[i];
		}
		{
			f[k++] = dU[Np];
			f[k++] = dG[Np];
			f[k++] = dH[Np];
			f[k++] = dT[Np];
			for(int j=1;j<=NC;j++)
				f[k++] = dW[Np][j];
		}
	}
	else
		ErrorMessage("Error in Recover Residuals (Reduced), Wrong Type: " + string_kind);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								SUMMARIES											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::printSummaryForOpposed()
{
	int i;

	cout << "-------------------------------------------------------" << endl;
	cout << "              COUNTER DIFFUSION FLAME                  " << endl;
	cout << "-------------------------------------------------------" << endl;
	cout << endl;

	cout << "FUEL SIDE" << endl;
	cout << "--------------------------------" << endl;
	cout << " Velocity:    "<< data->VC*100. << " cm/s" << endl;
	cout << " Temperature: "<< data->TC << " K" << endl;
	cout << " Density:     "<< rho[1] << " kg/m3" << endl;
	cout << " Composition: ";
	for (i=1;i<=data->xC.Size();i++)
		cout  << setw(20) << left << data->nameC[i] << data->xC[i] << endl << "              ";
	cout << endl;

	cout << "OXYDIZER SIDE" << endl;
	cout << "--------------------------------" << endl;
	cout << " Velocity:    "<< data->VO*100. << " cm/s" << endl;
	cout << " Temperature: "<< data->TO << " K" << endl;
	cout << " Density:     "<< rho[Np] << " kg/m3" << endl;
	cout << " Composition: ";
	for (i=1;i<=data->xO.Size();i++)
		cout  << setw(20) << left << data->nameO[i] << data->xO[i] << endl << "              ";
	cout << endl;

	cout << "FLAME PROPERTIES" << endl;
	cout << "--------------------------------" << endl;
	cout << " Inlet distance:     " << grid.L*100.			<< " cm" << endl;
	cout << " Stagnation plane:   " << xST*100.				<< " cm" << endl;
	cout << " Strain rate:        " << K					<< " 1/s" << endl;
	if (data->iPoolFire != POOL_FIRE_NONE)
	cout << " Strain rate (P.F.): " << 2.*data->VO/grid.L	<< " 1/s" << endl;

	cout << endl;

	cout << "SIMULATION" << endl;
	cout << "--------------------------------" << endl;
	cout << " Number of points:   " << Np << endl;
	cout << " G derivative:       " << data->iDerG << endl;
	cout << " T derivative:       " << data->iDerT << endl;
	cout << " W derivative:       " << data->iDerW << endl;
	cout << " Soret effect:       " << data->iSoretEffect << endl;
	cout << " Thermophoresis:     " << data->iThermophoreticEffect << endl;
	cout << endl << endl;

}

void OpenSMOKE_Flame1D::printSummaryForTwin()
{
	int i;

	cout << "-------------------------------------------------------" << endl;
	cout << "                 TWIN DIFFUSION FLAME                  " << endl;
	cout << "-------------------------------------------------------" << endl;
	cout << endl;

	cout << "INLET" << endl;
	cout << "--------------------------------" << endl;
	cout << " Velocity:    " << data->VC*100. << " cm/s" << endl;
	cout << " Temperature: " << data->TC << " K" << endl;
	cout << " Density:     " << rho[1] << " kg/m3" << endl;
	cout << " Composition: ";
	for (i=1;i<=data->xC.Size();i++)
		cout  << setw(20) << left << data->nameC[i] << data->xC[i] << endl << "              ";
	cout << endl;

	cout << "FLAME PROPERTIES" << endl;
	cout << "--------------------------------" << endl;
	cout << " Inlet distance: " << grid.L*100. << " cm" << endl;
	cout << " Strain rate:    "	<< K << " 1/s" << endl;
	cout << endl;

	cout << "SIMULATION" << endl;
	cout << "--------------------------------" << endl;
	cout << " Number of points: " << Np << endl;
	cout << " G derivative:     " << data->iDerG << endl;
	cout << " T derivative:     " << data->iDerT << endl;
	cout << " W derivative:     " << data->iDerW << endl;
	cout << endl << endl;

}

void OpenSMOKE_Flame1D::printSummaryForPremixed()
{
	int i;

	cout << "-------------------------------------------------------" << endl;
	cout << "              PREMIXED FLAME                  " << endl;
	cout << "-------------------------------------------------------" << endl;
	cout << endl;

	cout << "REACTANTS" << endl;
	cout << "--------------------------------" << endl;
	cout << " MassFlowRate: " << data->MassFlowRate*1.e3/(data->CrossSectionalArea*1.e4) << " g/s/cm2" << endl;
	cout << " Temperature:  " << data->TC << " K" << endl;
	cout << " Density:      " << rho[1] << " kg/m3" << endl;
	cout << " Composition:  ";
	for (i=1;i<=data->xC.Size();i++)
		cout  << setw(20) << left << data->nameC[i] << data->xC[i] << endl << "               ";
	cout << endl;

	cout << "SIMULATION" << endl;
	cout << "--------------------------------" << endl;
	cout << " Length:           " << grid.L*100. << " cm" << endl;
	cout << " Number of points: " << Np << endl;
	cout << " T derivative:     " << data->iDerT << endl;
	cout << " W derivative:     " << data->iDerW << endl;
	cout << endl << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								PROPERTIES											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenSMOKE_Flame1D::properties(int ReactionON, int jacobianIndex, BzzVectorInt &jacobianVariables, int dimBlock, int indexT)
{
	// jacobianIndex = -2 || -1 : valutazione normale delle proprieta
	// jacobianIndex =  0 : memorizzazione delle proprieta
	
	int i, j;
		
	int choice;
	BzzVectorInt pointToUpdate;

	// In this case the temperature dependencies have been solved and saved; 
	// therefore only need to recover all this information
	if (jacobianIndex == OnConstantT) 
		choice = 3;

	// In this case we don't solve the energy equation and therefore the temperature
	// profile is assigned and constant. However the first time all the dependencies
	// on T must be calculated
	else if (jacobianIndex == OffConstantT) 
	{
		choice = 2;
		tagConstantT = OnConstantT;
	}

	// If none of the two previous conditions is satisfied, it means that the energy
	// balance is necessary
	else
	{
		if (jacobianIndex < 0) choice = 1;			//	Calcolo di tutte le proprieta senza memorizzazione
		else if (jacobianIndex == 0) choice = 2;	//	Calcolo di tutte le proprieta con memorizzazione
		else //if (jacobianIndex > 0)
		{
			choice = 3;
			for(j=1;j<=jacobianVariables.Size();j++)
				if ( (jacobianVariables[j]%(dimBlock)) == indexT)
					pointToUpdate.Append(int(jacobianVariables[j]/dimBlock)+1);
		}
	}

	double cTot;

	molarFractionsAndPMtot();

	if ( data->iGasRadiation == true || data->iRadiativeSootModel != RADIATIVE_SOOT_MODEL_NONE )
		calculate_radiation();

	// --------------------------------------------------------------------------
	// PROPERTIES FOR DIFFERENT T 
	// --------------------------------------------------------------------------
	if (choice == 1)
	{
		for(i=1;i<=Np;i++)
		{
			// Estrazioni dei vettori delle omega e delle x
			W.GetRow(i,&wVector);
			X.GetRow(i,&xVector);

			// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
			cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			rho[i] = cTot * PMtot[i];
			urho[i] = 1./rho[i];
			cVector = cTot*xVector;


			// b. Calcolo dei calori specifici [J/kgK]
			mix->SpeciesCp(T[i]);
			Cp[i] = mix->MixCp_FromMassFractions(wVector);	
			Cpk.SetRow(i, mix->Cp);

			// c. Calcolo della conducibilitita termica [W/mK]
			mix->SpeciesConductivityFromFitting(T[i]);
			lambda[i] = mix->MixConductivity_FromMolarFractions(xVector);

			// c. Calcolo della viscosita dinamica [Pa.s]
			mix->SpeciesViscosityFromFitting(T[i]);
			mu[i] = mix->MixViscosity_FromMolarFractions(xVector);

			// d. Calcolo dei coefficienti di diffusione
			mix->SpeciesDiffusivityFromFitting(T[i], data->P_bar);
			mix->MixDiffusionCoefficient_FromMolarFractionsAndPM(DmixVector,xVector);
			Dm.SetRow(i, DmixVector);
			
			if (data->iTurbulentDiffusivity == true)
				Dm.SetRow(i, DmixVector[data->jINERT]*data->diffusivityEnhancingFactor);

			if (data->iSoretEffect == true)
			{
				mix->SpeciesThermalDiffusionRatiosFromFitting(T[i]);
				mix->MixThermalDiffusionRatios(TetaMixVector, xVector);
				Teta.SetRow(i, TetaMixVector);
			}

			// e. Reazioni chimiche
			if (ReactionON)
			{
				mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);
				mix->ComputeFromConcentrations( T[i], cVector, cTot, &RVector);// [kmol/m3/s]
				ElementByElementProduct(RVector, mix->M, &RVector);
				R.SetRow(i, RVector);			// [kg/m3/s]
				QReaction[i] = - mix->ComputeQReaction(T[i]);	// [J/m3/s]
			}

			if (data->iFocusReactionRates == true)	RR.SetRow(i, mix->r);

			if (data->iFocusElemetFluxAnalysis == true)
			{
				mix->ComputeFromConcentrations(T[i], cVector, cTot, rForwardVector, rBackwardVector);
				rForward.SetRow(i, rForwardVector);
				rBackward.SetRow(i, rBackwardVector);
			}
		}
	}

	// --------------------------------------------------------------------------
	// PROPERTIES FOR CONSTANT T + MEMORIZATION
	// --------------------------------------------------------------------------
	else if(choice==2)
	{
		// Table Creation (Memorization)
		// ----------------------------------------------------------
		for(i=1;i<=Np;i++)
		{
			// b. Calcolo dei calori specifici [J/kgK]
			mix->SpeciesCp(T[i]);
			CpMap.SetRow(i, mix->Cp);
			
			// c. Calcolo della conducibilitita termica [W/mK]
			mix->SpeciesConductivityFromFitting(T[i]);
			lambdaMap.SetRow(i, mix->lambda);

			// c. Calcolo della viscosita dinamica [Pa.s]
			mix->SpeciesViscosityFromFitting(T[i]);
			muMap.SetRow(i, mix->eta);

			// d. Calcolo dei coefficienti di diffusione
			mix->SpeciesDiffusivityFromFitting(T[i], data->P_bar);
			DjkMap[i] = mix->Djk;

			// d. Calcolo dei coefficienti di diffusione
			mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);
			k1Map.SetRow(i,mix->k1);
			k2Map.SetRow(i,mix->k2);
			uKeqMap.SetRow(i,mix->uKeq);
			logFcentMap.SetRow(i,mix->logFcent);
			reactionDHMap.SetRow(i,mix->reactionDH);
			reactionDSMap.SetRow(i,mix->reactionDS);

			if (data->iSoretEffect == true)
			{
				mix->SpeciesThermalDiffusionRatiosFromFitting(T[i]);
				TetakjMap[i] = mix->Tetakj;
			}
		}

		choice = 3;
	}


	if(choice == 3) // Da chiamare solo se indexJacobian>0 e la T non modificata
	{
		for(i=1;i<=Np;i++)
		{
			// Estrazioni dei vettori delle omega e delle x
			W.GetRow(i,&wVector);
			X.GetRow(i,&xVector);

			// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
			cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			rho[i] = cTot * PMtot[i];
			urho[i] = 1./rho[i];
			cVector = cTot*xVector;

			// b. Calcolo dei calori specifici [J/kgK]
			CpMap.GetRow(i, &mix->Cp);
			Cpk.SetRow(i, CpMap.GetRow(i));
			Cp[i] = mix->MixCp_FromMassFractions(wVector);	
			
			// c. Calcolo della conducibilitita termica [W/mK]
			lambdaMap.GetRow(i, &mix->lambda);
			lambda[i] = mix->MixConductivity_FromMolarFractions(xVector);

			// c. Calcolo della viscosita dinamica [Pa.s]
			muMap.GetRow(i, &mix->eta);
			mu[i] = mix->MixViscosity_FromMolarFractions(xVector);

			// d. Calcolo dei coefficienti di diffusione
			mix->Djk = DjkMap[i];
			mix->MixDiffusionCoefficient_FromMolarFractionsAndPM(DmixVector,xVector);
			Dm.SetRow(i, DmixVector);

			if (data->iTurbulentDiffusivity == true)
				Dm.SetRow(i, DmixVector[data->jINERT]*data->diffusivityEnhancingFactor);

			if (data->iSoretEffect == true)
			{
				mix->Tetakj = TetakjMap[i];
				mix->MixThermalDiffusionRatios(TetaMixVector, xVector);
				Teta.SetRow(i, TetaMixVector);
			}

			// e. Reazioni chimiche
			if (ReactionON)
			{
				k1Map.GetRow(i,&mix->k1);
				k2Map.GetRow(i,&mix->k2);
				uKeqMap.GetRow(i,&mix->uKeq);
				logFcentMap.GetRow(i,&mix->logFcent);
				reactionDHMap.GetRow(i,&mix->reactionDH);
				reactionDSMap.GetRow(i,&mix->reactionDS);
		
				mix->ComputeFromConcentrations( T[i], cVector, cTot, &RVector); // [kmol/m3/s]
				ElementByElementProduct(RVector, mix->M, &RVector);
				R.SetRow(i, RVector);									 // [kg/m3/s]
				QReaction[i] = - mix->ComputeQReaction(T[i]);			 // [J/m3.s]	
			}

			if (data->iFocusReactionRates == true)	RR.SetRow(i, mix->r);
		}
	}
	
	for(j=1;j<=pointToUpdate.Size();j++)
		{
			i = pointToUpdate[j];

			// Estrazioni dei vettori delle omega e delle x
			W.GetRow(i,&wVector);
			X.GetRow(i,&xVector);

			// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
			cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			rho[i] = cTot * PMtot[i];
			urho[i] = 1./rho[i];
			cVector = cTot*xVector;

			// b. Calcolo dei calori specifici [J/kgK]
			mix->SpeciesCp(T[i]);
			Cp[i] = mix->MixCp_FromMassFractions(wVector);	
			Cpk.SetRow(i, mix->Cp);

			// c. Calcolo della conducibilitita termica [W/mK]
			mix->SpeciesConductivityFromFitting(T[i]);
			lambda[i] = mix->MixConductivity_FromMolarFractions(xVector);

			// c. Calcolo della viscosita dinamica [Pa.s]
			mix->SpeciesViscosityFromFitting(T[i]);
			mu[i] = mix->MixViscosity_FromMolarFractions(xVector);

			// d. Calcolo dei coefficienti di diffusione
			mix->SpeciesDiffusivityFromFitting(T[i], data->P_bar);
			mix->MixDiffusionCoefficient_FromMolarFractionsAndPM(DmixVector,xVector);
			Dm.SetRow(i, DmixVector);

			if (data->iTurbulentDiffusivity == true)
				Dm.SetRow(i, DmixVector[data->jINERT]*data->diffusivityEnhancingFactor);

			if (data->iSoretEffect == true)
			{
				mix->SpeciesThermalDiffusionRatiosFromFitting(T[i]);
				mix->MixThermalDiffusionRatios(TetaMixVector,xVector);
				Teta.SetRow(i, TetaMixVector);
			}

			// e. Reazioni chimiche
			if (ReactionON)
			{
				mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);
				mix->ComputeFromConcentrations( T[i], cVector, cTot, &RVector);// [kmol/m3/s]
				ElementByElementProduct(RVector, mix->M, &RVector);
				R.SetRow(i, RVector);			// [kg/m3/s]
				QReaction[i] = - mix->ComputeQReaction(T[i]);	// [J/m3.s]
			}

			if (data->iFocusReactionRates == true)	RR.SetRow(i, mix->r);
		}

		// Global Kinetics
	if (data->iGlobalKinetics == 1)
	{
		for(i=1;i<=Np;i++)
		{
			// Estrazioni dei vettori delle omega e delle x
			W.GetRow(i,&wVector);
			X.GetRow(i,&xVector);

			// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
			cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			rho[i] = cTot * PMtot[i];
			urho[i] = 1./rho[i];
			cVector = cTot*xVector;

			global->GiveMeFormationRates(T[i],cVector,RVector);
			global->GiveMeReactionHeat(T[i], RVector, QReaction[i]);
			R.SetRow(i, RVector);			// [kg/m3/s]
		}
	}

	if (data->iCorrectionReactionRates == true)
	{ 
		QReaction	*= data->correctionReactionRates;
		R			*= data->correctionReactionRates;
		if (data->iFocusReactionRates == true)	
			RR *= data->correctionReactionRates;
	}

	if (data->iFakeTemperatureThermalConductivity == true)
		for(i=1;i<=Np;i++)
		{
			X.GetRow(i,&xVector);
			mix->SpeciesConductivityFromFitting(T[i]+data->fakeTemperatureThermalConductivityIncrement);
			lambda[i] = mix->MixConductivity_FromMolarFractions(xVector);
		}

	if (data->iUnityLewisNumbers == true)
	{
		for(i=1;i<=Np;i++)
			Dm.SetRow(i, lambda[i]/rho[i]/Cp[i]);
	}
	else if (data->iUserDefinedLewisNumbers == true)
	{
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				Dm[i][j] = lambda[i]/rho[i]/Cp[i] / data->user_defined_lewis_numbers[j];
	}

	if (data->iPhysicalSootDiffusionCoefficients == true)
	{
		unsigned int jBIN1A = 0;
		double MWBIN1A = 0.;

		for (int j = 1; j <= NC; j++)
		{
			if (mix->names[j] == "BIN1A")
			{
				jBIN1A = j;
				MWBIN1A = mix->M[j];
				break;
			}
		}

		if (jBIN1A > 0)
		{
			for (int j = 1; j <= NC; j++)
			{
				if (mix->names[j].compare(0, 3, "BIN") == 0)
				{
					const double MWratio = mix->M[j] / MWBIN1A;
					const double correctionCoefficient = pow(MWratio, -0.681);
					for (i = 1; i <= Np; i++)
						Dm[i][j] = Dm[i][jBIN1A] * correctionCoefficient;
				}
			}
		}
	}
}

void OpenSMOKE_Flame1D::properties(int ReactionON)
{
	molarFractionsAndPMtot();

	if ( data->iGasRadiation == true || data->iRadiativeSootModel != RADIATIVE_SOOT_MODEL_NONE )
		calculate_radiation();

	for(int i=1;i<=Np;i++)
	{
		// Estrazioni dei vettori delle omega e delle x
		W.GetRow(i,&wVector);
		X.GetRow(i,&xVector);

		// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		const double cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
		rho[i] = cTot * PMtot[i];
		urho[i] = 1./rho[i];
		cVector = cTot*xVector;

		// c. Calcolo della viscosita dinamica [Pa.s]
		mix->SpeciesViscosityFromFitting(T[i]);
		mu[i] = mix->MixViscosity_FromMolarFractions(xVector);

		if (i==1 || i==2 || i==Ni || i==Np)
		{
			// b. Calcolo dei calori specifici [J/kgK]
			mix->SpeciesCp(T[i]);
			Cp[i] = mix->MixCp_FromMassFractions(wVector);	
			Cpk.SetRow(i, mix->Cp);

			// c. Calcolo della conducibilitita termica [W/mK]
			mix->SpeciesConductivityFromFitting(T[i]);
			lambda[i] = mix->MixConductivity_FromMolarFractions(xVector);

			// d. Calcolo dei coefficienti di diffusione
			mix->SpeciesDiffusivityFromFitting(T[i], data->P_bar);
			mix->MixDiffusionCoefficient_FromMolarFractionsAndPM(DmixVector,xVector);
			Dm.SetRow(i, DmixVector);
			
			if (data->iTurbulentDiffusivity == true)
				Dm.SetRow(i, DmixVector[data->jINERT]*data->diffusivityEnhancingFactor);

			if (data->iSoretEffect == true)
			{
				mix->SpeciesThermalDiffusionRatiosFromFitting(T[i]);
				mix->MixThermalDiffusionRatios(TetaMixVector, xVector);
				Teta.SetRow(i, TetaMixVector);
			}

			// e. Reazioni chimiche
			if (ReactionON)
			{
				mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);
				mix->ComputeFromConcentrations( T[i], cVector, cTot, &RVector);// [kmol/m3/s]
				ElementByElementProduct(RVector, mix->M, &RVector);
				R.SetRow(i, RVector);			// [kg/m3/s]
				QReaction[i] = - mix->ComputeQReaction(T[i]);	// [J/m3/s]
			}

			if (data->iFocusReactionRates == true)	RR.SetRow(i, mix->r);

			if (data->iFocusElemetFluxAnalysis == true)
			{
				mix->ComputeFromConcentrations(T[i], cVector, cTot, rForwardVector, rBackwardVector);
				rForward.SetRow(i, rForwardVector);
				rBackward.SetRow(i, rBackwardVector);
			}
		}
	}

	// Additional

	if (data->iCorrectionReactionRates == true)
	{ 
		QReaction	*= data->correctionReactionRates;
		R			*= data->correctionReactionRates;
		if (data->iFocusReactionRates == true)	
			RR *= data->correctionReactionRates;
	}

	if (data->iFakeTemperatureThermalConductivity == true)
		for(int i=1;i<=Np;i++)
		{
			X.GetRow(i,&xVector);
			mix->SpeciesConductivityFromFitting(T[i]+data->fakeTemperatureThermalConductivityIncrement);
			lambda[i] = mix->MixConductivity_FromMolarFractions(xVector);
		}

	if (data->iUnityLewisNumbers == true)
	{
		for(int i=1;i<=Np;i++)
			Dm.SetRow(i, lambda[i]/rho[i]/Cp[i]);
	}
	else if (data->iUserDefinedLewisNumbers == true)
	{
		for(int i=1;i<=Np;i++)
			for(int j=1;j<=NC;j++)
				Dm[i][j] = lambda[i]/rho[i]/Cp[i] / data->user_defined_lewis_numbers[j];
	}
}

void OpenSMOKE_Flame1D::molarFractionsAndPMtot()
{
	int i, j;

	for(i=1;i<=Np;i++)
	{
		uPMtot[i] = 0.;
		for(j=1;j<=NC;j++)
			uPMtot[i] += W[i][j]*mix->uM[j];
		PMtot[i] = 1./uPMtot[i];
		for(j=1;j<=NC;j++)
			X[i][j] = W[i][j] * PMtot[i] * mix->uM[j];
	}
}

void OpenSMOKE_Flame1D::massFractionsAndPMtot()
{
	int i, j;
	for(i=1;i<=Np;i++)
	{
		PMtot[i] = 0.;
		for(j=1;j<=NC;j++)
			PMtot[i] += X[i][j]*mix->M[j];
		uPMtot[i] = 1./PMtot[i];
		for(j=1;j<=NC;j++)
			W[i][j] = X[i][j] * uPMtot[i] * mix->M[j];
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								DIFFUSION VELOCITIES											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void OpenSMOKE_Flame1D::compute_vStarOpposed()
{
	int i, j;

	vc = 0.;

	// No Soret Effect
	if (data->iSoretEffect == false)
	{
		for(i=1;i<=Ni;i++)
			for(j=1;j<=NC;j++)
			{
				vStar[i][j] = -Dm[i][j] * mix->M[j]*uPMtot[i] * (X[i+1][j]-X[i][j])*grid.udxe[i];
				vc[i] -= vStar[i][j];
			}
			
			for(j=1;j<=NC;j++)
			{
				vStar[Np][j] = -Dm[Np][j] * mix->M[j]*uPMtot[Np] * (X[Np][j]-X[Np-1][j])*grid.udxw[Np];
				vc[Np] -= vStar[Np][j];
			}
	}

	else	// Soret effect: inclusion of thermal diffusion velocities
	{
		for(i=1;i<=Ni;i++)
			for(j=1;j<=NC;j++)
			{
				vStar[i][j] = -Dm[i][j] * mix->M[j]*uPMtot[i] *( (X[i+1][j]-X[i][j])*grid.udxe[i] + Teta[i][j]/T[i]*(T[i+1]-T[i])*grid.udxe[i] );
				vc[i] -= vStar[i][j];
			}
			
			for(j=1;j<=NC;j++)
			{
				vStar[Np][j] = -Dm[Np][j]*mix->M[j]*uPMtot[Np]*( (X[Np][j]-X[Np-1][j])*grid.udxw[Np] + Teta[Np][j]/T[Np]*(T[Np]-T[Np-1])*grid.udxw[Np] );
				vc[Np] -= vStar[Np][j];
			}
	}


	// Correzione delle velocita diffusive
	for(i=1;i<=Np;i++)
		for(j=1;j<=NC;j++)
			vStar[i][j]+=vc[i]*W[i][j];

	// TODO
	for(i=1;i<=Np;i++)
	{
		double sum=0.;
		for(j=1;j<=NC;j++)
			sum += vStar[i][j];
		if (sum >1e-10 || sum <-1e-10)
		{	
	//		cout << "Sum: " << sum << " point: " << i << " x: " << grid.x[i] << " sumX: " << X.GetRow(i).GetSumElements() << " sum W: " << W.GetRow(i).GetSumElements() << endl;
	//		getchar();
		}
	}
}
*/

void OpenSMOKE_Flame1D::compute_vStarOpposed()
{
	// The diffusion contribution are calculated as: omega*V = -D*omega/X*grad(X), i.e. omega*V = -D*MW/MWmix*grad(X)
	// This is the correct formulation, but it is different from whatwe use in OpenFOAM
	if (data->iCorrectDiffusionFormulation == true)
	{
		int i, j;

		// Fake mass and mole fractions
		{
			BzzVector aux_c(Np);
			for (i = 1; i <= Np; i++)
			{
				aux_c = W.GetRow(i);
				double sum = aux_c.GetSumElements();
				for (j = 1; j <= NC; j++)
					W_c[i][j] = W[i][j] / sum;
			}
			mix->GetMWAndMoleFractionsFromMassFractions(MW_c, X_c, W_c);
		}

		// Correction velocity
		vc = 0.;

		bool correction_velocity_only_fick = true;

		// Correction applied only to Fick contribution
		if (correction_velocity_only_fick == true)
		{
			// No Soret Effect
			if (data->iSoretEffect == false)
			{
				if (data->iThermophoreticEffect == false)
				{
					for (i = 1; i <= Ni; i++)
					for (j = 1; j <= NC; j++)
					{
						vFick[i][j] = -Dm[i][j] * mix->M[j] / MW_c[i] * (X_c[i + 1][j] - X_c[i][j])*grid.udxe[i];
						vc[i] -= vFick[i][j];
					}

					for (j = 1; j <= NC; j++)
					{
						vFick[Np][j] = -Dm[Np][j] * mix->M[j] / MW_c[Np] * (X_c[Np][j] - X_c[Np - 1][j])*grid.udxw[Np];
						vc[Np] -= vFick[Np][j];
					}
				}
				else
				{
					BzzVector dT_over_dx(Np);
					grid.FirstDerivative('B', T, dT_over_dx);

					for (i = 1; i <= Ni; i++)
					{

						const double Vt = -0.55*mu[i] / rho[i] / T[i] * dT_over_dx[i];
						for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
							vThermophoretic[i][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[i][mix->polimiSoot->bin_indices()[j]];

						for (j = 1; j <= NC; j++)
						{
							vFick[i][j] = -Dm[i][j] * mix->M[j] / MW_c[i] * (X_c[i + 1][j] - X_c[i][j])*grid.udxe[i];
							vc[i] -= vFick[i][j];
						}
					}

					const double Vt = -0.55*mu[Np] / rho[Np] / T[Np] * dT_over_dx[Np];
					for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
						vThermophoretic[Np][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[Np][mix->polimiSoot->bin_indices()[j]];;

					for (j = 1; j <= NC; j++)
					{
						vFick[Np][j] = -Dm[Np][j] * mix->M[j] / MW_c[Np] * (X_c[Np][j] - X_c[Np - 1][j])*grid.udxw[Np];
						vc[Np] -= vFick[Np][j];
					}
				}
			}

			else	// Soret effect: inclusion of thermal diffusion velocities
			{
				if (data->iThermophoreticEffect == false)
				{
					for (i = 1; i <= Ni; i++)
					for (j = 1; j <= NC; j++)
					{
						vFick[i][j] = -Dm[i][j] * mix->M[j] / MW_c[i] * ((X_c[i + 1][j] - X_c[i][j])*grid.udxe[i] + Teta[i][j] / T[i] * (T[i + 1] - T[i])*grid.udxe[i]);
						vc[i] -= vFick[i][j];
					}

					for (j = 1; j <= NC; j++)
					{
						vFick[Np][j] = -Dm[Np][j] * mix->M[j] / MW_c[Np] * ((X_c[Np][j] - X_c[Np - 1][j])*grid.udxw[Np] + Teta[Np][j] / T[Np] * (T[Np] - T[Np - 1])*grid.udxw[Np]);
						vc[Np] -= vFick[Np][j];
					}
				}
				else
				{
					BzzVector dT_over_dx(Np);
					grid.FirstDerivative('B', T, dT_over_dx);

					for (i = 1; i <= Ni; i++)
					{
						const double Vt = -0.55*mu[i] / rho[i] / T[i] * dT_over_dx[i];
						for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
							vThermophoretic[i][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[i][mix->polimiSoot->bin_indices()[j]];;

						for (j = 1; j <= NC; j++)
						{
							vFick[i][j] = -Dm[i][j] * mix->M[j] / MW_c[i] * ((X_c[i + 1][j] - X_c[i][j])*grid.udxe[i] + Teta[i][j] / T[i] * (T[i + 1] - T[i])*grid.udxe[i]);
							vc[i] -= vFick[i][j];
						}
					}

					const double Vt = -0.55*mu[Np] / rho[Np] / T[Np] * dT_over_dx[Np];
					for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
						vThermophoretic[Np][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[Np][mix->polimiSoot->bin_indices()[j]];

					for (j = 1; j <= NC; j++)
					{
						vFick[Np][j] = -Dm[Np][j] * mix->M[j] / MW_c[Np] * ((X_c[Np][j] - X_c[Np - 1][j])*grid.udxw[Np] + Teta[Np][j] / T[Np] * (T[Np] - T[Np - 1])*grid.udxw[Np]);
						vc[Np] -= vFick[Np][j];
					}
				}
			}


			// Correction applied only to Fick contribution
			{
				for (i = 1; i <= Np; i++)
				for (j = 1; j <= NC; j++)
					vFick[i][j] += vc[i] * W_c[i][j];

				for (i = 1; i <= Np; i++)
				for (j = 1; j <= NC; j++)
					vStar[i][j] = vFick[i][j] + vThermophoretic[i][j];
			}
		}

		// Correction applied to both the contributions
		else
		{
			// No Soret Effect
			if (data->iSoretEffect == false)
			{
				if (data->iThermophoreticEffect == false)
				{
					for (i = 1; i <= Ni; i++)
					for (j = 1; j <= NC; j++)
					{
						vStar[i][j] = -Dm[i][j] * mix->M[j] / MW_c[i] * (X_c[i + 1][j] - X_c[i][j])*grid.udxe[i];
						vc[i] -= vStar[i][j];
					}

					for (j = 1; j <= NC; j++)
					{
						vStar[Np][j] = -Dm[Np][j] * mix->M[j] / MW_c[Np] * (X_c[Np][j] - X_c[Np - 1][j])*grid.udxw[Np];
						vc[Np] -= vStar[Np][j];
					}
				}
				else
				{
					BzzVector dT_over_dx(Np);
					grid.FirstDerivative('C', T, dT_over_dx);

					for (i = 1; i <= Ni; i++)
					{

						const double Vt = -0.55*mu[i] / rho[i] / T[i] * dT_over_dx[i];
						for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
							vThermophoretic[i][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[i][mix->polimiSoot->bin_indices()[j]];

						for (j = 1; j <= NC; j++)
						{
							vStar[i][j] = -Dm[i][j] * mix->M[j] / MW_c[i] * (X_c[i + 1][j] - X_c[i][j])*grid.udxe[i] + vThermophoretic[i][j];
							vc[i] -= vStar[i][j];
						}
					}

					const double Vt = -0.55*mu[Np] / rho[Np] / T[Np] * dT_over_dx[Np];
					for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
						vThermophoretic[Np][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[Np][mix->polimiSoot->bin_indices()[j]];;

					for (j = 1; j <= NC; j++)
					{
						vStar[Np][j] = -Dm[Np][j] * mix->M[j] / MW_c[Np] * (X_c[Np][j] - X_c[Np - 1][j])*grid.udxw[Np] + vThermophoretic[Np][j];
						vc[Np] -= vStar[Np][j];
					}
				}
			}

			else	// Soret effect: inclusion of thermal diffusion velocities
			{
				if (data->iThermophoreticEffect == false)
				{
					for (i = 1; i <= Ni; i++)
					for (j = 1; j <= NC; j++)
					{
						vStar[i][j] = -Dm[i][j] * mix->M[j] / MW_c[i] * ((X_c[i + 1][j] - X_c[i][j])*grid.udxe[i] + Teta[i][j] / T[i] * (T[i + 1] - T[i])*grid.udxe[i]);
						vc[i] -= vStar[i][j];
					}

					for (j = 1; j <= NC; j++)
					{
						vStar[Np][j] = -Dm[Np][j] * mix->M[j] / MW_c[Np] * ((X_c[Np][j] - X_c[Np - 1][j])*grid.udxw[Np] + Teta[Np][j] / T[Np] * (T[Np] - T[Np - 1])*grid.udxw[Np]);
						vc[Np] -= vStar[Np][j];
					}
				}
				else
				{
					BzzVector dT_over_dx(Np);
					grid.FirstDerivative('C', T, dT_over_dx);

					for (i = 1; i <= Ni; i++)
					{
						const double Vt = -0.55*mu[i] / rho[i] / T[i] * dT_over_dx[i];
						for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
							vThermophoretic[i][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[i][mix->polimiSoot->bin_indices()[j]];;

						for (j = 1; j <= NC; j++)
						{
							vStar[i][j] = -Dm[i][j] * mix->M[j] / MW_c[i] * ((X_c[i + 1][j] - X_c[i][j])*grid.udxe[i] + Teta[i][j] / T[i] * (T[i + 1] - T[i])*grid.udxe[i]) + vThermophoretic[i][j];
							vc[i] -= vStar[i][j];
						}
					}

					const double Vt = -0.55*mu[Np] / rho[Np] / T[Np] * dT_over_dx[Np];
					for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
						vThermophoretic[Np][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[Np][mix->polimiSoot->bin_indices()[j]];

					for (j = 1; j <= NC; j++)
					{
						vStar[Np][j] = -Dm[Np][j] * mix->M[j] / MW_c[Np] * ((X_c[Np][j] - X_c[Np - 1][j])*grid.udxw[Np] + Teta[Np][j] / T[Np] * (T[Np] - T[Np - 1])*grid.udxw[Np]) + vThermophoretic[Np][j];
						vc[Np] -= vStar[Np][j];
					}
				}
			}


			// Correction applied only to Fick contribution
			{
				for (i = 1; i <= Np; i++)
				for (j = 1; j <= NC; j++)
					vStar[i][j] += vc[i] * W_c[i][j];
			}
		}
	}

	// The diffusion contribution are calculated as: omega*V = -D*grad(omega)
	// In principle this is not correct, but it is what is done in OpenFOAM
	else
	{
		int i, j;

		// Fake mass and mole fractions
		{
			BzzVector aux_c(Np);
			for (i = 1; i <= Np; i++)
			{
				aux_c = W.GetRow(i);
				double sum = aux_c.GetSumElements();
				for (j = 1; j <= NC; j++)
					W_c[i][j] = W[i][j] / sum;
			}
			mix->GetMWAndMoleFractionsFromMassFractions(MW_c, X_c, W_c);
		}

		// Correction velocity
		vc = 0.;

		bool correction_velocity_only_fick = true;

		// Correction applied only to Fick contribution
		if (correction_velocity_only_fick == true)
		{
			// No Soret Effect
			if (data->iSoretEffect == false)
			{
				if (data->iThermophoreticEffect == false)
				{
					for (i = 1; i <= Ni; i++)
					for (j = 1; j <= NC; j++)
					{
						vFick[i][j] = -Dm[i][j] * (W_c[i + 1][j] - W_c[i][j])*grid.udxe[i];
						vc[i] -= vFick[i][j];
					}

					for (j = 1; j <= NC; j++)
					{
						vFick[Np][j] = -Dm[Np][j] * (W_c[Np][j] - W_c[Np - 1][j])*grid.udxw[Np];
						vc[Np] -= vFick[Np][j];
					}
				}
				else
				{
					BzzVector dT_over_dx(Np);
					grid.FirstDerivative('B', T, dT_over_dx);

					for (i = 1; i <= Ni; i++)
					{

						const double Vt = -0.55*mu[i] / rho[i] / T[i] * dT_over_dx[i];
						for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
							vThermophoretic[i][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[i][mix->polimiSoot->bin_indices()[j]];

						for (j = 1; j <= NC; j++)
						{
							vFick[i][j] = -Dm[i][j] * (W_c[i + 1][j] - W_c[i][j])*grid.udxe[i];
							vc[i] -= vFick[i][j];
						}
					}

					const double Vt = -0.55*mu[Np] / rho[Np] / T[Np] * dT_over_dx[Np];
					for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
						vThermophoretic[Np][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[Np][mix->polimiSoot->bin_indices()[j]];;

					for (j = 1; j <= NC; j++)
					{
						vFick[Np][j] = -Dm[Np][j] * (W_c[Np][j] - W_c[Np - 1][j])*grid.udxw[Np];
						vc[Np] -= vFick[Np][j];
					}
				}
			}

			else	// Soret effect: inclusion of thermal diffusion velocities
			{
				if (data->iThermophoreticEffect == false)
				{
					for (i = 1; i <= Ni; i++)
					for (j = 1; j <= NC; j++)
					{
						vFick[i][j] = -Dm[i][j] * ((W_c[i + 1][j] - W_c[i][j])*grid.udxe[i] + Teta[i][j] / T[i] * (T[i + 1] - T[i])*grid.udxe[i]);
						vc[i] -= vFick[i][j];
					}

					for (j = 1; j <= NC; j++)
					{
						vFick[Np][j] = -Dm[Np][j] * ((W_c[Np][j] - W_c[Np - 1][j])*grid.udxw[Np] + Teta[Np][j] / T[Np] * (T[Np] - T[Np - 1])*grid.udxw[Np]);
						vc[Np] -= vFick[Np][j];
					}
				}
				else
				{
					BzzVector dT_over_dx(Np);
					grid.FirstDerivative('B', T, dT_over_dx);

					for (i = 1; i <= Ni; i++)
					{
						const double Vt = -0.55*mu[i] / rho[i] / T[i] * dT_over_dx[i];
						for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
							vThermophoretic[i][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[i][mix->polimiSoot->bin_indices()[j]];;

						for (j = 1; j <= NC; j++)
						{
							vFick[i][j] = -Dm[i][j] * ((W_c[i + 1][j] - W_c[i][j])*grid.udxe[i] + Teta[i][j] / T[i] * (T[i + 1] - T[i])*grid.udxe[i]);
							vc[i] -= vFick[i][j];
						}
					}

					const double Vt = -0.55*mu[Np] / rho[Np] / T[Np] * dT_over_dx[Np];
					for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
						vThermophoretic[Np][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[Np][mix->polimiSoot->bin_indices()[j]];

					for (j = 1; j <= NC; j++)
					{
						vFick[Np][j] = -Dm[Np][j] * ((W_c[Np][j] - W_c[Np - 1][j])*grid.udxw[Np] + Teta[Np][j] / T[Np] * (T[Np] - T[Np - 1])*grid.udxw[Np]);
						vc[Np] -= vFick[Np][j];
					}
				}
			}


			// Correction applied only to Fick contribution
			{
				for (i = 1; i <= Np; i++)
				for (j = 1; j <= NC; j++)
					vFick[i][j] += vc[i] * W_c[i][j];

				for (i = 1; i <= Np; i++)
				for (j = 1; j <= NC; j++)
					vStar[i][j] = vFick[i][j] + vThermophoretic[i][j];
			}
		}

		// Correction applied to both the contributions
		else
		{
			// No Soret Effect
			if (data->iSoretEffect == false)
			{
				if (data->iThermophoreticEffect == false)
				{
					for (i = 1; i <= Ni; i++)
					for (j = 1; j <= NC; j++)
					{
						vStar[i][j] = -Dm[i][j] * (W_c[i + 1][j] - W_c[i][j])*grid.udxe[i];
						vc[i] -= vStar[i][j];
					}

					for (j = 1; j <= NC; j++)
					{
						vStar[Np][j] = -Dm[Np][j] * (W_c[Np][j] - W_c[Np - 1][j])*grid.udxw[Np];
						vc[Np] -= vStar[Np][j];
					}
				}
				else
				{
					BzzVector dT_over_dx(Np);
					grid.FirstDerivative('C', T, dT_over_dx);

					for (i = 1; i <= Ni; i++)
					{

						const double Vt = -0.55*mu[i] / rho[i] / T[i] * dT_over_dx[i];
						for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
							vThermophoretic[i][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[i][mix->polimiSoot->bin_indices()[j]];

						for (j = 1; j <= NC; j++)
						{
							vStar[i][j] = -Dm[i][j] * (W_c[i + 1][j] - W_c[i][j])*grid.udxe[i] + vThermophoretic[i][j];
							vc[i] -= vStar[i][j];
						}
					}

					const double Vt = -0.55*mu[Np] / rho[Np] / T[Np] * dT_over_dx[Np];
					for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
						vThermophoretic[Np][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[Np][mix->polimiSoot->bin_indices()[j]];;

					for (j = 1; j <= NC; j++)
					{
						vStar[Np][j] = -Dm[Np][j] * (W_c[Np][j] - W_c[Np - 1][j])*grid.udxw[Np] + vThermophoretic[Np][j];
						vc[Np] -= vStar[Np][j];
					}
				}
			}

			else	// Soret effect: inclusion of thermal diffusion velocities
			{
				if (data->iThermophoreticEffect == false)
				{
					for (i = 1; i <= Ni; i++)
					for (j = 1; j <= NC; j++)
					{
						vStar[i][j] = -Dm[i][j] * ((W_c[i + 1][j] - W_c[i][j])*grid.udxe[i] + Teta[i][j] / T[i] * (T[i + 1] - T[i])*grid.udxe[i]);
						vc[i] -= vStar[i][j];
					}

					for (j = 1; j <= NC; j++)
					{
						vStar[Np][j] = -Dm[Np][j] * ((W_c[Np][j] - W_c[Np - 1][j])*grid.udxw[Np] + Teta[Np][j] / T[Np] * (T[Np] - T[Np - 1])*grid.udxw[Np]);
						vc[Np] -= vStar[Np][j];
					}
				}
				else
				{
					BzzVector dT_over_dx(Np);
					grid.FirstDerivative('C', T, dT_over_dx);

					for (i = 1; i <= Ni; i++)
					{
						const double Vt = -0.55*mu[i] / rho[i] / T[i] * dT_over_dx[i];
						for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
							vThermophoretic[i][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[i][mix->polimiSoot->bin_indices()[j]];;

						for (j = 1; j <= NC; j++)
						{
							vStar[i][j] = -Dm[i][j] * ((W_c[i + 1][j] - W_c[i][j])*grid.udxe[i] + Teta[i][j] / T[i] * (T[i + 1] - T[i])*grid.udxe[i]) + vThermophoretic[i][j];
							vc[i] -= vStar[i][j];
						}
					}

					const double Vt = -0.55*mu[Np] / rho[Np] / T[Np] * dT_over_dx[Np];
					for (j = 1; j <= mix->polimiSoot->bin_indices().Size(); j++)
						vThermophoretic[Np][mix->polimiSoot->bin_indices()[j]] = Vt*W_c[Np][mix->polimiSoot->bin_indices()[j]];

					for (j = 1; j <= NC; j++)
					{
						vStar[Np][j] = -Dm[Np][j] * ((W_c[Np][j] - W_c[Np - 1][j])*grid.udxw[Np] + Teta[Np][j] / T[Np] * (T[Np] - T[Np - 1])*grid.udxw[Np]) + vThermophoretic[Np][j];
						vc[Np] -= vStar[Np][j];
					}
				}
			}


			// Correction applied only to Fick contribution
			{
				for (i = 1; i <= Np; i++)
				for (j = 1; j <= NC; j++)
					vStar[i][j] += vc[i] * W_c[i][j];
			}
		}
	}
}



void OpenSMOKE_Flame1D::compute_vStarPremixed()
{
	int i, j;

	vc_e = 0.;

	if (data->iSoretEffect == false)
	{
		if (data->iThermophoreticEffect == false)
		{
			for(j=1;j<=NC;j++)
			{
				for(i=1;i<=Ni;i++)
					coeff_e[i][j]	= - 0.50 * mix->M[j]*(Dm[i][j]*uPMtot[i]+Dm[i+1][j]*uPMtot[i+1]);
			}

			for(i=1;i<=Ni;i++)
				for(j=1;j<=NC;j++)
				{
					vStar_e[i][j]	 = coeff_e[i][j] * (X[i+1][j]-X[i][j])*grid.udxe[i];
					vc_e[i]			-= vStar_e[i][j];
				}
		}
		else
			ErrorMessage("ThermophoreticEffect not yet implemented for premixed flames!");
	}

	else
	{
		if (data->iThermophoreticEffect == false)
		{
			for(j=1;j<=NC;j++)
			{
				for(i=1;i<=Ni;i++)
					coeff_e[i][j]	= - 0.50 * mix->M[j]*(Dm[i][j]*uPMtot[i]+Dm[i+1][j]*uPMtot[i+1]);
			}

			for(i=1;i<=Ni;i++)
				for(j=1;j<=NC;j++)
				{
					vStar_e[i][j]	 = coeff_e[i][j] * ( (X[i+1][j]-X[i][j])*grid.udxe[i] + Teta[i][j]*(T[i+1]-T[i])/T[i]*grid.udxe[i] );
					vc_e[i]			-= vStar_e[i][j];
				}
		}
		else
			ErrorMessage("ThermophoreticEffect not yet implemented for premixed flames!");
	}

	
	// Correzione delle velocita diffusive
	for(i=1;i<=Ni;i++)
		for(j=1;j<=NC;j++)
		{
			vStar_e[i][j]   += vc_e[i] * 0.50*(W[i][j]+W[i+1][j]);
			vStar_w[i+1][j]  = vStar_e[i][j];
		}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									INITIALIZATION												   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::setup()
{
	ptFlame = this;
	tagConstantT = ResetConstantT;

	setupGasMixture();

	if		(data->gridKind == "EQUISPACED")				grid.Construct(data->Np, data->L, 0.);
	else if (data->gridKind == "STRETCHED")				grid.Construct(data->Np, data->L, data->alfa, 0.);
	else if (data->gridKind == "STRETCHED_STAGNATION")	grid.ConstructStretchedStagnation(data->Np, data->L, data->alfa, 0.);
	else if (data->gridKind == "STRETCHED_POOL_FIRE")	grid.Construct(data->Np, 0., data->L, data->poolfire_grid_alfa_fuel, data->poolfire_grid_alfa_oxidizer, data->poolfire_grid_point_fraction, data->poolfire_grid_distance_fraction);
	else if (data->gridKind == "CENTERED")				grid.Construct(data->Np, data->L, data->xcen, data->wmix, 0.);
	else if (data->gridKind == "USER")					grid.Construct("Grid.inp");
	else ErrorMessage("Error in Grid Geometry, Wrong Type: " + data->gridKind);

	Np = grid.Np;
	Ni = grid.Ni;

	// Setup of QMOM Module
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)
	{
		if (mix->iSootMode == false)
			ErrorMessage("This kinetic scheme was interpreted without the SootMode option and cannot be used together the QMOM model.\nMore details on the OpenSMOKE's User Guide.");

		assign_qmom_module_from_file(*mix, data->qmom_file_name);
	}

	// Setup of TwoEquation Module
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)
	{
		if (mix->iSootMode == false)
			ErrorMessage("This kinetic scheme was interpreted without the SootMode option and cannot be used together the 2E model.\nMore details on the OpenSMOKE's User Guide.");

		soot2EModel.setupFromFile(data->twoEquation_file_name);
		soot2EModel.assign_mixture(*mix);
	}

	// Memory allocation
	cout << "Allocating memory... ";
	allocate_all();
	cout << "Succesfully DONE!" << endl;

	// Radiation
	if ( data->iGasRadiation == true || data->iRadiativeSootModel != RADIATIVE_SOOT_MODEL_NONE )
		prepare_radiation();
	
	// Initial and Boundary conditions for conventional variables
	setupBoundaryConditions();

	// Initial and Boundary conditions for QMOM
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)
		initial_conditions_qmom_module();

	// Initial and Boundary conditions for TwoEquation Module
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)
		initial_conditions_soot_module();

	// Check List Of sensitivity species
//	int j;
//	if (data->names_sensitivity[1] == "ALL")
//	{
//		data->listSensitivitySpecies.resize(NC+1);
//		for(j=1;j<=NC;j++)
//			data->listSensitivitySpecies[j] = mix->names[j];
//	}
//	for(j=1;j<=data->nSensitivity;j++)
//		mix->recognize_species(data->listSensitivitySpecies[j]);
}

void OpenSMOKE_Flame1D::setupGasMixture()
{
	NC		= mix->NumberOfSpecies();
	NR		= mix->NumberOfReactions();
	nBlock	= 4 + NC;

//	mix->pah_manager.recognizeSpecies(NC, mix->names, mix->M);

	if (data->bin_index_zero > 0)
		mix->SetPolimiSoot(data->bin_index_zero, data->bin_density_A, data->bin_index_final, data->bin_density_B);
	else
		mix->SetPolimiSoot();
}

void OpenSMOKE_Flame1D::setupBoundaryConditions()
{
	if ( data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED	)
	{
		int i;
		double csiST;

		massFractionsAndPM(data->XC, WC, PMtot[1], mix->M);
		massFractionsAndPM(data->XO, WO, PMtot[Np], mix->M);

		setTemperatureAndMassFractionProfiles(FLAME1D_IGNITION_START);
		BzzVectorInt dummy;
		properties(0, -1, dummy, 0, 0);

		if (data->geometry == "AXIS")		nGeometry = 3.; 
		else if (data->geometry == "PLANAR") nGeometry = 2.; 
		else ErrorMessage("Error in Grid Geometry, Wrong Type " + data->geometry);

		UC =  rho[1]*data->VC / (nGeometry-1.);
		UO = -rho[Np]*data->VO / (nGeometry-1.);

		GC = -rho[1]*data->radialGradientC;
		GO = -rho[Np]*data->radialGradientO;

		for(i=1;i<=NC;i++)
		{
			BCW_C[i] =  rho[1] *data->VC*WC[i];
			BCW_O[i] = -rho[Np]*data->VO*WO[i];
		}

		csiST = 1./(1.+rho[Np]*data->VO*data->VO/(rho[1]*data->VC*data->VC));
		xST = grid.L*csiST;
		K = 2.*data->VO/grid.L*(1.+data->VC/data->VO*sqrt(rho[1]/rho[Np]));

		printSummaryForOpposed();
	}

  
	else if( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED)
	{
		int i;
		double sum;

		sum = 0.;
		for(i=1;i<=NC;i++)
			if(i!=data->jINERT) sum+=data->XO[i];
		data->XO[data->jINERT] = 1.-sum;

		massFractionsAndPM(data->XC, WC, PMtot[1], mix->M);	// INLET
		massFractionsAndPM(data->XO, WO, PMtot[Np], mix->M);	// OUTLET (estimated)
		setTemperatureAndMassFractionProfiles(FLAME1D_IGNITION_NONE);

		BzzVectorInt dummy;
		properties(0, -1, dummy, 0, 0);

		for(i=1;i<=NC;i++)
		{
			BCW_C[i] =  WC[i];		// INLET
			BCW_O[i] =  WO[i];		// OUTLET (estimated)
		}

		printSummaryForPremixed();
	}

	else if (	data->kind_of_flame == FLAME1D_PHYSICS_TWIN	)
	{
		int i;

		massFractionsAndPM(data->XC, WC, PMtot[1], mix->M);

		setTemperatureAndMassFractionProfiles(FLAME1D_IGNITION_START);
		BzzVectorInt dummy;
		properties(0, -1, dummy, 0, 0);

		if (data->geometry == "AXIS")			nGeometry = 3.; 
		else if (data->geometry == "PLANAR")	nGeometry = 2.; 
		else ErrorMessage("Error in Grid Geometry, Wrong Type " + data->geometry);

		UC =  rho[1]*data->VC / (nGeometry-1.);
		GC = -rho[1]*data->radialGradientC;

		for(i=1;i<=NC;i++)
			BCW_C[i] =  rho[1] *data->VC*WC[i];

		xST = grid.L;
		K = 4.*data->VC/(2.*grid.L);

		printSummaryForTwin();
	}
}

void OpenSMOKE_Flame1D::initializeVariables(const flame1d_ignition string_kind)
{
	int i;

	if ( data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED	)
	{
		double Gmin = -10.;
		double Hstar = -100.;

		for(i=1;i<=Np;i++)
		{
			if (grid.x[i]<=xST) 
			{
				U[i] = UC - UC/(xST-grid.xA)*(grid.x[i]-grid.xA);
				G[i] = GC + (Gmin-GC)/(xST-grid.xA)*(grid.x[i]-grid.xA);
			}
			if (grid.x[i] >xST) 
			{
				U[i] = UO/(grid.xB-xST)*(grid.x[i]-xST);
				G[i] = GO + (Gmin-GO)*(1. - 1./(grid.xB-xST)*(grid.x[i]-xST));
			}
		}

		H = Hstar;

		// In case of NO_SPARK ignition
		setTemperatureAndMassFractionProfiles(string_kind);
	}

	else if( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED )
	{

		H = data->MassFlowRate;					// portata massiva [kg/s]

		if (data->iAreaProfile==0)
			G = data->CrossSectionalArea;		// area di passaggio [m2]

		if (data->iAreaProfile==1)
			for(i=1;i<=Np;i++)
				G[i] = data->ud_cross_section_profile.GiveMeValue(0., grid.x[i]);

		// Fitting formula from  Bittner PhD Thesis (MIT)
		if (data->iAreaProfile==2)
			for(i=1;i<=Np;i++)
				G[i] = data->CrossSectionalArea * (1.+0.114*pow((grid.x[i]*1.e2), 1.44));		// area di passaggio [m2]

		// Calcolare proprieta
		BzzVectorInt dummy;
		properties(data->iHOT, -1, dummy, 0, 0);

		for(i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
	}

	else if ( data->kind_of_flame == FLAME1D_PHYSICS_TWIN )
	{
		double Gmin = -10.;
		double Hstar = -100.;

		for(i=1;i<=Np;i++)
		{
			U[i] = UC - UC/(grid.L-grid.xA)*(grid.x[i]-grid.xA);
			G[i] = GC + (Gmin-GC)/(grid.L-grid.xA)*(grid.x[i]-grid.xA);
		}

		H = Hstar;

		// In case of NO_SPARK ignition
		setTemperatureAndMassFractionProfiles(string_kind);
	}
}

void OpenSMOKE_Flame1D::setTemperatureAndMassFractionProfiles(const flame1d_ignition string_kind)
{
	int i, j;
	double sum;

	// ---------------------------------------------------------------------------
	// OPPOSED
	// ---------------------------------------------------------------------------
	if ( data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED	) 
	{
		// ---------------------------------------------------------------------------
		// Profili iniziali per il calcolo delle proprieta (START)
		// Profili iniziali nel caso di ignizione con scintilla (SPARK)
		// ---------------------------------------------------------------------------
		if (string_kind == FLAME1D_IGNITION_SPARK) 
			ErrorMessage("Spark ignition is not yet available");

		if (string_kind == FLAME1D_IGNITION_START) 
		{
			for(i=1;i<=Np;i++)
				T[i] = data->TC + (data->TO-data->TC)*(grid.x[i]-grid.xA)/grid.L;

			for(i=1;i<=Np;i++)
				for(j=1;j<=NC;j++)
					W[i][j] = WC[j] + (WO[j]-WC[j])*(grid.x[i]-grid.xA)/grid.L;

			for(i=1;i<=Np;i++)
			{
				sum = 0.;
				for(j=1;j<=NC;j++)
				if(j!=data->jINERT) sum+=W[i][j];
				W[i][data->jINERT] = 1.-sum;
			}
		}

		// ---------------------------------------------------------------------------
		// Profili iniziali nel caso di ignizione (USER DEFINED)
		// ---------------------------------------------------------------------------
		else if (string_kind == FLAME1D_IGNITION_NO_SPARK)
		{
			BzzVector sumX(Np); sumX = 0.;

			// Grid zones
			gridZones grid_for_interpolation;
			grid_for_interpolation.setup(100, grid.L, data->xcen, data->wmix);

			// 1. Temperatura
			// ---------------------------------------------------------------------------		
			// a. Profilo assegnato dall'utente
			if (data->iTemperatureProfile==1 || data->iTemperatureProfile==2)
				for(i=1;i<=Np;i++)
					T[i] = data->ud_temperature_profile.GiveMeValue(0., grid.x[i]);;

			// b. Profilo stimato
			if (data->iTemperatureProfile==0 && data->iPoolFire == POOL_FIRE_NONE)	// Profilo Peak
			{
				LinearInterpolation TinitialProfile;
				functionStep(TinitialProfile, grid_for_interpolation, data->TC, data->Tpeak, data->TO);

				for(i=1;i<=Np;i++)
					T[i] = TinitialProfile(grid.x[i]);
			}
			
			else if (data->iTemperatureProfile==0 && data->iPoolFire != POOL_FIRE_NONE)	// Profilo Peak
			{
				double xPeak = grid.x[Np]*0.2;
				double xSupport = grid.x[Np]*0.5;
				for(i=1;i<=Np;i++)
				{
					if (grid.x[i]<=xPeak)			T[i] = data->TC + (data->Tpeak-data->TC)/(xPeak)*(grid.x[i]);
					else if (grid.x[i]>=xSupport)	T[i] = data->TO;
					else 							T[i] = data->Tpeak + (data->TO-data->Tpeak)/(xSupport-xPeak)*(grid.x[i]-xPeak);
				}
			}

			// 2. Frazioni Molari
			// ---------------------------------------------------------------------------		
			// a. Intermediate species
			for(j=1;j<=data->iPeaks.Size();j++)
			{
				LinearInterpolation XinitialProfile;
				functionGaussian(XinitialProfile, grid_for_interpolation, data->xPeaks[j], 10.);

				for(i=1;i<=Np;i++)
				{	
					X[i][data->iPeaks[j]] = XinitialProfile(grid.x[i]);
					sumX[i] += X[i][data->iPeaks[j]];
				}
			}
				
			// b. Reactants and main products
			if (data->iPoolFire == POOL_FIRE_NONE)
			{
				for(i=1;i<=Np;i++)
				{
					for(j=1;j<=data->iXC.Size();j++)
					{
						if (data->iXC[j]!=data->jINERT)
						{
							if(grid.x[i]< grid_for_interpolation.xLeft) X[i][data->iXC[j]] = data->XC[data->iXC[j]] - sumX[i];
							else X[i][data->iXC[j]] = 0.;
						}
					}
			
					for(j=1;j<=data->iXO.Size();j++)
					{
						if (data->iXO[j]!=data->jINERT)
						{
							if(grid.x[i]> grid_for_interpolation.xRight) X[i][data->iXO[j]] = data->XO[data->iXO[j]] - sumX[i];
							else X[i][data->iXO[j]] = 0.;
						}
					}		
				}

				// c. Inert
				for(i=1;i<=Np;i++)
				{
					sum = 0.;
					for(j=1;j<=NC;j++)
						if(j!=data->jINERT) sum+=X[i][j];
					X[i][data->jINERT] = 1.-sum;
				}


			}

			else
			{
				double xPeak = grid.x[Np]*0.2;
				double xSupport = grid.x[Np]*0.5;

				for(i=1;i<=Np;i++)
				{
					for(j=1;j<=data->iXC.Size();j++)
					{
						if(grid.x[i]<=xPeak) X[i][data->iXC[j]] = data->XC[data->iXC[j]] - data->XC[data->iXC[j]]/xPeak*grid.x[i];
						else				X[i][data->iXC[j]] = 0.;
					}
			
					for(j=1;j<=data->iXO.Size();j++)
					{
						if (grid.x[i]>=xSupport)			X[i][data->iXO[j]] = data->XO[data->iXO[j]];
						else if (grid.x[i]<xPeak)			X[i][data->iXO[j]] = 0.;
						else								X[i][data->iXO[j]] = data->XO[data->iXO[j]]/(xSupport-xPeak)*(grid.x[i]-xPeak);
					}

					double sum = 1.-X.GetRow(i).GetSumElements();
					if (sum < 0.)
						ErrorMessage("Wrong initial conditions on mass fractions (internal error)");
					X[i][data->jINERT] += sum;
				}
			}


			// 3. Frazioni Massive
			// --------------------------------------------------------------------------
			massFractionsAndPMtot();
		}
	}


	// ---------------------------------------------------------------------------
	// PREMIXED
	// ---------------------------------------------------------------------------

	else if ( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED ) 
	{
		BzzVector sumX(Np); sumX = 0.;

		// Grid zones
		gridZones grid_for_interpolation;
		grid_for_interpolation.setup(100, grid.L, data->xcen, data->wmix);

		// 1. Temperatura			
		if (data->iTemperatureProfile==1 || data->iTemperatureProfile==2 || data->iAssignedFixedTemperatureProfileProvisional == true)
			for(i=1;i<=Np;i++)
				T[i] = data->ud_temperature_profile.GiveMeValue(0., grid.x[i]);

		if (data->iTemperatureProfile==0 && data->iAssignedFixedTemperatureProfileProvisional == false)
		{
			LinearInterpolation TinitialProfile;
			functionRampa(TinitialProfile, grid_for_interpolation, data->TC, data->TO);

			for(i=1;i<=Np;i++)
				T[i] = TinitialProfile(grid.x[i]);
		}

		// 2. Frazioni Molari
		// ---------------------------------------------------------------------------
		// a. Intermediate species
		for(j=1;j<=data->iPeaks.Size();j++)
		{
			LinearInterpolation XinitialProfile;
			functionGaussian(XinitialProfile, grid_for_interpolation, data->xPeaks[j], 10.);

			for(i=1;i<=Np;i++)
			{	
				X[i][data->iPeaks[j]] = XinitialProfile(grid.x[i]);
				sumX[i] += X[i][data->iPeaks[j]];
			}
		}
		// b. Reactants and main products
		for(j=1;j<=NC;j++)
		{
			LinearInterpolation XinitialProfile;
			functionRampa(XinitialProfile, grid_for_interpolation, data->XC[j], data->XO[j]);

			for(i=1;i<=Np;i++)
				X[i][j] = XinitialProfile(grid.x[i]) * (1.-sumX[i]);
		}
	
		// c. Intermediate species
		for(j=1;j<=data->iPeaks.Size();j++)
		{
			LinearInterpolation XinitialProfile;
			functionGaussian(XinitialProfile, grid_for_interpolation, data->xPeaks[j], 10.);

			for(i=1;i<=Np;i++)
				X[i][data->iPeaks[j]] = XinitialProfile(grid.x[i]);
		}
	
		// c. Inert
		for(i=1;i<=Ni;i++)
		{
			X[i][data->jINERT] = 1.;
			for(j=1;j<data->jINERT;j++)
				X[i][data->jINERT] = 1.-X[i][j];
			for(j=data->jINERT+1;j<=NC;j++)
				X[i][data->jINERT] = 1.-X[i][j];
				
			// Round Off Errors
			if (X[i][data->jINERT]<0.e0) 
			{
				cout << "WARNING: Correction on initial conditions due to round-off errors:!! " << endl;
				cout << "         Point: " << i << "\t" << "Error: " << X[i][data->jINERT] << endl;
			
				cout.setf(ios::scientific);
				for(j=1;j<=NC;j++)
					if(j!=data->jINERT) X[i][j]/=(1.-X[i][data->jINERT]);
				X[i][data->jINERT] = 0.;
			}
		}

		for(j=1;j<=NC;j++)
			X[Np][j] = X[Ni][j];
			
		// 2. Frazioni Massive
		// --------------------------------------------------------------------------
		massFractionsAndPMtot();
		cout << "Initialization OK!" << endl;
	}

	// ---------------------------------------------------------------------------
	// TWIN
	// ---------------------------------------------------------------------------
	else if ( data->kind_of_flame == FLAME1D_PHYSICS_TWIN ) 
	{
		// ---------------------------------------------------------------------------
		// Profili iniziali per il calcolo delle proprieta (START)
		// Profili iniziali nel caso di ignizione con scintilla (SPARK)
		// ---------------------------------------------------------------------------
		if (string_kind == FLAME1D_IGNITION_SPARK) 
			ErrorMessage("Spark ignition is not yet available");

		if (string_kind == FLAME1D_IGNITION_START) 
		{

				for(i=1;i<=Np;i++)
						T[i] = data->TC;

				for(i=1;i<=Np;i++)
					for(j=1;j<=NC;j++)
						W[i][j] = WC[j];

				for(i=1;i<=Np;i++)
				{
					sum = 0.;
					for(j=1;j<=NC;j++)
					if(j!=data->jINERT) sum+=W[i][j];
					W[i][data->jINERT] = 1.-sum;
				}

		}

		// ---------------------------------------------------------------------------
		// Profili iniziali nel caso di ignizione (USER DEFINED)
		// ---------------------------------------------------------------------------
		else if (string_kind == FLAME1D_IGNITION_NO_SPARK)
		{
			{
				BzzVector sumX(Np); sumX = 0.;

				// Grid zones
				gridZones grid_for_interpolation;
				grid_for_interpolation.setup(100, grid.L, data->wmix);

				// 1. Temperatura
				// ---------------------------------------------------------------------------		
				// a. Profilo assegnato dall'utente
				if (data->iTemperatureProfile==1 || data->iTemperatureProfile==2)
					for(i=1;i<=Np;i++)
						T[i] = data->ud_temperature_profile.GiveMeValue(0., grid.x[i]);

				// b. Profilo stimato
				if (data->iTemperatureProfile==0)	// Profilo Peak
				{
					LinearInterpolation TinitialProfile;
					functionStepTwin(TinitialProfile, grid_for_interpolation, data->TC, data->Tpeak);

					for(i=1;i<=Np;i++)
						T[i] = TinitialProfile(grid.x[i]);
				}

				// 2. Frazioni Molari
				// ---------------------------------------------------------------------------		
				// a. Intermediate species
				for(j=1;j<=data->iPeaks.Size();j++)
				{
					// TODO
					LinearInterpolation XinitialProfile;
					functionGaussian(XinitialProfile, grid_for_interpolation, data->xPeaks[j], 10.);

					for(i=1;i<=Np;i++)
					{	
						X[i][data->iPeaks[j]] = XinitialProfile(grid.x[i]);
						sumX[i] += X[i][data->iPeaks[j]];
					}
				}
				
				// b. Reactants and main products
				for(i=1;i<=Np;i++)
				{
					for(j=1;j<=data->iXC.Size();j++)
					{
						if (data->iXC[j]!=data->jINERT)
						{
							if(grid.x[i]< grid_for_interpolation.xLeft) X[i][data->iXC[j]] = data->XC[data->iXC[j]] - sumX[i];
							else X[i][data->iXC[j]] = 0.;
						}
					}		
				}		

				// c. Inert
				for(i=1;i<=Np;i++)
				{
					sum = 0.;
					for(j=1;j<=NC;j++)
						if(j!=data->jINERT) sum+=X[i][j];
					X[i][data->jINERT] = 1.-sum;
				}

				// 3. Frazioni Massive
				// --------------------------------------------------------------------------
				massFractionsAndPMtot();
			}			
		}
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								PRINT functions											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::printOnFile(const std::string fileNameOutput)
{
	int i, j;

	openOutputFileAndControl(fOutput, fileNameOutput);
	fOutput.setf(ios::scientific);

	GnuPlotInterface(fOutput);

	if ( data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED || data->kind_of_flame == FLAME1D_PHYSICS_TWIN )
	{
		ofstream fSoot;
		if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )
		{
			std::string fileName;
			fileName  = nameFolderAdditionalData + "/Soot_2E.out";
			openOutputFileAndControl(fSoot, fileName);
			fSoot.setf(ios::scientific);

            fSoot	<< setw(20) << left << "x[cm](1)";
            soot2EModel.GnuPlotInterface(fSoot, 2);
			fSoot << endl << endl;
		}

		for(i=1;i<=Np;i++)
		{
			fOutput << setw(20) << left << 1.e2*grid.x[i]; 
			fOutput << setw(20) << left << U[i];
			fOutput << setw(20) << left << 100.*(nGeometry-1.)*U[i]*urho[i];
			fOutput << setw(20) << left << G[i];
			fOutput << setw(20) << left << H[i];
			fOutput << setw(20) << left << T[i];
			
			fOutput << setw(20) << left << rho[i];
			fOutput << setw(20) << left << mu[i];
			fOutput << setw(20) << left << lambda[i];
			fOutput << setw(20) << left << Cp[i];
			fOutput << setw(20) << left << QReaction[i];
			fOutput << setw(20) << left << Qrad[i];

			// Elemental
			ElementalAnalysis();
			fOutput << setw(20) << left << Z[i];
			for(j=1;j<=mix->NumberOfElements();j++)
				fOutput << setw(20) << left << x_elemental[i][j];
			for(j=1;j<=mix->NumberOfElements();j++)
				fOutput << setw(20) << left << omega_elemental[i][j];

			// Species
			for(j=1;j<=data->iOut.Size();j++)
				fOutput << setw(20) << left << X[i][data->iOut[j]];
			for(j=1;j<=data->iOut.Size();j++)
				fOutput << setw(20) << left << W[i][data->iOut[j]];

			fOutput << endl;
			
			if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )
			{
				X.GetRow(i, &xVector);
				soot2EModel.update(T[i], data->P_atm, rho[i], 0., xVector, phiN[i], phiM[i]);
				soot2EModel.formation_rates();

				fSoot	<< setw(20) << left << 1.e2*grid.x[i];
				soot2EModel.write_on_file(fSoot, phiN[i], phiM[i]);
			}	

			if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM )
			{
				BzzVector momentsVector = moments.GetRow(i);
				qmom.updateMoments(momentsVector);
				
				double muDummy	=1.0;
				double trDummy	=1.0;
				double fDummy	=0.0;
				qmom.updateData(	T[i], data->P_Pascal, rho[i], muDummy, trDummy, 
									W[i][qmom.jO2], W[i][qmom.jC2H2], W[i][qmom.jOH], fDummy );
				qmom_sources = qmom.calculateSources();
				moments_source.SetRow(i, qmom_sources);
				qmom.soot.formation_rates();

				qmom.soot.write_on_file(fSootQMOM, 1.e2*grid.x[i], momentsVector);
			}
		}
		
		if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM )
			fSoot.close();
	}

	else if ( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED )
	{
		ofstream fSoot;
		if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )
		{
			std::string fileName;
			fileName = nameFolderAdditionalData + "/Soot_2E.out";
			openOutputFileAndControl(fSoot, fileName);
			fSoot.setf(ios::scientific);

            fSoot	<< setw(20) << left << "x[cm](1)";
            soot2EModel.GnuPlotInterface(fSoot, 2);
			fSoot << endl << endl;
		}
		
		for(i=1;i<=Np;i++)
		{
			fOutput << setw(20) << left << 1.e2*grid.x[i];
			fOutput << setw(20) << left << U[i]*100.;
			fOutput << setw(20) << left << H[i];
			fOutput << setw(20) << left << H[i]/G[i];
			fOutput << setw(20) << left << G[i];
			fOutput << setw(20) << left << T[i];

			fOutput << setw(20) << left << rho[i];
			fOutput << setw(20) << left << mu[i];
			fOutput << setw(20) << left << lambda[i];
			fOutput << setw(20) << left << Cp[i];
			fOutput << setw(20) << left << QReaction[i];
			fOutput << setw(20) << left << Qrad[i];

			// Elemental
			ElementalAnalysis();
			fOutput << setw(20) << left << Z[i];
			for(j=1;j<=mix->NumberOfElements();j++)
				fOutput << setw(20) << left << x_elemental[i][j];
			for(j=1;j<=mix->NumberOfElements();j++)
				fOutput << setw(20) << left << omega_elemental[i][j];

			// Species
			for(j=1;j<=data->iOut.Size();j++)
				fOutput << setw(20) << left << X[i][data->iOut[j]];
			for(j=1;j<=data->iOut.Size();j++)
				fOutput << setw(20) << left << W[i][data->iOut[j]];

			fOutput << endl;

			
			if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM ) 
			{
				BzzVector momentsVector = moments.GetRow(i);
				qmom.updateMoments(momentsVector);
				
				double muDummy	=1.0;
				double trDummy	=1.0;
				double fDummy	=0.0;
				qmom.updateData(	T[i], data->P_Pascal, rho[i], muDummy, trDummy, 
									W[i][qmom.jO2], W[i][qmom.jC2H2], W[i][qmom.jOH], fDummy );
				qmom_sources = qmom.calculateSources();
				moments_source.SetRow(i, qmom_sources);
				qmom.soot.formation_rates();

				qmom.soot.write_on_file(fSootQMOM, 1.e2*grid.x[i], momentsVector);
			}
			
			if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )
			{
				X.GetRow(i, &xVector);
				soot2EModel.update(T[i], data->P_atm, rho[i], 0., xVector, phiN[i], phiM[i]);
				soot2EModel.formation_rates();
				
				fSoot	<< setw(20) << left << 1.e2*grid.x[i];
				soot2EModel.write_on_file(fSoot, phiN[i], phiM[i]);
			}	
		}
		
		if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )
			fSoot.close();
	}
	fOutput.close();
	if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM ) 
		fSootQMOM.close();

	if (mix->polimiSoot->IsSoot() == true)
	{
		std::string fileName;

		ofstream fSoot;
		ofstream fSootDistribution;
	//	ofstream fPAH;	
		
		fileName = nameFolderAdditionalData + "/SootDetailed.out";
		openOutputFileAndControl(fSoot, fileName);
		fSoot.setf(ios::scientific);
/*
		fileName = nameFolderAdditionalData + "/PAHDetailed.out";
		openOutputFileAndControl(fPAH, fileName);
		fPAH.setf(ios::scientific);
*/
		fileName = nameFolderAdditionalData + "/SootDistribution.out";
		openOutputFileAndControl(fSootDistribution, fileName);
		fSootDistribution.setf(ios::scientific);

		GnuPlotSootInterface(fSoot);
		GnuPlotSootDistributionInterface(fSootDistribution);
//		GnuPlotPAHInterface(fPAH);

		for(i=1;i<=Np;i++)
		{
			X.GetRow(i, &cVector);
			double cTotForSoot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			cVector *= cTotForSoot;

			double MWgas=0.;
			for(int k=1;k<=NC;k++)
				MWgas += X[i][k]*mix->M[k];
			double rhoGas = cTotForSoot*MWgas;

			// PAH Distribution data
			W.GetRow(i, &wVector);
			X.GetRow(i, &xVector);
//			mix->pah_manager.pah_calculateAll(	 wVector, xVector, cVector, rhoGas);

			// Soot Distribution data
			mix->polimiSoot->Analysis(*mix, data->P_Pascal, T[i], wVector);
			mix->polimiSoot->Distribution();
			mix->polimiSoot->ProcessDistribution();

			// Print On File - PAH data
//			fPAH << setw(20) << left << 1.e2*grid.x[i];
//			fPAH << setw(20) << left << T[i];
//            mix->pah_manager.print_main_data_on_file(fPAH);
 //           fPAH << endl;

			// Print On File - Soot data
            fSoot << setw(20) << left << 1.e2*grid.x[i];
			fSoot << setw(20) << left << T[i];

			fSoot << setw(20) << left << mix->polimiSoot->fv_large();
			fSoot << setw(20) << left << mix->polimiSoot->x_large();
			fSoot << setw(20) << left << mix->polimiSoot->omega_large();
			fSoot << setw(20) << left << mix->polimiSoot->rho_large();
			fSoot << setw(20) << left << mix->polimiSoot->N_large();
			fSoot << setw(20) << left << mix->polimiSoot->h_over_c_large();
			fSoot << setw(20) << left << mix->polimiSoot->dmean_N_large();
			fSoot << setw(20) << left << mix->polimiSoot->d32_N_large();
			fSoot << setw(20) << left << mix->polimiSoot->dstd_N();

			fSoot << setw(20) << left << mix->polimiSoot->fv_small();
			fSoot << setw(20) << left << mix->polimiSoot->x_small();
			fSoot << setw(20) << left << mix->polimiSoot->omega_small();
			fSoot << setw(20) << left << mix->polimiSoot->rho_small();
			fSoot << setw(20) << left << mix->polimiSoot->N_small();

            fSoot << endl;

			// Distribution
			{
				for (int k=1;k<=mix->polimiSoot->baskets_d().Size();k++)
				{
					fSootDistribution << setw(20) << left << 1.e2*grid.x[i];
					fSootDistribution << setw(20) << left << k;
					
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_d()[k];
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_m()[k];
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_V()[k];

					fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_d()[k]);
					fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_m()[k]);
					fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_V()[k]);

					fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_dlog10d()[k]);
					fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_dlog10m()[k]);
					fSootDistribution << setw(20) << left << log10(mix->polimiSoot->baskets_dlog10V()[k]);

					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_fv()[k];
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_x()[k];
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_omega()[k];
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_rho()[k];
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_N()[k];
					
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_fv()[k]/(mix->polimiSoot->fv_large()+mix->polimiSoot->fv_small());
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_x()[k]/(mix->polimiSoot->x_large()+mix->polimiSoot->x_small());
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_omega()[k]/(mix->polimiSoot->omega_large()+mix->polimiSoot->omega_small());
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_rho()[k]/(mix->polimiSoot->rho_large()+mix->polimiSoot->rho_small());
					fSootDistribution << setw(20) << left << mix->polimiSoot->baskets_N()[k]/(mix->polimiSoot->N_large()+mix->polimiSoot->N_small());

					fSootDistribution << endl;
				}
			}
			fSootDistribution << endl;
		}

		fSoot.close();
//		fPAH.close();
		fSootDistribution.close();
	}

	if (data->iLewisNumbers == true)
	{
		std::string fileName;

		ofstream fLewis;
		fileName = nameFolderAdditionalData + "/Lewis.out";
		openOutputFileAndControl(fLewis, fileName);
		fLewis.setf(ios::scientific);

		GnuPlotLewisInterface(fLewis);

		for(i=1;i<=Np;i++)
		{
			int k;
			BzzVector Le(NC);
			BzzVector DVector(NC);

			// a. Calcolo delle frazioni massive e molari
			W.GetRow(i,&wVector);
			X.GetRow(i,&xVector);
			Dm.GetRow(i,&DVector);

			// b. Calcolo dei numeri di Lewis
			double alfa = lambda[i] / (rho[i]*Cp[i]);
			for(k=1;k<=NC;k++)
				Le[k] = alfa / Dm[i][k];

			double Le_MeanMass	= Dot(Le, wVector);
			double Le_MeanMole	= Dot(Le, xVector);
			double Le_TildeMass = alfa/Dot(DVector, wVector);
			double Le_TildeMole = alfa/Dot(DVector, xVector);


			// Print On File - Main data
			fLewis << setw(20) << left << 1.e2*grid.x[i];	// 1
			fLewis << setw(20) << left << T[i];	// 2
			fLewis << setw(20) << left << alfa;	// 3
			fLewis << setw(20) << left << Le_MeanMass;	// 4
			fLewis << setw(20) << left << Le_MeanMole;	// 5
			fLewis << setw(20) << left << Le_TildeMass;	// 6
			fLewis << setw(20) << left << Le_TildeMole;	// 7
			for(k=1;k<=NC;k++)
				fLewis << setw(20) << left << Le[k];	// 8...
			fLewis << endl;
		}

		fLewis.close();
	}

	if (data->iSoretEffect == true)
	{
		std::string fileName;

		ofstream fSoret;
		fileName = nameFolderAdditionalData + "/Soret.out";
		openOutputFileAndControl(fSoret, fileName);
		fSoret.setf(ios::scientific);

		GnuPlotSoretInterface(fSoret);
		
		for(i=1;i<=Np;i++)
		{
			fSoret << setw(20) << left << 1.e2*grid.x[i];	// 1
			fSoret << setw(20) << left << T[i];				// 2

			for(j=1;j<=NC;j++)	
				if (mix->M[j] <= 5.)
					fSoret << setw(20) << left << Teta[i][j];
			fSoret << endl;
		}
		fSoret.close();
	}

	if (data->iAssignedFormationRates == true)
	{
		std::string fileName;

		ofstream fFormationRates;
		fileName = nameFolderAdditionalData + "/FormationRates.out";
		openOutputFileAndControl(fFormationRates, fileName);
		fFormationRates.setf(ios::scientific);

		GnuPlotFormationRatesInterface(fFormationRates);

		for(i=1;i<=Np;i++)
		{
			// Print On File - Main data
			fFormationRates << setw(20) << left << 1.e2*grid.x[i];	// 1
			fFormationRates << setw(20) << left << T[i];	// 2
	
			for(int k=1;k<=data->index_formation_rates.Size();k++)
				fFormationRates << setw(20) << left << R[i][data->index_formation_rates[k]]/mix->M[data->index_formation_rates[k]];	// 3
			fFormationRates << endl;
		}
		fFormationRates.close();
	}

	if (data->iAssignedExperiment == true)
	{
		std::string fileName;
		double threshold = 1.e-6;

		ofstream fExperiment;
		fileName = nameFolderAdditionalData + "/Experiment.out";
		openOutputFileAndControl(fExperiment, fileName);
		fExperiment.setf(ios::scientific);

		fExperiment << "REACTOR ";
		if (data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED)	fExperiment << "PREMIXED"	<< endl;
		if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED)		fExperiment << "CFDF"		<< endl;
		if (data->kind_of_flame == FLAME1D_PHYSICS_TWIN)		fExperiment << "CFDF-TWIN"	<< endl;
		
		fExperiment << "SPACE   cm" << endl;
		fExperiment << "DATA    "   << 1+data->names_Experiment.size() << endl;
		fExperiment << endl;
		
		fExperiment << "T TEMPERATURE 1.0 CROSS" << endl;
		for(i=3;i<=Ni-1;i++)
			fExperiment << 1.e2*grid.x[i] << "\t" << T[i] << endl;
		fExperiment << "//" << endl << endl;
		
		for(int k=0;k<int(data->names_Experiment.size());k++)
		{
			fExperiment << data->names_Experiment[k] << " MOLE_FRACTION 1.0 CROSS" << endl;
			for(i=3;i<=Ni-1;i++)
				if ( X[i][data->index_Experiment[k+1]] > threshold)
					fExperiment << 1.e2*grid.x[i] << "\t" << X[i][data->index_Experiment[k+1]] << endl;
			fExperiment << "//" << endl << endl;
		}

		fExperiment.close();
	}

	if (data->iAssignedReactionRates == true)
	{
		std::string fileName;

		ofstream fReactionRates;
		fileName = nameFolderAdditionalData + "/ReactionRates.out";
		openOutputFileAndControl(fReactionRates, fileName);
		fReactionRates.setf(ios::scientific);

		GnuPlotReactionRatesInterface(fReactionRates);

		{
			data->iFocusReactionRates = true;
			BzzVectorInt dummy;
			properties(data->iHOT, -1, dummy, 0, 0);
			data->iFocusReactionRates = false;
		}

		for(i=1;i<=Np;i++)
		{
			// Print On File - Main data
			fReactionRates << setw(20) << left<< 1.e2*grid.x[i];	// 1
			fReactionRates << setw(20) << left<< T[i];	// 2
	
			for(int k=1;k<=data->index_reaction_rates.Size();k++)
				fReactionRates << setw(20) << left<< RR[i][data->index_reaction_rates[k]];	// 3
			fReactionRates << endl;
		}
		fReactionRates.close();
	}

	if (data->iVerboseMixtureProperties == true)
	{
		std::string fileName;

		ofstream fProperties;
		fileName = nameFolderAdditionalData + "/MixtureProperties.out";
		openOutputFileAndControl(fProperties, fileName);
		fProperties.setf(ios::scientific);

		GnuPlotMixturePropertiesInterface(fProperties);

		for(i=1;i<=Np;i++)
		{
			// Print On File - Main data
			fProperties << setw(20) << left << 1.e2*grid.x[i];	// 1
			fProperties << setw(20) << left << T[i];			// 2
			fProperties << setw(20) << left << PMtot[i];		// 3
			fProperties << setw(20) << left << rho[i];			// 4
			fProperties << setw(20) << left << Cp[i];			// 5
			fProperties << setw(20) << left << QReaction[i];	// 6
			fProperties << setw(20) << left << mu[i];			// 7
			fProperties << setw(20) << left << lambda[i];		// 8

			for(int k=1;k<=mix->NumberOfSpecies();k++)
				fProperties << setw(20) << left << Dm[i][k];	// 3
			fProperties << endl;
		}

		fProperties.close();
	}

	if (data->iSingleContributions == true)
	{
		std::string fileName;

		ofstream fSingleContributions;
		fileName = nameFolderAdditionalData + "/SingleContributions.out";
		openOutputFileAndControl(fSingleContributions, fileName);
		fSingleContributions.setf(ios::scientific);

		GnuPlotSingleContributionsInterface(fSingleContributions);

		for(i=2;i<=Np-1;i++)
		{
			// Print On File - Main data
			fSingleContributions << setw(20) << left<< 1.e2*grid.x[i];	// 1
			fSingleContributions << setw(20) << left<< T[i];			// 2

			for(int k=1;k<=single_contributions.Columns();k++)
				fSingleContributions << setw(20) << left << single_contributions[i][k];	// 3
			fSingleContributions << endl;
		}

		fSingleContributions.close();
	}

	bool iPrintVerbose = true;
	if (iPrintVerbose == true)
	{
		std::string fileName;

		ofstream fVerbose;
		fileName = nameFolderAdditionalData + "/Verbose.out";
		openOutputFileAndControl(fVerbose, fileName);
		fVerbose.setf(ios::scientific);

	//	GnuPlotVerboseInterface(fVerbose);

		for(i=1;i<=Np;i++)
		{
			// a. Calcolo delle frazioni massive e molari
			W.GetRow(i,&wVector);
			X.GetRow(i,&xVector);

			// Print On File - Main data
			fVerbose << setw(20) << left << 1.e2*grid.x[i];	// 1
			fVerbose << setw(20) << left << T[i];	// 2
			fVerbose << setw(20) << left << wVector.GetSumElements();	// 3
			fVerbose << setw(20) << left << xVector.GetSumElements();	// 4
			fVerbose << setw(20) << left << PMtot[i];	// 5
			fVerbose << setw(20) << left << rho[i];	// 6
			fVerbose << setw(20) << left << Cp[i];	// 7
			for(int k=1;k<=NC;k++)
				fVerbose << setw(20) << left << Cpk[i][k];	// 8...
			fVerbose << endl;
		}

		fVerbose.close();
	}
}
/*
void OpenSMOKE_Flame1D::GnuPlotPAHInterface(ofstream &fPAH)
{
	fPAH	<< setw(20) << left << "x[cm]"
			<< setw(20) << left << "T[K]";

	mix->pah_manager.GnuPlotInterface(fPAH, 3);
	fPAH << endl << endl;
}
*/
void OpenSMOKE_Flame1D::GnuPlotSootInterface(ofstream &fSoot)
{
	fSoot	<< setw(20) << left << "x[cm](1)"
			<< setw(20) << left << "T[K](2)";

    fSoot	<< setw(20) << left << "L_fv(3)"
			<< setw(20) << left << "L_x(4)"
			<< setw(20) << left << "L_y(5)"
			<< setw(20) << left << "L_rho[kg/m3](6)"
			<< setw(20) << left << "L_N[#/m3](7)"
			<< setw(20) << left << "L_H/C[-](8)"
			<< setw(20) << left << "L_d10[m](9)"
			<< setw(20) << left << "L_d32[m](10)"
			<< setw(20) << left << "L_dstd[m](11)";

    fSoot	<< setw(20) << left << "S_fv(12)"
			<< setw(20) << left << "S_x(13)"
			<< setw(20) << left << "S_y(14)"
			<< setw(20) << left << "S_rho[kg/m3](15)"
			<< setw(20) << left << "S_N[#/m3](16)";

	fSoot << endl << endl;
}

void OpenSMOKE_Flame1D::GnuPlotSootDistributionInterface(ofstream &fSoot)
{
	fSoot	<< setw(20) << left << "x[cm](1)"
			<< setw(20) << left << "bin(2)";

    fSoot	<< setw(20) << left << "d[m](3)"
			<< setw(20) << left << "m[?](4)"
			<< setw(20) << left << "V[?](5)"
			<< setw(20) << left << "log10d(6)"
			<< setw(20) << left << "log10m(7)"
			<< setw(20) << left << "log10V(8)"
			<< setw(20) << left << "dlog10d(9)"
			<< setw(20) << left << "dlog10m(10)"
			<< setw(20) << left << "dlog10V(11)";

    fSoot	<< setw(20) << left << "fv(12)"
			<< setw(20) << left << "x(13)"
			<< setw(20) << left << "y(14)"
			<< setw(20) << left << "rho[kg/m3](15)"
			<< setw(20) << left << "N[#/m3](16)";

    fSoot	<< setw(20) << left << "fvN(17)"
			<< setw(20) << left << "xN(18)"
			<< setw(20) << left << "yN(19)"
			<< setw(20) << left << "rhoN[](20)"
			<< setw(20) << left << "NN[](21)";

	fSoot << endl << endl;
}

void OpenSMOKE_Flame1D::GnuPlotLewisInterface(ofstream &fLewis)
{
	int fOutputCount = 1;

	PrintTagOnGnuplotLabel(20, fLewis, "x[cm]",			fOutputCount);
	PrintTagOnGnuplotLabel(20, fLewis, "T[K]",			fOutputCount);
	PrintTagOnGnuplotLabel(20, fLewis, "alfa[m2/s]",	fOutputCount);
	PrintTagOnGnuplotLabel(20, fLewis, "LeMeanMass",	fOutputCount);
	PrintTagOnGnuplotLabel(20, fLewis, "LeMeanMole",	fOutputCount);
	PrintTagOnGnuplotLabel(20, fLewis, "LeTildeMass",	fOutputCount);
	PrintTagOnGnuplotLabel(20, fLewis, "LeTildeMole",	fOutputCount);
	
	for(int j=1;j<=NC;j++)
		PrintTagOnGnuplotLabel(20, fLewis, mix->names[j],fOutputCount);
	fLewis << endl << endl;
}

void OpenSMOKE_Flame1D::GnuPlotSoretInterface(ofstream &fSoret)
{
	int fOutputCount = 1;

	PrintTagOnGnuplotLabel(20, fSoret, "x[cm]",			fOutputCount);
	PrintTagOnGnuplotLabel(20, fSoret, "T[K]",			fOutputCount);
	
	for(int j=1;j<=NC;j++)
		if (mix->M[j] <= 5.)	
			PrintTagOnGnuplotLabel(20, fSoret, mix->names[j],fOutputCount);
	fSoret << endl << endl;
}

void OpenSMOKE_Flame1D::GnuPlotSingleContributionsInterface(ofstream &fSingleContributions)
{
	int fOutputCount = 1;

	PrintTagOnGnuplotLabel(20,fSingleContributions, "x[cm]",			fOutputCount);
	PrintTagOnGnuplotLabel(20, fSingleContributions, "T[K]",			fOutputCount);

	for(int j=1;j<=data->index_SingleContributions.Size();j++)
	{
		PrintTagOnGnuplotLabel(20, fSingleContributions, data->names_SingleContributions[j-1] + "_D", fOutputCount);
		PrintTagOnGnuplotLabel(20, fSingleContributions, data->names_SingleContributions[j-1] + "_C", fOutputCount);
		PrintTagOnGnuplotLabel(20, fSingleContributions, data->names_SingleContributions[j-1] + "_R", fOutputCount);
	}
	fSingleContributions << endl << endl;
}

void OpenSMOKE_Flame1D::GnuPlotMixturePropertiesInterface(ofstream &fProperties)
{
	int fOutputCount = 1;

	PrintTagOnGnuplotLabel(20, fProperties, "x[cm]",			fOutputCount);
	PrintTagOnGnuplotLabel(20, fProperties, "T[K]",				fOutputCount);
	PrintTagOnGnuplotLabel(20, fProperties, "MW[kg/kmol]",		fOutputCount);
	PrintTagOnGnuplotLabel(20, fProperties, "rho[kg/m3]",		fOutputCount);
	PrintTagOnGnuplotLabel(20, fProperties, "Cp[J/kg/K]",		fOutputCount);
	PrintTagOnGnuplotLabel(20, fProperties, "QR[J/m3/s]",		fOutputCount);
	PrintTagOnGnuplotLabel(20, fProperties, "mu[Pa.s]",			fOutputCount);
	PrintTagOnGnuplotLabel(20, fProperties, "lambda[W/m/K]",	fOutputCount);

	for(int j=1;j<=mix->NumberOfSpecies();j++)
		PrintTagOnGnuplotLabel(20, fProperties, "Dm_" + mix->names[j]+"_[m2/s]", fOutputCount);
	fProperties << endl << endl;
}

void OpenSMOKE_Flame1D::GnuPlotFormationRatesInterface(ofstream &fFormationRates)
{
	int fOutputCount = 1;

	PrintTagOnGnuplotLabel(20, fFormationRates, "x[cm]", fOutputCount);
	PrintTagOnGnuplotLabel(20, fFormationRates, "T[K]",	 fOutputCount);

	for(int j=1;j<=data->index_formation_rates.Size();j++)
		PrintTagOnGnuplotLabel(20, fFormationRates, mix->names[data->index_formation_rates[j]], fOutputCount);
	fFormationRates << endl << endl;
}

void OpenSMOKE_Flame1D::GnuPlotReactionRatesInterface(ofstream &fReactionRates)
{
	int fOutputCount = 1;

	PrintTagOnGnuplotLabel(20, fReactionRates, "x[cm]", fOutputCount);
	PrintTagOnGnuplotLabel(20, fReactionRates, "T[K]",	fOutputCount);

	for(int j=1;j<=data->index_reaction_rates.Size();j++)
	{
		stringstream number_reaction;
		number_reaction << data->index_reaction_rates[j];
		PrintTagOnGnuplotLabel(20, fReactionRates, "r_" + number_reaction.str(),	fOutputCount);
	}
	fReactionRates << endl << endl;
}


void OpenSMOKE_Flame1D::GnuPlotInterfaceUnsteady()
{
	int j;

	{
		int fOutputCount = 1;
		PrintTagOnGnuplotLabel(20, fUnsteady, "t[s]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "t[-]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "K[Hz]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "x[cm]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "U[kg/m2/s]", fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "u[cm/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "H[Pa/m2]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "G[kg/m3/s]", fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "T[K]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "rho[kg/m3]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "MW[kg/kmol]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteady, "Qreac[W/m3]",	fOutputCount);

		for(j=1;j<=data->iOut.Size();j++)
			PrintTagOnGnuplotLabel(20, fUnsteady, data->nameOutput[j] + "_x", fOutputCount);
		fUnsteady << endl;
	}

	{
		int fOutputCount = 1;
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "t[s]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "t[-]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "K[Hz]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "KS[Hz]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "KJ[Hz]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "uF[cm/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "uO[cm/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "T[K]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "dummy",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "dummy",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "alfa[m2/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyMax, "Qreac[W/m3]",	fOutputCount);

		for(j=1;j<=data->iOut.Size();j++)
			PrintTagOnGnuplotLabel(20, fUnsteadyMax, data->nameOutput[j] + "_x", fOutputCount);
		fUnsteadyMax << endl;
	}

	// Unsteady QMOM Max
	{
		int fOutputCount = 1;
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "t[s]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "t[-]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "K[Hz]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "KS[Hz]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "KJ[Hz]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "uF[cm/s]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "uO[cm/s]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "T[K]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "fv[-]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "As[1/m]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fUnsteadyQMOMMax, "Lm[nm]",		fOutputCount);
		fUnsteadyQMOMMax << endl;
	}
}

void OpenSMOKE_Flame1D::GnuPlotInterface(ofstream &fOutput)
{
	int j;

	if ( data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED || data->kind_of_flame == FLAME1D_PHYSICS_TWIN )
	{
		int fOutputCount = 1;
		PrintTagOnGnuplotLabel(20, fOutput, "x[cm]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "U[kg/m2/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "u[cm/s]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "G[kg/m3/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "H[Pa/m2]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "T[K]",			fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "rho[kg/m3]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "mu[kg/m/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "k[W/m/K]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "Cp[J/kg/K]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "Qreac[W/m3]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "Qrad[W/m3]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "Z[-]",			fOutputCount);
	
		// Elements
		for(j=1;j<=mix->NumberOfElements();j++)
			PrintTagOnGnuplotLabel(20, fOutput,  mix->list_of_elements[j-1] + "_x",	fOutputCount);
		for(j=1;j<=mix->NumberOfElements();j++)
			PrintTagOnGnuplotLabel(20, fOutput,  mix->list_of_elements[j-1] + "_w",	fOutputCount);

		// Species
		for(j=1;j<=data->iOut.Size();j++)
			PrintTagOnGnuplotLabel(20, fOutput,  data->nameOutput[j] + "_x",	fOutputCount);
		for(j=1;j<=data->iOut.Size();j++)
			PrintTagOnGnuplotLabel(20, fOutput,  data->nameOutput[j] + "_w",	fOutputCount);
	
		if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM ) 
		{
			std::string fileName;
			fileName = nameFolderAdditionalData + "/Soot_QMOM.out";
			openOutputFileAndControl(fSootQMOM, fileName);
			fSootQMOM.setf(ios::scientific);

			qmom.soot.write_label_on_file(fSootQMOM, qmom.N*2);
		}
	}

	else if (  data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED ) 
	{
		int fOutputCount = 1;
		PrintTagOnGnuplotLabel(20, fOutput, "x[cm]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "u[cm/s]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "m[kg/s]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "m[kg/m2/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "A[m2]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "T[K]",			fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "rho[kg/m3]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "mu[kg/m/s]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "k[W/m/K]",		fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "Cp[J/kg/K]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "Qreac[W/m3]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "Qrad[W/m3]",	fOutputCount);
		PrintTagOnGnuplotLabel(20, fOutput, "Z[-]",			fOutputCount);

		// Elements
		for(j=1;j<=mix->NumberOfElements();j++)
			PrintTagOnGnuplotLabel(20, fOutput,  mix->list_of_elements[j-1] + "_x",	fOutputCount);
		for(j=1;j<=mix->NumberOfElements();j++)
			PrintTagOnGnuplotLabel(20, fOutput,  mix->list_of_elements[j-1] + "_w",	fOutputCount);

		// Species
		for(j=1;j<=data->iOut.Size();j++)
			PrintTagOnGnuplotLabel(20, fOutput,  data->nameOutput[j] + "_x",	fOutputCount);
		for(j=1;j<=data->iOut.Size();j++)
			PrintTagOnGnuplotLabel(20, fOutput,  data->nameOutput[j] + "_w",	fOutputCount);
	
		if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM )
		{
			std::string fileName;
			fileName = nameFolderAdditionalData + "/Soot_QMOM.out";
			openOutputFileAndControl(fSootQMOM, fileName);
			fSootQMOM.setf(ios::scientific);

			qmom.soot.write_label_on_file(fSootQMOM, qmom.N*2);
		}
	}

	fOutput << endl;
}

void OpenSMOKE_Flame1D::DAE_ODE_myPrint(BzzVector &v, double time)
{
	if(data->iterationVideoCounter == data->nStepsVideo)
	{
		cout.setf(ios::scientific);
		cout << data->iteration << "\t";
		
		if (data->iUnsteady == true)
		{
			cout << time								<< "\t";
			cout << unsteady.nPeriods					<< "\t";
			cout << unsteady.phase						<< "\t";
			cout << Tmax								<< "\t";
			cout <<  (nGeometry-1.)*U[1]/rho[1]  *100.	<< "\t";
			cout << -(nGeometry-1.)*U[Np]/rho[Np]*100.	<< "\t";
			cout << unsteady.K							<< "\t";
			cout << unsteady.KSeshadri					<< "\t";
			cout << endl;
		}

		else
		{
			cout << time << "\t";
			cout << Tmax << "\t";

			if(data->iDepositionWall == true)
			{
				double sum_max = -1e16;
				
				for(int i=1;i<=Np;i++)
				{ 
					double local_sum = 0.;
					for(int j=1;j<=mix->polimiSoot->bin_indices().Size();j++)
					{
						const int jj = mix->polimiSoot->bin_indices()[j];
						local_sum += W[i][jj];
					}
					if (local_sum > sum_max) sum_max = local_sum;
				}

				cout << sum_max << "\t";
				cout << soot_deposition << "\t";
				cout << -data->VO*100. << "\t";
			}

		}
		
		if (data->iFixedTemperature > 0 )
		{
			cout << U[1]*100.	<< "\t";
			cout << H[1]		<< "\t";
		}
		
		cout << endl;

		data->iterationVideoCounter = 0;
	}

	if(unsteady.onPrint == 1 && data->iUnsteady == true && data->iterationVideoCounter == 0)
	{
		for(int i=1;i<=Np;i++)
		{
			fUnsteady << setw(20) << left << time;								// 1. time
			fUnsteady << setw(20) << left << time/unsteady.T;					// 2. time'			
			fUnsteady << setw(20) << left << unsteady.K;						// 3. K
			fUnsteady << setw(20) << left << grid.x[i];							// 4. x
			fUnsteady << setw(20) << left << U[i];								// 5. U
			fUnsteady << setw(20) << left << 100.*(nGeometry-1.)*U[i]/rho[i];	// 6. velocity
			fUnsteady << setw(20) << left << H[i];								// 7. eigenvalue
			fUnsteady << setw(20) << left << G[i];								// 8. mass flow rate
			fUnsteady << setw(20) << left << T[i];								// 9. temperature
			fUnsteady << setw(20) << left << rho[i];							// density
			fUnsteady << setw(20) << left << PMtot[i];							// molecular weight
			fUnsteady << setw(20) << left << QReaction[i];						// heat release
			
			for(int k=1;k<=NC;k++)
				fUnsteady << setw(20) << left << X[i][k];						// 10... mole fractions
			fUnsteady << endl;
		}

		fUnsteady << endl;
		data->iterationFileCounter = 0;

		int iMax; T.Max(&iMax);
		{
			fUnsteadyMax << setw(20) << left << time;								// 1. time
			fUnsteadyMax << setw(20) << left << time/unsteady.T;					// 2. time			
			fUnsteadyMax << setw(20) << left << unsteady.K;							// 3. K
			fUnsteadyMax << setw(20) << left << unsteady.KSeshadri;					// 4. KSeshadri
			fUnsteadyMax << setw(20) << left << sqrt(-4.*H[1]/rho[1]);				// 5. KJackson
			fUnsteadyMax << setw(20) << left << 100.*(nGeometry-1.)*U[1]/rho[1];	// 6. velocity
			fUnsteadyMax << setw(20) << left << 100.*(nGeometry-1.)*U[Np]/rho[Np];	// 7. velocity
			fUnsteadyMax << setw(20) << left << T.Max();							// 8. temperature
			fUnsteadyMax << setw(20) << left << 0;									// 9. dummy
			fUnsteadyMax << setw(20) << left << 0;									// 10. dummy
			fUnsteadyMax << setw(20) << left << lambda[iMax]/Cp[iMax]/rho[iMax];	// 11. Thermal diffusivity [m2/s]
			fUnsteadyMax << setw(20) << left << QReaction[iMax];					// 12. Maximum heat release [W/m3]
				
			for(int k=1;k<=NC;k++)
			{
				int jMax;
				BzzVector row = X.GetColumn(k);
				row.Max(&jMax);
				fUnsteadyMax << setw(20) << left << X[jMax][k];						// 9. ecc mole fractions
			}
			fUnsteadyMax << endl;
		}

		// Unsteady QMOM max
		if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM )
		{
			double fvMax	= 0.;
			double AsMax	= 0.;
			double LMeanMax = 0.;
			for(int i=2;i<=Ni;i++)
			{
				BzzVector momentsVector = moments.GetRow(i);
				qmom.updateMoments(momentsVector);

				double valueFluctuations = 0;
				double mu=1.0;
				double tr=1.0;
				qmom.updateData(T[i], data->P_Pascal, rho[i], mu, tr, W[i][qmom.jO2], W[i][qmom.jC2H2], W[i][qmom.jOH], valueFluctuations );
	
				if (fvMax < qmom.soot.fv)	
				{
					fvMax		= qmom.soot.fv;
					AsMax		= qmom.soot.As;
					LMeanMax	= qmom.soot.Lmean;
				}
			}

			fUnsteadyQMOMMax << setw(20) << left << time;								// 1. time
			fUnsteadyQMOMMax << setw(20) << left << time/unsteady.T;					// 2. time			
			fUnsteadyQMOMMax << setw(20) << left << unsteady.K;							// 3. K
			fUnsteadyQMOMMax << setw(20) << left << unsteady.KSeshadri;					// 4. KSeshadri
			fUnsteadyQMOMMax << setw(20) << left << sqrt(-4.*H[1]/rho[1]);				// 5. KJackson
			fUnsteadyQMOMMax << setw(20) << left << 100.*(nGeometry-1.)*U[1]/rho[1];	// 6. velocity
			fUnsteadyQMOMMax << setw(20) << left << 100.*(nGeometry-1.)*U[Np]/rho[Np];	// 7. velocity
			fUnsteadyQMOMMax << setw(20) << left << T.Max();							// 8. temperature
			fUnsteadyQMOMMax << setw(20) << left << fvMax;								// 9. max fv soot
			fUnsteadyQMOMMax << setw(20) << left << AsMax;								// 10. max As soot
			fUnsteadyQMOMMax << setw(20) << left << LMeanMax*1.e9;						// 11. max LMean soot
			fUnsteadyQMOMMax << endl;
		}
	}

	if(data->iterationBackUpCounter == data->nStepsBackUp)
	{
		stringstream iteration;
		std::string iterationComplete;
		std::string fileName_BackUp;
		std::string fileName_BackUpInputData;
		
		if		(data->iteration > 9999)	iterationComplete = "";
		else if (data->iteration > 999)	iterationComplete = "0";
		else if (data->iteration > 99)	iterationComplete = "00";
		else if (data->iteration > 9)	iterationComplete = "000";
		else							iterationComplete = "0000";
		
		iteration << data->iteration;
		iterationComplete += iteration.str();

		fileName_BackUp				= nameFileBackupInputData + "_" + iterationComplete;
		fileName_BackUpInputData	= nameFileBackupInputData + "_" + iterationComplete;
		
		cout << "Printing BackUp Data on file: " << fileName_BackUp << ".out" << endl;
		printBackUpOnlyInputData(fileName_BackUpInputData);
		printBackUpOnlyData(fileName_BackUp);

		data->iterationBackUpCounter = 0;
	}

	data->iteration++;
	data->iterationVideoCounter++;
	data->iterationFileCounter++;
	data->iterationBackUpCounter++;

	if (data->iUnsteady == true)
	{
		unsteady_boundary_conditions(time);
		unsteady.update_time_target(time);
	}

	double TMax = T.Max();
	if (TMax <= Constants::TMinExtinction-500)
	{
		cout << "Minimum temperature reached: the flame is extinguishing..." << endl;
		cout << "Simulation is being interrupted..." << endl;
		bzzStop = 1;
	}

	if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED && data->iRobustTemperature == true)
	{
		cout << data->counterUnchanged << endl;
		cout << TMax << " " << data->TMaxOld << " " << fabs(TMax-data->TMaxOld)/TMax << endl;
		getchar();

		if ( fabs(TMax-data->TMaxOld)/TMax < 1.e-6)	
			data->counterUnchanged++;
		else
		{
			data->counterUnchanged = 0;
			data->TMaxOld = TMax;
		}

		if (data->counterUnchanged >= 1000)
		{
			cout << "Maximum number of iterations..." << endl;
			cout << "Simulation is skipped..." << endl;
			cout << "TMax Actual: " << TMax << endl;
			cout << "TMax Old:    " << data->TMaxOld << endl;
			cout << "Error(%):    " << fabs(TMax-data->TMaxOld)/TMax*100 << endl;
			bzzStop = 1;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								NLS Systems - Opposed											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::nonLinearSystem_Opposed_All(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, -1, dummy, 0, 0);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------

		give_Opposed_DU_DG_DH();
		give_Opposed_DT(0.);
		give_Opposed_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_ALL, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Twin_All(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, -1, dummy, 0, 0);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------

		give_Twin_DU_DG_DH();
		give_Twin_DT(0.);
		give_Twin_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_ALL, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Reduced_Opposed_All(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables_Reduced(OPPOSED_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		//properties(data->iHOT, -1, dummy, 0, 0);
		properties(data->iHOT);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------

		give_Opposed_DU_DG_DH();
		give_Opposed_DT(0.);
		give_Opposed_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals_Reduced(OPPOSED_ALL, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Reduced_Twin_All(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables_Reduced(TWIN_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		//properties(data->iHOT, -1, dummy, 0, 0);
		properties(data->iHOT);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------

		give_Twin_DU_DG_DH();
		give_Twin_DT(0.);
		give_Twin_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals_Reduced(TWIN_ALL, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Opposed_NoEnergy(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_NO_ENERGY, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		// TODO for correction kinetics
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DU_DG_DH();
		give_Opposed_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_NO_ENERGY, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Twin_NoEnergy(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_NO_ENERGY, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		// TODO for correction kinetics
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_DU_DG_DH();
		give_Twin_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_NO_ENERGY, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Opposed_OnlyUGH(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_ONLY_MOMENTUM, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DU_DG_DH();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_ONLY_MOMENTUM, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Twin_OnlyUGH(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_ONLY_MOMENTUM, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_DU_DG_DH();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_ONLY_MOMENTUM, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Opposed_OnlyMassFractions(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_ONLY_MASS_FRACTIONS, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 6. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_ONLY_MASS_FRACTIONS, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Twin_OnlyMassFractions(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_ONLY_MASS_FRACTIONS, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 6. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_ONLY_MASS_FRACTIONS, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Opposed_ColdReduced(BzzVector &x, BzzVector &f)
{	
	int i, j;

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_COLD_REDUCED, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		
	// -------------------------------------------------------------------------------------		
	// 3. Variabili Ausiliarie
	// -------------------------------------------------------------------------------------
		rhow[1] = rho[1];
		for(i=1;i<=Np-1;i++)
		{
			rhoe[i] = 0.50*(rho[i]+rho[i+1]);
			rhow[i+1] = rhoe[i];
		}
		rhoe[Np] = rho[Np];

	// -------------------------------------------------------------------------------------		
	// 4. Termini convettivi
	// -------------------------------------------------------------------------------------
		for(j=1;j<=data->nCold-1;j++)
		{ 
			W.GetColumn(data->iX[j], &auxNp);
			grid.FirstDerivative(data->iDerW, U, auxNp, diffWscalar);
			diffW.SetColumn(data->iX[j], diffWscalar);
		}

	// --------------------------------------------------------------------------------------
	// 5. Equazioni
	// --------------------------------------------------------------------------------------
		for(j=1;j<=data->nCold-1;j++)
			dW[1][data->iX[j]] = rho[1]*((nGeometry-1.)*U[1]*urho[1]*W[1][data->iX[j]] + vStar[1][data->iX[j]]) - BCW_C[data->iX[j]];
		dW[1][data->iX[data->nCold]] = W[1][data->iX[data->nCold]] - sumW[1];

		for(i=2;i<=Ni;i++)
		{
			for(j=1;j<=data->nCold-1;j++)
				dW[i][data->iX[j]] = -((nGeometry-1.)*U[i] * diffW[i][data->iX[j]]
				               - R[i][data->iX[j]]
						       + (rhoe[i]*vStar[i][data->iX[j]] - rhow[i]*vStar[i-1][data->iX[j]])*grid.udxc_over_2[i]
							 ) * urho[i];
			dW[i][data->iX[data->nCold]] = W[i][data->iX[data->nCold]] - sumW[i];
		}

		for(j=1;j<=data->nCold-1;j++)
			dW[Np][data->iX[j]] = rho[Np] * ((nGeometry-1.)*U[Np] * urho[Np] * W[Np][data->iX[j]] + vStar[Np][data->iX[j]]) - BCW_O[data->iX[j]];
		dW[Np][data->iX[data->nCold]] = W[Np][data->iX[data->nCold]] - sumW[Np];


	// --------------------------------------------------------------------------------------
	// 6. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_COLD_REDUCED, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Twin_ColdReduced(BzzVector &x, BzzVector &f)
{	
	ErrorMessage("TWIN_COLD_REDUCED not yet implemented!");
}

void OpenSMOKE_Flame1D::nonLinearSystem_Opposed_OnlyT(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		T = x;
	
	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, -1, dummy, 0, 0);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DT(0.);

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		f = dT;
}

void OpenSMOKE_Flame1D::nonLinearSystem_Twin_OnlyT(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		T = x;
	
	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, -1, dummy, 0, 0);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_DT(0.);

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		f = dT;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								DAE SYSTEMs	- Opposed											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenSMOKE_Flame1D::unsteady_boundary_conditions(double &time)
{
	if (data->iUnsteady == true)
	{
		unsteady.update_boundary_conditions(time, UC, UO, data->TO, rho[1], rho[Np], T.Max(), WC, WO);
		for(int i=1;i<=NC;i++)
		{
			BCW_C[i] =  (nGeometry-1.)*UC*WC[i];
			BCW_O[i] =  (nGeometry-1.)*UO*WO[i];
		}
	}
}

void OpenSMOKE_Flame1D::DAESystem_Opposed_All(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();
		unsteady_boundary_conditions(t);

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		// TODO
		properties(data->iHOT, DAE_Opposed_ALL.jacobianIndex, DAE_Opposed_ALL.jacobianVariables, NC+4, 4);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DU_DG_DH();
		give_Opposed_DT(t);
		give_Opposed_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_ALL, f);

}

void OpenSMOKE_Flame1D::DAESystem_Opposed_OnlyMomentum(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_ONLY_MOMENTUM, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DU_DG_DH();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_ONLY_MOMENTUM, f);
}

void OpenSMOKE_Flame1D::DAESystem_Twin_All(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();

		unsteady_boundary_conditions(t);

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		properties(data->iHOT, DAE_Twin_ALL.jacobianIndex, DAE_Twin_ALL.jacobianVariables, NC+4, 4);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_DU_DG_DH();
		give_Twin_DT(t);
		give_Twin_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_ALL, f);
}

void OpenSMOKE_Flame1D::DAESystem_Opposed_SOOT_ALL(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();

		unsteady_boundary_conditions(t);

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_SOOT_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		properties(data->iHOT, DAE_Opposed_SOOT_ALL.jacobianIndex, DAE_Opposed_SOOT_ALL.jacobianVariables, NC+6, 4);
		compute_vStarOpposed(); 
		
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_SOOT();
		give_Opposed_DU_DG_DH();
		give_Opposed_DT(t);
		give_Opposed_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_SOOT_ALL, f);
}

void OpenSMOKE_Flame1D::DAESystem_Twin_SOOT_ALL(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();

		unsteady_boundary_conditions(t);

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_SOOT_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		properties(data->iHOT, DAE_Twin_SOOT_ALL.jacobianIndex, DAE_Twin_SOOT_ALL.jacobianVariables, NC+6, 4);
		compute_vStarOpposed(); 
		
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_SOOT();
		give_Twin_DU_DG_DH();
		give_Twin_DT(t);
		give_Twin_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_SOOT_ALL, f);
}

void OpenSMOKE_Flame1D::DAESystem_Opposed_NoMomentum(BzzVector &x, double t, BzzVector &f)
{

	Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_NO_MOMENTUM, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)
		properties(data->iHOT,DAE_Opposed_NOMOMENTUM.jacobianIndex, DAE_Opposed_NOMOMENTUM.jacobianVariables, NC+4, 1);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DT(t);
		give_Opposed_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_NO_MOMENTUM, f);
}

void OpenSMOKE_Flame1D::DAESystem_Twin_NoMomentum(BzzVector &x, double t, BzzVector &f)
{

	Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_NO_MOMENTUM, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)
		properties(data->iHOT,DAE_Twin_NOMOMENTUM.jacobianIndex, DAE_Twin_NOMOMENTUM.jacobianVariables, NC+4, 1);
		compute_vStarOpposed();
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_DT(t);
		give_Twin_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_NO_MOMENTUM, f);
}

void OpenSMOKE_Flame1D::DAESystem_Opposed_NoEnergy(BzzVector &x, double t, BzzVector &f)
{
	Tmax=T.Max();

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_NO_ENERGY, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le propriet(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le propriet(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le propriet(bisogna memorizzare)

		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed(); 
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DU_DG_DH();
		give_Opposed_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_NO_ENERGY, f);
}

void OpenSMOKE_Flame1D::DAESystem_Twin_NoEnergy(BzzVector &x, double t, BzzVector &f)
{
	Tmax=T.Max();

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_NO_ENERGY, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le propriet(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le propriet(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le propriet(bisogna memorizzare)
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed(); 
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_DU_DG_DH();
		give_Twin_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_NO_ENERGY, f);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									GRID REFINING												   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_recursions(int Np, BzzMatrixInt &A, BzzVectorInt &V, BzzVector &weight)
{
	ChangeDimensions(Np, &V);
	for(int j=1;j<=A.Rows();j++)
		for(int i=1;i<=A.Columns();i++)
		{
			if(A[j][i]>0)
				V[A[j][i]] = V[A[j][i]] + weight[j];
		}
}

bool OpenSMOKE_Flame1D::newPoints(const std::string string_kind, char index)
{
	int  i;
	int  nDim;
	bool iOnlyTemperature;
	bool iAddedPoints = false;

	if		(string_kind == "TEMPERATURE")	iOnlyTemperature = true;
	else if (string_kind == "ALL")			iOnlyTemperature = false;
	else	ErrorMessage("Only TEMPERATURE || ALL options are available for grid refinement");

	if (iOnlyTemperature == true)	nDim = 1;
	if (iOnlyTemperature == false)	nDim = NC+1;

	BzzVector weight_D(nDim);
	BzzVector weight_G(nDim);
	
	weight_D = 1.;
	weight_G = 1.;
	
	cout << "Adding Points ("<< index << ") using " << string_kind << " profile!"<< endl;

	// Difference
	if (index == 'D')
	{
		BzzMatrixInt pointMatrix(nDim, data->nDiff);
		BzzVectorInt pointVector;

		if (iOnlyTemperature == false)
		{
			for(int j=1;j<=NC;j++)
			{
				W.GetColumn(j, &auxNp);
				BzzVectorInt pointList = grid.QueryNewPointsDifference(data->nDiff, data->deltaDiff, auxNp);

				if (pointList.Size()>0)
					cout << "Requiring new points (D) - Species: " << j << "\t - New points: " << pointList.Size() << endl;

				for(i=1;i<=pointList.Size();i++)
					pointMatrix[j][i] = pointList[i];
			}		
		}

		BzzVectorInt pointList = grid.QueryNewPointsDifference(data->nDiff, data->deltaDiff, T);
		if (pointList.Size()>0)
				cout << "Requiring new points (D) - Temperature:\t - New points: " << pointList.Size() << endl;
		for(i=1;i<=pointList.Size();i++)
			pointMatrix[nDim][i] = pointList[i];

		calculate_recursions(grid.Np, pointMatrix, pointVector, weight_D);

		BzzVectorInt list_of_new_points;
		for(i=1;i<=data->nDiff;i++)
		{
			int index;
			int maxValue = pointVector.Max(&index);
			if(maxValue > 0)
			{
				list_of_new_points.Append(index);
				pointVector[index] = 0;
			}

			cout << " ** Grid point: " << index << "\t - Number of requests: " << maxValue << endl;
		}

		if (list_of_new_points.Size() > 0)	
		{
			AddPoints(list_of_new_points);
			iAddedPoints = true;
		}
		
		cout << " New points added: " << pointList.Size() << endl;
	}

	// Gradient
	if (index == 'G')
	{
		BzzMatrixInt pointMatrix(nDim, data->nGrad);
		BzzVectorInt pointVector;
		
		if (iOnlyTemperature == false)
		{
			for(int j=1;j<=NC;j++)
			{
				W.GetColumn(j, &auxNp);
				BzzVectorInt pointList = grid.QueryNewPointsGradient(data->nGrad, data->deltaGrad, auxNp);
				
				if (pointList.Size()>0)
					cout << "Requiring new points (G) - Species: " << j << "\t - New points: " << pointList.Size() << endl;

				for(i=1;i<=pointList.Size();i++)
					pointMatrix[j][i] = pointList[i];
			}		
		}

		BzzVectorInt pointList = grid.QueryNewPointsGradient(data->nGrad, data->deltaGrad, T);
		if (pointList.Size()>0)
				cout << "Requiring new points (G) - Temperature:\t - New points: " << pointList.Size() << endl;
		for(i=1;i<=pointList.Size();i++)
			pointMatrix[nDim][i] = pointList[i];

		calculate_recursions(grid.Np, pointMatrix, pointVector, weight_G);

		BzzVectorInt list_of_new_points;
		for(i=1;i<=data->nGrad;i++)
		{
			int index;
			int maxValue = pointVector.Max(&index);
			if(maxValue > 0)
			{
				list_of_new_points.Append(index);
				pointVector[index] = 0;
			}

			cout << " ** Grid point: " << index << "\t - Number of requests: " << maxValue << endl;
		}

		if (list_of_new_points.Size() > 0)	
		{
			AddPoints(list_of_new_points);
			iAddedPoints = true;
		}

		cout << " New points added: " << pointList.Size() << endl;
	}

	return iAddedPoints;
}


void OpenSMOKE_Flame1D::doubleTheGrid()
{
	// For Flame Speed Problems
	data->iFixedTemperature += (data->iFixedTemperature-1);

	// Refine the grid
	grid.RefineDouble();
	Np = grid.Np;
	Ni = grid.Ni;

	// Linear Interpolation
	grid.DoubleField(U);
	grid.DoubleField(G);
	grid.DoubleField(H);
	grid.DoubleField(T);
	grid.DoubleField(W);

	if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM )
		grid.DoubleField(moments);
		
	if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )	
	{
		grid.DoubleField(phiN);
		grid.DoubleField(phiM);
	}

	// Updating Properties
	allocate_only_Np();
	BzzVectorInt dummy;
	properties(1, -1, dummy, 0, 0);
	updatingProfiles();

	cout << "Grid correctly refined" << endl; 
}


void OpenSMOKE_Flame1D::AddPoints(BzzVectorInt &listPoints)
{
	int i;

	// For Flame Speed Problems
	int count_shift = 0;
	for (i=1;i<=listPoints.Size();i++)
		if (listPoints[i] < data->iFixedTemperature)	count_shift++;
	data->iFixedTemperature += count_shift;

	// Refine the grid
	Sort(&listPoints);
	for (i=1;i<=listPoints.Size();i++)		// This is correct only if listPoints is sorted
		grid.Refine(listPoints[i]+(i-1));	// Min --> Max

	Np = grid.Np;
	Ni = grid.Ni;

	// Linear Interpolation
	grid.AddPointsField(U, listPoints);
	grid.AddPointsField(G, listPoints);
	grid.AddPointsField(H, listPoints);
	grid.AddPointsField(T, listPoints);
	grid.AddPointsField(W, listPoints);
	
	if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM )
		grid.AddPointsField(moments, listPoints);
		
	if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )	
	{
		grid.AddPointsField(phiN, listPoints);
		grid.AddPointsField(phiM, listPoints);
	}

	// Updating Properties
	allocate_only_Np();
	BzzVectorInt dummy;
	properties(1, -1, dummy, 0, 0);
	updatingProfiles();
}

void OpenSMOKE_Flame1D::refineGridPeak(double fraction)
{
	int i;
	BzzVectorInt listAddPoints;

	Tmax = T.Max();
	for(i=1;i<=Np;i++)
		if (T[i]>=(fraction*Tmax)) listAddPoints.Append(i);		

//	for(i=1;i<=Np;i++)
//		if (grid.x[i]>=0.55e-2 && grid.x[i]<=0.90e-2) listAddPoints.Append(i);		

	AddPoints(listAddPoints);
}

void OpenSMOKE_Flame1D::refineFlameBase()
{
	int i;
	int iCentered;
	BzzVectorInt listAddPoints;

	double Tmax = T.Max();
	double Tmin = T.Min();

	for(i=1;i<=Np;i++)
		if (T[i] >=  Tmin+0.05*(Tmax-Tmin)) 
		{
			iCentered = i;
			break;
		}

	cout << "Centered point: " << iCentered << endl;

	listAddPoints.Append(iCentered-3);		
	listAddPoints.Append(iCentered-2);		
	listAddPoints.Append(iCentered-1);		
	listAddPoints.Append(iCentered);		
	listAddPoints.Append(iCentered+1);		
	listAddPoints.Append(iCentered+2);		

	AddPoints(listAddPoints);
}

void OpenSMOKE_Flame1D::refineBase()
{
	int i;
	int iA, iB;
	BzzVectorInt listAddPoints;

	double Tmax = T.Max();
	double Tmin = T.Min();

	for(i=1;i<=Np;i++)
		if (T[i] >=  T[1]+0.05*(data->fixedTemperature-T[1])  )
		{
			iA = i>1 ? (i-1) : 1;
			break;
		}

	for(i=iA+1;i<=Np;i++)
		if (T[i] >=  Tmax-0.08*(Tmax-Tmin)) 
		{
			iB = i;
			break;
		}

	for(i=iA;i<=iB;i++)
		listAddPoints.Append(i);		

	AddPoints(listAddPoints);
}

void OpenSMOKE_Flame1D::refineAttachPoint()
{
	int i;
	int iA, iB;
	BzzVectorInt listAddPoints;

	double Tmax = T.Max();
	double Tmin = T.Min();

	for(i=1;i<=Np;i++)
		if (T[i] >=  T[1]+0.05*(data->fixedTemperature-T[1])  )
		{
			iA = i>1 ? (i-1) : 1;
			break;
		}

	for(i=iA+1;i<=Np;i++)
		if (T[i] >=  data->fixedTemperature + 0.05*(Tmax-Tmin)) 
		{
			iB = i;
			break;
		}

	for(i=iA;i<=iB;i++)
		listAddPoints.Append(i);		

	AddPoints(listAddPoints);
}

void OpenSMOKE_Flame1D::refineStagnationPlane(const double fraction)
{
	for(int i=1;i<=Np;i++)
		if (U[i] <= 0.)
		{
			int iA = i-1;
			int iB = i;
			double m = (U[iB]-U[iA])/(grid.x[iB]-grid.x[iA]);
			double q = U[iB]-m*grid.x[iB];
			double xstagnation = -q/m;
			double xA = xstagnation - 0.50*fraction*grid.L;
			double xB = xstagnation + 0.50*fraction*grid.L;
			xA = max(0., xA);
			xB = min(grid.L, xB);

			cout << "Stagnation plane: " << xstagnation*100. << " cm" << endl;
			cout << "Boundary A:       " << xA*100. << " cm" << endl;
			cout << "Boundary B:       " << xB*100. << " cm" << endl;

			refineFlameBase(xA*100, xB*100);

			break;
		}
}


OpenSMOKE_AdaptiveGrid *ptAdaptiveGrid;
void OdeAdaptiveGrid(BzzVector &x, double t, BzzVector &f)
{
	ptAdaptiveGrid->Ode(x,t,f);
}

void OpenSMOKE_Flame1D::adaptGrid(const int iOption)
{
	int i;
	BzzVector xAdapted = grid.x;

	cout << "Adapting grid..." << endl;

	// 1 solo prima parte (tutto)
	// 2 solo seconda parte (tutto)
	// 3 prima e seconda parte (tutto)

	// 11 solo prima parte (parziale)
	// 22 solo seconda parte (parziale)
	// 33 prima e seconda parte (parziale)
	
	if (data->iFixedTemperature>0)
	{
		if (iOption == 0)
		{
			OpenSMOKE_AdaptiveGrid adaptive_grid;
			ptAdaptiveGrid = &adaptive_grid;
			BzzVector xOld(Np);
			BzzVector xNew(Np);
			BzzVector TOld(Np);
			
			xOld = grid.x;
			TOld = T;
			
			adaptive_grid.Setup(xOld, TOld);
			adaptive_grid.SetModel(data->kind_of_adaptive_grid);
			if (data->iAssignedAdaptiveGridCoefficients == true)
				adaptive_grid.SetConstants(data->adaptive_grid_alfa, data->adaptive_grid_beta, data->adaptive_grid_gamma);
			
			BzzVector xMin(Np); xMin = xOld[1];
			BzzVector xMax(Np); xMax = xOld[xOld.Size()];
			BzzOdeStiff o(xOld, 0., OdeAdaptiveGrid);
			o.SetMinimumConstraints(xMin);
			o.SetMaximumConstraints(xMax);
			
			bzzStop = 0;
			xNew = o(1.e6);
			
			xAdapted = xNew; 
		}

		if (iOption == 1 || iOption == 11 || iOption == 3 || iOption == 33)
		{
			OpenSMOKE_AdaptiveGrid adaptive_grid;
			ptAdaptiveGrid = &adaptive_grid;
			BzzVector xOld(data->iFixedTemperature);
			BzzVector xNew(data->iFixedTemperature);
			BzzVector TOld(data->iFixedTemperature);
			
			for(i=1;i<=data->iFixedTemperature;i++)
			{
				xOld[i] = grid.x[i];
				TOld[i] = T[i];
			}
			
			adaptive_grid.Setup(xOld, TOld);
			adaptive_grid.SetModel(data->kind_of_adaptive_grid);
			if (data->iAssignedAdaptiveGridCoefficients == true)
				adaptive_grid.SetConstants(data->adaptive_grid_alfa, data->adaptive_grid_beta, data->adaptive_grid_gamma);
			
cout << "I am adapting the grid (I)..." << endl;
cout << "Option: " << iOption << endl;
cout << "alfa:   " << data->adaptive_grid_alfa << endl;
cout << "beta:   " << data->adaptive_grid_beta << endl;
cout << "gamma:  " << data->adaptive_grid_gamma << endl;
cout << "Fixed point: " << grid.x[data->iFixedTemperature] << endl;
cout << "Fixed point: " << T[data->iFixedTemperature] << endl;

			BzzVector xMin(xOld.Size()); xMin = xOld[1];
			BzzVector xMax(xOld.Size()); xMax = xOld[xOld.Size()];
			BzzOdeStiff o(xOld, 0., OdeAdaptiveGrid);
			o.SetMinimumConstraints(xMin);
			o.SetMaximumConstraints(xMax);

			bzzStop = 0;
			xNew = o(1.e6);
	//		adaptive_grid.PrintOnFile(xNew);
	//		getchar();
			
			for(i=1;i<=data->iFixedTemperature;i++)
				xAdapted[i] = xNew[i];

cout << "Number of points: " << xNew.Size() << endl;
cout << "xNew: " << xNew.Min() << " " << xNew.Max() << endl;

		}
	
		if (iOption == 2 || iOption == 22 || iOption == 3 || iOption == 33)
		{
			for(int count = 1;count<=1;count++)
			{
				OpenSMOKE_AdaptiveGrid adaptive_grid;
				ptAdaptiveGrid = &adaptive_grid;
				BzzVector xOld(Np+1-data->iFixedTemperature);
				BzzVector xNew(Np+1-data->iFixedTemperature);
				BzzVector TOld(Np+1-data->iFixedTemperature);
			
				for(i=data->iFixedTemperature;i<=Np;i++)
				{
					xOld[i+1-data->iFixedTemperature] = grid.x[i];
					TOld[i+1-data->iFixedTemperature] = T[i];
				}
			
				adaptive_grid.Setup(xOld, TOld);
				adaptive_grid.SetModel(data->kind_of_adaptive_grid);
				if (data->iAssignedAdaptiveGridCoefficients == true)
					adaptive_grid.SetConstants(data->adaptive_grid_alfa, data->adaptive_grid_beta, data->adaptive_grid_gamma);

				cout << "I am adapting the grid (II)..." << endl;
				cout << "Option: " << iOption << endl;
				cout << "alfa:   " << data->adaptive_grid_alfa << endl;
				cout << "beta:   " << data->adaptive_grid_beta << endl;
				cout << "gamma:  " << data->adaptive_grid_gamma << endl;
				cout << "Fixed point: " << grid.x[data->iFixedTemperature] << endl;
				cout << "Fixed point: " << T[data->iFixedTemperature] << endl;
			
				BzzVector xMin(xOld.Size()); xMin = xOld[1];
				BzzVector xMax(xOld.Size()); xMax = xOld[xOld.Size()];
				BzzOdeStiff o(xOld, 0., OdeAdaptiveGrid);
				o.SetMinimumConstraints(xMin);
				o.SetMaximumConstraints(xMax);
			
				bzzStop = 0;
				xNew = o(1.e6);
			//	adaptive_grid.PrintOnFile(xNew);
			//	getchar();
			
				if (fabs(xNew.Max()-xOld.Max())/xOld.Max() < 0.00001 )
				{
					for(i=data->iFixedTemperature;i<=Np;i++)
						xAdapted[i] = xNew[i+1-data->iFixedTemperature];
					
					cout << "Succesfully applied!" << endl;
					cout << "Number of points: " << xNew.Size() << endl;
					cout << "xNew: " << xNew.Min() << " " << xNew.Max() << endl;
					break;
				}
				else
				{
					for(i=data->iFixedTemperature;i<=Np;i++)
						xAdapted[i] = xOld[i+1-data->iFixedTemperature];
					cout << "Unsuccesfully applied!" << endl;
				}
			}
		}

	}
	cout << "I am adapting the solution..." << endl;

	// Linear Interpolation
	grid.AdaptField(U, xAdapted);
	grid.AdaptField(G, xAdapted);
	grid.AdaptField(H, xAdapted);
	grid.AdaptField(T, xAdapted);
	grid.AdaptField(W, xAdapted);
			
	
	if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM )
		grid.AdaptField(moments, xAdapted);
		
	if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )	
	{
		grid.AdaptField(phiN, xAdapted);
		grid.AdaptField(phiM, xAdapted);
	}

	// Construct New grid
	grid.Construct(xAdapted);
			
	// Updating Properties
	BzzVectorInt dummy;
	properties(1, -1, dummy, 0, 0);
	updatingProfiles();

	cout << "Grid correctly adapted" << endl; 
}

void OpenSMOKE_Flame1D::refineFlameBase(const double xA, const double xB)
{
	int i;
	BzzVectorInt listAddPoints;

	for(i=1;i<=Np;i++)
		if (grid.x[i] >= xA/100. && grid.x[i] <= xB/100.) 
			listAddPoints.Append(i);	

	cout << "Base flame refined "<< endl;	

	AddPoints(listAddPoints);
}


void OpenSMOKE_Flame1D::updatingProfiles()
{
	int i;

	// In case of assigned temperature profile is more accurate update the new temperature 
	// profile using the experimental data; therefore:
	if (data->iTemperatureProfile==1 || data->iAssignedFixedTemperatureProfileProvisional==true)
		for(i=1;i<=Np;i++)
			T[i] = data->ud_temperature_profile.GiveMeValue(0., grid.x[i]);

	// In case of assigned cross sectional area profile is more accurate update the new 
	// profile using the experimental data; therefore:
	if (data->iAreaProfile==1)
		for(i=1;i<=Np;i++)
			G[i] = data->ud_cross_section_profile.GiveMeValue(0., grid.x[i]);	// area di passaggio [m2]

	// Fitting formula from  Bittner PhD Thesis (MIT)
	if (data->iAreaProfile==2)
		for(i=1;i<=Np;i++)
			G[i] = data->CrossSectionalArea * (1.+0.114*pow((grid.x[i]*1.e2), 1.44));		// area di passaggio [m2]
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								NLS INTERFACE - Opposed												   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

int OpenSMOKE_Flame1D::solveNLS_Opposed(int Hot_or_Cold, const flame1d_model NLS_KIND)
{
	int i, j, k, dimensionBlock;
	char control;

	// Reaction On/Off
	// ---------------------------------------------------------------------------------------------
	data->iHOT = Hot_or_Cold;		// 0=Cold 1=Hot

	// Kind Of Problem
	// ---------------------------------------------------------------------------------------------
		 if (NLS_KIND == OPPOSED_COLD_REDUCED)			dimensionBlock = data->iX.Size();
	else if (NLS_KIND == OPPOSED_ONLY_MASS_FRACTIONS)	dimensionBlock = NC;
	else if (NLS_KIND == OPPOSED_ONLY_MOMENTUM)			dimensionBlock = 3;
	else if (NLS_KIND == OPPOSED_ONLY_TEMPERATURE)		dimensionBlock = 1;
	else if (NLS_KIND == OPPOSED_ALL)					dimensionBlock = nBlock;
	else if (NLS_KIND == OPPOSED_NO_ENERGY)				dimensionBlock = nBlock;
	
	else ErrorMessage("Error: unknow Reduced NLS!!");

	cout << "-----------------------------------------------------------------------" << endl;
	cout << "NLS - " << NLS_KIND << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(Np*dimensionBlock, &fNLS);
	ChangeDimensions(Np*dimensionBlock, &xFirstGuess);
	ChangeDimensions(Np*dimensionBlock, &xMin);
	ChangeDimensions(Np*dimensionBlock, &xMax);

	setMinimumAndMaximumValues(NLS_KIND);	

	if (NLS_KIND == OPPOSED_COLD_REDUCED)
	{
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=data->iX.Size();j++)
				xFirstGuess[k++] = W[i][data->iX[j]];

		tagConstantT = OffConstantT;

		NLS_Opposed_COLDREDUCED.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Opposed_COLDREDUCED, dimensionBlock);		
		control = nonLinearSystemSolution(o, dimensionBlock);
		tagConstantT = ResetConstantT;
	}


	if (NLS_KIND == OPPOSED_ONLY_MASS_FRACTIONS)
	{
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];

		tagConstantT = OffConstantT;

		NLS_Opposed_ONLYMASSFRACTIONS.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Opposed_ONLYMASSFRACTIONS, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);
		
		tagConstantT = ResetConstantT;
	}

	if (NLS_KIND == OPPOSED_ONLY_MOMENTUM)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
		}

		tagConstantT = OffConstantT;
 
		NLS_Opposed_ONLYUGH.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Opposed_ONLYUGH, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);

		tagConstantT = ResetConstantT;	
  }

	if (NLS_KIND == OPPOSED_ONLY_TEMPERATURE)
	{
		k=1;
		for(i=1;i<=Np;i++)
			xFirstGuess[k++] = T[i];
	
		NLS_Opposed_ONLYT.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Opposed_ONLYT, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);		
	}

	if (NLS_KIND == OPPOSED_ALL)
	{
		ErrorMessage("OPPOSED_ALL is not longer supported");
		/*
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		NLS_Opposed_ALL.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Opposed_ALL, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);	
		
		// Sensitivity Analysis
		final_sensitivity_analysis("OPPOSED", o);*/

	/*	{	
			BzzVector xSolution(xFirstGuess.Size());
			BzzVector fSolution(xFirstGuess.Size());
			o.GetSolution(&xSolution, &fSolution);
			final_sensitivity_analysis_diffusion("OPPOSED", o, xSolution, fSolution);
		}*/
	}

	if (NLS_KIND == OPPOSED_NO_ENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		tagConstantT = OffConstantT;

		NLS_Opposed_NOENERGY.assignFlame(this);

		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Opposed_NOENERGY, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);
		
		tagConstantT = ResetConstantT;
	}

	return int(control);
}

int OpenSMOKE_Flame1D::solveSensitivity(const flame1d_model NLS_KIND, BzzSave &fBinary)
{
	int dimensionBlock;
	char control;

	// Reaction On/Off
	data->iHOT = 1;		// 0=Cold 1=Hot

	// Kind Of Problem
		 if (NLS_KIND == OPPOSED_ALL)			dimensionBlock = nBlock;
	else if (NLS_KIND == TWIN_ALL)				dimensionBlock = nBlock;
	else if (NLS_KIND == PREMIXED_FLAMESPEED)	dimensionBlock = NC+2;
	else if (NLS_KIND == PREMIXED_ALL)			dimensionBlock = NC+1;
	else if (NLS_KIND == PREMIXED_NOENERGY)		dimensionBlock = NC;
	else ErrorMessage("Wrong sensitivity analysis choice");

	cout << "-----------------------------------------------------------------------" << endl;
	cout << " Sensitivity - Npoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(Np*dimensionBlock, &fNLS);
	ChangeDimensions(Np*dimensionBlock, &xFirstGuess);
	ChangeDimensions(Np*dimensionBlock, &xMin);
	ChangeDimensions(Np*dimensionBlock, &xMax);

	setMinimumAndMaximumValues(NLS_KIND);	

	BzzNonLinearSystemSparseObject o;

	if (NLS_KIND == OPPOSED_ALL)
	{
		int k=1;
		for(int i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		NLS_Opposed_ALL.assignFlame(this);
		o(xFirstGuess, &NLS_Opposed_ALL, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);	
	}

	if (NLS_KIND == TWIN_ALL)
	{
		int k=1;
		for(int i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		NLS_Twin_ALL.assignFlame(this);
		o(xFirstGuess, &NLS_Twin_ALL, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);			
	}
	
	if (NLS_KIND == PREMIXED_FLAMESPEED)
	{
		int k=1;
		for(int i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			xFirstGuess[k++] = H[i];
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		xFirstGuess( dimensionBlock*(data->iFixedTemperature-1) + 1) = data->fixedTemperature;

		NLS_Premixed_FLAMESPEED.assignFlame(this);
		o(xFirstGuess, &NLS_Premixed_FLAMESPEED, dimensionBlock);
		nonLinearSystemSolution(o, dimensionBlock);
	}

	if (NLS_KIND == PREMIXED_ALL)
	{
		int k=1;
		for(int i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}
		
		NLS_Premixed_ALL.assignFlame(this);
		o(xFirstGuess, &NLS_Premixed_ALL, dimensionBlock);
		nonLinearSystemSolution(o, dimensionBlock);
	}

	if (NLS_KIND == PREMIXED_NOENERGY)
	{
		int k=1;
		for(int i=1;i<=Np;i++)
		{
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}
		
		NLS_Premixed_NOENERGY.assignFlame(this);
		o(xFirstGuess, &NLS_Premixed_NOENERGY, dimensionBlock);
		nonLinearSystemSolution(o, dimensionBlock);
	}

	// Sensitivity Analysis
	if (kindOfSensitivity == FREQUENCY_FACTOR)
	{
		if (NLS_KIND == OPPOSED_ALL)			final_sensitivity_analysis("OPPOSED_ALL", o, fBinary);
		if (NLS_KIND == TWIN_ALL)				final_sensitivity_analysis("TWIN_ALL", o, fBinary);
		if (NLS_KIND == PREMIXED_FLAMESPEED)	final_sensitivity_analysis("PREMIXED_FLAMESPEED", o, fBinary);
		if (NLS_KIND == PREMIXED_ALL)			final_sensitivity_analysis("PREMIXED_ALL", o, fBinary);
		if (NLS_KIND == PREMIXED_NOENERGY)		final_sensitivity_analysis("PREMIXED_NOENERGY", o, fBinary);
	}
	
	else if (kindOfSensitivity == TRANSPORT_PROPERTIES)
	{
		BzzVector xSolution(xFirstGuess.Size());
		BzzVector fSolution(xFirstGuess.Size());
		o.GetSolution(&xSolution, &fSolution);

		if (NLS_KIND == OPPOSED_ALL)			final_sensitivity_analysis_transport_properties("OPPOSED_ALL", o, xSolution, fSolution, fBinary);
		if (NLS_KIND == TWIN_ALL)				final_sensitivity_analysis_transport_properties("TWIN_ALL", o, xSolution, fSolution, fBinary);
		if (NLS_KIND == PREMIXED_FLAMESPEED)	final_sensitivity_analysis_transport_properties("PREMIXED_FLAMESPEED", o, xSolution, fSolution, fBinary);
		if (NLS_KIND == PREMIXED_ALL)			final_sensitivity_analysis_transport_properties("PREMIXED_ALL", o, xSolution, fSolution, fBinary);
		if (NLS_KIND == PREMIXED_NOENERGY)		final_sensitivity_analysis_transport_properties("PREMIXED_NOENERGY", o, xSolution, fSolution, fBinary);
	}

	else if (kindOfSensitivity == FREQUENCY_FACTOR_AND_TRANSPORT_PROPERTIES)
	{
		BzzVector xSolution(xFirstGuess.Size());
		BzzVector fSolution(xFirstGuess.Size());
		o.GetSolution(&xSolution, &fSolution);

		if (NLS_KIND == OPPOSED_ALL)			final_sensitivity_analysis("OPPOSED_ALL", o, fBinary);
		if (NLS_KIND == TWIN_ALL)				final_sensitivity_analysis("TWIN_ALL", o, fBinary);
		if (NLS_KIND == PREMIXED_FLAMESPEED)	final_sensitivity_analysis("PREMIXED_FLAMESPEED", o, fBinary);
		if (NLS_KIND == PREMIXED_ALL)			final_sensitivity_analysis("PREMIXED_ALL", o, fBinary);
		if (NLS_KIND == PREMIXED_NOENERGY)		final_sensitivity_analysis("PREMIXED_NOENERGY", o, fBinary);


		if (NLS_KIND == OPPOSED_ALL)			final_sensitivity_analysis_transport_properties("OPPOSED_ALL", o, xSolution, fSolution, fBinary);
		if (NLS_KIND == TWIN_ALL)				final_sensitivity_analysis_transport_properties("TWIN_ALL", o, xSolution, fSolution, fBinary);
		if (NLS_KIND == PREMIXED_FLAMESPEED)	final_sensitivity_analysis_transport_properties("PREMIXED_FLAMESPEED", o, xSolution, fSolution, fBinary);
		if (NLS_KIND == PREMIXED_ALL)			final_sensitivity_analysis_transport_properties("PREMIXED_ALL", o, xSolution, fSolution, fBinary);
		if (NLS_KIND == PREMIXED_NOENERGY)		final_sensitivity_analysis_transport_properties("PREMIXED_NOENERGY", o, xSolution, fSolution, fBinary);
	}

	return int(control);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								NLS INTERFACE - Twin												   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

int OpenSMOKE_Flame1D::solveNLS_Twin(int Hot_or_Cold, const flame1d_model NLS_KIND)
{
	int i, j, k, dimensionBlock;
	char control;

	// Reaction On/Off
	// ---------------------------------------------------------------------------------------------
	data->iHOT = Hot_or_Cold;		// 0=Cold 1=Hot

	// Kind Of Problem
	// ---------------------------------------------------------------------------------------------
		 if (NLS_KIND == TWIN_COLD_REDUCED)			dimensionBlock = data->iX.Size();
	else if (NLS_KIND == TWIN_ONLY_MASS_FRACTIONS)	dimensionBlock = NC;
	else if (NLS_KIND == TWIN_ONLY_MOMENTUM)		dimensionBlock = 3;
	else if (NLS_KIND == TWIN_ONLY_TEMPERATURE)		dimensionBlock = 1;
	else if (NLS_KIND == TWIN_ALL)					dimensionBlock = nBlock;
	else if (NLS_KIND == TWIN_NO_ENERGY)			dimensionBlock = nBlock;
	
	else ErrorMessage("Error: unknow Reduced NLS!!");

	cout << "-----------------------------------------------------------------------" << endl;
	cout << "NLS - " << NLS_KIND << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(Np*dimensionBlock, &fNLS);
	ChangeDimensions(Np*dimensionBlock, &xFirstGuess);
	ChangeDimensions(Np*dimensionBlock, &xMin);
	ChangeDimensions(Np*dimensionBlock, &xMax);

	setMinimumAndMaximumValues(NLS_KIND);	

	if (NLS_KIND == TWIN_COLD_REDUCED)
	{
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=data->iX.Size();j++)
				xFirstGuess[k++] = W[i][data->iX[j]];

		tagConstantT = OffConstantT;

		NLS_Twin_COLDREDUCED.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Twin_COLDREDUCED, dimensionBlock);		
		control = nonLinearSystemSolution(o, dimensionBlock);
		
		tagConstantT = ResetConstantT;
	}


	if (NLS_KIND == TWIN_ONLY_MASS_FRACTIONS)
	{
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];

		tagConstantT = OffConstantT;

		NLS_Twin_ONLYMASSFRACTIONS.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Twin_ONLYMASSFRACTIONS, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);
		
		tagConstantT = ResetConstantT;
	}

	if (NLS_KIND == TWIN_ONLY_MOMENTUM)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
		}

		tagConstantT = OffConstantT;
 
		NLS_Twin_ONLYUGH.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Twin_ONLYUGH, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);

		tagConstantT = ResetConstantT;	
  }

	if (NLS_KIND == TWIN_ONLY_TEMPERATURE)
	{
		k=1;
		for(i=1;i<=Np;i++)
			xFirstGuess[k++] = T[i];
	
		NLS_Twin_ONLYT.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Twin_ONLYT, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);		
	}

	if (NLS_KIND == TWIN_ALL)
	{
		ErrorMessage("TWIN ALL is not longer supported");
		/*
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		NLS_Twin_ALL.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Twin_ALL, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);			

		// Sensitivity Analysis
		final_sensitivity_analysis("TWIN", o);*/
	}

	if (NLS_KIND == TWIN_NO_ENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		tagConstantT = OffConstantT;

		NLS_Twin_NOENERGY.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Twin_NOENERGY, dimensionBlock);
		control = nonLinearSystemSolution(o, dimensionBlock);
		
		tagConstantT = ResetConstantT;
	}

	return int(control);
}

int OpenSMOKE_Flame1D::solveNLS_Reduced_Opposed(int Hot_or_Cold, const flame1d_model NLS_KIND)
{
	// Reaction On/Off
	// ---------------------------------------------------------------------------------------------
	data->iHOT = Hot_or_Cold;		// 0=Cold 1=Hot

	// Kind Of Problem
	// ---------------------------------------------------------------------------------------------
	int number_of_equations = 0;
	if (NLS_KIND == OPPOSED_ALL)					number_of_equations = 2+2*NC+3*Np;
	else ErrorMessage("Error: unknow Reduced NLS!!");

	cout << "-----------------------------------------------------------------------" << endl;
	cout << "NLS Reduced - " << NLS_KIND << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(number_of_equations, &fNLS);
	ChangeDimensions(number_of_equations, &xFirstGuess);
	ChangeDimensions(number_of_equations, &xMin);
	ChangeDimensions(number_of_equations, &xMax);

	setMinimumAndMaximumValuesReduced(NLS_KIND);	

	char control;
	if (NLS_KIND == OPPOSED_ALL)
	{
		int k=1;
		{
			xFirstGuess[k++] = U[1];
			xFirstGuess[k++] = G[1];
			xFirstGuess[k++] = H[1];
			xFirstGuess[k++] = T[1];
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[1][j];
		}

		for(int i=2;i<=Ni;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
		}

		{
			xFirstGuess[k++] = U[Np];
			xFirstGuess[k++] = G[Np];
			xFirstGuess[k++] = H[Np];
			xFirstGuess[k++] = T[Np];
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[Np][j];
		}

		NLS_Reduced_Opposed_ALL.assignFlame(this);
		BzzNonLinearSystemObject o(xFirstGuess, &NLS_Reduced_Opposed_ALL);
		control = nonLinearSystemSolution(o);			
	}

	return int(control);
}

int OpenSMOKE_Flame1D::solveNLS_Reduced_Twin(int Hot_or_Cold, const flame1d_model NLS_KIND)
{
	// Reaction On/Off
	// ---------------------------------------------------------------------------------------------
	data->iHOT = Hot_or_Cold;		// 0=Cold 1=Hot

	// Kind Of Problem
	// ---------------------------------------------------------------------------------------------
	int number_of_equations = 0;
	if (NLS_KIND == TWIN_ALL)					number_of_equations = 2+2*NC+3*Np;
	else ErrorMessage("Error: unknow Reduced NLS!!");

	cout << "-----------------------------------------------------------------------" << endl;
	cout << "NLS Reduced - " << NLS_KIND << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(number_of_equations, &fNLS);
	ChangeDimensions(number_of_equations, &xFirstGuess);
	ChangeDimensions(number_of_equations, &xMin);
	ChangeDimensions(number_of_equations, &xMax);

	setMinimumAndMaximumValuesReduced(NLS_KIND);	

	char control;
	if (NLS_KIND == TWIN_ALL)
	{
		int k=1;
		{
			xFirstGuess[k++] = U[1];
			xFirstGuess[k++] = G[1];
			xFirstGuess[k++] = H[1];
			xFirstGuess[k++] = T[1];
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[1][j];
		}

		for(int i=2;i<=Ni;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
		}

		{
			xFirstGuess[k++] = U[Np];
			xFirstGuess[k++] = G[Np];
			xFirstGuess[k++] = H[Np];
			xFirstGuess[k++] = T[Np];
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = W[Np][j];
		}

		NLS_Reduced_Twin_ALL.assignFlame(this);
		BzzNonLinearSystemObject o(xFirstGuess, &NLS_Reduced_Twin_ALL);
		control = nonLinearSystemSolution(o);			
	}

	return int(control);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								NLS INTERFACE - Premixed												   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::solveNLS_Premixed(const flame1d_model string_kind)
{
	int i, j, k, dimensionBlock;

	if (string_kind == PREMIXED_ALL)				{ data->iHOT=1; dimensionBlock = NC+1;}
	else if (string_kind == PREMIXED_NOENERGY)	{ data->iHOT=1; dimensionBlock = NC;}
	else if (string_kind == PREMIXED_FLAMESPEED)	{ data->iHOT=1; dimensionBlock = NC+2;}	

	else ErrorMessage("Error in SolveNLS, Wrong Type:" + string_kind);

	cout << "-----------------------------------------------------------------------" << endl;
	cout << "NLS - " << string_kind << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(Np*dimensionBlock, &fNLS);
	ChangeDimensions(Np*dimensionBlock, &xFirstGuess);
	ChangeDimensions(Np*dimensionBlock, &xMin);
	ChangeDimensions(Np*dimensionBlock, &xMax);

	setMinimumAndMaximumValues(string_kind);

	if (string_kind == PREMIXED_ALL)
	{
		ErrorMessage("PREMIXED_ALL is not longer supported");
/*
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}
		
		NLS_Premixed_ALL.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Premixed_ALL, dimensionBlock);
		nonLinearSystemSolution(o, dimensionBlock);

		// Sensitivity Analysis
		final_sensitivity_analysis("PREMIXED", o);*/
	}

	else if (string_kind == PREMIXED_NOENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];

		tagConstantT = OffConstantT;

		NLS_Premixed_NOENERGY.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Premixed_NOENERGY, dimensionBlock);
		nonLinearSystemSolution(o, dimensionBlock);		

		tagConstantT = ResetConstantT;
	}

	if (string_kind == PREMIXED_FLAMESPEED)
	{
		ErrorMessage("PREMIXED_FLAMESPEED is not longer supported");
	/*
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			xFirstGuess[k++] = H[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		xFirstGuess( dimensionBlock*(data->iFixedTemperature-1) + 1) = data->fixedTemperature;

		NLS_Premixed_FLAMESPEED.assignFlame(this);
		BzzNonLinearSystemSparseObject o(xFirstGuess, &NLS_Premixed_FLAMESPEED, dimensionBlock);
		nonLinearSystemSolution(o, dimensionBlock);
				
		// Sensitivity Analysis
		final_sensitivity_analysis("PREMIXED_FLAMESPEED", o);*/
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DAE INTERFACE - Opposed										   //
//																								   //
////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::solveDAE_Opposed(int Hot_Or_Cold, const flame1d_model string_kind, double tEnd)
{
	int i, j, k;
	int dimensionBlock;
	//data->iUnsteady = true;		// TODO

	data->iHOT = Hot_Or_Cold;

		 if (string_kind == OPPOSED_ALL)				dimensionBlock=nBlock;
	else if (string_kind == OPPOSED_NO_MOMENTUM)		dimensionBlock=nBlock-3;
	else if (string_kind == OPPOSED_NO_ENERGY)			dimensionBlock=nBlock-1;
	else if (string_kind == OPPOSED_ONLY_MOMENTUM)		dimensionBlock=3;
	else if (string_kind == OPPOSED_SOOT_ALL)			dimensionBlock=nBlock+2;
	else if (string_kind == OPPOSED_QMOM_ALL)			dimensionBlock=nBlock+2*qmom.N;
	
	else ErrorMessage("Error in SolveDAE, Wrong Type" + string_kind);

	cout << "-----------------------------------------------------------------------" << endl;
	cout << "DAE solution: " << string_kind << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(dimensionBlock*Np, &xFirstGuess);
	ChangeDimensions(dimensionBlock*Np, &xMin);
	ChangeDimensions(dimensionBlock*Np, &xMax);
	ChangeDimensions(dimensionBlock*Np, &inDerAlg);
	
	std::string nameUnsteadyFile;

	if (data->iUnsteady == true)
	{
		nameUnsteadyFile = nameFolderUnsteadyData + "/Unsteady.out";
		openOutputFileAndControl(fUnsteady, nameUnsteadyFile);
		fUnsteady.setf(ios::scientific);
		
		nameUnsteadyFile = nameFolderUnsteadyData + "/Unsteady_Max.out";
		openOutputFileAndControl(fUnsteadyMax, nameUnsteadyFile);
		fUnsteadyMax.setf(ios::scientific);

		nameUnsteadyFile = nameFolderUnsteadyData + "/Unsteady_QMOMMax.out";
		openOutputFileAndControl(fUnsteadyQMOMMax, nameUnsteadyFile);
		fUnsteadyQMOMMax.setf(ios::scientific);
	}

	setDifferentialAndAlgebraic(string_kind);
	setMinimumAndMaximumValues(string_kind);
	
	if (string_kind == OPPOSED_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		if (iUnsteadyFromBackUp!=-1 && data->iUnsteady == true)
			unsteady.setup(	U[1], U[Np], rho[1], rho[Np], PMtot[1], PMtot[Np], 
							T[1], T[Np], grid.L, data->P_Pascal, mix, data->unsteady_flame_file_name);
		GnuPlotInterfaceUnsteady();


		OpenSMOKE_Flame1D_Solution previous_solution;
		previous_solution.PasteFromExternalSolution(*this, *data);
		
		bzzStop = 0;
		DAE_Opposed_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Opposed_ALL, dimensionBlock);	
		DAESystemSolution(&o, tEnd);  
		bzzStop = 0;

		bool iSuccess;
		if (o.GetHUsed() <= tEnd/1000.)	iSuccess = false;
		else														iSuccess = true;

		if (iSuccess == false)
		{
		}

	}

	if (string_kind == OPPOSED_NO_MOMENTUM)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		DAE_Opposed_NOMOMENTUM.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Opposed_NOMOMENTUM, dimensionBlock);
		DAESystemSolution(&o, tEnd); 
	}

	if (string_kind == OPPOSED_ONLY_MOMENTUM)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
		}

		DAE_Opposed_ONLYMOMENTUM.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Opposed_ONLYMOMENTUM, dimensionBlock);
		DAESystemSolution(&o, tEnd); 
	}

	if (string_kind == OPPOSED_NO_ENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		tagConstantT = OffConstantT;

		DAE_Opposed_NOENERGY.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Opposed_NOENERGY, dimensionBlock);		
		DAESystemSolution(&o, tEnd); 

		tagConstantT = ResetConstantT;
	}

	if (string_kind == OPPOSED_SOOT_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
			xFirstGuess[k++] = phiN[i];
			xFirstGuess[k++] = phiM[i];
		}

		if (iUnsteadyFromBackUp!=-1 && data->iUnsteady == true)
			unsteady.setup(	U[1], U[Np], rho[1], rho[Np], PMtot[1], PMtot[Np], 
							T[1], T[Np], grid.L, data->P_Pascal, mix, data->unsteady_flame_file_name);

		GnuPlotInterfaceUnsteady();

		DAE_Opposed_SOOT_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Opposed_SOOT_ALL, dimensionBlock);
		DAESystemSolution(&o, tEnd); 
	}

	if (string_kind == OPPOSED_QMOM_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
			for(j=1;j<=2*qmom.N;j++)
				xFirstGuess[k++] = moments[i][j];
		}

		if (iUnsteadyFromBackUp!=-1 && data->iUnsteady == true)
			unsteady.setup(	U[1], U[Np], rho[1], rho[Np], PMtot[1], PMtot[Np], 
							T[1], T[Np], grid.L, data->P_Pascal, mix, data->unsteady_flame_file_name);

		GnuPlotInterfaceUnsteady();

		DAE_Opposed_QMOM_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Opposed_QMOM_ALL, dimensionBlock);		
		DAESystemSolution(&o, tEnd); 
	}

	fUnsteady.close();
	fUnsteadyMax.close();
	fUnsteadyQMOMMax.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DAE INTERFACE - Twin									   //
//																								   //
////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::solveDAE_Twin(int Hot_Or_Cold, const flame1d_model string_kind, double tEnd)
{
	int i, j, k;
	int dimensionBlock;
	//data->iUnsteady = true;		// TODO

	data->iHOT = Hot_Or_Cold;

		 if (string_kind == TWIN_ALL)										dimensionBlock=nBlock;
	else if (string_kind == TWIN_NO_MOMENTUM)		dimensionBlock=nBlock-3;
	else if (string_kind == TWIN_NO_ENERGY)					dimensionBlock=nBlock-1;
	else if (string_kind == TWIN_SOOT_ALL)						dimensionBlock=nBlock+2;
	else if (string_kind == TWIN_QMOM_ALL)					dimensionBlock=nBlock+2*qmom.N;
	
	else ErrorMessage("Error in SolveDAE, Wrong Type" + string_kind);

	cout << "-----------------------------------------------------------------------" << endl;
	cout << "DAE solution: " << string_kind << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(dimensionBlock*Np, &xFirstGuess);
	ChangeDimensions(dimensionBlock*Np, &xMin);
	ChangeDimensions(dimensionBlock*Np, &xMax);
	ChangeDimensions(dimensionBlock*Np, &inDerAlg);
	
	std::string nameUnsteadyFile;

	if (data->iUnsteady == true)
	{
		nameUnsteadyFile = nameFolderUnsteadyData + "/Unsteady.out";
		openOutputFileAndControl(fUnsteady, nameUnsteadyFile);
		fUnsteady.setf(ios::scientific);
		
		nameUnsteadyFile = nameFolderUnsteadyData + "/Unsteady_Max.out";
		openOutputFileAndControl(fUnsteadyMax, nameUnsteadyFile);
		fUnsteadyMax.setf(ios::scientific);

		nameUnsteadyFile = nameFolderUnsteadyData + "/Unsteady_QMOMMax.out";
		openOutputFileAndControl(fUnsteadyQMOMMax, nameUnsteadyFile);
		fUnsteadyQMOMMax.setf(ios::scientific);
	}

	setDifferentialAndAlgebraic(string_kind);
	setMinimumAndMaximumValues(string_kind);
	
	if (string_kind == TWIN_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		if (iUnsteadyFromBackUp!=-1 && data->iUnsteady == true)
			unsteady.setup(	U[1], U[Np], rho[1], rho[Np], PMtot[1], PMtot[Np], 
							T[1], T[Np], grid.L, data->P_Pascal, mix, data->unsteady_flame_file_name);
		GnuPlotInterfaceUnsteady();


		OpenSMOKE_Flame1D_Solution previous_solution;
		previous_solution.PasteFromExternalSolution(*this, *data);
		
		DAE_Twin_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Twin_ALL, dimensionBlock);	
		DAESystemSolution(&o, tEnd);  
		
		bool iSuccess;
		if (o.GetHUsed() <= tEnd/1000.)	iSuccess = false;
		else														iSuccess = true;

		if (iSuccess == false)
		{
		}

	}

	if (string_kind == TWIN_NO_MOMENTUM)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		DAE_Twin_NOMOMENTUM.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Twin_NOMOMENTUM, dimensionBlock);
		DAESystemSolution(&o, tEnd); 
	}

	if (string_kind == TWIN_NO_ENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		tagConstantT = OffConstantT;

		DAE_Twin_NOENERGY.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Twin_NOENERGY, dimensionBlock);		
		DAESystemSolution(&o, tEnd); 

		tagConstantT = ResetConstantT;
	}

	if (string_kind == TWIN_SOOT_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
			xFirstGuess[k++] = phiN[i];
			xFirstGuess[k++] = phiM[i];
		}

		if (iUnsteadyFromBackUp!=-1 && data->iUnsteady == true)
			unsteady.setup(	U[1], U[Np], rho[1], rho[Np], PMtot[1], PMtot[Np], 
							T[1], T[Np], grid.L, data->P_Pascal, mix, data->unsteady_flame_file_name);

		GnuPlotInterfaceUnsteady();

		DAE_Twin_SOOT_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Twin_SOOT_ALL, dimensionBlock);
		DAESystemSolution(&o, tEnd); 
	}

	if (string_kind == TWIN_QMOM_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = U[i];
			xFirstGuess[k++] = G[i];
			xFirstGuess[k++] = H[i];
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
			for(j=1;j<=2*qmom.N;j++)
				xFirstGuess[k++] = moments[i][j];
		}

		if (iUnsteadyFromBackUp!=-1 && data->iUnsteady == true)
			unsteady.setup(	U[1], U[Np], rho[1], rho[Np], PMtot[1], PMtot[Np], 
							T[1], T[Np], grid.L, data->P_Pascal, mix, data->unsteady_flame_file_name);

		GnuPlotInterfaceUnsteady();

		DAE_Twin_QMOM_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Twin_QMOM_ALL, dimensionBlock);		
		DAESystemSolution(&o, tEnd); 
	}

	fUnsteady.close();
	fUnsteadyMax.close();
	fUnsteadyQMOMMax.close();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DAE INTERFACE - Premixed										   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::solveDAE_Premixed(const flame1d_model string_kind, double tEnd)
{
	int i,j,k;
	int dimensionBlock;

	data->iHOT = 1;

		 if (string_kind == PREMIXED_ALL)			dimensionBlock=NC+1;	
	else if (string_kind == PREMIXED_FLAMESPEED)	dimensionBlock=NC+2;	
	else if (string_kind == PREMIXED_NOENERGY)		dimensionBlock=NC;	
	else if (string_kind == PREMIXED_QMOM_ALL)		dimensionBlock=NC+1+2*qmom.N;	
	else if (string_kind == PREMIXED_QMOM_NOENERGY)	dimensionBlock=NC  +2*qmom.N;	
	else if (string_kind == PREMIXED_SOOT_ALL)		dimensionBlock=NC+1+2;	
	else if (string_kind == PREMIXED_SOOT_NOENERGY)	dimensionBlock=NC  +2;	
	else ErrorMessage("Error in SolveDAE, Wrong Type " + string_kind);
	
	cout << "-----------------------------------------------------------------------" << endl;
	cout << "DAE solution: " << string_kind << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(dimensionBlock*Np, &xFirstGuess);
	ChangeDimensions(dimensionBlock*Np, &xMin);
	ChangeDimensions(dimensionBlock*Np, &xMax);
	ChangeDimensions(dimensionBlock*Np, &inDerAlg);
	
	if (data->iUnsteady == true)
	{
		std::string nameUnsteadyFile;
		nameUnsteadyFile = nameFolderUnsteadyData + "/Unsteady.out";
		openOutputFileAndControl(fUnsteady, nameUnsteadyFile);
		fUnsteady.setf(ios::scientific);
	}

	setDifferentialAndAlgebraic(string_kind);
	setMinimumAndMaximumValues(string_kind);
	
	if (string_kind == PREMIXED_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		DAE_Premixed_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Premixed_ALL, dimensionBlock);
		DAESystemSolution(&o, tEnd);
	}

	else if (string_kind == PREMIXED_NOENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		tagConstantT = OffConstantT;

		DAE_Premixed_NOENERGY.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Premixed_NOENERGY, dimensionBlock);		
		DAESystemSolution(&o, tEnd);

		tagConstantT = ResetConstantT;
	}

	if (string_kind == PREMIXED_FLAMESPEED)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			xFirstGuess[k++] = H[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		xFirstGuess( dimensionBlock*(data->iFixedTemperature-1) + 1) = data->fixedTemperature;
				
		DAE_Premixed_FLAMESPEED.assignFlame(this);
		BzzDaeSparseObject o;
		o(xFirstGuess, 0., inDerAlg, &DAE_Premixed_FLAMESPEED, dimensionBlock);
		DAESystemSolution(&o, tEnd);
	}

	if (string_kind == PREMIXED_QMOM_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
			for(j=1;j<=2*qmom.N;j++)
				xFirstGuess[k++] = moments[i][j];
		}

		DAE_Premixed_QMOM_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Premixed_QMOM_ALL, dimensionBlock);
		DAESystemSolution(&o, tEnd);
	}

	if (string_kind == PREMIXED_QMOM_NOENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
			for(j=1;j<=2*qmom.N;j++)
				xFirstGuess[k++] = moments[i][j];
		}

		tagConstantT = OffConstantT;

		DAE_Premixed_QMOM_NOENERGY.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Premixed_QMOM_NOENERGY, dimensionBlock);
		DAESystemSolution(&o, tEnd);
	}

	if (string_kind == PREMIXED_SOOT_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
			xFirstGuess[k++] = phiN[i];
			xFirstGuess[k++] = phiM[i];
		}

		DAE_Premixed_SOOT_ALL.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Premixed_SOOT_ALL, dimensionBlock);
		DAESystemSolution(&o, tEnd);
	}

	if (string_kind == PREMIXED_SOOT_NOENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
			xFirstGuess[k++] = phiN[i];
			xFirstGuess[k++] = phiM[i];
		}

		tagConstantT = OffConstantT;

		DAE_Premixed_SOOT_NOENERGY.assignFlame(this);
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &DAE_Premixed_SOOT_NOENERGY, dimensionBlock);
		DAESystemSolution(&o, tEnd);
	}

	fUnsteady.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE INTERFACE - Premixed										   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::solveODE_Premixed(const flame1d_model string_kind, double tEnd)
{
	int i,j,k;
	int dimensionBlock;

	data->iHOT = 1;

		 if (string_kind == PREMIXED_ALL)				dimensionBlock=NC+1;	
	else if (string_kind == PREMIXED_NOENERGY)		dimensionBlock=NC;	
	else ErrorMessage("Error in SolveODE, Wrong Type: " + string_kind);

	cout << "-----------------------------------------------------------------------" << endl;
	cout << "ODE solution: " << string_kind << "\tNpoints: " << Np << endl;
	cout << "-----------------------------------------------------------------------" << endl;

	ChangeDimensions(dimensionBlock*Np, &xFirstGuess);
	ChangeDimensions(dimensionBlock*Np, &xMin);
	ChangeDimensions(dimensionBlock*Np, &xMax);
	
	if (data->iUnsteady == true)
	{
		std::string nameUnsteadyFile;
		nameUnsteadyFile = nameFolderUnsteadyData + "/Unsteady.out";
		openOutputFileAndControl(fUnsteady, nameUnsteadyFile);
		fUnsteady.setf(ios::scientific);
	}

	setMinimumAndMaximumValues(string_kind);
	
	if (string_kind == PREMIXED_ALL)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			xFirstGuess[k++] = T[i];
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		ODE_Premixed_ALL.assignFlame(this);
		BzzOdeSparseStiffObject o(xFirstGuess, 0., &ODE_Premixed_ALL, dimensionBlock);
		ODESystemSolution(o, tEnd);
	}

	else if (string_kind == PREMIXED_NOENERGY)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			for(j=1;j<=NC;j++)
				xFirstGuess[k++] = W[i][j];
		}

		tagConstantT = OffConstantT;

		ODE_Premixed_NOENERGY.assignFlame(this);
		BzzOdeSparseStiffObject o(xFirstGuess, 0., &ODE_Premixed_NOENERGY, dimensionBlock);
		ODESystemSolution(o, tEnd);

		tagConstantT = ResetConstantT;
	}

	fUnsteady.close();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								GIVE RESIDUALS - Opposed											//
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::give_Opposed_DU_DG_DH()
{
	int i;

	for(i=1;i<=Np;i++)
		M[i] = U[i]*G[i]*urho[i];

	grid.FirstDerivative(data->iDerG, U, M, diffM);

	// Equations
	// --------------------------------------------------------------------------------------------
	dU[1] = U[1] - UC;
	dG[1] = G[1] - GC;
	dH[1] = H[1] - H[2];

	if (data->iPoolFire == POOL_FIRE_TASSIGNED || data->iPoolFire == POOL_FIRE_EQUILIBRIUM)
	{
		double DHvap = data->pool_fire_liquid_species->Hv(T[1]) * data->correctionFactorVaporizationHeat;
		dU[1] = (nGeometry-1.)*U[1] * DHvap - lambda[1] * (T[2] - T[1])*grid.udxe[1];
		data->VC = (nGeometry-1.)*U[1] / rho[1];
	}

	if (data->iPoolFire == POOL_FIRE_LIQUIDPOOL)
	{
		double lambdaLiquid = data->pool_fire_liquid_species->lambda(0.50*(T[1]+data->pool_fire_feed_temperature));
		double DHvap = data->pool_fire_liquid_species->Hv(T[1]) * data->correctionFactorVaporizationHeat;
		dU[1] = (nGeometry-1.)*U[1] * DHvap - lambda[1] * (T[2] - T[1])*grid.udxe[1] + lambdaLiquid*(T[1] - data->pool_fire_feed_temperature) / data->pool_fire_depth;
		data->VC = (nGeometry-1.)*U[1] / rho[1];
	}
	
	for(i=2;i<=Ni;i++)
	{
		dU[i] = -0.50*(G[i-1]+G[i]) + (-U[i-1]+U[i])*grid.udxw[i];

		dG[i] = - (	+ (nGeometry-1.)*diffM[i] + 
					- nGeometry*G[i]*G[i]*urho[i] +
					- (  mu[i]*(G[i+1]*urho[i+1]-G[i]*urho[i])*grid.udxe[i] -
						 mu[i-1]*(G[i]*urho[i]-G[i-1]*urho[i-1])*grid.udxw[i])*grid.udxc_over_2[i] +
					- H[i] );

		dH[i] = H[i] - H[i+1]; 
	}

	// Soot deposition (start)
	if (data->iDepositionWall == true)
	{
		// Units: kg/m2/s
		soot_deposition = 0.;
		for(int j=1;j<=mix->polimiSoot->bin_indices().Size();j++)
		{
			const int jj = mix->polimiSoot->bin_indices()[j];
			soot_deposition += vThermophoretic[Np][jj];
		}
		soot_deposition *= rho[Np];
	}
	// Soot deposition (end)

	dU[Np] = -0.50*(G[Np-1]+G[Np]) + (-U[Np-1]+U[Np])*grid.udxw[Np];
	dG[Np] = G[Np] - GO;
	dH[Np] = U[Np] - UO;
}	

void OpenSMOKE_Flame1D::give_Twin_DU_DG_DH()
{
	give_Opposed_DU_DG_DH();

	// Original set boundary conditions (twin flames)
	dU[Np] = - (+ (nGeometry-1.)*diffM[Np] +
				- nGeometry*G[Np]*G[Np]*urho[Np] +
				- (  mu[Np]*(G[Np-1]*urho[Np-1]-G[Np]*urho[Np])*grid.udxw[Np] -
						mu[Np-1]*(G[Np]*urho[Np]-G[Np-1]*urho[Np-1])*grid.udxw[Np])*(grid.udxw[Np]/2.) +
				- H[Np] );

	dG[Np] = 0.50*(G[Np-1]+G[Np]) - (-U[Np-1]+U[Np])*grid.udxw[Np];
	dH[Np] = U[Np]; 
}			

void OpenSMOKE_Flame1D::give_Opposed_DW(const std::string string_kind)
{
	int i, j;

	rhow[1] = rho[1];
	for(i=1;i<=Np-1;i++)
	{
		rhoe[i] = 0.50*(rho[i]+rho[i+1]);
		rhow[i+1] = rhoe[i];
	}
	rhoe[Np] = rho[Np];

	for(j=1;j<=NC;j++)
	{
		W.GetColumn(j, &auxNp);
		grid.FirstDerivative(data->iDerW, U, auxNp, diffWscalar);
		diffW.SetColumn(j, diffWscalar);
	}

	// Corrections for soot
	// --------------------------------------------------------------------------------------
	if ( data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT )
		R += SootGasCorrection;
	
	// Equazioni
	// --------------------------------------------------------------------------------------
	for(j=1;j<=NC;j++)
		dW[1][j] = rho[1] * ((nGeometry-1.)*U[1] * urho[1] * W[1][j] + vStar[1][j]) - BCW_C[j];
	
	for(i=2;i<=Ni;i++)
		for(j=1;j<=NC;j++)
			dW[i][j] = -((nGeometry-1.)*U[i] * diffW[i][j]
			               - R[i][j]
					       + (rhoe[i]*vStar[i][j] - rhow[i]*vStar[i-1][j])*grid.udxc_over_2[i]
						 ) * urho[i];

	for(j=1;j<=NC;j++)
		dW[Np][j] = rho[Np] * ((nGeometry-1.)*U[Np] * urho[Np] * W[Np][j] + vStar[Np][j]) - BCW_O[j];

	if(data->iDepositionWall == true)
	{	
		// Gas phase species
		for(j=1;j<=NC;j++)
			dW[Np][j] = (nGeometry-1.)*U[Np] * W[Np][j] + rho[Np] * vStar[Np][j] - BCW_O[j];

		// Soot species
		for(int j=1;j<=mix->polimiSoot->bin_indices().Size();j++)
		{
			const int jj = mix->polimiSoot->bin_indices()[j];
			dW[Np][jj] = (nGeometry-1.)*U[Np] * W[Np][jj] + rho[Np] * vStar[Np][jj] - rho[Np] * vThermophoretic[Np][jj];
		}
	}
	
	if (data->iPoolFire == POOL_FIRE_TASSIGNED || data->iPoolFire == POOL_FIRE_EQUILIBRIUM || data->iPoolFire == POOL_FIRE_LIQUIDPOOL)
	{
		for(j=1;j<=NC;j++)
			dW[1][j] = rho[1] * ((nGeometry-1.)*U[1] * urho[1] * W[1][j] + vStar[1][j]);
		dW[1][data->jFUEL] = (nGeometry-1.)*U[1] * (1. - W[1][data->jFUEL]) - rho[1] * vStar[1][data->jFUEL];
	}

	if (data->iSingleContributions == true)
	{
		for(int j=1;j<=data->index_SingleContributions.Size();j++)
		{
			int k = (j-1)*3;
			for(int i=2;i<=Ni;i++)
			{
				single_contributions[i][k+1]  = - (rhoe[i]*vStar[i][data->index_SingleContributions[j]] - rhow[i]*vStar[i-1][data->index_SingleContributions[j]])*grid.udxc_over_2[i];
				single_contributions[i][k + 2] = -(nGeometry-1.)*U[i] * diffW[i][data->index_SingleContributions[j]];
				single_contributions[i][k+3] =   R[i][data->index_SingleContributions[j]];
			}
		}
	}
}

void OpenSMOKE_Flame1D::give_Twin_DW(const std::string string_kind)
{
	give_Opposed_DW(string_kind);

	// Original set boundary conditions (twin flames)
	for(int j=1;j<=NC;j++)
		dW[Np][j] = -W[Np-1][j] + W[Np][j] ;
}


void OpenSMOKE_Flame1D::give_Opposed_DT(double t)
{
	int i, j;
	
	// -------------------------------------------------------------------------------------		
	// Variabili Ausiliarie
	// -------------------------------------------------------------------------------------
	sumCpDiffusive = 0.;
	for(i=2;i<=Ni;i++)
		for(j=1;j<=NC;j++)
//			sumCpDiffusive[i] += Cpk[i][j]*vStar[i][j];
			sumCpDiffusive[i] += 0.50*Cpk[i][j]*(vStar[i-1][j]+vStar[i][j]);
	
	// -------------------------------------------------------------------------------------		
	// Termini convettivi
	// -------------------------------------------------------------------------------------
	grid.FirstDerivative(data->iDerT, U, T, diffT);
	grid.FirstDerivative('C', T, diffTcentral);


	// --------------------------------------------------------------------------------------
	// Equazioni
	// --------------------------------------------------------------------------------------

	if (data->iTemperatureProfile==1)
	{
		dT[1] = T[1] - data->TC;
		for(i=2;i<=Ni;i++)
			dT[i] = 0.;
		dT[Np] = T[Np] - data->TO; 
	}
	else
	{
		dT[1] = T[1] - data->TC;

		if (data->iPoolFire == POOL_FIRE_EQUILIBRIUM || data->iPoolFire == POOL_FIRE_LIQUIDPOOL)
		{
			double Pvap = data->pool_fire_liquid_species->Pv(T[1]) * data->correctionFactorVaporPressure;
			dT[1] = Pvap/data->P_Pascal - X[1][data->jFUEL];
			data->TC = T[1];
		}

		for(i=2;i<=Ni;i++)
		
			dT[i] = -(+(nGeometry-1.)*U[i] * diffT[i] +
					    - (  lambda[i]*(T[i+1]-T[i])*grid.udxe[i] - lambda[i-1]*(T[i]-T[i-1])*grid.udxw[i] ) 
					      *grid.udxc_over_2[i] / Cp[i]
						- QReaction[i] / Cp[i] 
					    + rho[i]/Cp[i]*sumCpDiffusive[i]*diffTcentral[i] - Qrad[i]/Cp[i]
					  ) * urho[i] ;
		
		dT[Np] = T[Np] - data->TO; 
	}

	
}

void OpenSMOKE_Flame1D::give_Twin_DT(double t)
{
	give_Opposed_DT(t);

	dT[Np] = -(+(nGeometry-1.)*U[Np] * diffT[Np] * 0. +	// La derivata prima e' nulla 
				    - (  lambda[Np]*(T[Np-1]-T[Np])*grid.udxw[Np] - lambda[Np-1]*(T[Np]-T[Np-1])*grid.udxw[Np] ) 
				      *grid.udxc_over_2[Np/2] / Cp[Np]
					- QReaction[Np] / Cp[Np] 
				    + rho[Np]/Cp[Np]*sumCpDiffusive[Np]*diffTcentral[Np]*0 - 0*Qrad[Np]/Cp[Np]
				  ) * urho[Np] ;

	// Adiabatic condition
	dT[Np] = T[Np] - T[Np-1];
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									BACKUP												   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::recoverFromBackUp(const std::string fileName)
{
	char binary_string[40];
	BzzVector x_grid;

	std::string fileNameCase = fileName + ".inp";
	std::string fileNameData = fileName + ".bin";

	ptFlame = this;

	// Recover the kind of problem from backup data
	BzzLoad fInput('*', fileNameData);
		fInput.fileLoad.read((char*) binary_string, sizeof(binary_string));
		std::string kindOfProblem = binary_string;
		fInput >> Np;	ChangeDimensions(Np, &x_grid);
		fInput >> NC;
		fInput >> x_grid;
		data->Setup(kindOfProblem);
	fInput.End();

	setupGasMixture();

	if ( data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED || data->kind_of_flame == FLAME1D_PHYSICS_TWIN)	data->readFromFileForOpposed(fileNameCase);
	if ( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED)												data->readFromFileForPremixed(fileNameCase);
	

	grid.RecoverFromBackUp(x_grid);
	Np = grid.Np;
	Ni = grid.Ni;

	cout << "Backup not ready" << endl;

	// Setup of QMOM Module
	if (data->iQMOM == true)
		assign_qmom_module_from_file(*mix, data->qmom_file_name);

	// Setup of TwoEquations Module
	if (data->i2E == true)
	{
		soot2EModel.setupFromFile(data->twoEquation_file_name);
		soot2EModel.assign_mixture(*mix);
	}

	allocate_all();

	cout << "Radiation Setup" << endl;
	if ( data->iGasRadiation == true || data->iRadiativeSootModel != RADIATIVE_SOOT_MODEL_NONE )
		prepare_radiation();

	cout << "Boundary Conditions Setup" << endl;
	setupBoundaryConditions();
	
	cout << "Initial Conditions Setup" << endl;
	// Initial and Boundary conditions for QMOM
	if (data->iQMOM == true)
		initial_conditions_qmom_module();
		
	// Initial and Boundary conditions for TwoEquations Module
	if (data->i2E == true)
		initial_conditions_soot_module();

	cout << "Recovering variables from backup" << endl;	
	recoverVariables(fileNameData);
	cout << "Variables correctly recovered..." << endl;	
}

void OpenSMOKE_Flame1D::recoverVariables(const std::string fileName)
{
	char binary_string[40];
	
	BzzLoad fInput('*', fileName);

	fInput.fileLoad.read((char*) binary_string, sizeof(binary_string));
	
	fInput >> Np;
	fInput >> NC;
	fInput >> grid.x;

	fInput >> U;
	fInput >> G;
	fInput >> H;
	fInput >> T;
	fInput >> W;
		
	if (data->iQMOM == true)
		fInput >> moments;
			
	if (data->i2E == true)
	{
		fInput >> phiN;
		fInput >> phiM;
	}

	molarFractionsAndPMtot();

	if (iUnsteadyFromBackUp==-1)
	{
	}

	fInput.End();
}

void OpenSMOKE_Flame1D::printBackUpOnlyInputData(const std::string fileName)
{
	std::string fileNameCase = fileName + ".inp";

	if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED || 
		data->kind_of_flame == FLAME1D_PHYSICS_TWIN)
	{
		double vFuel =  (nGeometry-1.)*U[1] / rho[1];
		double vAir  = -(nGeometry-1.)*U[Np]/ rho[Np];

		data->printFileForOpposed(fileNameCase, vFuel, vAir);
	}

	if (data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED)
		data->printFileForPremixed(fileNameCase, H[1]);
}

void OpenSMOKE_Flame1D::printBackUpOnlyData(const std::string fileName)
{
	char binary_string[40];
	std::string kindOfProblem;

	if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED)
	{
		     if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)	kindOfProblem = "OPPOSED";
		else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)	kindOfProblem = "OPPOSED_SOOT";
		else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)	kindOfProblem = "OPPOSED_QMOM";
	}
	else if (data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED)
	{
		     if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)	kindOfProblem = "PREMIXED";
		else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)	kindOfProblem = "PREMIXED_SOOT";
		else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)	kindOfProblem = "PREMIXED_QMOM";
	}
	else if (data->kind_of_flame == FLAME1D_PHYSICS_TWIN)
	{
		     if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)	kindOfProblem = "TWIN";
		else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)	kindOfProblem = "TWIN_SOOT";
		else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)	kindOfProblem = "TWIN_QMOM";
	}

	std::string fileNameData = fileName + ".bin";

	BzzSave fOutput('*', fileNameData);
	strcpy(binary_string, kindOfProblem.c_str());
	fOutput.fileSave.write((char*) binary_string, sizeof(binary_string));
	fOutput << Np;
	fOutput << NC;
	fOutput << grid.x;
	fOutput << U;
	fOutput << G;
	fOutput << H;
	fOutput << T;
	fOutput << W;
		
	if (data->iQMOM == true)
		fOutput << moments;
		
	if (data->i2E == true)
	{
		fOutput << phiN;
		fOutput << phiM;
	}
	
	if (data->iUnsteady == true)
	{
	}

	fOutput.End();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								NLS SYSTEMs - Premixed											//
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::nonLinearSystem_Premixed_All(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, -1, dummy, 0, 0);
		compute_vStarOpposed();
		compute_vStarPremixed();

	for(int i=1;i<=Np;i++)
		U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DT();
		give_Premixed_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_ALL, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Premixed_NoEnergy(BzzVector &x, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_NOENERGY, x);
	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT , tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DW("NLS");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_NOENERGY, f);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								ODE SYSTEMs - Premixed											//
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::ODESystem_Premixed_All(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		properties(data->iHOT, ODE_Premixed_ALL.jacobianIndex, ODE_Premixed_ALL.jacobianVariables, NC+1, 1);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DT();
		give_Premixed_DW("ODE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_ALL, f);
}


void OpenSMOKE_Flame1D::ODESystem_Premixed_NoEnergy(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_NOENERGY, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DW("ODE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_NOENERGY, f);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								DAE SYSTEMs - Premixed											//
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::DAESystem_Premixed_All(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		properties(data->iHOT, DAE_Premixed_ALL.jacobianIndex, DAE_Premixed_ALL.jacobianVariables, NC+1, 1);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DT();
		give_Premixed_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_ALL, f);
}

void OpenSMOKE_Flame1D::DAESystem_Premixed_FlameSpeed(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_FLAMESPEED, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		properties(data->iHOT, DAE_Premixed_FLAMESPEED.jacobianIndex, DAE_Premixed_FLAMESPEED.jacobianVariables, NC+2, 1);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DT();
		give_Premixed_DW("DAE");
		give_Premixed_DH();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_FLAMESPEED, f);
}

void OpenSMOKE_Flame1D::nonLinearSystem_Premixed_FlameSpeed(BzzVector &x, BzzVector &f)
{
		Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_FLAMESPEED, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		BzzVectorInt dummy;
		properties(data->iHOT, -1, dummy, 0, 0);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DT();
		give_Premixed_DW("NLS");
		give_Premixed_DH();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_FLAMESPEED, f);
}

void OpenSMOKE_Flame1D::DAESystem_Premixed_NoEnergy(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_NOENERGY, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)
		
		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_NOENERGY, f);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								GIVE RESIDUALS - Premixed											//
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Flame1D::give_Premixed_DW(const std::string string_kind)
{
	int i, j;

	// 1. Computation: A x rho
	// ------------------------------------------------------------
	A_x_rhow[1] = G[1]*rho[1];
	for(i=1;i<=Np-1;i++)
	{
		A_x_rhoe[i] = 0.50*(G[i]*rho[i]+G[i+1]*rho[i+1]);
		A_x_rhow[i+1] = A_x_rhoe[i];
	}
	A_x_rhoe[Np] = G[Np]*rho[Np];


	// 2. Convective terms
	// ------------------------------------------------------------
	for(j=1;j<=NC;j++)
	{
		W.GetColumn(j, &auxNp);
		grid.FirstDerivative(data->iDerW, U, auxNp, diffWscalar);
		diffW.SetColumn(j, diffWscalar);
	}
	
	// Corrections for soot
	// --------------------------------------------------------------------------------------
	if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)
		R += SootGasCorrection;

	// Equations
	// --------------------------------------------------------------------------------------

	for(j=1;j<=NC;j++)
		dW[1][j] = W[1][j] + A_x_rhoe[1]/H[1]*vStar_e[1][j] - BCW_C[j];
//	{
	//	cout << "Boundary without retrodiffusion" << endl;
	//	dW[1][j] = W[1][j] - BCW_C[j];
//	}

	for(i=2;i<=Ni;i++)
	{
		for(j=1;j<=NC;j++)
			dW[i][j] = - (   H[i]*diffW[i][j]
			               - G[i]*R[i][j]
					       + (A_x_rhoe[i]*vStar_e[i][j] - A_x_rhow[i]*vStar_w[i][j])*grid.udxc_over_2[i]
						 ) / (rho[i]*G[i]);
	}

	for(j=1;j<=NC;j++)
		dW[Np][j] = -W[Np-1][j] + W[Np][j];	

	// Boundary equations correction for ODE system
	// ---------------------------------------------------------------------
	if (string_kind == "ODE")
		for(j=1;j<=NC;j++)
		{
			dW[1][j] = 0.;
			dW[Np][j] = 0.;
		}

	if (data->iSingleContributions == true)
	{
		for(int j=1;j<=data->index_SingleContributions.Size();j++)
		{
			int k = (j-1)*3;
			for(int i=2;i<=Ni;i++)
			{
				single_contributions[i][k+1] = - (A_x_rhoe[i]*vStar_e[i][data->index_SingleContributions[j]] - A_x_rhow[i]*vStar_w[i][data->index_SingleContributions[j]])*grid.udxc_over_2[i]/G[i];
				single_contributions[i][k+2] = - H[i]*diffW[i][data->index_SingleContributions[j]]/G[i];
				single_contributions[i][k+3] =   R[i][data->index_SingleContributions[j]];
			}
		}
	}
}

void OpenSMOKE_Flame1D::give_Premixed_DT()
{
	int i, j;

	if (data->iAssignedFixedTemperatureProfileProvisional == true)
	{
		dT = 0.;
		dT[1] = T[1] - data->TC;
		dT[Np] = T[Np] - data->ud_temperature_profile.GiveMeValue(0., grid.x[Np]);;
	}
	else
	{

	// 1. Computation: A x lambda
	// ------------------------------------------------------------
	A_x_lambdaw[1] = G[1]*lambda[1];
	for(i=1;i<=Np-1;i++)
	{
		A_x_lambdae[i] = 0.50 * (G[i]*lambda[i] + G[i+1]*lambda[i+1]);
		A_x_lambdaw[i+1] = A_x_lambdae[i];
	}
	A_x_lambdae[Np] = G[Np]*lambda[Np];


	// 2. Computation: mass diffusional fluxes
	// ------------------------------------------------------------
	sumCpDiffusive = 0.;
	for(i=2;i<=Ni;i++)
		for(j=1;j<=NC;j++)
//			sumCpDiffusive[i] += Cpk[i][j]*vStar[i][j];
			sumCpDiffusive[i] += 0.50*Cpk[i][j]*(vStar[i-1][j]+vStar[i][j]);
	

	// 3. Computation: convective terms
	// ------------------------------------------------------------
	grid.FirstDerivative(data->iDerT, U, T, diffT);
	grid.FirstDerivative('C', T, diffTcentral);


	// 4. Equations
	// ------------------------------------------------------------
	dT[1] = T[1] - data->TC;

	for(i=2;i<=Ni;i++)
	{
		dT[i] = - ( + H[i] * diffT[i] + 
				    - (A_x_lambdae[i]*(T[i+1]-T[i])*grid.udxe[i] - A_x_lambdaw[i]*(T[i]-T[i-1])*grid.udxw[i] ) 
				      *grid.udxc_over_2[i] / Cp[i]
					- G[i]*QReaction[i] / Cp[i] 
				    + G[i]*rho[i]/Cp[i]*sumCpDiffusive[i]*diffTcentral[i]
				  ) / (rho[i]*G[i]);
		
		dT[i]  += Qrad[i]/Cp[i]/rho[i] ;
	}

	dT[Np] = T[Np] - T[Np-1];
	}
}

void OpenSMOKE_Flame1D::give_Premixed_DH()
{
	int i;

	for(i=1;i<data->iFixedTemperature+1;i++)
		dH[i] = H[i] - H[i+1];

	dH[data->iFixedTemperature] = T[data->iFixedTemperature] - data->fixedTemperature;
//	dH[data->iFixedTemperature] = T[data->iFixedTemperature] - data->fixedTemperature + 
//									 2.*H[data->iFixedTemperature] - H[data->iFixedTemperature-1] - H[data->iFixedTemperature+1];
	
	for(i=data->iFixedTemperature+1;i<=Np;i++)
		dH[i] = H[i-1] - H[i];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//								UTILITIES											//
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

int OpenSMOKE_Flame1D::nonLinearSystemSolution(BzzNonLinearSystemSparseObject &o, int dimensionBlock)
{
	int i;
	int iMaxfNLS;
	double sum, maxfNLS;
	double startTime, endTime;

	o.SetMinimumConstraints(xMin);
	o.SetMaximumConstraints(xMax);

	// Default values: (A) 1e-10      (R) sqrt(MachEpsFloat())
	o.SetTolerance(data->abs_nls_Tolerances, data->rel_nls_Tolerances);
	cout << "Absolute tolerance: " << data->abs_nls_Tolerances << endl;
	cout << "Relative tolerance: " << data->rel_nls_Tolerances << endl;

	startTime = BzzGetCpuTime();
	bzzStop = 0;
	char control = o();
	endTime = BzzGetCpuTime() - startTime;
	
	nonLinearSystemLabel(o, control);

	// Information
	// ---------------------------------------------
	o.GetSolution(&xFirstGuess, &fNLS);
	
	maxfNLS = fNLS.MaxAbs(&iMaxfNLS);
	
	sum = 0.;
	for (i=1;i<=(Np*dimensionBlock);i++)
		sum += fabs(fNLS[i]);
	sum=sum/(Np*dimensionBlock);

	cout << "Mean residual:" << sum << endl;
	cout << "Max. residual:" << maxfNLS << endl;
	cout << "Time NLS solution: " <<  endTime << " s" << endl;
	{
		int iPoint = iMaxfNLS/dimensionBlock+1;

		cout << "Index Max Residual: " << iMaxfNLS << " - Point: " << iPoint << " Variable: " << iMaxfNLS%dimensionBlock << endl;
		cout << "x[cm]: " << grid.x[iPoint]*100. << endl;
		cout << "U: " << U[iPoint-1] << " " << U[iPoint] << " " << U[iPoint+1] << endl;
		cout << "G: " << G[iPoint-1] << " " << G[iPoint] << " " << G[iPoint+1] << endl;
		cout << "H: " << H[iPoint-1] << " " << H[iPoint] << " " << H[iPoint+1] << endl;
	}

	cout << endl;

	return int(control);
}

int OpenSMOKE_Flame1D::nonLinearSystemSolution(BzzNonLinearSystemObject &o)
{
	int i;
	int iMaxfNLS;
	double sum, maxfNLS;
	double startTime, endTime;

	o.SetMinimumConstraints(xMin);
	o.SetMaximumConstraints(xMax);

	// Default values: (A) 1e-10      (R) sqrt(MachEpsFloat())
	o.SetTolerance(data->abs_nls_Tolerances, data->rel_nls_Tolerances);
	cout << "Absolute tolerance: " << data->abs_nls_Tolerances << endl;
	cout << "Relative tolerance: " << data->rel_nls_Tolerances << endl;

	startTime = BzzGetCpuTime();
	bzzStop = 0;
	char control = o();
	endTime = BzzGetCpuTime() - startTime;
	
	nonLinearSystemLabel(o, control);

	// Information
	// ---------------------------------------------
	o.GetSolution(&xFirstGuess, &fNLS);
	
	maxfNLS = fNLS.MaxAbs(&iMaxfNLS);
	
	sum = 0.;
	for (i=1;i<=fNLS.Size();i++)
		sum += fabs(fNLS[i]);
	sum=sum/double(fNLS.Size());

	cout << "Mean residual:" << sum << endl;
	cout << "Max. residual:" << maxfNLS << endl;
	cout << "Time NLS solution: " <<  endTime << " s" << endl;
	cout << endl;

	return int(control);
}

//TODO
void OpenSMOKE_Flame1D::nonLinearSystemSolutionWithControl_FlameSpeed(const double time_first)
{
	// Check
	if (data->iUnsteady == true)
		ErrorMessage("NLS system cannot be solved for unsteady problems...");

	double tEnd					= 1.e-4;
	int dimensionBlock			= NC+2;
	int max_number_of_guesses	= 6;

	// Memory allocation
	ChangeDimensions(Np*dimensionBlock, &fNLS);
	ChangeDimensions(Np*dimensionBlock, &xFirstGuess);
	ChangeDimensions(Np*dimensionBlock, &xMin);
	ChangeDimensions(Np*dimensionBlock, &xMax);
	ChangeDimensions(Np*dimensionBlock, &inDerAlg);

	// Set minimum values and differential-algebraic solutions
	setMinimumAndMaximumValues(PREMIXED_FLAMESPEED);
	setDifferentialAndAlgebraic(PREMIXED_FLAMESPEED);

	// Store original solution
	OpenSMOKE_Flame1D_Solution original_solution;
	original_solution.PasteFromExternalSolution(*this, *data);

	data->iHOT = 1;
	if (time_first > 0.)	tEnd = time_first;
	int number_of_guesses = 1;
	for(;;)
	{
		if (number_of_guesses > max_number_of_guesses)
			ErrorMessage("Maximum number of guesses for non linear system solution...");

		// Differential-Algebraic system solution
		if (time_first == 0. && number_of_guesses == 1)	
		{
		}
		else
		{
			cout << "-----------------------------------------------------------------------" << endl;
			cout << " DAE Solution" << endl;
			cout << "  Npoints:   " << Np	<< endl;
			cout << "  Guess:     " << number_of_guesses	<< endl;
			cout << "  Time:      " << tEnd << " s"			<< endl;
			cout << "-----------------------------------------------------------------------" << endl;

			PasteFromExternalSolution(original_solution);

			int k=1;
			for(int i=1;i<=Np;i++)
			{
				xFirstGuess[k++] = T[i];
				xFirstGuess[k++] = H[i];
				for(int j=1;j<=NC;j++)
					xFirstGuess[k++] = W[i][j];
			}

			xFirstGuess( dimensionBlock*(data->iFixedTemperature-1) + 1) = data->fixedTemperature;
					
			DAE_Premixed_FLAMESPEED.assignFlame(this);
			BzzDaeSparseObject dae_object;
			dae_object(xFirstGuess, 0., inDerAlg, &DAE_Premixed_FLAMESPEED, dimensionBlock);
			DAESystemSolution(&dae_object, tEnd);

			original_solution.PasteFromExternalSolution(*this, *data);
		}

		// Non linear system solution block
		{
			cout << "-----------------------------------------------------------------------" << endl;
			cout << " NLS Solution" << endl;
			cout << "  Npoints:   " << Np	<< endl;
			cout << "  Guess:     " << number_of_guesses	<< endl;
			cout << "-----------------------------------------------------------------------" << endl;
	
			int k=1;
			for(int i=1;i<=Np;i++)
			{
				xFirstGuess[k++] = T[i];
				xFirstGuess[k++] = H[i];
				for(int j=1;j<=NC;j++)
					xFirstGuess[k++] = W[i][j];
			}

			xFirstGuess( dimensionBlock*(data->iFixedTemperature-1) + 1) = data->fixedTemperature;

			NLS_Premixed_FLAMESPEED.assignFlame(this);
			BzzNonLinearSystemSparseObject nls(xFirstGuess, &NLS_Premixed_FLAMESPEED, dimensionBlock);
			char control = nonLinearSystemSolution(nls, dimensionBlock);
			
			if (T.Min()>0. && T.Max()<6000.)
			{
				if (int(control)>=2 && int(control)<=6)	
				{
					nls.GetSolution(&xFirstGuess, &fNLS);
					break;
				}
			}
		}
		
		cout << endl;
		cout << "It was impossibile to reach a solution!" << endl;
		cout << "A better first guess solution is required!" << endl;
		cout << endl;

		// In case of failure
		number_of_guesses++;
		tEnd *= 10.;
	}

	// Information
	// ---------------------------------------------
	double sum = 0.;
	for (int i=1;i<=(Np*dimensionBlock);i++)
		sum += fabs(fNLS[i]);
	sum=sum/(Np*dimensionBlock);

	cout << "Mean residual:" << sum << endl;
	cout << "Max. residual:" << fNLS.MaxAbs() << endl;
	cout << endl;
}

//TODO
void OpenSMOKE_Flame1D::nonLinearSystemSolutionWithControl_Opposed(const double time_first)
{
	// Check
	if (data->iUnsteady == true)
		ErrorMessage("NLS system cannot be solved for unsteady problems...");

	double tEnd					= 1.e-4;
	int dimensionBlock			= NC+4;
	int max_number_of_guesses	= 6;

	// Memory allocation
	ChangeDimensions(Np*dimensionBlock, &fNLS);
	ChangeDimensions(Np*dimensionBlock, &xFirstGuess);
	ChangeDimensions(Np*dimensionBlock, &xMin);
	ChangeDimensions(Np*dimensionBlock, &xMax);
	ChangeDimensions(Np*dimensionBlock, &inDerAlg);

	// Set minimum values and differential-algebraic solutions
	setMinimumAndMaximumValues(OPPOSED_ALL);
	setDifferentialAndAlgebraic(OPPOSED_ALL);

	// Store original solution
	OpenSMOKE_Flame1D_Solution original_solution;
	original_solution.PasteFromExternalSolution(*this, *data);

	data->iHOT = 1;
	if (time_first > 0.)	tEnd = time_first;
	int number_of_guesses = 1;
	for(;;)
	{
		if (number_of_guesses > max_number_of_guesses)
			ErrorMessage("Maximum number of guesses for non linear system solution...");

		// Differential-Algebraic system solution
		if (time_first == 0. && number_of_guesses == 1)	
		{
		}
		else
		{
			cout << "-----------------------------------------------------------------------" << endl;
			cout << " DAE Solution" << endl;
			cout << "  Npoints:   " << Np	<< endl;
			cout << "  Guess:     " << number_of_guesses	<< endl;
			cout << "  Time:      " << tEnd << " s"			<< endl;
			cout << "-----------------------------------------------------------------------" << endl;

			PasteFromExternalSolution(original_solution);

			int k=1;
			for(int i=1;i<=Np;i++)
			{
				xFirstGuess[k++] = U[i];
				xFirstGuess[k++] = G[i];
				xFirstGuess[k++] = H[i];
				xFirstGuess[k++] = T[i];
				for(int j=1;j<=NC;j++)
					xFirstGuess[k++] = W[i][j];
			}

			DAE_Opposed_ALL.assignFlame(this);
			BzzDaeSparseObject dae_object;
			dae_object(xFirstGuess, 0., inDerAlg, &DAE_Opposed_ALL, dimensionBlock);
			DAESystemSolution(&dae_object, tEnd);

			original_solution.PasteFromExternalSolution(*this, *data);
		}

		// Non linear system solution block
		{
			cout << "-----------------------------------------------------------------------" << endl;
			cout << " NLS Solution" << endl;
			cout << "  Npoints:   " << Np	<< endl;
			cout << "  Guess:     " << number_of_guesses	<< endl;
			cout << "-----------------------------------------------------------------------" << endl;
	
			int k=1;
			for(int i=1;i<=Np;i++)
			{
				xFirstGuess[k++] = U[i];
				xFirstGuess[k++] = G[i];
				xFirstGuess[k++] = H[i];
				xFirstGuess[k++] = T[i];
				for(int j=1;j<=NC;j++)
					xFirstGuess[k++] = W[i][j];
			}

			NLS_Opposed_ALL.assignFlame(this);
			BzzNonLinearSystemSparseObject nls(xFirstGuess, &NLS_Opposed_ALL, dimensionBlock);
			char control = nonLinearSystemSolution(nls, dimensionBlock);
			
			if (T.Min()>0. && T.Max()<6000.)
			{
				if (int(control)>=2 && int(control)<=6)	
				{
					nls.GetSolution(&xFirstGuess, &fNLS);
					break;
				}
			}
		}
		
		cout << endl;
		cout << "It was impossibile to reach a solution!" << endl;
		cout << "A better first guess solution is required!" << endl;
		cout << endl;

		// In case of failure
		number_of_guesses++;
		tEnd *= 10.;
	}

	// Information
	// ---------------------------------------------
	double sum = 0.;
	for (int i=1;i<=(Np*dimensionBlock);i++)
		sum += fabs(fNLS[i]);
	sum=sum/(Np*dimensionBlock);

	cout << "Mean residual:" << sum << endl;
	cout << "Max. residual:" << fNLS.MaxAbs() << endl;
	cout << endl;
}

//TODO
void OpenSMOKE_Flame1D::nonLinearSystemSolutionWithControl_Twin(const double time_first)
{
	// Check
	if (data->iUnsteady == true)
		ErrorMessage("NLS system cannot be solved for unsteady problems...");

	double tEnd					= 1.e-4;
	int dimensionBlock			= NC+4;
	int max_number_of_guesses	= 6;

	// Memory allocation
	ChangeDimensions(Np*dimensionBlock, &fNLS);
	ChangeDimensions(Np*dimensionBlock, &xFirstGuess);
	ChangeDimensions(Np*dimensionBlock, &xMin);
	ChangeDimensions(Np*dimensionBlock, &xMax);
	ChangeDimensions(Np*dimensionBlock, &inDerAlg);

	// Set minimum values and differential-algebraic solutions
	setMinimumAndMaximumValues(TWIN_ALL);
	setDifferentialAndAlgebraic(TWIN_ALL);

	// Store original solution
	OpenSMOKE_Flame1D_Solution original_solution;
	original_solution.PasteFromExternalSolution(*this, *data);

	data->iHOT = 1;
	if (time_first > 0.)	tEnd = time_first;
	int number_of_guesses = 1;
	for(;;)
	{
		if (number_of_guesses > max_number_of_guesses)
			ErrorMessage("Maximum number of guesses for non linear system solution...");

		// Differential-Algebraic system solution
		if (time_first == 0. && number_of_guesses == 1)	
		{
		}
		else
		{
			cout << "-----------------------------------------------------------------------" << endl;
			cout << " DAE Solution" << endl;
			cout << "  Npoints:   " << Np	<< endl;
			cout << "  Guess:     " << number_of_guesses	<< endl;
			cout << "  Time:      " << tEnd << " s"			<< endl;
			cout << "-----------------------------------------------------------------------" << endl;

			PasteFromExternalSolution(original_solution);

			int k=1;
			for(int i=1;i<=Np;i++)
			{
				xFirstGuess[k++] = U[i];
				xFirstGuess[k++] = G[i];
				xFirstGuess[k++] = H[i];
				xFirstGuess[k++] = T[i];
				for(int j=1;j<=NC;j++)
					xFirstGuess[k++] = W[i][j];
			}

			DAE_Twin_ALL.assignFlame(this);
			BzzDaeSparseObject dae_object;
			dae_object(xFirstGuess, 0., inDerAlg, &DAE_Twin_ALL, dimensionBlock);
			DAESystemSolution(&dae_object, tEnd);

			original_solution.PasteFromExternalSolution(*this, *data);
		}

		// Non linear system solution block
		{
			cout << "-----------------------------------------------------------------------" << endl;
			cout << " NLS Solution" << endl;
			cout << "  Npoints:   " << Np	<< endl;
			cout << "  Guess:     " << number_of_guesses	<< endl;
			cout << "-----------------------------------------------------------------------" << endl;
	
			int k=1;
			for(int i=1;i<=Np;i++)
			{
				xFirstGuess[k++] = U[i];
				xFirstGuess[k++] = G[i];
				xFirstGuess[k++] = H[i];
				xFirstGuess[k++] = T[i];
				for(int j=1;j<=NC;j++)
					xFirstGuess[k++] = W[i][j];
			}

			NLS_Twin_ALL.assignFlame(this);
			BzzNonLinearSystemSparseObject nls(xFirstGuess, &NLS_Twin_ALL, dimensionBlock);
			char control = nonLinearSystemSolution(nls, dimensionBlock);
			
			if (T.Min()>0. && T.Max()<6000.)
			{
				if (int(control)>=2 && int(control)<=6)	
				{
					nls.GetSolution(&xFirstGuess, &fNLS);
					break;
				}
			}
		}
		
		cout << endl;
		cout << "It was impossibile to reach a solution!" << endl;
		cout << "A better first guess solution is required!" << endl;
		cout << endl;

		// In case of failure
		number_of_guesses++;
		tEnd *= 10.;
	}

	// Information
	// ---------------------------------------------
	double sum = 0.;
	for (int i=1;i<=(Np*dimensionBlock);i++)
		sum += fabs(fNLS[i]);
	sum=sum/(Np*dimensionBlock);

	cout << "Mean residual:" << sum << endl;
	cout << "Max. residual:" << fNLS.MaxAbs() << endl;
	cout << endl;
}


void OpenSMOKE_Flame1D::nonLinearSystemLabel(BzzNonLinearSystemSparseObject &o, char &control)
{
	cout << "Maximum temperature: " << T.Max() << " K" << endl;
	cout << "Number of functions for Numerical Jacobian: " << o.NumFunctionsForNumericalJacobian() << endl;
	cout << "Numerical Jacobians: " << o.NumNumericalJacobians() << endl;
	cout << "Number of Newtons: " << o.NumNewtons() << endl;
	cout << "Number of Quasi Newton: " << o.NumQuasiNewtons() << endl;
	
	cout << "State " << int(control) << ": " << endl;
	if (control==1)			cout << " ATTENTION! The maximum number of functions calls has been performed!" << endl;
	else if (control==2)	cout << " The Newton correction has reached the required precision!" << endl;
	else if (control==3)	cout << " The Quasi Newton correction has reached the required precision!" << endl;
	else if (control==4)	cout << " The Gradient has reached the required precision!" << endl;
	else if (control==5)	cout << " The Objective Function phiNew has reached the required precision!" << endl;
	else if (control==6)	cout << " The Objective Function phiW has reached the required precision!" << endl;
	else if (control==7)	cout << " ATTENTION! The Objective Function phiNew has reached the required precision but the solution is dubious!" << endl;
	else if (control==8)	cout << " ATTENTION! Reached the assigned max value for Newton calls!" << endl;
	else if (control==-1)	cout << " ATTENTION! Impossible to reach the solution!" << endl;
	else if (control==-2)	cout << " ATTENTION! The search has been stopped!" << endl;
	else if (control==-3)	cout << " ATTENTION! The object is not initialized!" << endl;
	else if (control==-4)	cout << " ATTENTION! Impossible to reach the solution in Restart!" << endl;

	cout << endl;
}

void OpenSMOKE_Flame1D::nonLinearSystemLabel(BzzNonLinearSystemObject &o, char &control)
{
	cout << "Maximum temperature: " << T.Max() << " K" << endl;
	cout << "Number of functions for Numerical Jacobian: " << o.NumFunctionsForNumericalJacobian() << endl;
	cout << "Numerical Jacobians: " << o.NumNumericalJacobians() << endl;
	cout << "Number of Newtons: " << o.NumNewtons() << endl;
	cout << "Number of Quasi Newton: " << o.NumQuasiNewtons() << endl;
	
	cout << "State " << int(control) << ": " << endl;
	if (control==1)			cout << " ATTENTION! The maximum number of functions calls has been performed!" << endl;
	else if (control==2)	cout << " The Newton correction has reached the required precision!" << endl;
	else if (control==3)	cout << " The Quasi Newton correction has reached the required precision!" << endl;
	else if (control==4)	cout << " The Gradient has reached the required precision!" << endl;
	else if (control==5)	cout << " The Objective Function phiNew has reached the required precision!" << endl;
	else if (control==6)	cout << " The Objective Function phiW has reached the required precision!" << endl;
	else if (control==7)	cout << " ATTENTION! The Objective Function phiNew has reached the required precision but the solution is dubious!" << endl;
	else if (control==8)	cout << " ATTENTION! Reached the assigned max value for Newton calls!" << endl;
	else if (control==-1)	cout << " ATTENTION! Impossible to reach the solution!" << endl;
	else if (control==-2)	cout << " ATTENTION! The search has been stopped!" << endl;
	else if (control==-3)	cout << " ATTENTION! The object is not initialized!" << endl;
	else if (control==-4)	cout << " ATTENTION! Impossible to reach the solution in Restart!" << endl;

	cout << endl;
}

void OpenSMOKE_Flame1D::DAESystemSolution(BzzDaeSparseObject *o, double tEnd)
{	
	o->StepPrint(DAE_ODE_Print);
	o->SetMinimumConstraints(xMin);
	o->SetMaximumConstraints(xMax);

	// Default values: (A) 1e-10      (R) 100*MachEpsFloat()
	o->SetTolRel(data->rel_dae_Tolerances);
	o->SetTolAbs(data->abs_dae_Tolerances);
	cout << "Absolute tolerance: " << data->abs_dae_Tolerances << endl;
	cout << "Relative tolerance: " << data->rel_dae_Tolerances << endl;

	if (data->initial_time_step > 0)
	{
		cout << "Initial time step: " << data->initial_time_step << endl;
		o->SetH0(data->initial_time_step);
	}

	if (data->max_integration_order > 0)
	{
		cout << "Maximum integration order: " << data->max_integration_order << endl;
		o->SetMaxOrder(data->max_integration_order);
	}

	double timeStart = BzzGetCpuTime();
	bzzStop = 0;
	xFirstGuess = (*o)(tEnd, tEnd);
	
	cout << endl;
	cout << "Number of steps: "					<< o->GetNumStep() << endl;
	cout << "Number of function for Jacobian: " << o->GetNumFunctionForJacobian() << endl;
	cout << "Numerical Jacobians: "				<< o->GetNumNumericalJacobian() << endl;
	cout << "Time DAE solution: "				<< BzzGetCpuTime() - timeStart << " s" << endl << endl;

	if (data->iUnsteady == true)
		fUnsteady.close();
}

void OpenSMOKE_Flame1D::ODESystemSolution(BzzOdeSparseStiffObject &o, double tEnd)
{
	o.StepPrint(DAE_ODE_Print);
	o.SetMinimumConstraints(xMin);
	o.SetMaximumConstraints(xMax);
	
	double timeStart = BzzGetCpuTime();
	bzzStop = 0;
	xFirstGuess = o(tEnd);
	
	cout << endl;
	cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian() << endl;
	cout << "Numerical Jacobians: " << o.GetNumNumericalJacobian() << endl;
	cout << "Time DAE solution: " << BzzGetCpuTime() - timeStart << " s" << endl << endl;

	if (data->iUnsteady == true)
		fUnsteady.close();
}

void OpenSMOKE_Flame1D::assign_qmom_module_from_file(OpenSMOKE_ReactingGas &mix, const std::string fileName)
{
	qmom.setupGasMixture(mix);
	qmom.readFromFile(fileName);
	allocate_QMOM_N();

	qmom.soot.assign_mixture(mix);
}

void OpenSMOKE_Flame1D::initial_conditions_qmom_module()
{
	int i;

	// Definition of Inlet-Initial Conditions [units: m]
	momentsC[1] = qmom.seed_density;
	for(i=2;i<=2*qmom.N;i++)
		momentsC[i] = qmom.seed_density*pow(qmom.epsilon, i-1)/double(i);

	// Normalization of moments
	double coeff = 1.;
	for (i=1;i<=2*qmom.N;i++)
	{
		momentsC[i]/=(qmom.seed_density*coeff);
		coeff *= qmom.epsilon;
	}

	for(i=1;i<=Np;i++)
		moments.SetRow(i, momentsC);
}

void OpenSMOKE_Flame1D::initial_conditions_soot_module()
{
	int i;

	soot2EModel.initial_values(rho[1]);
	
	phiNC = soot2EModel.phiNStart;
	phiMC = soot2EModel.phiMStart;
	phiNO = soot2EModel.phiNStart;
	phiMO = soot2EModel.phiMStart;

	for (i=1;i<=Np;i++)
	{
		soot2EModel.initial_values(rho[i]);
		phiN = soot2EModel.phiNStart;
		phiM = soot2EModel.phiMStart;
	}
}

void OpenSMOKE_Flame1D::give_Premixed_QMOM()
{
	int i, j;

	// 0. Diffusion coefficients
	// ------------------------------------------------------------
	DiffusionMoments = 0.;

	// 1. Computation: A x lambda
	// ------------------------------------------------------------
	for(j=1;j<=2*qmom.N;j++)
	{
		A_x_rho_x_Diffw[1][j] = G[1]*rho[1]*DiffusionMoments[1][j];
		for(i=1;i<=Ni;i++)
		{
			A_x_rho_x_Diffe[i][j]	= 0.50 * (	G[i]*rho[i]*DiffusionMoments[i][j] + 
												G[i+1]*rho[i+1]*DiffusionMoments[i+1][j]);
			A_x_rho_x_Diffw[i+1][j] = A_x_rho_x_Diffe[i][j];
		}
		A_x_rho_x_Diffe[Np][j] = G[Np]*rho[Np]*DiffusionMoments[Np][j];
	}	

	// 2. Computation: convective terms
	// ------------------------------------------------------------
	grid.FirstDerivative('B', moments, diffMoments);


	// 3. Moment Source Terms
	// ------------------------------------------------------------
	double mu=1.0;double tr=1.0;	// they are used only if the fractal dimension is enabled
	for(i=2;i<=Np;i++)
	{
		BzzVector momentsVector = moments.GetRow(i);
		qmom.updateMoments(momentsVector);

		double valueFluctuations = 0;
		qmom.updateData( T[i], data->P_Pascal, rho[i], mu, tr, 
						 W[i][qmom.jO2], W[i][qmom.jC2H2], W[i][qmom.jOH], valueFluctuations );

		qmom_sources = qmom.calculateSources();
		moments_source.SetRow(i, qmom_sources);
		qmom.soot.formation_rates();

		if (i!=Np)
		{
			for(j=1;j<=NC-1;j++)
				dW[i][j] += qmom.soot.SGas[j]/rho[i];
			dW[i][NC] += qmom.soot.S/rho[i];
		}
	}
	
	// 3. Equations
	// ------------------------------------------------------------
	for(j=1;j<=2*qmom.N;j++)
		dmoments[1][j] = moments[1][j] - momentsC[j];

	for(j=1;j<=2*qmom.N;j++)
		for(i=2;i<=Ni;i++)
			dmoments[i][j] =  - ( H[i] * diffMoments[i][j] 
		                     - (A_x_rho_x_Diffe[i][j]*(moments[i+1][j]-moments[i][j])*grid.udxe[i] - 
						        A_x_rho_x_Diffw[i][j]*(moments[i][j]-moments[i-1][j])*grid.udxw[i] ) 
						       *grid.udxc_over_2[i] )
						  /(rho[i]*G[i]) + 
						  moments_source[i][j];
	

	// This is the correct boundary condition if the diffusivity is set to zero (Differential equation)
	for(j=1;j<=2*qmom.N;j++)
		dmoments[Np][j] =  - ( H[Np] * diffMoments[Np][j] ) /(rho[Np]*G[Np]) + 
											moments_source[Np][j];

}

void OpenSMOKE_Flame1D::give_Opposed_QMOM()
{
	int i,j;

	// 0. Diffusion coefficients
	// ------------------------------------------------------------
	DiffusionMoments = 1.e-7;

	// 1. Moment Source Terms
	// ------------------------------------------------------------
	double mu=1.0;double tr=1.0;	// they are used only if the fractal dimension is enabled
	for(i=2;i<=Ni;i++)
	{
		BzzVector momentsVector = moments.GetRow(i);
		qmom.updateMoments(momentsVector);

		double valueFluctuations = 0;
		qmom.updateData( T[i], data->P_Pascal, rho[i], mu, tr, 
						 W[i][qmom.jO2], W[i][qmom.jC2H2], W[i][qmom.jOH], valueFluctuations );

		qmom_sources = qmom.calculateSources();
		moments_source.SetRow(i, qmom_sources);
		qmom.soot.formation_rates();

		for(j=1;j<=NC-1;j++)
			dW[i][j] +=  qmom.soot.SGas[j]/rho[i];
		dW[i][NC] +=  qmom.soot.S/rho[i];
	}
	
	// 2. Computation: convective terms
	// ------------------------------------------------------------
	grid.FirstDerivative('U', U, moments, diffMoments);
	
	// 3. Equations
	// ------------------------------------------------------------
	for(j=1;j<=2*qmom.N;j++)
		dmoments[1][j] = moments[1][j] - momentsC[j];

	for(i=2;i<=Ni;i++)
		for(j=1;j<=2*qmom.N;j++)
			dmoments[i][j] = (-(nGeometry-1.)*U[i] * diffMoments[i][j]) * urho[i] +
								moments_source[i][j] +
								( DiffusionMoments[i][j]  *(moments[i+1][j]-moments[i][j])*grid.udxe[i] - 
								  DiffusionMoments[i-1][j]*(moments[i][j]-moments[i-1][j])*grid.udxw[i] ) * grid.udxc_over_2[i];
	for(j=1;j<=2*qmom.N;j++)
		dmoments[Np][j] = moments[Np][j] - momentsC[j];
}

void OpenSMOKE_Flame1D::give_Twin_QMOM()
{
	give_Opposed_QMOM();
	for(int j=1;j<=2*qmom.N;j++)
		dmoments[Np][j] = moments[Np][j] - moments[Np-1][j];
}

void OpenSMOKE_Flame1D::DAESystem_Premixed_QMOM_All(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_QMOM_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)
 
		properties(data->iHOT,  DAE_Premixed_QMOM_ALL.jacobianIndex, DAE_Premixed_QMOM_ALL.jacobianVariables, NC+1 + 2*qmom.N, 1);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DT();
		give_Premixed_DW("DAE");
		give_Premixed_QMOM();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_QMOM_ALL, f);
}

void OpenSMOKE_Flame1D::DAESystem_Opposed_QMOM_ALL(BzzVector &x, double t, BzzVector &f)
{
	Tmax=T.Max();

	unsteady_boundary_conditions(t);

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(OPPOSED_QMOM_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		properties(data->iHOT, DAE_Opposed_QMOM_ALL.jacobianIndex, DAE_Opposed_QMOM_ALL.jacobianVariables, NC+4+2*qmom.N, 4);
		compute_vStarOpposed(); 
		
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DU_DG_DH();
		give_Opposed_DT(t);
		give_Opposed_DW("DAE");
		give_Opposed_QMOM();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(OPPOSED_QMOM_ALL, f);
}

void OpenSMOKE_Flame1D::DAESystem_Twin_QMOM_ALL(BzzVector &x, double t, BzzVector &f)
{
	Tmax=T.Max();

	unsteady_boundary_conditions(t);

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(TWIN_QMOM_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		properties(data->iHOT, DAE_Twin_QMOM_ALL.jacobianIndex, DAE_Twin_QMOM_ALL.jacobianVariables, NC+4+2*qmom.N, 4);
		compute_vStarOpposed(); 
		
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Twin_DU_DG_DH();
		give_Twin_DT(t);
		give_Twin_DW("DAE");
		give_Twin_QMOM();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(TWIN_QMOM_ALL, f);
}

void OpenSMOKE_Flame1D::DAESystem_Premixed_QMOM_NoEnergy(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_QMOM_NOENERGY, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_DW("DAE");
		give_Premixed_QMOM();

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_QMOM_NOENERGY, f);
}

void OpenSMOKE_Flame1D::give_Premixed_SOOT()
{
	int i;

	// 0. Moment Source Terms
	// ------------------------------------------------------------
	for(i=2;i<=Np;i++)
	{
		double DiffC = 0.;

		X.GetRow(i, &xVector);
		soot2EModel.update(T[i], data->P_atm, rho[i], DiffC, xVector, phiN[i], phiM[i]);

		soot2EModel.formation_rates();
		source_phiN[i] = soot2EModel.s;
		source_phiM[i] = soot2EModel.S;

		for (int j=1;j<=NC-1;j++)
			SootGasCorrection[i][j] = soot2EModel.SGas[j];	// Gas species
		SootGasCorrection[i][NC]	= soot2EModel.S;		// Soot

		DiffusionSoot[i] = soot2EModel.Diff;
	}

	// 1. Computation: A x lambda
	// ------------------------------------------------------------
	A_x_rho_x_Diffw[1][1] = G[1]*rho[1]*DiffusionSoot[1];
	for(i=1;i<=Ni;i++)
	{
		A_x_rho_x_Diffe[i][1]	= 0.50 * (	G[i]*rho[i]*DiffusionSoot[i] + G[i+1]*rho[i+1]*DiffusionSoot[i+1]);
		A_x_rho_x_Diffw[i+1][1] = A_x_rho_x_Diffe[i][1];
	}
	A_x_rho_x_Diffe[Np][1] = G[Np]*rho[Np]*DiffusionSoot[Np];
	
	// 2. Computation: convective terms
	// ------------------------------------------------------------
	grid.FirstDerivative('B', phiN, diff_phiN);
	grid.FirstDerivative('B', phiM, diff_phiM);

	
	// 3. Equations
	// ------------------------------------------------------------
	dphiN[1] = phiN[1] - phiNC;
	dphiM[1] = phiM[1] - phiMC;

	for(i=2;i<=Ni;i++)
	{
		dphiN[i] =  - ( H[i] * diff_phiN[i] 
		                     - (A_x_rho_x_Diffe[i][1]*(phiN[i+1]-phiN[i])*grid.udxe[i] - 
						        A_x_rho_x_Diffw[i][1]*(phiN[i]-phiN[i-1])*grid.udxw[i] ) 
						       *grid.udxc_over_2[i] )
						  /(rho[i]*G[i]) + 
						  source_phiN[i] / rho[i];

		dphiM[i] =  - ( H[i] * diff_phiM[i] 
		                     - (A_x_rho_x_Diffe[i][1]*(phiM[i+1]-phiM[i])*grid.udxe[i] - 
						        A_x_rho_x_Diffw[i][1]*(phiM[i]-phiM[i-1])*grid.udxw[i] ) 
						       *grid.udxc_over_2[i] )
						  /(rho[i]*G[i]) + 
						  source_phiM[i] / rho[i];
	}
	
	// Gradient zero
	dphiN[Np] = phiN[Np] - phiN[Np-1];
	dphiM[Np] = phiM[Np] - phiM[Np-1];
}

void OpenSMOKE_Flame1D::give_Opposed_SOOT()
{
	int i;

	// 1. Computation: A x lambda
	// ------------------------------------------------------------
	for(i=2;i<=Np;i++)
	{
		double DiffC = 1e-4;
		
		X.GetRow(i, &xVector);
		soot2EModel.update(T[i], data->P_atm, rho[i], DiffC, xVector, phiN[i], phiM[i]);	

		soot2EModel.formation_rates();
		source_phiN[i] = soot2EModel.s;
		source_phiM[i] = soot2EModel.S;
		DiffusionSoot[i] = 1e-6;
		
		for (int j=1;j<=NC-1;j++)
			SootGasCorrection[i][j] = soot2EModel.SGas[j];	// Gas species
		SootGasCorrection[i][NC]	= soot2EModel.S;		// Soot
	}
	
	// 2. Computation: convective terms
	// ------------------------------------------------------------
	grid.FirstDerivative('U', U, phiN, diff_phiN);
	grid.FirstDerivative('U', U, phiM, diff_phiM);
	
	// 3. Equations
	// ------------------------------------------------------------
	dphiN[1] = phiN[1] - phiNC;
	dphiM[1] = phiM[1] - phiMC;

	for(i=2;i<=Ni;i++)
	{
		dphiN[i] = (-(nGeometry-1.)*U[i] * diff_phiN[i] + source_phiN[i]) * urho[i] +
					 ( DiffusionSoot[i]*(phiN[i+1]-phiN[i])*grid.udxe[i] - DiffusionSoot[i-1]*(phiN[i]-phiN[i-1])*grid.udxw[i] ) *grid.udxc_over_2[i];
			
		dphiM[i] = (-(nGeometry-1.)*U[i] * diff_phiM[i] + source_phiM[i]) * urho[i] +
					 ( DiffusionSoot[i]*(phiM[i+1]-phiM[i])*grid.udxe[i] - DiffusionSoot[i-1]*(phiM[i]-phiM[i-1])*grid.udxw[i] ) *grid.udxc_over_2[i];
	}

	dphiN[Np] = phiN[Np] - phiNO;
	dphiM[Np] = phiM[Np] - phiMO;
}

void OpenSMOKE_Flame1D::give_Twin_SOOT()
{
	give_Opposed_SOOT();
	dphiN[Np] = phiN[Np] - phiN[Np-1];
	dphiM[Np] = phiM[Np] - phiM[Np-1];
}

void OpenSMOKE_Flame1D::DAESystem_Premixed_SOOT_ALL(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();
		
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_SOOT_ALL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		properties(data->iHOT, DAE_Premixed_SOOT_ALL.jacobianIndex, DAE_Premixed_SOOT_ALL.jacobianVariables, NC+1 + 2, 1);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_SOOT();
		give_Premixed_DT();
		give_Premixed_DW("DAE");

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_SOOT_ALL, f);
}

void OpenSMOKE_Flame1D::DAESystem_Premixed_SOOT_NoEnergy(BzzVector &x, double t, BzzVector &f)
{
		Tmax=T.Max();

	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(PREMIXED_SOOT_NOENERGY, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		// jacobianIndex = -2: Calcolo di tutte le proprieta(prima volta) 
		// jacobianIndex = -1: Calcolo di tutte le proprieta(valutazione del sistema)
		// jacobianIndex =  0: Calcolo di tutte le proprieta(bisogna memorizzare)

		BzzVectorInt dummy;
		properties(data->iHOT, tagConstantT, dummy, 0, 0);
		compute_vStarOpposed();
		compute_vStarPremixed();

		for(int i=1;i<=Np;i++)
			U[i] = H[i]/(G[i]*rho[i]);	// velocita [m/s]
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Premixed_SOOT();
		give_Premixed_DW("DAE");


	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(PREMIXED_SOOT_NOENERGY, f);
}

void OpenSMOKE_Flame1D::prepare_radiation()
{
	sigmaSB = 5.67e-8;									// [W/m2K4]
	Tenv4   = BzzPow4(data->environmentTemperature);	// [K4]

	if ( data->iGasRadiation == true )
	{
		iH2O = mix->recognize_species_without_exit("H2O");
		iCO2 = mix->recognize_species_without_exit("CO2");
		iCO  = mix->recognize_species_without_exit("CO");
		iCH4 = mix->recognize_species_without_exit("CH4");
	}
}

void OpenSMOKE_Flame1D::calculate_radiation()
{
	double uT;
	double K_H2O, K_CO2, K_CO, K_CH4;
	BzzVector as(4);
	double asTot;

	for(int i=2; i<=Ni; i++)
	{
		Qrad[i] = 0.;

		if ( data->iGasRadiation == true )
		{
			uT = 1000./T[i];
	
			// 1. Water [1/m.bar]
			K_H2O = -0.23093 +uT*(-1.1239+uT*(9.4153 +uT*(-2.9988 +uT*( 0.51382 + uT*(-1.8684e-5)))));
	
			// 2. Carbon Dioxide [1/m.bar]
			K_CO2 =  18.741  +uT*(-121.31+uT*(273.5  +uT*(-194.05 +uT*( 56.31   + uT*(-5.8169)))));

			// 3. Carbon monoxide [1/m.bar]
			if( T[i] < 750. )	
				K_CO = 4.7869+T[i]*(-0.06953 + T[i]*(2.95775e-4 + T[i]*(-4.25732e-7 + T[i]*2.02894e-10)));
			else
				K_CO = 10.09+T[i]*(-0.01183 + T[i]*(4.7753e-6 + T[i]*(-5.87209e-10 + T[i]*-2.5334e-14)));

			// 4. Methane [1/m.bar]
			K_CH4 = 6.6334 +T[i]*(- 0.0035686+T[i]*(1.6682e-08+T[i]*(2.5611e-10-2.6558e-14*T[i])));

			// Absorption coefficients // Secondo me devono poi essere moltiplicati per la pressione in bar
			if(iH2O!=0)	as[1] = K_H2O*X[i][iH2O];		// [1/m]
			if(iCO2!=0)	as[2] = K_CO2*X[i][iCO2];		// [1/m]
			if(iCO!=0)	as[3] = K_CO*X[i][iCO];			// [1/m]
			if(iCH4!=0)	as[4] = K_CH4*X[i][iCH4];		// [1/m]

			asTot = as.GetSumElements();				// Absorption Coefficient


			// Source term
			Qrad[i] += - 4.*sigmaSB * asTot * (BzzPow4(T[i]) - Tenv4); // Source term [W/m3]
		}

		if ( mix->polimiSoot->IsSoot() && data->iRadiativeSootModel != RADIATIVE_SOOT_MODEL_NONE )
		{
			double as_Soot, fv_Soot;

			X.GetRow(i, &cVector);
			double cTotForSoot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			cVector *= cTotForSoot;
			
			double MWgas=0.;
			for(int k=1;k<=NC;k++)
				MWgas += X[i][k]*mix->M[k];
			double rhoGas = cTotForSoot*MWgas;
		

			W.GetRow(i, &wVector);
			X.GetRow(i, &xVector);
			mix->polimiSoot->Analysis(*mix, data->P_Pascal, T[i], wVector);
			fv_Soot = mix->polimiSoot->fv_large();

			// a. FLUENT Model
			if (data->iRadiativeSootModel == RADIATIVE_SOOT_MODEL_FLUENT)
				as_Soot = 1232.4 * (fv_Soot*1927.5946) *( 1. + 4.8e-4*(T[i]-2000.) );	// [1/m]

			// b. Bilger Model
			else if (data->iRadiativeSootModel == RADIATIVE_SOOT_MODEL_BILGER)
				as_Soot = 1307. * fv_Soot * T[i];										// [1/m]

			// c. Widmann Model
			else if (data->iRadiativeSootModel == RADIATIVE_SOOT_MODEL_WIDMANN)
				as_Soot = 2370. * fv_Soot * T[i];										// [1/m]

			// Source term
			Qrad[i] -= 4.*sigmaSB * as_Soot * (BzzPow4(T[i]) - Tenv4);					// Source term [W/m3]
		}
	}
}

void OpenSMOKE_Flame1D::final_sensitivity_analysis(std::string kindOfProblem,  BzzNonLinearSystemSparseObject &nls, BzzSave &fBinary)
{
	OpenSMOKE_SensitivityAnalysis_Fast_Flame1D sensitivity;
	int kind_of_flame;

	     if (kindOfProblem == "PREMIXED_FLAMESPEED")	kind_of_flame			= 1;
	else if (kindOfProblem == "PREMIXED_ALL")			kind_of_flame			= 2;
	else if (kindOfProblem == "OPPOSED_ALL")			kind_of_flame			= 3;
	else if (kindOfProblem == "TWIN_ALL")				kind_of_flame			= 3;
	else if (kindOfProblem == "PREMIXED_NOENERGY")		kind_of_flame			= 4;
	else ErrorMessage("Wrong kind of sensitivity analysis...");

	cout << "------------------------------------------------------------------"	<< endl;
	cout << "     Final Sensitivity Analysis: " << kindOfProblem					<< endl;
	cout << "------------------------------------------------------------------"	<< endl;

	cout << "1. Initialize sensitivity..." << endl;
	cout << kindOfProblem << " " << endl;
	sensitivity.Initialize(kind_of_flame, mix, data->index_sensitivity, Np, kindOfSensitivity);
	sensitivity.VideoSummary();

	cout << "2. Sensitivity Analysis..." << endl;
	sensitivity.AFactorized = &nls.JTrid;
	cout << "1" << endl;
	if (data->iAssignedFixedTemperatureProfileProvisional==true)
	{
		BzzVector CpProvisional(Np); CpProvisional=-1.;
		sensitivity.BuildJAlfaMatrix(T, rho, CpProvisional, data->names_sensitivity, data->P_Pascal, X);
	}
	else
	{
		sensitivity.BuildJAlfaMatrix(T, rho, Cp, data->names_sensitivity, data->P_Pascal, X);
	}
	cout << "2" << endl;
	sensitivity.SaveOnBinaryFile(fBinary);
	cout << "3" << endl;
	sensitivity.PrintOnFile(nameFolderAdditionalData, grid.x, T, H, W, U, G, rho, PMtot, data->names_sensitivity);
	cout << "4" << endl;

	std::string nameXMLFolder = nameFolderAdditionalData;
	sensitivity.SaveOnXMLFile(nameXMLFolder);
}

void OpenSMOKE_Flame1D::final_sensitivity_analysis_transport_properties(std::string kindOfProblem,  BzzNonLinearSystemSparseObject &nls, BzzVector &xSolution,  BzzVector &fSolution, BzzSave &fBinary)
{
	OpenSMOKE_LennardJonesSensitivityCoefficientsManager lennard_jones_manager;
	lennard_jones_manager.Setup(mix->path_complete + "/TransportSensitivity.out", data->sensitivityLennardJonesMode);

	int kind_of_flame;
	     if (kindOfProblem == "PREMIXED_FLAMESPEED")	kind_of_flame			= 1;
	else if (kindOfProblem == "PREMIXED_ALL")			kind_of_flame			= 2;
	else if (kindOfProblem == "OPPOSED_ALL")			kind_of_flame			= 3;
	else if (kindOfProblem == "TWIN_ALL")				kind_of_flame			= 3;
	else if (kindOfProblem == "PREMIXED_NOENERGY")		kind_of_flame			= 4;
	else ErrorMessage("Wrong kind of sensitivity analysis...");

	cout << "------------------------------------------------------------------"	<< endl;
	cout << "  Final Sensitivity Analysis  (Transport Properties):             "    << kindOfProblem << endl;
	cout << "------------------------------------------------------------------"	<< endl;

//	BzzVectorInt indices;
//	sensitivity.Initialize(kind_of_flame, mix, indices, Np);
//	sensitivity.VideoSummary();

//	cout << "0. Recovering Jacobian..." << endl;
//	sensitivity.AFactorized = &nls.JTrid;

	cout << "1. Assembling Jalfa..." << endl;
	int NPtransport     = 5;
	int dimEquations	= fSolution.Size();
	int dimBlock		= fSolution.Size()/Np;

	BzzVector fPlus(fSolution.Size());
	BzzVector fMinus(fSolution.Size());
	BzzMatrix JAlfa(dimEquations, NC*NPtransport);

	for(int kk=1;kk<=NPtransport;kk++)
	{
		BzzVector denominator;
		fittingCoefficientExtraction fittingPlus;
		fittingCoefficientExtraction fittingMinus;

		if (kk==1)	{	fittingPlus  = OPENSMOKE_FITTING_KE_PLUS;		fittingMinus = OPENSMOKE_FITTING_KE_MINUS;		denominator = lennard_jones_manager.delta_epsylon_over_kb; }
		if (kk==2)	{	fittingPlus  = OPENSMOKE_FITTING_SIGMA_PLUS;	fittingMinus = OPENSMOKE_FITTING_SIGMA_MINUS;	denominator = lennard_jones_manager.delta_sigma; }
		if (kk==3)	{	fittingPlus  = OPENSMOKE_FITTING_MU_PLUS;		fittingMinus = OPENSMOKE_FITTING_MU_MINUS;		denominator = lennard_jones_manager.delta_mu; }
		if (kk==4)	{	fittingPlus  = OPENSMOKE_FITTING_ALFA_PLUS;		fittingMinus = OPENSMOKE_FITTING_ALFA_MINUS;	denominator = lennard_jones_manager.delta_alfa; }
		if (kk==5)	{	fittingPlus  = OPENSMOKE_FITTING_ZROT_PLUS;		fittingMinus = OPENSMOKE_FITTING_ZROT_MINUS;	denominator = lennard_jones_manager.delta_zRot298; }

		for(int j=1;j<=NC;j++)
		{
			cout << "Species: " << j << " " << dimEquations << " " << dimBlock << endl;

			cout << "Before " << fSolution[dimEquations/2] << " " <<  mu[Np/2] << endl;

			// Positive increment
			{
				lennard_jones_manager.GetFittingCoefficients(fittingPlus, j, *mix);
				if (kindOfProblem == "PREMIXED_FLAMESPEED")		nonLinearSystem_Premixed_FlameSpeed(xSolution, fPlus);
				if (kindOfProblem == "PREMIXED_ALL")			nonLinearSystem_Premixed_All(xSolution, fPlus);
				if (kindOfProblem == "OPPOSED_ALL")				nonLinearSystem_Opposed_All(xSolution, fPlus);
				if (kindOfProblem == "TWIN_ALL")				nonLinearSystem_Twin_All(xSolution, fPlus);
				if (kindOfProblem == "PREMIXED_NOENERGY")		nonLinearSystem_Premixed_NoEnergy(xSolution, fPlus);
			}

			cout << "Plus   " << fPlus[dimEquations/2] << " " << mu[Np/2] << endl;
			
			// Negative increment
			{
				lennard_jones_manager.GetFittingCoefficients(fittingMinus, j, *mix);
				if (kindOfProblem == "PREMIXED_FLAMESPEED")		nonLinearSystem_Premixed_FlameSpeed(xSolution, fMinus);
				if (kindOfProblem == "PREMIXED_ALL")			nonLinearSystem_Premixed_All(xSolution, fMinus);
				if (kindOfProblem == "OPPOSED_ALL")				nonLinearSystem_Opposed_All(xSolution, fMinus);
				if (kindOfProblem == "TWIN_ALL")				nonLinearSystem_Twin_All(xSolution, fMinus);
				if (kindOfProblem == "PREMIXED_NOENERGY")		nonLinearSystem_Premixed_NoEnergy(xSolution, fMinus);
			}
			cout << "Minus   " << fMinus[dimEquations/2] << " " << mu[Np/2] << endl;
			
			
			for(int i=1;i<=Np;i++)
				for(int k=1;k<=dimBlock;k++)
				{	
					int index  = (i-1)*dimBlock+k;
					JAlfa[index][j+(kk-1)*NC] = (fPlus[index]-fMinus[index]) / denominator[j];
				}

			// Go back to original values
			lennard_jones_manager.GetFittingCoefficients(OPENSMOKE_FITTING_BASE, j, *mix);
		}
	}

	cout << "2. Sensitivity Coefficients..." << endl;
	BzzMatrix S(dimEquations, NC);
	JAlfa *= -1.;
	Solve(&nls.JTrid, JAlfa, &S);

	cout << "3. Save on binary file..." << endl;
	{
		BzzVectorInt indices_print_species(NC);
		for(int j=1;j<=NC;j++)
			indices_print_species[j] = j;

		BzzVector parameters(NC*NPtransport);
		for(int kk=1;kk<=NPtransport;kk++)
		{
			BzzVector numerator;

			if (kk==1)	{	numerator = lennard_jones_manager.epsylon_over_kb; }
			if (kk==2)	{	numerator = lennard_jones_manager.sigma; }
			if (kk==3)	{	numerator = lennard_jones_manager.mu; }
			if (kk==4)	{	numerator = lennard_jones_manager.alfa; }
			if (kk==5)	{	numerator = lennard_jones_manager.zRot298; }

			for(int j=1;j<=NC;j++)
				parameters[(kk-1)*NC+j] = numerator[j];
		}

		std::string dummy;
		char name[Constants::NAME_SIZE];

		dummy = "SENSITIVITY-DIFF";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));

		dummy = "V20100417";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));

		if (kindOfProblem == "PREMIXED_ALL")		dummy = "PREMIXED";
		if (kindOfProblem == "PREMIXED_FLAMESPEED")	dummy = "FLAMESPEED";
		if (kindOfProblem == "OPPOSED_ALL")			dummy = "OPPOSED";
		if (kindOfProblem == "TWIN_ALL")			dummy = "TWIN";
		if (kindOfProblem == "PREMIXED_NOENERGY")	dummy = "PREMIXED";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));

		dummy = "INDICES";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));
		fBinary << indices_print_species;

		dummy = "NP";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));
		fBinary << NC*NPtransport;

		dummy = "NV";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));
		fBinary << dimBlock;

		dummy = "N";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));
		fBinary << Np;

		dummy = "PARAMETERS";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));
		fBinary << parameters;

		dummy = "S_S";
		strcpy(name, dummy.c_str());
		fBinary.fileSave.write((char*) name, sizeof(name));
		fBinary << S;
	}

//	cout << "4. Print on file..." << endl;
//	OpenSMOKE_SensitivityAnalysis_Diffusion_Flame1D post_processing;
//	post_processing.Initialize(kind_of_flame, mix, data->index_sensitivity, Np);
//	post_processing.PrintOnFile(grid.x, T, H, W, U, G, rho, PMtot, data->names_sensitivity, S, Diffusion);
}

void OpenSMOKE_Flame1D::FoldersAndFilesManager()
{
	std::string MSDOScommand;

	if (data->iUserDefinedFolderName == false)
	{
		// 1. Output Folder name
		if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED)
		{
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)	nameOutputFolder = "Output_Opposed_";
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)	nameOutputFolder = "Output_Opposed_QMOM_";
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)	nameOutputFolder = "Output_Opposed_SOOT_";
		}

		if (data->kind_of_flame == FLAME1D_PHYSICS_TWIN)
		{
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)	nameOutputFolder = "Output_Twin_";
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)	nameOutputFolder = "Output_Twin_QMOM_";
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)	nameOutputFolder = "Output_Twin_SOOT_";
		}

		if (data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED)
		{
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)	nameOutputFolder = "Output_Premixed_";
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)	nameOutputFolder = "Output_Premixed_QMOM_";
			if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)	nameOutputFolder = "Output_Premixed_SOOT_";
		}
		nameOutputFolder += mix->name_kinetic_scheme;
	}
	else
		nameOutputFolder = data->userDefinedFolderName;

	MSDOScommand = "mkdir " + nameOutputFolder;
	system(MSDOScommand.c_str());

	// 2. Backup Input data file name
	nameFolderBackupData		= nameOutputFolder;
	nameFolderUnsteadyData		= nameOutputFolder;
	nameFolderSteadyData		= nameOutputFolder;
	nameFolderAdditionalData	= nameOutputFolder;
	#if LINUX_SO==1
		nameFolderBackupData     += "/BackUp";
		nameFolderUnsteadyData   += "/Unsteady";
		nameFolderSteadyData     += "/Steady";
		nameFolderAdditionalData += "/Additional";
	#else
		nameFolderBackupData     += "\\BackUp";
		nameFolderUnsteadyData   += "\\Unsteady";
		nameFolderSteadyData     += "\\Steady";
		nameFolderAdditionalData += "\\Additional";
	#endif

	MSDOScommand = "mkdir " + nameFolderBackupData;
	system(MSDOScommand.c_str());

	MSDOScommand = "mkdir " + nameFolderUnsteadyData;
	system(MSDOScommand.c_str());

	MSDOScommand = "mkdir " + nameFolderSteadyData;
	system(MSDOScommand.c_str());

	MSDOScommand = "mkdir " + nameFolderAdditionalData;
	system(MSDOScommand.c_str());
	
	nameFileBackupInputData = nameFolderBackupData + "/BackUp";
}

void OpenSMOKE_Flame1D::FoldersAndFilesManager(const std::string backupFolder)
{
	// 1. Output Folder name
	nameOutputFolder = backupFolder;

	// 2. BackUp Folder name
	nameFolderBackupData		=	nameOutputFolder;
	nameFolderUnsteadyData		=	nameOutputFolder;
	nameFolderSteadyData		=	nameOutputFolder;
	nameFolderAdditionalData	=	nameOutputFolder;

	#if LINUX_SO==1
		nameFolderBackupData		+= "/BackUp";
		nameFileBackupInputData		 =	nameFolderBackupData;
		nameFileBackupInputData		+= "/BackUp";
		nameFolderUnsteadyData		+= "/Unsteady";
		nameFolderSteadyData		+= "/Steady";
		nameFolderAdditionalData	+= "/Additional";
	#else
		nameFolderBackupData		+= "\\BackUp";
		nameFileBackupInputData		 = nameFolderBackupData;
		nameFileBackupInputData		+= "\\BackUp";
		nameFolderUnsteadyData		+= "\\Unsteady";
		nameFolderSteadyData		+= "\\Steady";
		nameFolderAdditionalData	+= "\\Additional";
	#endif
}

void OpenSMOKE_Flame1D::Run()
{
	const int Hot  = 1;
	const int Cold = 0;

	std::string nameFileSolution;
	std::string nameFileSolutionComplete;

	nameFileSolution = nameFolderSteadyData + "/Solution_";


	if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED)
	{
		OpenSMOKE_Flame1D_OpposedFlameManager opposed_flame_manager;

		if (data->iOpposedFlameAnalysis == true)
			opposed_flame_manager.SetupFromFile(this, data);

		if (data->iBackUp == false)
		{
			// COLD solution (0) - NoSpark (2)
			// ---------------------------------------------------------------
			initializeVariables(FLAME1D_IGNITION_NO_SPARK);
		}

		//for(int indexOperation=1;indexOperation<=operations->nOperations;indexOperation++)
		int indexOperation = 1;
		while (indexOperation <= operations->nOperations)
		{
			bool iCalculate = true;
			bool iSensitivityNow = false;

			// NLS HOT SOLUTION (1)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==1)
			{
				if (operations->iOptionA[indexOperation] == 1) solveNLS_Opposed(Hot, OPPOSED_ONLY_MOMENTUM);
				if (operations->iOptionA[indexOperation] == 2) solveNLS_Opposed(Hot, OPPOSED_ALL);
				if (operations->iOptionA[indexOperation] == 3) solveNLS_Opposed(Hot, OPPOSED_ONLY_MASS_FRACTIONS);
			}

			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==111)
			{
				if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED) nonLinearSystemSolutionWithControl_Opposed(operations->iOptionB[indexOperation]);
				if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED)	nonLinearSystemSolutionWithControl_Twin(operations->iOptionB[indexOperation]);
			}

			// SENSITIVITY ANALYSIS
			if (operations->iOperation[indexOperation]==6)
			{
				if (operations->iOptionA[indexOperation] == 1) kindOfSensitivity = FREQUENCY_FACTOR;
				if (operations->iOptionA[indexOperation] == 2) kindOfSensitivity = TRANSPORT_PROPERTIES;
				if (operations->iOptionA[indexOperation] == 3) kindOfSensitivity = FREQUENCY_FACTOR_AND_TRANSPORT_PROPERTIES;
				iSensitivityNow = true;
				//solveNLS_Opposed(Hot, OPPOSED_ALL);
			}

			// DAE HOT SOLUTION (3)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==3)
			{
				if (iUnsteadyFromBackUp!=-1)
				{	
				//	solveNLS_Opposed(Hot, OPPOSED_ONLY_MOMENTUM);
				//	if (data->iUnsteady != true)	solveDAE_Opposed(Hot, OPPOSED_ONLY_MOMENTUM, 1e5);
					solveNLS_Opposed(Hot, OPPOSED_ONLY_MOMENTUM);
				}

				if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)
				{
					solveNLS_Reduced_Opposed(Hot, OPPOSED_ALL);
					solveDAE_Opposed(Hot, OPPOSED_ALL, operations->iOptionB[indexOperation]);
				}
				else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)
					solveDAE_Opposed(Hot, OPPOSED_SOOT_ALL, operations->iOptionB[indexOperation]);
				else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)
					solveDAE_Opposed(Hot, OPPOSED_QMOM_ALL, operations->iOptionB[indexOperation]);
			}

			if (operations->iOperation[indexOperation]==61)
				SolveODE_SingleReactors(operations->iOptionB[indexOperation]);

			// NEW POINTS (40)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==41)
			{
				iCalculate = newPoints("TEMPERATURE", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==42)
			{
				iCalculate = newPoints("QREACTION", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==43)
			{
				iCalculate = newPoints("ALL", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==44)
			{
				char option[200];
				strcpy(option, operations->species[indexOperation].c_str());
				iCalculate = newPoints(option, operations->iOptionA[indexOperation]);
			}
			
			// DOUBLE - REFINE the grid (50)
			// ---------------------------------------------------------------
			if (operations->iOperation[indexOperation]==51)
			{
				doubleTheGrid();
			}

			if (operations->iOperation[indexOperation]==52)
			{
				refineGridPeak(0.80);
			}

			if (operations->iOperation[indexOperation]==54)
			{
				refineFlameBase(operations->iOptionA[indexOperation], operations->iOptionB[indexOperation]);
			}

			if (operations->iOperation[indexOperation]==55)
			{
				adaptGrid(operations->iOptionA[indexOperation]);
			}

			if (operations->iOperation[indexOperation]==58)
			{
				refineStagnationPlane(operations->iOptionA[indexOperation]);
			}

			// ------------------------------------------------------------------
			// Stampa su file
			// ------------------------------------------------------------------
			if (iCalculate == true)
			{
				stringstream operation;
				operation << indexOperation;
				nameFileSolutionComplete = nameFileSolution + operation.str() + ".out";
				printOnFile(nameFileSolutionComplete);

				#if LINUX_SO==1
					nameFileSolutionComplete = nameFolderAdditionalData + "/opposed.osm";
					std::string instruction = "cp " + mix->path_complete + "idealgas.bin" + " nameFolderAdditionalData/idealgas.bin";
					system(instruction.c_str());
				#else
					nameFileSolutionComplete = nameFolderAdditionalData + "\\opposed.osm";
					std::string instruction;
					instruction = "copy " + mix->path_complete + "\\idealgas.bin " + nameFolderAdditionalData + "\\idealgas.bin /Y";
					system(instruction.c_str());
					instruction = "copy " + mix->path_complete + "\\reactions.bin " + nameFolderAdditionalData + "\\reactions.bin /Y";
					system(instruction.c_str());
				#endif
				SaveOnBinaryFile(nameFileSolutionComplete, iSensitivityNow, OPPOSED_ALL);
			}
			else
				indexOperation++;

			indexOperation++;
		}

		if (data->iOpposedFlameAnalysis == true)
			opposed_flame_manager.Run();
	}

	else if (data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED)
	{
		flame1d_model KIND;
		OpenSMOKE_Flame1D_FlameSpeedManager flame_speed_manager;


		if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)
		{
			if (data->iTemperatureProfile==0 || data->iTemperatureProfile==2)
				KIND = PREMIXED_ALL;
			if (data->iTemperatureProfile==1)
				KIND = PREMIXED_NOENERGY;
			if (data->iFixedTemperature>0)
				KIND = PREMIXED_FLAMESPEED;
		}

		if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)
		{
			if (data->iTemperatureProfile==0 || data->iTemperatureProfile==2)
				KIND = PREMIXED_QMOM_ALL;
			if (data->iTemperatureProfile==1)
				KIND = PREMIXED_QMOM_NOENERGY;
		}
		
		if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)
		{
			if (data->iTemperatureProfile==0 || data->iTemperatureProfile==2)
				KIND = PREMIXED_SOOT_ALL;
			if (data->iTemperatureProfile==1)
				KIND = PREMIXED_SOOT_NOENERGY;
		}

		if (data->iBackUp == false)
			initializeVariables(FLAME1D_IGNITION_NONE);
	
		if (data->iFlameSpeedAnalysis == true)
			flame_speed_manager.SetupFromFile(this, data);

		//for(int indexOperation=1;indexOperation<=operations->nOperations;indexOperation++)
		int indexOperation = 1;
		while (indexOperation <= operations->nOperations)
		{
			bool iCalculate = true;
			bool iSensitivityNow = false;

			// ------------------------------------------------------------------------
			
			// For premixed problems we have three different kinds of calculations
			// 1. Premixed flame with mass flow rate assigned
			// 2. Premixed flame with mass flow rate and temperature assigned
			// 3. Primixed flame with mass flow rate unknown
			

			// ------------------------------------------------------------------
			// START CALCULATIONS
			// ------------------------------------------------------------------

			// NLS SOLUTION
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==1)
			{
				solveNLS_Premixed(KIND);
			}

			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==111)
			{
				nonLinearSystemSolutionWithControl_FlameSpeed(operations->iOptionB[indexOperation]);
			}

			// SENSITIVITY (6)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==6)
			{
				if (operations->iOptionA[indexOperation] == 1) kindOfSensitivity = FREQUENCY_FACTOR;
				if (operations->iOptionA[indexOperation] == 2) kindOfSensitivity = TRANSPORT_PROPERTIES;
				if (operations->iOptionA[indexOperation] == 3) kindOfSensitivity = FREQUENCY_FACTOR_AND_TRANSPORT_PROPERTIES;
				iSensitivityNow = true;
			//	solveNLS_Premixed(KIND);
			}
			
			if (operations->iOperation[indexOperation]==4)
			{
				ErrorMessage("Premixed Flame Speed Operation 4 not allowed");
				solveNLS_Premixed(PREMIXED_FLAMESPEED);
			}

			// ODE SOLUTION
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==2)
			{
				solveODE_Premixed(KIND, operations->iOptionB[indexOperation]);
			}

			// DAE SOLUTION
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==3)
			{
				solveDAE_Premixed(KIND, operations->iOptionB[indexOperation]);
			}

			// NEW POINTS (40)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==41)
			{
				iCalculate = newPoints("TEMPERATURE", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==42)
			{
				iCalculate = newPoints("QREACTION", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==43)
			{
				iCalculate = newPoints("ALL", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==44)
			{
				char option[200];
				strcpy(option, operations->species[indexOperation].c_str());
				iCalculate = newPoints(option, operations->iOptionA[indexOperation]);
			}


			// REFINE (50)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==51)
			{
				doubleTheGrid();
			}
			if (operations->iOperation[indexOperation]==53)
			{
				refineFlameBase();
			}
			if (operations->iOperation[indexOperation]==54)
			{
				refineFlameBase(operations->iOptionA[indexOperation], operations->iOptionB[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==55)
			{
				adaptGrid(operations->iOptionA[indexOperation]);
			}

			if (operations->iOperation[indexOperation]==56)
			{
				refineBase();
			}

			if (operations->iOperation[indexOperation]==57)
			{
				refineAttachPoint();
			}

			if (operations->iOperation[indexOperation]==58)
			{
				refineStagnationPlane(operations->iOptionA[indexOperation]);
			}

			// ------------------------------------------------------------------
			// Stampa su file
			// ------------------------------------------------------------------
			if (iCalculate == true)
			{
				stringstream operation;
				operation << indexOperation;
				nameFileSolutionComplete = nameFileSolution + operation.str() + ".out";
				printOnFile(nameFileSolutionComplete);

				#if LINUX_SO==1
				nameFileSolutionComplete = nameFolderAdditionalData + "/premixed.osm";
				#else
					nameFileSolutionComplete = nameFolderAdditionalData + "\\premixed.osm";
					std::string instruction;
					instruction = "copy " + mix->path_complete + "\\idealgas.bin " + nameFolderAdditionalData + "\\idealgas.bin /Y";
					system(instruction.c_str());
					instruction = "copy " + mix->path_complete + "\\reactions.bin " + nameFolderAdditionalData + "\\reactions.bin /Y";
					system(instruction.c_str());
				#endif
				SaveOnBinaryFile(nameFileSolutionComplete, iSensitivityNow, KIND);
			}
			else
				indexOperation++;
			indexOperation++;
		}

		if (data->iFlameSpeedAnalysis == true)
			flame_speed_manager.Run();
	}

	else if (data->kind_of_flame == FLAME1D_PHYSICS_TWIN)
	{
		if (data->iBackUp == false)
		{
			// COLD solution (0) - NoSpark (2)
			// ---------------------------------------------------------------
			initializeVariables(FLAME1D_IGNITION_NO_SPARK);
		}

		int indexOperation = 1;
		while (indexOperation <= operations->nOperations)
		{
			bool iCalculate = true;
			bool iSensitivityNow = false;

			// NLS HOT SOLUTION (1)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==1)
			{
				if (operations->iOptionA[indexOperation] == 1) solveNLS_Twin(Hot, TWIN_ONLY_MOMENTUM);
				if (operations->iOptionA[indexOperation] == 2) solveNLS_Twin(Hot, TWIN_ALL);
				if (operations->iOptionA[indexOperation] == 3) solveNLS_Twin(Hot, TWIN_ONLY_MASS_FRACTIONS);
			}


			// SENSITIVITY (6)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==6)
			{
				if (operations->iOptionA[indexOperation] == 1) kindOfSensitivity = FREQUENCY_FACTOR;
				if (operations->iOptionA[indexOperation] == 2) kindOfSensitivity = TRANSPORT_PROPERTIES;
				if (operations->iOptionA[indexOperation] == 3) kindOfSensitivity = FREQUENCY_FACTOR_AND_TRANSPORT_PROPERTIES;
				iSensitivityNow = true;
			//	solveNLS_Opposed(Hot, TWIN_ALL);
			}

			// DAE HOT SOLUTION (3)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==3)
			{
				if (iUnsteadyFromBackUp!=-1)
					solveNLS_Twin(Hot, TWIN_ONLY_MOMENTUM);

				if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_NONE)
				{
					//solveDAE_Twin(Hot, "TWIN_NO_ENERGY", operations->iOptionB[indexOperation]);
					// TODO
					solveNLS_Reduced_Twin(Hot, TWIN_ALL);
					solveDAE_Twin(Hot, TWIN_ALL, operations->iOptionB[indexOperation]);
				}
				else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_SOOT)
					solveDAE_Twin(Hot, TWIN_SOOT_ALL, operations->iOptionB[indexOperation]);
				else if (data->kind_of_subphysics == FLAME1D_SUBPHYSICS_QMOM)
					solveDAE_Twin(Hot, TWIN_QMOM_ALL, operations->iOptionB[indexOperation]);
			}

			// NEW POINTS (40)
			// ------------------------------------------------------------------
			if (operations->iOperation[indexOperation]==41)
			{
				iCalculate = newPoints("TEMPERATURE", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==42)
			{
				iCalculate = newPoints("QREACTION", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==43)
			{
				iCalculate = newPoints("ALL", operations->iOptionA[indexOperation]);
			}
			if (operations->iOperation[indexOperation]==44)
			{
				char option[200];
				strcpy(option, operations->species[indexOperation].c_str());
				iCalculate = newPoints(option, operations->iOptionA[indexOperation]);
			}
			
			// DOUBLE - REFINE the grid (50)
			// ---------------------------------------------------------------
			if (operations->iOperation[indexOperation]==51)
			{
				doubleTheGrid();
			}

			if (operations->iOperation[indexOperation]==52)
			{
				refineGridPeak(0.80);
			}

			if (operations->iOperation[indexOperation]==54)
			{
				refineFlameBase(operations->iOptionA[indexOperation], operations->iOptionB[indexOperation]);
			}

			if (operations->iOperation[indexOperation]==55)
			{
				adaptGrid(operations->iOptionA[indexOperation]);
			}

			if (operations->iOperation[indexOperation]==58)
			{
				refineStagnationPlane(operations->iOptionA[indexOperation]);
			}

			// ------------------------------------------------------------------
			// Stampa su file
			// ------------------------------------------------------------------
			if (iCalculate == true)
			{
				stringstream operation;
				operation << indexOperation;
				nameFileSolutionComplete = nameFileSolution + operation.str() + ".out";
				printOnFile(nameFileSolutionComplete);

				#if LINUX_SO==1
				nameFileSolutionComplete = nameFolderAdditionalData + "/twin.osm";
				#else
					nameFileSolutionComplete = nameFolderAdditionalData + "\\twin.osm";
					std::string instruction;
					instruction = "copy " + mix->path_complete + "\\idealgas.bin " + nameFolderAdditionalData + "\\idealgas.bin /Y";
					system(instruction.c_str());
					instruction = "copy " + mix->path_complete + "\\reactions.bin " + nameFolderAdditionalData + "\\reactions.bin /Y";
					system(instruction.c_str());
				#endif

				SaveOnBinaryFile(nameFileSolutionComplete, iSensitivityNow, TWIN_ALL);
			}
			else
				indexOperation++;

			indexOperation++;
		}
	}

	// ------------------------------------------------------------------
	// Stampa su file
	// ------------------------------------------------------------------
	
	nameFileSolutionComplete = nameFileSolution + "Final.out";
	printOnFile(nameFileSolutionComplete);

	nameFileSolution = nameFolderBackupData + "/BackUp_FINAL";
	printBackUpOnlyInputData(nameFileSolution);

	nameFileSolution = nameFolderBackupData + "/BackUp_FINAL";
	printBackUpOnlyData(nameFileSolution);

	std::string nameXMLFile = nameFolderAdditionalData + "\\Output.xml";
	PrintXMLFile(nameXMLFile);
}

// Print XML Files
void OpenSMOKE_Flame1D::PrintXMLFile(const std::string file_name)
{
	unsigned int n_additional = 7;
	if ( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED )
		n_additional++;
	ofstream fXML;
	fXML.open(file_name.c_str(), ios::out);
	fXML.setf(ios::scientific);

	fXML << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
	fXML << "<opensmoke version=\"0.1a\">" << endl;
		
	fXML << "<Type> Flame1D </Type>" << endl;
		
	fXML << "<additional>" << endl;
	fXML << n_additional << endl;
	fXML << "axial-coordinate [cm] 2" << endl;
	fXML << "temperature [K] 3" << endl;
	fXML << "pressure [Pa] 4" << endl;
	fXML << "mol-weight [kg/kmol] 5" << endl;
	fXML << "density [kg/m3] 6" << endl;
	fXML << "heat-release [W/m3] 7" << endl;
	fXML << "axial-velocity [m/s] 8" << endl;
	if ( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED )
		fXML << "mass-flow-rate [kg/m2/s] 9" << endl;
	fXML << "</additional>" << endl;
		
	fXML << "<t-p-mw>" << endl;
	fXML << "1 2 3" << endl;
	fXML << "</t-p-mw>" << endl;
		
	fXML << "<mass-fractions>" << endl;
	fXML << NC << endl;
	for(int i=1;i<=NC;i++)
		fXML << mix->names[i] << " " << mix->M[i] << " " << n_additional+i << std::endl;
	fXML << "</mass-fractions>" << endl;

	fXML << "<profiles>" << endl;
	for(int i=1;i<=Np;i++)
	{
		fXML << 1.e2*grid.x[i] << " ";
		fXML << T[i] << " ";
		fXML << data->P_Pascal << " ";
		fXML << PMtot[i] << " ";
		fXML << rho[i] << " ";
		fXML << QReaction[i] << " ";

		if ( data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED || data->kind_of_flame == FLAME1D_PHYSICS_TWIN )
			fXML << (nGeometry-1.)*U[i] / rho[i] << " ";
		else if ( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED )
			fXML << U[i] << " ";

		if ( data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED )
			fXML << H[i] << " ";

		for(int j=1;j<=NC;j++)
			fXML << W[i][j] << " ";
		fXML << std::endl;
	}
	fXML << "</profiles>" << endl;

	fXML << "<profiles-size>" << endl;
	fXML << Np << " " << NC+n_additional << std::endl;
	fXML << "</profiles-size>" << endl;
	fXML << "</opensmoke>" << endl;
	fXML.close();
}

void OpenSMOKE_Flame1D::ReRun()
{
	data->iBackUp = true;

	Run();
}

void OpenSMOKE_Flame1D::PasteFromExternalSolution(OpenSMOKE_Flame1D_Solution &solution)
{
	ptFlame = this;
	
	if (Np != solution.Np)
	{	
		grid.Construct(solution.coordinates);
		Np = grid.Np;
		Ni = grid.Ni;
		allocate_only_Np();
	
	//	ErrorMessage("External Solution must have the same dimension of destination solution");
	}

	U = solution.U;
	G = solution.G;
	H = solution.H;
	T = solution.T;
	W = solution.W;
		
	if (data->iQMOM == true)
		moments = solution.moments;
			
	if (data->i2E == true)
	{
		phiN = solution.phiN;
		phiM = solution.phiM;
	}

	molarFractionsAndPMtot();

	// Data manager
	data->VC = solution.VC;
	data->VO = solution.VO;
	data->xO = solution.xO;
	data->xC = solution.xC;
	data->XO = solution.XO;
	data->XC = solution.XC;
	data->MassFlowRate = solution.MassFlowRate;
	data->CrossSectionalArea = solution.CrossSectionalArea;
	data->iFixedTemperature = solution.iFixedTemperature;
	data->fixedTemperature = solution.fixedTemperature;
}

void OpenSMOKE_Flame1D_Solution::PasteFromExternalSolution(OpenSMOKE_Flame1D &flame, OpenSMOKE_Flame1D_DataManager &data_manager)
{
	Np = flame.Np;
	x  = flame.grid.x;
	U  = flame.U;
	G  = flame.G;
	H  = flame.H;
	T  = flame.T;
	W  = flame.W;
	if (flame.data->iQMOM == true)
		moments = flame.moments;
			
	if (flame.data->i2E == true)
	{
		phiN = flame.phiN;
		phiM = flame.phiM;
	}

	// Data manager
	VC = data_manager.VC;
	VO = data_manager.VO;
	xO = data_manager.xO;
	xC = data_manager.xC;
	XO = data_manager.XO;
	XC = data_manager.XC;
	MassFlowRate = data_manager.MassFlowRate;
	CrossSectionalArea = data_manager.CrossSectionalArea;
	iFixedTemperature = data_manager.iFixedTemperature;
	fixedTemperature = data_manager.fixedTemperature;

	coordinates = flame.grid.x;
}

void OpenSMOKE_Flame1D::ElementalAnalysis()
{
	BzzVector omega_elemental_vector;
	BzzVector omega_elemental_fuel_vector;
	BzzVector omega_elemental_air_vector;
	
	// Elemental analysis
	mix->GetElementalMoleFractionsFromSpeciesMoleFractions(x_elemental, X);
	mix->GetElementalMassFractionsFromSpeciesMassFractions(omega_elemental, W);

	for(int i=1;i<=Np;i++)
	{
		omega_elemental_vector = omega_elemental.GetRow(i);
		omega_elemental_fuel_vector = omega_elemental.GetRow(1);
		omega_elemental_air_vector = omega_elemental.GetRow(Np);
		Z[i] = mix->GetMixtureFraction(omega_elemental_vector, omega_elemental_fuel_vector, omega_elemental_air_vector);
	}
}

void OpenSMOKE_Flame1D::SaveOnBinaryFile(BzzSave &fOutput)
{
	std::string dummy;
	char name[Constants::NAME_SIZE];

	dummy = "FLAME1D";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "V20100417";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	if (data->kind_of_flame == FLAME1D_PHYSICS_PREMIXED)	dummy = "PREMIXED";
	if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED)		dummy = "OPPOSED";
	if (data->kind_of_flame == FLAME1D_PHYSICS_TWIN)		dummy = "TWIN";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	// Grid
	dummy = "GRID";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << grid.x;

	// Velocity
	if (data->kind_of_flame == FLAME1D_PHYSICS_OPPOSED || data->kind_of_flame == FLAME1D_PHYSICS_TWIN)
	{
		dummy = "MIXTUREFRACTION";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
		fOutput << Z;
	}
	else
	{
		dummy = "MIXTUREFRACTION";	// TODO
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
		fOutput << grid.x;
	}

	// Velocity
	dummy = "VELOCITY";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << U;

	// Temperature
	dummy = "TEMPERATURE";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << T;

	// Pressure
	BzzVector P_Pa(T.Size()); P_Pa = data->P_Pascal;
	dummy = "PRESSURE";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << P_Pa;

	// Molecular weights
	dummy = "MW";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << PMtot;

	// Molecular weights
	dummy = "G";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << G;

	// Molecular weights
	dummy = "H";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << H;

	// Indices of species for local analysis
	dummy = "INDICES";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << data->iOut;

	// Indices of species for local analysis
	dummy = "MOLEFRACTIONS";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	for(int j=1;j<=NC;j++)
	{
		BzzVector aux = X.GetColumn(j);
		fOutput << aux;
	}
}

void OpenSMOKE_Flame1D::SaveOnBinaryFile(const std::string filename, const bool iSensitivity, const flame1d_model NLS_KIND)
{

	cout << "Save on binary file..." << endl;
	BzzSave fOutput;
	fOutput('*', filename);
	PrintHeaderOnBinaryFile(fOutput);

	cout << "  -- mixture details..." << endl;
	mix->SaveOnBinaryFile(fOutput);

	cout << "  -- flame details..." << endl;
	this->SaveOnBinaryFile(fOutput);
	

	// ROPA
	if (data->iAssignedROPA == true)
	{
		cout << "  -- rate of production analysis..." << endl;
		{
			data->iFocusReactionRates = true;
			BzzVectorInt dummy;
			properties(data->iHOT, -1, dummy, 0, 0);
			data->iFocusReactionRates = false;
		}

		if (data->iVerboseAssignedROPA == true)
		{
			ropa.SetNumberOfPoints(Np);
			ropa.Run(RR, grid.x);
			ropa.PrintRateOfProductionAnalyses(nameFolderAdditionalData + "/ROPA.out", grid.x, T);
			ropa.PrintIntegralRateOfProductionAnalyses(nameFolderAdditionalData + "/ROPA_Integral.out");
		}

		std::string dummy;
		char name[Constants::NAME_SIZE];

		dummy = "ROPA";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));

		dummy = "V20100417";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));

		// Indices of species for rate of production analysis
		BzzVectorInt aux_int(NC);
		for (int j=1;j<=NC;j++)	aux_int[j] = j;
		dummy = "INDICES";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
		fOutput << aux_int;

		// Reaction rates
		dummy = "REACTIONRATES";
		strcpy(name, dummy.c_str());
		fOutput.fileSave.write((char*) name, sizeof(name));
		fOutput << RR;
	}

	// Sensitivity
	if (iSensitivity == true)
	{
		cout << "  -- sensitivity analysis..." << endl;
		solveSensitivity(NLS_KIND, fOutput);
	}

	
	PrintEndOnBinaryFile(fOutput);
	fOutput.End();
			
			
	// Post processor test
	//OpenSMOKE_PostProcessor post_processor;
	//post_processor.ReadFromBinaryFile(filename);

	// Element Flux Analysis
	if (data->iAssignedGlobalElementFluxAnalysis == true)
	{
		cout << "  -- element flux analysis..." << endl;
		{
			ChangeDimensions(grid.Np, mix->NumberOfReactions(), &rForward);
			ChangeDimensions(grid.Np, mix->NumberOfReactions(), &rBackward);
			ChangeDimensions(mix->NumberOfReactions(), &rForwardVector);
			ChangeDimensions(mix->NumberOfReactions(), &rBackwardVector);

			data->iFocusElemetFluxAnalysis = true;
			BzzVectorInt dummy;
			properties(data->iHOT, -1, dummy, 0, 0);
			data->iFocusElemetFluxAnalysis = false;

			for(int j=1;j<=data->element_flux_analysis_manager->n;j++)
				data->element_flux_analysis->RunGlobal(nameFolderAdditionalData + "\\" + data->element_flux_analysis_manager->file_names[j-1], grid.x, data->element_flux_analysis_manager->xA[j-1], data->element_flux_analysis_manager->xB[j-1], rForward, rBackward); 
		}
	}
}


void OpenSMOKE_Flame1D::propertiesSingleReactors(const int i)
{
	molarFractionsAndPMtot();

	// --------------------------------------------------------------------------
	// PROPERTIES FOR DIFFERENT T 
	// --------------------------------------------------------------------------
	{
		// Estrazioni dei vettori delle omega e delle x
		W.GetRow(i,&wVector);
		X.GetRow(i,&xVector);

		// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		double cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
		rho[i] = cTot * PMtot[i];
		urho[i] = 1./rho[i];
		cVector = cTot*xVector;

		// b. Calcolo dei calori specifici [J/kgK]
		//mix->SpeciesCp(T[i]);
		//Cp[i] = mix->MixCp_FromMassFractions(wVector);	
		//Cpk.SetRow(i, mix->Cp);

		// e. Reazioni chimiche
		mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);
		mix->ComputeFromConcentrations( T[i], cVector, cTot, &RVector);					// [kmol/m3/s]
		ElementByElementProduct(RVector, mix->M, &RVector);
		R.SetRow(i, RVector);															// [kg/m3/s]
		//QReaction[i] = - mix->ComputeQReaction(T[i]);									// [J/m3/s]
	}
}

void OpenSMOKE_Flame1D::ODESystem_SingleReactor_Isothermal(BzzVector &x, double t, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recupero variabili fisiche
	// --------------------------------------------------------------------------------------
		recoverPhysicalVariables(SINGLEREACTOR_ISOTHERMAL, x);

	// -------------------------------------------------------------------------------------		
	// 2. Calcolo delle proprieta
	// -------------------------------------------------------------------------------------
		propertiesSingleReactors(indexReactor);
		
	// --------------------------------------------------------------------------------------
	// 3. Equazioni
	// --------------------------------------------------------------------------------------
		give_Opposed_DW(indexReactor);

	// --------------------------------------------------------------------------------------
	// 4. Recupero residui
	// --------------------------------------------------------------------------------------
		recoverResiduals(SINGLEREACTOR_ISOTHERMAL, f);

}

void OpenSMOKE_Flame1D::give_Opposed_DW(const int i)
{
	for(int j=1;j<=NC;j++)
		dW[i][j] =  R[i][j]/rho[i];
}

void OpenSMOKE_Flame1D::SolveODE_SingleReactors(const double tEnd)
{
	ChangeDimensions(NC, &xMin);
	ChangeDimensions(NC, &xMax);
	ChangeDimensions(NC, &xFirstGuess);

	for(int j=1;j<=NC;j++)	xMin[j] = 0.;
	for(int j=1;j<=NC;j++)	xMax[j] = 1.;

	for(indexReactor=1;indexReactor<=Np;indexReactor++)
	{
		cout << "Reactor: " << indexReactor << endl;
		cout << "T: " << T[indexReactor] << " K" << endl;
		cout << "x: " << grid.x[indexReactor]*100. << " cm" << endl;
		
		for(int j=1;j<=NC;j++)
			xFirstGuess[j] = W[indexReactor][j];

		ODE_SingleReactor_Isothermal.assignFlame(this);
		BzzOdeStiffObject o(xFirstGuess, 0., &ODE_SingleReactor_Isothermal);	
		
		//o.StepPrint(DAE_ODE_Print);
		o.SetMinimumConstraints(xMin);
		o.SetMaximumConstraints(xMax);
		//o.SetTolRel(data->rel_dae_Tolerances);
		//o.SetTolAbs(data->abs_dae_Tolerances);
	
		double timeStart = BzzGetCpuTime();
		xFirstGuess = o(tEnd, tEnd);
		cout << "Time DAE solution: " << BzzGetCpuTime() - timeStart << " s" << endl << endl;
	}
}
