/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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
#include "idealreactors/flamelet/OpenSMOKE_Flamelet.h"

OpenSMOKE_Flamelet *ptFlamelet;

const double max_residual = 1.e-10;
const int max_newton_iterations = 30;
const int max_searching_iterations = 100;
const double max_temperature_increment = 10.;
	
 
void ODE_Print(BzzVector &x, double t)
{
	ptFlamelet->MyODE_Print(x, t);
}

void OpenSMOKE_Flamelet::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flamelet"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Flamelet::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Flamelet"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

OpenSMOKE_Flamelet::OpenSMOKE_Flamelet()
{
	name_object = "[Name not assigned]";
}

void OpenSMOKE_Flamelet::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Flamelet::Assign(OpenSMOKE_Flamelet_DataManager *_data)
{
	data = _data;
}

void OpenSMOKE_Flamelet::Assign(OpenSMOKE_Flamelet_ScheduleClass *_operations)
{
	operations = _operations;
}

void OpenSMOKE_Flamelet::Assign(OpenSMOKE_ReactingGas *_mix)
{
	mix = _mix;
}

void OpenSMOKE_Flamelet::Assign(OpenSMOKE_GlobalKinetics *_global)
{
	global = _global;
}

void OpenSMOKE_Flamelet::Allocate_Master_Variables()
{
	cout << Np << endl;
	cout << NC << endl;
	ChangeDimensions(Np,	&T);
	ChangeDimensions(Np,NC, &w);

	if (data->i2E == true)
	{
		ChangeDimensions(Np, &phiN);
		ChangeDimensions(Np, &phiM);
	}
}
void OpenSMOKE_Flamelet::Allocate_Np_Variables()
{
	dimTot	 = dimBlock*Np;
	
	// Time derivatives
	ChangeDimensions(Np,	&dT_over_dt);
	ChangeDimensions(Np,NC, &dw_over_dt);
	
	// Spatial derivatives: First Order
	ChangeDimensions(Np,	&dT_over_dx);
	ChangeDimensions(Np,	&dCp_over_dx);
	ChangeDimensions(Np,NC, &dw_over_dx);

	// Spatial derivatives: Second Order
	ChangeDimensions(Np,	&d2T_over_dx2);
	ChangeDimensions(Np,NC, &d2w_over_dx2);

	// Ode Solver initialization
	ChangeDimensions(dimTot, &initialValues);
	ChangeDimensions(dimTot, &xMin);
	ChangeDimensions(dimTot, &xMax);

	// Properties and kinetic scheme
	ChangeDimensions(Np, NC, &X);
	ChangeDimensions(Np, NC, &R);
	ChangeDimensions(Np, NC, &Cpk);
	ChangeDimensions(Np, &PMtot);
	ChangeDimensions(Np, &uPMtot);
	ChangeDimensions(Np, &rho);
	ChangeDimensions(Np, &Cp);
	ChangeDimensions(Np, &enthalpy);
	ChangeDimensions(Np, &QReaction);
	ChangeDimensions(Np, &chi);

	// LookUp tables
	ChangeDimensions(Np, NC, &CpMap);
	ChangeDimensions(Np, NR, &k1Map);
	ChangeDimensions(Np, NR, &k2Map);
	ChangeDimensions(Np, NR, &uKeqMap);
	ChangeDimensions(Np, NR, &logFcentMap);
	ChangeDimensions(Np, NR, &reactionDSMap);
	ChangeDimensions(Np, NR, &reactionDHMap);

	ChangeDimensions(Np, &hProfile);
	ChangeDimensions(Np, &hDefectProfile);
	ChangeDimensions(Np, mix->NumberOfElements(), &omega_elemental_profile);

	ChangeDimensions(Np, &Qrad);
	ChangeDimensions(Np, &asTot);

	ChangeDimensions(Np,NC, &x_elemental);
	ChangeDimensions(Np,NC, &omega_elemental);
	

	if (data->i2E == true)
	{
		ChangeDimensions(Np, &dphiN_over_dt);
		ChangeDimensions(Np, &dphiM_over_dt);
		ChangeDimensions(Np, &source_phiN);
		ChangeDimensions(Np, &source_phiM);
		ChangeDimensions(Np, &d2phiN_over_dx2);
		ChangeDimensions(Np, &d2phiM_over_dx2);
		ChangeDimensions(Np, NC, &SootGasCorrection);
		ChangeDimensions(Np, &DiffusionSoot);
	}
}

void OpenSMOKE_Flamelet::Allocate()
{
	Allocate_Master_Variables();

	Allocate_Np_Variables();

	ChangeDimensions(NC, &xVector);
	ChangeDimensions(NC, &wVector);
	ChangeDimensions(NC, &RVector);
	ChangeDimensions(NC, &c);
}

void OpenSMOKE_Flamelet::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	int i, j, k;
	double sumCp;

	// -----------------------------------------------------------------
	// Unknowns Recovering
	// -----------------------------------------------------------------
	k=1;
	for(i=1;i<=Np;i++)
	{
		T[i] = x[k++];
		for(j=1;j<=NC;j++)
			w[i][j] = x[k++];
	}

	// -----------------------------------------------------------------
	// Properties evaluation
	// -----------------------------------------------------------------
	Properties(odeFlamelet.jacobianIndex, odeFlamelet.jacobianVariables, dimBlock);
	ChiEvaluation(data->chiSt/data->correctionChi);			


	// -----------------------------------------------------------------
	// Derivatives
	// -----------------------------------------------------------------
	grid.FirstDerivative('C',  T,  dT_over_dx);
	grid.FirstDerivative('C', Cp, dCp_over_dx);
	grid.FirstDerivative('C',  w,  dw_over_dx);
	grid.SecondDerivative(T, d2T_over_dx2);
	grid.SecondDerivative(w, d2w_over_dx2);


	// -----------------------------------------------------------------
	// Energy balance
	// -----------------------------------------------------------------
	dT_over_dt[1] = 0.;
	for(i=2;i<=Ni;i++)
	{	
		dT_over_dt[i] =  0.50*chi[i] * d2T_over_dx2[i] + QReaction[i] / (rho[i]*Cp[i]); 
						
		sumCp = 0.;
		for(j=1;j<=NC;j++)
			sumCp += Cpk[i][j]*dw_over_dx[i][j];

		dT_over_dt[i] += 0.50*chi[i]/Cp[i]*(dCp_over_dx[i] + sumCp)*dT_over_dx[i];
	}
	dT_over_dt[Np] = 0.;


	// -----------------------------------------------------------------
	// Mass balances
	// -----------------------------------------------------------------
	for(j=1;j<=NC;j++)
		dw_over_dt[1][j] = 0.;

	for(i=2;i<=Ni;i++)
		for(j=1;j<=NC;j++)
			dw_over_dt[i][j] = 0.50*chi[i]*d2w_over_dx2[i][j] + R[i][j]/rho[i];

	for(j=1;j<=NC;j++)
		dw_over_dt[Np][j] = 0.;


	// -----------------------------------------------------------------
	// Residuals Recovering
	// -----------------------------------------------------------------
	k=1;
	for(i=1;i<=Np;i++)
	{
		f[k++] = dT_over_dt[i];
		for(j=1;j<=NC;j++)
			f[k++] = dw_over_dt[i][j];
	}
}

void OpenSMOKE_Flamelet::GetSystemFunctions_Enthalpy(BzzVector &x,double t,BzzVector &f)
{
	int i, j, k;
	// -----------------------------------------------------------------
	// Unknowns Recovering
	// -----------------------------------------------------------------
	k=1;
	for(i=1;i<=Np;i++)
		for(j=1;j<=NC;j++)
			w[i][j] = x[k++];

	// -----------------------------------------------------------------
	// Recovering temperatures
	// -----------------------------------------------------------------
	#if LINUX_SO==1
		BzzVector aux(NC);
		for(i=2;i<=Ni;i++)
		{
			aux = w.GetRow(i);
			//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hProfile[i], aux);
			T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hProfile[i], aux, 
				   max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
		}
	#else
		for(i=2;i<=Ni;i++)
			//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hProfile[i], w.GetRow(i));
			T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions_NewVersion(T[i], hProfile[i], w.GetRow(i), 
				   max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
	#endif
	// -----------------------------------------------------------------
	// Properties evaluation
	// -----------------------------------------------------------------
	BzzVectorInt dummy;
	Properties(-1, dummy, dimBlock);
	ChiEvaluation(data->chiSt/data->correctionChi);			

	// -----------------------------------------------------------------
	// Derivatives
	// -----------------------------------------------------------------
	grid.FirstDerivative('C',  w,  dw_over_dx);
	grid.SecondDerivative(w, d2w_over_dx2);

	// -----------------------------------------------------------------
	// Mass balances
	// -----------------------------------------------------------------
	for(j=1;j<=NC;j++)
		dw_over_dt[1][j] = 0.;

	for(i=2;i<=Ni;i++)
		for(j=1;j<=NC;j++)
			dw_over_dt[i][j] = 0.50*chi[i]*d2w_over_dx2[i][j] + R[i][j]/rho[i];

	for(j=1;j<=NC;j++)
		dw_over_dt[Np][j] = 0.;


	// -----------------------------------------------------------------
	// Residuals Recovering
	// -----------------------------------------------------------------
	k=1;
	for(i=1;i<=Np;i++)
		for(j=1;j<=NC;j++)
			f[k++] = dw_over_dt[i][j];
}

void OpenSMOKE_Flamelet::GetSystemFunctions_EnthalpyDefect(BzzVector &x,double t,BzzVector &f)
{
	int i, j, k;
	// -----------------------------------------------------------------
	// Unknowns Recovering
	// -----------------------------------------------------------------
	k=1;
	for(i=1;i<=Np;i++)
		for(j=1;j<=NC;j++)
			w[i][j] = x[k++];

	// -----------------------------------------------------------------
	// Recovering temperatures
	// -----------------------------------------------------------------
	#if LINUX_SO==1


		for(i=Np/2;i<=Np;i++)
		{
			BzzVector aux = w.GetRow(i);
			//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], aux);
			T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], aux, 
				   max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
			if (T[i]<=data->minimumTemperature)
			{
				for(int j=i;j<=Np;j++)
					T[j] = data->minimumTemperature;
				break;
			}
		}
		for(i=Np/2-1;i>=1;i--)
		{
			BzzVector aux = w.GetRow(i);
			//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], aux);
			T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], aux, 
				   max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
			if (T[i]<=data->minimumTemperature)
			{
				for(int j=i;j>=1;j--)
					T[j] = data->minimumTemperature;
				break;
			}
		}
	
	#else
		
		for(i=Np/2;i<=Np;i++)
		{
			//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], w.GetRow(i));
			T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions_NewVersion(T[i], hDefectProfile[i], w.GetRow(i), 
				   max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
			if (T[i]<=data->minimumTemperature)
			{
				for(int j=i;j<=Np;j++)
					T[j] = data->minimumTemperature;
				break;
			}
		}
		for(i=Np/2-1;i>=1;i--)
		{
			//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], w.GetRow(i));
			T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions_NewVersion(T[i], hDefectProfile[i], w.GetRow(i), 
				   max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
			if (T[i]<=data->minimumTemperature)
			{
				for(int j=i;j>=1;j--)
					T[j] = data->minimumTemperature;
				break;
			}
		}

	#endif
	// -----------------------------------------------------------------
	// Properties evaluation
	// -----------------------------------------------------------------
	BzzVectorInt dummy;
	Properties(-1, dummy, dimBlock);
	ChiEvaluation(data->chiSt/data->correctionChi);			

	// -----------------------------------------------------------------
	// Derivatives
	// -----------------------------------------------------------------
	grid.FirstDerivative('C',  w,  dw_over_dx);
	grid.SecondDerivative(w, d2w_over_dx2);

	// -----------------------------------------------------------------
	// Mass balances
	// -----------------------------------------------------------------
	for(j=1;j<=NC;j++)
		dw_over_dt[1][j] = 0.;

	for(i=2;i<=Ni;i++)
		for(j=1;j<=NC;j++)
			dw_over_dt[i][j] = 0.50*chi[i]*d2w_over_dx2[i][j] + R[i][j]/rho[i];

	for(j=1;j<=NC;j++)
		dw_over_dt[Np][j] = 0.;


	// -----------------------------------------------------------------
	// Residuals Recovering
	// -----------------------------------------------------------------
	k=1;
	for(i=1;i<=Np;i++)
		for(j=1;j<=NC;j++)
			f[k++] = dw_over_dt[i][j];
}

void OpenSMOKE_Flamelet::GetSystemFunctions_Soot(BzzVector &x,double t,BzzVector &f)
{
	int i, j, k;
	double sumCp;

	// -----------------------------------------------------------------
	// Unknowns Recovering
	// -----------------------------------------------------------------
	k=1;
	for(i=1;i<=Np;i++)
	{
		T[i]			=	x[k++];
		for(j=1;j<=NC;j++)
			w[i][j]		=	x[k++];
		
		phiN[i]	= x[k++];
		phiM[i]	= x[k++];
	}

	// -----------------------------------------------------------------
	// Properties evaluation
	// -----------------------------------------------------------------
	Properties(odeFlamelet_Soot.jacobianIndex, odeFlamelet_Soot.jacobianVariables, dimBlock);
	ChiEvaluation(data->chiSt/data->correctionChi);

	R += SootGasCorrection;
	
	// -----------------------------------------------------------------
	// Derivatives
	// -----------------------------------------------------------------
	// Energy transport equation
	grid.FirstDerivative('C',  T,  dT_over_dx);
	grid.FirstDerivative('C', Cp, dCp_over_dx);
	grid.SecondDerivative(T, d2T_over_dx2);

	// Species Transport Equations
	grid.FirstDerivative('C',  w,  dw_over_dx);
	grid.SecondDerivative(w, d2w_over_dx2);

	// -----------------------------------------------------------------
	// Soot Properties Transport Equations
	// -----------------------------------------------------------------
	grid.SecondDerivative(phiN, d2phiN_over_dx2);
	grid.SecondDerivative(phiM, d2phiM_over_dx2);

	// Soot Particle Density Equation
	dphiN_over_dt[1]	=	0.;
	for(i=2;i<=Ni;i++)		
		dphiN_over_dt[i] =	0.50*chi[i]*d2phiN_over_dx2[i]	+	// Molecular diffusion 
							source_phiN[i] / rho[i];		; 	// Reaction						
	dphiN_over_dt[Np]	=	0.;

	// Soot Mass Fraction Equation
	dphiM_over_dt[1]	=	0.;
	for(i=2;i<=Ni;i++)		
		dphiM_over_dt[i] =	0.50*chi[i]*d2phiM_over_dx2[i] +		// Molecular diffusion 
							source_phiM[i] / rho[i]; ;				// Reaction						
	dphiM_over_dt[Np]		=	0.;

	// -----------------------------------------------------------------
	// Energy balance
	// -----------------------------------------------------------------

	dT_over_dt[1] = 0.;
	for(i=2;i<=Ni;i++)
	{	
		dT_over_dt[i] =  0.50*chi[i] * d2T_over_dx2[i]  + QReaction[i] / (rho[i]*Cp[i]); 
						
		sumCp = 0.;
		for(j=1;j<=NC;j++)
			sumCp += Cpk[i][j]*dw_over_dx[i][j];

		dT_over_dt[i] += 0.50*chi[i]/Cp[i]*(dCp_over_dx[i] + sumCp)*dT_over_dx[i];
	}
	dT_over_dt[Np] = 0.;

	// -----------------------------------------------------------------
	// Mass balances
	// -----------------------------------------------------------------
	for(j=1;j<=NC;j++)
		dw_over_dt[1][j] = 0.;

	for(i=2;i<=Ni;i++)
		for(j=1;j<=NC;j++)
			dw_over_dt[i][j] = 0.50*chi[i]*d2w_over_dx2[i][j] + R[i][j]/rho[i];

	for(j=1;j<=NC;j++)
		dw_over_dt[Np][j] = 0.;

	// -----------------------------------------------------------------
	// Residuals Recovering
	// -----------------------------------------------------------------
	k=1;
	for(i=1;i<=Np;i++)
	{
		f[k++] = dT_over_dt[i];
		for(j=1;j<=NC;j++)
			f[k++] = dw_over_dt[i][j];
		
		f[k++] = dphiN_over_dt[i];
		f[k++] = dphiM_over_dt[i];
	}
}


void OpenSMOKE_Flamelet::Setup()
{
	NC = mix->NumberOfSpecies();
	NR = mix->NumberOfReactions();

	grid.Construct(data->Np, 1.0, 0.);

	if (data->i2E == false && data->iEnthalpy == false)		dimBlock = 1 + NC;
	if (data->i2E == false && data->iEnthalpy == true)		dimBlock = NC;
	if (data->i2E == true  && data->iEnthalpy == false)		dimBlock = 3 + NC;

	Np = grid.Np;
	Ni = grid.Ni;

	Allocate();

	InitialTemperatureProfile();
	InitialMassFractionProfiles();

	ChiEvaluation(data->chiSt/data->correctionChi);
	prepare_radiation();

	if (data->i2E == true)	// 2EModel
	{
		sootModel.setupFromFile(data->twoEquation_file_name);
		sootModel.assign_mixture(*mix);
		Initial_conditions_soot_module();
	}
}

void OpenSMOKE_Flamelet::Prepare()
{
	int k, i, j;
	double ZERO = 0.0;
	double ONE  = 1.0;

	double TMIN = 250.;
	double TMAX = 6000.;

	if (data->i2E == false && data->iEnthalpy == false)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			initialValues[k++] = T[i];
			for(j=1;j<=NC;j++)
				initialValues[k++] = w[i][j];
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++] = TMIN;
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++] = TMAX;
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;
		}
	}

	else if (data->i2E == false && data->iEnthalpy == true)
	{
		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				initialValues[k++] = w[i][j];

		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xMin[k++] = ZERO;

		k=1;
		for(i=1;i<=Np;i++)
			for(j=1;j<=NC;j++)
				xMax[k++] = ONE;
	}

	else if (data->i2E == true && data->iEnthalpy == false)
	{
		k=1;
		for(i=1;i<=Np;i++)
		{
			initialValues[k++]		= T[i];
			for(j=1;j<=NC;j++)
				initialValues[k++]	= w[i][j];
			initialValues[k++]		= phiN[i];
			initialValues[k++]		= phiM[i];
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMin[k++]		= TMIN;
			for(j=1;j<=NC;j++)
				xMin[k++]	= ZERO;
			xMin[k++]		= ZERO;		// phi_N
			xMin[k++]		= ZERO;		// phi_M
		}

		k=1;
		for(i=1;i<=Np;i++)
		{
			xMax[k++]		= TMAX;
			for(j=1;j<=NC;j++)
				xMax[k++]	= ONE;
			xMax[k++]		= ONE;		// phi_N
			xMax[k++]		= ONE;		// phi_M
		}
	}
}

void OpenSMOKE_Flamelet::Solve(const double tEnd)
{
	double timeStart;

	// -----------------------------------------------------------------
	// File Output Setup
	// -----------------------------------------------------------------
	GnuPlotODE_label(Np);
	GnuPlot_label("Solution_", Np);
	
	// -----------------------------------------------------------------
	// Preparation
	// -----------------------------------------------------------------
	Prepare();

	// -----------------------------------------------------------------
	// Solution 
	// -----------------------------------------------------------------
	bool jEquilibrium = false;
	if (data->chiSt <= 1.e-6) jEquilibrium = true;

	// 1
	if (jEquilibrium == true)
	{
		BzzVector  x_elemental_profile(mix->NumberOfElements());
		BzzVector  x_equilibrium(mix->NumberOfSpecies());
		BzzVector  w_equilibrium(mix->NumberOfSpecies());

		int count_warning = 0;
		double time_start = BzzGetCpuTime();
		
		for(int j=2;j<=Np-1;j++)
		{
			int flag;
			double mw;
			double N_equilibrium=1.;
			double T_equilibrium = T[j-1];
			if (j==2)	x_equilibrium = data->XC;
			else		X.GetRow(j-1, &x_equilibrium);

			

			#if LINUX_SO==1
				BzzVector aux = omega_elemental_profile.GetRow(j);
				mix->GetMWAndElementalMoleFractionsFromElementalMassFractions(mw, x_elemental_profile, aux);
			#else
				mix->GetMWAndElementalMoleFractionsFromElementalMassFractions(mw, x_elemental_profile, omega_elemental_profile.GetRow(j));
			#endif
			mix->Equilibrium_HP(T_equilibrium, x_equilibrium, N_equilibrium, hProfile[j], data->P_Pascal, x_elemental_profile, false, flag);
			
			// Solution
			T[j] = T_equilibrium;
			mix->GetMWAndMassFractionsFromMoleFractions(mw, w_equilibrium, x_equilibrium);
			X.SetRow(j, x_equilibrium);
			w.SetRow(j, w_equilibrium);

			if (flag!=1)
				count_warning++;

			cout << "Warning points: " << count_warning << endl;
		}

		if (data->iEnthalpyDefect == true)
		{
			int jStart = 0;
			for(int j=2;j<=Np;j++)
				if (grid.x[j]>=(1.-data->zStoichiometric))
				{
					jStart = j;
					break;
				}

			cout << "Stochiometric csi: " << 1.-data->zStoichiometric << endl;
			cout << "Stochiometric j:   " << jStart << endl;

			int i;

			for(i=1;i<=Np;i++)
				hDefectProfile[i] = hProfile[i] + data->enthalpyDefect;

			for(i=jStart;i<=Np-1;i++)
			{
					#if LINUX_SO==1
						BzzVector aux = w.GetRow(i);
						//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], aux);
						T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions_NewVersion(T[i], hDefectProfile[i], aux, 
								 max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);

					#else
						//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], w.GetRow(i));
						T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions_NewVersion(T[i], hDefectProfile[i], w.GetRow(i), 
								 max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
					#endif

				
				if (T[i]<=data->minimumTemperature)
				{
					for(int j=i;j<=Np-1;j++)
						T[j] = data->minimumTemperature;
					break;
				}
			}

			for(i=jStart-1;i>=2;i--)
			{
				#if LINUX_SO==1
					BzzVector aux = w.GetRow(i);
					T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], aux);
					T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], aux, 
						 max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
				#else
					//T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions(T[i], hDefectProfile[i], w.GetRow(i));
					T[i] = mix->GetTemperatureFromMassEnthalpyAndMassFractions_NewVersion(T[i], hDefectProfile[i], w.GetRow(i), 
							max_residual, max_newton_iterations, max_searching_iterations, max_temperature_increment);
				#endif
				
				if (T[i]<=data->minimumTemperature)
				{
					for(int j=i;j>=2;j--)
						T[j] = data->minimumTemperature;
					break;
				}
			}

			for(i=jStart;i<=Np-1;i++)
			{
				if (T[i] > data->minimumTemperature)
				{
					int flag;
					double mw;
					double N_equilibrium=1.;
					double T_equilibrium = T[i];
					X.GetRow(i, &x_equilibrium);

					#if LINUX_SO==1
						BzzVector aux = omega_elemental_profile.GetRow(i);
						mix->GetMWAndElementalMoleFractionsFromElementalMassFractions(mw, x_elemental_profile, aux);
					#else
						mix->GetMWAndElementalMoleFractionsFromElementalMassFractions(mw, x_elemental_profile, omega_elemental_profile.GetRow(i));
					#endif
					mix->Equilibrium_HP(T_equilibrium, x_equilibrium, N_equilibrium, hDefectProfile[i], data->P_Pascal, x_elemental_profile, false, flag);
					
					// Solution
					if (T_equilibrium>data->minimumTemperature)
					{
						T[i] = T_equilibrium;
						mix->GetMWAndMassFractionsFromMoleFractions(mw, w_equilibrium, x_equilibrium);
						X.SetRow(i, x_equilibrium);
						w.SetRow(i, w_equilibrium);
					}
					else
					{
						for(int j=i;j<=Np-1;j++)
							T[j] = data->minimumTemperature;
						break;
					}
				}
			}

			for(i=jStart-1;i>=2;i--)
			{
				if (T[i] > data->minimumTemperature)
				{
					int flag;
					double mw;
					double N_equilibrium=1.;
					double T_equilibrium = T[i];
					X.GetRow(i, &x_equilibrium);

					#if LINUX_SO==1
						BzzVector aux = omega_elemental_profile.GetRow(i);
						mix->GetMWAndElementalMoleFractionsFromElementalMassFractions(mw, x_elemental_profile, aux);
					#else
						mix->GetMWAndElementalMoleFractionsFromElementalMassFractions(mw, x_elemental_profile, omega_elemental_profile.GetRow(i));
					#endif
					mix->Equilibrium_HP(T_equilibrium, x_equilibrium, N_equilibrium, hDefectProfile[i], data->P_Pascal, x_elemental_profile, false, flag);
					
					// Solution
					if (T_equilibrium>data->minimumTemperature)
					{
						T[i] = T_equilibrium;
						mix->GetMWAndMassFractionsFromMoleFractions(mw, w_equilibrium, x_equilibrium);
						X.SetRow(i, x_equilibrium);
						w.SetRow(i, w_equilibrium);
					}
					else
					{
						for(int j=i;j>=2;j--)
							T[j] = data->minimumTemperature;
						break;
					}
				}
			}
		}
		
		double time_end = BzzGetCpuTime();
		cout << "Equilibrium time:   " << -time_start+time_end << endl;
		cout << "Number of warnings: " << count_warning << endl;

		PrintGnuPlot();		
	}

	// 2
	if (data->i2E == false && jEquilibrium == false && data->iEnthalpy == false)
	{
		ptFlamelet = this;
		odeFlamelet.assignFlamelet(this);
		BzzOdeSparseStiffObject o(initialValues, 0., &odeFlamelet, dimBlock);
		o.SetMinimumConstraints(xMin);
		o.SetMaximumConstraints(xMax);
		
		o.SetTolRel(data->relTolerances*MachEps());
		o.SetTolAbs(data->absTolerances);

		o.StepPrint(ODE_Print);
		Video_label();

		timeStart = BzzGetCpuTime();
		o(tEnd, tEnd);
		fGnuPlotODE.close();
		PrintGnuPlot();

		cout << endl;
		cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian() << endl;
		cout << "Numerical Jacobians:             "	<< o.GetNumNumericalJacobian() << endl;
		cout << "Time ODE solution:               "	<< BzzGetCpuTime() - timeStart << " s" << endl << endl;
	}

	// 3
	if (data->i2E == false && jEquilibrium == false && data->iEnthalpy == true && data->iEnthalpyDefect == false)
	{
		ptFlamelet = this;
		odeFlamelet_Enthalpy.assignFlamelet(this);
		
		BzzOdeSparseStiffObject o(initialValues, 0., &odeFlamelet_Enthalpy, dimBlock);

		o.SetMinimumConstraints(xMin);
		o.SetMaximumConstraints(xMax);

		o.SetTolRel(data->relTolerances*MachEps());
		o.SetTolAbs(data->absTolerances);

		o.StepPrint(ODE_Print);
		Video_label();

		timeStart = BzzGetCpuTime();
		o(tEnd, tEnd);
		fGnuPlotODE.close();
		PrintGnuPlot();

		cout << endl;
		cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian() << endl;
		cout << "Numerical Jacobians:             "	<< o.GetNumNumericalJacobian() << endl;
		cout << "Time ODE solution:               "	<< BzzGetCpuTime() - timeStart << " s" << endl << endl;
	}

	// 3 - Enthalpy Defect
	if (data->i2E == false && jEquilibrium == false && data->iEnthalpy == true && data->iEnthalpyDefect == true)
	{
		for(int j=1;j<=Np;j++)
			hDefectProfile[j] = hProfile[j] + data->enthalpyDefect;

		ptFlamelet = this;
		odeFlamelet_EnthalpyDefect.assignFlamelet(this);
		BzzOdeSparseStiffObject o(initialValues, 0., &odeFlamelet_EnthalpyDefect, dimBlock);
		o.SetMinimumConstraints(xMin);
		o.SetMaximumConstraints(xMax);
		
		o.SetTolRel(data->relTolerances*MachEps());
		o.SetTolAbs(data->absTolerances);


		o.StepPrint(ODE_Print);
		Video_label();
		timeStart = BzzGetCpuTime();
		o(tEnd, tEnd);
		fGnuPlotODE.close();
		PrintGnuPlot();

		cout << endl;
		cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian() << endl;
		cout << "Numerical Jacobians:             "	<< o.GetNumNumericalJacobian() << endl;
		cout << "Time ODE solution:               "	<< BzzGetCpuTime() - timeStart << " s" << endl << endl;
	}

	// 4
	if (data->i2E == true && jEquilibrium == false)
	{
		ptFlamelet = this;
		odeFlamelet_Soot.assignFlamelet(this);
		BzzOdeSparseStiffObject o(initialValues, 0., &odeFlamelet_Soot, dimBlock);
		o.SetMinimumConstraints(xMin);
		o.SetMaximumConstraints(xMax);

		o.SetTolRel(data->relTolerances*MachEps());
		o.SetTolAbs(data->absTolerances);

		o.StepPrint(ODE_Print);
	
		Video_label();

		timeStart = BzzGetCpuTime();
		o(tEnd, tEnd);
	
		fGnuPlotODE.close();
		PrintGnuPlot();

		cout << endl;
		cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian() << endl;
		cout << "Numerical Jacobians:             "	<< o.GetNumNumericalJacobian() << endl;
		cout << "Time ODE solution:               "	<< BzzGetCpuTime() - timeStart << " s" << endl << endl;
	}
}

 
void OpenSMOKE_Flamelet::MyODE_Print(BzzVector &x, double t)
{
	if(data->iterationVideoCounter == data->nStepsVideo)
	{
		cout << data->iteration	<< "\t";
		cout << t				<< "\t";
		if (data->i2E == false)
		{
			cout << T.Max()		<< "\t";
			cout << T.Min()		<< "\t";
		}

		if (data->i2E == true)
		{
			cout << T.Max()		<< "\t";
			cout << phiN.Max()	<< "\t";
			cout << phiM.Max()	<< "\t";
		}
		cout << endl;

		data->iterationVideoCounter = 0;
	}

	if(data->iterationFileCounter == data->nStepsFile)
	{
		PrintGnuPlotODE(t);
		data->iterationFileCounter = 0;
	}

	data->iteration++;
	data->iterationVideoCounter++;
	data->iterationFileCounter++;
}

void OpenSMOKE_Flamelet::Video_label()
{
	cout.setf(ios::scientific);

	cout << endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "ODE System Solution: NC = " << NC << " NP = " << Np << endl;
	cout << "-------------------------------------------------------------" << endl;

	cout << "-------------------------------------------------------------" << endl;
	cout << "#Iter."		<< "\t"
		 << "Time [s]"		<< "\t"
		 << "Tmax [K]"		<< "\t"
		 << endl;
	cout << "-------------------------------------------------------------" << endl;
}

void OpenSMOKE_Flamelet::Video_label_soot()
{
	cout.setf(ios::scientific);

	cout << endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "ODE System Solution - Soot: NC = " << NC << " NP = " << Np << endl;
	cout << "-------------------------------------------------------------" << endl;

	cout << "-------------------------------------------------------------" << endl;
	cout << "#Iter."			<< "\t"
		 << "Time [s]"			<< "\t"
		 << "T_Max [K]"			<< "\t"
		 << "N_Max [kmol/kg]"	<< "\t"
		 << "wSoot_Max [-]"		<< "\t"
		 << endl;
	cout << "-------------------------------------------------------------" << endl;
}

void OpenSMOKE_Flamelet::GnuPlotODE_label(const int N)
{
	int j;

	std::string fileName;
	char numberOfPoints[4];

	my_itoa(N, numberOfPoints, 10);

	fileName  = nameFolderUnsteadyData + "/ODESolution_";
	fileName += numberOfPoints;
	fileName += ".out";

	openOutputFileAndControl(fGnuPlotODE, fileName);
	fGnuPlotODE.setf(ios::scientific);

	fGnuPlotODE << "Time[s]"	<< "\t";
	fGnuPlotODE << "f[-]"		<< "\t\t";
	fGnuPlotODE << "T[K]"		<< "\t\t";

	int count = 4;
	for(j=1;j<=data->iOut.Size();j++)
		fGnuPlotODE << data->nameOutput[j] << " [X" << count++ << "]  \t";
	fGnuPlotODE << endl << endl;
}

void OpenSMOKE_Flamelet::GnuPlot_label(const std::string name, const int N)
{
	int j;

	std::string fileName;

	if (N>0)
	{
		char numberOfPoints[4];
		my_itoa(N, numberOfPoints, 10);
		fileName  = nameFolderSteadyData + "/" + name;
		fileName += numberOfPoints;
		fileName += ".out";
	}
	else if (N<0)
		fileName = nameFolderSteadyData + "/" + name;
	else if (N==0)
		fileName = nameFolderProfilesData + "/" + name;

	openOutputFileAndControl(fGnuPlot, fileName);
	fGnuPlot.setf(ios::scientific);

	fGnuPlot << "#Point"	<< "\t\t";
	fGnuPlot << "f[-]"		<< "\t\t";
	fGnuPlot << "T[K]"		<< "\t\t";
	fGnuPlot << "H[J/kmol]"	<< "\t";
	fGnuPlot << "H[J/kg]  "	<< "\t";
	fGnuPlot << "P[Pa]	  "	<< "\t";
	fGnuPlot << "chi[Hz]  "	<< "\t";

	int count = 8;
	for(j=1;j<=NC;j++)
		fGnuPlot << mix->names[j] << " [X" << count++ << "]  \t";

	for(j=1;j<=NC;j++)
		fGnuPlot << mix->names[j] << " [W" << count++ << "]  \t";

	fGnuPlot << "chi[1/s]("		<< count++ << ")"	<< "\t";
	fGnuPlot << "Cp[J/kgK]("	<< count++ << ")"	<< "\t";
	fGnuPlot << "MW("			<< count++ << ")"	<< "\t";
	fGnuPlot << "rho[kg/m3]("	<< count++ << ")"	<< "\t";
	fGnuPlot << "as[1/m]("		<< count++ << ")"	<< "\t";
	fGnuPlot << "Qrad[W/m3]("	<< count++ << ")"	<< "\t";

	// Elemental
	for(j=1;j<=mix->NumberOfElements();j++)
		fGnuPlot << mix->list_of_elements[j-1] << " [X" << count++ << "]  \t";
	for(j=1;j<=mix->NumberOfElements();j++)
		fGnuPlot << mix->list_of_elements[j-1] << " [W" << count++ << "]  \t";

	fGnuPlot << endl << endl;

	if (data->i2E == true)
	{
		openOutputFileAndControl(fGnuPlotSoot, nameFolderAdditionalData + "/Soot_2E.out");
		fGnuPlotSoot.setf(ios::scientific);

        fGnuPlotSoot << setw(20) << left << "f[-]";
        sootModel.GnuPlotInterface(fGnuPlotSoot, 2);
		fGnuPlotSoot << endl << endl;
	}
}

void OpenSMOKE_Flamelet::PrintGnuPlotODE(const double t)
{
	int i, j;
	for(i=1;i<=Np;i++)
	{
		fGnuPlotODE	<< t <<			"\t"
					<< grid.x[i] << "\t"
					<< T[i]   <<	"\t";
		
		for(j=1;j<=data->iOut.Size();j++)
			fGnuPlotODE << X[i][data->iOut[j]] << "\t";
		
		fGnuPlotODE << endl;
	}

	fGnuPlotODE << endl;		
}

void OpenSMOKE_Flamelet::PrintGnuPlot()
{
	int i, j;
	
	calculate_radiation();

	for(i=1;i<=Np;i++)
	{
		fGnuPlot	<< i							<< "\t\t"
					<< grid.x[i]					<< "\t"
					<< T[i]							<< "\t"
					<< enthalpy[i]*PMtot[i]			<< "\t"		// J/kmol
					<< enthalpy[i]					<< "\t"		// J/kg
					<< data->P_Pascal				<< "\t"		// Pa
					<< data->chiSt					<< "\t";	// Hz
		
		for(j=1;j<=NC;j++)
			fGnuPlot << X[i][j] << "\t";
		for(j=1;j<=NC;j++)
			fGnuPlot << w[i][j] << "\t";

		fGnuPlot	<< chi[i]		<<	"\t";
		fGnuPlot	<< Cp[i]		<<	"\t";
		fGnuPlot	<< PMtot[i]		<<	"\t";
		fGnuPlot	<< rho[i]		<<	"\t";
		fGnuPlot	<< asTot[i]		<<	"\t";
		fGnuPlot	<< -Qrad[i]		<<	"\t";

		// Elemental
		ElementalAnalysis();
		for(j=1;j<=mix->NumberOfElements();j++)
			fGnuPlot << x_elemental[i][j]		<< "\t";
		for(j=1;j<=mix->NumberOfElements();j++)
			fGnuPlot << omega_elemental[i][j]	<< "\t";

		fGnuPlot	<< endl;

		if (data->i2E == true)
		{
			X.GetRow(i, &xVector);
			sootModel.update(T[i], data->P_atm, rho[i], 0., xVector, phiN[i], phiM[i]);
			sootModel.formation_rates();

			fGnuPlotSoot << setw(20) << left << grid.x[i];
			sootModel.write_on_file(fGnuPlotSoot, phiN[i], phiM[i]);
		}	
	}

	fGnuPlot.close();	
	fGnuPlotSoot.close();	
}

void OpenSMOKE_Flamelet::InitialTemperatureProfile()
{
	int i;

	for(i=1;i<=Np;i++)
	{
		if(grid.x[i]<=data->xcen) T[i] = data->TC+(data->Tpeak-data->TC)/data->xcen*grid.x[i];
		if(grid.x[i] >data->xcen) T[i] = data->Tpeak-(data->Tpeak-data->TO)/(1.-data->xcen)*(grid.x[i]-data->xcen);
	}
}

void OpenSMOKE_Flamelet::InitialMassFractionProfiles()
{
	int i, j;
	double sum;

	// b. Reactants and main products
	for(i=1;i<=Np;i++)
	{
//		for(j=1;j<=data->iXC.Size();j++)
//			X[i][data->iXC[j]] = data->XC[data->iXC[j]]*(1.-grid.x[i]);
		
//		for(j=1;j<=data->iXO.Size();j++)
//			X[i][data->iXO[j]] = data->XO[data->iXO[j]]*grid.x[i];

		for(j=1;j<=NC;j++)
			X[i][j] = data->XC[j] + (data->XO[j]-data->XC[j])*grid.x[i];
		
//		for(j=1;j<=data->iXO.Size();j++)
//			X[i][data->iXO[j]] = data->XO[data->iXO[j]]*grid.x[i];

	}		

	// c. Inert
	for(i=1;i<=Np;i++)
	{
		sum = 0.;
		for(j=1;j<=NC;j++)
			if(j!=data->jINERT) sum+=X[i][j];
		X[i][data->jINERT] = 1.-sum;

		if (sum > 1.0) ErrorMessage("WARNING: wrong initial composition!! ");
	}

	// 3. Frazioni Massive
	// --------------------------------------------------------------------------
	MassFractionsAndPMtot();

	// 4. Enthalpy profile
	for(i=1;i<=Np;i++)
		hProfile[i] = data->enthalpyC + (data->enthalpyO-data->enthalpyC)*grid.x[i];	// Linear profile

	for(int j=1;j<=mix->NumberOfElements();j++)
		for(int i=1;i<=Np;i++)
			omega_elemental_profile[i][j] = data->omega_elemental_C[j] + (data->omega_elemental_O[j]-data->omega_elemental_C[j])*grid.x[i];
}

void OpenSMOKE_Flamelet::MoleFractionsAndPMtot()
{
	int i, j;

	for(i=1;i<=Np;i++)
	{
		uPMtot[i] = 0.;
		for(j=1;j<=NC;j++)
			uPMtot[i] += w[i][j]*mix->uM[j];
		PMtot[i] = 1./uPMtot[i];
		for(j=1;j<=NC;j++)
			X[i][j] = w[i][j] * PMtot[i] * mix->uM[j];
	}
}

void OpenSMOKE_Flamelet::MassFractionsAndPMtot()
{
	int i, j;
	for(i=1;i<=Np;i++)
	{
		PMtot[i] = 0.;
		for(j=1;j<=NC;j++)
			PMtot[i] += X[i][j]*mix->M[j];
		uPMtot[i] = 1./PMtot[i];
		for(j=1;j<=NC;j++)
			w[i][j] = X[i][j] * uPMtot[i] * mix->M[j];
	}
}

void OpenSMOKE_Flamelet::ChiEvaluation(const double chi0)
{
	int i;

	for(i=1;i<=Np;i++)
		chi[i] = chi0*exp(-2.*BzzPow2(ErfInv.at(2.*grid.x[i]-1.)));

	if (data->iDensityCorrection == true)
		for(i=1;i<=Np;i++)
			chi[i] *= 3./4. * BzzPow2(sqrt(rho[Np]/rho[i])+1.) / (2.*sqrt(rho[Np]/rho[i]) + 1.);
}

double OpenSMOKE_Flamelet::StoichiometricChiEvaluation(const double chi0)
{
	double as;
	double rhoSt;
	double chiSt;
	double zStar = 1.-data->zStoichiometric;
	double Cc = 1.0;

	if (data->iDensityCorrection == true)
	{
		int j=0;
		for(int i=1;i<=Np;i++)
			if (grid.x[i]>=zStar)	
			{
				j = i-1;
				break;
			}

			rhoSt = rho[j] + (rho[j+1]-rho[j])/(grid.x[j+1]-grid.x[j])*(zStar-grid.x[j]);
			Cc	  = 3./4. * BzzPow2(sqrt(rho[Np]/rhoSt)+1.) / (2.*sqrt(rho[Np]/rhoSt) + 1.);
	}

	chiSt	= chi0*exp(-2.*BzzPow2(ErfInv.at(2.*zStar-1.))) * Cc;
	as		= Constants::pi*chiSt;									// [????]

	cout << " Correction coefficient:    " << Cc		<< endl;
	cout << " Maximum Chi:               " << chi0	<< " 1/s" << endl;
	cout << " Stoichiometric Chi:        " << chiSt	<< " 1/s" << endl;
	cout << " Caracteristic Strain rate: " << as		<< " 1/s" << endl;
	
	return chiSt;
}


void OpenSMOKE_Flamelet::Properties()
{
	BzzVectorInt dummy;
	Properties(0, dummy, 0);
}

void OpenSMOKE_Flamelet::Properties(const int jjacobianIndex, BzzVectorInt &jjacobianVariables, const int nBlock)
{
	// jacobianIndex = -2 || -1 : valutazione normale delle proprieta
	// jacobianIndex =  0 : memorizzazione delle proprieta
	
	int i, j;
	int indexT = 1;
	int choice;
	BzzVectorInt pointToUpdate;

	if (jjacobianIndex < 0) choice = 1;			//	Calcolo di tutte le proprieta senza memorizzazione
	else if (jjacobianIndex == 0) choice = 2;	//	Calcolo di tutte le proprieta con memorizzazione
	else //if (jjacobianIndex > 0)
	{
		choice = 3;
		for(j=1;j<=jjacobianVariables.Size();j++)
			if ( (jjacobianVariables[j]%(nBlock)) == indexT)
			{
				pointToUpdate.Append(int(jjacobianVariables[j]/(nBlock))+1);
			}
	}
	

	double cTot;

	MoleFractionsAndPMtot();

	// --------------------------------------------------------------------------
	// PROPERTIES FOR DIFFERENT T 
	// --------------------------------------------------------------------------
	if (choice == 1)
	{
		for(i=1;i<=Np;i++)
		{
			// Estrazioni dei vettori delle omega e delle x
			w.GetRow(i,&wVector);
			X.GetRow(i,&xVector);

			// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
			cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			rho[i] = cTot * PMtot[i];
			c = cTot*xVector;

			// b. Calcolo dei calori specifici [J/kgK]
			mix->SpeciesCp(T[i]);
			Cp[i] = mix->MixCp_FromMassFractions(wVector);	
			Cpk.SetRow(i, mix->Cp);

			// e. Reazioni chimiche
			{
				mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);	
				mix->ComputeFromConcentrations( T[i], c, cTot, &RVector);	// [kmol/m3/s]
				ElementByElementProduct(RVector, mix->M, &RVector);
				R.SetRow(i, RVector);										// [kg/m3/s]
				QReaction[i] = - mix->ComputeQReaction(T[i]);				// [J/m3.s]
			}
			
			// This enthalpy is in [J/kg] 
			enthalpy[i] = mix->GetMixEnthalpy_Mass(T[i], wVector);		// [J/kg]
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
			
			// d. Calcolo dei coefficienti di diffusione
			mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);
			k1Map.SetRow(i,mix->k1);
			k2Map.SetRow(i,mix->k2);
			uKeqMap.SetRow(i,mix->uKeq);
			logFcentMap.SetRow(i,mix->logFcent);
			reactionDHMap.SetRow(i,mix->reactionDH);
			reactionDSMap.SetRow(i,mix->reactionDS);
		}

		choice = 3;
	}


	if(choice == 3) // Da chiamare solo se indexJacobian>0 e la Tnon modificata
	{
		for(i=1;i<=Np;i++)
		{
			// Estrazioni dei vettori delle omega e delle x
			w.GetRow(i,&wVector);
			X.GetRow(i,&xVector);

			// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
			cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			rho[i] = cTot * PMtot[i];
			c = cTot*xVector;

			// b. Calcolo dei calori specifici [J/kgK]
			CpMap.GetRow(i, &mix->Cp);
			Cpk.SetRow(i, CpMap.GetRow(i));
			Cp[i] = mix->MixCp_FromMassFractions(wVector);	

			// e. Reazioni chimiche
			{
				k1Map.GetRow(i,&mix->k1);
				k2Map.GetRow(i,&mix->k2);
				uKeqMap.GetRow(i,&mix->uKeq);
				logFcentMap.GetRow(i,&mix->logFcent);
				reactionDHMap.GetRow(i,&mix->reactionDH);
				reactionDSMap.GetRow(i,&mix->reactionDS);
		
				mix->ComputeFromConcentrations( T[i], c, cTot, &RVector);// [kmol/m3/s]
				ElementByElementProduct(RVector, mix->M, &RVector);
				R.SetRow(i, RVector);			// [kg/m3/s]
				QReaction[i] = - mix->ComputeQReaction(T[i]);	// [J/m3.s]
			}

			// This enthalpy is in [J/kg] 
			enthalpy[i] = mix->GetMixEnthalpy_Mass(T[i], wVector);		// [J/kg]
		}
	}

	for(j=1;j<=pointToUpdate.Size();j++)
	{
		i = pointToUpdate[j];

		// Estrazioni dei vettori delle omega e delle x
		w.GetRow(i,&wVector);
		X.GetRow(i,&xVector);

		// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
		rho[i] = cTot * PMtot[i];
		c = cTot*xVector;

		// b. Calcolo dei calori specifici [J/kgK]
		mix->SpeciesCp(T[i]);
		Cp[i] = mix->MixCp_FromMassFractions(wVector);	
		Cpk.SetRow(i, mix->Cp);

		// e. Reazioni chimiche
		{
			mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);
			mix->ComputeFromConcentrations( T[i], c, cTot, &RVector);	// [kmol/m3/s]
			ElementByElementProduct(RVector, mix->M, &RVector);
			R.SetRow(i, RVector);										// [kg/m3/s]
			QReaction[i] = - mix->ComputeQReaction(T[i]);				// [J/m3.s]
		}
			
		// This enthalpy is in [J/kg] 
		enthalpy[i] = mix->GetMixEnthalpy_Mass(T[i], wVector);		// [J/kg]
	}

	// Global Kinetics
	if (data->iGlobalKinetics == 1)
	{
		for(i=1;i<=Np;i++)
		{
			// Estrazioni dei vettori delle omega e delle x
			w.GetRow(i,&wVector);
			X.GetRow(i,&xVector);

			// a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
			cTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);
			rho[i] = cTot * PMtot[i];
			c = cTot*xVector;

			global->GiveMeFormationRates(T[i],c,RVector);
			global->GiveMeReactionHeat(T[i], RVector, QReaction[i]);
			R.SetRow(i, RVector);			// [kg/m3/s]
		}
	}

	if (data->i2E == true)
	{
		for(i=2;i<=Np;i++)
		{
			double DiffC = 0.;

			X.GetRow(i, &xVector);
			sootModel.update(T[i], data->P_atm, rho[i], DiffC, xVector, phiN[i], phiM[i]);
			sootModel.formation_rates();
			source_phiN[i] = sootModel.s;
			source_phiM[i] = sootModel.S;
			for (int j=1;j<=NC-1;j++)
				SootGasCorrection[i][j] = sootModel.SGas[j];	// Gas species
			SootGasCorrection[i][NC]	= sootModel.S;			// Soot
	
			DiffusionSoot[i] = sootModel.Diff;
		}
	}
}

void OpenSMOKE_Flamelet::RefineGridPeak(const double fraction)
{
	BzzVectorInt listAddPoints;

	double Tmax = T.Max();
	for(int i=1;i<=Np;i++)
		if (T[i]>=(fraction*Tmax)) listAddPoints.Append(i);

	AddPoints(listAddPoints);
}

void OpenSMOKE_Flamelet::RefineLeanSide(const double fraction)
{
	bool iMax = false;
	double Tmax = T.Max();
	BzzVectorInt listAddPoints;

	for(int i=Ni;i>=1;i--)
	{
		if (iMax == false)
		{
			if (T[i]<=T[i+1])	iMax = true;
			listAddPoints.Append(i);
		}
		else
		{
			if (T[i]>=(fraction*Tmax)) listAddPoints.Append(i);
		}
	}

	AddPoints(listAddPoints);
}

void OpenSMOKE_Flamelet::DoubleTheGrid()
{
	// Refine the grid
	grid.RefineDouble();
	Np = grid.Np;
	Ni = grid.Ni;

	// Linear Interpolation
	grid.DoubleField(T);
	grid.DoubleField(w);
	
	if (data->i2E == true)
	{
		grid.DoubleField(phiN);
		grid.DoubleField(phiM);
	}

	// Properties updating
	Allocate_Np_Variables();
	Properties();
	ChiEvaluation(data->chiSt/data->correctionChi);

	// Update enthalpy profile
	for(int i=1;i<=Np;i++)
		hProfile[i] = data->enthalpyC + (data->enthalpyO-data->enthalpyC)*grid.x[i];	// Linear profile

	for(int j=1;j<=mix->NumberOfElements();j++)
		for(int i=1;i<=Np;i++)
			omega_elemental_profile[i][j] = data->omega_elemental_C[j] + (data->omega_elemental_O[j]-data->omega_elemental_C[j])*grid.x[i];
}

void OpenSMOKE_Flamelet::AddPoints(BzzVectorInt &listPoints)
{
	int i;

	// Refine the grid
	Sort(&listPoints);
	for (i=1;i<=listPoints.Size();i++)	// This is correct only if listPoints is sorted
		grid.Refine(listPoints[i]+(i-1));	// Min --> Max

	Np = grid.Np;
	Ni = grid.Ni;

	grid.AddPointsField(T, listPoints);
	grid.AddPointsField(w, listPoints);

	if (data->i2E == true)
	{
		grid.AddPointsField(phiN, listPoints);
		grid.AddPointsField(phiM, listPoints);
	}

	// Properties updating
	Allocate_Np_Variables();
	Properties();
	ChiEvaluation(data->chiSt/data->correctionChi);

	// Update enthalpy profile
	for(i=1;i<=Np;i++)
		hProfile[i] = data->enthalpyC + (data->enthalpyO-data->enthalpyC)*grid.x[i];	// Linear profile

	for(int j=1;j<=mix->NumberOfElements();j++)
		for(int i=1;i<=Np;i++)
			omega_elemental_profile[i][j] = data->omega_elemental_C[j] + (data->omega_elemental_O[j]-data->omega_elemental_C[j])*grid.x[i];
}

void OpenSMOKE_Flamelet::NewPoints(const std::string KIND, const char index)
{
	BzzVector phi;
	BzzVectorInt pointList;

	int count = 0;

	if		("TEMPERATURE" == KIND)	phi = T;
	else if ("QREACTION" == KIND)	phi = QReaction;
	else	ErrorMessage("You can adapt the grid only on the temperature or Qreaction profile!!");
	
	cout << "Adding Points ("<< index << ") using " << KIND << " profile!"<< endl;

	// Difference
	if (index == 'D')
		pointList = grid.QueryNewPointsDifference(data->nDiff, data->deltaDiff, phi);

	// Gradient
	if (index == 'G')
		pointList = grid.QueryNewPointsGradient(data->nGrad, data->deltaGrad, phi);

	AddPoints(pointList);
	
	cout << " New points added: "		<< pointList.Size() << endl;
	cout << " Total number of points: " << Np				<< endl;
}

void OpenSMOKE_Flamelet::RecoverFromBackUp(const std::string fileName)
{
	int i, j;

	std::string stringa;
	std::string firstSpecies;
	std::string *nameSpecies;
	int numberOfSpecies;
	int MaxNumberOfSpecies = 500;
	int nFound;
	BzzVector coordinate;
	BzzVector temperature;
	BzzMatrix molarfractions;
	BzzMatrix massfractions;
	BzzVector vector;
	double dummy;

	LinearInterpolation dummylinear;


	nameSpecies = new std::string[MaxNumberOfSpecies+1]; 
	
	ifstream fInput;
	openInputFileAndControl(fInput, fileName);
	fInput.setf(ios::scientific);

	fInput >> stringa;	// Point
	fInput >> stringa;	// Radial Coordinate
	fInput >> stringa;  // Temperature
	fInput >> stringa;  // Enthalpy (molar)
	fInput >> stringa;	// Enthalpy (mass)
	fInput >> stringa;  // pressure
	fInput >> stringa;	// chi reference

	fInput >> firstSpecies;					// First species name
	fInput >> stringa;
	nameSpecies[1] = firstSpecies;
	
	numberOfSpecies = 1;
	for(;;)
	{
		fInput >> stringa;
		if (stringa == firstSpecies)
			break;
		else
		{
			numberOfSpecies++;
			nameSpecies[numberOfSpecies] = stringa;
		}
		fInput >> stringa;
	}

	for(j=1;j<=2*numberOfSpecies;j++)
		fInput >> stringa;

	char commento[Constants::COMMENT_SIZE];
	fInput.getline(commento, Constants::COMMENT_SIZE);


	ChangeDimensions(numberOfSpecies, &vector);

	i=0;
	do 
	{
		i++;
		fInput >> dummy;
		fInput >> dummy; coordinate.Append(dummy);
		fInput >> dummy; temperature.Append(dummy);
		fInput >> dummy;
		fInput >> dummy;
		fInput >> dummy;
		fInput >> dummy;
		
		for(j=1;j<=numberOfSpecies;j++)
			fInput >> vector[j];
		molarfractions.AppendRow(vector);
		
		for(j=1;j<=numberOfSpecies;j++)
			fInput >> vector[j];
		massfractions.AppendRow(vector);

		fInput.getline(commento, Constants::COMMENT_SIZE);

	} while (coordinate[i] < 1.00);

	fInput.close();

	cout << "Initialize profile from file " << fileName << "!" << endl;
	cout << "Number of species: "			<< numberOfSpecies << endl;

	// Costruzione della griglia non equispaziata
	//if (coordinate.Size()!=Np)
	//	ErrorMessage("The number of points in the backup file doesn't match with the number in Data.dat\n");
	Np = coordinate.Size();
	Ni = Np-1;
	NC = mix->NumberOfSpecies();
	NR = mix->NumberOfReactions();
	
	if (data->i2E == false && data->iEnthalpy == false)		dimBlock = 1 + NC;
	if (data->i2E == false && data->iEnthalpy == true)		dimBlock =     NC;
	if (data->i2E == true  && data->iEnthalpy == false)		dimBlock = 3 + NC;


	Allocate();

	BzzVectorInt iFound(NC);
	nFound=0;
	for(i=1;i<=NC;i++)
	{
		for(j=1;j<=numberOfSpecies;j++)
			if (mix->names[i] == nameSpecies[j])
			{
				iFound[i]=j;
				nFound++;
				break;
			}
	}
	cout << "Number of species found: "	<< nFound << endl;

	// -----------------------------------------------------------------
	// File Output Setup
	// -----------------------------------------------------------------
	GnuPlotODE_label(Np);
	GnuPlot_label("Solution_", Np);

	grid.Construct(coordinate);
	InitialTemperatureProfile();
	InitialMassFractionProfiles();
	ChiEvaluation(data->chiSt/data->correctionChi);

	if (data->i2E == true)	// 2EModel
	{
		sootModel.setupFromFile(data->twoEquation_file_name);
		sootModel.assign_mixture(*mix);
		Initial_conditions_soot_module();
	}


	// 0. Temperature
	// --------------------------------------------------------------------------
	dummylinear(coordinate, temperature);
	for(i=2;i<=Ni;i++)
		T[i] = dummylinear(grid.x[i]);

	// 1. MolarFractions
	// --------------------------------------------------------------------------
	for(j=1;j<=NC;j++)
	{
		if (iFound[j]!=0)
		{
			BzzVector aux = molarfractions.GetColumn(iFound[j]);
			dummylinear(coordinate, aux);
			for(i=2;i<=Ni;i++)
				X[i][j] = dummylinear(grid.x[i]);
		}
		else if (iFound[j]==0)
		{
			cout << "Species not available: " << mix->names[j] << endl;
		}
	}

	// 2. Inert Species
	// --------------------------------------------------------------------------
	for(i=1;i<=Np;i++)
	{
		double sum = 0.;
		for(j=1;j<=NC;j++)
			if(j!=data->jINERT) sum+=X[i][j];
		X[i][data->jINERT] = 1.-sum;
	}

	// 3. Mass Fractions
	// --------------------------------------------------------------------------
	MassFractionsAndPMtot();

	cout << T.Max() << endl;
}

void OpenSMOKE_Flamelet::Initial_conditions_soot_module()
{
	for (int i=1;i<=Np;i++)
	{
		sootModel.initial_values(1.0);	// TODO
		phiN[i] = sootModel.phiNStart;
		phiM[i] = sootModel.phiMStart;
	}
}

void OpenSMOKE_Flamelet::prepare_radiation()
{
	Tenv4   = BzzPow4(295.0);								// [K4]

	iH2O = mix->recognize_species_without_exit("H2O");		// [-]
	iCO2 = mix->recognize_species_without_exit("CO2");		// [-]
	iCO  = mix->recognize_species_without_exit("CO");		// [-]
	iCH4 = mix->recognize_species_without_exit("CH4");		// [-]
}

void OpenSMOKE_Flamelet::calculate_radiation()
{
	double uT;
	double K_H2O, K_CO2, K_CO, K_CH4;
	BzzVector as(4);

	for(int i=1; i<=Np; i++)
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

		// Absorption coefficients
		if(iH2O!=0)	as[1] = K_H2O*X[i][iH2O];		// [1/m]
		if(iCO2!=0)	as[2] = K_CO2*X[i][iCO2];		// [1/m]
		if(iCO!=0)	as[3] = K_CO*X[i][iCO];			// [1/m]
		if(iCH4!=0)	as[4] = K_CH4*X[i][iCH4];		// [1/m]

		asTot[i] = as.GetSumElements();				// Absorption Coefficient

		// Source term
		Qrad[i] = - 4.*Constants::sigma * asTot[i] * (BzzPow4(T[i]) - Tenv4); // Source term [W/m3]
	}
}


void OpenSMOKE_Flamelet::Run()
{
	if (data->iBackUp == false)										// Start from scratch
		Setup();

	for(int k=1;k<=operations->nOperations;k++)
	{
		if (operations->iOperation[k]==1)							// ODE System
			Solve(operations->iOptionB[k]);

		else if (operations->iOperation[k]==111)					// Add points
			NewPoints("TEMPERATURE", operations->iOptionA[k]);
		
		else if (operations->iOperation[k]==112)					// Add points
			NewPoints("QREACTION", operations->iOptionA[k]);

		else if (operations->iOperation[k]==12)						// Double grid
			DoubleTheGrid();
	
		else if (operations->iOperation[k]==13)						// Refine peak
			RefineGridPeak(operations->iOptionB[k]);

		else if (operations->iOperation[k]==14)						// Refine lean side
			RefineLeanSide(operations->iOptionB[k]);

	}

	// Print Final Solution
	GnuPlot_label(nameFileOutput, 0);
	PrintGnuPlot();
}

void OpenSMOKE_Flamelet::RunMemorySaved()
{
	Solve(1.e0);
	Solve(1.e2);

	double Told = T.Max();
	Solve(1.e4);
	for(int j=1;j<=3;j++)
	{
		if (fabs(T.Max()-Told) > 0.001)
		{
			Told = T.Max();
			Solve(1.e5 * double(j));
		}
		else
			break;
	}
	
	// Print Final Solution
	GnuPlot_label(nameFileOutput, 0);
	PrintGnuPlot();
}

void OpenSMOKE_Flamelet::ReRun()
{
	data->iBackUp = true;

	Run();
}

void OpenSMOKE_Flamelet::FoldersAndFilesManager()
{
	std::string MSDOScommand;

	// 1. Output Folder name
	nameOutputFolder = "Flamelet_";	
	nameOutputFolder += mix->name_kinetic_scheme;

	MSDOScommand = "mkdir " + nameOutputFolder;
	system(MSDOScommand.c_str());

	// 2. Backup Input data file name
	nameFolderBackupData		= nameOutputFolder;
	nameFolderUnsteadyData		= nameOutputFolder;
	nameFolderSteadyData		= nameOutputFolder;
	nameFolderAdditionalData	= nameOutputFolder;
	nameFolderProfilesData		= nameOutputFolder;
	#if LINUX_SO==1
		nameFolderBackupData     += "/BackUp";
		nameFolderUnsteadyData   += "/Unsteady";
		nameFolderSteadyData     += "/Steady";
		nameFolderAdditionalData += "/Additional";
		nameFolderProfilesData   += "/Profiles";
		nameFileOutput			  = nameFolderProfilesData + "/Final_Solution.out";
	#else
		nameFolderBackupData     += "\\BackUp";
		nameFolderUnsteadyData   += "\\Unsteady";
		nameFolderSteadyData     += "\\Steady";
		nameFolderAdditionalData += "\\Additional";
		nameFolderProfilesData   += "\\Profiles";
		nameFileOutput			  = nameFolderProfilesData + "\\Final_Solution.out";
	#endif

	MSDOScommand = "mkdir " + nameFolderBackupData;
	system(MSDOScommand.c_str());

	MSDOScommand = "mkdir " + nameFolderUnsteadyData;
	system(MSDOScommand.c_str());

	MSDOScommand = "mkdir " + nameFolderSteadyData;
	system(MSDOScommand.c_str());

	MSDOScommand = "mkdir " + nameFolderAdditionalData;
	system(MSDOScommand.c_str());

	MSDOScommand = "mkdir " + nameFolderProfilesData;
	system(MSDOScommand.c_str());
	
	nameFileBackupInputData = nameFolderBackupData + "/BackUp.inp";
}

void OpenSMOKE_Flamelet::FoldersAndFilesManager(const std::string backupFolder)
{
	// 1. Output Folder name
	nameOutputFolder = backupFolder;

	// 2. BackUp Folder name
	nameFolderBackupData		=	nameOutputFolder;
	nameFolderUnsteadyData		=	nameOutputFolder;
	nameFolderSteadyData		=	nameOutputFolder;
	nameFolderProfilesData		=	nameOutputFolder;
	nameFolderAdditionalData	=	nameOutputFolder;

	#if LINUX_SO==1
		nameFolderBackupData		+= "/BackUp";
		nameFileBackupInputData		 =	nameFolderBackupData;
		nameFileBackupInputData		+= "/BackUp.inp";
		nameFolderUnsteadyData		+= "/Unsteady";
		nameFolderSteadyData		+= "/Steady";
		nameFolderAdditionalData	+= "/Additional";
		nameFolderProfilesData	    += "/Profiles";
		nameFileOutput			     = "Final_Solution.out";
	#else
		nameFolderBackupData		+= "\\BackUp";
		nameFileBackupInputData		 = nameFolderBackupData;
		nameFileBackupInputData		+= "\\BackUp.inp";
		nameFolderUnsteadyData		+= "\\Unsteady";
		nameFolderSteadyData		+= "\\Steady";
		nameFolderAdditionalData	+= "\\Additional";
		nameFolderProfilesData	    += "\\Profiles";
		nameFileOutput			     = "\\Final_Solution.out";
	#endif
}


void OpenSMOKE_Flamelet::PasteFromExternalSolution(OpenSMOKE_Flamelet_Solution &solution)
{
	ptFlamelet = this;
	
	if (Np != solution.Np)
		ErrorMessage("External Solution must have the same dimension of destination solution");

	T = solution.T;
	w = solution.W;
			
	if (data->i2E == true)
	{
		phiN = solution.phiN;
		phiM = solution.phiM;
	}

	MoleFractionsAndPMtot();
}

void OpenSMOKE_Flamelet_Solution::PasteFromExternalSolution(OpenSMOKE_Flamelet &flame)
{
	Np = flame.Np;
	x  = flame.grid.x;
	T  = flame.T;
	W  = flame.w;
			
	if (flame.data->i2E == true)
	{
		phiN = flame.phiN;
		phiM = flame.phiM;
	}
}

void OpenSMOKE_Flamelet::ElementalAnalysis()
{	
	// Elemental analysis
	mix->GetElementalMoleFractionsFromSpeciesMoleFractions(x_elemental, X);
	mix->GetElementalMassFractionsFromSpeciesMassFractions(omega_elemental, w);
}
