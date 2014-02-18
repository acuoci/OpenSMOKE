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

#include <iomanip>
#include "droplet/OpenSMOKE_DropletMicrogravity.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "droplet/OpenSMOKE_DropletMicrogravity_DataManager.h"
#include "OpenSMOKE_DropletMicrogravity_Utilities.hpp"
#include "droplet/OpenSMOKE_Droplet.h"

void OpenSMOKE_DropletMicrogravity::PrintBackup()
{
	cout << " * Writing Backup file..." << endl;
	BzzSave fBackup("Backup.out");
	fBackup << NCGas;
	for (int j=1;j<=NCGas;j++)
		fBackup << mix->names[j].c_str();
	fBackup << grid.x;
	fBackup << m;
	fBackup << T;
	fBackup << OmegaGas;
	fBackup.End();
}

void OpenSMOKE_DropletMicrogravity::DAE_myPrint(BzzVector &x, double t)
{
	if(data->iterationVideoCounter == data->nStepsVideo)
	{
		if (data->iMode == EIGENVALUE)
		{
			cout << data->iteration													<< "\t";
			cout << t																<< "\t";
			cout << 4./Constants::pi*m[1]/InterfaceLiquidDensity()/(2.*grid.x[1]) * 1e6		<< "\t";
			cout << T[1]															<< "\t";
			cout << OmegaGas[1][data->jFuel[1]]										<< "\t";
			if (NCLiquid==2)
				cout << OmegaGas[1][data->jFuel[2]]									<< "\t";
			cout << T.Max()															<< "\t";
			cout << mass[1]*1000.													<< "\t";
			cout << endl;
		}
		else if (data->iMode == UNSTEADY_BATCH)
		{
			cout << data->iteration													<< "\t";
			cout << t																<< "\t";
			cout << 4./Constants::pi*m[1]/InterfaceLiquidDensity()/(2.*grid.x[1]) * 1e6		<< "\t";
			cout << T[1]															<< "\t";
			cout << OmegaGas[1][data->jFuel[1]]										<< "\t";
			if (NCLiquid == 2)
				cout << OmegaGas[1][data->jFuel[2]]									<< "\t";
			cout << T.Max()															<< "\t";
			cout << grid.x[1]*2.*1000.												<< "\t";

			if (data->enhancingFactor > 0.)
				cout << YGlycerol(t, data->enhancingFactor);

			cout << endl;
		}

		data->iterationVideoCounter = 0;
	}

	if (data->iterationUnsteadyFileCounter == 50)
	{
		// Unsteady maps
		{
			for(int i=1;i<=N;i++)
			{
				fUnsteady << setw(16) << left << t;
				fUnsteady << setw(16) << left << grid.x[i]*1000.;
				fUnsteady << setw(16) << left << grid.x[i]/grid.x[1];
				fUnsteady << setw(16) << left << T[i];
				fUnsteady << setw(16) << left << u[i];
				for(int j=1;j<=NCLiquid;j++)
					fUnsteady << setw(16) << left << XGas[i][data->jFuel[j]];
				fUnsteady << setw(16) << left << XGas[i][data->jOxidizer];
				fUnsteady << setw(16) << left << XGas[i][data->jInert];
				fUnsteady << setw(16) << left << XGas[i][data->jH2O]	;
				fUnsteady << setw(16) << left << XGas[i][data->jCO2]	;
				fUnsteady << setw(16) << left << XGas[i][data->jCO];
				fUnsteady << setw(16) << left << T_reaction[i];
				fUnsteady << setw(16) << left << endl;
			}
			fUnsteady << endl;

			if (data->iRadiation == 1)
			{
				for(int i=1;i<=N;i++)
				{
					fUnsteadyRadiation << setw(16) << left <<  t;
					fUnsteadyRadiation << setw(16) << left << grid.x[i]*1000.;
					fUnsteadyRadiation << setw(16) << left << grid.x[i]/grid.x[1];
					fUnsteadyRadiation << setw(16) << left << T[i]; // 4
					fUnsteadyRadiation << setw(16) << left << Kp[i]; // 5
					fUnsteadyRadiation << setw(16) << left << radiation.divq[i]; // 6
					fUnsteadyRadiation << setw(16) << left << radiation.G0[i]; // 7
					fUnsteadyRadiation << setw(16) << left << radiation.G1[i]; // 8
					fUnsteadyRadiation << setw(16) << left << radiation.G2[i]; // 9
					fUnsteadyRadiation << setw(16) << left << radiation.G3; // 10

					fUnsteadyRadiation << setw(16) << left << AbsorptionCoefficient_H2O(T[i])*data->P_atm*XGas[i][data->jH2O];
					fUnsteadyRadiation << setw(16) << left << AbsorptionCoefficient_CO2(T[i])*data->P_atm*XGas[i][data->jH2O];
					fUnsteadyRadiation << setw(16) << left << AbsorptionCoefficient_CO(T[i]) *data->P_atm*XGas[i][data->jH2O];
					fUnsteadyRadiation << setw(16) << left << AbsorptionCoefficient_CH4(T[i])*data->P_atm*XGas[i][data->jH2O];

					fUnsteadyRadiation << setw(16) << left << AbsorptionCoefficient_H2O(T[i])*data->P_atm;
					fUnsteadyRadiation << setw(16) << left << AbsorptionCoefficient_CO2(T[i])*data->P_atm;
					fUnsteadyRadiation << setw(16) << left << AbsorptionCoefficient_CO(T[i]) *data->P_atm;
					fUnsteadyRadiation << setw(16) << left << AbsorptionCoefficient_CH4(T[i])*data->P_atm;

					fUnsteadyRadiation << endl;
				}
				fUnsteadyRadiation << endl;
			}
		}

		data->iterationUnsteadyFileCounter = 0;
	}
	
	if(data->iterationFileCounter == data->nStepsFile)
	{
		// Interface
		{
			if (data->iteration == 1)		
				PrintTagInterfaceFile();

			fInterface << setw(16) << left << t; 
			fInterface << setw(16) << left << T[1]; 
			fInterface << setw(16) << left << u[1]*1000.;
			fInterface << setw(16) << left << grid.x[1]*2.*1000.;
			fInterface << setw(16) << left << mass[1]*1000.;
			fInterface << setw(16) << left << BzzPow2(grid.x[1]*2.*1000.);
			fInterface << setw(16) << left << 2./Constants::pi*m[1]/InterfaceLiquidDensity()/BzzPow2(2.*grid.x[1]) * 1e3;
			fInterface << setw(16) << left << 4./Constants::pi*m[1]/InterfaceLiquidDensity()/(2.*grid.x[1]) * 1e6;
			fInterface << setw(16) << left << m[1]*1000.;
			fInterface << setw(16) << left << radiation.qSurface;
			fInterface << setw(16) << left << radiation.qSurface*grid.A[1];
			fInterface << setw(16) << left << radiation.G3;

			for(int j=1;j<=NCLiquid;j++)
				fInterface << setw(16) << left << data->OmegaDroplet[data->jFuel[j]];
			for(int j=1;j<=NCGas;j++)
				fInterface << setw(16) << left << XGas[1][j];

			fInterface << endl;
		}

		// Additional
		{
			if (data->iteration == 1)		
				PrintTagFlameFile();

			double xTMax, TMax;
			FindMaximumUsingSecondOrderPolynomial(grid.x, T, xTMax, TMax);

			double xOHMax, OHMax;
			int jOH = mix->recognize_species("OH");
			BzzVector auxOH; auxOH= XGas.GetColumn(jOH);
			FindMaximumUsingSecondOrderPolynomial(grid.x, auxOH, xOHMax, OHMax);

			double xC2H2Max, C2H2Max;
			int jC2H2 = mix->recognize_species("C2H2");
			BzzVector auxC2H2 = XGas.GetColumn(jC2H2);
			FindMaximumUsingSecondOrderPolynomial(grid.x, auxC2H2, xC2H2Max, C2H2Max);

			double flameWidthT  = FindWidth(grid.x, T);
			double flameWidthOH = FindWidth(grid.x, auxOH);
			double mSoot = SootPropensity();

			fFlame << setw(20) << left << t; 
			fFlame << setw(20) << left << T[1]; 

			fFlame << setw(20) << left << TMax;
			fFlame << setw(20) << left << xTMax*1000.*2;
			fFlame << setw(20) << left << xTMax/grid.x[1];

			fFlame << setw(20) << left << OHMax;
			fFlame << setw(20) << left << xOHMax*1000.*2.;
			fFlame << setw(20) << left << xOHMax/grid.x[1];

			fFlame << setw(20) << left << C2H2Max;
			fFlame << setw(20) << left << xC2H2Max*1000.*2.;
			fFlame << setw(20) << left << xC2H2Max/grid.x[1];

			fFlame << setw(20) << left << flameWidthT;
			fFlame << setw(20) << left << flameWidthT/(2.*grid.x[1]);
			fFlame << setw(20) << left << flameWidthOH;
			fFlame << setw(20) << left << flameWidthOH/(2.*grid.x[1]);

			fFlame << setw(20) << left << mSoot/m[1];

			fFlame << endl;
		}

		data->iterationFileCounter = 0;
	}

	if(data->iterationBackUpCounter == data->nStepsBackUp)
	{
		PrintBackup();

		data->iterationBackUpCounter = 0;
	}

	if (data->iLiquidPhase != LIQUID_DROPLET_PERFECTLY_STIRRED)
	{
		if (data->iteration <= 100)	droplet->T = T[1];
		droplet->UpdateFromGasPhase(tOld, t);
//		droplet->UpdateFromGasPhase(m[1], lambda[1] * dT_over_dr[1] * grid.A[1], radiation.qSurface* grid.A[1]);
		droplet->UpdateFromGasPhase(T[1]);
		droplet->Run();
	//	getchar();
	}

	// Update Grid
	if (data->iMode == UNSTEADY_BATCH)
		grid.Construct(N, data->rEnvironment - DropletDiameter()/2., data->stretchingFactor, DropletDiameter()/2.);

	// Update radiative heat transfer [W/m3]
	if (data->iRadiation == 1)
		radiation.Calculate(grid.x, T, Kp);

	// Update counters and time
	data->iteration++;
	data->iterationVideoCounter++;
	data->iterationFileCounter++;
	data->iterationBackUpCounter++;
	data->iterationUnsteadyFileCounter++;
	tOld = t;
}

void OpenSMOKE_DropletMicrogravity::PrintFinalSolution()
{
	ofstream fFinalSolution;
	openOutputFileAndControl(fFinalSolution, data->nameOutputFolder + "/FinalSolution.out");
	fFinalSolution.setf(ios::scientific);

	int count = 1;
	PrintTagOnGnuplotLabel(20, fFinalSolution, "r[mm]",			count);
	PrintTagOnGnuplotLabel(20, fFinalSolution, "csi[-]",		count);
	PrintTagOnGnuplotLabel(20, fFinalSolution, "m[g/s]",		count);
	PrintTagOnGnuplotLabel(20, fFinalSolution, "dD[mm/s]",		count);
	PrintTagOnGnuplotLabel(20, fFinalSolution, "dD2[mm2/s]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolution, "T[K]",			count);
	PrintTagOnGnuplotLabel(20, fFinalSolution, "u[mm/s]",		count);
	PrintTagOnGnuplotLabel(20, fFinalSolution, "divq[W/m3]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolution, "Kp[1/m]",		count);
	for(int j=1;j<=NCGas;j++)
		PrintTagOnGnuplotLabel(20, fFinalSolution, "x_" + mix->names[j],count);
	for(int j=1;j<=NCGas;j++)
		PrintTagOnGnuplotLabel(20, fFinalSolution, "w_" + mix->names[j],count);
	fFinalSolution << endl;
	fFinalSolution << endl;

	for(int i=1;i<=N;i++)
	{
		fFinalSolution	<< setw(20) << left << grid.x[i]*1.e3;
		fFinalSolution	<< setw(20) << left << grid.x[i]/grid.x[1];
		fFinalSolution	<< setw(20) << left << m[i]*1e3;
		fFinalSolution	<< setw(20) << left << 2.*(m[i]*1e3)/Constants::pi/InterfaceLiquidDensity()/BzzPow2(2.*grid.x[1]);
		fFinalSolution	<< setw(20) << left << 4.*(m[i]*1e6)/Constants::pi/InterfaceLiquidDensity()/(2.*grid.x[1]);
		fFinalSolution	<< setw(20) << left << T[i];
		fFinalSolution	<< setw(20) << left << u[i]*1e3;
		fFinalSolution	<< setw(20) << left << Kp[i];
		if (data->iRadiation == 1)
			fFinalSolution	<< setw(20) << left << radiation.divq[i];
		else
			fFinalSolution	<< setw(20) << left << 0;

		for(int j=1;j<=NCGas;j++)
			fFinalSolution	<< setw(20) << left << XGas[i][j];
		for(int j=1;j<=NCGas;j++)
			fFinalSolution	<< setw(20) << left << OmegaGas[i][j];
		fFinalSolution << endl;
	}

	fFinalSolution.close();

	// Print Backup
	PrintBackup();
}


void OpenSMOKE_DropletMicrogravity::PrintFinalSummary()
{
	ofstream fFinalSummary;
	openOutputFileAndControl(fFinalSummary, data->nameOutputFolder + "/FinalSummary.out");
	fFinalSummary.setf(ios::scientific);

	double xTMax, TMax;
	FindMaximumUsingSecondOrderPolynomial(grid.x, T, xTMax, TMax);

	double xOHMax, OHMax;
	int jOH = mix->recognize_species("OH");
	BzzVector auxOH; auxOH= XGas.GetColumn(jOH);
	FindMaximumUsingSecondOrderPolynomial(grid.x, auxOH, xOHMax, OHMax);

	double xC2H2Max, C2H2Max;
	int jC2H2 = mix->recognize_species("C2H2");
	BzzVector auxC2H2 = XGas.GetColumn(jC2H2);
	FindMaximumUsingSecondOrderPolynomial(grid.x, auxC2H2, xC2H2Max, C2H2Max);

	double flameWidthT  = FindWidth(grid.x, T);
	double flameWidthOH = FindWidth(grid.x, auxOH);

	double mSoot = SootPropensity();
	

	fFinalSummary << "Droplet_diameter_[mm]:           " << 2.*grid.x[1]*1000. << endl;
	fFinalSummary << "Droplet_surface_[mm2]:           " << grid.A[1]*1.e6 << endl;
	fFinalSummary << "Droplet_volume_[mm3]:            " << Constants::pi/6.*BzzPow3(2.*grid.x[1])*1.e9 << endl;
	fFinalSummary << "Droplet_mass_[g]:                " << DropletMass()*1000. << endl;
	fFinalSummary << "Droplet_surface_temperature_[K]: " << T[1] << endl;
	fFinalSummary << "Environment_O2_mass_fraction:    " << data->OmegaEnvironment[data->jOxidizer] << endl;
	fFinalSummary << "Environment_O2_temperature_[K]:  " << data->TEnvironment << endl;

	fFinalSummary << "Vaporization_rate_[g/s]:         " << m[1]*1000. << endl;
	fFinalSummary << "Vaporization_rate_[mm/s]:        " << 2.*(m[1]*1e3)/Constants::pi/InterfaceLiquidDensity()/BzzPow2(2.*grid.x[1]) << endl;
	fFinalSummary << "Vaporization_rate_[mm2/s]:       " << 4.*(m[1]*1e6)/Constants::pi/InterfaceLiquidDensity()/(2.*grid.x[1]) << endl;

	fFinalSummary << "Flame_temperature_[K]:           " << TMax << endl;
	fFinalSummary << "Flame_position_[mm]:             " << xTMax*1000. << endl;
	fFinalSummary << "Flame_stand-off_ratio_[-]:       " << xTMax/grid.x[1] << endl;
	fFinalSummary << "Flame_width_[mm]:                " << flameWidthT*1000. << endl;
	fFinalSummary << "Flame_width_(width/D)[-]:        " << flameWidthT/(2.*grid.x[1]) << endl;

	fFinalSummary << "Flame_OH_peak_(mole_fraction):   " << OHMax << endl;
	fFinalSummary << "Flame_position_(OH):             " << xOHMax*1000. << endl;
	fFinalSummary << "Flame_stand-off_ratio_(OH):      " << xOHMax/grid.x[1] << endl;
	fFinalSummary << "Flame_width_[mm]:                " << flameWidthOH*1000. << endl;
	fFinalSummary << "Flame_width_(width/D)[-]:        " << flameWidthOH/(2.*grid.x[1]) << endl;

	fFinalSummary << "Flame_C2H2_peak_(mole_fraction): " << C2H2Max << endl;
	fFinalSummary << "Mass_C2H2_[g]:                   " << MassIntegral("C2H2") << endl;
	fFinalSummary << "Mass_C2H2_over_mass_droplet[-]:  " << MassIntegral("C2H2")/DropletMass() << endl;
	fFinalSummary << "Soot_propensity[-]:              " << mSoot/m[1] << endl;
	fFinalSummary << "Volume[mm3]:                     " << MassIntegral("volume")*1e9 << endl;
	
	
	fFinalSummary.close();
}
void OpenSMOKE_DropletMicrogravity::PrintFinalSolutionProperties()
{
	BzzVector sumDiffusionFluxes(N);
	for(int i=1;i<=N;i++)
		for(int j=1;j<=NCGas;j++)
			sumDiffusionFluxes[i] += vStar_e[i][j];

	ofstream fFinalSolutionProperties;
	openOutputFileAndControl(fFinalSolutionProperties, data->nameOutputFolder + "/FinalSolution.prop.out");
	fFinalSolutionProperties.setf(ios::scientific);

	int count = 1;
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "r[mm]",			count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "csi[-]",		count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "rho[kg/m3]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "MW[kg/kmol]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "Cp[J/kg/K]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "lamb[W/m/K]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "T_Con[W/m3]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "T_Dif[W/m3]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "T_Cp[W/m3]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "T_Rea[W/m3]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "Sum_V[W/m3]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "V*O2[m/s]",		count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "DmO2[m2/s]",	count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "SumOmega",		count);
	PrintTagOnGnuplotLabel(20, fFinalSolutionProperties, "SumDOmegaDt",	count);
	fFinalSolutionProperties << endl;
	fFinalSolutionProperties << endl;

	for(int i=1;i<=N;i++)
	{
		fFinalSolutionProperties << setw(20) << left << grid.x[i]*1.e3;
		fFinalSolutionProperties << setw(20) << left << grid.x[i]/grid.x[1];
		fFinalSolutionProperties << setw(20) << left << rho[i];
		fFinalSolutionProperties << setw(20) << left << MW[i];
		fFinalSolutionProperties << setw(20) << left << Cp[i];
		fFinalSolutionProperties << setw(20) << left << lambda[i];
		fFinalSolutionProperties << setw(20) << left << T_convection[i];
		fFinalSolutionProperties << setw(20) << left << T_conduction[i];
		fFinalSolutionProperties << setw(20) << left << T_diffusion_fluxes[i];
		fFinalSolutionProperties << setw(20) << left << T_reaction[i];
		fFinalSolutionProperties << setw(20) << left << sumDiffusionFluxes[i];
		fFinalSolutionProperties << setw(20) << left << vStar_e[i][data->jOxidizer];
		fFinalSolutionProperties << setw(20) << left << Dm[i][data->jOxidizer];
		fFinalSolutionProperties << setw(20) << left << OmegaGas.GetRow(i).GetSumElements();
		fFinalSolutionProperties << setw(20) << left << dOmegaGas_over_dt.GetRow(i).GetSumElements();
		fFinalSolutionProperties << endl;
	}

	fFinalSolutionProperties.close();
}

void OpenSMOKE_DropletMicrogravity::PrintTagInterfaceFile()
{
	int count = 1;
	PrintTagOnGnuplotLabel(16, fInterface, "time[s]",		count);
	PrintTagOnGnuplotLabel(16, fInterface, "T[K]",			count);
	PrintTagOnGnuplotLabel(16, fInterface, "u[mm/s]",		count);
	PrintTagOnGnuplotLabel(16, fInterface, "D[mm]",			count);
	PrintTagOnGnuplotLabel(16, fInterface, "m[g]",			count);
	PrintTagOnGnuplotLabel(16, fInterface, "D2[mm2]",		count);
	PrintTagOnGnuplotLabel(16, fInterface, "dD[mm/s]",		count);
	PrintTagOnGnuplotLabel(16, fInterface, "dD2[mm2/s]",	count);
	PrintTagOnGnuplotLabel(16, fInterface, "dm[g/s]",		count);
	PrintTagOnGnuplotLabel(16, fInterface, "qSur[W/m2]",	count);
	PrintTagOnGnuplotLabel(16, fInterface, "QSur[W]",		count);
	PrintTagOnGnuplotLabel(16, fInterface, "G3[-]",			count);
	for(int j=1;j<=NCLiquid;j++)
		PrintTagOnGnuplotLabel(16, fInterface, "wL_" + data->nameFuels[j-1],count);
	for(int j=1;j<=NCGas;j++)
		PrintTagOnGnuplotLabel(16, fInterface, "xG_" + mix->names[j],count);
	fInterface << endl;
	fInterface << endl;
}

void OpenSMOKE_DropletMicrogravity::PrintTagFlameFile()
{
	int count = 1;

	PrintTagOnGnuplotLabel(20, fFlame, "time[s]",			count);
	PrintTagOnGnuplotLabel(20, fFlame, "Ti[K]",			count);

	PrintTagOnGnuplotLabel(20, fFlame, "Tmax[K]",			count);
	PrintTagOnGnuplotLabel(20, fFlame, "Df(T)[mm]",		count);
	PrintTagOnGnuplotLabel(20, fFlame, "Df(T)[-]",			count);

	PrintTagOnGnuplotLabel(20, fFlame, "OHmax[K]",			count);
	PrintTagOnGnuplotLabel(20, fFlame, "Df(T)[mm]",		count);
	PrintTagOnGnuplotLabel(20, fFlame, "Df(T)[-]",			count);

	PrintTagOnGnuplotLabel(20, fFlame, "C2H2max[K]",		count);
	PrintTagOnGnuplotLabel(20, fFlame, "Df(C2H2)[mm]",		count);
	PrintTagOnGnuplotLabel(20, fFlame, "Df(OH)[-]",		count);

	PrintTagOnGnuplotLabel(20, fFlame, "width(T)[mm]",		count);
	PrintTagOnGnuplotLabel(20, fFlame, "width(T)[-]",		count);

	PrintTagOnGnuplotLabel(20, fFlame, "width(OH)[mm]",	count);
	PrintTagOnGnuplotLabel(20, fFlame, "width(OH)[-]",		count);

	PrintTagOnGnuplotLabel(20, fFlame, "soot_prop",		count);
	fFlame << endl;
	fFlame << endl;
}
