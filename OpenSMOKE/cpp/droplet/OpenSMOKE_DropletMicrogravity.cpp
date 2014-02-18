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
double YGlycerol(const double t, const double enhancing_factor);
double XGlycerol(const double t, const double enhancing_factor);
#include "OpenSMOKE_DropletMicrogravity_Utilities.hpp"
#include "OpenSMOKE_DropletMicrogravity_Output.hpp"

double YGlycerol(const double t, const double enhancing_factor)
{
	BzzVector x1(15, 0., 1.52E-03, 1.03E-01,2.35E-01,2.97E-01,3.26E-01,3.50E-01,3.86E-01,4.41E-01,5.12E-01,6.27E-01, 1.02E+00,1.69E+00, 1.99E+00,2.00E+00);
	BzzVector y1(15, 0.2, 2.08E-01, 4.47E-01 ,7.65E-01, 9.21E-01, 9.60E-01, 9.78E-01, 9.88E-01, 9.93E-01, 9.96E-01, 9.96E-01, 9.97E-01, 9.99E-01, 1.00E+00, 1.00E+00);

	BzzVector x2(21, 0., 3.03E-03, 6.52E-02, 1.58E-01, 2.33E-01, 3.21E-01, 4.29E-01, 5.33E-01, 6.29E-01, 7.32E-01, 8.14E-01, 8.95E-01, 9.52E-01, 1.01E+00, 1.05E+00, 1.10E+00, 1.17E+00, 1.25E+00, 1.38E+00, 1.75E+00, 2.00E+00);
	BzzVector y2(21, 0.2, 2.06E-01, 2.31E-01, 2.73E-01, 3.21E-01, 3.85E-01, 4.61E-01, 5.39E-01, 6.20E-01, 7.12E-01, 7.94E-01, 8.80E-01, 9.34E-01, 9.63E-01, 9.78E-01, 9.89E-01, 9.93E-01, 9.96E-01, 9.99E-01, 9.99E-01, 1.00E+00);

	LinearInterpolation lin1;
	LinearInterpolation lin2;
	lin1(x1,y1);
	lin2(x2,y2);

	double z1 = lin1(t);
	double z2 = lin2(t);
	double w = 1./9*enhancing_factor - 1./9.;

	return (1.-w)*z1 + w*z2;
}

double XGlycerol(const double t, const double enhancing_factor)
{
	double z = YGlycerol(t, enhancing_factor);
	double pm = 1./(z/92.093+(1.-z)/60.10);

	return z*pm/92.093;
}



OpenSMOKE_DropletMicrogravity *ptDroplet;

void DAE_Print(BzzVector &x, double t)
{
	ptDroplet->DAE_myPrint(x, t);
}

const double OpenSMOKE_DropletMicrogravity::MMIN  =  -1.e6;

void OpenSMOKE_DropletMicrogravity::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_DropletMicrogravity"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_DropletMicrogravity::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_DropletMicrogravity"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

void OpenSMOKE_DropletMicrogravity::MemoryAllocation()
{
	ChangeDimensions(N, &m);
	ChangeDimensions(N, &dm_over_dt);
	ChangeDimensions(N, &du_over_dr);

	ChangeDimensions(N, &mass);
	ChangeDimensions(N, &dmass_over_dt);

	ChangeDimensions(N, &T);
	ChangeDimensions(N, &dT_over_dt);
	ChangeDimensions(N, &dT_over_dr);
	ChangeDimensions(N, &dT_over_dr_central);

	ChangeDimensions(N, &u);

	ChangeDimensions(N, &Kp);
	ChangeDimensions(N, &lambda);
	ChangeDimensions(N, &Cp);
	ChangeDimensions(N, &rho);
	ChangeDimensions(N, &QReaction);
	ChangeDimensions(N, NCGas, &RGas);

	ChangeDimensions(N, &A_x_rhoe);
	ChangeDimensions(N, &A_x_rhow);
	ChangeDimensions(N, &A_x_lambdae);
	ChangeDimensions(N, &A_x_lambdaw);

	ChangeDimensions(N, &T_convection);
	ChangeDimensions(N, &T_conduction);
	ChangeDimensions(N, &T_diffusion_fluxes);
	ChangeDimensions(N, &T_reaction);

	ChangeDimensions(N, &MW);
	ChangeDimensions(N, &uMW);
	
	ChangeDimensions(N, NCGas, &OmegaGas_diffusion);
	ChangeDimensions(N, NCGas, &OmegaGas_convection);
	ChangeDimensions(N, NCGas, &OmegaGas_reaction);

	ChangeDimensions(N, NCGas, &OmegaGas);
	ChangeDimensions(N, NCGas, &dOmegaGas_over_dt);
	ChangeDimensions(N, NCGas, &dOmegaGas_over_dr);

	ChangeDimensions(N, NCGas, &XGas);
	ChangeDimensions(N, NCGas, &Dm);
	ChangeDimensions(N, NCGas, &Teta);
	ChangeDimensions(N, NCGas, &Cpk);

	ChangeDimensions(N, &vc_e);
	ChangeDimensions(N, NCGas, &vStar_e);
	ChangeDimensions(N, NCGas, &vStar_w);
	ChangeDimensions(N, NCGas, &coeff_e);
}

void OpenSMOKE_DropletMicrogravity::MemoryAllocationKinetics()
{
	ChangeDimensions(NCGas, &OmegaGasVector);
	ChangeDimensions(NCGas, &XGasVector);
	ChangeDimensions(NCGas, &CGasVector);
	ChangeDimensions(NCGas, &RGasVector);
	ChangeDimensions(NCGas, &DmixVector);
	ChangeDimensions(NCGas, &TetaMixVector);
}

void OpenSMOKE_DropletMicrogravity::PrepareSystem()
{
	if (data->iMode == EIGENVALUE)			dimensionBlock	= 2 + NCGas;
	else if (data->iMode == UNSTEADY_BATCH)	dimensionBlock	= 3 + NCGas;
	
	// Number of equations
	NE = N*dimensionBlock;

	// Memory allocation
	ChangeDimensions(NE, &xMin);
	ChangeDimensions(NE, &xMax);
	ChangeDimensions(NE, &inDerAlg);
	ChangeDimensions(NE, &xFirstGuess);

	// DAE system preparation
	setMinimumAndMaximumValues();
	setDifferentialAndAlgebraic();
	setFirstGuess();
}	

OpenSMOKE_DropletMicrogravity::OpenSMOKE_DropletMicrogravity()
{
	cout.setf(ios::scientific);

	ptDroplet = this;
	grid.SetSpherical();
	droplet = new OpenSMOKE_Droplet();

	tOld = 0.;
}

void OpenSMOKE_DropletMicrogravity::SetName(const string name)
{
	name_object = name;
}

void OpenSMOKE_DropletMicrogravity::Assign(OpenSMOKE_ReactingGas *_mix)
{
	mix = _mix;
}

void OpenSMOKE_DropletMicrogravity::Assign(OpenSMOKE_DropletMicrogravity_DataManager *_data)
{
	data = _data;
}

void OpenSMOKE_DropletMicrogravity::Assign(OpenSMOKE_GlobalKinetics *_global)
{
	global = _global;
}

double OpenSMOKE_DropletMicrogravity::InterfaceLiquidSpecificHeat()
{
	double sum=0.;
	for(int j=1;j<=NCLiquid;j++)
		sum += cpDroplet[j] * data->OmegaDroplet[data->jFuel[j]];
	return sum;
}

double OpenSMOKE_DropletMicrogravity::InterfaceLiquidDensity()
{
	double sum=0.;
	for(int j=1;j<=NCLiquid;j++)
		sum += rhoDroplet[j] * data->OmegaDroplet[data->jFuel[j]];
	return sum;
}

double OpenSMOKE_DropletMicrogravity::InterfaceLiquidVaporizationHeat(const double t)
{
	if (data->enhancingFactor > 0.)
	{
		double yglycerol = YGlycerol(t, data->enhancingFactor);
		return dHVaporization[1]*(1.-yglycerol) + dHVaporization[2]*yglycerol;
	}
	else
	{
		double sum=0.;
		for(int j=1;j<=NCLiquid;j++)
			sum += dHVaporization[j] * data->OmegaDroplet[data->jFuel[j]];
		return sum;
	}
}

double OpenSMOKE_DropletMicrogravity::InterfaceLiquidMolecularWeight()
{
	double sum=0.;
	for(int j=1;j<=NCLiquid;j++)
		sum += data->OmegaDroplet[data->jFuel[j]] / mix->M[data->jFuel[j]];
	return 1./sum;
}

double OpenSMOKE_DropletMicrogravity::DropletMass()
{
	if (data->iMode == EIGENVALUE)
		return InterfaceLiquidDensity()*DropletVolume();
	else if (data->iMode == UNSTEADY_BATCH)
		return InterfaceLiquidDensity()*DropletVolume();
	else
		return 0;
}

double OpenSMOKE_DropletMicrogravity::DropletDiameter()
{
	if (data->iMode == EIGENVALUE)
		return data->diameterDroplet0;
	else if (data->iMode == UNSTEADY_BATCH)
		return pow(6.*mass[1]/Constants::pi/InterfaceLiquidDensity(), 1./3.);
	else
		return 0;
}

double OpenSMOKE_DropletMicrogravity::DropletVolume()
{
	return Constants::pi/6.*BzzPow3(2.*grid.x[1]);
}

void OpenSMOKE_DropletMicrogravity::Setup()
{
	// Output files
	{
		string MSDOScommand = "mkdir " + data->nameOutputFolder;
		system(MSDOScommand.c_str());

		openOutputFileAndControl(fUnsteady, data->nameOutputFolder + "/SolutionUnsteady.out");
		fUnsteady.setf(ios::scientific);
	
		openOutputFileAndControl(fInterface, data->nameOutputFolder + "/Interface.out");
		fInterface.setf(ios::scientific);
	
		openOutputFileAndControl(fFlame, data->nameOutputFolder + "/Flame.out");
		fFlame.setf(ios::scientific);

		if (data->iRadiation == 1)
		{
			openOutputFileAndControl(fUnsteadyRadiation, data->nameOutputFolder + "/SolutionUnsteadyRadiation.out");
			fUnsteadyRadiation.setf(ios::scientific);
		}
	}

	// Liquid droplet
	if (data->iLiquidPhase != LIQUID_DROPLET_PERFECTLY_STIRRED)
	{
		droplet->Assign(mix);
		droplet->Assign(data);
		droplet->Setup();
	}

	NCGas			= mix->NumberOfSpecies();
	NCLiquid		= data->nameFuels.size();
	N				= data->N;

	ChangeDimensions(NCLiquid, &rhoDroplet);
	ChangeDimensions(NCLiquid, &dHVaporization);
	ChangeDimensions(NCLiquid, &pvDroplet);
	ChangeDimensions(NCLiquid, &cpDroplet);
	
	ChangeDimensions(NCLiquid, &data->jFuel);
	for(int j=1;j<=NCLiquid;j++)
		data->jFuel[j] = mix->recognize_species(data->nameFuels[j-1]);
	data->jOxidizer = mix->recognize_species("O2");
	data->jInert    = mix->recognize_species("N2");
	data->jH2O    = mix->recognize_species_without_exit("H2O");
	data->jCO2    = mix->recognize_species_without_exit("CO2");
	data->jCO     = mix->recognize_species_without_exit("CO");
	data->jCH4    = mix->recognize_species_without_exit("CH4");


	data->rEnvironment = data->ratioRadii*(data->diameterDroplet0/2.);
	
	double Tinterface=data->TDroplet0;
	if (data->dropletInterfaceTemperature == INTERFACE_USER_DEFINED)		Tinterface = data->interfaceTemperature;

	for(int j=1;j<=NCLiquid;j++)
	{
		rhoDroplet[j]      = data->liquid_species[j].rho(Tinterface);
		dHVaporization[j]  = data->liquid_species[j].Hv(Tinterface);
	}

	// Adiabatic temperature and reaction heat estimation
	{
		OpenSMOKE_GasStream gas_stream;
		gas_stream.AssignKineticScheme(*mix);
		gas_stream.AssignTemperature(Tinterface, "K");
		gas_stream.AssignPressure(data->P_Pascal, "Pa");
		gas_stream.AssignMassFlowRate(1., "kg/s");
		gas_stream.AssignFuelMassFractions(data->OmegaDroplet);
		gas_stream.AssignOxidizerMassFractions(data->OmegaEnvironment);
		gas_stream.AssignEquivalenceRatio(1.0);
		gas_stream.lock();

		for(int i=1;i<=mix->NumberOfSpecies();i++)
			if (gas_stream.x[i] > 0.)	cout << mix->names[i] << " " << gas_stream.x[i] << endl;

		int flag;
		double TFlame = 2200.;
		double NFlame = 1.;
		BzzVector xFlame(mix->NumberOfSpecies());
		BzzVector xElementalFlame(mix->NumberOfElements());
		mix->GetElementalMoleFractionsFromSpeciesMoleFractions(xElementalFlame, gas_stream.x);

		mix->Equilibrium_HP(TFlame, xFlame, NFlame, gas_stream.enthalpy, data->P_Pascal, xElementalFlame, false, flag);
		
		cout << "Equilibrium composition (mole fractions)" << endl;
		cout << " * Temperature [K]: " << TFlame << endl;
		for(int i=1;i<=mix->NumberOfSpecies();i++)
			if (xFlame[i] > 1.e-8)	
				cout << " * " << setw(14) << left << mix->names[i] << setw(16) << xFlame[i] << endl;
	}

	MemoryAllocation();
	MemoryAllocationKinetics();

	// Grid construction
	grid.Construct(N, data->rEnvironment - data->diameterDroplet0/2., data->stretchingFactor, data->diameterDroplet0/2.);
		
	InitialEstimations();

	cout << "Droplet density [kg/m3]:             "	<< InterfaceLiquidDensity() << endl;
	cout << "Droplet heat of vaporization [J/kg]: "	<< InterfaceLiquidVaporizationHeat(0) << endl;
	cout << "Droplet mass [kg]:                   "	<< DropletMass() << endl;
}	

void OpenSMOKE_DropletMicrogravity::Summary()
{
	cout << "Droplet  diameter [mm]:    "	<< data->diameterDroplet0*1000. << endl;
	cout << "Droplet  surface [mm2]:    "	<< Constants::pi*BzzPow2(data->diameterDroplet0*1000.) << endl;
	cout << "Droplet  volume  [mm3]:    "	<< Constants::pi/6.*BzzPow3(data->diameterDroplet0*1000.) << endl;
	
	cout << "Environ. diameter [mm]:    "	<< data->rEnvironment*2.*1000. << endl;
	cout << "Environ. surface [mm2]:    "	<< Constants::pi*BzzPow2(data->rEnvironment*2.*1000.) << endl;
	cout << "Environ. volume  [mm3]:    "	<< Constants::pi/6.*BzzPow3(data->rEnvironment*2.*1000.) << endl;

	cout << "Radii ratio [-]:           "	<< (data->rEnvironment)/(data->diameterDroplet0/2.) << endl;

	cout << "Droplet  temperature [K]:  "	<< data->TDroplet0 << endl;
	cout << "Environ. temperature [K]:  "	<< data->TEnvironment << endl;

	cout << "Integration time [s]:      "	<< data->tEnd << endl;
}

void OpenSMOKE_DropletMicrogravity::GiveMeSpatialDerivatives()
{
	if (data->iMode == EIGENVALUE)
	{
		grid.FirstDerivative('C', u, T, dT_over_dr_central);
		grid.FirstDerivative(data->iDerT, u, T, dT_over_dr);
		grid.FirstDerivative(data->iDerW, u, OmegaGas, dOmegaGas_over_dr);
		grid.FirstDerivative(data->iDerU, u, du_over_dr);
	}
	else if (data->iMode == UNSTEADY_BATCH)
	{
		grid.FirstDerivative('C', u, T, dT_over_dr_central);
		grid.FirstDerivative(data->iDerT, u, T, dT_over_dr);
		grid.FirstDerivative(data->iDerW, u, OmegaGas, dOmegaGas_over_dr);
		grid.FirstDerivative(data->iDerU, u, du_over_dr);
	}
	else
		ErrorMessage("GiveMeSpatialDerivatives: Wrong option!");
}


void OpenSMOKE_DropletMicrogravity::GiveMeContinuityEquation(const double t)	// ConEq
{
	// Interface equations	
	if (data->iMode == EIGENVALUE)
	{
		if (data->dropletInterfaceTemperature == INTERFACE_USER_DEFINED)
			dm_over_dt[1] = T[1] - data->interfaceTemperature;
	
		else if (data->dropletInterfaceTemperature == INTERFACE_EQUILIBRIUM)
		{
		/*	if (NCLiquid == 1)
				dm_over_dt[1] = pvDroplet[1]/data->P_Pascal - OmegaGas[1][data->jFuel[1]]*MW[1]/mix->M[data->jFuel[1]];
			else if (NCLiquid == 2)
				dm_over_dt[1] = XGas[1][data->jFuel[1]]*pvDroplet[1] + XGas[1][data->jFuel[2]]*pvDroplet[2] - data->P_Pascal;
			else
				ErrorMessage("Only 2 liquid species can be considered up to now!");
				*/

			double sum=0.;
			for(int j=1;j<=NCLiquid;j++)
				sum += XGas[1][data->jFuel[j]];
			dm_over_dt[1] = 0.;
			for(int j=1;j<=NCLiquid;j++)
				dm_over_dt[1] += data->XDroplet[data->jFuel[j]]*pvDroplet[j];
			dm_over_dt[1] -= data->P_Pascal*sum;
		}
	}
	else if (data->iMode == UNSTEADY_BATCH)
	{
		if (NCLiquid == 1)
			dm_over_dt[1] = pvDroplet[1]/data->P_Pascal - OmegaGas[1][data->jFuel[1]]*MW[1]/mix->M[data->jFuel[1]];
		else if (NCLiquid == 2)
		{
		
			double sum=0.;
			for(int j=1;j<=NCLiquid;j++)
				sum += XGas[1][data->jFuel[j]];
			dm_over_dt[1] = 0.;
			for(int j=1;j<=NCLiquid;j++)
				dm_over_dt[1] += data->XDroplet[data->jFuel[j]]*pvDroplet[j];
			dm_over_dt[1] -= data->P_Pascal*sum;
		//	dm_over_dt[1] = data->XDroplet[data->jFuel[1]]*pvDroplet[1] + data->XDroplet[data->jFuel[2]]*pvDroplet[2] - data->P_Pascal*sum;
		
			// TODO (very provisional)
			if (data->enhancingFactor > 0.)
			{
				double sum=0.;
				for(int j=1;j<=NCLiquid;j++)
					sum += XGas[1][data->jFuel[j]];
				dm_over_dt[1] = 0.;
				double xglycerol = XGlycerol(t, data->enhancingFactor);
				dm_over_dt[1] += (1.-xglycerol)*pvDroplet[1];
				dm_over_dt[1] +=      xglycerol*pvDroplet[2];
				dm_over_dt[1] -= data->P_Pascal*sum;	
			}
		}
		else
			ErrorMessage("Only 2 liquid species can be considered up to now!");
	}

	// Internal points
	if (data->iVelocity == false)
	{
		for(int i=2;i<=N;i++)
			dm_over_dt[i] = m[i] - m[i-1];
	}
	else if (data->iVelocity == true)
	{
		for(int i=2;i<=N-1;i++)
		{
			double sum_dt = 0.;
			for(int k=1;k<=NCGas;k++)
				sum_dt += dOmegaGas_over_dt[i][k]/mix->M[k];
			sum_dt *= MW[i];

			double sum_dr = 0.;
			for(int k=1;k<=NCGas;k++)
				sum_dr += dOmegaGas_over_dr[i][k]/mix->M[k];
			sum_dr *= MW[i];

			dm_over_dt[i] = u[i]*(2./grid.x[i]-sum_dr-dT_over_dr[i]/T[i]) + du_over_dr[i] - sum_dt - dT_over_dt[i]/T[i];
		}
		dm_over_dt[N] = m[N] - m[N-1];
	}

}

void OpenSMOKE_DropletMicrogravity::GiveMeDropletMassEquation(const double t)	// MassEq
{
	// Interface
	dmass_over_dt[1] = -m[1];

	// Internal points
	for(int i=2;i<=N;i++)
		dmass_over_dt[i] = mass[i] - mass[i-1];
}

void OpenSMOKE_DropletMicrogravity::GiveMeDT(const double t)
{
	// a. Convective term [W/m]
	for(int i=1;i<=N;i++)
		T_convection[i] = -rho[i]*u[i]*Cp[i]*grid.A[i]*dT_over_dr[i];

	// b. Conduction term [W/m]
	for(int i=2;i<N;i++)
		T_conduction[i] = (  A_x_lambdae[i]*(T[i+1]-T[i])*grid.udxe[i] - A_x_lambdaw[i]*(T[i]-T[i-1])*grid.udxw[i] ) 
															*grid.udxc_over_2[i];;

	// c. Species diffusion fluxes [W/m]
	for(int i=2;i<N;i++)
	{
		double sumCpDiffusive = 0.;
		for(int k=1;k<=NCGas;k++)
			sumCpDiffusive += Cpk[i][k] * 0.50*(vStar_e[i][k]+vStar_w[i][k]);
		T_diffusion_fluxes[i] = -grid.A[i]*rho[i]*sumCpDiffusive*dT_over_dr_central[i];
	}

	// d. Reactions [W/m]
	for(int i=1;i<=N;i++)
		T_reaction[i] = grid.A[i]*QReaction[i];

	// e. Radiative heat transfer [W/m]
//	if (data->iRadiation == 1)
//		radiation.Calculate(grid.x, T, Kp);
	
	// Equations (Teq)
	if (data->iLiquidPhase == LIQUID_DROPLET_PERFECTLY_STIRRED)
	{
		if (data->iMode == EIGENVALUE)
			dT_over_dt[1] = InterfaceLiquidVaporizationHeat(t)*m[1] - lambda[1] * dT_over_dr[1] * grid.A[1];
		else if (data->iMode == UNSTEADY_BATCH)
			dT_over_dt[1] = (lambda[1] * dT_over_dr[1] * grid.A[1] - InterfaceLiquidVaporizationHeat(t)*m[1]) / 
											(InterfaceLiquidSpecificHeat()*mass[1]);

		if (data->iRadiation == 1)
			dT_over_dt[1] += - (radiation.qSurface* grid.A[1]) / (InterfaceLiquidSpecificHeat()*mass[1]);
	}
	else
	{
	//	cout << droplet->conductionFluxInterface() << endl;
		dT_over_dt[1] = lambda[1] * dT_over_dr[1] * grid.A[1] - InterfaceLiquidVaporizationHeat(t)*m[1] - droplet->conductionFluxInterface();
		if (data->iRadiation == 1)
			dT_over_dt[1] += - radiation.qSurface* grid.A[1];
	}
	
	for(int i=2;i<N;i++)
		dT_over_dt[i] = (T_conduction[i] + T_convection[i] + T_diffusion_fluxes[i] + T_reaction[i])  /
												(rho[i]*Cp[i]*grid.A[i]);

	if (data->iRadiation == 1)
		for(int i=2;i<N;i++)
			dT_over_dt[i] -= radiation.divq[i] / (rho[i]*Cp[i]);

	
	if (data->boundaryConditions == BCS_NEUMANN)
			dT_over_dt[N] = T[N] - T[N-1];

	else if (data->boundaryConditions == BCS_DIRICHLET)
			dT_over_dt[N] = T[N] - data->TEnvironment;
}

void OpenSMOKE_DropletMicrogravity::GiveMeDOmegaGas(double t)
{
	// a. Convective term [kg/s/m]
	for(int i=2;i<N;i++)
		for(int k=1;k<=NCGas;k++)
			OmegaGas_convection[i][k] = -rho[i]*u[i]*grid.A[i]*dOmegaGas_over_dr[i][k];

	// b. Diffusion terms
	for(int i=2;i<N;i++)
		for(int k=1;k<=NCGas;k++)
			OmegaGas_diffusion[i][k] = 	- (A_x_rhoe[i]*vStar_e[i][k] - A_x_rhow[i]*vStar_w[i][k] )*grid.udxc_over_2[i];

	// c. Reaction terms
	for(int i=1;i<=N;i++)
		for(int k=1;k<=NCGas;k++)
			OmegaGas_reaction[i][k] = RGas[i][k]*grid.A[i];


	// Equations (omegaEq)
	for(int k=1;k<=NCGas;k++)
		dOmegaGas_over_dt[1][k] = OmegaGas[1][k] - data->OmegaDroplet[k];

	// TODO (very provisional)
	if (data->enhancingFactor > 0.)
	{
		double yglycerol = YGlycerol(t, data->enhancingFactor);
		for(int k=1;k<=NCGas;k++)
			dOmegaGas_over_dt[1][k] = OmegaGas[1][k] - data->OmegaDroplet[k];
		dOmegaGas_over_dt[1][data->jFuel[1]] = OmegaGas[1][data->jFuel[1]] - (1.-yglycerol);
		dOmegaGas_over_dt[1][data->jFuel[2]] = OmegaGas[1][data->jFuel[2]] - yglycerol;
	}
	
	for(int k=1;k<=NCGas;k++)
		dOmegaGas_over_dt[1][k] =   (OmegaGas[1][k] - data->OmegaDroplet[k])*m[1]
     					          + rho[1]*grid.A[1]*vStar_w[1][k] ;

	for(int i=2;i<N;i++)
		for(int k=1;k<=NCGas;k++)
			dOmegaGas_over_dt[i][k] = (OmegaGas_diffusion[i][k]+OmegaGas_convection[i][k] + OmegaGas_reaction[i][k]) /
														(rho[i]*grid.A[i]);


	if (data->boundaryConditions == BCS_NEUMANN)
		for(int k=1;k<=NCGas;k++)
			dOmegaGas_over_dt[N][k] = OmegaGas[N][k] - OmegaGas[N-1][k];

	else if (data->boundaryConditions == BCS_DIRICHLET)
		for(int k=1;k<=NCGas;k++)
			dOmegaGas_over_dt[N][k] = OmegaGas[N][k] - data->OmegaEnvironment[k];
}

void OpenSMOKE_DropletMicrogravity::setMinimumAndMaximumValues()
{
	int k;
	double ZERO  =  0.;
	double ONE   =  1.;
	double TMIN  =  200.;
	double TMAX  =  6000.;
	double UMAX  =  1.e6;
	double MMAX  =  1.e0;
	double MASSMIN  =  1.e-14;
	double MASSMAX  =  1.e0;

	if (data->iMode == EIGENVALUE)
	{
		k=1;
		for(int i=1;i<=N;i++)
		{
			xMin[k++] =  TMIN;			// Temperature
			xMin[k++] =  MMIN;			// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				xMin[k++] = ZERO;		// Mass fractions
		}

		k=1;
		for(int i=1;i<=N;i++)
		{
			xMax[k++] = TMAX;			// Temperature
			xMax[k++] = MMAX;			// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				xMax[k++] = ONE;		// Mass fractions
		}
	}

	else if (data->iMode == UNSTEADY_BATCH)
	{
		k=1;
		for(int i=1;i<=N;i++)
		{
			xMin[k++] =  MASSMIN;		// Droplet mass
			xMin[k++] =  TMIN;			// Temperature
			xMin[k++] =  MMIN;			// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				xMin[k++] = ZERO;		// Mass fractions
		}

		k=1;
		for(int i=1;i<=N;i++)
		{
			xMax[k++] = MASSMAX;		// Droplet mass
			xMax[k++] = TMAX;			// Temperature
			xMax[k++] = MMAX;			// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				xMax[k++] = ONE;		// Mass fractions
		}
	}
}

void OpenSMOKE_DropletMicrogravity::setDifferentialAndAlgebraic()
{
	const int DIFFERENTIAL	= 1;
	const int ALGEBRAIC		= 0;

	if (data->iMode == EIGENVALUE)
	{
		int k=1;
		inDerAlg[k++] = ALGEBRAIC;				// Temperature
		inDerAlg[k++] = ALGEBRAIC;				// Mass flow rate
		for(int j=1;j<=NCGas;j++)
			inDerAlg[k++] = ALGEBRAIC;			// Mass fractions

		for(int i=2;i<N;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;		// Temperature
			inDerAlg[k++] = ALGEBRAIC;			// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				inDerAlg[k++] = DIFFERENTIAL;	// Mass fractions
		}

		inDerAlg[k++] = ALGEBRAIC;				// Temperature
		inDerAlg[k++] = ALGEBRAIC;				// Mass flow rate
		for(int j=1;j<=NCGas;j++)
			inDerAlg[k++] = ALGEBRAIC;			// Mass fractions
	}
	else if (data->iMode == UNSTEADY_BATCH)
	{
		int k=1;
		inDerAlg[k++] = DIFFERENTIAL;			// Mass droplet
		if (data->iLiquidPhase == LIQUID_DROPLET_PERFECTLY_STIRRED)
			inDerAlg[k++] = DIFFERENTIAL;		// Temperature
		else
			inDerAlg[k++] = ALGEBRAIC;		// Temperature
		inDerAlg[k++] = ALGEBRAIC;				// Mass flow rate
		for(int j=1;j<=NCGas;j++)
			inDerAlg[k++] = ALGEBRAIC;			// Mass fractions

		for(int i=2;i<N;i++)
		{
			inDerAlg[k++] = ALGEBRAIC;			// Mass droplet
			inDerAlg[k++] = DIFFERENTIAL;		// Temperature
			inDerAlg[k++] = ALGEBRAIC;			// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				inDerAlg[k++] = DIFFERENTIAL;	// Mass fractions
		}

		inDerAlg[k++] = ALGEBRAIC;				// Mass droplet
		inDerAlg[k++] = ALGEBRAIC;				// Temperature
		inDerAlg[k++] = ALGEBRAIC;				// Mass flow rate
		for(int j=1;j<=NCGas;j++)
			inDerAlg[k++] = ALGEBRAIC;			// Mass fractions
	}
}

void OpenSMOKE_DropletMicrogravity::setFirstGuess()
{
	if (data->iMode == EIGENVALUE)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			xFirstGuess[k++] = T[i];					// Temperature
			xFirstGuess[k++] = m[i];					// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				xFirstGuess[k++] = OmegaGas[i][j];		// Mass fractions
		}
	}
	else if (data->iMode == UNSTEADY_BATCH)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			xFirstGuess[k++] = mass[i];					// Droplet mass
			xFirstGuess[k++] = T[i];					// Temperature
			xFirstGuess[k++] = m[i];					// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				xFirstGuess[k++] = OmegaGas[i][j];		// Mass fractions
		}
	}
}

void OpenSMOKE_DropletMicrogravity::recoverUnknowns(BzzVector &x)
{	
	if (data->iMode == EIGENVALUE)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			T[i] = x[k++];					// Temperature
			m[i] = x[k++];					// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				OmegaGas[i][j] = x[k++];	// Mass fractions
		}
	}
	else if (data->iMode == UNSTEADY_BATCH)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			mass[i] = x[k++];				// Droplet mass
			T[i] = x[k++];					// Temperature
			m[i] = x[k++];					// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				OmegaGas[i][j] = x[k++];	// Mass fractions
		}
	}
}

void OpenSMOKE_DropletMicrogravity::recoverResiduals(BzzVector &f)
{
	if (data->iMode == EIGENVALUE)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			f[k++] = dT_over_dt[i];					// Temperature
			f[k++] = dm_over_dt[i];					// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				f[k++] = dOmegaGas_over_dt[i][j];	// Mass fractions
		}
	}
	else if (data->iMode == UNSTEADY_BATCH)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			f[k++] = dmass_over_dt[i];				// Droplet mass
			f[k++] = dT_over_dt[i];					// Temperature
			f[k++] = dm_over_dt[i];					// Mass flow rate
			for(int j=1;j<=NCGas;j++)
				f[k++] = dOmegaGas_over_dt[i][j];	// Mass fractions
		}
	}
}


void OpenSMOKE_DropletMicrogravity::Velocities()
{
	for(int i=1;i<=N;i++)
		u[i] = m[i]/(rho[i]*grid.A[i]);	
}

void OpenSMOKE_DropletMicrogravity::MassFlow()
{
	for(int i=1;i<=N;i++)
		m[i] = u[i]*rho[i]*grid.A[i];	
}

void OpenSMOKE_DropletMicrogravity::setInitialConditions()
{	
		// Resetting initial conditions
		{
			T = data->TEnvironment;
			m = 0.;
			OmegaGas = 0.;
		}

		// From zero
		if (data->iBackupFromBinaryFile == false)
		{
			// Interface initial conditions
			{
				// Temperature
				T[1] = data->TDroplet0;
				if (data->dropletInterfaceTemperature == INTERFACE_USER_DEFINED)		
					T[1] = data->interfaceTemperature;

				// Update liquid properties
				UpdateLiquidProperties();

				// Composition
				{
					BzzVector xInterface(NCGas);
					BzzVector omegaInterface(NCGas);
					double MWInterface;
					for (int j=1;j<=NCLiquid;j++)
						xInterface[data->jFuel[j]] = pvDroplet[j] / data->P_Pascal;
					double sum = xInterface.GetSumElements();
			
					BzzVector xEnvironment(NCGas);
					double MWEnvironment;
					mix->GetMWAndMoleFractionsFromMassFractions(MWEnvironment, xEnvironment, data->OmegaEnvironment);
					xInterface[data->jOxidizer] = xEnvironment[data->jOxidizer]*(1.-sum);
					xInterface[data->jInert]    = xEnvironment[data->jInert]*(1.-sum);

					// Mass fractions at the interface
					mix->GetMWAndMassFractionsFromMoleFractions(MWInterface, omegaInterface, xInterface);
					if ((omegaInterface.GetSumElements() > 1.+1.e-5) || (omegaInterface.GetSumElements() < 1.-1.e-5))
						ErrorMessage("The sum of mass fractions at the interface is not equal to 1...");
					else
					{
						double sum = omegaInterface.GetSumElements();
						omegaInterface /= sum;
					}
					OmegaGas.SetRow(1, omegaInterface);
				}
			}

			// Air
			for(int i=2;i<=N;i++)
			{
				T[i] = data->TEnvironment;
				m[i] = m[1];

				OmegaGas[i][data->jOxidizer] = data->OmegaEnvironment[data->jOxidizer];
				OmegaGas[i][data->jInert] = data->OmegaEnvironment[data->jInert];
			}

			// Spark
			if (data->iSpark == true)
				SparkProfileComplete(grid.x, T, data->SparkRatio, T[1], data->Tpeak, data->TEnvironment);

			// Mass flow rate
			UpdateProperties();
			double mdot = lambda[1]*(T[2]-T[1])*grid.udxe[1]*grid.A[1] / InterfaceLiquidVaporizationHeat(0);
			m = mdot;

			// Droplet mass
			mass = DropletMass();
		}

		// From backup file
		else if (data->iBackupFromBinaryFile == true)
		{
			ReadFromBackupFile(*mix, data->nameBackupFile,grid.x, T, m, OmegaGas);

			if (data->iSpark == true)
				SparkProfilePartial(grid.x, T, data->SparkRatio, data->Tpeak);

			UpdateLiquidProperties();
			mass = DropletMass();
		}

		// Summary on video
		cout << endl;
		cout << "Interface **" << endl;
		cout << " * Temperature [K]:       " << T[1] << endl;
		cout << " * Mass flow rate [kg/s]: " << m[1] << endl;
		cout << " * Mass fractions" << endl;
		for (int j=1;j<=NCGas;j++)
			if (OmegaGas[1][j] != 0.)	
				cout << "     Omega " << mix->names[j] << ": " << OmegaGas[1][j] << endl;
		cout << endl;
	
		// Properties
		UpdateProperties();
		UpdateLiquidProperties();
		Velocities();
		
		// Update radiative heat transfer [W/m3]
		if (data->iRadiation == 1)
		{
			radiation.Initialize(N);
			radiation.SetOptions(data->radiationOptions);
			radiation.Calculate(grid.x, T, Kp);
		}

		// Liquid droplet
		if (data->iLiquidPhase != LIQUID_DROPLET_PERFECTLY_STIRRED)
		{
			droplet->SetInitialConditions();
		}
}

void OpenSMOKE_DropletMicrogravity::InitialEstimations()
{
	double TFlame		= 2200.;		// Flame temperature [K]
	double Hc			= 44926.e3;		// Combustion heat [J/kg]
	double nu			= 15.08;		// Stoichiometry

	double Tinterface=data->TDroplet0;
	if (data->dropletInterfaceTemperature == INTERFACE_USER_DEFINED)		Tinterface = data->interfaceTemperature;
	double TMean = 0.50*(Tinterface+TFlame);

	double mdot;
	double Boq;
	double Kappa;
	
	// Mole fractions
	double MWFuel, MWOxidizer;
	BzzVector xOxidizer(NCGas), xFuel(NCGas);
	mix->GetMWAndMoleFractionsFromMassFractions(MWFuel, xFuel, data->OmegaDroplet);
	mix->GetMWAndMassFractionsFromMoleFractions(MWOxidizer, xOxidizer, data->OmegaEnvironment);
	
	// Thermal conductivity
	mix->SpeciesConductivityFromFitting(TMean);
	double lambdaFuel		= mix->MixConductivity_FromMolarFractions(xFuel);	// [W/m/K]
	double lambdaOxidizer	= mix->MixConductivity_FromMolarFractions(xOxidizer);	// [W/m/K]
	double lambdaG			= 0.40*lambdaFuel+0.60*lambdaOxidizer;

	// Specific heat
	mix->SpeciesCp(TMean);
	double CpFuel = mix->MixCp_FromMassFractions(data->OmegaDroplet);			// [J/kg/K]	
	double CpOxidizer = mix->MixCp_FromMassFractions(data->OmegaEnvironment);	// [J/kg/K]	

	// Evaporation rate
	Boq = (Hc/nu + CpFuel*(data->TEnvironment-Tinterface)) / InterfaceLiquidVaporizationHeat(0);
	mdot = 4.*Constants::pi*lambdaG*data->diameterDroplet0/2./CpFuel*log(Boq+1.);
	Kappa = 8.*lambdaG/InterfaceLiquidDensity()/CpFuel*log(Boq+1.);

	// Flame temperature and radius
	double Tf = InterfaceLiquidVaporizationHeat(0)/(CpOxidizer*(1.+nu))*(nu*Boq-1.) + Tinterface;
	double csif = log(1.+Boq)/log((nu+1.)/nu);
	double rf = data->diameterDroplet0/2.*csif;
	double omegaf = (Boq-1./nu)/(Boq+1.);

	// Droplet lifetime
	double tDroplet = BzzPow2(data->diameterDroplet0)/Kappa;

	double ZT = CpOxidizer/(4.*Constants::pi*lambdaG);
	double ZF = ZT; 
	double c1 = ZT*mdot/rf;
	double c2 = ZT*mdot/(data->diameterDroplet0/2.);

	ChangeDimensions(N, &TInitialEstimation);
	ChangeDimensions(N, NCGas, &OmegaInitialEstimation);
	
	for(int i=1;i<=N;i++)
		if (grid.x[i] <= rf)
		{
			TInitialEstimation[i] = ( (Tinterface-TFlame)*exp(-ZT*mdot/grid.x[i])+TFlame*exp(-c2)-Tinterface*exp(-c1) ) /
															( exp(-c2)-exp(-c1));
			OmegaInitialEstimation[i][data->jFuel[1]] = (1.-((1.-omegaf)*exp(-ZF*mdot/grid.x[i]))/exp(-c2))*.90;	
			OmegaInitialEstimation[i][data->jOxidizer] = (1.-((1.-omegaf)*exp(-ZF*mdot/grid.x[i]))/exp(-c2))*.10;	
			OmegaInitialEstimation[i][data->jInert] = 1.-OmegaInitialEstimation[i][data->jFuel[1]]-OmegaInitialEstimation[i][data->jOxidizer];	
		}
		else
		{
			TInitialEstimation[i] = ( (Tf-data->TEnvironment)*exp(-ZT*mdot/grid.x[i])+data->TEnvironment*exp(-c1)-Tf ) /
																	( exp(-c1)-1.);
			OmegaInitialEstimation[i][data->jOxidizer] = nu*(exp(-ZF*mdot/grid.x[i])/exp(-c1)-1.);
			OmegaInitialEstimation[i][data->jInert] = 1.-OmegaInitialEstimation[i][data->jOxidizer];
		}

	// Summary on video
	cout << "Effective thermal conductivity [W/m/K]: " << lambdaG << endl;
	cout << "Specific heat (fuel) [J/kg/K]:          " << CpFuel  << endl;
	cout << "Specific heat (oxid.) [J/kg/K]:         " << CpOxidizer  << endl;
	cout << "Transfer number B [-]:                  " << Boq  << endl;
	cout << "Evaporation mass flow rate [kg/s]:      " << mdot  << endl;
	cout << "Evaporation constant [mm2/s]:           " << Kappa*1e6  << endl;
	cout << "Flame temperature [K]:                  " << Tf  << endl;
	cout << "Flame radius [-]:                       " << csif  << endl;
	cout << "Fuel mass fraction at interface[-]:     " << omegaf  << endl;
	cout << "Ox mass fraction at oo[-]:              " << OmegaInitialEstimation[N][data->jOxidizer]  << endl;
	cout << "Droplet lifetime [ms]:                  " << tDroplet*1e3  << endl;
	cout << "Vaporization rate [kg/s]:               "	<< mdot << endl;
	cout << "Vaporization rate [mm2/s]:              "	<< 4.*mdot/InterfaceLiquidDensity()/data->diameterDroplet0*1e6 << endl;

}

void OpenSMOKE_DropletMicrogravity::MoleFractionsAndMolecularWeights()
{
	int i, j;

	for(i=1;i<=N;i++)
	{
		uMW[i] = 0.;
		for(j=1;j<=NCGas;j++)
			uMW[i] += OmegaGas[i][j]*mix->uM[j];
		MW[i] = 1./uMW[i];
		for(j=1;j<=NCGas;j++)
			XGas[i][j] = OmegaGas[i][j] * MW[i] * mix->uM[j];
	}
}

void OpenSMOKE_DropletMicrogravity::DiffusionVelocities()
{
	vc_e = 0.;
	
	for(int i=1;i<=N-1;i++)
		for(int j=1;j<=NCGas;j++)
		{
			coeff_e[i][j]	= - 0.50 * mix->M[j]*(Dm[i][j]/MW[i]+Dm[i+1][j]/MW[i+1]);
			vStar_e[i][j]	 = coeff_e[i][j] * (XGas[i+1][j]-XGas[i][j])*grid.udxe[i];
			vc_e[i]			-= vStar_e[i][j];
		}

	// Correzione delle velocita diffusive
	for(int i=1;i<=N-1;i++)
		for(int j=1;j<=NCGas;j++)
		{
			vStar_e[i][j]   += vc_e[i] * 0.50*(OmegaGas[i][j]+OmegaGas[i+1][j]);
			vStar_w[i+1][j]  = vStar_e[i][j];
		}

	// Fake mass and mole fractions
	BzzMatrix W_c(N, NCGas);
	BzzMatrix X_c(N, NCGas);
	BzzVector MW_c(N);
	{
		BzzVector aux_c(N);
		for(int i=1;i<=N;i++)
		{
			aux_c = OmegaGas.GetRow(i);
			double sum = aux_c.GetSumElements();
			for(int j=1;j<=NCGas;j++)	
				W_c[i][j] = OmegaGas[i][j]/sum;
		}
		mix->GetMWAndMoleFractionsFromMassFractions(MW_c, X_c, W_c);
	}

	double vc_w = 0.;
	for(int j=1;j<=NCGas;j++)
	{
			vStar_w[1][j]	 = - mix->M[j]*Dm[1][j]/MW[1] * (XGas[2][j]-XGas[1][j])/(grid.x[2]-grid.x[1]);
			vc_w			-= vStar_w[1][j];
	}

	for(int j=1;j<=NCGas;j++)
		vStar_w[1][j]   += vc_w * W_c[1][j];
}

void OpenSMOKE_DropletMicrogravity::UpdateProperties()
{
	double cGasTot;

	MoleFractionsAndMolecularWeights();

	for(int i=1;i<=N;i++)
	{
		// Mass and mole fractions
		OmegaGas.GetRow(i,&OmegaGasVector);
		XGas.GetRow(i,&XGasVector);

		// a. Concentration and density
		cGasTot = data->P_Pascal  / (Constants::R_J_kmol*T[i]);		// [kmol/m3]
		rho[i] = cGasTot * MW[i];									// [kg/m3]
		CGasVector = cGasTot*XGasVector;							// [kmol/m3]

		// b. Specific heat
		mix->SpeciesCp(T[i]);
		Cp[i] = mix->MixCp_FromMassFractions(OmegaGasVector);	// [J/kg/K]	
		Cpk.SetRow(i, mix->Cp);									// [J/kg/K]

		// c. Thermal conductivity
		mix->SpeciesConductivityFromFitting(T[i]);
		lambda[i] = mix->MixConductivity_FromMolarFractions(XGasVector);	// [W/m/K]

		// e. Diffusion coefficients
		mix->SpeciesDiffusivityFromFitting(T[i], data->P_bar);
		mix->MixDiffusionCoefficient_FromMolarFractionsAndPM(DmixVector,XGasVector);
		Dm.SetRow(i, DmixVector);											// [m2/s]
		
		// f. Radiative heat transfer
		if (data->iRadiation == 1)
			AbsorptionCoefficients();

		// f. Soret coefficients
		if (data->iSoretEffect == true)
		{
			mix->SpeciesThermalDiffusionRatiosFromFitting(T[i]);
			mix->MixThermalDiffusionRatios(TetaMixVector, XGasVector);
			Teta.SetRow(i, TetaMixVector);
		}

		// g. Chemical reactions
		if (data->iReactions == true)
		{
			mix->ComputeKineticParameters( T[i], log(T[i]), 1./T[i], data->P_Pascal);
			mix->ComputeFromConcentrations( T[i], CGasVector, cGasTot, &RGasVector);	// [kmol/m3/s]
			ElementByElementProduct(RGasVector, mix->M, &RGasVector);					// [kg/m3/s]
			RGas.SetRow(i, RGasVector);													// [kg/m3/s]
			QReaction[i] = - mix->ComputeQReaction(T[i]);								// [J/m3/s]
		}

		// Global Kinetics
		if (data->iGlobalKinetics == true)
		{
			global->GiveMeFormationRates(T[i],CGasVector,RGasVector);
			global->GiveMeReactionHeat(T[i], RGasVector, QReaction[i]);
			RGas.SetRow(i, RGasVector);			// [kg/m3/s]
		}
	}

	// Density at centers
	A_x_rhow[1] = grid.A[1]*rho[1];
	for(int i=1;i<N;i++)
	{
		A_x_rhoe[i] = 0.50*(grid.A[i]*rho[i]+grid.A[i+1]*rho[i+1]);
		A_x_rhow[i+1] = A_x_rhoe[i];
	}
	A_x_rhoe[N] = grid.A[N]*rho[N];


	// Thermal conductivity at centers
	A_x_lambdaw[1] = grid.A[1]*lambda[1];
	for(int i=1;i<N;i++)
	{
		A_x_lambdae[i] = 0.50 * (grid.A[i]*lambda[i] + grid.A[i+1]*lambda[i+1]);
		A_x_lambdaw[i+1] = A_x_lambdae[i];
	}
	A_x_lambdae[N] = grid.A[N]*lambda[N];
}

void OpenSMOKE_DropletMicrogravity::UpdateLiquidProperties()
{
	for(int j=1;j<=NCLiquid;j++)
	{
		rhoDroplet[j]		= data->liquid_species[j].rho(T[1]);
		dHVaporization[j]   = data->liquid_species[j].Hv(T[1]);
		pvDroplet[j]		= data->liquid_species[j].Pv(T[1]);
		cpDroplet[j]		= data->liquid_species[j].Cp(T[1]);
	}
}


void OpenSMOKE_DropletMicrogravity::Run()
{
	if (data->iMode == EIGENVALUE)
	{
		setInitialConditions();
		Summary();	
		PrepareSystem();
		daeSystemEigenValue.assignDroplet(this);
		cout << " a. Initial solution..." << endl;
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &daeSystemEigenValue, dimensionBlock);	
		cout << " b. Time integration..." << endl;
		DAESystemSolution(&o, data->tEnd);  
		PrintFinalSolution();
		PrintFinalSolutionProperties();
		PrintFinalSummary();
	}
	else if (data->iMode == UNSTEADY_BATCH)
	{
		setInitialConditions();
		Summary();	
		PrepareSystem();
		daeSystemUnsteadyBatch.assignDroplet(this);
		cout << " a. Initial solution..." << endl;
		BzzDaeSparseObject o(xFirstGuess, 0., inDerAlg, &daeSystemUnsteadyBatch, dimensionBlock);	
		cout << " b. Time integration..." << endl;
		DAESystemSolution(&o, data->tEnd);  
		PrintFinalSolution();
		PrintFinalSolutionProperties();
		PrintFinalSummary();
	}
}

void OpenSMOKE_DropletMicrogravity::DAESystemSolution(BzzDaeSparseObject *o, double tEnd)
{	
	o->StepPrint(DAE_Print);
	o->SetMinimumConstraints(xMin);
	o->SetMaximumConstraints(xMax);

	// Default values: (A) 1e-10      (R) 100*MachEps()
	o->SetTollRel(data->relTolerances);
	o->SetTollAbs(data->absTolerances);

	cout << "Absolute tolerance: " << data->absTolerances << endl;
	cout << "Relative tolerance: " << data->relTolerances << endl;

	double timeStart = BzzGetCpuTime();
	bzzStop = 0;
	xFirstGuess = (*o)(tEnd, tEnd);

	cout << endl;
	cout << "Number of steps: "					<< o->GetNumStep() << endl;
	cout << "Number of function for Jacobian: " << o->GetNumFunctionForJacobian() << endl;
	cout << "Numerical Jacobians: "				<< o->GetNumNumericalJacobian() << endl;
	
	cout << "Time DAE solution: "				<< BzzGetCpuTime() - timeStart << " s" << endl << endl;
}

int OpenSMOKE_DropletMicrogravity::nonLinearSystemSolution(BzzNonLinearSystemSparseObject &o, int dimensionBlock)
{
	o.SetMinimumConstraints(xMin);
	o.SetMaximumConstraints(xMax);

	o.SetTolerance(data->absTolerances, data->relTolerances);
	double startTime = BzzGetCpuTime();
	bzzStop = 0;
	char control = o();
	double endTime = BzzGetCpuTime() - startTime;
	
	BzzVector fNLS(NE);
	o.GetSolution(&xFirstGuess, &fNLS);
	
	double maxfNLS = fNLS.MaxAbs();
	
	double sum = 0.;
	for (int i=1;i<=(N*dimensionBlock);i++)
		sum += fabs(fNLS[i]);
	sum=sum/(N*dimensionBlock);

	cout << "Mean residual:" << sum << endl;
	cout << "Max. residual:" << maxfNLS << endl;

	cout << "Time NLS solution: " <<  endTime << " s" << endl;

	cout << endl;

	return int(control);
}

void OpenSMOKE_DropletMicrogravity::DAESystemEigenValue(BzzVector &x, double t, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recovering
	// --------------------------------------------------------------------------------------
	recoverUnknowns(x);

	// -------------------------------------------------------------------------------------		
	// 2. Properties
	// -------------------------------------------------------------------------------------
	UpdateProperties();
	UpdateLiquidProperties();
	Velocities();
	DiffusionVelocities();
	GiveMeSpatialDerivatives();

	// --------------------------------------------------------------------------------------
	// 3. Equations
	// --------------------------------------------------------------------------------------
	GiveMeDT(t);
	GiveMeDOmegaGas(t);
	GiveMeContinuityEquation(t);
	
	// --------------------------------------------------------------------------------------
	// 4. Recovering
	// --------------------------------------------------------------------------------------
	recoverResiduals(f);	
}

void OpenSMOKE_DropletMicrogravity::DAESystemUnsteadyBatch(BzzVector &x, double t, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recovering
	// --------------------------------------------------------------------------------------
	recoverUnknowns(x);

	// -------------------------------------------------------------------------------------		
	// 2. Properties
	// -------------------------------------------------------------------------------------
	UpdateProperties();
	UpdateLiquidProperties();
	Velocities();
	DiffusionVelocities();
	GiveMeSpatialDerivatives();

	// --------------------------------------------------------------------------------------
	// 3. Equations
	// --------------------------------------------------------------------------------------
	GiveMeDropletMassEquation(t);
	GiveMeDT(t);
	GiveMeDOmegaGas(t);
	GiveMeContinuityEquation(t);
	
	// --------------------------------------------------------------------------------------
	// 4. Recovering
	// --------------------------------------------------------------------------------------
	recoverResiduals(f);	
}

double OpenSMOKE_DropletMicrogravity::SootPropensity(const int j)
{
	const int jC2H2 = mix->recognize_species("C2H2");
	const double MWsoot = 12000.;						// [kg/kmol]
	const double A = 10.36e6;							// [1/s]
	const double Ea = 172e6;							// [J/kmol]

	double cC2H2 = (rho[j]/MW[j])*XGas[j][jC2H2];		// [kmol/m3]
	double k = A*exp(-Ea/Constants::R_J_kmol/T[j]);		// [1/s]
	return k*cC2H2 * MWsoot;							// [kg/m3/s]
}

double OpenSMOKE_DropletMicrogravity::SpeciesIntegralMass(const int j, const int indexSpecies)
{
	if (indexSpecies == 0)
		return 1.;

	return rho[j] * OmegaGas[j][indexSpecies];							// [kg/m3/s]
}

double OpenSMOKE_DropletMicrogravity::SootPropensity()
{
	double sum = 0.;
	for(int j=2;j<=N;j++)
	{
		double y1 = SootPropensity(j-1);
		double y2 = SootPropensity(j);
		sum += 0.50*(y1+y2) * Constants::pi/6. * ( BzzPow3(2.*grid.x[j])-BzzPow3(2.*grid.x[j-1]) );
	}
	return sum;
}

double OpenSMOKE_DropletMicrogravity::MassIntegral(const string names)
{
	int indexSpecies = 0;
	if (names != "volume")
		indexSpecies = mix->recognize_species(names);

	double sum = 0.;
	for(int j=2;j<=N;j++)
	{
		double y1 = SpeciesIntegralMass(j-1, indexSpecies);
		double y2 = SpeciesIntegralMass(j, indexSpecies);
		sum += 0.50*(y1+y2) * Constants::pi/6. * ( BzzPow3(2.*grid.x[j])-BzzPow3(2.*grid.x[j-1]) );
	}
	return sum;
}

void OpenSMOKE_DropletMicrogravity::AbsorptionCoefficients()
{
	for(int i=1; i<=N; i++)
	{
		Kp[i] = 0.;

		if(data->jH2O > 0)	Kp[i] += AbsorptionCoefficient_H2O(T[i])*XGas[i][data->jH2O];		// [1/m/atm]
		if(data->jCO2 > 0)	Kp[i] += AbsorptionCoefficient_CO2(T[i])*XGas[i][data->jCO2];		// [1/m/atm]
		if(data->jCO  > 0)	Kp[i] += AbsorptionCoefficient_CO(T[i])*XGas[i][data->jCO];			// [1/m/atm]
		if(data->jCH4 > 0)	Kp[i] += AbsorptionCoefficient_CH4(T[i])*XGas[i][data->jCH4];		// [1/m/atm]

		Kp[i] *= data->P_atm;		// [1/m]
	}
}
