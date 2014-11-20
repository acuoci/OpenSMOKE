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
#include "droplet/OpenSMOKE_Droplet.h"
#include "engine/OpenSMOKE_ReactingGas.h"

OpenSMOKE_Droplet *ptDroplet;

void DAE_Print_Droplet(BzzVector &x, double t)
{
	ptDroplet->DAE_myPrint(x, t);
}

void OpenSMOKE_Droplet::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Droplet"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Droplet::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_DropletMicrogravity"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

void OpenSMOKE_Droplet::MemoryAllocation()
{
	ChangeDimensions(N, &T);
	ChangeDimensions(N, &dT_over_dt);
	ChangeDimensions(N, &dT_over_dr);

	ChangeDimensions(N, &lambda);
	ChangeDimensions(N, &cp);
	ChangeDimensions(N, &rho);

	ChangeDimensions(N, &A_x_rhoe);
	ChangeDimensions(N, &A_x_rhow);
	ChangeDimensions(N, &A_x_lambdae);
	ChangeDimensions(N, &A_x_lambdaw);

	ChangeDimensions(N, &T_conduction);

	ChangeDimensions(N, &MW);
	ChangeDimensions(N, &uMW);
	
	ChangeDimensions(N, NC, &Omega_diffusion);

	ChangeDimensions(N, NC, &Omega);
	ChangeDimensions(N, NC, &dOmega_over_dt);
	ChangeDimensions(N, NC, &dOmega_over_dr);

	ChangeDimensions(NC, &OmegaVector);
	ChangeDimensions(NC, &XVector);
//	ChangeDimensions(NC, &mGasInterface);

	ChangeDimensions(N, NC, &X);
	ChangeDimensions(N, NC, &Dm);

	ChangeDimensions(N, &vc_e);
	ChangeDimensions(N, NC, &vStar_e);
	ChangeDimensions(N, NC, &vStar_w);
	ChangeDimensions(N, NC, &coeff_e);
}

void OpenSMOKE_Droplet::PrepareSystem()
{
	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_NOMOMENTUM)			dimensionBlock	= 1 + NC;
	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_ONLYT)	dimensionBlock	= 1;
	
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

OpenSMOKE_Droplet::OpenSMOKE_Droplet()
{
	cout.setf(ios::scientific);

	ptDroplet = this;
	grid.SetSpherical();
}

void OpenSMOKE_Droplet::SetName(const std::string name)
{
	name_object = name;
}

void OpenSMOKE_Droplet::Assign(OpenSMOKE_ReactingGas *_mix)
{
	mix = _mix;
}

void OpenSMOKE_Droplet::Assign(OpenSMOKE_DropletMicrogravity_DataManager *_data)
{
	data = _data;
}

void OpenSMOKE_Droplet::Setup()
{
	NC					= data->nameFuels.size();
	N					= data->NDroplet;
	stretchingFactor	= data->dropletStretchingFactor;
	if (data->dropletGridMode == "fixedRatio")
		gridMode		= LIQUID_DROPLET_GRID_FIXED_RATIO;
	else if (data->dropletGridMode == "exponential")
		gridMode		= LIQUID_DROPLET_GRID_EXPONENTIAL;
	
	MemoryAllocation();

	if (N<10)	ErrorMessage("Minimum number of points for liquid phase is 10!");

	// Grid construction
	if (gridMode == LIQUID_DROPLET_GRID_FIXED_RATIO)
		grid.Construct(N, data->diameterDroplet0/2., stretchingFactor, 0.);
	else if (gridMode == LIQUID_DROPLET_GRID_EXPONENTIAL)
		grid.ConstructExponential(N, data->diameterDroplet0/2., stretchingFactor);

	openOutputFileAndControl(fUnsteady, data->nameOutputFolder + "/DropletUnsteady.out");
	fUnsteady.setf(ios::scientific);
}	

void OpenSMOKE_Droplet::UpdateFromGasPhase(const double tStart_, const double tEnd_)
{
	tStart = tStart_;
	tEnd = tEnd_;
}

void OpenSMOKE_Droplet::UpdateFromGasPhase(const double TGasInterface_)
{
	TGasInterface = TGasInterface_;
}

/*
void OpenSMOKE_Droplet::UpdateFromGasPhase(const double mGas_, const double conductionGas_, const double radiationFlux_)
{
	mGas = mGas_;
	conductionGas = conductionGas_;
	radiationFlux = radiationFlux_;

	//	for(int k=1;k<=NC;k++)
	//		mGasInterface[k] = mGas_ * omegaGasInterface_[data->jFuel[k]];
}
*/

void OpenSMOKE_Droplet::GiveMeSpatialDerivatives()
{
	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_NOMOMENTUM)
	{
		grid.FirstDerivative('C', T, dT_over_dr);
		grid.FirstDerivative('C', Omega, dOmega_over_dr);
	}
	else if (data->iLiquidPhase == LIQUID_DROPLET_PDE_ONLYT)
	{
		grid.FirstDerivative('C', T, dT_over_dr);
	}
}

void OpenSMOKE_Droplet::GiveMeDT(const double t)
{
	// a. Convective term [W/m]

	// b. Conduction term [W/m]
	for(int i=2;i<N;i++)
		T_conduction[i] = (  A_x_lambdae[i]*(T[i+1]-T[i])*grid.udxe[i] - A_x_lambdaw[i]*(T[i]-T[i-1])*grid.udxw[i] ) 
															*grid.udxc_over_2[i];

	// d. Reactions [W/m]
	
	// Equations (Teq)
	dT_over_dt[1] = T[1] - T[2];

	for(int i=2;i<N;i++)
		dT_over_dt[i] = (T_conduction[i] )  / (rho[i]*cp[i]*grid.A[i]);

	dT_over_dt[N] = T[N] - TGasInterface;
/*
	double VaporizationFlux = 0.;
	for(int j=1;j<=NC;j++)
		VaporizationFlux += (Omega[N][j]*mGas) *data->liquid_species[j].Hv(T[N]);	// [kg/s*J/kg] = [J/s]
	
	double conductionLiquid = lambda[N] * dT_over_dr[N] * grid.A[N];
	dT_over_dt[N] = conductionLiquid + VaporizationFlux + radiationFlux - conductionGas;
	*/
/*	cout << lambda[N] << " " << dT_over_dr[N] << grid.A[N] << endl;
	cout << mGas << " " << Omega[N][1] << " " << Omega[N][2] << " " << VaporizationFlux/mGas << endl;
	cout << conductionLiquid << " " << radiationFlux << " " << conductionGas << " " << VaporizationFlux << endl;
	
	getchar();*/
}

void OpenSMOKE_Droplet::GiveMeDOmega(double t)
{/*
	// a. Convective term [kg/s/m]

	// b. Diffusion terms
	for(int i=2;i<N;i++)
		for(int k=1;k<=NC;k++)
			Omega_diffusion[i][k] = 	- (A_x_rhoe[i]*vStar_e[i][k] - A_x_rhow[i]*vStar_w[i][k] )*grid.udxc_over_2[i];

	// c. Reaction terms

	// Equations (omegaEq)
	for(int k=1;k<=NC;k++)
		dOmega_over_dt[1][k] = Omega[1][k] - Omega[2][k];

	for(int i=2;i<=N-1;i++)
		for(int k=1;k<=NC;k++)
			dOmega_over_dt[i][k] = (Omega_diffusion[i][k]) / (rho[i]*grid.A[i]);

	for(int k=1;k<=NC;k++)
		dOmega_over_dt[N][k] = grid.A[N]*rho[N]*vStar_w[N][k] + mGasInterface[k]*0;
*/}

void OpenSMOKE_Droplet::setMinimumAndMaximumValues()
{
	double ZERO  =  0.;
	double ONE   =  1.;
	double TMIN  =  200.;
	double TMAX  =  1000.;

	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_NOMOMENTUM)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			xMin[k++] =  TMIN;			// Temperature
			for(int j=1;j<=NC;j++)
				xMin[k++] = ZERO;		// Mass fractions
		}

		k=1;
		for(int i=1;i<=N;i++)
		{
			xMax[k++] = TMAX;			// Temperature
			for(int j=1;j<=NC;j++)
				xMax[k++] = ONE;		// Mass fractions
		}
	}

	else if (data->iLiquidPhase == LIQUID_DROPLET_PDE_ONLYT)
	{
		int k=1;
		for(int i=1;i<=N;i++)
			xMin[k++] =  TMIN;			// Temperature

		k=1;
		for(int i=1;i<=N;i++)
			xMax[k++] = TMAX;			// Temperature
	}
}

void OpenSMOKE_Droplet::setDifferentialAndAlgebraic()
{
	const int DIFFERENTIAL	= 1;
	const int ALGEBRAIC		= 0;

	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_NOMOMENTUM)
	{
		int k=1;
		inDerAlg[k++] = ALGEBRAIC;				// Temperature
		for(int j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;			// Mass fractions

		for(int i=2;i<N;i++)
		{
			inDerAlg[k++] = DIFFERENTIAL;		// Temperature
			for(int j=1;j<=NC;j++)
				inDerAlg[k++] = DIFFERENTIAL;	// Mass fractions
		}

		inDerAlg[k++] = ALGEBRAIC;				// Temperature
		for(int j=1;j<=NC;j++)
			inDerAlg[k++] = ALGEBRAIC;			// Mass fractions
	}

	else if (data->iLiquidPhase == LIQUID_DROPLET_PDE_ONLYT)
	{
		int k=1;
		inDerAlg[k++] = ALGEBRAIC;				// Temperature

		for(int i=2;i<N;i++)
			inDerAlg[k++] = DIFFERENTIAL;		// Temperature

		inDerAlg[k++] = ALGEBRAIC;				// Temperature
	}
}

void OpenSMOKE_Droplet::setFirstGuess()
{
	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_NOMOMENTUM)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			xFirstGuess[k++] = T[i];				// Temperature
			for(int j=1;j<=NC;j++)
				xFirstGuess[k++] = Omega[i][j];		// Mass fractions
		}
	}

	else if (data->iLiquidPhase == LIQUID_DROPLET_PDE_ONLYT)
	{
		int k=1;
		for(int i=1;i<=N;i++)
			xFirstGuess[k++] = T[i];				// Temperature
	}
}

void OpenSMOKE_Droplet::recoverUnknowns(BzzVector &x)
{	
	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_NOMOMENTUM)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			T[i] = x[k++];					// Temperature
			for(int j=1;j<=NC;j++)
				Omega[i][j] = x[k++];	// Mass fractions
		}
	}

	else if (data->iLiquidPhase == LIQUID_DROPLET_PDE_ONLYT)
	{
		int k=1;
		for(int i=1;i<=N;i++)
			T[i] = x[k++];					// Temperature
	}
}

void OpenSMOKE_Droplet::recoverResiduals(BzzVector &f)
{
	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_NOMOMENTUM)
	{
		int k=1;
		for(int i=1;i<=N;i++)
		{
			f[k++] = dT_over_dt[i];				// Temperature
			for(int j=1;j<=NC;j++)
				f[k++] = dOmega_over_dt[i][j];	// Mass fractions
		}
	}

	else if (data->iLiquidPhase == LIQUID_DROPLET_PDE_ONLYT)
	{
		int k=1;
		for(int i=1;i<=N;i++)
			f[k++] = dT_over_dt[i];				// Temperature
	}
}

void OpenSMOKE_Droplet::SetInitialConditions()
{	
	// Resetting initial conditions
	{
		T = data->TDroplet0;
		Omega = 0.;
	}

	// From zero
	if (data->iBackupFromBinaryFile == false)
	{
		T = data->TDroplet0;
		for(int i=1;i<=N;i++)
			for(int j=1;j<=NC;j++)
				Omega[i][j] = data->OmegaDroplet[data->jFuel[j]];

		UpdateProperties();
	}

	// From backup file
	else if (data->iBackupFromBinaryFile == true)
	{
		ErrorMessage("Backup not yet implemented...");
	//	ReadFromBackupFile(*mix, data->nameBackupFile, grid.x, T, Omega);
	//	UpdateProperties();
	}
	
	// Properties
	UpdateProperties();
}

void OpenSMOKE_Droplet::MoleFractionsAndMolecularWeights()
{
	int i, j;

	for(i=1;i<=N;i++)
	{
		uMW[i] = 0.;
		for(j=1;j<=NC;j++)
			uMW[i] += Omega[i][j]*mix->uM[data->jFuel[j]];
		MW[i] = 1./uMW[i];
		for(j=1;j<=NC;j++)
			X[i][j] = Omega[i][j] * MW[i] * mix->uM[data->jFuel[j]];
	}
}

void OpenSMOKE_Droplet::DiffusionVelocities()
{
	vc_e = 0.;
	
	for(int i=1;i<=N-1;i++)
		for(int j=1;j<=NC;j++)
		{
			coeff_e[i][j]	= - 0.50 * mix->M[data->jFuel[j]]*(Dm[i][j]/MW[i]+Dm[i+1][j]/MW[i+1]);
			vStar_e[i][j]	 = coeff_e[i][j] * (X[i+1][j]-X[i][j])*grid.udxe[i];
			vc_e[i]			-= vStar_e[i][j];
		}

	// Correzione delle velocita diffusive
	for(int i=1;i<=N-1;i++)
		for(int j=1;j<=NC;j++)
		{
			vStar_e[i][j]   += vc_e[i] * 0.50*(Omega[i][j]+Omega[i+1][j]);
			vStar_w[i+1][j]  = vStar_e[i][j];
		}

	double vc_w = 0.;
	for(int j=1;j<=NC;j++)
	{
		vStar_w[N][j]	 = - mix->M[data->jFuel[j]]*Dm[N][j]/MW[N] * (X[N][j]-X[N-1][j])/(grid.x[N]-grid.x[N-1]);
		vc_w			-= vStar_w[N][j];
	}

	for(int j=1;j<=NC;j++)
		vStar_w[N][j]   += vc_w * Omega[N][j];
}

void OpenSMOKE_Droplet::UpdateProperties()
{
	double cGasTot;

	MoleFractionsAndMolecularWeights();

	for(int i=1;i<=N;i++)
	{
		// Mass and mole fractions
		Omega.GetRow(i,&OmegaVector);
		X.GetRow(i,&XVector);

		// a. Concentration and density
		rho[i] = GetDensity(T[i], OmegaVector);						// [kg/m3]

		// b. Specific heat
		cp[i] = GetSpecificHeat(T[i], OmegaVector);					// [J/kg/K]	

		// c. Thermal conductivity
		lambda[i] = GetThermalConductivity(T[i], OmegaVector);		// [W/m/K]

		// e. Diffusion coefficients
		Dm.SetRow(i, GetDiffusionCoefficients(T[i], OmegaVector));	// [m2/s]
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

void OpenSMOKE_Droplet::Run()
{
	if (data->iLiquidPhase == LIQUID_DROPLET_PDE_NOMOMENTUM)
	{
		PrepareSystem();
		daeSystemNoMomentum.assignDroplet(this);
		cout << " a. Initial solution..." << endl;
		BzzDaeSparseObject o(xFirstGuess, tStart, inDerAlg, &daeSystemNoMomentum, dimensionBlock);	
		cout << " b. Time integration..." << endl;
		DAESystemSolution(&o);
	}

	else if (data->iLiquidPhase == LIQUID_DROPLET_PDE_ONLYT)
	{
		PrepareSystem();
		daeSystemOnlyTemperature.assignDroplet(this);
		cout << " a. Initial solution..." << endl;
		BzzDaeSparseObject o(xFirstGuess, tStart, inDerAlg, &daeSystemOnlyTemperature, dimensionBlock);	
		cout << " b. Time integration..." << endl;
		DAESystemSolution(&o);
	}
}

void OpenSMOKE_Droplet::DAESystemSolution(BzzDaeSparseObject *o)
{	
	o->StepPrint(DAE_Print_Droplet);
	o->SetMinimumConstraints(xMin);
	o->SetMaximumConstraints(xMax);

	// Default values: (A) 1e-10      (R) 100*MachEps()
	o->SetTolRel(data->relTolerances);
	o->SetTolAbs(data->absTolerances);

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

void OpenSMOKE_Droplet::DAESystemNoMomentum(BzzVector &x, double t, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recovering
	// --------------------------------------------------------------------------------------
	recoverUnknowns(x);

	// -------------------------------------------------------------------------------------		
	// 2. Properties
	// -------------------------------------------------------------------------------------
	UpdateProperties();
	DiffusionVelocities();
	GiveMeSpatialDerivatives();

	// --------------------------------------------------------------------------------------
	// 3. Equations
	// --------------------------------------------------------------------------------------
	GiveMeDT(t);
	GiveMeDOmega(t);
	
	// --------------------------------------------------------------------------------------
	// 4. Recovering
	// --------------------------------------------------------------------------------------
	recoverResiduals(f);	
}

void OpenSMOKE_Droplet::DAESystemOnlyTemperature(BzzVector &x, double t, BzzVector &f)
{
	// --------------------------------------------------------------------------------------
	// 1. Recovering
	// --------------------------------------------------------------------------------------
	recoverUnknowns(x);

	// -------------------------------------------------------------------------------------		
	// 2. Properties
	// -------------------------------------------------------------------------------------
	UpdateProperties();
	GiveMeSpatialDerivatives();

	// --------------------------------------------------------------------------------------
	// 3. Equations
	// --------------------------------------------------------------------------------------
	GiveMeDT(t);
	
	// --------------------------------------------------------------------------------------
	// 4. Recovering
	// --------------------------------------------------------------------------------------
	recoverResiduals(f);	
}

double OpenSMOKE_Droplet::GetDensity(const double t, const BzzVector &omega)
{
	double sum=0.;
	for (int j=1;j<=NC;j++)
		sum += data->liquid_species[j].rho(t) * omega[j];
	return sum;
}

double OpenSMOKE_Droplet::GetSpecificHeat(const double t, const BzzVector &omega)
{
	double sum=0.;
	for (int j=1;j<=NC;j++)
		sum += data->liquid_species[j].Cp(t) * omega[j];
	return sum;
}

double OpenSMOKE_Droplet::GetThermalConductivity(const double t, const BzzVector &omega)
{
	return 0.50e3;	// [W/m/K]
}

BzzVector OpenSMOKE_Droplet::GetDiffusionCoefficients(const double t, const BzzVector &omega)
{
	BzzVector D(NC);
	D = 1e-8;
	return D;
}

void OpenSMOKE_Droplet::DAE_myPrint(BzzVector &x, double t)
{
	cout << t << " " << T[N/2] << " " << T[N] << endl;

	double dummy = 0.;
	for (int i=1;i<=N;i++)
	{
		fUnsteady << setw(16) << left << t;						// 1 time [s]
		fUnsteady << setw(16) << left << grid.x[i]*1000.;		// 2 radial coordinate [mm]
		fUnsteady << setw(16) << left << grid.x[i]/grid.x[N];	// 3 radial coordinate [-]
		fUnsteady << setw(16) << left << T[i];					// 4 temperature [K]
		fUnsteady << setw(16) << left << dummy*1000.;			// 5 velocity [mm/s]
		fUnsteady << setw(16) << left << rho[i];				// 6 density [kg/m3]
		fUnsteady << setw(16) << left << cp[i];					// 7 specific heat [J/kg/K]
		fUnsteady << setw(16) << left << lambda[i];				// 8 thermal conductivity [W/m/K]
		fUnsteady << setw(16) << left << MW[i];					// 9 Molecular weight [kg/kmol]
		for(int j=1;j<=NC;j++)
			fUnsteady << setw(16) << left << Omega[i][j];
		fUnsteady << endl;
	}
	fUnsteady << endl;
}