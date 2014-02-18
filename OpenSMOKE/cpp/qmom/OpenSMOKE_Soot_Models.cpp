/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci   	                               *
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

#include "basic/OpenSMOKE_Constants.h"
#include "qmom/OpenSMOKE_Soot_Models.h"

double calculate_correction_coefficient(double T, double EoverR, double beta, double qan);

//////////////////////////////////////////////////////////////////////
// Non linear system to obtain the fractal dimension
//////////////////////////////////////////////////////////////////////
void MyNonLinearSystem_Fractal_Dimension::ObjectBzzPrint(void)
{
}
void MyNonLinearSystem_Fractal_Dimension::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptSoot->nls_fractal_dimension(x, f);
}

void MyNonLinearSystem_Fractal_Dimension::assignSoot(OpenSMOKE_Soot_Models *soot)
{
	ptSoot = soot;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void OpenSMOKE_Soot_Models::initialize(	const int _N, const int _iFractalDimension, const int _iNucleation,
								const int _iGrowth, const int _iOxidation, const int _iAggregation,
								const double _L0, const int _iTfluctuations)
{
	N	= _N;					// dimension of multi dirac distribution

	iNucleation			= _iNucleation;				// nucleation model
	iGrowth				= _iGrowth;					// growth model
	iOxidation			= _iOxidation;				// oxidation model
	iAggregation		= _iAggregation;			// aggregation model

	iFractalDimension	= _iFractalDimension;		// 1 = calculate fractal dimension 
													// 0 = spherical particles

	iTfluctuations		= _iTfluctuations;			// 1 = T fluctuations on
													// 0 = T fluctuations off

	L0Soot = _L0*1.e-9;			// primary soot particle size	[nm] --> [m]
	
	allocate_memory();			// allocation of memory

	define_constants();			// definition of physical constants and properties of soot

	user_defined_variables();	// user-defined varibles

	nls.assignSoot(this);		// setup of non linear system
}

void OpenSMOKE_Soot_Models::allocate_memory()
{
	ChangeDimensions(N,N, &Beta);
	ChangeDimensions(N, &Dp);
	ChangeDimensions(N, &Diff);
	ChangeDimensions(N, &c);
	ChangeDimensions(N, &g);

	ChangeDimensions(2, &xFirstGuess);
	ChangeDimensions(2, &fnls);
}

void OpenSMOKE_Soot_Models::define_constants()
{
	Df    = 2.40;				// fractal dimension (guess value)
	Dfmin = 1.80;				// minimum fractal dimension
	Dfmax = 3.00;				// maximum fractal dimension
	Df0	  = 2.40;				// mean fractal dimension

	rhoSoot =	1927.5946;			// soot density [kg/m3]
	PMSoot	=	12.010999679565430;	// soot molecular weight [kg/kmol]

	pi = acos(-1.0);			// [-]
	pi_over_6 = pi/6.;			// [-]
	_8_over_pi = 8./pi;			// [-]
	_2pi = 2.*pi;				// [-]
	_3pi = 3.*pi;				// [-]
	_pi_plus_8 = 3.*pi;			// [-]
	pi_squared = sqrt(pi);		// [-]
	_36pi_to_1_over_3 
		       = pow(36.*pi,1./3.); // [-]
}

void OpenSMOKE_Soot_Models::user_defined_variables()
{	
	sFractal = 1.;					// slope of fractal dimension variation [0.5-1.5]

	// Derived Variables
	L0Soot_3 = BzzPow3(L0Soot); // cube primary particle size [m3]
}

void OpenSMOKE_Soot_Models::update(BzzVector &_w, BzzVector &_L)
{
	int i;

	w	= _w;		// particle density [1/m3]
	L   = _L;		// particle sizes   [m]

	// Number density of soot particles [1/m3]
	m0 = w.GetSumElements();

	// Soot volume fraction [-]
	fv = 0.;
	for(i=1;i<=N;i++)
		fv += w[i]*(L[i]*L[i]*L[i]);
	fv *= pi_over_6;

	// Mean soot particle size - Linear [m]
	Lmean = 0.;
	for(i=1;i<=N;i++)
		Lmean += w[i]*L[i];
	Lmean /= m0;

	// Moment order 2
	m2 = 0.;
	for(i=1;i<=N;i++)
		m2 += w[i]*(L[i]*L[i]);

	// Mean soot particle size - Cubic [m]
	Lmean_Volume = (fv/pi_over_6)/m0;
	Lmean_Volume = pow(Lmean_Volume, (1./3.));

	// Specific surface [1/m]
	As = _36pi_to_1_over_3 * pow(m0, 1./3.) * pow(fv, 2./3.);
}

void OpenSMOKE_Soot_Models::assign_physical_properties(double _T, double _P, double _rho, 
											 double _mu, double _tr,
											 double _omega_O2, double _omega_C2H2,
											 double _omega_OH, double _STNV)
{
	T_gas = _T;							// temperature [K]
	P_gas = _P;							// pressure [Pa]
	rho_gas = _rho;						// gas density [kg/m3]
	mu_gas = _mu;						// viscosity [kg/ms]
	tr_gas = _tr;						// turbulence micro-scale time [s]

	cTot_gas = P_gas/(8314.*T_gas);		// total concentration [kmol/m3]
	PM_gas   = rho_gas/cTot_gas;		// molecular weight [kg/kmol]			

	omega_O2   = _omega_O2;				// oxygen mass fraction
	omega_C2H2 = _omega_C2H2;			// acetylene mass fraction
	omega_OH   = _omega_OH;				// OH mass fraction

	x_O2	 = omega_O2*PM_gas/PM_O2;		// mole fraction of oxygen
	x_C2H2	 = omega_C2H2*PM_gas/PM_C2H2;	// mole fraction of acetylene
	x_OH	 = omega_OH*PM_gas/PM_OH;		// mole fraction of OH

	c_O2	 = x_O2*cTot_gas;			// oxygen concentration [kmol/m3]
	c_C2H2	 = x_C2H2*cTot_gas;			// acetylene concentration [kmol/m3]
	c_OH	 = x_OH*cTot_gas;			// OH concentration [kmol/m3]

	p_O2	 = x_O2*P_gas/101325.;		// oxygen pressure [atm]
	p_C2H2	 = x_C2H2*P_gas/101325.;	// acetylene pressure [atm]
	p_OH	 = x_OH*P_gas/101325.;		// OH pressure [atm]

	STNV	 = _STNV;					// squared temperature normalized variance [-]

	if (iFractalDimension==0)
		Df = 3.0;
	else
		solve_for_fractal_dimension();		// the fractal dimension must be evaluated as soon
											// as possible
	// Collision diameters [m]
	for(int i=1; i<=N; i++)
		Dp[i] = collision_diameter(L[i]);

	omegaSoot	= fv*rhoSoot/rho_gas;
	xSoot		= omegaSoot * PM_gas/PMSoot;
	cSoot		= xSoot * P_gas/(8314.*T_gas);
	MSoot		= rhoSoot * fv;							// soot density [kg/m3]
	dp			= Lmean;								// soot particle diameter [m]
}


void OpenSMOKE_Soot_Models::coefficient_for_aggregation_kernel(const double Dp, double &Diff, double &c, double &g)
{
	double lambda;
	double Kn;
	double Cc;
	double l;
	double volume;
	double mass;

	
	// Calculation of particle volume	[m3]
	// Input:	particle diameter		[m]
	// Output:	particle volume			[m3]

		volume = pi_over_6 * (Dp*Dp*Dp);

	// Calculation of particle mass		[kg]
	// Input:	particle volume			[m3]
	//			particle density		[kg/m3]
	// Output:	particle mass			[kg]

		mass = volume * rhoSoot;


	// Calculation of mean free path	[m]
	// Input:	gas viscosity			[Pa.s]
	//			temperature				[K]
	//			molecular weight		[kg/kmol]
	//			pressure				[Pa]
	// Output:	mean free path			[m]

		lambda = 2.*mu_gas / (P_gas*pow(8.*PM_gas/(pi*8314.*T_gas), 0.50));


	// Calculation of Knudsen number	[-]
	// Input:	mean free path			[m]
	//			particle diameter		[m]
	// Output:	Knudsen number			[-]
	
		Kn = 2*lambda / Dp;

	// Calculation of correction coefficient for diffusivity	[-]
	// Input:	Knudsen number			[-]
	// Output:	Correction coefficient	[-]

		// Cc = (5+Kn*(4+Kn*(6+18*Kn)))/(5+Kn*(-1+(_pi_plus_8*Kn)));	// Artlet, 2003
		Cc = 1. + Kn*(1.257+0.40*exp(-1.10/Kn));						// Seinfeld, 2006

	// Calculation of particle diffusivity	[m2/s]
	// Input:	gas viscosity				[Pa.s]
	//			temperature					[K]
	//			particle diameter			[m]
	// Output:	particle diffusivity		[m2/s]

		Diff = Constants::kBoltzmann * T_gas /(_3pi*mu_gas*Dp) *Cc;
		
	// Calculation of collisional velocity		[m/s]
	// Input:	particle mass					[kg]
	//			temperature						[K]
	// Output:	collisional velocity			[m/s]

		c = pow(_8_over_pi*Constants::kBoltzmann*T_gas/mass, 0.50);

	// Calculation of modified mean free path	[m]
	// Input:	particle diffusivity			[m2/s]
	//			collisional velocity			[m/s]
	// Output:	modified mean free path			[m]

		l = _8_over_pi*Diff/c;

	// Calculation of modified collisional diameter [m]
	// Input:	particle diameter					[m]
	//			modified mean free path				[m]
	// Output:	modified collisional diameter		[m]
	
		g = 1./(3.*Dp*l)*(BzzPow3(Dp+l) - pow(Dp*Dp+l*l,1.50)) - Dp;
}

void OpenSMOKE_Soot_Models::aggregation_kernel(BzzMatrix &Kernel)
{
	int i, j;

	// Collision diameters [m]
	// This function is called when all the properties have been updated and therefore
	// it is not necessary to recall it

	for(i=1; i<=N; i++)
		Dp[i] = collision_diameter(L[i]);

	// These coeffients are used in the next function to estimate the aggregation kernel
	for(i=1; i<=N; i++)
		coefficient_for_aggregation_kernel(Dp[i], Diff[i], c[i], g[i]);

	// Aggregation kernel [m3/s]
	for(i=1; i<=N; i++)
		for(j=1; j<=N; j++)
			Beta[i][j] = evaluate_aggregation_kernel(Dp[i],Dp[j],Diff[i],Diff[j],c[i],c[j],g[i],g[j]);
	
	Kernel = Beta;

	Kernel *= CC_Aggregation;
}


double OpenSMOKE_Soot_Models::collision_diameter(double Lenght)
{
	// This is the collision diameter [m]
	double Dc = (iFractalDimension==0) ? Lenght : L0Soot*pow(BzzPow3(Lenght)/L0Soot_3, 1./Df);

	return Dc;
}

double OpenSMOKE_Soot_Models::evaluate_aggregation_kernel(const double Dp1,   const double Dp2, 
											    const double Diff1,	const double Diff2,
											    const double c1,	const double c2,
											    const double g1,	const double g2)
{
	// This is the evaluation of correction by Fuchs (1964)
	double A = (Dp1 + Dp2)/(Dp1+Dp2+2.*pow(g1*g1 + g2*g2, 0.50));
	double B = 8.*(Diff1+Diff2)/pow(c1*c1 + c2*c2, 0.50) / (Dp1 + Dp2);
	double correction_coefficient = 1./(A+B);

	// Aggregation Kernel
	double Kernel0 = _2pi*(Dp1+Dp2)*(Diff1+Diff2);

	// Aggregation Kernel corrected by Fuchs
	return Kernel0*correction_coefficient;
}


// **************************************************************************************** //
//							FRACTAL DIMENSION EVALUATION									//
// **************************************************************************************** //

void OpenSMOKE_Soot_Models::mean_aggregation_kernel()
{
	double DiffMean;
	double cMean;
	double gMean;
	double DpMean;

	// The definition of this mean particle lenght should be revised, it is not so clear
	// how to calculate it. I think the second choice is the correct one: it is more logic 
	// using the mean diameter evaluated on the volume (or mass) of soot particles. The
	// difference in the resulting fractal dimension is however not so large
	DpMean = collision_diameter(Lmean);
	DpMean = collision_diameter(Lmean_Volume);		// WARNING

	coefficient_for_aggregation_kernel(DpMean, DiffMean, cMean, gMean);

	// Mean Aggregation Kernel for mean particle size [m3/s]
	Beta_mean_aggregation_kernel = 
		evaluate_aggregation_kernel( DpMean,DpMean,	DiffMean,DiffMean,
									 cMean,cMean,	gMean,gMean			);
}

void OpenSMOKE_Soot_Models::evaluate_fractal_dimension()
{
	// Characteristic collision time [s]
	// It is higher for low soot concentrations and in this case the fractal dimension
	// tends to the maximum value (which is 3)
	tc  = 1./(Beta_mean_aggregation_kernel*m0);		// [s] 
	
	// Characteristic collision number [-]
	// Higher values tends to give a higher value of fractal dimension
	Tau = tc / tr_gas;											// [-]

	// Fractal dimension [-]
	// This function is 1.8 for Tau below 0.01 and 3 for Tau over 100
	// Different values are possible only if the two characteristic times differ 
	// at maximum 2 order of magnitude; the slope of the curve is not a very important
	// parameter: lower values (0.5) give a soft transition between the two 
	// asymptotic values
	if (Tau<=1)
		Df = Dfmin + pow( Df0-Dfmin, 1./pow(Tau, sFractal) );
	else
		Df = Dfmax - pow( Dfmax-Df0,    pow(Tau, sFractal)  );
}

void OpenSMOKE_Soot_Models::nls_fractal_dimension(BzzVector &x, BzzVector &f)
{
	// Recovering variables
	Df								= x[1];
	Beta_mean_aggregation_kernel	= x[2];

	// Evaluation of mean aggregation kernel and fractal dimension
	mean_aggregation_kernel();
	evaluate_fractal_dimension();
	
	// Equations to solve
	f[1] = x[1] - Df;
	f[2] = x[2] - Beta_mean_aggregation_kernel;

	// The solution must satisfy this constraints
	if (Df<1.80 || Df>3.00)
	{
		cout << "Error in solving the non linear system. Fractal Dimension = " << Df << endl;
		exit(1);
	}
}

void OpenSMOKE_Soot_Models::solve_for_fractal_dimension()
{
	// Initial values
	mean_aggregation_kernel();
	evaluate_fractal_dimension();
	
	xFirstGuess[1] = Df;
	xFirstGuess[2] = Beta_mean_aggregation_kernel;

	// Non Linear System Solution
	BzzNonLinearSystemObject o(xFirstGuess, &nls);
	o();
}


// **************************************************************************************** //
//								NUCLEATION MODELS										//
// **************************************************************************************** //

void OpenSMOKE_Soot_Models::nucleation_rate(double &J, double &epsilon)
{
	double mNucleation;

	// Maximum particle size diameter [nm]
	epsilon = 2.*L0Soot;

	// Nucleation model 1 - Liu 2001 (Combustion and Flame, 2006)
	// -----------------------------------------------------------------------------------
	if (iNucleation == 1)
	{
		mNucleation = 2.70e6*exp(-20643./T_gas)*c_C2H2;		// [kg/m3s]
		J = mNucleation * Constants::Nav_kmol / 1.08e6;		// [1/m3s]
		Sn = mNucleation / 1.08e6;							// [kmol/m3/s]
	}

	// Nucleation model 2 - Liu 2006 (Combustion and Flame, 2006)
	// -----------------------------------------------------------------------------------
	// This model has been applied to methane-air diffusion flames in laminar conditions
	// but at pressures between 5 and 40 atm
	// The particle size is 2.4 nm and the number of carbon atoms for particle is 700
	if (iNucleation == 2)
	{
		mNucleation = 24000.*exp(-16103./T_gas)*c_C2H2;		// [kg/m3s]
		J = mNucleation * Constants::Nav_kmol / 8400.;		// [1/m3s]
		Sn = mNucleation / 8400;							// [kmol/m3/s]
	}


	// Nucleation model 3 - Liu 2003 (Combustion theory and modelling, 2003)
	// -----------------------------------------------------------------------------------
	// This model has been applied to ethylene-air diffusion flames in laminar conditions
	// and atmospheric pressure
	// The particle size is 2.4 nm and the number of carbon atoms for particle is 700
	else if (iNucleation == 3)
	{
		mNucleation = 40.80*exp(-7548./T_gas)*c_C2H2;		// [kg/m3s]
		J = mNucleation * Constants::Nav_kmol / 8400.;		// [1/m3s]
		Sn = mNucleation / 8400.;							// [kmol/m3/s]
	}


	// Nucleation model 4 - Moss 1999 (Combustion and flame, 1999)
	// -----------------------------------------------------------------------------------
	// This model has been applied to turbulent methane-air jet flames
	// and atmospheric pressure and 3 atm
	// The particle size is between 0.50247 and 0.620 nm and the number of carbon 
	// atoms for particle is 12
	else if (iNucleation == 4)
	{
		mNucleation = 7776.*exp(-21100./T_gas)*c_C2H2;		// [kg/m3s]
		J = mNucleation * Constants::Nav_kmol / 144.;		// [1/m3s]
		Sn = mNucleation / 144.;							// [kmol/m3/s]
	}


	// Nucleation model 5 - Wen 2003 (Combustion and flame, 2003)
	// -----------------------------------------------------------------------------------
	// This model has been applied to turbulent kerosene-air jet flames
	// and atmospheric pressure 
	// The particle size is between 1.20 and 1.25 nm and the number of carbon 
	// atoms for particle is 100
	else if (iNucleation == 5)
	{
		mNucleation =  7776.*exp(-21100./T_gas)*c_C2H2;	// [kg/m3s]
		J =  mNucleation * Constants::Nav_kmol / 144.;	// [1/m3s]
		Sn = mNucleation / 144.;						// [kmol/m3/s]
	}


	// Nucleation model 6 - Lindstedt 1994
	// -----------------------------------------------------------------------------------
	else if (iNucleation == 6)
	{	
		mNucleation =  151200.*exp(-21100./T_gas)*c_C2H2;	// [kg/m3s]
		J = mNucleation * Constants::Nav_kmol / 720.;		// [1/m3s]
		Sn = mNucleation / 720.;							// [kmol/m3/s]
	}

		// Nucleation model 7 - Leung
	// -----------------------------------------------------------------------------------
	else if (iNucleation == 7)
	{	
		mNucleation =  243840.*exp(-21100./T_gas)*c_C2H2;	// [kg/m3s]
		J = mNucleation * Constants::Nav_kmol / 384.;		// [1/m3s]
		Sn = mNucleation / 384.;							// [kmol/m3/s]
	}

	mNucleation *= CC_Nucleation;	// [kg/m3s]
	J			*= CC_Nucleation;	// [1/m3s]
	Sn			*= CC_Nucleation;	// [kmol/m3/s]

	SN				= mNucleation;						// [kg/m3.s]
	SNGas_C2H2		= - SN/(2.*PMSoot) * PM_C2H2;		// [kg/m3.s] Acetylene consumption
	SNGas_H2		=   SN/(2.*PMSoot) * PM_H2;			// [kg/m3.s] Hydrogen formation
}


// **************************************************************************************** //
//									GROWTH MODELS										//
// **************************************************************************************** //

void OpenSMOKE_Soot_Models::growth_rate(double &alfa, double &exponent)
{	
	// Growth model 1 - Liu 2001 
	// -----------------------------------------------------------------------------------
	if (iGrowth == 1)
	{
		// mGrowth = 2.*Msoot*6.00 * exp(-6038./T) * c_C2H2 * As;					// [kg/m3]	
		// mGrowth =	       144 * exp(-6038./T) * c_C2H2 * As;					// [kg/m3]	
		
		// Perfectly spherical particles
		if (iFractalDimension == 0)
		{
			exponent = 0.;											// [-]
			alfa = 2./rhoSoot * 12000.*exp(-12083./T_gas) * c_C2H2;
		}
		// Non spherical particles
		else
		{
			// TODO
		}

		SG = 12000. * exp(-12083/T_gas) * As* c_C2H2;	// [kg/m3.s] Soot growth rate
	}

	// Growth model  - Liu 2006 (Combustion and Flame, 2006)
	// -----------------------------------------------------------------------------------
	// This model has been applied to methane-air diffusion flames in laminar conditions
	// but at pressures between 5 and 40 atm
	// This model is very different from the next ones because the dependence on the
	// surface is not linear, but the growth rate depends on the A SQUARED
	if (iGrowth == 2)
	{
		// mGrowth = 2.*Msoot * 1750. * exp(-10064./T) * c_C2H2 * sqrt(As);		// [kg/m3s]
		//         =            42000 * exp(-10064./T) * c_C2H2 * sqrt(As);		// [kg/m3s]
		
		// Perfectly spherical particles
		if (iFractalDimension == 0)
		{
			alfa =  2./rhoSoot * 42000.*exp(-10064./T_gas) * c_C2H2;
			exponent = 0.;
			alfa /= pow(pi*m2, 0.50);		// this is sqrt(As)...  As=pi*m2

			// Correction for temperature turbulent fluctuation
			if (iTfluctuations==1)
			{
				double correction = calculate_correction_coefficient(T_gas, 10064., 0., STNV);
				alfa *= correction;
			}
		}
		// Non spherical particles
		else
		{
			// TODO
		}

		SG = 42000. * exp(-10064./T_gas) * sqrt(As)* c_C2H2;		// [kg/m3.s]			
	}	
	
	// Growth model 3 - Liu 2003 (Combustion theory and modelling, 2003)
	// -----------------------------------------------------------------------------------
	// This model has been applied to ethylene-air diffusion flames in laminar conditions
	// and atmospheric pressure
	// The dependence on the surface is linear
	else if (iGrowth == 3)
	{
		// mGrowth = 2.*Msoot*6.00 * exp(-6038./T) * c_C2H2 * As;					// [kg/m3]	
		// mGrowth =	       144 * exp(-6038./T) * c_C2H2 * As;					// [kg/m3]	
		
		// Perfectly spherical particles
		if (iFractalDimension == 0)
		{
			exponent = 0.;											// [-]
			alfa = 2./rhoSoot * 144.*exp(-6038./T_gas) * c_C2H2;
		}
		// Non spherical particles
		else
		{
			// TODO
		}

		SG = 144. * exp(-6038./T_gas) * As * c_C2H2;		// [kg/m3.s]	
	}


	// Growth model 4 - Moss 1999 (Combustion and flame, 1999)
	// -----------------------------------------------------------------------------------
	// This model has been applied to turbulent methane-air jet flames
	// and atmospheric pressure and 3 atm
	// The model is always the same; the difference is the value of preexponential
	// factor and the exponent of the acetylene concentration
	else if (iGrowth == 4)
	{
		// Mass Growth Rate		
		// mGrowth = b * exp(-12100./T_gas) * pow(c_C2H2, m) * As;		// [kg/m3s]
		
		// Perfectly spherical particles
		if (iFractalDimension == 0)
		{
			exponent = 0.;						// [-]
			alfa = 2./rhoSoot *  (11700.*exp(-12100./T_gas)*c_C2H2);
		}
		// Non spherical particles
		else
		{
			// TODO
		}

		SG = 11700. * exp(-12100./T_gas) * As * c_C2H2;		// [kg/m3.s]	
	}

	// Growth model 5 - Wen 2003 (Combustion and Flame, 2003)
	// -----------------------------------------------------------------------------------
	// This model has been applied to turbulent kerosene-air jet flames
	// and atmospheric pressure 
	// The particle size is between 1.20 and 1.25 nm and the number of carbon 
	// atoms for particle is 100
	// It is very similar to the model 3 but with a different value for the constants
	else if (iGrowth == 5)
	{
		// mGrowth = 9000.6 * exp(-12100./T_gas) * c_C2H2 * As;				// [kg/m3]	
		
		// Perfectly spherical particles
		if (iFractalDimension == 0)
		{
			exponent = 0.;													// [-]
			alfa = 2./rhoSoot * 9000.6*exp(-12100./T_gas) * c_C2H2;			// [m/s]
		}

		// Non spherical particles
		else
		{
			// TODO
		}

		SG =  9000.6 * exp(-12100./T_gas) * As * c_C2H2;		// [kg/m3.s]	
	}

	// Growth model 6 - Lindstedt 1994
	// -----------------------------------------------------------------------------------
	// This model has been applied to turbulent kerosene-air jet flames
	// and atmospheric pressure 
	// The particle size is between 1.20 and 1.25 nm and the number of carbon 
	// atoms for particle is 100
	// It is very similar to the model 3 but with a different value for the constants
	else if (iGrowth == 6)
	{
		// mGrowth = 9000.6 * exp(-12100./T_gas) * c_C2H2 * As;				// [kg/m3]	
		
		// Perfectly spherical particles
		if (iFractalDimension == 0)
		{
			exponent = 0.;													// [-]
			alfa = 2./rhoSoot * 18000.*exp(-12100./T_gas) * c_C2H2;			// [m/s]
		}

		// Non spherical particles
		else
		{
			// TODO
		}

		SG =  18000. * exp(-12100./T_gas) * As * c_C2H2;		// [kg/m3.s]	
	}

	// Growth model 7 - Leung 1991 (Combustion and Flame, 1991)
	// -----------------------------------------------------------------------------------
	// The dependence on the surface is root square
	else if (iGrowth == 7)
	{
		// mGrowth = 144000.*exp(-12100./T) * c_C2H2 * sqrt(As);		// [kg/m3]	
		
		// Perfectly spherical particles
		if (iFractalDimension == 0)
		{
			alfa =  2./rhoSoot * 144000.*exp(-12100./T_gas) * c_C2H2;
			exponent = 0.;
			alfa /= pow(pi*m2, 0.50);
		}
		// Non spherical particles
		else
		{
			// TODO
		}

		SG =  144000. * exp(-12100./T_gas) * sqrt(As)* c_C2H2;		// [kg/m3.s]	
	}

	SG		*= CC_Growth;	// [kg/m3s]
	alfa	*= CC_Growth;	// [1/m3s]
	
	SGGas_C2H2	= -SG / (2.*PMSoot) * PM_C2H2;			// [kg/m3.s] Acetylene consumptioM

	SGGas_H2	=  SG / (2.*PMSoot) * PM_H2;			// [kg/m3.s] Hydrogen formation
}

// **************************************************************************************** //
//									OXIDATION MODELS										//
// **************************************************************************************** //
	
void OpenSMOKE_Soot_Models::oxydation_rate(double &alfa, double &exponent)
{
	double C_O2 = p_O2*101325./(Constants::R_J_kmol*T_gas);		// [kmol/m3]
	double C_OH = p_OH*101325./(Constants::R_J_kmol*T_gas);		// [kmol/m3]

	alfa		= 0.;
	exponent	= 0.;
	SO_O2		= 0.;
	SO_OH		= 0.;

	if (fv>=1e-9 && p_O2>=1e-5)
	{
		// Oxidation model 1 - Lee
		// -----------------------------------------------------------------------------------
		// This is a very simple model which accounts just for the oxidation by the O2
		if (iOxidation == 1)
		{
			//	mOxidation = 8903.51*exp(-19778./T)*sqrt(T)*C_O2*As;	// [kg/m3s]
			
			exponent = 0.;
			alfa =   8903.51 * exp(-19778./T_gas) * sqrt(T_gas) * C_O2;	// [kg/m2/s]
			alfa *= -2./rhoSoot;										// [m/s]		

			SO_O2 = - rhoSoot/2. * alfa * As;							// [kg/m3/s]
			SO_OH = 0.;													// [kg/m3/s]
		}

		// Oxidation model 3 - Nagle Strickland-Constable  
		// -----------------------------------------------------------------------------------
		// This is a very simple model which accounts just for the oxidation by the O2
		else if (iOxidation == 3)
		{
			double kA = 20.	  * exp(-15098./T_gas);			// [g/cm2.s.atm]		
			double kB = 4.46e-3  * exp(-7650./T_gas);		// [g/cm2.s.atm]
			double kT = 1.510e5  * exp(-48817./T_gas);		// [g/cm2.s]
			double kZ = 21.3	  * exp(2063./T_gas);		// [1/atm]
			double chi = 1./(1.+kT/(kB*p_O2));				// [-]

			// mOxidation = 120.*(kA*p_O2*chi/(1.+kZ*p_O2) + kB*p_O2*(1.-chi)) * As;
			
			exponent = 0.;
			alfa  = 12. * (kA*p_O2/(1.+kZ*p_O2)*chi + kB*p_O2*(1.-chi));	// [g/cm2.s]
			alfa *= 10.;													// [kg/m2.s]
			alfa *= 1./(1+exp(-(T_gas-1650.)/80));							// correction in Liu CST 2003
			alfa *= - 2./rhoSoot;											// [m/s] 

			SO_O2 = - rhoSoot/2. * alfa * As;								// [kg/m3.s]
			SO_OH = 0.;														// [kg/m3.s]													
		}
	}

	if (fv>=1e-9 && p_OH>=1e-7)
	{
		// Oxidation model 2 - Fenimore Jones 
		// -----------------------------------------------------------------------------------
		// This is a very simple model which accounts just for the oxidation by the OH
		// molecule; there is no activation energy
		// The efficiency is set to 0.04
		if (iOxidation == 2)
		{		
			// mOxidation = etaNeoh * 105.81 * sqrt(T) * C_OH * As;		// [kg/m3/s]
		
			exponent = 0.;
			alfa  =  0.04 * 105.81 * sqrt(T_gas) * C_OH;				// [kg/m2/s]
			alfa *= -2./rhoSoot;										// [m/s]

			SO_OH = - rhoSoot/2. * alfa * As;							// [kg/m3.s]
			SO_O2 = 0.;													// [kg/m3.s]
		}
	}

/*
	{
		const double fvstar = 1.e-9;
		const double ALFA  = 1.e-5;
		const double H = 1.50*log(ALFA/(1.-ALFA));
		const double K = 2.00*log((1.-ALFA)/ALFA) / fvstar;
		double delta = 1.e10;

		double m = ( tanh(K*fv + H) + 1. )/2.;
		double gammafv	=	m*pow(fv + m/delta, 2./3.) + 
							(1-m)*pow(fvstar, 2./3.-1.)*fv;

		double ASootCorrected = pow(36.*Constants::pi, 1./3.) * pow(m0, 1./3.) * gammafv;
		SO_O2 *= ASootCorrected/As;
		SO_OH *= ASootCorrected/As;
	}
*/

	alfa  *= CC_Oxidation;			// [m/s]		
	SO_O2 *= CC_Oxidation;			// [kg/m3/s]
	SO_OH *= CC_Oxidation;			// 

	SOGas_O2 = - SO_O2/PMSoot * 0.50 * PM_O2;						// [kg/m3/s]	O2 consumption
	SOGas_OH = - SO_OH/PMSoot * 1.00 * PM_OH;						// [kg/m3/s]	OH consumption
	SOGas_CO =   (SO_O2+SO_OH)/PMSoot * PM_CO;						// [kg/m3/s]	CO formation
	SOGas_H  =   SO_OH/PMSoot * PM_H;								// [kg/m3/s]	H  formation
}


double calculate_correction_coefficient(double T, double EoverR, double beta, double qan)
{
	double CoeffCorr;
	double n,n2,n3,n4,n5,n6;
	double EsuRT,EsuRT2,EsuRT3,EsuRT4,EsuRT5,EsuRT6;
	double sommaTer2,sommaTer4,sommaTer6;
	double qan2, qan3;
				
	EsuRT = EoverR/T;
	n = beta;
				
	n2=n*n;
	n3=n2*n;
	n4=n3*n;
	n5=n4*n;
	n6=n5*n;

	EsuRT2=EsuRT*EsuRT;
	EsuRT3=EsuRT2*EsuRT;
	EsuRT4=EsuRT3*EsuRT;
	EsuRT5=EsuRT4*EsuRT;
	EsuRT6=EsuRT5*EsuRT;
		
	qan2=qan*qan;
	qan3=qan2*qan;
	sommaTer2=-n+n2+EsuRT*(-2+2*n)+EsuRT2;
	sommaTer4=-6*n+11*n2-6*n3+n4+EsuRT*(-24+44*n-24*n2+4*n3)+EsuRT2*(36-30*n+6*n2)+EsuRT3*(-12+4*n)+EsuRT4;
	sommaTer6=-120*n+274*n2-225*n3+85*n4-15*n5+n6+EsuRT*(-720+1644*n-1350*n2+510*n3-90*n4+6*n5);
	sommaTer6+=EsuRT2*(1800-2310*n+1065*n2-210*n3+15*n4)+EsuRT3*(-1200+940*n-240*n2+20*n3);
	sommaTer6+=EsuRT4*(300-135*n+15*n2)+EsuRT5*(-30+6*n)+EsuRT6;
	CoeffCorr=1+0.25*sommaTer2*qan+0.015625*sommaTer4*qan2+4.340278E-4*sommaTer6*qan3;

	return CoeffCorr;
}

void OpenSMOKE_Soot_Models::assign_mixture(OpenSMOKE_ReactingGas &mix)
{
	iC2H2	= mix.recognize_species("C2H2");
	iO2		= mix.recognize_species("O2");
	iOH		= mix.recognize_species("OH");
	iH2		= mix.recognize_species("H2");
	iH		= mix.recognize_species("H");
	iCO		= mix.recognize_species("CO");

	PM_C2H2 = mix.M[iC2H2];
	PM_O2	= mix.M[iO2];
	PM_OH	= mix.M[iOH];
	PM_CO	= mix.M[iCO];
	PM_H2	= mix.M[iH2];
	PM_H	= mix.M[iH];

	ChangeDimensions(mix.NumberOfSpecies(), &SGas);

	SNGas_C2H2	= 0.;
	SNGas_H2	= 0.;
	SGGas_C2H2	= 0.;
	SGGas_H2	= 0.;
	SOGas_O2	= 0.;
	SOGas_OH	= 0.;
	SOGas_CO	= 0.;
	SOGas_H		= 0.;
	SN			= 0.;
	SG			= 0.;
	SO_O2		= 0.;
	SO_OH		= 0.;

}

void OpenSMOKE_Soot_Models::formation_rates()
{
	SGas = 0.;

	S = SN + SG - (SO_O2 + SO_OH);		// [kg/m3.s]

	// Gas species source terms
	// -----------------------------------------------------------------------
	S_C2H2	= SNGas_C2H2 + SGGas_C2H2;	// [kg/m3/s] Acetylene consumption
	S_H2	= SNGas_H2 + SGGas_H2;		// [kg/m3/s] Hydrogen formation
	S_OH	= SOGas_OH;					// [kg/m3/s] OH consumption
	S_O2	= SOGas_O2;					// [kg/m3/s] O2 consumption
	S_CO	= SOGas_CO;					// [kg/m3/s] CO formation
	S_H		= SOGas_H;					// [kg/m3/s] H formation

	SGas[iC2H2] = 	S_C2H2;				// [kg/m3/s] Acetylene consumption
	SGas[iH2]	= 	S_H2;				// [kg/m3/s] Hydrogen formation
	SGas[iOH]	= 	S_OH;				// [kg/m3/s] OH consumption
	SGas[iCO]	= 	S_CO;				// [kg/m3/s] CO formation
	SGas[iO2]	= 	S_O2;				// [kg/m3/s] O2 consumption
	SGas[iH]	= 	S_H;				// [kg/m3/s] H formation
}

void OpenSMOKE_Soot_Models::write_on_file(ofstream &fOutput, double csi, BzzVector &moments)
{
	fOutput << csi			<< "\t";		// soot particle number density [#/m3]
	fOutput << m0			<< "\t";		// soot particle number density [#/m3]
	fOutput << fv			<< "\t";		// soot volume fraction [-]
	fOutput << omegaSoot	<< "\t";		// soot mass fraction [-]
	fOutput << xSoot		<< "\t";		// soot mole fraction [-]
	fOutput << cSoot		<< "\t";		// soot concentration [kmol/m3]
	fOutput << MSoot		<< "\t";		// soot mass density [kg/m3]
	fOutput << dp*1.e9		<< "\t";		// soot particle diameter [nm]
	fOutput << As			<< "\t";		// soot specific area [1/m]
	
	fOutput << S			<< "\t";		// soot mass fraction source term [kg/m3.s]
	fOutput << Sn			<< "\t";		// soot particle number density source term: nucleation [kmol/m3.s]
	fOutput << SN			<< "\t";		// soot mass fraction source term: nucleation [kg/m3.s]
	fOutput << SG			<< "\t";		// soot mass fraction source term: growth [kg/m3.s]
	fOutput << SO_O2		<< "\t";		// soot mass fraction source term: O2 oxidation [kg/m3.s]
	fOutput << SO_OH		<< "\t";		// soot mass fraction source term: OH oxidation [kg/m3.s]
	
	fOutput << S_C2H2/PM_C2H2	<< "\t";	// gas source term: acetylene [kmol/m3.s]
	fOutput << S_H2/PM_H2		<< "\t";	// gas source term: hydrogen [kmol/m3.s]
	fOutput << S_OH/PM_OH		<< "\t";	// gas source term: OH [kmol/m3.s]
	fOutput << S_O2/PM_O2		<< "\t";	// gas source term: oxygen [kmol/m3.s]
	fOutput << S_CO/PM_CO		<< "\t";	// gas source term: carbon monoxide [kmol/m3.s]
	fOutput << S_H/PM_H			<< "\t";	// gas source term: H [kmol/m3.s]


	int j;

	// Normalized Moments
	for(j=1;j<=moments.Size();j++)	fOutput << moments[j]	<< "\t";
	
	// Abscissas
	for(j=1;j<=moments.Size()/2;j++)	fOutput	<< L[j]	<< "\t";

	// Weigths
	for(j=1;j<=moments.Size()/2;j++)	fOutput << w[j]	<< "\t";

	fOutput << endl;
}

void OpenSMOKE_Soot_Models::write_label_on_file(ofstream &fOutput, int N)
{
	int count = 1;

	fOutput << "csi ("			<< count++ << ")\t" << endl; 
	fOutput << "m0[#/m3] ("		<< count++ << ")\t" << endl; 
	fOutput << "fv[-] ("		<< count++ << ")\t" << endl; 
	fOutput << "omegaSoot[-] ("	<< count++ << ")\t" << endl; 
	fOutput << "xSoot[-] ("		<< count++ << ")\t" << endl; 
	fOutput << "cSoot[-] ("		<< count++ << ")\t" << endl; 
	fOutput << "MSoot[-] ("		<< count++ << ")\t" << endl; 
	fOutput << "dp [nm] ("		<< count++ << ")\t" << endl; 
	fOutput << "ASoot[1/m] ("	<< count++ << ")\t" << endl; 

	fOutput << "S[kg/m3.s] ("		<< count++ << ")\t" << endl; 
	fOutput << "Sn[kmol/m3.s] ("	<< count++ << ")\t" << endl; 
	fOutput << "SN[kg/m3.s] ("		<< count++ << ")\t" << endl; 
	fOutput << "SG[kg/m3.s] ("		<< count++ << ")\t" << endl; 
	fOutput << "SO_O2[kg/m3.s] ("	<< count++ << ")\t" << endl; 
	fOutput << "SO_OH[kg/m3.s] ("	<< count++ << ")\t" << endl; 

	fOutput << "S_C2H2[kmol/m3.s] ("	<< count++ << ")\t" << endl; 
	fOutput << "S_H2[kmol/m3.s] ("		<< count++ << ")\t" << endl; 
	fOutput << "S_OH[kmol/m3.s] ("		<< count++ << ")\t" << endl; 
	fOutput << "S_O2[kmol/m3.s] ("		<< count++ << ")\t" << endl; 
	fOutput << "S_CO[kmol/m3.s] ("		<< count++ << ")\t" << endl; 
	fOutput << "S_H[kmol/m3.s] ("		<< count++ << ")\t" << endl; 

	int j;

	for(j=1;j<=N;j++)	fOutput << "mN"		<< j << "(" << count++ << ")\t" << endl; 

	for(j=1;j<=N/2;j++)	fOutput << "csi_"	<< j << "(" << count++ << ")\t" << endl;

	for(j=1;j<=N/2;j++)	fOutput << "w_"		<< j << "(" << count++ << ")\t" << endl;

	fOutput << endl;
}
