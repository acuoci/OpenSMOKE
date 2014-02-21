/***************************************************************************
 *   Copyright (C) 2006-2009 by Alberto Cuoci						       *
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
#include <sstream>
#include "flamelet.h"
#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Utilities.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "distributions/OpenSMOKE_BetaDistribution.h"
#include "distributions/OpenSMOKE_ClippedGaussianAccurateDistribution.h"

void flamelet::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_FlameletSingle"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void flamelet::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LookUp_Table_FlameletSingle"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void flamelet::assign_flamelet(OpenSMOKE_ReactingGas &mix, vector<string> _names, const double _pressure_Pa, const double _chi, BzzVector _csi, BzzVector _temperature, BzzMatrix _massfractions)
{
	chi = _chi;
	species.names = _names;
	pressure_Pa = _pressure_Pa;
	csi = _csi;
	temperature = _temperature;
	massfractions = _massfractions;

	NP = csi.Size();
	NC = massfractions.Columns();
	clean_massfractions();

	ChangeDimensions(NC, &species.pm);
	ChangeDimensions(NP, &density);
	ChangeDimensions(NP, &udensity);
	ChangeDimensions(NP, &auxiliary_NP);
	ChangeDimensions(NP, &z_scaled);
	ChangeDimensions(NP, &zv2_scaled);

	if (temperature.Size() != NP || massfractions.Rows() != NP)
		ErrorMessage("The vector dimensions are wrong");

	// Species for soot nucleation and growth
	species.pm = mix.M;
	
	// Species for soot 
	species.iC2H2	= mix.recognize_species_without_exit("C2H2");
	species.iO2		= mix.recognize_species("O2");
	species.iOH		= mix.recognize_species("OH");
	species.iO		= mix.recognize_species("O");
	species.iH2		= mix.recognize_species("H2");
	
	species.iC6H6	= mix.recognize_species_without_exit("C6H6");
	species.iC6H5	= mix.recognize_species_without_exit("C6H5");

	// Calculate density
	for (int i=1;i<=NP;i++)
	{
		BzzVector x(NC);
		double MW;
		BzzVector aux = massfractions.GetRow(i);
		mix.GetMWAndMoleFractionsFromMassFractions(MW,x,aux);
		double cTot			= pressure_Pa / (Constants::R_J_kmol*temperature[i]);		// [kmol/m3]
			   density[i]	= cTot * MW;												// [kg/m3]
			   udensity[i]	= 1./density[i];											// [m3/kg]
	}

	for (int i=1;i<=NP;i++)
	{
		z_scaled[i]		= csi[i]*udensity[i];
		zv2_scaled[i]	= csi[i]*csi[i]*udensity[i];
	}
}

void flamelet::clean_massfractions()
{
	for(int i=1;i<=NP;i++)
		for(int j=1;j<=NC;j++)
			if (massfractions[i][j]<0.) massfractions[i][j] = 0.; 
}

void flamelet::print_on_file_appending(ofstream &fOutput)
{
	for(int j=1;j<=NP;j++)
	{
		fOutput << chi						<< "\t";
		fOutput << csi[j]					<< "\t";
		fOutput << density[j]				<< "\t";
		fOutput << temperature[j]			<< "\t";
		for(int i=1;i<=NC;i++)
			fOutput << massfractions[j][i]	<< "\t";
		fOutput << endl;
	}
	fOutput << endl;
}

void calculate_radiation(OpenSMOKE_ReactingGas &mix, const double TCell, BzzVector &mole_fractions, double &asTot, double &Qrad)
{
	double uT;
	double K_H2O, K_CO2, K_CO, K_CH4;
	BzzVector as(4);

	double	Tenv4	= BzzPow4(295.0);								// [K4]
	int		iH2O	= mix.recognize_species_without_exit("H2O");	// [-]
	int		iCO2	= mix.recognize_species_without_exit("CO2");	// [-]
	int		iCO		= mix.recognize_species_without_exit("CO");		// [-]
	int		iCH4	= mix.recognize_species_without_exit("CH4");	// [-]

	
	uT = 1000./TCell;
	
	// 1. Water [1/m.bar]
	K_H2O = -0.23093 +uT*(-1.1239+uT*(9.4153 +uT*(-2.9988 +uT*( 0.51382 + uT*(-1.8684e-5)))));
	
	// 2. Carbon Dioxide [1/m.bar]
	K_CO2 =  18.741  +uT*(-121.31+uT*(273.5  +uT*(-194.05 +uT*( 56.31   + uT*(-5.8169)))));

	// 3. Carbon monoxide [1/m.bar]
	if( TCell < 750. )	
		K_CO = 4.7869+TCell*(-0.06953 + TCell*(2.95775e-4 + TCell*(-4.25732e-7 + TCell*2.02894e-10)));
	else
		K_CO = 10.09+TCell*(-0.01183 + TCell*(4.7753e-6 + TCell*(-5.87209e-10 + TCell*-2.5334e-14)));

	// 4. Methane [1/m.bar]
	K_CH4 = 6.6334 +TCell*(- 0.0035686+TCell*(1.6682e-08+TCell*(2.5611e-10-2.6558e-14*TCell)));

	// Absorption coefficients
	if( iH2O !=0 )	as[1] = K_H2O*mole_fractions[iH2O];		// [1/m]
	if( iCO2 !=0 )	as[2] = K_CO2*mole_fractions[iCO2];		// [1/m]
	if( iCO  !=0 )	as[3] = K_CO*mole_fractions[iCO];		// [1/m]
	if( iCH4 !=0 )	as[4] = K_CH4*mole_fractions[iCH4];		// [1/m]

	asTot = as.GetSumElements();				// Absorption Coefficient

	// Source term
	Qrad = - 4.*Constants::sigma * asTot * (BzzPow4(TCell) - Tenv4); // Source term [W/m3]
}

void flamelet::print_on_file_appending_FLUENT(ofstream &fOutput, const pdf_write_type index, OpenSMOKE_ReactingGas &mix)
{
	if (index == PDF_WRITE_MF_REYNOLDS)		// 3
	{
		cout << " - Writing mixture fraction (reynolds)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << z_r_library[i][j] << endl;
	}
	
	if (index == PDF_WRITE_MFV_REYNOLDS)	// 4
	{
		cout << " - Writing variance of mixture fraction (reynolds)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << zv2_r_library[i][j] << endl;
	}

	if (index == PDF_WRITE_TEMPERATURE)		// 5
	{
		cout << " - Writing temperature (favre)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << T_f_library[i][j] << endl;
	}

	if (index == PDF_WRITE_DENSITY)			// 6
	{
		cout << " - Writing density (reynolds)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << density_r_library[i][j] << endl;
	}

	if (index == PDF_WRITE_ENTHALPY)		// 7
	{
		cout << " - Writing enthalpy (favre)..." << endl;

		BzzVector mole_fractions(NC);
		BzzVector enthalpy(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				for(int k=1;k<=NC;k++)
					mole_fractions[k] = x_f_library[i][j][k];

				enthalpy[count++] = mix.GetMixEnthalpy_Mole(T_f_library[i][j], mole_fractions)/mw_f_library[i][j];		// [J/kg]
			}
	
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
			fOutput << enthalpy[k] << endl;
	}

	if (index == PDF_WRITE_MOLECULAR_WEIGHT)	// 8
	{
		cout << " - Writing molecular weight (favre)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << mw_f_library[i][j] << endl;
	}

	if (index == PDF_WRITE_SPECIFIC_HEAT)		// 9
	{
		cout << " - Writing specific heat (favre)..." << endl;

		BzzVector mass_fractions(NC);
		BzzVector cp(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				for(int k=1;k<=NC;k++)
					mass_fractions[k] = w_f_library[i][j][k];
				mix.SpeciesCp(T_f_library[i][j]);
				cp[count++] = mix.MixCp_FromMassFractions(mass_fractions);
			}
	
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
				fOutput << cp[k] << endl;
	}

	if (index == PDF_WRITE_THERMAL_CONDUCTIVITY)	// 10
	{
		cout << " - Writing thermal conductivity (favre)..." << endl;

		BzzVector mole_fractions(NC);
		BzzVector lambda(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				for(int k=1;k<=NC;k++)
					mole_fractions[k] = x_f_library[i][j][k];
				mix.SpeciesConductivityFromFitting(T_f_library[i][j]);

				lambda[count++] = mix.MixConductivity_FromMolarFractions(mole_fractions);
			}
	
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
			fOutput << lambda[k] << endl;
	}

	if (index == PDF_WRITE_DYNAMIC_VISCOSITY)	// 11
	{
		cout << " - Writing dynamic viscosity (favre)..." << endl;

		BzzVector mole_fractions(NC);
		BzzVector mu(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				for(int k=1;k<=NC;k++)
					mole_fractions[k] = x_f_library[i][j][k];

				mix.SpeciesViscosityFromFitting(T_f_library[i][j]);

				mu[count++] = mix.MixViscosity_FromMolarFractions(mole_fractions);
			}
	
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
				fOutput << mu[k] << endl;
	}

	if (index == PDF_WRITE_ABSORPTION_COEFFICIENT)	// 12
	{
		cout << " - Writing absorption coefficient (favre)..." << endl;

		BzzVector mole_fractions(NC);
		BzzVector as(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				double asTot, Qrad;

				for(int k=1;k<=NC;k++)
					mole_fractions[k] = x_f_library[i][j][k];
				
				calculate_radiation(mix, T_f_library[i][j], mole_fractions, asTot, Qrad);
				as[count++] = asTot;
			}
		
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
			fOutput << as[k] << endl;
	}
}

void flamelet::print_on_file_appending_FLUENT(ofstream &fOutput, const int index, OpenSMOKE_ReactingGas &mix)
 {
	for(int j=1;j<=Nslices_Variance;j++)
		for(int i=1;i<=Nslices_MixtureFraction;i++)
			fOutput << w_f_library[i][j][index] << endl;
}

void flamelet::print_on_file_appending_FLUENT(BzzSave &fOutput, const pdf_write_type index, OpenSMOKE_ReactingGas &mix)
{
	if (index == PDF_WRITE_MF_REYNOLDS)		// 3
	{
		cout << " - Writing mixture fraction (reynolds)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << z_r_library[i][j];
	}
	
	if (index == PDF_WRITE_MFV_REYNOLDS)	// 4
	{
		cout << " - Writing variance of mixture fraction (reynolds)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << zv2_r_library[i][j];
	}

	if (index == PDF_WRITE_TEMPERATURE)		// 5
	{
		cout << " - Writing temperature (favre)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << T_f_library[i][j];
	}

	if (index == PDF_WRITE_DENSITY)			// 6
	{
		cout << " - Writing density (reynolds)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << density_r_library[i][j];
	}

	if (index == PDF_WRITE_ENTHALPY)		// 7
	{
		cout << " - Writing enthalpy (favre)..." << endl;

		BzzVector mole_fractions(NC);
		BzzVector enthalpy(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				for(int k=1;k<=NC;k++)
					mole_fractions[k] = x_f_library[i][j][k];

				enthalpy[count++] = mix.GetMixEnthalpy_Mole(T_f_library[i][j], mole_fractions)/mw_f_library[i][j];		// [J/kg]
			}
	
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
			fOutput << enthalpy[k];
	}

	if (index == PDF_WRITE_MOLECULAR_WEIGHT)	// 8
	{
		cout << " - Writing molecular weight (favre)..." << endl;

		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
				fOutput << mw_f_library[i][j];
	}

	if (index == PDF_WRITE_SPECIFIC_HEAT)		// 9
	{
		cout << " - Writing specific heat (favre)..." << endl;

		BzzVector mass_fractions(NC);
		BzzVector cp(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				for(int k=1;k<=NC;k++)
					mass_fractions[k] = w_f_library[i][j][k];
				mix.SpeciesCp(T_f_library[i][j]);
				cp[count++] = mix.MixCp_FromMassFractions(mass_fractions);
			}
	
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
				fOutput << cp[k];
	}

	if (index == PDF_WRITE_THERMAL_CONDUCTIVITY)	// 10
	{
		cout << " - Writing thermal conductivity (favre)..." << endl;

		BzzVector mole_fractions(NC);
		BzzVector lambda(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				for(int k=1;k<=NC;k++)
					mole_fractions[k] = x_f_library[i][j][k];
				mix.SpeciesConductivityFromFitting(T_f_library[i][j]);

				lambda[count++] = mix.MixConductivity_FromMolarFractions(mole_fractions);
			}
	
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
			fOutput << lambda[k];
	}

	if (index == PDF_WRITE_DYNAMIC_VISCOSITY)	// 11
	{
		cout << " - Writing dynamic viscosity (favre)..." << endl;

		BzzVector mole_fractions(NC);
		BzzVector mu(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				for(int k=1;k<=NC;k++)
					mole_fractions[k] = x_f_library[i][j][k];

				mix.SpeciesViscosityFromFitting(T_f_library[i][j]);

				mu[count++] = mix.MixViscosity_FromMolarFractions(mole_fractions);
			}
	
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
				fOutput << mu[k];
	}

	if (index == PDF_WRITE_ABSORPTION_COEFFICIENT)	// 12
	{
		cout << " - Writing absorption coefficient (favre)..." << endl;

		BzzVector mole_fractions(NC);
		BzzVector as(Nslices_Variance*Nslices_MixtureFraction);

		int count = 1;
		for(int j=1;j<=Nslices_Variance;j++)
			for(int i=1;i<=Nslices_MixtureFraction;i++)
			{
				double asTot, Qrad;

				for(int k=1;k<=NC;k++)
					mole_fractions[k] = x_f_library[i][j][k];
				
				calculate_radiation(mix, T_f_library[i][j], mole_fractions, asTot, Qrad);
				as[count++] = asTot;
			}
		
		for(int k=1;k<=Nslices_Variance*Nslices_MixtureFraction;k++)
			fOutput << as[k];
	}
}

void flamelet::print_on_file_appending_FLUENT(BzzSave &fOutput, const int index, OpenSMOKE_ReactingGas &mix)
 {
	for(int j=1;j<=Nslices_Variance;j++)
		for(int i=1;i<=Nslices_MixtureFraction;i++)
			fOutput << w_f_library[i][j][index];
}

void flamelet::update_mole_fractions(OpenSMOKE_ReactingGas &mix)
 {
	BzzVector mole_fractions(NC);
	BzzVector mass_fractions(NC);

	for(int j=1;j<=Nslices_Variance;j++)
		for(int i=1;i<=Nslices_MixtureFraction;i++)
		{
			for(int k=1;k<=NC;k++)
				mass_fractions[k] = w_f_library[i][j][k];
			mix.GetMWAndMoleFractionsFromMassFractions(mw_f_library[i][j], mole_fractions, mass_fractions);
		
			for(int k=1;k<=NC;k++)
				x_f_library[i][j][k] = mole_fractions[k];
		}
}


void flamelet::prepare_library(int _Nslices_Variance, double alfa)
{
	int j, k;
	maxNormalVariance = 0.98;

	// Number Of Mixture Fraction Slices
	Nslices_MixtureFraction = csi.Size();

	// Number Of Variance Slices
	Nslices_Variance = _Nslices_Variance;

	// Memory Allocations
	density_r_library	= new BzzVector[Nslices_MixtureFraction + 1];
	z_r_library			= new BzzVector[Nslices_MixtureFraction + 1];
	zv2_r_library		= new BzzVector[Nslices_MixtureFraction + 1];
	T_f_library			= new BzzVector[Nslices_MixtureFraction + 1];
	mw_f_library		= new BzzVector[Nslices_MixtureFraction + 1];
	w_f_library			= new BzzMatrix[Nslices_MixtureFraction + 1];
	x_f_library			= new BzzMatrix[Nslices_MixtureFraction + 1];
	

	for(j=1;j<=Nslices_MixtureFraction;j++)
	{
		ChangeDimensions(Nslices_Variance, &density_r_library[j]);
		ChangeDimensions(Nslices_Variance, &z_r_library[j]);
		ChangeDimensions(Nslices_Variance, &zv2_r_library[j]);
		ChangeDimensions(Nslices_Variance, &T_f_library[j]);
		ChangeDimensions(Nslices_Variance, &mw_f_library[j]);
		ChangeDimensions(Nslices_Variance, NC, &w_f_library[j]);
		ChangeDimensions(Nslices_Variance, NC, &x_f_library[j]);
	}

	ChangeDimensions(Nslices_MixtureFraction,					&mixture_fraction_space);
	ChangeDimensions(Nslices_Variance,							&normal_variance_space);
	ChangeDimensions(Nslices_MixtureFraction, Nslices_Variance, &variance_space);
	ChangeDimensions(Nslices_MixtureFraction, Nslices_Variance, &scaled_variance_space);

	// Mixture Fraction Space
	mixture_fraction_space = csi;
	
	// Variance Space
	double sum = 1.;
	for(j=1;j<=Nslices_Variance-2;j++)
		sum+=pow(alfa, j);
	double dh = maxNormalVariance/sum;

	normal_variance_space[1] = 0.;
	for(j=2;j<=Nslices_Variance;j++)
		normal_variance_space[j] = normal_variance_space[j-1] + dh*pow(alfa, j-2);

	for(k=1;k<=Nslices_MixtureFraction;k++)
		for(j=1;j<=Nslices_Variance;j++)
		{
			variance_space[k][j] = mixture_fraction_space[k]*(1.-mixture_fraction_space[k]) *
								   normal_variance_space[j];
			scaled_variance_space[k][j] = 0.25 *variance_space[k][j];
		}
}

		
void flamelet::import_pdf(BzzVector &_normal_variance, BzzMatrix &_density, BzzMatrix &_temperature, BzzMatrix *_w,
						  BzzMatrix &_z_r, BzzMatrix &_zv2_r)
{
	int j, k;

	// Number Of Mixture Fraction Slices
	Nslices_MixtureFraction = csi.Size();

	// Number Of Variance Slices
	Nslices_Variance = _normal_variance.Size();

	// Memory Allocations
	density_r_library	= new BzzVector[Nslices_MixtureFraction + 1];
	z_r_library			= new BzzVector[Nslices_MixtureFraction + 1];
	zv2_r_library		= new BzzVector[Nslices_MixtureFraction + 1];
	T_f_library			= new BzzVector[Nslices_MixtureFraction + 1];
	mw_f_library		= new BzzVector[Nslices_MixtureFraction + 1];
	w_f_library			= new BzzMatrix[Nslices_MixtureFraction + 1];
	x_f_library			= new BzzMatrix[Nslices_MixtureFraction + 1];

	for(j=1;j<=Nslices_MixtureFraction;j++)
	{
		ChangeDimensions(Nslices_Variance,		&density_r_library[j]);
		ChangeDimensions(Nslices_Variance,		&z_r_library[j]);
		ChangeDimensions(Nslices_Variance,		&zv2_r_library[j]);
		ChangeDimensions(Nslices_Variance,		&T_f_library[j]);
		ChangeDimensions(Nslices_Variance,		&mw_f_library[j]);
		ChangeDimensions(Nslices_Variance, NC,	&w_f_library[j]);
		ChangeDimensions(Nslices_Variance, NC,	&x_f_library[j]);
	}

	ChangeDimensions(Nslices_MixtureFraction,					&mixture_fraction_space);
	ChangeDimensions(Nslices_Variance,							&normal_variance_space);
	ChangeDimensions(Nslices_MixtureFraction, Nslices_Variance, &variance_space);
	ChangeDimensions(Nslices_MixtureFraction, Nslices_Variance, &scaled_variance_space);

	// Mixture Fraction Space
	mixture_fraction_space = csi;
				
	// Variance Space
	normal_variance_space = _normal_variance;
	for(k=1;k<=Nslices_MixtureFraction;k++)
		for(j=1;j<=Nslices_Variance;j++)
		{
			variance_space[k][j] = mixture_fraction_space[k]*(1.-mixture_fraction_space[k]) *
								   normal_variance_space[j];
			scaled_variance_space[k][j] = 0.25 *variance_space[k][j];
		}

	// Import Library
	for(j=1;j<=Nslices_MixtureFraction;j++)
		T_f_library[j] = _temperature.GetRow(j);
	for(j=1;j<=Nslices_MixtureFraction;j++)
		density_r_library[j] = _density.GetRow(j);
	for(j=1;j<=Nslices_MixtureFraction;j++)
		z_r_library[j] = _z_r.GetRow(j);
	for(j=1;j<=Nslices_MixtureFraction;j++)
		zv2_r_library[j] = _zv2_r.GetRow(j);
	w_f_library = _w;
}


void flamelet::calculate_betaPDF_and_mean_values(	double mean, double variance,
													double &density_mean, double &temperature_mean, 
													BzzVector &massfractions_mean,
													double &z_r_mean, double &zv2_r_mean)
{
	// Case A: AIR - No variance information is needed
	if (mean==0.)
	{
		density_mean = density[1];
		temperature_mean = temperature[1];
		z_r_mean = density_mean*z_scaled[1];
		zv2_r_mean = density_mean*zv2_scaled[1] - BzzPow2(z_r_mean);
		for(int j=1;j<=NC;j++)
			massfractions_mean[j]= massfractions[1][j];
	}

	// Case B: FUEL - No variance information is needed
	else if (mean==1.)
	{
		density_mean = density[NP];
		temperature_mean = temperature[NP];
		z_r_mean = density_mean*z_scaled[NP];
		zv2_r_mean = density_mean*zv2_scaled[NP] - BzzPow2(z_r_mean);
		for(int j=1;j<=NC;j++)
			massfractions_mean[j]= massfractions[NP][j];
	}

	// CASE C: AIR + FUEL
	else
	{
		// CASE C1: Variance=0 - BetaPDF=Delta Dirac
		if (variance==0.)
		{
			int iMean;
			double interpolationFactor;
			search_index_for_linear_interpolation(csi, mean, iMean, interpolationFactor);

			density_mean = interpolate_function(density, iMean, interpolationFactor);
			z_r_mean = density_mean*interpolate_function(z_scaled, iMean, interpolationFactor);
			zv2_r_mean = density_mean*interpolate_function(zv2_scaled, iMean, interpolationFactor) - BzzPow2(z_r_mean);
			temperature_mean = interpolate_function(temperature, iMean, interpolationFactor);
			for(int j=1;j<=NC;j++)
			{
				BzzVector aux = massfractions.GetColumn(j);	
				massfractions_mean[j]=interpolate_function(aux, iMean, interpolationFactor);
			}
		}

		// CASE C3: Calculation of betaPDF Distribution
		else
		{
			OpenSMOKE_BetaDistribution	BETA;
			double alfa = mean*(mean*(1.-mean)/variance-1.);
			double beta = (1.-mean)*(mean*(1.-mean)/variance-1.);
			BETA.Set(alfa, beta);
			
			// Mean Value of Density
			{
				LinearInterpolation interpolation;
				interpolation(csi, udensity);
				density_mean = 1./BETA.GetValue(interpolation);
			}

			// Mean Value of Csi
			{
				LinearInterpolation interpolation;
				interpolation(csi, z_scaled);
				z_r_mean = density_mean*BETA.GetValue(interpolation);
			}

			// Mean Value of CsiV2
			{
				LinearInterpolation interpolation;
				interpolation(csi, zv2_scaled);
				zv2_r_mean = density_mean*BETA.GetValue(interpolation) - BzzPow2(z_r_mean);
			}

			// Mean Value of Temperature
			{
				LinearInterpolation interpolation;
				interpolation(csi, temperature);
				temperature_mean = BETA.GetValue(interpolation);
			}

			// Mean Value of Mass Fractions
			for(int j=1;j<=NC;j++)
			{
				LinearInterpolation interpolation;
				BzzVector aux = massfractions.GetColumn(j);	
				interpolation(csi, aux);
				massfractions_mean[j] = BETA.GetValue(interpolation);
			}
		}
	}
}

void flamelet::calculate_clippedGaussian_and_mean_values
			   ( double mean, double variance, double &density_mean, 
			     double &temperature_mean, BzzVector &massfractions_mean,
				 double &z_r_mean, double &zv2_r_mean)
{
	// Case A: AIR - No variance information is needed
	if (mean==0.)
	{
		density_mean = density[1];
		z_r_mean = density_mean*z_scaled[1];
		zv2_r_mean = density_mean*zv2_scaled[1] - BzzPow2(z_r_mean);
		temperature_mean = temperature[1];
		for(int j=1;j<=NC;j++)
			massfractions_mean[j]= massfractions[1][j];
	}

	// Case B: FUEL - No variance information is needed
	else if (mean==1.)
	{
		density_mean = density[NP];
		z_r_mean = density_mean*z_scaled[NP];
		zv2_r_mean = density_mean*zv2_scaled[NP] - BzzPow2(z_r_mean);
		temperature_mean = temperature[NP];
		for(int j=1;j<=NC;j++)
			massfractions_mean[j]= massfractions[NP][j];
	}

	// CASE C: AIR + FUEL
	else
	{
		// CASE C1: Variance=0 - BetaPDF=Delta Dirac
		if (variance==0.)
		{
			int iMean;
			double interpolationFactor;
			search_index_for_linear_interpolation(csi, mean, iMean, interpolationFactor);

			density_mean = interpolate_function(density, iMean, interpolationFactor);
			z_r_mean = density_mean*interpolate_function(z_scaled, iMean, interpolationFactor);
			zv2_r_mean = density_mean*interpolate_function(zv2_scaled, iMean, interpolationFactor) - BzzPow2(z_r_mean);
			temperature_mean = interpolate_function(temperature, iMean, interpolationFactor);
			for(int j=1;j<=NC;j++)
			{
				BzzVector aux = massfractions.GetColumn(j);	
				massfractions_mean[j]=interpolate_function(aux, iMean, interpolationFactor);
			}
		}

		// CASE C2: Calculation of Clipped Gaussian Distribution 
		else
		{
			OpenSMOKE_ClippedGaussianAccurateDistribution GAUSS;
			GAUSS.Set(mean, variance);
			
			// Mean Value of Density
			{
				LinearInterpolation interpolation;
				interpolation(csi, udensity);
				density_mean = 1./GAUSS.GetValue(interpolation);
			}

			// Mean Value of Csi
			{
				LinearInterpolation interpolation;
				interpolation(csi, z_scaled);
				z_r_mean = density_mean*GAUSS.GetValue(interpolation);
			}

			// Mean Value of CsiV2
			{
				LinearInterpolation interpolation;
				interpolation(csi, zv2_scaled);
				zv2_r_mean = density_mean*GAUSS.GetValue(interpolation) - BzzPow2(z_r_mean);
			}

			// Mean Value of Temperature
			{
				LinearInterpolation interpolation;
				interpolation(csi, temperature);
				temperature_mean = GAUSS.GetValue(interpolation);
			}

			// Mean Value of Mass Fractions
			for(int j=1;j<=NC;j++)
			{
				LinearInterpolation interpolation;
				BzzVector aux = massfractions.GetColumn(j);
				interpolation(csi, aux);
				massfractions_mean[j] = GAUSS.GetValue(interpolation);
			}
		}	
	}
}

void flamelet::apply_clippedGaussian()
{
	// Apply the Clipped Gaussian PDF for each mixture fraction slice 
	for(int k=1;k<=Nslices_MixtureFraction;k++)	
	{
		cout << "Slice no. " << k << " - Mixture fraction = " << mixture_fraction_space[k] << endl;
		apply_clippedGaussian(k, mixture_fraction_space[k]);
	}
}

void flamelet::apply_betaPDF(OpenSMOKE_ReactingGas &mix)
{
	// Apply the BetaPDF for each mixture fraction slice 
	for(int k=1;k<=Nslices_MixtureFraction;k++)	
	{
		cout << "Slice no. " << k << " - Mixture fraction = " << mixture_fraction_space[k] << endl;
		apply_betaPDF(k, mixture_fraction_space[k]);
	}
/*
	bool iNO = true;
	if (iNO == true)
	{
		int index_no = mix.recognize_species_without_exit("NO");
		BzzVector omega_no(temperature.Size());
		BzzVector y(temperature.Size());
		BzzVector x(temperature.Size());
		BzzVector c(temperature.Size());
		BzzVector R(mix.NumberOfSpecies());
		
		for(int i=1;i<=omega_no.Size();i++)
		{
			double cTot, MW;

			for(int j=1;j<=NC;j++)
				y[j]= massfractions[i][j];
			
			mix.GetMWAndMoleFractionsFromMassFractions(MW,x,y);
			cTot = pressure_Pa/Constants::R_J_kmol/temperature[i];

			for(int j=1;j<=NC;j++)
				c[j]= cTot*x[j];
			

			mix.ComputeKineticParameters( temperature[i], log(temperature[i]), 1./temperature[i], pressure_Pa);
			mix.ComputeFromConcentrations( temperature[i], c, cTot, &R );	// [kmol/m3/s]
			ElementByElementProduct(R, mix.M, &R);							// [kg/m3]
			omega_no[i] = R[index_no]/density[i];
		}

	}
	*/
}

void flamelet::apply_clippedGaussian(int k, double mean)
{
	BzzVector w_f_library_auxiliary(NC);

	for(int j=1;j<=Nslices_Variance;j++)
	{
		calculate_clippedGaussian_and_mean_values(mean, variance_space[k][j], density_r_library[k][j], T_f_library[k][j], w_f_library_auxiliary, z_r_library[k][j],zv2_r_library[k][j]);
		w_f_library[k].SetRow(j, w_f_library_auxiliary); 
	}
}

void flamelet::apply_betaPDF(int k, double mean)
{
	BzzVector w_f_library_auxiliary(NC);

	for(int j=1;j<=Nslices_Variance;j++)
	{
		calculate_betaPDF_and_mean_values(mean, variance_space[k][j], density_r_library[k][j], T_f_library[k][j], w_f_library_auxiliary, z_r_library[k][j],zv2_r_library[k][j]);
		w_f_library[k].SetRow(j, w_f_library_auxiliary); 
	}
}



void flamelet::print_library(char* fileName)
{
	ofstream fOut;
	fOut.open(fileName, ios::out);
	fOut.setf(ios::scientific);
	double ZERO = 0.;

	fOut << "Csi[-](1)     ";
	fOut << "Var[-](2)     ";
	fOut << "VarN[-](3)    ";
	fOut << "VarS[-](4)    ";
	fOut << "rho[kg/m3](5) ";
	fOut << "T[K](6)       ";

	int count = 7;
	for(int j=1;j<=NC;j++)
		fOut << species.names[j] << "(" << count++ << ")\t\t";
	fOut << endl;
		

	for(int k=1;k<=Nslices_MixtureFraction;k++)
	{
		for(int j=1;j<=Nslices_Variance;j++)
		{
			fOut << (mixture_fraction_space[k]		> 1.e-16 ? mixture_fraction_space[k]	: ZERO) << "\t";
			fOut << (variance_space[k][j]			> 1.e-16 ? variance_space[k][j]			: ZERO)	<< "\t";
			fOut << (normal_variance_space[j]		> 1.e-16 ? normal_variance_space[j]		: ZERO)	<< "\t";
			fOut << (scaled_variance_space[k][j]	> 1.e-16 ? scaled_variance_space[k][j]	: ZERO)	<< "\t";
			fOut << (density_r_library[k][j]		> 1.e-16 ? density_r_library[k][j]		: ZERO)	<< "\t";
			fOut << (T_f_library[k][j]				> 1.e-16 ? T_f_library[k][j]			: ZERO)	<< "\t";
			for(int i=1;i<=NC;i++)
				fOut << (w_f_library[k][j][i]		> 1.e-16 ? w_f_library[k][j][i]			: ZERO)	<< "\t";
			fOut << endl;
		}
		fOut << endl;
	}

	fOut.close();
}

void flamelet::print_library(ofstream &fOut)
{
	int i, j, k;
	fOut << "pressure "		<< pressure_Pa				<< endl;
	fOut << "chi "			<< chi						<< endl;
	fOut << "nCsi "			<< Nslices_MixtureFraction	<< endl;
	fOut << "nVariances "	<< Nslices_Variance			<< endl;
	
	fOut << "nSpecies "		<< NC						<< endl;
	for(j=1;j<=NC;j++)
		fOut << species.names[j] << endl;
	
	fOut << "csi" << endl;
	for(k=1;k<=Nslices_MixtureFraction;k++)
		fOut << mixture_fraction_space[k]	<< endl;

	fOut << "variance" << endl;
	for(j=1;j<=Nslices_Variance;j++)
		fOut << normal_variance_space[j]	<< endl;

	fOut << "density" << endl;
	for(k=1;k<=Nslices_MixtureFraction;k++)
		for(j=1;j<=Nslices_Variance;j++)
			fOut << density_r_library[k][j]	<< endl;

	fOut << "z_reynolds" << endl;
	for(k=1;k<=Nslices_MixtureFraction;k++)
		for(j=1;j<=Nslices_Variance;j++)
			fOut << z_r_library[k][j]	<< endl;

	fOut << "zv2_reynolds" << endl;
	for(k=1;k<=Nslices_MixtureFraction;k++)
		for(j=1;j<=Nslices_Variance;j++)
			fOut << zv2_r_library[k][j]	<< endl;

	fOut << "temperature" << endl;
	for(k=1;k<=Nslices_MixtureFraction;k++)
		for(j=1;j<=Nslices_Variance;j++)
			fOut << T_f_library[k][j]			<< endl;

	fOut << "massfractions" << endl;
	for(k=1;k<=Nslices_MixtureFraction;k++)
	{
		fOut << k << endl;
		for(j=1;j<=Nslices_Variance;j++)
		{
			for(i=1;i<=NC;i++)
				fOut << w_f_library[k][j][i]	<< " ";
			fOut << endl;
		}
		fOut << endl;
	}
	fOut << "end" << endl;
	fOut << endl;
}

void flamelet::apply_source_betaPDF(OpenSMOKE_ReactingGas &mix, bool iPerfectlyCorrelatedApproach, int index_in_flamelet_group, int nSources, 
									ofstream *file_list, string *file_names, BzzVector &m0_N, BzzVector &fv_N,
									nucleation_models nucleation_model, growth_models growth_model, aggregation_models aggregation_model, oxidation_models oxidation_model)
{
	int i,k;
	double N0, Mp0, dp0;

	if (species.iC2H2 <= 0)
		ErrorMessage("Acetylene is not included in the kinetic scheme...");

	if (iPerfectlyCorrelatedApproach == true)
	{
		cout << "Perfectly Correlated Approach!" << endl;
		int a = m0_N.Size();
		int b = fv_N.Size();
				
		if (a!=b || a!=Nslices_MixtureFraction || b!=Nslices_MixtureFraction)
			ErrorMessage("The dimensions of Soot Profiles Vector and Mixture fraction space don't match...");

		if (m0_N.Max()!=1.0 || fv_N.Max()!=1.0)
			ErrorMessage("Error: the soot profiles have not been normalized...");

		cout << fv_N.Max() << endl;
	}

	bool	iDQMOM = false;
	double	mGas_H2;
	double	mGas_C2H2;
	double	mw_C2H2	= mix.M[mix.recognize_species("C2H2")];
	double	mw_H2	= mix.M[mix.recognize_species("H2")];

	sources = new sources_in_flamelet_library[nSources + 1];
	for(i=1;i<=nSources;i++)	
		sources[i].initialize(Nslices_MixtureFraction, Nslices_Variance);

	     if (nucleation_model == NUCLEATION_LIU_2001) 				N0 = 90000.;
	else if (nucleation_model == NUCLEATION_LIU_2002) 				N0 = 700.;
	else if (nucleation_model == NUCLEATION_LIU_2003) 				N0 = 700.;
	else if (nucleation_model == NUCLEATION_MOSS_1999) 				N0 = 12.;
	else if (nucleation_model == NUCLEATION_WEN_2003) 				N0 = 100.;
	else if (nucleation_model == NUCLEATION_LINDSTEDT_1994) 		N0 = 60.;
	else if (nucleation_model == NUCLEATION_LEUNG_1991) 			N0 = 32.;
	else if (nucleation_model == NUCLEATION_HALL_1997) 				N0 = 100.;
		
	dp0 = 0.27029949e-9*pow(N0, 1./3.);		// [m]
	Mp0 = Constants::Msoot * N0;			// [kg/kmol]

	// Constants
	BzzVector x(NC);								
	BzzVector c(NC);
	BzzVector p(NC);

	for(k=1;k<=Nslices_MixtureFraction;k++)	
		for(i=1;i<=Nslices_Variance;i++)
		{
			//int j;
			double MW;
			BzzVector aux = w_f_library[k].GetRow(i);
			mix.GetMWAndMoleFractionsFromMassFractions(MW,x,aux);

			double cTot =	pressure_Pa / 
							(Constants::R_J_kmol*T_f_library[k][i]);	// [kmol/m3]
			double rho = cTot * MW;										// [kg/m3]
		
			c = cTot*x;													// [kmol/m3]
			p = (pressure_Pa/Constants::P_Reference)*x;					// [atm]

			// -----------------------------------------------------------------------
			// Nucleation
			// -----------------------------------------------------------------------
			{				
				double nNucleation, mNucleation;
				double c_C2H2	= c[species.iC2H2];
				double T		= T_f_library[k][i];

				if (nucleation_model == NUCLEATION_LIU_2001)
					nNucleation = 30. * exp(-20643./T) * c_C2H2;			// [kmol/m3/s]

				else if (nucleation_model == NUCLEATION_LIU_2002)		
					nNucleation = 2.85 * exp(-16103./T) * c_C2H2;			// [kmol/m3/s]

				else if (nucleation_model == NUCLEATION_LIU_2003)		
					nNucleation = 0.004857 * exp(-7548./T) * c_C2H2;		// [kmol/m3/s]
			
				else if (nucleation_model == NUCLEATION_MOSS_1999)		
					nNucleation = 54. * exp(-21100./T) * c_C2H2;			// [kmol/m3/s]
			
				else if (nucleation_model == NUCLEATION_WEN_2003)		
					nNucleation = 54. * exp(-21100./T) * c_C2H2;			// [kmol/m3/s]

				else if (nucleation_model == NUCLEATION_LINDSTEDT_1994)		
					nNucleation = 210. * exp(-21100./T) * c_C2H2;			// [kmol/m3/s]

				else if (nucleation_model == NUCLEATION_LEUNG_1991)
					nNucleation = 635. * exp(-21100./T) * c_C2H2;			// [kmol/m3/s]

				else if (nucleation_model == NUCLEATION_HALL_1997)
				{
					if (species.iC6H6 <= 0 || species.iC6H5 <= 0)
						ErrorMessage("C6H6 and/or C6H5 are not available in the kinetic scheme...");

					double cTot		= c.GetSumElements();
					double x_C6H6	= c[species.iC6H6]/cTot;
					double x_C6H5	= c[species.iC6H5]/cTot;
					double x_H2		= c[species.iH2]/cTot;
					double x_C2H2	= c_C2H2/cTot;

					if (x_H2 > 1.e-14)
						nNucleation = 8. * cTot*cTot/Mp0*x_C6H5/x_H2*x_C2H2 * (9.63e10*x_C2H2*exp(-4378./T) + 5.62e11*x_C6H6*exp(-6390./T));
					else
						nNucleation = 0.;		// [kmol/m3/s]
				}

				else
					nNucleation = 0.;


				if (iDQMOM == false)
				{
					mNucleation =   Mp0 * nNucleation;										// [kg/m3/s] 
					mGas_C2H2	= - mNucleation/(2.*Constants::Msoot) * mw_C2H2;			// [kg/m3/s] Acetylene consumption
					mGas_H2		=   mNucleation/(2.*Constants::Msoot) * mw_H2;				// [kg/m3/s] Hydrogen formation
				}
				else
				{
					mNucleation = nNucleation * Constants::Nav_kmol;						// [1/m3/s]	
					mGas_C2H2	= - Mp0 * nNucleation/(2.*Constants::Msoot) * mw_C2H2;		// [kg/m3/s] Acetylene consumption
					mGas_H2		=   Mp0 * nNucleation/(2.*Constants::Msoot) * mw_H2;		// [kg/m3/s] Hydrogen formation
				}

				sources[1].source_mean_nofluctuations[k][i] = mNucleation;					// [kg/m3/s] 
			}				

			// -----------------------------------------------------------------------
			// Molecular Growth
			// -----------------------------------------------------------------------
			{
				double mGrowth;
				double c_C2H2	= c[species.iC2H2];
				double T		= T_f_library[k][i];


				if (growth_model == GROWTH_LIU_2001)						// linear
					mGrowth = 12000. * exp(-12083./T) * c_C2H2 ;			// [kg/m3/s] / [1/m] = [kg/m2/s]

				else if (growth_model == GROWTH_LIU_2002)					// sqrt
					mGrowth = 42000. * exp(-10064./T) * c_C2H2 ;			// [kg/m3/s] / [1/m]^0.5
	
				else if (growth_model == GROWTH_LIU_2003)					// linear
					mGrowth = 144. * exp(-6038./T) * c_C2H2 ;				// [kg/m3/s] / [1/m] = [kg/m2/s]
				
				else if (growth_model == GROWTH_MOSS_1999)					// linear
					mGrowth = 11700. * exp(-12100./T) * c_C2H2 ;			// [kg/m3/s] / [1/m] = [kg/m2/s] 

				else if (growth_model == GROWTH_WEN_2003)					// linear
					mGrowth = 9000.6 * exp(-12100./T) * c_C2H2 ;			// [kg/m3/s] / [1/m] = [kg/m2/s] 

				else if (growth_model == GROWTH_LINDSTEDT_1994)				// linear
					mGrowth = 18000. * exp(-12100./T) * c_C2H2 ;			// [kg/m3/s] / [1/m] = [kg/m2/s] 

				else if (growth_model == GROWTH_LEUNG_1991)					// sqrt
					mGrowth = 144000. * exp(-12100./T) * c_C2H2 ;			// [kg/m3/s] / [1/m]^0.5
					
				else
					mGrowth = 0.;


				if (iDQMOM == false)
				{
					mGas_C2H2 += -mGrowth / (2.*Constants::Msoot) * mw_C2H2;		// [kg/m3/s] 	Acetylene consumption
					mGas_H2   +=  mGrowth / (2.*Constants::Msoot) * mw_H2;			// [kg/m3/s] 	Hydrogen formation
				}
				else
				{
					ErrorMessage("DQMOM non yet implemented...");
				}

				sources[2].source_mean_nofluctuations[k][i] = mGrowth;

				// Perfectly Correlated Approach
				if (iPerfectlyCorrelatedApproach == true)
					sources[2].source_mean_nofluctuations[k][i] *= pow(fv_N[k], 1./3.);
			}


			// -----------------------------------------------------------------------
			// Aggregation
			// -----------------------------------------------------------------------
			{
				double nAggregation;
				double T = T_f_library[k][i];

				if (aggregation_model == AGGREGATION_SMOLUCHOWSKI)
					nAggregation =            sqrt(T)                                    ;	// [kmol/m3/s]						
					//			 = 7.94e-40 * sqrt(T) * pow(fv, 0.16667)* pow(m0, 11./6.);				
				
				else if (aggregation_model == AGGREGATION_MOSS)
					nAggregation =                     sqrt(T)          ;					// [kmol/m3/s]
					//			 = 1e9/(Nav*Nav) * sqrt(T) * (m0*m0);					
				
				else
					nAggregation = 0.;

				if (iDQMOM == false)
				{
				}
				else
				{
					ErrorMessage("DQMOM non yet implemented...");
				}

				sources[3].source_mean_nofluctuations[k][i] = nAggregation;	
				
				// Perfectly Correlated Approach
				if (iPerfectlyCorrelatedApproach == true)
					sources[3].source_mean_nofluctuations[k][i] *= pow(fv_N[k], 1./6.);

			}

			// -----------------------------------------------------------------------
			// Oxidation
			// -----------------------------------------------------------------------
			{
				double mOxidation;
				double T		= T_f_library[k][i];
				double c_O2		= c[species.iO2];
				double p_O2		= p[species.iO2];
				double c_OH		= c[species.iOH];

				if (oxidation_model == OXIDATION_LEE)			
					mOxidation = 8903.51 * exp(-19778./T) * sqrt(T) * c_O2;					// [kg/m3/s]
			
				else if (oxidation_model == OXIDATION_NSC)
				{
					double kA = 20.*exp(-15098./T);				// [g/cm2/s/atm]
					double kB = 4.46e-3*exp(-7650./T);			// [g/cm2/s/atm]
					double kT = 1.51e5*exp(-48817./T);			// [g/cm2/s]
					double kZ = 21.3*exp(2063./T);				// [1/atm]

					double chi = 1./(1.+kT/(kB*p_O2));			// [-]

					mOxidation = 120.*(kA*p_O2*chi/(1.+kZ*p_O2) + kB*p_O2*(1.-chi));		// [kg/m3/s]
				}
			
				else if (oxidation_model == OXIDATION_NEOH)	
					mOxidation = 0.04 * 105.81 * sqrt(T) * c_OH;							// [kg/m3/s]
					
				else
					mOxidation = 0.;

				if (iDQMOM == false)
				{
				}
				else
				{
					ErrorMessage("DQMOM non yet implemented...");
				}

				sources[4].source_mean_nofluctuations[k][i] = mOxidation;	
				
				// Perfectly Correlated Approach
				if (iPerfectlyCorrelatedApproach == true)
					sources[4].source_mean_nofluctuations[k][i] *= pow(fv_N[k], 2./3.);
			}			
		}

	// Source profile for variance equal to zero	
	// --------------------------------------------------------------------------------------
	for(i=1;i<=nSources;i++)	
		sources[i].source = sources[i].source_mean_nofluctuations.GetColumn(1);

	// Beta PDF Calculation of source term profiles
	// --------------------------------------------------------------------------------------
	for(i=1;i<=nSources;i++)
	{
		cout << endl;
		cout << "Source PDF No. " << i << endl;
		cout << "------------------------------------------------------------------" << endl;

		for(k=1;k<=Nslices_MixtureFraction;k++)	
		{
			cout << "Slice no. " << k << " - Mixture fraction = " << mixture_fraction_space[k] << endl;
			apply_source_betaPDF(k, mixture_fraction_space[k], sources[i].source, sources[i].source_mean);
		}
	}

	// Writing Formatted File for Comparison
	// --------------------------------------------------------------------------------------
	for(i=1;i<=nSources;i++)
	{
		char nameFile[60];
		char number[3];
		system("mkdir Output");
		strcpy(nameFile, "Output/");
		strcat(nameFile, file_names[i].c_str());
		my_itoa(index_in_flamelet_group, number, 10);
		strcat(nameFile, "_");
		strcat(nameFile, number);
		strcat(nameFile, ".out");

		sources[i].print_on_file(file_list[i], nameFile, mixture_fraction_space, normal_variance_space, chi, pressure_Pa);
	}

	cout << "End" << endl;
}

void flamelet::apply_source_betaPDF(int k, double mean, BzzVector &bzz_source, BzzMatrix &bzz_source_mean)
{
	for(int j=1;j<=Nslices_Variance;j++)
		calculate_source_betaPDF(mean, variance_space[k][j], bzz_source_mean[k][j], bzz_source);
}

void flamelet::calculate_source_betaPDF( double mean, double variance, double &source_mean, BzzVector &bzz_source)
{
	// Case A: AIR - No variance information is needed
	if (mean==0.)
		source_mean = bzz_source[1];

	// Case B: FUEL - No variance information is needed
	else if (mean==1.)
		source_mean = bzz_source[NP];

	// CASE C: AIR + FUEL
	else
	{
		// CASE C1: Variance=0 - BetaPDF=Delta Dirac
		if (variance==0.)
		{
			int iMean;
			double interpolationFactor;
			search_index_for_linear_interpolation(csi, mean, iMean, interpolationFactor);
			source_mean = interpolate_function(bzz_source, iMean, interpolationFactor);
		}

		// CASE C2: Calculation of betaPDF Distribution
		else
		{
			OpenSMOKE_BetaDistribution BETA;
			double alfa = mean*(mean*(1.-mean)/variance-1.);
			double beta = (1.-mean)*(mean*(1.-mean)/variance-1.);
			BETA.Set(alfa, beta);
			
			// Mean Value of Density
			LinearInterpolation interpolation;
			interpolation(csi, bzz_source);
			source_mean = BETA.GetValue(interpolation);
		}
	}
}


void flamelet::requests_for_soot(	OpenSMOKE_ReactingGas &mix, int csiRequest, int varianceRequest,
									BzzVector &wRequest, BzzVector &xRequest, 
									BzzVector &cRequest, double &rhoRequest, double &temperatureRequest)
{
	int k=csiRequest;	
	int i=varianceRequest;

	wRequest = w_f_library[k].GetRow(i);

	double MW;
	mix.GetMWAndMoleFractionsFromMassFractions(MW,xRequest, wRequest);
	double cTot = pressure_Pa / (Constants::R_J_kmol*T_f_library[k][i]);		// [kmol/m3]
	rhoRequest = cTot * MW;														// [kg/m3]
	cRequest = cTot*xRequest;													// [kmol/m3]
	temperatureRequest = T_f_library[k][i];										// [K]
}

void flamelet::apply_source_betaPDF_SootProperties(int index_in_flamelet_group, int nSources, ofstream *file_list, string *file_names,
												   BzzVector &m0, BzzVector &fv)
{
	int i,k;

	sources = new sources_in_flamelet_library[nSources + 1];
	for(i=1;i<=nSources;i++)	
		sources[i].initialize(Nslices_MixtureFraction, Nslices_Variance);

	for(k=1;k<=Nslices_MixtureFraction;k++)	
		for(i=1;i<=Nslices_Variance;i++)
		{
			// -----------------------------------------------------------------------
			// Soot Particle Density
			// -----------------------------------------------------------------------
			sources[1].source_mean_nofluctuations[k][i] = m0[k];

			// -----------------------------------------------------------------------
			// Soot Volume fraction
			// -----------------------------------------------------------------------	
			sources[2].source_mean_nofluctuations[k][i] = fv[k];
		}

	// Source profile for variance equal to zero	
	// --------------------------------------------------------------------------------------
	for(i=1;i<=nSources;i++)	
		sources[i].source = sources[i].source_mean_nofluctuations.GetColumn(1);

	// Beta PDF Calculation of source term profiles
	// --------------------------------------------------------------------------------------
	for(i=1;i<=nSources;i++)
	{
		cout << endl;
		cout << "Source PDF No. " << i << endl;
		cout << "------------------------------------------------------------------" << endl;

		for(k=1;k<=Nslices_MixtureFraction;k++)	
		{
			cout << "Slice no. " << k << " - Mixture fraction = " << mixture_fraction_space[k] << endl;
			apply_source_betaPDF(k, mixture_fraction_space[k], sources[i].source, sources[i].source_mean);
		}
	}

	// Writing Formatted File for Comparison
	// --------------------------------------------------------------------------------------
	for(i=1;i<=nSources;i++)
	{
		char nameFile[30];
		char number[3];
		strcpy(nameFile, "Output/");
		strcat(nameFile, file_names[i].c_str());
		my_itoa(index_in_flamelet_group, number, 10);
		strcat(nameFile, number);
		strcat(nameFile, ".out");

		sources[i].print_on_file(file_list[i], nameFile, mixture_fraction_space, normal_variance_space, chi, pressure_Pa);
	}
}

void flamelet::write_progress_variables(OpenSMOKE_ReactingGas &mix, const int nProgressVariable, vector<string>& progress_variable_name, BzzVectorInt* progress_variable_index, BzzVector* progress_variable_value)
{
	stringstream index;
	index << chi;
	string file_name = "Flamelet_" + index.str() + ".out";

	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);
	fOutput.setf(ios::scientific);

	int fOutputCount = 1;
	PrintTagOnGnuplotLabel(20, fOutput, "csi[-]", fOutputCount);
	PrintTagOnGnuplotLabel(20, fOutput, "rho[kg/m3]", fOutputCount);
	PrintTagOnGnuplotLabel(20, fOutput, "T[K]", fOutputCount);
	for(int j=1;j<=nProgressVariable;j++)
		PrintTagOnGnuplotLabel(20, fOutput, progress_variable_name[j], fOutputCount);
	for(int j=1;j<=nProgressVariable;j++)
		PrintTagOnGnuplotLabel(20, fOutput, "S-"+progress_variable_name[j], fOutputCount);
	for(int j=1;j<=NC;j++)
		PrintTagOnGnuplotLabel(20, fOutput, mix.names[j]+"_w", fOutputCount);
	fOutput << endl;
	
	for (int i=1;i<=NP;i++)
	{
		// Memory allocation
		BzzVector x(mix.NumberOfSpecies());
		BzzVector c(mix.NumberOfSpecies());
		BzzVector R(mix.NumberOfSpecies());
		
		// Reconstruct source terms
		double MW;
		BzzVector aux=massfractions.GetRow(i);
		mix.GetMWAndMoleFractionsFromMassFractions(MW,x,aux);
		double cTot = pressure_Pa / (Constants::R_J_kmol*temperature[i]);		// [kmol/m3]
		c = cTot*x;																// [kmol/m3]
		mix.ComputeKineticParameters( temperature[i], log(temperature[i]), 1./temperature[i], pressure_Pa);
		mix.ComputeFromConcentrations( temperature[i], c, cTot, &R);			// [kmol/m3/s]
		ElementByElementProduct(R, mix.M, &R);									// [kg/m3/s]
			
		// Write on file
		fOutput << setw(20) << left << csi[i];
		fOutput << setw(20) << left << density[i];
		fOutput << setw(20) << left << temperature[i];
		for(int j=1;j<=nProgressVariable;j++)
		{
			double sum = 0.;
			for(int k=1;k<=progress_variable_index[j].Size();k++)
				sum += massfractions[i][progress_variable_index[j][k]]*progress_variable_value[j][k];
			fOutput << setw(20) << left << sum;
		}
		for(int j=1;j<=nProgressVariable;j++)
		{
			double sum = 0.;
			for(int k=1;k<=progress_variable_index[j].Size();k++)
				sum += R[progress_variable_index[j][k]]*progress_variable_value[j][k];
			fOutput << setw(20) << left << sum;
		}
		for(int j=1;j<=NC;j++)
			fOutput << setw(20) << left << massfractions[i][j];

		fOutput << endl;
	}
}
