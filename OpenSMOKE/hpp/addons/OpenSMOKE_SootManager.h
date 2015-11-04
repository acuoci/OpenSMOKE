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

#ifndef OPENSMOKE_SOOTMANAGER
#define OPENSMOKE_SOOTMANAGER

#include "BzzMath.hpp"
using namespace std;

class OpenSMOKE_ReactingGas;

/*
class OpenSMOKE_SootManager
{

private:

	// Small BIN data
	BzzVectorInt smallBins;
	BzzVector small_MW;
	BzzVector small_m;
	BzzVector small_d;
	BzzVector small_N;
	BzzVector small_omega;
	BzzVector small_fv;
	BzzVector small_p;


	// Large BIN data
	BzzVectorInt largeBins;
	BzzVector large_MW;
	BzzVector large_m;
	BzzVector large_d;
	BzzVector large_N;
	BzzVector large_omega;
	BzzVector large_fv;
	BzzVector large_p;


	BzzVectorInt *classBin;
	int nClassBin;
	BzzVector dClassBin;
	BzzVector mClassBin;
	BzzVector mwClassBin;



	// Central Moments
	// -------------------------------------------------------------------------------------------------------
	double giveNOrderStandardMomentForLargeBins(BzzVector &x, BzzVector &f, double N);
	double giveNOrderStandardMomentForSmallBins(BzzVector &x, BzzVector &f, double N);

	// Central Moments
	// -------------------------------------------------------------------------------------------------------
	double giveNOrderCentralMomentForLargeBins(BzzVector &x, BzzVector &f, double N);
	double giveNOrderCentralMomentForSmallBins(BzzVector &x, BzzVector &f, double N);

	// Normalized Moments
	// -------------------------------------------------------------------------------------------------------
	double giveNOrderNormalizedMomentForLargeBins(BzzVector &x, BzzVector &f, double N);
	double giveNOrderNormalizedMomentForSmallBins(BzzVector &x, BzzVector &f, double N);

	// Mean
	// -------------------------------------------------------------------------------------------------------
	double giveMeanForLargeBins(BzzVector &x, BzzVector &f);
	double giveMeanForSmallBins(BzzVector &x, BzzVector &f);

	// Standard Deviation
	// -------------------------------------------------------------------------------------------------------
	double giveStandardDeviationForLargeBins(BzzVector &x, BzzVector &f);
	double giveStandardDeviationForSmallBins(BzzVector &x, BzzVector &f);

	// Variance
	// -------------------------------------------------------------------------------------------------------
	double giveVarianceForLargeBins(BzzVector &x, BzzVector &f);
	double giveVarianceForSmallBins(BzzVector &x, BzzVector &f);

	// Skewness
	// -------------------------------------------------------------------------------------------------------
	double giveSkewnessForLargeBins(BzzVector &x, BzzVector &f);
	double giveSkewnessForSmallBins(BzzVector &x, BzzVector &f);

	// Kurtoisis
	// -------------------------------------------------------------------------------------------------------
	double giveKurtoisisForLargeBins(BzzVector &x, BzzVector &f);
	double giveKurtoisisForSmallBins(BzzVector &x, BzzVector &f);

	double mean_large_d;
	double mean_large_d32;
	double std_large_d;
	double variance_large_d;
	double skewness_large_d;
	double kurtoisis_large_d;

	double mean_small_d;
	double mean_small_d32;
	double std_small_d;
	double variance_small_d;
	double skewness_small_d;
	double kurtoisis_small_d;

	double mean_large_m;
	double mean_large_m32;
	double std_large_m;
	double variance_large_m;
	double skewness_large_m;
	double kurtoisis_large_m;

	double mean_small_m;
	double mean_small_m32;
	double std_small_m;
	double variance_small_m;
	double skewness_small_m;
	double kurtoisis_small_m;

	double mean_large_MW;
	double mean_large_MW32;
	double std_large_MW;
	double variance_large_MW;
	double skewness_large_MW;
	double kurtoisis_large_MW;

	double mean_small_MW;
	double mean_small_MW32;
	double std_small_MW;
	double variance_small_MW;
	double skewness_small_MW;
	double kurtoisis_small_MW;

	double giveMeanDiameterForLargeBins(BzzVector &xGas);
	double giveMeanDiameterForSmallBins(BzzVector &xGas);


public:

	// Index to check if the soot bins are available in the current kinetic scheme
	bool iSoot;

	// Soot CRECK manager setup
	void	recognizeSpecies(int NC, std::string *names, BzzVector &M);

	// General purpose functions
	double	giveSumSmallBins(BzzVector &vector);
	double	giveSumLargeBins(BzzVector &vector);

	// Soot volume fraction functions
	double	giveVolumeFractionFromMassFractionForLargeBins(double large_omega_tot_Gas, double rhoGas);
	double  giveVolumeFractionFromMassFractionForSmallBins(double small_omega_tot_Gas, double rhoGas);

	// Soot particle number density
	BzzVector giveParticleNumberDensityFromMassFractionForSmallBins(BzzVector &omega, double rhoGas);
	BzzVector giveParticleNumberDensityFromMassFractionForLargeBins(BzzVector &omega, double rhoGas);

	// Soot particle mass fraction
	BzzVector giveParticleMassFractionFromMassFractionForSmallBins(BzzVector &omegaGas, double rhoGas);
	BzzVector giveParticleMassFractionFromMassFractionForLargeBins(BzzVector &omegaGas, double rhoGas);

	// Soot Volume Fraction
	BzzVector giveParticleVolumeFractionFromMassFractionForSmallBins(BzzVector &omegaGas, double rhoGas);
	BzzVector giveParticleVolumeFractionFromMassFractionForLargeBins(BzzVector &omegaGas, double rhoGas);

	// Soot calculate every information
	void large_calculateAll(BzzVector &omegaGas, BzzVector &xGas, BzzVector &cGas, double rhoGas);
	void small_calculateAll(BzzVector &omegaGas, BzzVector &xGas, BzzVector &cGas, double rhoGas);

	// Soot Statistical Analysis
	void large_statistical_analysis();
	void small_statistical_analysis();

	// Soot Particle Size Distribution
	BzzVector giveParticleSizeDistributionForSmallBins();
	BzzVector giveParticleSizeDistributionForLargeBins();



	// Soot Distribution properties - Large particles
	// --------------------------------------------------------------------------------------------------------------
	double	large_omega_tot;		// soot mass fraction
	double	large_x_tot;			// soot mole fraction
	double	large_c_tot;			// soot concentration
	double	large_fv_tot;			// soot volume fraction
	double	large_N_tot;			// soot particle number density



	// Soot Distribution properties - Small particles
	// --------------------------------------------------------------------------------------------------------------
	double	small_omega_tot;		// soot mass fraction
	double	small_x_tot;			// soot mole fraction
	double	small_c_tot;			// soot concentration
	double	small_fv_tot;			// soot volume fraction
	double	small_N_tot;			// soot particle number density



	void print_main_data_on_file(ofstream &fSoot);
	void GnuPlotInterface(ofstream &fSoot, int count);

	void print_distribution_data_on_file_large(double csi, ofstream &fSoot);

};
*/
class OpenSMOKE_PAHManager
{

private:

	// PAH 340nm
	BzzVectorInt pah340nm;
	BzzVector pah340nm_MW;

	// PAH 400nm
	BzzVectorInt pah400nm;
	BzzVector pah400nm_MW;

	// PAH 500nm
	BzzVectorInt pah500nm;
	BzzVector pah500nm_MW;

public:

	// Index to check if the soot bins are available in the current kinetic scheme
	bool iPAH;

	// Soot CRECK manager setup
	void	recognizeSpecies(int NC, std::string *names, BzzVector &M);

	// General purpose functions
	double	giveSum340nm(BzzVector &vector);
	double	giveSum400nm(BzzVector &vector);
	double	giveSum500nm(BzzVector &vector);

	// Calculate every information
	void pah_calculateAll(BzzVector &omegaGas, BzzVector &xGas, BzzVector &cGas, double rhoGas);


	// PAH 340nm
	// --------------------------------------------------------------------------------------------------------------
	double	pah340nm_omega_tot;		// pah mass fraction
	double	pah340nm_x_tot;			// pah mole fraction
	double	pah340nm_c_tot;			// pah concentration

	// PAH 400nm
	// --------------------------------------------------------------------------------------------------------------
	double	pah400nm_omega_tot;		// pah mass fraction
	double	pah400nm_x_tot;			// pah mole fraction
	double	pah400nm_c_tot;			// pah concentration

	// PAH 500nm
	// --------------------------------------------------------------------------------------------------------------
	double	pah500nm_omega_tot;		// pah mass fraction
	double	pah500nm_x_tot;			// pah mole fraction
	double	pah500nm_c_tot;			// pah concentration


	void print_main_data_on_file(ofstream &fPAH);
	void GnuPlotInterface(ofstream &fPAH, int count);
};


#endif // OPENSMOKE_SOOTMANAGER
