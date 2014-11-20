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

#include <iomanip>
#include "addons/OpenSMOKE_SootManager.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"

/*
void OpenSMOKE_SootManager::recognizeSpecies(int numComponents, std::string *names, BzzVector &M)
{
	int i;
	double pi = acos(-1.);
	nClassBin = 16;

	classBin = new BzzVectorInt[nClassBin+1];


	// Recognizing soot
	int indexLarge = 0;
	int indexSmall = 0;
	for(i=1;i<=numComponents;i++)
		if (names[i].compare(0, 3, "BIN") == 0)
		{
			if ( 
				 (names[i]=="BIN1")   ||
				 (names[i]=="BIN1A")  ||
				 (names[i]=="BIN1B")  ||
				 (names[i]=="BIN1C")  ||
				 (names[i]=="BIN2A")  ||
				 (names[i]=="BIN2B")  ||
				 (names[i]=="BIN2C")  ||
				 (names[i]=="BIN3A")  ||
				 (names[i]=="BIN3B")  ||
				 (names[i]=="BIN3C")  ||
				 (names[i]=="BIN4A")  ||
				 (names[i]=="BIN4B")  ||
				 (names[i]=="BIN4C")  ||
				 (names[i]=="BIN1AJ")  ||
				 (names[i]=="BIN1BJ")  ||
				 (names[i]=="BIN1CJ")  ||
				 (names[i]=="BIN2AJ")  ||
				 (names[i]=="BIN2BJ")  ||
				 (names[i]=="BIN2CJ")  ||
				 (names[i]=="BIN3AJ")  ||
				 (names[i]=="BIN3BJ")  ||
				 (names[i]=="BIN3CJ")  ||
				 (names[i]=="BIN4AJ")  ||
				 (names[i]=="BIN4BJ")  ||
				 (names[i]=="BIN4CJ")  
			   )
			{
				indexSmall++;
				smallBins.Append(i);															// Bin Global Index
				small_MW.Append(M[i]);															// Molecular Weight [kg/kmol]
				small_m.Append(M[i]/Constants::Nav_kmol);										// Mass of particle [kg]
				small_d.Append( pow(6./pi*M[i]/1.8/(Constants::Nav_kmol/1000.), 1./3.) *1e7 );	// Diameter of particle [nm] 
			}
			else
			{
				indexLarge++;
				largeBins.Append(i);															// Bin Global Index
				large_MW.Append(M[i]);															// Molecular Weight [kg/kmol]
				large_m.Append(M[i]/Constants::Nav_kmol);										// Mass of particle [kg]
				large_d.Append( pow(6./pi*M[i]/1.8/(Constants::Nav_kmol/1000.), 1./3.) *1e7 );	// Diameter of particle [nm] 

				//if (indexLarge<=39)
				{
				if		(names[i].compare(0, 4, "BIN5") == 0)	classBin[1].Append(indexLarge);
				else if (names[i].compare(0, 4, "BIN6") == 0)	classBin[2].Append(indexLarge);
				else if (names[i].compare(0, 4, "BIN7") == 0)	classBin[3].Append(indexLarge);
				else if (names[i].compare(0, 4, "BIN8") == 0)	classBin[4].Append(indexLarge);
				else if (names[i].compare(0, 4, "BIN9") == 0)	classBin[5].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN10") == 0)	classBin[6].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN11") == 0)	classBin[7].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN12") == 0)	classBin[8].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN13") == 0)	classBin[9].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN14") == 0)	classBin[10].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN15") == 0)	classBin[11].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN16") == 0)	classBin[12].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN17") == 0)	classBin[13].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN18") == 0)	classBin[14].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN19") == 0)	classBin[15].Append(indexLarge);
				else if (names[i].compare(0, 5, "BIN20") == 0)	classBin[16].Append(indexLarge);
				else
				{
					cout << "This BIN is not expected: " << names[i] << endl;
					system("pause");
					exit(1);
				}
				}

			}
		}

		ChangeDimensions(nClassBin, &dClassBin);
		ChangeDimensions(nClassBin, &mClassBin);
		ChangeDimensions(nClassBin, &mwClassBin);
		for(i=1;i<=nClassBin;i++)
		{
			for(int j=1;j<=classBin[i].Size();j++)
			{
				dClassBin[i] += large_d[ classBin[i][j]];
				mClassBin[i] += large_m[ classBin[i][j]];
				mwClassBin[i] += large_MW[ classBin[i][j]];
			}

			if (classBin[i].Size()>0)
			{
				dClassBin[i] /= classBin[i].Size();
				mClassBin[i] /= classBin[i].Size();
				mwClassBin[i] /= classBin[i].Size();
			}
			else
			{
				dClassBin[i]  = 0.;
				mClassBin[i]  = 0.;
				mwClassBin[i] = 0.;
			}

		}

		bool iBINsExistency = false;
		for(i=1;i<=nClassBin;i++)
			if (classBin[i].Size() != 0)
			{
				iBINsExistency = true;
				break;
			}

		if(iBINsExistency == true)
		{
			cout << "-------------------------------------------------" << endl;
			cout << " BIN Classes" << endl;
			cout << "-------------------------------------------------" << endl;
			for(i=1;i<=nClassBin;i++)
				if (classBin[i].Size()) cout << " Class " << i << "(" << i+4 << "): " << classBin[i].Size() << " BINs " << dClassBin[i] << endl;
			cout << endl << endl;
		}

		if (smallBins.Size()==0 && largeBins.Size()==0)
			iSoot = false;										// no soot in kinetic scheme
		else iSoot = true;										//    soot in kinetic scheme

		if (iSoot == true)
		{
			cout << "Number of Small BIN Classes: " << smallBins.Size() << endl;
			cout << "Number of Large BIN Classes: " << largeBins.Size() << endl;

			for(i=1;i<=smallBins.Size();i++)
			{
				for(int j=1;j<=smallBins.Size();j++)
					if (i!=j)
						if (small_MW[i]==small_MW[j])
							cout << "WARNING: Two Small BIN Classes have the same molecular weight! " <<
									names[smallBins[i]] << "\t" << names[smallBins[j]] << endl; 
			}	

			for(i=1;i<=largeBins.Size();i++)
			{
				for(int j=1;j<=largeBins.Size();j++)
					if (i!=j)
						if (large_MW[i]==large_MW[j])
							cout << "WARNING: Two Large BIN Classes have the same molecular weight! " <<
									 names[largeBins[i]] << "\t" << names[largeBins[j]] << endl; 
			}

			cout.setf(ios::scientific);
			cout << "Small BINs" << endl;
			cout << "-----------------------------------------------------------------------------------------" << endl;
			for(i=1;i<=smallBins.Size();i++)
				cout << i << "\t" << names[smallBins[i]]		<< "\t" 
								  << smallBins[i]				<< "\t" 
								  << small_MW[i]				<< "\t" 
								  << small_m[i]					<< "\t" 
								  << small_d[i]					<< endl; 
			
			cout << endl;
			cout << "Large BINs" << endl;
			cout << "-----------------------------------------------------------------------------------------" << endl;
			for(i=1;i<=largeBins.Size();i++)
				cout << i << "\t"	<< names[largeBins[i]]		<< "\t" 
									<< largeBins[i]				<< "\t" 
									<< large_MW[i]				<< "\t" 
									<< large_m[i]				<< "\t" 
									<< large_d[i]				<< endl; 
		}
}

double OpenSMOKE_SootManager::giveSumSmallBins(BzzVector &vector)
{
	double sum = 0.;
	for(int i=1;i<=smallBins.Size();i++)
		sum += vector[smallBins[i]];
	return sum;
}

double OpenSMOKE_SootManager::giveSumLargeBins(BzzVector &vector)
{
	double sum = 0.;
	for(int i=1;i<=largeBins.Size();i++)
		sum += vector[largeBins[i]];
	return sum;
}

double OpenSMOKE_SootManager::giveVolumeFractionFromMassFractionForLargeBins(double large_omega_tot_Gas, double rhoGas)
{
	return large_omega_tot_Gas * rhoGas/Constants::rhoSoot;
}

double OpenSMOKE_SootManager::giveVolumeFractionFromMassFractionForSmallBins(double small_omega_tot_Gas, double rhoGas)
{
	return small_omega_tot_Gas * rhoGas/Constants::rhoSoot;
}

BzzVector OpenSMOKE_SootManager::giveParticleNumberDensityFromMassFractionForSmallBins(BzzVector &omega, double rhoGas)
{
	BzzVector N(smallBins.Size());
	for(int i=1;i<=smallBins.Size();i++)
		N[i] = omega[smallBins[i]] / small_m[i] * rhoGas;
	return N;
}

BzzVector OpenSMOKE_SootManager::giveParticleNumberDensityFromMassFractionForLargeBins(BzzVector &omega, double rhoGas)
{
	BzzVector N(largeBins.Size());
	for(int i=1;i<=largeBins.Size();i++)
		N[i] = omega[largeBins[i]] / large_m[i] * rhoGas;
	return N;
}

BzzVector OpenSMOKE_SootManager::giveParticleMassFractionFromMassFractionForSmallBins(BzzVector &omegaGas, double rhoGas)
{
	BzzVector w(smallBins.Size());
	for(int i=1;i<=smallBins.Size();i++)
		w[i] = omegaGas[smallBins[i]];
	return w;
}

BzzVector OpenSMOKE_SootManager::giveParticleMassFractionFromMassFractionForLargeBins(BzzVector &omegaGas, double rhoGas)
{
	BzzVector w(largeBins.Size());
	for(int i=1;i<=largeBins.Size();i++)
		w[i] = omegaGas[largeBins[i]];
	return w;
}

BzzVector OpenSMOKE_SootManager::giveParticleVolumeFractionFromMassFractionForSmallBins(BzzVector &omegaGas, double rhoGas)
{
	BzzVector w(smallBins.Size());
	for(int i=1;i<=smallBins.Size();i++)
		w[i]= omegaGas[smallBins[i]] * rhoGas / Constants::rhoSoot;
	return w;
}

BzzVector OpenSMOKE_SootManager::giveParticleVolumeFractionFromMassFractionForLargeBins(BzzVector &omegaGas, double rhoGas)
{
	BzzVector w(largeBins.Size());
	for(int i=1;i<=largeBins.Size();i++)
		w[i]= omegaGas[largeBins[i]] * rhoGas / Constants::rhoSoot;
	return w;
}

BzzVector OpenSMOKE_SootManager::giveParticleSizeDistributionForSmallBins()
{
	double Ntot = small_N.GetSumElements();

	BzzVector p(smallBins.Size());
	for(int i=1;i<=smallBins.Size();i++)
		p[i] = small_N[i] / Ntot;
	return p;
}

BzzVector OpenSMOKE_SootManager::giveParticleSizeDistributionForLargeBins()
{
	double Ntot = large_N.GetSumElements();

	BzzVector p(largeBins.Size());
	for(int i=1;i<=largeBins.Size();i++)
		p[i] = large_N[i] / Ntot;
	return p;
}



// Standard Moments
// -------------------------------------------------------------------------------------------------------
double OpenSMOKE_SootManager::giveNOrderStandardMomentForLargeBins(BzzVector &x, BzzVector &f, double N)
{
	double fTot = f.GetSumElements();
	double sum = 0;
	for(int i=1;i<=largeBins.Size();i++)
		sum += pow(x[i], N)*f[i] / fTot;
	return sum;
}

double OpenSMOKE_SootManager::giveNOrderStandardMomentForSmallBins(BzzVector &x, BzzVector &f, double N)
{
	double fTot = f.GetSumElements();
	double sum = 0;
	for(int i=1;i<=smallBins.Size();i++)
		sum += pow(x[i], N)*f[i] / fTot;
	return sum;
}


// Central Moments
// -------------------------------------------------------------------------------------------------------
double OpenSMOKE_SootManager::giveNOrderCentralMomentForLargeBins(BzzVector &x, BzzVector &f, double N)
{
	int i;
	
	double fTot = f.GetSumElements();

	double mu = 0.;
	for(i=1;i<=largeBins.Size();i++)
		mu += x[i]*f[i]/fTot;
	
	double sum = 0;
	for(i=1;i<=largeBins.Size();i++)
		sum += pow(x[i] - mu, N)*f[i]/fTot;
	return sum;
}

double OpenSMOKE_SootManager::giveNOrderCentralMomentForSmallBins(BzzVector &x, BzzVector &f, double N)
{
	int i;

	double fTot = f.GetSumElements();

	double mu = 0.;
	for(i=1;i<=smallBins.Size();i++)
		mu += x[i]*f[i]/fTot;
	
	double sum = 0;
	for(i=1;i<=smallBins.Size();i++)
		sum += pow(x[i] - mu, N)*f[i]/fTot;
	return sum;
}

// Normalized Moments
// -------------------------------------------------------------------------------------------------------
double OpenSMOKE_SootManager::giveNOrderNormalizedMomentForLargeBins(BzzVector &x, BzzVector &f, double N)
{
	double mN = giveNOrderCentralMomentForLargeBins(x, f, N);
	double sigma = giveStandardDeviationForLargeBins(x, f);

	return mN / pow(sigma,N);
}

double OpenSMOKE_SootManager::giveNOrderNormalizedMomentForSmallBins(BzzVector &x, BzzVector &f, double N)
{
	double mN = giveNOrderCentralMomentForSmallBins(x, f, N);
	double sigma = giveStandardDeviationForSmallBins(x, f);

	return mN / pow(sigma,N);
}

// Mean
// -------------------------------------------------------------------------------------------------------
double OpenSMOKE_SootManager::giveMeanForLargeBins(BzzVector &x, BzzVector &f)
{
	double sum = 0;
	for(int i=1;i<=largeBins.Size();i++)
		sum += x[i]*f[i];
	return sum / f.GetSumElements();
}

double OpenSMOKE_SootManager::giveMeanForSmallBins(BzzVector &x, BzzVector &f)
{
	double sum = 0;
	for(int i=1;i<=smallBins.Size();i++)
		sum += x[i]*f[i];
	return sum / f.GetSumElements();
}

// Standard Deviation
// -------------------------------------------------------------------------------------------------------
double OpenSMOKE_SootManager::giveStandardDeviationForLargeBins(BzzVector &x, BzzVector &f)
{
	return sqrt(giveNOrderCentralMomentForLargeBins(x, f, 2.));
}

double OpenSMOKE_SootManager::giveStandardDeviationForSmallBins(BzzVector &x, BzzVector &f)
{
	return sqrt(giveNOrderCentralMomentForSmallBins(x, f, 2.));
}

// Variance
// -------------------------------------------------------------------------------------------------------
double OpenSMOKE_SootManager::giveVarianceForLargeBins(BzzVector &x, BzzVector &f)
{
	return giveNOrderCentralMomentForLargeBins(x, f, 2.);
}

double OpenSMOKE_SootManager::giveVarianceForSmallBins(BzzVector &x, BzzVector &f)
{
	return giveNOrderCentralMomentForSmallBins(x, f, 2.);
}

// Skewness
// -------------------------------------------------------------------------------------------------------
double OpenSMOKE_SootManager::giveSkewnessForLargeBins(BzzVector &x, BzzVector &f)
{
	return giveNOrderCentralMomentForLargeBins(x, f, 3.);
}

double OpenSMOKE_SootManager::giveSkewnessForSmallBins(BzzVector &x, BzzVector &f)
{
	return giveNOrderCentralMomentForSmallBins(x, f, 3.);
}

// Kurtoisis
// -------------------------------------------------------------------------------------------------------
double OpenSMOKE_SootManager::giveKurtoisisForLargeBins(BzzVector &x, BzzVector &f)
{
	return giveNOrderNormalizedMomentForLargeBins(x, f, 4.) - 3.;
}

double OpenSMOKE_SootManager::giveKurtoisisForSmallBins(BzzVector &x, BzzVector &f)
{
	return giveNOrderNormalizedMomentForSmallBins(x, f, 4.) - 3.;
}


double OpenSMOKE_SootManager::giveMeanDiameterForLargeBins(BzzVector &xGas)
{
	double sum = 0.;
	double dmean = 0.;
	
	for(int i=1;i<=largeBins.Size();i++)
		if (xGas[largeBins[i]] > 1e-100)
		{
			dmean	+= xGas[largeBins[i]] * large_d[i];
			sum		+= xGas[largeBins[i]];
		}

	if (sum <= 1.e-12) 
		return 0.;
	else
		return dmean / sum;
}

double OpenSMOKE_SootManager::giveMeanDiameterForSmallBins(BzzVector &xGas)
{
	double sum = 0.;
	double dmean = 0.;
	
	for(int i=1;i<=smallBins.Size();i++)
		if (xGas[smallBins[i]] > 1e-100)
		{
			dmean	+= xGas[smallBins[i]] * small_d[i];
			sum		+= xGas[smallBins[i]];
		}

	if (sum <= 1.e-12) 
		return 0.;
	else
		return dmean / sum;
}

void OpenSMOKE_SootManager::large_calculateAll(BzzVector &omegaGas, BzzVector &xGas, 
											BzzVector &cGas, double rhoGas)
{
	large_omega_tot		= giveSumLargeBins(omegaGas);
	large_x_tot			= giveSumLargeBins(xGas);
	large_c_tot			= giveSumLargeBins(cGas);

	large_fv_tot		= giveVolumeFractionFromMassFractionForLargeBins(large_omega_tot, rhoGas);

	large_N				= giveParticleNumberDensityFromMassFractionForLargeBins(omegaGas, rhoGas);
	large_N_tot			= large_N.GetSumElements();
	
	large_omega			= giveParticleMassFractionFromMassFractionForLargeBins(omegaGas, rhoGas);
	large_fv			= giveParticleVolumeFractionFromMassFractionForLargeBins(omegaGas, rhoGas);
	mean_large_d		= giveMeanDiameterForLargeBins(xGas);


//	if (large_fv_tot > 1.e-12)
//		large_statistical_analysis();
}

void OpenSMOKE_SootManager::small_calculateAll(BzzVector &omegaGas, BzzVector &xGas, 
											BzzVector &cGas, double rhoGas)
{
	small_omega_tot		= giveSumSmallBins(omegaGas);
	small_x_tot			= giveSumSmallBins(xGas);
	small_c_tot			= giveSumSmallBins(cGas);

	small_fv_tot		= giveVolumeFractionFromMassFractionForSmallBins(small_omega_tot, rhoGas);

	small_N				= giveParticleNumberDensityFromMassFractionForSmallBins(omegaGas, rhoGas);
	small_N_tot			= small_N.GetSumElements();

	small_omega			= giveParticleMassFractionFromMassFractionForSmallBins(omegaGas, rhoGas);
	small_fv			= giveParticleVolumeFractionFromMassFractionForSmallBins(omegaGas, rhoGas);

	mean_small_d		= giveMeanDiameterForSmallBins(xGas);

//	if (small_fv_tot > 1.e-12)
//		small_statistical_analysis();
}

void OpenSMOKE_SootManager::large_statistical_analysis()
{
	large_p = giveParticleSizeDistributionForLargeBins();


	mean_large_d		= giveMeanForLargeBins(large_d, large_fv);
	mean_large_d32		= giveNOrderStandardMomentForLargeBins(large_d, large_fv, 3.) / giveNOrderStandardMomentForLargeBins(large_d, large_fv, 2.);
	std_large_d			= giveStandardDeviationForLargeBins(large_d, large_fv);
	variance_large_d	= giveVarianceForLargeBins(large_d, large_fv);
	skewness_large_d	= giveSkewnessForLargeBins(large_d, large_fv);
	kurtoisis_large_d	= giveKurtoisisForLargeBins(large_d, large_fv);

	mean_large_m		= giveMeanForLargeBins(large_m, large_fv);
	mean_large_m32		= giveNOrderStandardMomentForLargeBins(large_m, large_fv, 3.) / giveNOrderStandardMomentForLargeBins(large_m, large_fv, 2.);
	std_large_m			= giveStandardDeviationForLargeBins(large_m, large_fv);
	variance_large_m	= giveVarianceForLargeBins(large_m, large_fv);
	skewness_large_m	= giveSkewnessForLargeBins(large_m, large_fv);
	kurtoisis_large_m	= giveKurtoisisForLargeBins(large_m, large_fv);

	mean_large_MW		= giveMeanForLargeBins(large_MW, large_fv);
	mean_large_MW32		= giveNOrderStandardMomentForLargeBins(large_MW, large_fv, 3.) / giveNOrderStandardMomentForLargeBins(large_MW, large_fv, 2.);
	std_large_MW		= giveStandardDeviationForLargeBins(large_MW, large_fv);
	variance_large_MW	= giveVarianceForLargeBins(large_MW, large_fv);
	skewness_large_MW	= giveSkewnessForLargeBins(large_MW, large_fv);
	kurtoisis_large_MW	= giveKurtoisisForLargeBins(large_MW, large_fv);
}

void OpenSMOKE_SootManager::small_statistical_analysis()
{
	small_p = giveParticleSizeDistributionForSmallBins();

	mean_small_d		= giveMeanForSmallBins(small_d, small_fv);
	mean_small_d32		= giveNOrderStandardMomentForSmallBins(small_d, small_fv, 3.) / giveNOrderStandardMomentForSmallBins(small_d, small_fv, 2.);
	std_small_d			= giveStandardDeviationForSmallBins(small_d, small_fv);
	variance_small_d	= giveVarianceForSmallBins(small_d, small_fv);
	skewness_small_d	= giveSkewnessForSmallBins(small_d, small_fv);
	kurtoisis_small_d	= giveKurtoisisForSmallBins(small_d, small_fv);

	mean_small_m		= giveMeanForSmallBins(small_m, small_fv);
	mean_small_m32		= giveNOrderStandardMomentForSmallBins(small_m, small_fv, 3.) / giveNOrderStandardMomentForSmallBins(small_m, small_fv, 2.);
	std_small_m			= giveStandardDeviationForSmallBins(small_m, small_fv);
	variance_small_m	= giveVarianceForSmallBins(small_m, small_fv);
	skewness_small_m	= giveSkewnessForSmallBins(small_m, small_fv);
	kurtoisis_small_m	= giveKurtoisisForSmallBins(small_m, small_fv);

	mean_small_MW		= giveMeanForSmallBins(small_MW, small_fv);
	mean_small_MW32		= giveNOrderStandardMomentForSmallBins(small_MW, small_fv, 3.) / giveNOrderStandardMomentForSmallBins(small_MW, small_fv, 2.);
	std_small_MW		= giveStandardDeviationForSmallBins(small_MW, small_fv);
	variance_small_MW	= giveVarianceForSmallBins(small_MW, small_fv);
	skewness_small_MW	= giveSkewnessForSmallBins(small_MW, small_fv);
	kurtoisis_small_MW	= giveKurtoisisForSmallBins(small_MW, small_fv);
}

void OpenSMOKE_SootManager::print_main_data_on_file(ofstream &fSoot)
{			
	// LARGE SOOT
	fSoot << setw(20) << left << large_fv_tot;			// fv
	fSoot << setw(20) << left << large_x_tot;			// xSoot
	fSoot << setw(20) << left << large_omega_tot;		// omegaSoot
	fSoot << setw(20) << left << large_c_tot;			// cSoot [kmol/m3]
	fSoot << setw(20) << left << large_N_tot;			// soot particle number density [#/m3]

	fSoot << setw(20) << left << mean_large_d;			// soot particle mean diameter [nm]
	fSoot << setw(20) << left << mean_large_d32;		// soot particle mean diameter [nm]
	fSoot << setw(20) << left << std_large_d;			// soot particle standard deviation [nm]
	fSoot << setw(20) << left << variance_large_d;		// soot particle variance  [nm2]
	fSoot << setw(20) << left << skewness_large_d;		// soot particle skewness  [nm3]
	fSoot << setw(20) << left << kurtoisis_large_d;		// soot particle kurtoisis [-]

	fSoot << setw(20) << left << mean_large_m;			// soot particle mean mass [kg]
	fSoot << setw(20) << left << mean_large_m32;		// soot particle mean mass [kg]
	fSoot << setw(20) << left << std_large_m;			// soot particle standard deviation [kg]
	fSoot << setw(20) << left << variance_large_m;		// soot particle variance  [kg2]
	fSoot << setw(20) << left << skewness_large_m;		// soot particle skewness  [kg3]
	fSoot << setw(20) << left << kurtoisis_large_m;		// soot particle kurtoisis [-]

	fSoot << setw(20) << left << mean_large_MW;			// soot particle mean molecular weight [kg/kmol]
	fSoot << setw(20) << left << mean_large_MW32;		// soot particle mean molecular weight [kg/kmol]
	fSoot << setw(20) << left << std_large_MW;			// soot particle standard deviation [kg/kmol]
	fSoot << setw(20) << left << variance_large_MW;		// soot particle variance  [kg2/kmol22]
	fSoot << setw(20) << left << skewness_large_MW;		// soot particle skewness  [kg3/kmol3]
	fSoot << setw(20) << left << kurtoisis_large_MW;		// soot particle kurtoisis [-]


	// SMALL SOOT
	fSoot << setw(20) << left << small_fv_tot;		// fv
	fSoot << setw(20) << left << small_x_tot;		// xSoot
	fSoot << setw(20) << left << small_omega_tot;	// omegaSoot
	fSoot << setw(20) << left << small_c_tot;		// cSoot [kmol/m3]
	fSoot << setw(20) << left << small_N_tot;		// soot particle number density [#/m3]
				
	fSoot << setw(20) << left << mean_small_d;			// soot particle mean diameter [nm]
	fSoot << setw(20) << left << mean_small_d32;		// soot particle mean diameter [nm]
	fSoot << setw(20) << left << std_small_d;			// soot particle standard deviation [nm]
	fSoot << setw(20) << left << variance_small_d;		// soot particle variance  [nm2]
	fSoot << setw(20) << left << skewness_small_d;		// soot particle skewness  [nm3]
	fSoot << setw(20) << left << kurtoisis_small_d;		// soot particle kurtoisis [-]

	fSoot << setw(20) << left << mean_small_m;			// soot particle mean mass [kg]
	fSoot << setw(20) << left << mean_small_m32;		// soot particle mean mass [kg]
	fSoot << setw(20) << left << std_small_m;			// soot particle standard deviation [kg]
	fSoot << setw(20) << left << variance_small_m;		// soot particle variance  [kg2]
	fSoot << setw(20) << left << skewness_small_m;		// soot particle skewness  [kg3]
	fSoot << setw(20) << left << kurtoisis_small_m;		// soot particle kurtoisis [-]

	fSoot << setw(20) << left << mean_small_MW;			// soot particle mean molecular weight [kg/kmol]
	fSoot << setw(20) << left << mean_small_MW32;		// soot particle mean molecular weight [kg/kmol]
	fSoot << setw(20) << left << std_small_MW;			// soot particle standard deviation [kg/kmol]
	fSoot << setw(20) << left << variance_small_MW;		// soot particle variance  [kg2/kmol22]
	fSoot << setw(20) << left << skewness_small_MW;		// soot particle skewness  [kg3/kmol3]
	fSoot << setw(20) << left << kurtoisis_small_MW;	// soot particle kurtoisis [-]
}

void OpenSMOKE_SootManager::GnuPlotInterface(ofstream &fSoot, int count)
{

	PrintTagOnGnuplotLabel(20, fSoot, "Lfv",			count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lx",			count);
	PrintTagOnGnuplotLabel(20, fSoot, "Ly",			count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lc[kmol/m3]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "LN[#/m3]",		count);

	PrintTagOnGnuplotLabel(20, fSoot, "Ld[nm]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Ld32[nm]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Ldstd[nm]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Ldvar[nm2]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Ldske[nm3]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Ldkur[-]",		count);

	PrintTagOnGnuplotLabel(20, fSoot, "Lm[kg]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lm32[kg]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lmstd[kg]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lmvar[kg2]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lmske[kg3]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lmkur[-]",		count);

	PrintTagOnGnuplotLabel(20, fSoot, "L[kg/kml]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "L32[g/ml]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lstd[g/ml]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lva[g2/ml2]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lsk[g3/ml3]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Lkur[-]",			count);



	PrintTagOnGnuplotLabel(20, fSoot, "Sfv",			count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sx",			count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sy",			count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sc[kmol/m3]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "SN[#/m3]",		count);

	PrintTagOnGnuplotLabel(20, fSoot, "Sd[nm]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sd32[nm]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sdstd[nm]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sdvar[nm2]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sdske[nm3]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sdkur[-]",		count);

	PrintTagOnGnuplotLabel(20, fSoot, "Sm[kg]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sm32[kg]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Smstd[kg]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Smvar[kg2]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Smske[kg3]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Smkur[-]",		count);

	PrintTagOnGnuplotLabel(20, fSoot, "S[g/ml]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "S32[g/ml]",		count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sstd[g/ml]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Sva[g2/ml2]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Ssk[g3/ml3]",	count);
	PrintTagOnGnuplotLabel(20, fSoot, "Skur[-]",			count);
}

void OpenSMOKE_SootManager::print_distribution_data_on_file_large(double csi, ofstream &fSoot)
{	
	int i, j;
	BzzVector p(nClassBin);
	BzzVector fv(nClassBin);
	BzzVector omega(nClassBin);
	BzzVector N(nClassBin);

	if (large_fv_tot > 1.e-12)
	{
	for(i=1;i<=nClassBin;i++)
	{
		for(j=1;j<=classBin[i].Size();j++)
		{
			p[i]		+= large_p[ classBin[i][j]  ];
			omega[i]	+= large_omega[ classBin[i][j]  ];
			fv[i]		+= large_fv[ classBin[i][j]  ];
			N[i]		+= large_N[ classBin[i][j]  ];
		}
	}

	double meanD_mole	= 0.;
	double sum_mole		= 0.;
	double meanD_mass	= 0.;
	double sum_mass		= 0.;
	double meanD_moleRatio	= 0.;

	for(i=1;i<=nClassBin;i++)
	{
		meanD_mole	+= N[i] * dClassBin[i];
		sum_mole	+= N[i];
		meanD_mass	+= omega[i] * dClassBin[i];
		sum_mass	+= omega[i];

		meanD_moleRatio	+= 3.14159/4 * N[i] * dClassBin[i] * dClassBin[i];

	}
	
	if (sum_mass > 1.e-12)
	{
		meanD_moleRatio = meanD_moleRatio / meanD_mole;

		meanD_mole		= meanD_mole/sum_mole;
		meanD_mass		= meanD_mass/sum_mass;
	}
	else
	{
		meanD_mole = 0.;
		meanD_mass = 0.;
		meanD_moleRatio = 0.;
	}


	for(i=1;i<=nClassBin;i++)
	{
		fSoot << setw(20) << left << csi;
		
		fSoot << setw(20) << left << dClassBin[i];
		fSoot << setw(20) << left << mClassBin[i];
		fSoot << setw(20) << left << mwClassBin[i];

		fSoot << setw(20) << left << p[i];
		fSoot << setw(20) << left << N[i];
		fSoot << setw(20) << left << fv[i];
		fSoot << setw(20) << left << omega[i];

		fSoot << setw(20) << left << meanD_mole;
		fSoot << setw(20) << left << meanD_mass;
		fSoot << setw(20) << left << meanD_moleRatio;
		
		fSoot << endl;
	}
	}
}
*/


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

void OpenSMOKE_PAHManager::recognizeSpecies(int numComponents, std::string *names, BzzVector &M)
{
	int i;

	// Recognizing PAH
	for(i=1;i<=numComponents;i++)
	{
		if (names[i].compare(0, 3, "BIN") == 0)
		{
			if ( 
				 (names[i]=="BIN1")   ||
				 (names[i]=="BIN1A")  ||
				 (names[i]=="BIN1B")  ||
				 (names[i]=="BIN1C")  ||
				 (names[i]=="BIN2A")  ||
				 (names[i]=="BIN2B")  ||
				 (names[i]=="BIN2C")  ||
				 (names[i]=="BIN3A")  ||
				 (names[i]=="BIN3B")  ||
				 (names[i]=="BIN3C")  ||
				 (names[i]=="BIN4A")  ||
				 (names[i]=="BIN4B")  ||
				 (names[i]=="BIN4C")  ||
				 (names[i]=="BIN1AJ")  ||
				 (names[i]=="BIN1BJ")  ||
				 (names[i]=="BIN1CJ")  ||
				 (names[i]=="BIN2AJ")  ||
				 (names[i]=="BIN2BJ")  ||
				 (names[i]=="BIN2CJ")  ||
				 (names[i]=="BIN3AJ")  ||
				 (names[i]=="BIN3BJ")  ||
				 (names[i]=="BIN3CJ")  ||
				 (names[i]=="BIN4AJ")  ||
				 (names[i]=="BIN4BJ")  ||
				 (names[i]=="BIN4CJ")  
			   )
			{
				pah500nm.Append(i);				// Bin Global Index
				pah500nm_MW.Append(M[i]);		// Molecular Weight [kg/kmol]
			}
		}
		
		if ( 
			 (names[i] == "FENA")  ||
			 (names[i] == "PYRE")  
		   )
		{
			pah400nm.Append(i);				// Bin Global Index
			pah400nm_MW.Append(M[i]);		// Molecular Weight [kg/kmol]
		}
		if ( 
			(names[i] == "BENZ")  ||
			(names[i] == "TOLUO") ||
			(names[i] == "FC2H")  ||
			(names[i] == "INDE")  ||
			(names[i] == "FUEL1") ||
			(names[i] == "C12H8") ||
			(names[i] == "FEN2")  ||
			(names[i] == "FLUO")  
		   )
		{
			pah340nm.Append(i);				// Bin Global Index
			pah340nm_MW.Append(M[i]);		// Molecular Weight [kg/kmol]
		}
	}

	if (pah340nm.Size()==0 && pah400nm.Size()==0 && pah500nm.Size()==0)
		 iPAH = false;										// no PAH in kinetic scheme
	else iPAH = true;										//    PAH in kinetic scheme

	if (iPAH == true)
	{
		cout << "Number of 340nm PAH Classes: " << pah340nm.Size() << endl;
		cout << "Number of 400nm PAH Classes: " << pah400nm.Size() << endl;
		cout << "Number of 500nm PAH Classes: " << pah500nm.Size() << endl;

		cout.setf(ios::scientific);
			
		cout << "340nm PAH" << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl;
		for(i=1;i<=pah340nm.Size();i++)
			cout << i << "\t" << names[pah340nm[i]]			<< "\t" 
							  << pah340nm[i]				<< "\t" 
							  << pah340nm_MW[i]				<< endl; 
		cout << endl;

		cout << "400nm PAH" << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl;
		for(i=1;i<=pah400nm.Size();i++)
			cout << i << "\t" << names[pah400nm[i]]			<< "\t" 
							  << pah400nm[i]				<< "\t" 
							  << pah400nm_MW[i]				<< endl; 
		cout << endl;

		cout << "500nm PAH" << endl;
		cout << "-----------------------------------------------------------------------------------------" << endl;
		for(i=1;i<=pah500nm.Size();i++)
			cout << i << "\t" << names[pah500nm[i]]			<< "\t" 
							  << pah500nm[i]				<< "\t" 
							  << pah500nm_MW[i]				<< endl; 
		cout << endl;
	}
}

double OpenSMOKE_PAHManager::giveSum340nm(BzzVector &vector)
{
	double sum = 0.;
	for(int i=1;i<=pah340nm.Size();i++)
		sum += vector[pah340nm[i]];
	return sum;
}

double OpenSMOKE_PAHManager::giveSum400nm(BzzVector &vector)
{
	double sum = 0.;
	for(int i=1;i<=pah400nm.Size();i++)
		sum += vector[pah400nm[i]];
	return sum;
}

double OpenSMOKE_PAHManager::giveSum500nm(BzzVector &vector)
{
	double sum = 0.;
	for(int i=1;i<=pah500nm.Size();i++)
		sum += vector[pah500nm[i]];
	return sum;
}

void OpenSMOKE_PAHManager::pah_calculateAll(BzzVector &omegaGas, BzzVector &xGas, BzzVector &cGas, double rhoGas)
{
	pah340nm_omega_tot	= giveSum340nm(omegaGas);
	pah340nm_x_tot		= giveSum340nm(xGas);
	pah340nm_c_tot		= giveSum340nm(cGas);

	pah400nm_omega_tot	= giveSum400nm(omegaGas);
	pah400nm_x_tot		= giveSum400nm(xGas);
	pah400nm_c_tot		= giveSum400nm(cGas);

	pah500nm_omega_tot	= giveSum500nm(omegaGas);
	pah500nm_x_tot		= giveSum500nm(xGas);
	pah500nm_c_tot		= giveSum500nm(cGas);

}

void OpenSMOKE_PAHManager::GnuPlotInterface(ofstream &fPAH, int count)
{
	PrintTagOnGnuplotLabel(20, fPAH, "x340",			count);
	PrintTagOnGnuplotLabel(20, fPAH, "x400",			count);
	PrintTagOnGnuplotLabel(20, fPAH, "x500",			count);
	PrintTagOnGnuplotLabel(20, fPAH, "w340",			count);
	PrintTagOnGnuplotLabel(20, fPAH, "w400",			count);
	PrintTagOnGnuplotLabel(20, fPAH, "w500",			count);
	PrintTagOnGnuplotLabel(20, fPAH, "c340[kmol/m3]",	count);
	PrintTagOnGnuplotLabel(20, fPAH, "c400[kmol/m3]",	count);
	PrintTagOnGnuplotLabel(20, fPAH, "c500[kmol/m3]",	count);
}


void OpenSMOKE_PAHManager::print_main_data_on_file(ofstream &fPAH)
{			
	fPAH << setw(20) << left << pah340nm_x_tot;		// mole fraction
	fPAH << setw(20) << left << pah400nm_x_tot;		// mole fraction
	fPAH << setw(20) << left << pah500nm_x_tot;		// mole fraction

	fPAH << setw(20) << left << pah340nm_omega_tot;		// mass fraction
	fPAH << setw(20) << left << pah400nm_omega_tot;		// mass fraction
	fPAH << setw(20) << left << pah500nm_omega_tot;		// mass fraction

	fPAH << setw(20) << left << pah340nm_c_tot;		// concentration [kmol/m3]
	fPAH << setw(20) << left << pah400nm_c_tot;		// concentration [kmol/m3]
	fPAH << setw(20) << left << pah500nm_c_tot;		// concentration [kmol/m3]
}