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
#include "addons/OpenSMOKE_PolimiSoot.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"

void OpenSMOKE_PolimiSoot::Setup(OpenSMOKE_IdealGas &gas, const std::string minimum_bin)
{
	const unsigned int bin_index_zero  = 10;
	const unsigned int bin_index_final = 20;
	const double bin_density_A   = 1500.;
	const double bin_density_B   = 1700.;
	const double Df = 1.8;
	Setup(gas, minimum_bin, bin_index_zero, bin_density_A, bin_index_final, bin_density_B, Df);	
}


void OpenSMOKE_PolimiSoot::Setup(OpenSMOKE_IdealGas &gas, const std::string minimum_bin, 
									const unsigned int bin_index_zero, const double bin_density_A, 
									const unsigned int bin_index_final, const double bin_density_B, const double Df)
{
	iBin_ = false;
	iSoot_ = false;
	bin_names_.push_back("bin_names");

	fv_small_		= 0.;
	rho_small_		= 0.;
	N_small_		= 0.;
	omega_small_	= 0.;
	x_small_		= 0.;
	fv_large_		= 0.;
	rho_large_		= 0.;
	N_large_		= 0.;
	omega_large_	= 0.;
	x_large_		= 0.;
	h_over_c_small_ = 0.;
	o_over_c_small_ = 0.;
	o_over_h_small_ = 0.;
	h_over_c_large_ = 0.;
	o_over_c_large_ = 0.;
	o_over_h_large_ = 0.;
	dmean_N_large_	= 0.;
	dmean_N_small_	= 0.;
	d32_N_large_	= 0.;
	d32_N_small_	= 0.;
	dvariance_N_	= 0.;
	dstd_N_			= 0.;

	bin_index_zero_ = bin_index_zero;
	bin_density_A_ = bin_density_A;
	bin_index_final_ = bin_index_final;
	bin_density_B_ = bin_density_B;
	Df_ = Df;

	for(int j=1;j<=gas.NumberOfSpecies();j++)
		if (gas.names[j].compare(0, 3, "BIN") == 0)
		{
			iBin_ = true;
			break;
		}

	double min_mw_soot = 1.e32;
	for(int j=1;j<=gas.NumberOfSpecies();j++)
		if (gas.names[j].compare(0, minimum_bin.size(), minimum_bin) == 0)
		{
			iSoot_ = true;
			if (gas.M[j] < min_mw_soot)	min_mw_soot = gas.M[j];
		}

	if (iBin_ == true)
	{
		int iC = gas.recognize_element("c");
		int iH = gas.recognize_element("h");
		int iO = gas.recognize_element("o");
		unsigned int BIN12_index = 0;

		for(int j=1;j<=gas.NumberOfSpecies();j++)
			if (gas.names[j].compare(0, 3, "BIN") == 0)
			{
				// BIN density
				{
					const int nc = gas.elements(iC,j);
					const double index = log(nc/24.)/log(2.) + 1.;
					if (index<=bin_index_zero_) bin_density_.Append(bin_density_A_);
					else
					{
						const double m = (bin_density_B_-bin_density_A_)/(bin_index_final_-bin_index_zero_);
						const double c = bin_density_A_-m*bin_index_zero_;
						bin_density_.Append(c+m*index);
					}
				}

				if (gas.names[j].compare(0, 5, "BIN12") == 0 && BIN12_index == 0)
					BIN12_index = bin_indices_.Size()+1;

				bin_indices_.Append(j);																								// Index of bin in the gas phase kinetic scheme [-]
				bin_names_.push_back(gas.names[j]);																					// Index of bin in the gas phase kinetic scheme [-]
				bin_mw_.Append(gas.M[j]);																							// Molecular weight [kg/kmol]
				bin_m_.Append(gas.M[j]/Constants::Nav_kmol);																		// Mass of particle [kg]
				bin_ds_.Append( pow(6./Constants::pi*gas.M[j]/(bin_density_[bin_density_.Size()]/1000.)
										/(Constants::Nav_mol), 1./3.)*1.e-2 );														// Diameter of particle [m]
				bin_V_.Append( Constants::pi/6.*BzzPow3(bin_ds_[bin_ds_.Size()]) );													// Volume of particle [m3]
				bin_c_.Append(gas.elements(iC,j));															// C
				bin_h_.Append(gas.elements(iH,j));															// H
				bin_o_.Append(gas.elements(iO,j));															// O
				bin_h_over_c_.Append(gas.elements(iH,j)/gas.elements(iC,j));								// Ratio H/C
				bin_o_over_c_.Append(gas.elements(iO,j)/gas.elements(iC,j));								// Ratio O/C
				if (gas.elements(iH,j)!=0)	bin_o_over_h_.Append(gas.elements(iO,j)/gas.elements(iH,j));	// Ratio O/H
				else                        bin_o_over_h_.Append(0.);										// Ratio O/H
				
				if ( gas.M[j] >=  min_mw_soot)	bin_indices_large_.Append(bin_indices_.Size());
				else                            bin_indices_small_.Append(bin_indices_.Size());

				// Collisional diameter and diameter
				{
					int iCollisional = false;
					const int nc = gas.elements(iC, j);
					const double index = log(nc / 24.) / log(2.) + 1;
					if (index > 12)
						iCollisional = true;
					
					if (iCollisional == true)
					{
						bin_np_.Append(bin_m_[bin_m_.Size()] / bin_m_[BIN12_index]);
						bin_dc_.Append(sqrt(5. / 3.)*bin_ds_[BIN12_index] * pow(bin_np_[bin_np_.Size()] / pow(1. + 2. / Df_, Df_ / 2.), 1. / Df_));
						bin_d_.Append(bin_dc_[bin_dc_.Size()]);
					}
					else
					{
						bin_np_.Append(0.);
						bin_dc_.Append(0.);
						bin_d_.Append(bin_ds_[bin_ds_.Size()]);
					}
				}
			}

		// Memory allocation
		ChangeDimensions(gas.NumberOfSpecies(), &xGas_);
		ChangeDimensions(bin_indices_.Size(), &bin_omega_);
		ChangeDimensions(bin_indices_.Size(), &bin_x_);
		ChangeDimensions(bin_indices_.Size(), &bin_fv_);
		ChangeDimensions(bin_indices_.Size(), &bin_rho_);
		ChangeDimensions(bin_indices_.Size(), &bin_N_);

		if (iBin_ == true)
		{
			for(int i=1;i<=bin_indices_.Size();i++)
			{
				bool iNew = true;
				for(int k=1;k<=bin_baskets_.Size();k++)
					if (bin_c_[i] == bin_baskets_[k])
					{
						iNew = false;
						break;
					}
				if (iNew == true)
					bin_baskets_.Append(bin_c_[i]);
			}

			Sort(&bin_baskets_);
			bin_baskets_indices_ = new BzzVectorInt[bin_baskets_.Size()+1];

			for(int i=1;i<=bin_indices_.Size();i++)
				for(int k=1;k<=bin_baskets_.Size();k++)
				{
					if (bin_c_[i] == bin_baskets_[k])
						bin_baskets_indices_[k].Append(i);
				}

			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_d_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_mw_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_log10d_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_dlog10d_);
			ChangeDimensions(bin_baskets_.Size(), &dN_over_dlog10d_);

			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_V_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_log10V_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_dlog10V_);
			ChangeDimensions(bin_baskets_.Size(), &dN_over_dlog10V_);

			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_m_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_log10m_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_dlog10m_);
			ChangeDimensions(bin_baskets_.Size(), &dN_over_dlog10m_);

			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_N_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_fv_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_rho_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_x_);
			ChangeDimensions(bin_baskets_.Size(), &bin_baskets_omega_);


			for(int k=1;k<=bin_baskets_.Size();k++)
			{
				bin_baskets_mw_[k] = 0.;
				bin_baskets_d_[k] = 0.;
				bin_baskets_m_[k] = 0.;
				bin_baskets_V_[k] = 0.;

				for(int i=1;i<=bin_baskets_indices_[k].Size();i++)
				{	
					int j = bin_baskets_indices_[k][i];
					bin_baskets_d_[k] += bin_d_[j];
					bin_baskets_mw_[k] += bin_mw_[j];
					bin_baskets_V_[k] += bin_V_[j];
					bin_baskets_m_[k] += bin_m_[j];
				}
				
				bin_baskets_d_[k] /= double(bin_baskets_indices_[k].Size());
				bin_baskets_mw_[k] /= double(bin_baskets_indices_[k].Size());
				bin_baskets_V_[k] /= double(bin_baskets_indices_[k].Size());
				bin_baskets_m_[k] /= double(bin_baskets_indices_[k].Size());
				
				bin_baskets_log10d_[k] = log10(bin_baskets_d_[k]);
				bin_baskets_log10V_[k] = log10(bin_baskets_V_[k]);
				bin_baskets_log10m_[k] = log10(bin_baskets_m_[k]);
			}

			// Intervals
			bin_baskets_dlog10d_[1] = ((bin_baskets_log10d_[2]+bin_baskets_log10d_[1])/2. - bin_baskets_log10d_[1])*2.;
			bin_baskets_dlog10V_[1] = ((bin_baskets_log10V_[2]+bin_baskets_log10V_[1])/2. - bin_baskets_log10V_[1])*2.;
			bin_baskets_dlog10m_[1] = ((bin_baskets_log10m_[2]+bin_baskets_log10m_[1])/2. - bin_baskets_log10m_[1])*2.;

			for(int k=2;k<=bin_baskets_.Size()-1;k++)
			{
				bin_baskets_dlog10d_[k] = (bin_baskets_log10d_[k+1]+bin_baskets_log10d_[k])/2. - (bin_baskets_log10d_[k]+bin_baskets_log10d_[k-1])/2.;
				bin_baskets_dlog10V_[k] = (bin_baskets_log10V_[k+1]+bin_baskets_log10V_[k])/2. - (bin_baskets_log10V_[k]+bin_baskets_log10V_[k-1])/2.;
				bin_baskets_dlog10m_[k] = (bin_baskets_log10m_[k+1]+bin_baskets_log10m_[k])/2. - (bin_baskets_log10m_[k]+bin_baskets_log10m_[k-1])/2.;
			}

			bin_baskets_dlog10d_[bin_baskets_.Size()] = ((bin_baskets_log10d_[bin_baskets_.Size()]+bin_baskets_log10d_[bin_baskets_.Size()-1])/2. - bin_baskets_log10d_[bin_baskets_.Size()-1])*2.;
			bin_baskets_dlog10V_[bin_baskets_.Size()] = ((bin_baskets_log10V_[bin_baskets_.Size()]+bin_baskets_log10V_[bin_baskets_.Size()-1])/2. - bin_baskets_log10V_[bin_baskets_.Size()-1])*2.;
			bin_baskets_dlog10m_[bin_baskets_.Size()] = ((bin_baskets_log10m_[bin_baskets_.Size()]+bin_baskets_log10m_[bin_baskets_.Size()-1])/2. - bin_baskets_log10m_[bin_baskets_.Size()-1])*2.;
 

		}
	}

	if (iBin_==true || iSoot_==true)
		WriteSummaryFiles();

}

void OpenSMOKE_PolimiSoot::Analysis(OpenSMOKE_IdealGas &gas, const double P_Pa, const double T, BzzVector &omegaGas)
{
	if (iBin_ == true)
	{
		double small_eps = 1e-20;

		double MWGas;

		gas.GetMWAndMoleFractionsFromMassFractions(MWGas, xGas_, omegaGas);
		double rhoGas = P_Pa*MWGas/(Constants::R_J_kmol*T);

		for(int i=1;i<=bin_indices_.Size();i++)
		{
			int j = bin_indices_[i];
			bin_omega_[i]	= omegaGas[j];							// mass fraction
			bin_x_[i]		= xGas_[j];								// mole fraction
			bin_rho_[i]		= rhoGas*omegaGas[j];					// density [kg_soot/m3]
			bin_fv_[i]		= bin_rho_[i]/bin_density_[i];			// volume fraction [m3_soot/m3]
			bin_N_[i]		= bin_fv_[i]/bin_V_[i];					// [1/m3]
		}

		fv_small_		= 0.;
		rho_small_		= 0.;
		N_small_		= 0.;
		omega_small_	= 0.;
		x_small_		= 0.;
		h_over_c_small_ = 0.;
		o_over_c_small_ = 0.;
		o_over_h_small_ = 0.;

		for(int i=1;i<=bin_indices_small_.Size();i++)
		{
			int j = bin_indices_small_[i];
			fv_small_		+= bin_fv_[j];
			rho_small_		+= bin_rho_[j];
			N_small_		+= bin_N_[j];
			omega_small_	+= bin_omega_[j];
			x_small_		+= bin_x_[j];
			h_over_c_small_ += bin_omega_[j]*bin_h_over_c_[j];
			o_over_c_small_ += bin_omega_[j]*bin_o_over_c_[j];
			o_over_h_small_ += bin_omega_[j]*bin_o_over_h_[j];
		}

		double denominator = omega_small_;
		if ( denominator < small_eps)	denominator = 1e32;
		h_over_c_small_ /= denominator;
		o_over_c_small_ /= denominator;
		o_over_h_small_ /= denominator;

		if (iSoot_ == true)
		{
			fv_large_		= 0.;
			rho_large_		= 0.;
			N_large_		= 0.;
			omega_large_	= 0.;
			x_large_		= 0.;
			h_over_c_large_ = 0.;
			o_over_c_large_ = 0.;
			o_over_h_large_ = 0.;

			for(int i=1;i<=bin_indices_large_.Size();i++)
			{
				int j = bin_indices_large_[i];
				fv_large_		+= bin_fv_[j];
				rho_large_		+= bin_rho_[j];
				N_large_		+= bin_N_[j];
				omega_large_	+= bin_omega_[j];
				x_large_		+= bin_x_[j];
				h_over_c_large_ += bin_omega_[j]*bin_h_over_c_[j];
				o_over_c_large_ += bin_omega_[j]*bin_o_over_c_[j];
				o_over_h_large_ += bin_omega_[j]*bin_o_over_h_[j];
			}
			
			double denominator = omega_large_;
			if ( denominator < small_eps)	denominator = 1e32;
			h_over_c_large_ /= denominator;
			o_over_c_large_ /= denominator;
			o_over_h_large_ /= denominator;
		}

	}
}

void OpenSMOKE_PolimiSoot::Distribution()
{
	if (iBin_ == true)
	{
		for(int k=1;k<=bin_baskets_.Size();k++)
		{
			bin_baskets_N_[k]   = 0.;
			bin_baskets_fv_[k]  = 0.;
			bin_baskets_rho_[k] = 0.;
			bin_baskets_omega_[k] = 0.;
			bin_baskets_x_[k] = 0.;
			for(int i=1;i<=bin_baskets_indices_[k].Size();i++)
			{	
				int j = bin_baskets_indices_[k][i];
				bin_baskets_N_[k]		+= bin_N_[j];
				bin_baskets_fv_[k]		+= bin_fv_[j];
				bin_baskets_rho_[k]		+= bin_rho_[j];
				bin_baskets_omega_[k]	+= bin_omega_[j];
				bin_baskets_x_[k]		+= bin_x_[j];
			}

			dN_over_dlog10d_[k] = bin_baskets_N_[k] / bin_baskets_dlog10d_[k];
			dN_over_dlog10V_[k] = bin_baskets_N_[k] / bin_baskets_dlog10V_[k];
			dN_over_dlog10m_[k] = bin_baskets_N_[k] / bin_baskets_dlog10m_[k];
		}
	}
}

void OpenSMOKE_PolimiSoot::ProcessDistribution()
{
	if (iBin_ == true)
	{
		double small_eps = 1e-16;

		double m0 = bin_baskets_N_.GetSumElements();
		double m1 = Dot(bin_baskets_N_, bin_baskets_d_);
		
		double m2 = 0.;
		double m3 = 0.;
		for(int k=1;k<=bin_baskets_.Size();k++)
		{
			m2 += bin_baskets_N_[k]*bin_baskets_d_[k]*bin_baskets_d_[k];
			m3 += bin_baskets_N_[k]*bin_baskets_d_[k]*bin_baskets_d_[k]*bin_baskets_d_[k];
		}

		if (m0 < small_eps)	m0 = 1e32;
		if (m2 < small_eps)	m2 = 1e32;

		// Mean diameters [m]
		dmean_N_large_	= m1/m0;
		d32_N_large_	= m3/m2;

		// Variance [m2]
		dvariance_N_ = 0.;
		for(int k=1;k<=bin_baskets_.Size();k++)
			dvariance_N_ += bin_baskets_N_[k]*BzzPow2(bin_baskets_d_[k]-dmean_N_large_);
		dvariance_N_	/= m0;

		// Standard deviation [m]
		dstd_N_ = sqrt(dvariance_N_);
	}
}

void OpenSMOKE_PolimiSoot::WriteSummaryFiles()
{
	// Write Bin Summary
	{
		ofstream fOutput;
		openOutputFileAndControl(fOutput, "BinSummary.out");
		fOutput.setf(ios::scientific);
		fOutput << setw(5)  << left << "#"; 
		fOutput << setw(5)  << left << "#";
		fOutput << setw(10) << left << "Name";	
		fOutput << setw(16) << fixed << left << setprecision(2) << "MW[kg/kmol]";	// [kg/kmol]
		fOutput << setw(12) << fixed << left << setprecision(5) << "d[nm]";			// [nm]
		fOutput << setw(12) << fixed << left << setprecision(5) << "dsph[nm]";		// [nm]
		fOutput << setw(12) << fixed << left << setprecision(5) << "dcol[nm]";		// [nm]
		fOutput << setw(16) << scientific << left << "m[mug]";						// [mug] 
		fOutput << setw(16) << scientific << left << "V[cm3]";						// [cm3]
		fOutput << setw(12) << scientific << left << "np[-]";						// [-]
		fOutput << setw(10) << fixed << left << setprecision(0) << "C";
		fOutput << setw(10) << fixed << left << setprecision(0) << "H"; 
		fOutput << setw(10) << fixed << left << setprecision(0) << "O";
		fOutput << setw(10) << fixed << left << setprecision(4) << "H/C"; 
		fOutput << setw(10) << fixed << left << setprecision(4) << "O/C"; 
		fOutput << setw(10) << fixed << left << setprecision(4) << "O/H";
		fOutput << endl;

		for(int i=1;i<=bin_indices_small_.Size();i++)
		{
			int j = bin_indices_small_[i];
			fOutput << setw(5)  << left << j; 
			fOutput << setw(5)  << left << i; 
			fOutput << setw(10) << left << bin_names_[j].c_str();	
			fOutput << setw(16) << fixed << left << setprecision(2) << bin_mw_[j];			// [kg/kmol]
			fOutput << setw(12) << fixed << left << setprecision(5) << bin_d_[j] * 1e9;		// [nm]
			fOutput << setw(12) << fixed << left << setprecision(5) << bin_ds_[j] * 1e9;	// [nm]
			fOutput << setw(12) << fixed << left << setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
			fOutput << setw(16) << scientific << left << bin_m_[j] * 1e9;					// [mug] 
			fOutput << setw(16) << scientific << left << bin_V_[j] * 1e6;					// [cm3]
			fOutput << setw(12) << fixed << left << setprecision(2) << bin_np_[j];			// [-]
			fOutput << setw(10) << fixed << left << setprecision(0) << bin_c_[j];
			fOutput << setw(10) << fixed << left << setprecision(0) << bin_h_[j]; 
			fOutput << setw(10) << fixed << left << setprecision(0) << bin_o_[j];
			fOutput << setw(10) << fixed << left << setprecision(4) << bin_h_over_c_[j]; 
			fOutput << setw(10) << fixed << left << setprecision(4) << bin_o_over_c_[j]; 
			fOutput << setw(10) << fixed << left << setprecision(4) << bin_o_over_h_[j];
			fOutput << endl;
		}
		fOutput << endl;
	

		if (iSoot_ == true)
		{
			for(int i=1;i<=bin_indices_large_.Size();i++)
			{
				int j = bin_indices_large_[i];
				fOutput << setw(5)  << left << j; 
				fOutput << setw(5)  << left << i; 
				fOutput << setw(10) << left << bin_names_[j].c_str();	
				fOutput << setw(16) << fixed << left << setprecision(2) << bin_mw_[j];			// [kg/kmol]
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_d_[j] * 1e9;		// [nm]
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_ds_[j] * 1e9;	// [nm]
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
				fOutput << setw(16) << scientific << left << bin_m_[j]*1e9;						// [mug] 
				fOutput << setw(16) << scientific << left << bin_V_[j]*1e6;						// [cm3]
				fOutput << setw(12) << fixed << left << setprecision(2) << bin_np_[j];			// [-]
				fOutput << setw(10) << fixed << left << setprecision(0) << bin_c_[j]; 
				fOutput << setw(10) << fixed << left << setprecision(0) << bin_h_[j]; 
				fOutput << setw(10) << fixed << left << setprecision(0) << bin_o_[j];
				fOutput << setw(10) << fixed << left << setprecision(4) << bin_h_over_c_[j]; 
				fOutput << setw(10) << fixed << left << setprecision(4) << bin_o_over_c_[j]; 
				fOutput << setw(10) << fixed << left << setprecision(4) << bin_o_over_h_[j];
				fOutput << endl;
			}
		}

		fOutput.close();
	}
	
	// Write Soot Summary
	{
		ofstream fOutput;
		openOutputFileAndControl(fOutput, "SootSummary.out");
		fOutput.setf(ios::scientific);

		fOutput << setw(5)  << left << "#"; 
		fOutput << setw(5)  << left << "#";
		fOutput << setw(10) << left << "Name";	
		fOutput << setw(16) << fixed << left << setprecision(2) << "MW[kg/kmol]";	// [kg/kmol]
		fOutput << setw(12) << fixed << left << setprecision(5) << "d[nm]";			// [nm]
		fOutput << setw(12) << fixed << left << setprecision(5) << "dsph[nm]";		// [nm]
		fOutput << setw(12) << fixed << left << setprecision(5) << "dcol[nm]";		// [nm]
		fOutput << setw(16) << scientific << left << "m[mug]";						// [mug] 
		fOutput << setw(16) << scientific << left << "V[cm3]";						// [cm3]
		fOutput << setw(12) << scientific << left << "np[-]";						// [-]
		fOutput << setw(10) << fixed << left << setprecision(0) << "C"; 
		fOutput << setw(10) << fixed << left << setprecision(0) << "H"; 
		fOutput << setw(10) << fixed << left << setprecision(0) << "O";
		fOutput << setw(10) << fixed << left << setprecision(4) << "H/C"; 
		fOutput << setw(10) << fixed << left << setprecision(4) << "O/C"; 
		fOutput << setw(10) << fixed << left << setprecision(4) << "O/H";
		fOutput << endl;
		
		for(int k=1;k<=bin_baskets_.Size();k++)
		{
			for(int i=1;i<=bin_baskets_indices_[k].Size();i++)
			{	
				int j = bin_baskets_indices_[k][i];
				fOutput << setw(5)  << left << j; 
				fOutput << setw(5)  << left << i; 
				fOutput << setw(10) << left << bin_names_[j].c_str();	
				fOutput << setw(16) << fixed << left << setprecision(2) << bin_mw_[j];			// [kg/kmol]
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_d_[j]*1e9;		// [nm]
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_ds_[j]*1e9;		// [nm]
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
				fOutput << setw(16) << scientific << left << bin_m_[j]*1e9;						// [mug] 
				fOutput << setw(16) << scientific << left << bin_V_[j]*1e6;						// [cm3]
				fOutput << setw(12) << fixed << left << setprecision(2) << bin_np_[j];			// [-]
				fOutput << setw(10) << fixed << left << setprecision(0) << bin_c_[j]; 
				fOutput << setw(10) << fixed << left << setprecision(0) << bin_h_[j]; 
				fOutput << setw(10) << fixed << left << setprecision(0) << bin_o_[j];
				fOutput << setw(10) << fixed << left << setprecision(4) << bin_h_over_c_[j]; 
				fOutput << setw(10) << fixed << left << setprecision(4) << bin_o_over_c_[j]; 
				fOutput << setw(10) << fixed << left << setprecision(4) << bin_o_over_h_[j];
				fOutput << endl;
			}

			fOutput << endl;
		}
		fOutput << endl;
		fOutput << endl;

		fOutput << setw(5)  << left << "#"; 
		fOutput << setw(10) << fixed << left << setprecision(0) << "C"; 
		fOutput << setw(12) << fixed << left << setprecision(5) << "d[nm]";			
		fOutput << setw(12) << fixed << left << setprecision(5) << "log10d[m]";		
		fOutput << setw(12) << fixed << left << setprecision(5) << "dlog10d[m]";	
		fOutput << setw(12) << fixed << left << setprecision(5) << "m[kg*1e21]";			
		fOutput << setw(12) << fixed << left << setprecision(5) << "log10m[kg]";	
		fOutput << setw(12) << fixed << left << setprecision(5) << "dlog10m[kg]";	
		fOutput << setw(12) << fixed << left << setprecision(5) << "V[m3*1e26]";			
		fOutput << setw(12) << fixed << left << setprecision(5) << "log10V[m3]";	
		fOutput << setw(12) << fixed << left << setprecision(5) << "dlog10V[m3]";	

		fOutput << endl;

		for(int k=1;k<=bin_baskets_.Size();k++)
			{
				cout << k << endl;
				fOutput << setw(5)  << left << k; 
				fOutput << setw(10) << fixed << left << setprecision(0) << bin_baskets_[k]; 

				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_d_[k]*1e9;		
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_log10d_[k];		
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_dlog10d_[k];		

				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_m_[k]*1e21;		
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_log10m_[k];		
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_dlog10m_[k];		

				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_V_[k]*1e26;		
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_log10V_[k];		
				fOutput << setw(12) << fixed << left << setprecision(5) << bin_baskets_dlog10V_[k];		

				fOutput << endl;
			}

		fOutput.close();
	}

}
