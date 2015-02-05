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

#ifndef OPENSMOKE_POLIMISOOT_H
#define OPENSMOKE_POLIMISOOT_H

#include "BzzMath.hpp"
#include "engine/OpenSMOKE_IdealGas.h"

class OpenSMOKE_PolimiSoot
{
public:

	void Setup(OpenSMOKE_IdealGas &gas, const std::string minimum_bin,
							const unsigned int bin_index_zero, const double bin_density_A, 
							const unsigned int bin_index_final, const double bin_density_B,
							const double Df);

	void Setup(OpenSMOKE_IdealGas &gas, const std::string minimum_bin);

	void Analysis(OpenSMOKE_IdealGas &gas, const double P_Pa, const double T, BzzVector &omegaGas);
	void Distribution();
	void ProcessDistribution();

	bool IsBin()	{ return iBin_; }
	bool IsSoot()	{ return iSoot_; }

	BzzVectorInt& bin_indices()		{ return bin_indices_; }

	double fv_small()		{return fv_small_; }
	double rho_small()		{return rho_small_; }
	double N_small()		{return N_small_; }
	double omega_small()	{return omega_small_; }
	double x_small()		{return x_small_; }
	
	double fv_large()		{return fv_large_; }
	double rho_large()		{return rho_large_; }
	double N_large()		{return N_large_; }
	double omega_large()	{return omega_large_; }
	double x_large()		{return x_large_; }

	BzzVector& baskets_d()					{ return bin_baskets_d_; }
	BzzVector& baskets_m()					{ return bin_baskets_m_; }
	BzzVector& baskets_V()					{ return bin_baskets_V_; }
	BzzVector& baskets_x()					{ return bin_baskets_x_; }
	BzzVector& baskets_omega()				{ return bin_baskets_omega_; }
	BzzVector& baskets_N()					{ return bin_baskets_N_; }
	BzzVector& baskets_fv()					{ return bin_baskets_fv_; }
	BzzVector& baskets_rho()				{ return bin_baskets_rho_; }
	BzzVector& baskets_dN_over_dlog10d()	{ return dN_over_dlog10d_; }
	BzzVector& baskets_dN_over_dlog10m()	{ return dN_over_dlog10m_; }
	BzzVector& baskets_dN_over_dlog10V()	{ return dN_over_dlog10V_; }
	BzzVector& baskets_dlog10d()			{ return bin_baskets_dlog10d_; }
	BzzVector& baskets_dlog10m()			{ return bin_baskets_dlog10m_; }
	BzzVector& baskets_dlog10V()			{ return bin_baskets_dlog10V_; }

	

	double h_over_c_small()		{ return h_over_c_small_; }
	double o_over_c_small()		{ return o_over_c_small_; }
	double o_over_h_small()		{ return o_over_h_small_; }
	double h_over_c_large()		{ return h_over_c_large_; }
	double o_over_c_large()		{ return o_over_c_large_; }
	double o_over_h_large()		{ return o_over_h_large_; }
	double dmean_N_small()		{ return dmean_N_small_; }
	double dmean_N_large()		{ return dmean_N_large_; }
	double dvariance_N()		{ return dvariance_N_; }
	double dstd_N()				{ return dstd_N_; }
	double d32_N_small()		{ return d32_N_small_; }
	double d32_N_large()		{ return d32_N_large_; }

	void SetIndexZero(const unsigned int bin_index_zero)
	{
		bin_index_zero_ = bin_index_zero;
	}

	void SetIndexFinal(const unsigned int bin_index_final)
	{
		bin_index_final_ = bin_index_final;
	}

	void SetBinDensityA(const double bin_density_A)
	{
		 bin_density_A_ =  bin_density_A;
	}

	void SetBinDensityB(const double bin_density_B)
	{
		 bin_density_B_ =  bin_density_B;
	}

private:

	BzzVector       bin_density_;
	BzzVectorInt	bin_indices_;
	BzzVectorInt    bin_indices_small_;
	BzzVectorInt    bin_indices_large_;
	BzzVector		bin_mw_;
	BzzVector		bin_m_;
	BzzVector		bin_d_;
	BzzVector		bin_V_;
	BzzVector		bin_c_;
	BzzVector		bin_h_;
	BzzVector		bin_o_;
	BzzVector		bin_h_over_c_;
	BzzVector		bin_o_over_c_;
	BzzVector		bin_o_over_h_;
	vector<string>	bin_names_;

	BzzVector		bin_dc_;
	BzzVector		bin_ds_;
	BzzVector		bin_np_;

	BzzVector		xGas_;
	BzzVector		bin_omega_;
	BzzVector		bin_x_;
	BzzVector		bin_fv_;
	BzzVector		bin_rho_;
	BzzVector		bin_N_;

	BzzVector     bin_baskets_;
	BzzVectorInt* bin_baskets_indices_;
	
	BzzVector     bin_baskets_d_;
	BzzVector     bin_baskets_mw_;
	BzzVector     bin_baskets_log10d_;
	BzzVector     bin_baskets_dlog10d_;
	BzzVector     dN_over_dlog10d_;

	BzzVector     bin_baskets_V_;
	BzzVector     bin_baskets_log10V_;
	BzzVector     bin_baskets_dlog10V_;
	BzzVector     dN_over_dlog10V_;

	BzzVector     bin_baskets_m_;
	BzzVector     bin_baskets_log10m_;
	BzzVector     bin_baskets_dlog10m_;
	BzzVector     dN_over_dlog10m_;

	BzzVector     bin_baskets_N_;
	BzzVector     bin_baskets_fv_;
	BzzVector     bin_baskets_rho_;
	BzzVector     bin_baskets_x_;
	BzzVector     bin_baskets_omega_;
	
	// Index to check if the soot bins are available in the current kinetic scheme
	bool iSoot_;
	bool iBin_;

	double fv_small_;
	double rho_small_;
	double N_small_;
	double omega_small_;
	double x_small_;
	double fv_large_;
	double rho_large_;
	double N_large_;
	double omega_large_;
	double x_large_;
	double h_over_c_small_;
	double o_over_c_small_;
	double o_over_h_small_;
	double h_over_c_large_;
	double o_over_c_large_;
	double o_over_h_large_;
	double dmean_N_small_;
	double dmean_N_large_;
	double dvariance_N_;
	double dstd_N_;
	double d32_N_small_;
	double d32_N_large_;
	double Df_;

	void WriteSummaryFiles();

	unsigned int bin_index_zero_;
	unsigned int bin_index_final_;
	double bin_density_A_;
	double bin_density_B_;
};


#endif // OPENSMOKE_POLIMISOOT_H
