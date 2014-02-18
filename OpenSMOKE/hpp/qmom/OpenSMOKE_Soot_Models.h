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

#ifndef		QMOM_SOOTMODULE
#define		QMOM_SOOTMODULE

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"
#include "engine/OpenSMOKE_ReactingGas.h"

class OpenSMOKE_Soot_Models;

class MyNonLinearSystem_Fractal_Dimension : public BzzMyNonLinearSystemObject
{
public:
	void assignSoot(OpenSMOKE_Soot_Models *soot);

	OpenSMOKE_Soot_Models *ptSoot;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Soot_Models  
{
public:

	// Function to initialize this object; must be called just once
	void initialize(const int _N, const int _iFractalDimension, const int _iNucleation,
		            const int _iGrowth, const int _iOxidation, const int _iAggregation,
					const double _L0, const int _iTfluctuations);

	// These two function must be called every for every computational cell
	void update(BzzVector &_w, BzzVector &_L);
	void assign_physical_properties(double _T,  double _P, double _rho, 
									double _mu, double _tr,
									double _omega_O2, double _omega_C2H2,
									double _omega_OH, double _STNV);
	
	// Public functions
	void nucleation_rate(double &J, double &epsilon);
	void oxydation_rate(double &alfa, double &exponent);
	void growth_rate(double &alfa, double &exponent);
	void aggregation_kernel(BzzMatrix &Kernel);
	
	// This function is not really public; it is used just by the non linear solver
	void nls_fractal_dimension(BzzVector &x, BzzVector &f);

	void assign_mixture(OpenSMOKE_ReactingGas &_mix);
	void formation_rates();
	void write_on_file(ofstream &fOutput, double csi, BzzVector &moments);
	void write_label_on_file(ofstream &fOutput, int N);

	BzzVector SGas;
	double S;
	
	// They are public just to easily plot them in FLUENT
	BzzVector Dp;
	double Df;
	double Tau;
	double m0, m2;
	double fv, As, Lmean;

	void read_correction_coefficients_from_file(ifstream &fInput)
	{
		char comment[Constants::COMMENT_SIZE];

		fInput >> CC_Nucleation;	fInput.getline(comment, Constants::COMMENT_SIZE);
		fInput >> CC_Growth;		fInput.getline(comment, Constants::COMMENT_SIZE);
		fInput >> CC_Oxidation;		fInput.getline(comment, Constants::COMMENT_SIZE);
		fInput >> CC_Aggregation;	fInput.getline(comment, Constants::COMMENT_SIZE);
	}

private:
	int N;

	BzzVector w;
	BzzVector L;

	double pi;
	double pi_over_6;
	double _8_over_pi;
	double _2pi;
	double _3pi;
	double _pi_plus_8;
	double pi_squared;
	double _36pi_to_1_over_3;
	
	// Soot properties
	double rhoSoot;
	double PMSoot;
	double L0Soot;
	double L0Soot_3;

	// Gas phase properties
	double T_gas;
	double P_gas;
	double PM_gas;
	double mu_gas;
	double tr_gas;
	double rho_gas;
	double cTot_gas;
	double omega_O2;
	double omega_C2H2;
	double omega_OH;
	double x_O2;
	double x_C2H2;
	double x_OH;
	double c_O2;
	double c_C2H2;
	double c_OH;
	double p_O2;
	double p_C2H2;
	double p_OH;
	double STNV;



	// Aggregation Kernel
	BzzMatrix Beta;
	
	BzzVector Diff;
	BzzVector c;
	BzzVector g;

	// Fractal Dimension
	double Dfmin, Dfmax, Df0;
	double tc;
	double sFractal;
	double Lmean_Volume;
	double Beta_mean_aggregation_kernel;
	MyNonLinearSystem_Fractal_Dimension nls;
	BzzVector xFirstGuess;
	BzzVector fnls;


	void allocate_memory();
	void define_constants();
	void user_defined_variables();	
	void solve_for_fractal_dimension();	
	void coefficient_for_aggregation_kernel(const double Dp, double &Diff, double &c, double &g);
	double collision_diameter(double Lenght);
	void evaluate_fractal_dimension();
	void mean_aggregation_kernel();
	double evaluate_aggregation_kernel(const double Dp1,    const double Dp2, 
									   const double Diff1,	const double Diff2,
									   const double c1,	    const double c2,
									   const double g1,	    const double g2);

	int iNucleation;
	int iGrowth;
	int iOxidation;
	int iAggregation;
	int iFractalDimension;
	int iTfluctuations;

	
	double SNGas_C2H2;
	double SNGas_H2;
	double SGGas_C2H2;
	double SGGas_H2;
	double SOGas_O2;
	double SOGas_OH;
	double SOGas_CO;
	double SOGas_H;;

	int iC2H2;
	int iO2;
	int iOH;
	int iH2;
	int iH;
	int iCO;

	double PM_C2H2;
	double PM_O2;
	double PM_OH;
	double PM_CO;
	double PM_H2;
	double PM_H;

	double omegaSoot;
	double xSoot;
	double cSoot;
	double MSoot;		// soot density [kg/m3]
	double dp;			// d10

	double S_C2H2;
	double S_H2;
	double S_O2;
	double S_OH;
	double S_CO;
	double S_H;;
	
	double Sn;
	double SG;
	double SN;
	double SO_O2;
	double SO_OH;

	double CC_Nucleation;
	double CC_Oxidation;
	double CC_Growth;
	double CC_Aggregation;
};

#endif // !defined(AFX_SOOT_MODELS_H__CA141DA3_02E5_45CC_A7F7_33E9A7BFFDC8__INCLUDED_)
