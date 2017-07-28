/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
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
#if !defined(OpenSMOKE_Flame1D_OscillatingBoundary_H)
#define OpenSMOKE_Flame1D_OscillatingBoundary_H

#include "OpenSMOKE.hpp"
#include "OpenSMOKE_Flame1D_DataManager.h"

enum unsteady_boundary_kinds {	OSCILLATING_BOUNDARY_NONE, OSCILLATING_BOUNDARY_SIN, 
				OSCILLATING_BOUNDARY_AIR_VELOCITY, OSCILLATING_BOUNDARY_FUEL_AIR_VELOCITIES, 
				OSCILLATING_BOUNDARY_AIR_TEMPERATURE, OSCILLATING_BOUNDARY_FUEL_COMPOSITION, OSCILLATING_BOUNDARY_SIN_EXP, 
				RELAXATION_BOUNDARY_FUEL_VELOCITY, RELAXATION_BOUNDARY_FUEL_TEMPERATURE, RELAXATION_BOUNDARY_FUEL_EQRATIO,
				DYNAMIC_TEMPERATURE, DYNAMIC_VELOCITY };

class OpenSMOKE_Flame1D_OscillatingBoundary
{
public:

	double tDelay;
	unsteady_boundary_kinds	unsteady_boundary_kind;
	
	int iShape;
	int nPeriods;
	double T;
	double phase;
	double time;
	
	int onPrint;
	int iPrint;
	int lastPrint;
	
	double	K;
	double  KSeshadri;
	double	Kold;
	double	slope;

	// Steady state values
	double USteadyFuel;
	double USteadyAir;
	double rhoSteadyFuel;
	double rhoSteadyAir;
	double MWSteadyFuel; 
	double MWSteadyAir;
	BzzVector XSteadyFuel;
	BzzVector XSteadyAir;
	BzzVector YSteadyFuel;
	BzzVector YSteadyAir;
	double TSteadyFuel; 
	double TSteadyAir;
	double vSteadyFuel; 
	double vSteadyAir;
	double phiSteadyFuel; 
	double phiSteadyAir;
	double KSteady;
	double P; 
	double L; 
	double xSt;
	double betaCoefficient;

	// Relaxation on the fuel side
	double vTargetFuel;
	double TTargetFuel;
	double phiTargetFuel;
	double relaxationTimeFuel;
	double relaxationFuelCoefficient;

	// Dynamics
	double slope_velocity_fuel_;
	double slope_velocity_oxidizer_;
	double slope_temperature_fuel_;
	double slope_temperature_oxidizer_;

	// Equivalent Strain rate
	int Ntimes;
	BzzVector Kvector;
	BzzVector timeVector;
	double Kequivalent;
	void update_K_equivalent(double _time);
	double integral_a(int j);
	double integral_b();	
	
	// Oscillation data
	double frequency;
	double semiAmplitudeFuel;
	double semiAmplitudeAir;
	double _2pi;
	double	InPhase;

	void setup(double _USteadyFuel, double _USteadyAir,   double _rhoSteadyFuel, double _rhoSteadyAir, 
			   double _MWSteadyFuel, double _MWSteadyAir, BzzVector& _YSteadyFuel, BzzVector& _YSteadyAir, double _TSteadyFuel, double _TSteadyAir,
			   double _L, double _P, OpenSMOKE_ReactingGas  *_mix, const std::string fileName);

	void update_boundary_conditions(double _time, 	double &UC,    double &UO, 
							double &TFuel, double &TAir,
							double &rhoC,  double &rhoO, double _Tmax,
							BzzVector &WC, BzzVector &WO,
							OpenSMOKE_Flame1D_DataManager& data);
	
	void update_time_target(const double time);

private:

	double give_oscillating_profile_Fuel(double time);
	double give_oscillating_profile_Air(double time);
	double give_oscillating_profile_Fuel(double time, double incresing_factor);
	double give_oscillating_profile_Air(double time, double incresing_factor);

	OpenSMOKE_ReactingGas  *mix;

	vector<double> list_target_fuel_velocities;
	vector<double> list_target_fuel_temperatures;
	vector<double> list_target_fuel_equivalence_ratios;
	int next_target;

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};

#endif
