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

#ifndef OPENSMOKE_CONVERSIONS
#define OPENSMOKE_CONVERSIONS

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"
#include <string>
using namespace std;

namespace OpenSMOKE_Conversions
{
    const double Pa_from_atm 		= 101325.0;
    const double Pa_from_bar  		= 100000.0;
    const double Pa_from_mbar 		= 100.0;
    const double Pa_from_torr 		= 133.32237;
    const double Pa_from_kPa  		= 1000.;
    const double Pa_from_psi  		= 6894.75728;

    const double m_from_mm 		    = 0.0010;
    const double m_from_cm 		    = 0.0100;
    const double m_from_in 		    = 0.0254;
    const double m_from_ft 		    = 0.3048;
    const double m_from_km 		    = 0.3048;

    const double m2_from_mm2 		= 1.e-6;
    const double m2_from_cm2 		= 1.e-4;
    const double m2_from_in2 		= 6.4516e-4;
    const double m2_from_ft2 		= 0.09290304;

    const double m3_from_mm3 		    = 1.e-9;
    const double m3_from_cm3 		    = 1.e-6;
    const double m3_from_l 		        = 1.e-3;
    const double m3_from_in3 		    = 1.6387064e-5;
    const double m3_from_ft3 		    = 0.028316846592;
    const double m3_from_gallon_UK 	    = 0.00454609;
    const double m3_from_oz_UK 		    = 0.0000284130625;
    const double m3_from_oz_USA 	    = 0.000029573529563;
    const double m3_from_gallon_dry_USA = 0.0044048838;
    const double m3_from_gallon_liq_USA = 0.003785411784;

    const double kg_from_g 		= 1.e-3;
    const double kg_from_lb 	= 0.454;
    const double kg_from_oz 	= 0.028349;

    const double s_from_hr 		= 3600.;
    const double s_from_min 	= 60.00;
    const double s_from_ms		= 1e-3;

    const double J_from_kJ 		= 1.e3;
    const double J_from_erg 	= 1.e-7;
    const double J_from_BTU		= 1055.0559;
    const double J_from_kWh		= 3600000.;
    const double J_from_cal		= Constants::R_J_mol/Constants::R_cal_mol;
    const double J_from_kcal	= 1.e3*J_from_cal;
    const double J_from_eV		= 1.6021765314e-19;

    const double kmol_from_mol	= 1.e-3;
    
    const double W_from_kW 		= 1.e3;

    double  conversion_length(const double value, const std::string unit);
    double  conversion_u_length(const double value, const std::string units);
    double  conversion_area(const double value, const std::string unit);
    double  conversion_volume(const double value, const std::string unit);
    double  conversion_specificVolume(const double value, const std::string unit);
    double  conversion_pressure(const double value, const std::string unit);
    double  conversion_time(const double value, const std::string unit);
    double  conversion_energy(const double value, const std::string unit);
    double  conversion_entropy(const double value, const std::string unit);
    double  conversion_specificEnergy(const double value, const std::string unit);
    double  conversion_specificEntropy(const double value, const std::string unit);
	double  conversion_mass(const double value, const std::string unit);
	double  conversion_temperature(const double value, const std::string unit);
    double  conversion_frequency(const double value, const std::string unit);
    double  conversion_velocity(const double value, const std::string unit);
    double  conversion_massFlowRate(const double value, const std::string unit);
    double  conversion_moleFlowRate(const double value, const std::string unit);
    double  conversion_volumetricFlowRate(const double value, const std::string unit);
    double  conversion_heat_flux(const double value, const std::string unit);
    double  conversion_heat_exchange_coefficient(const double value, const std::string unit);
	double  conversion_dynamic_viscosity(const double value, const std::string units);
	double  conversion_density(const double value, const std::string units);
	double  conversion_angle(const double value, const std::string units);
	double  conversion_angular_velocity(const double value, const std::string units);
	double  conversion_area_velocity(const double value, const std::string units);
	double  conversion_valve_flow_coefficient(const double value, const std::string units);
    double  conversion_specificEnergyMolar(const double value, const std::string units);

	void    ErrorMessage(const std::string message);
};

#endif // OPENSMOKE_CONVERSIONS
