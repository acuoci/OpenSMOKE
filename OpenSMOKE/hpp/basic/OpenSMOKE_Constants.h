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

#ifndef OPENSMOKE_CONSTANTS
#define OPENSMOKE_CONSTANTS

namespace Constants
{
    const double pi          		= 3.1415926535897932384626433;  // pi
    const double _2pi          		= 2.*pi;						// 2pi
    const double pi_over_2			= pi/2.;						// pi/2

									  
    const double R_J_mol     		= 8.3144621;			// Ideal gas constant [J/mol/K]
    const double R_J_kmol    		= 8314.4621;			// Ideal gas constant [J/kmol/K]
	const double R_cal_mol			= 8.3144621/4.18443;	// Ideal gas constant [cal/mol/K]
    const double T_Zero      		= 273.15;           // Freezing point [K]
    const double m_u         		= 1.660538782e-27;  // atomic mass constant [kg]
    const double Nav_mol         	= 6.0221417930e23;  // Avogadro constant [1/mol]
    const double Nav_kmol      		= 6.0221417930e26;  // Avogadro constant [1/kmol]
    const double kBoltzmann  		= 1.3806504e-23;    // Boltzmann constant [J/K]
    const double hPlanck     		= 6.62606896e-34;	// Planck constant [J.s]
    const double vLight      		= 299792458.0;		// Speed of light in vacuum [m/s]
    const double sigma       		= 5.670400e-8;		// Stefan-Boltzmann constant [W/m2/K4]
    const double g           		= 9.80665;			// gravitational acceleration [m/s2]
    const double T_Reference        = 298.10;			// std reference temperature [K]
    const double P_Reference        = 101325.;			// std reference pressure [Pa]
	const double SmallMoleFraction	= 1.e-32;			// minimum mole fraction [-]
	const double Msoot				= 12.010999;		// carbon molecular weight [kg/kmol]
	const double densitySoot		= 1800.00;			// soot density [kg/m3]

	const double Conversion_cal_J	= R_J_kmol/R_cal_mol;	// [J/cal]

	const double Kg_min				= 1.e-3;
	const double Kg_max				= 0.98;

	const double csi_min			= 1.e-3;
	const double csi_max			= 0.995;

	const double TMinExtinction		= 1000.;			// [K]

	const int	 COMMENT_SIZE		= 300;				// size of comments
	const int	 NAME_SIZE			= 40;				// size of names
	const int	 REACTION_NAME_SIZE	= 400;				// size of reaction names

	const double MaxMWThermalDiffusionRatios = 20.;		// kg/kmol
}

enum ProfileKind				 { NONE, CONSTANT, USERDEFINED };
enum kindOfSensitivityParameter  { FREQUENCY_FACTOR, TRANSPORT_PROPERTIES, FREQUENCY_FACTOR_AND_TRANSPORT_PROPERTIES };
enum SymbolicKinetics			 { NOJACOBIAN, POLIMI_C1C3HTNOX_0810, POLIMI_C1C3HTNOX_AVIO, GRI30, SANDIEGO_AVIO, FLUENT_GLARBORG152,
								   FLUENT_DRM22_POLIMI, FLUENT_DRM22_POLIMI_NOX, FLUENT_DRM22_POLIMI_THERMALNOX, GRI12 };
enum SolidReactorKind			 { NOSOLIDREACTOR, CONSTP, CONSTV, TG_ISO, TG_NONISO};
enum SolidParticleModel			 { NOSOLIDMODEL,	NORESISTANCE, SHRINKING, MODEL1D};
enum kindOfPressureDependence	 { PRESSURE_LINDEMANN, PRESSURE_TROE, PRESSURE_SRI};
enum versionKinetics			 { V090524, V090905, V101116};
enum equilibriumProblem			 { EQUILIBRIUM_TP, EQUILIBRIUM_HP, EQUILIBRIUM_SP, EQUILIBRIUM_TV, EQUILIBRIUM_UV, EQUILIBRIUM_SV, EQUILIBRIUM_TRHO, EQUILIBRIUM_URHO, EQUILIBRIUM_SRHO};
enum solid_regression_parameters {	SOLID_REGRESSION_ALL, SOLID_REGRESSION_A, SOLID_REGRESSION_BETA, SOLID_REGRESSION_E,
									SOLID_REGRESSION_A_E, SOLID_REGRESSION_A_BETA, SOLID_REGRESSION_BETA_E, SOLID_REGRESSION_USERDEFINED,
									SOLID_REGRESSION_ALL_PHI, SOLID_REGRESSION_A_E_PHI };
enum adaptive_grid_model		 { ADAPTIVE_GRADIENT, ADAPTIVE_CHEMKIN, ADAPTIVE_CHEMKIN_SQRT, ADAPTIVE_CHEMKIN_POW, ADAPTIVE_EISEMAN};

enum flame1d_model				 {	OPPOSED_ALL, OPPOSED_ONLY_MOMENTUM, OPPOSED_ONLY_TEMPERATURE, OPPOSED_ONLY_MASS_FRACTIONS,
									OPPOSED_NO_ENERGY, OPPOSED_NO_MOMENTUM, OPPOSED_COLD_REDUCED, 
									OPPOSED_SOOT_ALL, OPPOSED_QMOM_ALL,
									PREMIXED_ALL, PREMIXED_NOENERGY, PREMIXED_FLAMESPEED,
									PREMIXED_QMOM_ALL, PREMIXED_QMOM_NOENERGY, PREMIXED_SOOT_ALL, PREMIXED_SOOT_NOENERGY,
									TWIN_COLD_REDUCED, TWIN_ONLY_MASS_FRACTIONS, TWIN_ONLY_MOMENTUM, TWIN_ONLY_TEMPERATURE, TWIN_ALL, TWIN_NO_ENERGY,
									TWIN_NO_MOMENTUM, TWIN_SOOT_ALL, TWIN_QMOM_ALL,
									SINGLEREACTOR_ISOTHERMAL };

enum flame1d_physics			 {	FLAME1D_PHYSICS_PREMIXED, FLAME1D_PHYSICS_OPPOSED, FLAME1D_PHYSICS_TWIN};
enum flame1d_subphysics			 {	FLAME1D_SUBPHYSICS_NONE, FLAME1D_SUBPHYSICS_SOOT, FLAME1D_SUBPHYSICS_QMOM};

enum nucleation_models 			{	NUCLEATION_NONE, NUCLEATION_LIU_2001, NUCLEATION_LIU_2002, NUCLEATION_LIU_2003, NUCLEATION_MOSS_1999,
									NUCLEATION_WEN_2003, NUCLEATION_LINDSTEDT_1994, NUCLEATION_LEUNG_1991, NUCLEATION_HALL_1997};
enum growth_models 				{	GROWTH_NONE, GROWTH_LIU_2001, GROWTH_LIU_2002, GROWTH_LIU_2003, GROWTH_MOSS_1999,
									GROWTH_WEN_2003, GROWTH_LINDSTEDT_1994, GROWTH_LEUNG_1991};
enum aggregation_models 		{	AGGREGATION_NONE, AGGREGATION_SMOLUCHOWSKI, AGGREGATION_MOSS};												
enum oxidation_models 			{	OXIDATION_NONE, OXIDATION_LEE, OXIDATION_NSC, OXIDATION_NEOH};

enum interactions_modes 		{	INTERACTIONS_NONE, INTERACTIONS_UNCORRELATED, INTERACTIONS_CORRELATED, INTERACTIONS_BRUTAL };

enum surfaceReactionTags { SURFACE_REACTION_CONVENTIONAL, 
						   SURFACE_REACTION_IRREVERSIBLE, SURFACE_REACTION_REVERSIBLE,
						   SURFACE_REACTION_STICKY, SURFACE_REACTION_LH };


namespace liquid_specificheat_equation
{
	enum  equation { EQ1, EQ2 };
}

namespace liquid_density_equation 
{
	enum  equation {EQ0, EQ1, EQ7};
}

namespace liquid_thermalconductivity_equation
{	
	enum  equation { EQ0, EQ1, EQ2 };
}

namespace liquid_vaporizationheat_equation 
{
	enum  equation { EQ1 };
}

namespace liquid_vaporpressure_equation 
{
	enum  equation { EQ1 };
}
#endif // OPENSMOKE_CONSTANTS