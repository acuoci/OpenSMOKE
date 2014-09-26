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

#ifndef OPENSMOKE_ICEM_MULTIZONE_H
#define OPENSMOKE_ICEM_MULTIZONE_H

#include "OpenSMOKE.hpp"
#include "idealreactors/OpenSMOKE_0DReactor.h"

class  OpenSMOKE_ICEM_MultiZone;

enum icem_multizone_models { ICEM_KOMNINOS_SPLITTING, ICEM_CHEMKIN_MULTIZONE };
enum icem_heat_transfer_models {ICEM_HEAT_TRANSFER_NONE, ICEM_HEAT_TRANSFER_WOSCHNI, ICEM_HEAT_TRANSFER_ASSANIS, ICEM_HEAT_TRANSFER_HOHENMBERG, ICEM_HEAT_TRANSFER_ANNAND};
enum icem_exchange_area_models {ICEM_EXCHANGE_AREA_LATERAL, ICEM_EXCHANGE_AREA_TOTAL};

class MyOdeSystem_ICEM_MultiZone : public BzzOdeSystemObject
{
private:

public:
	OpenSMOKE_ICEM_MultiZone *ptICEM_MultiZone;
	void assignICEM_MultiZone(OpenSMOKE_ICEM_MultiZone *icem);
	virtual void MyODEPrint(BzzVector &y, double eta);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};


class MyOdeSystem_Reaction_ICEM_MultiZone : public BzzOdeSystemObject
{
private:

public:
	OpenSMOKE_ICEM_MultiZone *ptICEM_MultiZone;
	void assignICEM_MultiZone(OpenSMOKE_ICEM_MultiZone *icem);
	virtual void MyODEPrint(BzzVector &y, double eta);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class MyDaeSystem_ICEM_MultiZone : public BzzDaeSystemObject
{
private:

public:
	OpenSMOKE_ICEM_MultiZone *ptICEM_MultiZone;
	void assignICEM_MultiZone(OpenSMOKE_ICEM_MultiZone *icem);
	virtual void MyDAEPrint(BzzVector &y, double eta);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class MyNLSSystem_ICEM_MultiZone_Gaussian : public BzzMyNonLinearSystemObject
{
public:
	void assignICEM(OpenSMOKE_ICEM_MultiZone *icem);

	OpenSMOKE_ICEM_MultiZone *ptICEM;
	virtual void GetResiduals(BzzVector &x, BzzVector &f);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Dictionary_ICEM_MultiZone : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_ICEM_MultiZone();
};


class OpenSMOKE_ICEM_MultiZone : public OpenSMOKE_0DReactor
{

public:	// Virtual functions (compulsory)

	OpenSMOKE_ICEM_MultiZone();

	void Lock();
	void VideoFinalResult();
	void VideoGeneralInfo();
	void VideoSummary();
	void DefineFromFile(const std::string inputFile);
	void Solve();
	void AssignEnd(const std::string units, const double value);
	void EnergyAnalysis(OpenSMOKE_GasStream &outletStream);
	void EquationSystemPrint(BzzVector &y, double eta);
	void ODEPrint(BzzVector &y, double t);
	void SummaryOnFile();

public:	// Pseudo-Private Functions

	void ODESystem_ICEM_MultiZone(BzzVector &x, double t, BzzVector &f);
	void ODESystem_Reaction_ICEM_MultiZone(BzzVector &x, double t, BzzVector &f);
	void DAESystem_ICEM_MultiZone(BzzVector &x, double t, BzzVector &f);
	void NLSSystem_ICEM_Gaussian(BzzVector &x, BzzVector &f);

	MyOdeSystem_ICEM_MultiZone			ODE_ICEM_MultiZone_Object;
	MyOdeSystem_Reaction_ICEM_MultiZone ODE_Reaction_ICEM_MultiZone_Object;
	MyDaeSystem_ICEM_MultiZone			DAE_ICEM_MultiZone_Object;
	MyNLSSystem_ICEM_MultiZone_Gaussian NLS_Gaussian_MultiZone_Object;
	
private: // Virtual functions (compulsory)
	
	void LabelODEFile();
	void LabelEquationSystemFile();
	void Initialize();
	void UpdateProperties_isothermal(int memoIndex);
	void UpdateProperties(int jacobianIndex, int indexT);
	void UpdatePropertiesSingleReactor();
	void UpdatePropertiesMixing();
	void UpdateMeanMassFractions();
	
	void UpdateHeatFlux(const double tau, const double csi);
	void UpdateExchangeArea(const double tau, const double csi);
	void UpdateTwoEquationModel(BzzVector &y, BzzVector &dy);


public: // Specific functions 
	
	void AssignEngineModel(const icem_multizone_models value);
	void AssignClearanceVolume(const std::string units, const double value);
	void AssignTotalDisplacement(const std::string units, const double value);
	void AssignRotationRate(const std::string units, const double value);
	void AssignCompressionRatio(const double value);
	void AssignArmRatio(const double value);
	void AssignExhaustRatio(const double value);
	void AssignStartAngle(const std::string units, const double value);
	void AssignEndAngle(const std::string units, const double value);
	void AssignDiameter(const std::string units, const double value);
	void AssignNumberOfZones(const int value);
	void AssignCrevicesVolume(const std::string units, const double value);
	void AssignNumberOfCycles(const int value);
	
	void SetConstantExchangeArea(const double value, const std::string units);
	void SetUserDefinedExchangeArea(const std::string fileName);
	void UnsetUserDefinedExchangeArea();
	void SetExchangeAreaModel(const icem_exchange_area_models value);
	void SetRelaxation(const int value);


private: // Specific data

	int	   number_of_cycles;
	double teta;
//	double teta_deg;
	double rotation_rate;
	double volume_clearance;
	double volume_displacement_total;
	double volume_displacement_swept;
	double arm_ratio;
	double compression_ratio;
	double exhaust_ratio;
	double start_angle;
	double end_angle;
	double delta_step_reference;
	double volume_total_starting;

	// Geometry
	OpenSMOKE_UD_Profile	ud_A_profile;			// UD A profile

private:

	bool assignedClearanceVolume;
	bool assignedCompressionRatio;
	bool assignedArmRatio;
	bool assignedExhaustRatio;
	bool assignedStartAngle;
	bool assignedEndAngle;
	bool assignedRotationRate;
	bool assignedDiameter;
	bool assignedCrevicesVolume;
	bool assignedTotalDisplacement;
	bool assignedNumberOfZones;
	bool assignedNumberOfCycles;

	void MemoryAllocation();
	void GeometrySummary();
	void InitializeGeometry();
	void UpdateAngleAndVolume(const double teta);
	void UpdateVolumeDerivatives(const double teta);
	void UpdateProperties(const int jacobianIndex, BzzVectorInt &jacobianVariables, int dimBlock, const int indexT);
	void UpdateTimeDerivative(const double t);
	void UpdateMassFluxes();
	void UpdateInterfaces();
	void UpdateEnthalpyFluxes();

	void Equations_Species();
	void Equations_Energy();
	void Equations_Pressure();

private:

	// Cylinder data
	double height_cylinder;				// Cylinder (time-dependent)
	double height_cylinder_maximum;
	double height_cylinder_minimum;
	double volume_crevices;				// Crevices volume
	double thickness;					// Zone thickness
	double height_core_minimum;			// Core minimum height
	double diameter_cylinder;			// Cylinder (constant)
	double area_crevices;				// Crevices (constant)
	double area_cylinder_total;			// Cylinder	(time-dependent)
	double area_cylinder_base;			// Cylinder (time-dependent)
	double area_cylinder_lateral;		// Cylinder (time-dependent)
	double volume_total;				// Cylinder and crevices
	double volume_cylinder_total;		// Cylinder (time-dependent)
	double velocity_piston;
	double arm_La;
	double arm_Lc;

	// Zones data
	BzzVector zone_volume;				// Zone volumes
	BzzVector area_i;					// Zone internal areas
	BzzVector area_e;					// Zone external areas
	BzzVector diameter_i;				// Zone internal diameter	
	BzzVector diameter_e;				// Zone external diameter
	BzzVector height_i;					// Zone internal height
	BzzVector height_e;					// Zone external height
	BzzVector area_crown;				// Zone crown area

private:

	int N;
	int nEquations;
	int dimBlock;

	// Main unknowns
	BzzVector zone_temperature;			// Zone temperature
	BzzMatrix zone_omega;				// Zone mass fractions
	BzzVector zone_pressure;			// Zone total mole
	BzzVector zone_g;					// 
	BzzVector zone_heat_release;		// 

	// Equations
	BzzMatrix domega_over_dt;
	BzzVector dtemperature_over_dt;
	BzzVector dmass_over_dt;
	BzzVector p_equation;
	BzzVector g_equation;
	BzzVector dvolume_over_dt;

	// Zone properties (I)
	BzzVector zone_mass;				// Zone total mass
	BzzVector zone_mole;				// Zone total mole
	BzzVector zone_mw;					// Zone molecular weight
	BzzMatrix zone_x;					// Zone mole fractions
	BzzMatrix zone_c;					// Zone concentrations

	// Zone properties (II)
	BzzVector zone_rho;
	BzzVector zone_ctot;
	BzzVector zone_Cp;
	BzzVector zone_Cv;
	BzzVector zone_mu;
	BzzVector zone_lambda;
	BzzVector zone_QReaction;
	BzzVector zone_UReaction;
	BzzMatrix zone_hmass;
	BzzMatrix zone_umass;
	BzzMatrix zone_R;
	BzzVector zone_Qe;
	BzzVector zone_enthalpy;
	BzzVector zone_internal_energy;
	BzzVector dummy_vector_nc;
	BzzVector zone_area_effective;		// Zone area effective
	BzzVector zone_volume_old;			// Zone volume (old value)
	BzzVector zone_mass_old;			// Zone mass (old value)

	double work;						// [J]
	double volume_total_old;			// [m3]

	bool iRelaxation;
	int relaxation_cycles;
	int relaxation_index;
	BzzMatrix *relaxation_matrix;
	bool iMonteCarloMixing;
	double TauMixing;
	double Cphi;

private:

	// Fluxes
	BzzVector zone_m_in;				// Zone inlet mass fluxes
	BzzVector zone_h_in;				//
	BzzMatrix zone_hmass_in;			// Zone inlet enthalpy 
	BzzMatrix zone_omega_in;			// Zone inlet mass fractions

	// Interface variables
	BzzVector zone_lambda_in;
	BzzVector zone_temperature_gradient_in;
	BzzVector yPlus_in;
	BzzVector viscosity_ratio_in;
	BzzVector zone_Qdiffusion;

private:

	// Total variables
	double mass_total;
	double mass_total_initial;
	double enthalpy_total;
	double internal_energy_total;
	double QReaction_total;
	double HeatRelease_total;

	// Mean properties
	double mean_rho;
	double mean_mw;
	BzzVector mean_omega;
	BzzVector mean_x;

	BzzMatrix zone_omega_initial;
	BzzVector zone_temperature_initial;
	BzzVector zone_enthalpy_initial;
	BzzVector zone_g_initial;
	BzzVector zone_pressure_initial;
	BzzVector zone_heat_release_initial;
	
private:

	// Indices
	int jGlobal;
	int jCrevices;
	int jOutmost;

	BzzVectorInt	index_species_zones;
	vector<string>	names_species_zones;

private:

	// Files
	ofstream fZones;
	ofstream fGeometry;
	ofstream fMass;
	ofstream fRecycle;
	ofstream fExhaust;

	// File functions
	void PrepareAdditionalFiles();
	void LabelAdditionalFiles();
	void PrintAdditionalFiles(const double t, const double teta);

private:

	// Models
	icem_multizone_models		icem_multizone_model;
	icem_heat_transfer_models	icem_heat_transfer_model;
	icem_exchange_area_models   icem_exchange_area_model;

private:

	bool iHot;
	bool iVerboseZones;
	bool iVerboseGeometry;
	bool iVerboseMasses;
	bool iPlugFlowExhaustGases;
	bool iPlugFlowExhaustGasesInletTemperature;

	void SetVerboseGeometry();
	void SetVerboseMasses();
	void SetVerboseZones(const vector<string> _names);
	void SetColdSimulation();
	void SetZoneInitialConditions(const vector<string> string_vector);
	void SetGaussianInitialConditions(const vector<string> string_vector);
	void SetEquations();
	void SetMinimumValues(BzzVector &xMin);
	void SetMaximumValues(BzzVector &xMax);
	void SetInitialValues(BzzVector &xInitial);
	void ReSetInitialValues(BzzVector &xInitial);
	void SetDifferentialEquations(BzzVectorInt &fDifferential);
	void SetRelativeTolerances(BzzVector &fRelTol);
	void SetMonteCarloMixing();
	void SetMicromixingTime(const double value, const std::string units);
	void SetMicromixingConstant(const double value);
	void SetPlugFlowTime(const std::string units, const double value);
	void SetPlugFlowInletTemperature(const std::string units, const double value);
	void ResetSpecies(BzzMatrix &species);
	void SetResetSpecies(const vector<string> names);

	double gaussian_mean;
	double gaussian_sigma;
	double plugFlowExhaustGasesTime;
	double plugFlowExhaustGasesInletTemperature;

	BzzVector indicesResetSpecies;
	bool iResetSpecies;

private:

	// Kinetic Maps
	BzzMatrix zone_CpMap;
	BzzMatrix zone_k1Map;
	BzzMatrix zone_k2Map;
	BzzMatrix zone_uKeqMap;
	BzzMatrix zone_logFcentMap;
	BzzMatrix zone_reactionDHMap;
	BzzMatrix zone_reactionDSMap;

	// Properties maps
	BzzMatrix zone_lambdaMap;
	BzzMatrix zone_muMap;

private:

	double Re;						// [-]
	double heat_transfer_velocity;	// [m/s]
	double pressure_reference;		// [Pa]
	double volume_reference;		// [m3]
	double temperature_reference;	// [K]

	void UpdateHeatTransfer();
	void SetHeatTransferModel_Woschni();
	void SetHeatTransferModel_Assanis();
	void SetHeatTransferModel_Hohenberg();
	void SetHeatTransferModel_Annand();

	void SetHeatTransferModel_WoschniParameters(const vector<string> string_vector);
	void SetHeatTransferModel_AssanisParameters(const vector<string> string_vector);
	void SetHeatTransferModel_HohenbergParameters(const vector<string> string_vector);
	void SetHeatTransferModel_AnnandParameters(const vector<string> string_vector);

private:

	double woschni_alfa;		// [-]
	double woschni_p;			// [-]
	double woschni_l;			// [-]
	double woschni_t;			// [-]
	double woschni_v;			// [-]

	double woschni_c1_phase1;	// [-]
	double woschni_c1_phase2;	// [-]

	double woschni_c2_phase1;	// [m/s/K]
	double woschni_c2_phase2;	// [m/s/K]

private:

	double assanis_alfa;		// [-]
	double assanis_p;			// [-]
	double assanis_l;			// [-]
	double assanis_t;			// [-]
	double assanis_v;			// [-]

	double assanis_c1_phase1;	// [-]
	double assanis_c1_phase2;	// [-]

	double assanis_c2_phase1;	// [m/s/K]
	double assanis_c2_phase2;	// [m/s/K]

private:

	double hohenberg_alfa;

	double hohenberg_p;
	double hohenberg_V;
	double hohenberg_t;
	double hohenberg_b;
	double hohenberg_v;

private:

	double annand_alfa;
	double annand_Re;

private:

	int			zone_count;
	BzzVector	zone_list_temperature;
	BzzVector	zone_list_mass_fraction;
	BzzVector	zone_list_external_exchange;

private:
	
	static const double rad_to_degrees;
	static const double degrees_to_rad;

	static const double min_temperature;
	static const double max_temperature;

	static const double min_pressure;
	static const double max_pressure;

private:
	
	void SaveOnBinaryFile(BzzSave &fOutput);
	void SaveOnBinaryFile(const std::string filename);
	void UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr);

	BzzDaeObject		dae_single;
	BzzDaeSparseObject	dae_sparse;
};

#endif // OPENSMOKE_ICEM_MULTIZONE_H
