/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci   	   *
 *   alberto.cuoci@polimi.it   						   *
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

#ifndef OPENSMOKE_IDEALGAS_H
#define OPENSMOKE_IDEALGAS_H

#include "Linux_Windows.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
//#include "addons/OpenSMOKE_SootManager.h"
#include "engine/OpenSMOKE_EquilibriumStanjan.h"
using namespace std;

class OpenSMOKE_PolimiSoot;

class OpenSMOKE_IdealGas
{
protected:
	int NC;

private:


	// Variabili che indicano se la corrispondente proprieta deve essere calcolata (1)
	int iMixViscositySimplified;

	BzzVector *CpHT;	// calore specifico a pressione costante (High Temperature)
	BzzVector *CpLT;	// calore specifico a pressione costante (Low Temperature)

	BzzVector		sumK;
	BzzVector		sqrtEta;
	BzzVector		usqrtEta;
	BzzVector		PMRatio1su4;
	BzzVector		delta_phi;
	BzzVector		phi_eta_sup;
	BzzVector		phi_eta_inf;
	BzzVector		sqrtPMRatio_inf;
	BzzVector		sqrtPMRatio_sup;
		
	BzzVectorInt	iD;
	BzzVector		x_aux;
	BzzVector Tintervals;

	void Allocate(void);
	void PostProcessMeanTemperatureForThermodynamicProperties();

public:

	// Matrici delle costanti di fitting
	BzzMatrix fittingEta;
	BzzMatrix fittingLambda;
	BzzMatrix fittingDbinary;
	BzzMatrix fittingTetaBinary;

	// Variabili di ingresso per la lettura dei dati termodinamici
	BzzVector *aDH; // alta T
	BzzVector *bDH; // bassa T
	BzzVector *aDS; // alta T
	BzzVector *bDS; // bassa T
	BzzVector T1;
	BzzVector T2;
	BzzVector T3;


	// Names and molecular weigths	
	std::string				*names;
	BzzVector		M;		// Peso molecolare [g/mol]
	BzzVector		uM;		// Peso molecolare [g/mol]
	BzzVector		Cp;		// Species Specific Heat (Constant pressure)
	BzzVector		lambda;	// Transport properties of single components
	BzzMatrix		Djk;	// Transport properties of single components
	BzzVector		eta;	// Transport properties of single components
	BzzMatrix		Tetakj;	// Thermal diffusion ratios

	BzzVector		Cv;		// Species Specific Heat (Constant volume)
	BzzVector		h;		// Species Enthalpy (dimensionless)
	BzzVector		s;		// Species Enthalpy (dimensionless)

	BzzVectorInt iThermalDiffusionRatios;

	void ReadVersion(BzzLoad &binaryFile);
	void Setup(BzzLoad &binaryFile);
	void Fitting(BzzLoad &binaryFile);
	void SetPolimiSoot(	const unsigned int bin_index_zero, const double bin_density_A, 
						const unsigned int bin_index_final, const double bin_density_B,
						const double Df);
	void SetPolimiSoot(	);

	void CorrectSpeciesDiffusivityForSoot();

	// Utilities
	double	GetMWFromMassFractions(BzzVector &y);
	double	GetMWFromMoleFractions(BzzVector &x);
	void	GetMWAndMassFractionsFromMoleFractions(double &MWmix, BzzVector &y, BzzVector &x);
	void	GetMWAndMoleFractionsFromMassFractions(double &MWmix, BzzVector &x, BzzVector &y);
	void	GetMassFractionsFromMoleFractionsAndMW(BzzVector &y, BzzVector &x, const double MWmix);
	void	GetMoleFractionsFromMassFractionsAndMW(BzzVector &x, BzzVector &y, const double MWmix);

	void	GetMWAndMoleFractionsFromMassFractions(BzzVector &MWmix, BzzMatrix &x, BzzMatrix &y);
	void	GetMWAndMassFractionsFromMoleFractions(BzzVector &MWmix, BzzMatrix &y, BzzMatrix &x);

	// Proprieta singole specie (termodinamiche)
	void SpeciesCp(double T);
	void SpeciesCv(void);
	
	void SpeciesEnthalpyAndEntropy(double T);
	void SpeciesEnthalpy(double T);
	void SpeciesEntropy(double T);
	
	
	void	SpeciesViscosityFromFitting(double T);
	void	SpeciesConductivityFromFitting(double T);
	void	SpeciesDiffusivityFromFitting(double T, double P);
	void	SpeciesThermalDiffusionRatiosFromFitting(double T);

	// Calcolo delle proprieta di miscela
	void	MixDiffusionCoefficient_FromMolarFractionsAndPM(BzzVector &Dkm, BzzVector &x);
	double	MixViscosity_FromMolarFractions(BzzVector &x);
	double	MixConductivity_FromMolarFractions(BzzVector &x);
	void	MixThermalDiffusionRatios(BzzVector &Teta, BzzVector &x);

	// Thermodynamics: mixture properties
	double MixCp_FromMassFractions(BzzVector &y);
	double MixCp_FromMoleFractions(BzzVector &x);
	double MixCv_FromCpMix(const double CpMix, const double MWmix);
	
	int recognize_species(char* name);
	int recognize_species_without_exit(char* name);
    int recognize_species(const std::string name);
	int recognize_species_without_exit(const std::string name);

	int NumberOfSpecies();
	int	NumberOfElements();

	OpenSMOKE_IdealGas();

	void SetMinimumTemperature(const double tmin);
	void SetMaximumTemperature(const double tmax);
	void SetName(const std::string name);
	void SetAuthorName(const std::string name);

	std::string name_object;
	std::string name_author;
	std::string building_date;

	double TMIN;
	double TMAX;

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);


	// Managers
//	OpenSMOKE_PAHManager	pah_manager;	// PAH manager
	OpenSMOKE_PolimiSoot*	polimiSoot;		// Soot manager

	void CorrectBinaryDiffusionCoefficients(const int j, const int k, const double correction);
	void CorrectFormationEnthalpy(const int j, const double correction);

	bool iTransportMode;

	bool TransportMode();

	void SpeciesCp_ForcingCorrelation(BzzVector &T, const char low_or_high, BzzVector &cpi, BzzVector &d1cpi, BzzVector &d2cpi, BzzVector &d3cpi, BzzVector &d4cpi);
	bool binary_version_;
	void SetASCII();

public:

	versionKinetics	iVersion;

	vector<string>	list_of_elements;
	BzzMatrix	elements;
	BzzVector	m_elements;
	BzzVector	um_elements;

	int  recognize_element(const std::string name);
	int  recognize_element_without_error(const std::string name);
	void GetElementalMoleFractionsFromSpeciesMoleFractions(BzzVector &x_elemental, BzzVector &x);
	void GetElementalMassFractionsFromSpeciesMassFractions(BzzVector &omega_elemental, BzzVector &omega);
	void GetElementalMoleFractionsFromSpeciesMoleFractions(BzzMatrix &x_elemental, BzzMatrix &x);
	void GetElementalMassFractionsFromSpeciesMassFractions(BzzMatrix &omega_elemental, BzzMatrix &omega);
	void GetElementalMassFractionsFromElementalMoleFractionsAndMW(BzzVector &omega_elemental, BzzVector &x_elemental, const double MWmix);
	void GetElementalMoleFractionsFromElementalMassFractionsAndMW(BzzVector &x_elemental, BzzVector &omega_elemental, const double MWmix);
	void GetMWAndElementalMassFractionsFromElementalMoleFractions(double &MWmix, BzzVector &omega_elemental, BzzVector &x_elemental);
	void GetMWAndElementalMoleFractionsFromElementalMassFractions(double &MWmix,BzzVector &x_elemental, BzzVector &omega_elemental);
	void GetElementalMassFractionsFromSpeciesMoleFractions(BzzVector &omega_elemental, BzzVector &x);

	double GetMixtureFraction(BzzVector &omega, BzzVector &omegaFuel, BzzVector &omegaAir);

	void VerboseDataSpecies(ofstream &fOutput);
	void VerboseDataSpecies(ofstream &fOutput, const int k, BzzVector &T_Vector, BzzMatrix &Cp_Matrix, BzzMatrix &DH_Matrix, BzzMatrix &DS_Matrix, BzzMatrix &S_Matrix, BzzVector &H_Tref, BzzVector &G_Tref);


public:

	// Standard properties [J, kg, K]
	void GetStandardEnthalpy_Mass(BzzVector &h, const double T);
	void GetStandardInternalEnergy_Mass(BzzVector &u, const double T);
	void GetStandardEntropy_Mass(BzzVector &s, const double T);
	void GetStandardGibbsFreeEnergy_Mass(BzzVector &g, const double T);
	void GetStandardHelmotzFreeEnergy_Mass(BzzVector &a, const double T);

	// Mixture-averaged properties [J, kg, K]
	void GetMixAveragedEnthalpy_Mass(BzzVector &h, const double T);
	void GetMixAveragedInternalEnergy_Mass(BzzVector &u, const double T);
	void GetMixAveragedEntropy_Mass(BzzVector &s, const double P_Pa, const double T, BzzVector &omega);
	void GetMixAveragedGibbsFreeEnergy_Mass(BzzVector &g, const double P_Pa, const double T, BzzVector &omega);
	void GetMixAveragedHelmotzFreeEnergy_Mass(BzzVector &a, const double P_Pa, const double T, BzzVector &omega);

	// Mixture formation enthalpy [J, kmol, K]
	double GetMixFormationEnthalpy_Mass(BzzVector &omega);
	double GetMixFormationInternalEnergy_Mass(BzzVector &omega);
	double GetMixFormationEntropy_Mass(const double P_Pa, BzzVector &omega);
	double GetMixFormationGibbsFreeEnergy_Mass(const double P_Pa, BzzVector &omega);
	double GetMixFormationHelmotzFreeEnergy_Mass(const double P_Pa, BzzVector &omega);

	// Mixture formation enthalpy [J, kmol, K]	
	double GetMixEnthalpy_Mass(const double T, BzzVector &omega);								// [J/kg]
	double GetMixInternalEnergy_Mass(const double T, BzzVector &omega);						// [J/kg]
	double GetMixEntropy_Mass(const double P_Pa, const double T, BzzVector &omega);			// [J/kg]
	double GetMixGibbsFreeEnergy_Mass(const double P_Pa, const double T, BzzVector &omega);	// [J/kg]
	double GetMixHelmotzFreeEnergy_Mass(const double P_Pa, const double T, BzzVector &omega);	// [J/kg]
	
	// Standard properties [J, kg, K]
	void GetStandardEnthalpy_Mole(BzzVector &h, const double T);
	void GetStandardInternalEnergy_Mole(BzzVector &u, const double T);
	void GetStandardEntropy_Mole(BzzVector &s, const double T);
	void GetStandardGibbsFreeEnergy_Mole(BzzVector &g, const double T);
	void GetStandardHelmotzFreeEnergy_Mole(BzzVector &a, const double T);

	// Mixture-averaged properties [J, kg, K]
	void GetMixAveragedEnthalpy_Mole(BzzVector &h, const double T);
	void GetMixAveragedInternalEnergy_Mole(BzzVector &u, const double T);
	void GetMixAveragedEntropy_Mole(BzzVector &s, const double P_Pa, const double T, BzzVector &x);
	void GetMixAveragedGibbsFreeEnergy_Mole(BzzVector &g, const double P_Pa, const double T, BzzVector &x);
	void GetMixAveragedHelmotzFreeEnergy_Mole(BzzVector &a, const double P_Pa, const double T, BzzVector &x);

	// Mixture formation enthalpy [J, kmol, K]
	double GetMixFormationEnthalpy_Mole(BzzVector &x);
	double GetMixFormationInternalEnergy_Mole(BzzVector &x);
	double GetMixFormationEntropy_Mole(const double P_Pa, BzzVector &x);
	double GetMixFormationGibbsFreeEnergy_Mole(const double P_Pa, BzzVector &x);
	double GetMixFormationHelmotzFreeEnergy_Mole(const double P_Pa, BzzVector &x);

	// Mixture formation enthalpy [J, kmol, K]	
	double GetMixEnthalpy_Mole(const double T, BzzVector &x);								// [J/kmol]
	double GetMixInternalEnergy_Mole(const double T, BzzVector &x);						// [J/kmol]
	double GetMixEntropy_Mole(const double P_Pa, const double T, BzzVector &x);			// [?]
	double GetMixGibbsFreeEnergy_Mole(const double P_Pa, const double T, BzzVector &x);	// [?]
	double GetMixHelmotzFreeEnergy_Mole(const double P_Pa, const double T, BzzVector &x);	// [?]

	// Inverse Functions
	double GetTemperatureFromMassEnthalpyAndMoleFractions(const double TFirstGuess, const double Hmass, BzzVector &xmole);
	double GetTemperatureFromMassEnthalpyAndMassFractions(const double TFirstGuess, const double Hmass, BzzVector &omega);
	
	double GetTemperatureFromMassEnthalpyAndMoleFractions(const double TFirstGuess, const double Hmass, BzzVector &xmole, const double dtmax);
	double GetTemperatureFromMassEnthalpyAndMassFractions(const double TFirstGuess, const double Hmass, BzzVector &omega, const double dtmax);

	double GetTemperatureFromMassEnthalpyAndMoleFractions(const double TFirstGuess, const double Hmass, BzzVector &xmole, 
									const double max_residual, const unsigned int max_newton_iterations, 
									const unsigned int max_searching_iterations, const double max_temperature_increment);

	double GetTemperatureFromMassEnthalpyAndMoleFractions_NewVersion(const double TFirstGuess, const double Hmass, BzzVector &xmole, 
									const double max_residual, const unsigned int max_newton_iterations, 
									const unsigned int max_searching_iterations, const double max_temperature_increment);


	double GetTemperatureFromMassEnthalpyAndMassFractions(const double TFirstGuess, const double Hmass, BzzVector &omega, 
									const double max_residual, const unsigned int max_newton_iterations, 
									const unsigned int max_searching_iterations, const double max_temperature_increment);

	double GetTemperatureFromMassEnthalpyAndMassFractions_NewVersion(const double TFirstGuess, const double Hmass, BzzVector &omega, 
									const double max_residual, const unsigned int max_newton_iterations, 
									const unsigned int max_searching_iterations, const double max_temperature_increment);

	double GetTemperatureFromMassInternalEnergyAndMoleFractions(const double TFirstGuess, const double Umass, BzzVector &xmole);
	double GetTemperatureFromMassEntropyAndMoleFractions(const double TFirstGuess, const double P_Pascal, const double Smass, BzzVector &xmole);

	double GetPressureFromMassEntropyAndMoleFractions(const double T, const double Smass, BzzVector &xmole);

	// Get Composition
	BzzVector GetMoleFractionsFromEquivalenceRatio(const double equivalence_ratio, const std::string fuel_name);
	BzzVector GetMoleFractionsFromEquivalenceRatio(const double equivalence_ratio,	const vector<string> fuel_names,		BzzVector &moles_fuel,
																					const vector<string> oxidizer_names,	BzzVector &moles_oxidizer);


	// Equilibrium
	void Equilibrium_TP(BzzVector &xmole_E, double &N_E, const double T, const double P_Pa, BzzVector &x_elements, const bool iVerbose);
	void Equilibrium_HP(double &T_E, BzzVector &xmole_E, double &N_E, const double Hmass, const double P_Pa, BzzVector &x_elements, const bool iVerbose, int &flag);
	void Equilibrium_SP(double &T_E, BzzVector &xmole_E, double &N_E, const double Smass, const double P_Pa, BzzVector &x_elements, const bool iVerbose);
	void Equilibrium_TRHO(double &P_E_Pascal, BzzVector &x_E, double &N_E, const double T, const double rho, BzzVector &x_elements, const bool iVerbose);
	void Equilibrium_URHO(double &T_E, double &P_E_Pascal, BzzVector &xmole_E, double &N_E, const double UMass, const double rho, BzzVector &x_elements, const bool iVerbose);
	void Equilibrium_SRHO(double &T_E, double &P_E_Pascal, BzzVector &xmole_E, double &N_E, const double SMass, const double rho, BzzVector &x_elements, const bool iVerbose);
	
	void VerboseDataSpeciesCheckConsistency(ofstream &fOutput);

private:
	static const int	N_LOOKUP;

protected:
	static const double	logCATM;
	

};


class OpenSMOKE_Equilibrium_SRHO_MyNonLinearSystem: public BzzMyNonLinearSystemObject
{
public:
	void assignIdealGas(OpenSMOKE_IdealGas *mix, const double _Smass, const double _rho, BzzVector &_xmole, const bool iVerbose);

	OpenSMOKE_IdealGas *ptMix;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void) {};

public:
	double T_E;
	double P_E;
	double N_E;
	BzzVector xmole_E;

private:
	double Smass;
	double rho;
	BzzVector x_elements;

	OpenSMOKE_EquilibriumStanjan equilibrium;
};

class OpenSMOKE_Equilibrium_URHO_MyNonLinearSystem: public BzzMyNonLinearSystemObject
{
public:
	void assignIdealGas(OpenSMOKE_IdealGas *mix, const double _Umass, const double _rho, BzzVector &_xmole, const bool iVerbose);

	OpenSMOKE_IdealGas *ptMix;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void) {};

public:
	double T_E;
	double P_E;
	double N_E;
	BzzVector xmole_E;

private:
	double Umass;
	double rho;
	BzzVector x_elements;

	OpenSMOKE_EquilibriumStanjan equilibrium;
};

#endif // OPENSMOKE_IDEALGAS_H


