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

#ifndef OPENSMOKE_0DREACTOR_H
#define OPENSMOKE_0DREACTOR_H

#include "OpenSMOKE.hpp"
#include "addons/OpenSMOKE_RateOfProductionAnalysis.h"
#include "addons/OpenSMOKE_SensitivityAnalysis.h"
#include "addons/OpenSMOKE_ElementFluxAnalysis.h"

class OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D;

class OpenSMOKE_0DReactor
{
public:

    // 1. Assignements (compulsory)
	void AssignSoot2EModel(OpenSMOKE_2EModel &_soot2EModel);
    void AssignGlobalKineticScheme(OpenSMOKE_GlobalKinetics &_global);
	void AssignKineticScheme(OpenSMOKE_ReactingGas &_mix);
	void AssignEnergy(const bool iEnergy);
	void AssignInletFlows(OpenSMOKE_GasStream &stream);
	virtual void AssignEnd(const string units, const double value) = 0;
	virtual void Lock() = 0;

    // 2. Option settings
	void SetName(const string name);
    void SetOptions(string option);
    void SetVideoOptions(const int _nVideoSteps);
    void SetFileOptions(const int _nFileSteps);
    void SetOutputFolder(const string _outputFolderName);
    void SetGlobalKinetics();
	void SetTwoEquationModel();
	void SetVerbose();
	void SetVerboseEnergy();
	void SetHistory();
    void SetRelativeTolerance(const double relative);
    void SetAbsoluteTolerance(const double absolute);
    void SetUserDefinedTemperature(const string fileName);
	void SetUserDefinedTemperature(OpenSMOKE_UD_Profile &udp);
	void SetViscosity(const double value, const string units);
	void SetKeySpecies(const vector<string> key_species);
	
	void SetUserDefinedHeatFlux(const string fileName);
	void SetUserDefinedHeatExchangeCoefficient(const string fileName);
	void SetUserDefinedAmbientTemperature(const string fileName);
	void SetConstantHeatFlux(const double value, const string units);
	void SetConstantHeatExchangeCoefficient(const double value, const string units);
	void SetConstantAmbientTemperature(const double value, const string units);
	
	void SetReactionRatesOnFile(const vector<string> _names);
	void SetFormationRatesOnFile(const vector<string> _names);
	void SetSelectivitiesOnFile();
	void SetVerboseROPAOnFile(const vector<string> _names);
	void SetROPAOnFile();
	void SetSensitivityOnFile(const vector<string> _names);
	void SetSensitivityOptions(const vector<string> _names);
	void SetExperimentOnFile(const vector<string> _names);
	void SetElementFluxAnalysisOnFile(const vector<string> _names);
	void SetLocalElementFluxAnalysisOnFile(const vector<string> _names);
	void SetPostProcessing();
	void SetOutputSpecies(const vector<string> string_vector);

	// 3. Unset Options
    void UnsetGlobalKinetics();
    void UnsetTwoEquationModel();
    void UnsetVerbose();
    void UnsetHistory();
    void UnsetRelativeTolerance();
    void UnsetAbsoluteTolerance();
    void UnsetUserDefinedTemperature();
	void UnsetUserDefinedHeatFlux();
	void UnsetUserDefinedHeatExchangeCoefficient();
	void UnsetUserDefinedAmbientTemperature();
    void UnsetViscosity();
    void UnsetVerboseEnergy();
	void UnsetKeySpecies();

	// 4. Messages for user
	virtual void VideoFinalResult() = 0;
	virtual void VideoSummary() = 0;
	virtual void SummaryOnFile() = 0;

	// Constructor
    virtual void DefineFromFile(const string inputFile) = 0;

	// Solve Reactor
	virtual void Solve() = 0;

	// Final Analyses
	void OutletStream(OpenSMOKE_GasStream &outlet);
	void MassAnalysis(OpenSMOKE_GasStream &outletStream);
	virtual void EnergyAnalysis(OpenSMOKE_GasStream &outlet) = 0;

	// ODE print
	virtual void ODEPrint(BzzVector &y, double eta) = 0;
	
protected:

	// Error and warning messages
	void	ErrorMessage(const string message);
	void	WarningMessage(const string message);

	// Message functions
	virtual void VideoGeneralInfo() = 0;
	virtual void LabelODEFile() = 0;

	// Output preparation
	void PrepareFiles();
	void CloseFiles();

    // Updating properties at each ODE call
			void	Setup();
	virtual void    Initialize() = 0;
    virtual void    UpdateHeatFlux(const double tau, const double csi) = 0;
    virtual void    UpdateExchangeArea(const double tau, const double csi) = 0;
    virtual void    UpdateProperties(int indexJacobian, int indexT) = 0;
	virtual void    UpdateProperties_isothermal(int indexMemo) = 0;
	virtual void	UpdateTwoEquationModel(BzzVector &y, BzzVector &dy) = 0;

	// Properties
	double SimplifiedViscosity();

protected:

	// Error and warning messages
	string name_object;
	string class_name;
	string out_name;

	// Control Variables
	bool assignedKineticScheme;
	bool assignedGlobalKineticScheme;
	bool assignedSoot2EModel;
	bool assignedEnergy;
	bool assignedInletFlows;
	bool assignedEnd;

	// Options
	bool iEnergy;
	bool iGlobalKinetics;
	bool iTwoEquationModel;
    bool iHistory;
    bool iVerbose;
	bool iVerboseEnergy;
	bool iVerboseReactionRates;
	bool iVerboseFormationRates;
	bool iAssignedROPA;
	bool iVerboseAssignedROPA;
	bool iVerboseSelectivity;
	bool iVerboseSensitivity;
	bool iVerboseElementFluxAnalysis;
	bool iVerboseLocalElementFluxAnalysis;
	bool iRelativeTolerance;
	bool iAbsoluteTolerance;
	bool iUserDefinedTemperature;
	bool iUserDefinedViscosity;
    bool iKeySpecies;
	bool iVerboseSolid;
	bool iVerboseExperiment;

    ProfileKind iUserDefinedHeatFlux;
    ProfileKind iUserDefinedExchangeArea;
	ProfileKind iUserDefinedHeatExchangeCoefficient;
	ProfileKind iUserDefinedAmbientTemperature;

	// Properties evaluation
	int     indexProperties;
	int     indexTwoEquations;

	// Print Options
	string outputFolderName;
	string outputName;
	string outputSootName;
	string outputPAHName;
	string output2EName;
	string outputOSMName;
	string outputEnergyName;
	string outputReactionRatesName;
	string outputFormationRatesName;
	string outputSolidName;
	string outputConversionsName;
	string outputExperimentName;
	string outputSensitivityName;
	string outputSelectivityName;


    // Internal counters
    int countIterations;
	int countVideoSteps;
	int countFileSteps;

    // User defined options
	int     nVideoSteps;
	int     nFileSteps;
	double  relativeTolerance;
	double  absoluteTolerance;

    // Main objects
	OpenSMOKE_ReactingGas				*mix;						// Gas Mixture
	OpenSMOKE_GlobalKinetics			*global;					// Global Kinetics
	OpenSMOKE_2EModel					*soot2EModel;				// Two Equation Model
	OpenSMOKE_GasStream					*inletStream;				// Inlet stream
	OpenSMOKE_RateOfProductionAnalysis	 ropa;						// Rate of production analysis
	OpenSMOKE_SensitivityAnalysis		 Sensitivity;				// Sensitivity analysis
	OpenSMOKE_ElementFluxAnalysis		 ElementFluxAnalysis;		// Element Flux analysis
	
	OpenSMOKE_UD_Profile	ud_temperature_profile;	// UD temperature profile
	OpenSMOKE_UD_Profile	ud_Qe_profile;			// UD Qe profile
	OpenSMOKE_UD_Profile	ud_U_profile;			// UD Qe profile
	OpenSMOKE_UD_Profile	ud_Tambient_profile;	// UD Ae profile

	BzzVectorInt	iOutputSpecies;
	vector<string>	namesOutputSpecies;
	vector<string>	sensitivityOptions;

	// Kinetic scheme data
	int NC;							// Number of Species
	int NR;							// Number of Reactions
	
	
	double timeOld;
	BzzOdeStiffObject o;

	// Unknowns and residuals
	BzzVector	omega;		    // Mass fractions [-]
    BzzVector domega;
	double          T;				// Reactor temperature [K]
	double          dT;


    // Main data
	double TauTotal;				// Total time of simulation [s]
    double P;				        // Reactor pressure [Pa]
	double rho;						// Reactor density [kg/m3]
	double MWtot;					// Reactor molecular weight [kg/kmol]
	double Cp;						// Constant pressure specific heat [J/kg/K]
	double Cv;						// Constant volume specific heat [J/kg/K]
	double QReaction;				// Reactor reaction heat [J/m3/s]

	BzzVector		x;			// Mole fractions [-]
	BzzVector		c;			// Concentrations [kmol/m3]
	BzzVector		R;			// Formation Rates [kg/m3/s]
	BzzVector		Rtilde;		// Formation rates [kmol/m3/s]

    // Reactor geometry
    double Qe;						// External heat flux [J/m2/s]
    double U;						// Heat exchange coefficient [W/m2/K]
    double Tambient;				// Ambient temperature [K]
    double Viscosity;				// Inlet stream viscosity [K]

	// Kinetic Maps
	BzzVector CpMap;
	BzzVector k1Map;
	BzzVector k2Map;
	BzzVector uKeqMap;
	BzzVector logFcentMap;
	BzzVector reactionDHMap;
	BzzVector reactionDSMap;

	// Key species
	vector<string>	key_species_names;
	vector<int>		key_species_index;

    // Output files
	ofstream    fPAH;
	ofstream    fSoot;
    ofstream    fSoot2E;
    ofstream    fOutput;
    ofstream    fEnergy;
	ofstream	fReactionRates;
	ofstream	fFormationRates;
	ofstream	fSolid;
	ofstream	fSensitivity;
	ofstream	fConversions;
	ofstream	fSelectivity;

	
	// Additional variables 
	BzzVectorInt index_reaction_rates;
	BzzVectorInt index_formation_rates;
	BzzVectorInt index_sensitivity;
	BzzVectorInt index_ROPA;
	BzzVectorInt index_Experiment;
	BzzVector    index_local_ElementFluxAnalysis;
	BzzVector    selectivity_mass;
	vector<string> names_ROPA;
	vector<string> names_reaction_rates;
	vector<string> names_formation_rates;
	vector<string> names_sensitivity;
	vector<string> names_Experiment;

public:

    BzzVector P_History;
    BzzVector Tau_History;
    BzzVector Csi_History;
    BzzVector T_History;
    BzzMatrix mass_History;
    BzzMatrix mole_History;

    int countGlobalIterations;

	inline bool isVerbose() { return iVerbose; };

protected:

	// Sensitivity
	OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D *sensitivity_fast;


protected:

    static const double   ONE;
    static const double   ZERO;
    static const double   MAXNUMBER;
    static const int      MAX_TIME_STEPS;
    static const int      MINUSONE;
	static const double   MINPRESSURE;
	static const double   MAXERRORONMASSFRACTIONS;
};

#endif // OPENSMOKE_0DREACTOR_H
