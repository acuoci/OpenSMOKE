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

#ifndef OPENSMOKE_BATCH_H
#define OPENSMOKE_BATCH_H

#include "OpenSMOKE.hpp"

class   OpenSMOKE_Batch;
void    ODEPrintExternal(BzzVectorDouble &omega, double t);


class MyOdeSystem_Batch : public BzzOdeSystemDoubleObject
{
private:
	bool iEnergy;
	bool iConstantPressure;

public:
	OpenSMOKE_Batch *ptBatch;
	void assignBatch(OpenSMOKE_Batch *batch, const bool _iEnergy, const bool _iConstantPressure);
	virtual void GetSystemFunctions(BzzVectorDouble &y,double t,BzzVectorDouble &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Dictionary_Batch : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_Batch();
};

class OpenSMOKE_Batch
{
public:

    // Default constructor
    OpenSMOKE_Batch();

    // Assignements (compulsory)
    void AssignGlobalKineticScheme(OpenSMOKE_GlobalKinetics &_global);
	void AssignSoot2EModel(OpenSMOKE_2EModel &_soot2EModel);
	void AssignKineticScheme(OpenSMOKE_ReactingGas &_mix);
	void AssignEnergy(const bool iEnergy);
	void AssignConstantPressure();
	void AssignInletFlows(OpenSMOKE_GasStream &stream);
	void AssignVolume(const string units, const double value);
	void AssignEnd(const string units, const double value);
	void lock();

    // Option settings
	void SetName(const string name);
    void SetOptions(const string option);
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
	void SetViscosity(const double value, const string units);
	
	void SetUserDefinedHeatFlux(const string fileName);
	void SetUserDefinedExchangeArea(const string fileName);
	void SetUserDefinedHeatExchangeCoefficient(const string fileName);
	void SetUserDefinedAmbientTemperature(const string fileName);

	void SetConstantHeatFlux(const double value, const string units);
	void SetConstantExchangeArea(const double value, const string units);
	void SetConstantHeatExchangeCoefficient(const double value, const string units);
	void SetConstantAmbientTemperature(const double value, const string units);

	void SetReactionRatesOnFile(const vector<string> _names);
	void SetFormationRatesOnFile(const vector<string> _names);
	void SetROPAOnFile(const vector<string> _names);
	void SetSensitivityOnFile(const vector<string> _names);
	

    void UnsetGlobalKinetics();
    void UnsetTwoEquationModel();
    void UnsetVerbose();
    void UnsetHistory();
    void UnsetRelativeTolerance();
    void UnsetAbsoluteTolerance();
    void UnsetUserDefinedTemperature();
	void UnsetUserDefinedHeatFlux();
	void UnsetUserDefinedExchangeArea();
	void UnsetUserDefinedHeatExchangeCoefficient();
	void UnsetUserDefinedAmbientTemperature();
    void UnsetViscosity();
    void UnsetVerboseEnergy();

	// Messages for user
	void videoFinalResult();

	// Constructor
    void CreateReactor(OpenSMOKE_ReactingGas &mix, OpenSMOKE_GasStream &inletStream);
    void DefineFromFile(const string inputFile);

	// Solve Batch Reactor
	void Solve();

	// Video Summary
	void videoSummary();
	void OutletStream(OpenSMOKE_GasStream &outlet);
	void EnergyAnalysis(OpenSMOKE_GasStream &outlet);
	void MassAnalysis(OpenSMOKE_GasStream &outletStream);
	void SummaryOnFile(const string file_name);


private:

    // ODE System
	MyOdeSystem_Batch ODE_Batch_Object;

	// Error and warning messages
	string  name_object;
	void	ErrorMessage(const string message);
	void	WarningMessage(const string message);

	// Message functions
	void VideoGeneralInfo();
	void LabelODEFile();

	// Control Variables
	bool assignedKineticScheme;
	bool assignedGlobalKineticScheme;
	bool assignedSoot2EModel;
	bool assignedEnergy;
	bool assignedInletFlows;
	bool assignedVolume;
	bool assignedEnd;
	

	// Options
	bool iEnergy;
	bool iGlobalKinetics;
	bool iTwoEquationModel;
	bool iConstantPressure;
    bool iHistory;
    bool iVerbose;
	bool iVerboseEnergy;
	bool iVerboseReactionRates;
	bool iVerboseFormationRates;
	bool iVerboseROPA;
	bool iVerboseSensitivity;
    bool iRelativeTolerance;
	bool iAbsoluteTolerance;
	bool iUserDefinedTemperature;
	bool iUserDefinedViscosity;
	
    ProfileKind iUserDefinedHeatFlux;
    ProfileKind iUserDefinedExchangeArea;
	ProfileKind iUserDefinedHeatExchangeCoefficient;
	ProfileKind iUserDefinedAmbientTemperature;

	// Properties evaluation
	int     indexProperties;
	int     indexTwoEquations;

	// Print Options
	string outputFolderName;
	string outputBatchName;
	string outputSootName;
	string outputPAHName;
	string output2EName;
	string outputEnergyName;
	string outputReactionRatesName;
	string outputFormationRatesName;
	string outputROPAName;
	string outputSensitivityName;


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
	OpenSMOKE_ReactingGas		*mix;						// Gas Mixture
	OpenSMOKE_GlobalKinetics	*global;					// Global Kinetics
	OpenSMOKE_2EModel			*soot2EModel;				// Two Equation Model
	OpenSMOKE_GasStream			*inletStream;				// Inlet stream
	
	OpenSMOKE_UD_Profile	ud_temperature_profile;	// UD temperature profile
	OpenSMOKE_UD_Profile	ud_Qe_profile;			// UD Qe profile
	OpenSMOKE_UD_Profile	ud_A_profile;			// UD A profile
	OpenSMOKE_UD_Profile	ud_U_profile;			// UD Qe profile
	OpenSMOKE_UD_Profile	ud_Tambient_profile;	// UD Tambient profile

	// Kinetic scheme data
	int NC;							// Number of Species
	int NR;							// Number of Reactions

	// Unknowns
	BzzVectorDouble		omega;		    // Mass fractions [-]
	double              T;				// Reactor temperature [K]

    // Residuals
    BzzVectorDouble domega;
	double          dT;

	// Additional data
	double volume;					// Reactor volume	[m3]
	double P;				        // Reactor pressure [Pa]
	double mass;					// Reactor mass		[kg]
	double moles;					// Reactor moles	[kmol]

    // Main data
    
	double rho;						// Reactor density [kg/m3]
	double MWtot;					// Reactor molecular weight [kg/kmol]
	double Cp;						// Reactor specific heat [J/kg/K]
	double QReaction;				// Reactor reaction heat [J/m3/s]
	BzzVectorDouble		x;			// Mole fractions [-]
	BzzVectorDouble		c;			// Concentrations [kmol/m3]
	BzzVectorDouble		R;			// Formation Rates [kg/m3.s]

    // Reactor geometry
	double TauTotal;                // Reactor Contact Time [s]
    double Qe;						// External heat flux [W/m2]
    double A;						// Exchange surface [m2]
    double U;						// Heat exchange coefficient [W/m2/K]
    double Tambient;				// Ambient temperature [K]
    double Viscosity;				// Inlet stream viscosity [K]

	// Kinetic Maps
	BzzVectorDouble CpMap;
	BzzVectorDouble k1Map;
	BzzVectorDouble k2Map;
	BzzVectorDouble uKeqMap;
	BzzVectorDouble logFcentMap;
	BzzVectorDouble reactionDHMap;
	BzzVectorDouble reactionDSMap;

    // Output files
	ofstream    fPAH;
	ofstream    fSoot;
    ofstream    fSoot2E;
    ofstream    fOutput;
    ofstream    fEnergy;
	ofstream	fReactionRates;
	ofstream	fFormationRates;
	ofstream	fROPA;
	ofstream	fSensitivity;

    // Updating properties at each ODE call
    void    Initialize();
    void    UpdateHeatFlux(const double Tau);
	void    UpdateExchangeArea(const double Tau);
    void    UpdateProperties(const int indexJacobian, const int indexT);
	void    UpdateProperties_isothermal(int indexMemo);
	void    UpdateTwoEquationModel(BzzVectorDouble &y, BzzVectorDouble &dy);

    // Additional variables for momentum equation solution
    double  SimplifiedViscosity();

	// Additional variables 
	BzzVectorInt index_reaction_rates;
	BzzVectorInt index_formation_rates;
	BzzVectorInt index_ROPA;
	BzzVectorInt index_sensitivity;
	vector<string> names_ROPA;
	vector<string> names_reaction_rates;
	vector<string> names_formation_rates;
	vector<string> names_sensitivity;


public:

	// Pseudo-Private Functions
	void ODESystem_Isothermal_Batch(BzzVectorDouble &x, double t, BzzVectorDouble &f);
	void ODESystem_NonIsothermal_ConstantPressure_Batch(BzzVectorDouble &x, double t, BzzVectorDouble &f);
	void ODESystem_NonIsothermal_ConstantVolume_Batch(BzzVectorDouble &x, double t, BzzVectorDouble &f);
	void ODEPrint(BzzVectorDouble &omega, double t);


public:

    BzzVectorDouble Tau_History;
    BzzVectorDouble T_History;
    BzzMatrixDouble mass_History;
    BzzMatrixDouble mole_History;

    int countGlobalIterations;

private:

    static const double   ONE;
    static const double   ZERO;
    static const double   MAXNUMBER;
    static const int      MAX_TIME_STEPS;
    static const int      MINUSONE;
};

#endif // OPENSMOKE_BATCH_H
