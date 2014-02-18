/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci							   *
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
#include <string>
#include <sstream>
#include "idealreactors/batch/OpenSMOKE_Batch.h"
#include "basic/OpenSMOKE_Conversions.h"


OpenSMOKE_Batch	*ptBatch;

const double    OpenSMOKE_Batch::ONE             =  1.0;
const double    OpenSMOKE_Batch::ZERO            =  0.0;
const double    OpenSMOKE_Batch::MAXNUMBER       =  1.e16;
const int       OpenSMOKE_Batch::MAX_TIME_STEPS  =  50000;
const int       OpenSMOKE_Batch::MINUSONE        = -1;

OpenSMOKE_Batch::OpenSMOKE_Batch()
{
    ptBatch   =  this;

    // Control Variables
    assignedKineticScheme		= false;
    assignedGlobalKineticScheme	= false;
    assignedSoot2EModel			= false;
	assignedInletFlows			= false;
    assignedEnd					= false;
    assignedVolume				= false;
    assignedEnergy				= false;

    // Properties Index
    iGlobalKinetics             = false;    // Global Kinetics 0=OFF 1=ON
    iEnergy		                = false;    // 0 = Isothermal 1 = Energy Equation
	iConstantPressure			= false;	// 0 = Constant volume (default) 1=Constant Pressure
    iTwoEquationModel	        = false;    // 2E Model for soot predictions
    iRelativeTolerance	        = false;    // User defined relative tolerance
    iAbsoluteTolerance	        = false;    // User defined absolute tolerance
    iHistory                    = false;    // Save solution history
    iVerbose                    = true;     // Write on video and file
    iVerboseEnergy              = false;    // Report on energy
    iUserDefinedTemperature     = false;    // User defined temperature
	iUserDefinedViscosity		= false;	// User defined viscosity

	iVerboseReactionRates		= false;	// Reaction rates
	iVerboseFormationRates		= false;	// Formation rates
	iVerboseROPA				= false;	// Rate of production analysis
	iVerboseSensitivity			= false;	// Sensitivity analysis

	iUserDefinedHeatFlux				= NONE; 
	iUserDefinedExchangeArea			= NONE; 
	iUserDefinedHeatExchangeCoefficient	= NONE; 
	iUserDefinedAmbientTemperature		= NONE; 

    // Internal variables
    countGlobalIterations       =-1;    // Global index

    // Properties Index (It means that the Jacobian must be calculated)
    indexProperties			    = 0;

    // Default options
    nFileSteps		=   countFileSteps		= 10;
    nVideoSteps		=   countVideoSteps		= 10;
    countIterations = 0;

    // Output Folder
    outputFolderName			= "Output";
    outputBatchName				=  outputFolderName + "/" + "Batch.out";
    outputSootName				=  outputFolderName + "/" + "Soot.out";
    outputPAHName				=  outputFolderName + "/" + "PAH.out";
    output2EName				=  outputFolderName + "/" + "2E.out";
    outputEnergyName			=  outputFolderName + "/" + "Energy.out";
    outputReactionRatesName		=  outputFolderName + "/" + "ReactionRates.out";
    outputFormationRatesName	=  outputFolderName + "/" + "FormationRates.out";
    outputROPAName				=  outputFolderName + "/" + "ROPA.out";
    outputSensitivityName		=  outputFolderName + "/" + "Sensitivity.out";

	// Name Object
	name_object = "[Name not assigned]";

    // Default Values
    Qe						= 0.0;
	Viscosity				= 0.0;
}

void OpenSMOKE_Batch::AssignKineticScheme(OpenSMOKE_ReactingGas &_mix)
{
    ptBatch		=  this;
    mix			= &_mix;

    NC = mix->NumberOfSpecies();
    NR = mix->NumberOfReactions();

    // These vectors are used for the evaluation of gas mixture properties
    ChangeDimensions(NC, &omega);
    ChangeDimensions(NC, &domega);
    ChangeDimensions(NC, &x);
    ChangeDimensions(NC, &c);
    ChangeDimensions(NC, &R);

    ChangeDimensions(NC, &CpMap);
    ChangeDimensions(NR, &k1Map);
    ChangeDimensions(NR, &k2Map);
    ChangeDimensions(NR, &uKeqMap);
    ChangeDimensions(NR, &logFcentMap);
    ChangeDimensions(NR, &reactionDHMap);
    ChangeDimensions(NR, &reactionDSMap);

    // Control Variables
    assignedKineticScheme = true;
}

void OpenSMOKE_Batch::AssignInletFlows(OpenSMOKE_GasStream &_inletStream)
{
    if (assignedKineticScheme == false)
        ErrorMessage("You must define the kinetic scheme before defining the initial composition!");

    inletStream = &_inletStream;

    assignedInletFlows = true;
}

void OpenSMOKE_Batch::AssignSoot2EModel(OpenSMOKE_2EModel &_soot2EModel)
{
    soot2EModel = &_soot2EModel;

    assignedSoot2EModel = true;
}

void OpenSMOKE_Batch::AssignGlobalKineticScheme(OpenSMOKE_GlobalKinetics &_global)
{
    global = &_global;

    assignedGlobalKineticScheme = true;
}

void OpenSMOKE_Batch::AssignEnergy(const bool _iEnergy)
{
	iEnergy = _iEnergy;

    assignedEnergy = true;
}

void OpenSMOKE_Batch::AssignConstantPressure()
{
	iConstantPressure = true;
}

void OpenSMOKE_Batch::AssignEnd(const string units, const double value)
{
	TauTotal	= OpenSMOKE_Conversions::conversion_time(value, units);
    assignedEnd = true;
}

void OpenSMOKE_Batch::AssignVolume(const string units, const double value)
{
	volume			= OpenSMOKE_Conversions::conversion_volume(value, units);
    assignedVolume	= true;
}

void OpenSMOKE_Batch::SetName(const string name)
{
    name_object = name;
}

void OpenSMOKE_Batch::SetGlobalKinetics()
{
    iGlobalKinetics = true;
}

void OpenSMOKE_Batch::UnsetGlobalKinetics()
{
    iGlobalKinetics = false;
}

void OpenSMOKE_Batch::SetTwoEquationModel()
{
    iTwoEquationModel = true;
}

void OpenSMOKE_Batch::UnsetTwoEquationModel()
{
    iTwoEquationModel = false;
}

void OpenSMOKE_Batch::SetVerbose()
{
    iVerbose = true;
}

void OpenSMOKE_Batch::UnsetVerbose()
{
    iVerbose = false;
}

void OpenSMOKE_Batch::SetRelativeTolerance(const double relative)
{
    iRelativeTolerance	= true;
    relativeTolerance	= relative;
}

void OpenSMOKE_Batch::UnsetRelativeTolerance()
{
    iRelativeTolerance = false;
}

void OpenSMOKE_Batch::SetAbsoluteTolerance(const double absolute)
{
    iAbsoluteTolerance = true;
    absoluteTolerance = absolute;
}

void OpenSMOKE_Batch::UnsetAbsoluteTolerance()
{
    iAbsoluteTolerance = false;
}

void OpenSMOKE_Batch::SetUserDefinedTemperature(const string fileName)
{
    iUserDefinedTemperature = true;
    ud_temperature_profile.AssignFromFile(fileName, "TEMPERATURE");
	ud_temperature_profile.SetName(name_object + " - Temperature Profile");
}

void OpenSMOKE_Batch::UnsetUserDefinedTemperature()
{
    iUserDefinedTemperature = false;
}

void OpenSMOKE_Batch::SetUserDefinedHeatFlux(const string fileName)
{
    iUserDefinedHeatFlux = USERDEFINED;
    ud_Qe_profile.AssignFromFile(fileName, "HEAT_FLUX");
	ud_Qe_profile.SetName(name_object + " - Heat Flux Profile");
}

void OpenSMOKE_Batch::UnsetUserDefinedHeatFlux()
{
    iUserDefinedHeatFlux = NONE;
}

void OpenSMOKE_Batch::SetUserDefinedExchangeArea(const string fileName)
{
    iUserDefinedExchangeArea = USERDEFINED;
    ud_A_profile.AssignFromFile(fileName, "AREA");
	ud_A_profile.SetName(name_object + " - Exchange Area Profile");
}

void OpenSMOKE_Batch::UnsetUserDefinedExchangeArea()
{
    iUserDefinedExchangeArea = NONE;
}

void OpenSMOKE_Batch::SetUserDefinedHeatExchangeCoefficient(const string fileName)
{
    iUserDefinedHeatExchangeCoefficient = USERDEFINED;
    ud_U_profile.AssignFromFile(fileName, "HEAT_TRANSFER_COEFFICIENT");
	ud_U_profile.SetName(name_object + " - Heat Exchange Coefficient Profile");
}

void OpenSMOKE_Batch::UnsetUserDefinedHeatExchangeCoefficient()
{
    iUserDefinedHeatExchangeCoefficient = NONE;
}

void OpenSMOKE_Batch::SetUserDefinedAmbientTemperature(const string fileName)
{
    iUserDefinedAmbientTemperature = USERDEFINED;
    ud_Tambient_profile.AssignFromFile(fileName, "TEMPERATURE");
	ud_Tambient_profile.SetName(name_object + " - Ambient Temperature Profile");
}

void OpenSMOKE_Batch::UnsetUserDefinedAmbientTemperature()
{
    iUserDefinedAmbientTemperature = NONE;
}

void OpenSMOKE_Batch::SetVerboseEnergy()
{
    iVerboseEnergy = true;
}

void OpenSMOKE_Batch::UnsetVerboseEnergy()
{
    iVerboseEnergy = false;
}

void OpenSMOKE_Batch::SetHistory()
{
    iHistory = true;

    ChangeDimensions(MAX_TIME_STEPS, &Tau_History);
    ChangeDimensions(MAX_TIME_STEPS, &T_History);
    ChangeDimensions(MAX_TIME_STEPS, mix->NumberOfSpecies(), &mass_History);
    ChangeDimensions(MAX_TIME_STEPS, mix->NumberOfSpecies(), &mole_History);
}

void OpenSMOKE_Batch::UnsetHistory()
{
    iHistory = false;
}

void OpenSMOKE_Batch::SetConstantHeatFlux(const double value, const string units)
{
	iUserDefinedHeatFlux = CONSTANT;
	Qe = OpenSMOKE_Conversions::conversion_heat_flux(value, units);
}

void OpenSMOKE_Batch::SetConstantExchangeArea(const double value, const string units)
{
	iUserDefinedExchangeArea = CONSTANT;
	A = OpenSMOKE_Conversions::conversion_area(value, units);
}

void OpenSMOKE_Batch::SetConstantHeatExchangeCoefficient(const double value, const string units)
{
	iUserDefinedHeatExchangeCoefficient = CONSTANT;
	U = OpenSMOKE_Conversions::conversion_heat_exchange_coefficient(value, units);
}

void OpenSMOKE_Batch::SetConstantAmbientTemperature(const double value, const string units)
{
	iUserDefinedAmbientTemperature = CONSTANT;
	Tambient = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_Batch::SetViscosity(const double value, const string units)
{
	iUserDefinedViscosity = true;
	Viscosity = OpenSMOKE_Conversions::conversion_dynamic_viscosity(value, units);
}

void OpenSMOKE_Batch::UnsetViscosity()
{
	iUserDefinedViscosity = false;
	Viscosity = 0.;
}

void OpenSMOKE_Batch::SetOptions(const string option)
{
    if (option == "NOTHING")
    {
    }
    else
        ErrorMessage("This option is not allowed: " + option);
}

void OpenSMOKE_Batch::SetOutputFolder(const string _outputFolderName)
{
    outputFolderName			= _outputFolderName;
    outputBatchName				=  outputFolderName + "/" + "Batch.out";
    outputSootName				=  outputFolderName + "/" + "Soot.out";
    outputPAHName				=  outputFolderName + "/" + "PAH.out";
    output2EName				=  outputFolderName + "/" + "2E.out";
    outputEnergyName			=  outputFolderName + "/" + "Energy.out";
	outputReactionRatesName		=  outputFolderName + "/" + "ReactionRates.out";
    outputFormationRatesName	=  outputFolderName + "/" + "FormationRates.out";
    outputROPAName				=  outputFolderName + "/" + "ROPA.out";
    outputSensitivityName		=  outputFolderName + "/" + "Sensitivity.out";
}

void OpenSMOKE_Batch::SetVideoOptions(const int _nVideoSteps)
{
    nVideoSteps		=   countVideoSteps		= _nVideoSteps;
                        countIterations     = 0;
}

void OpenSMOKE_Batch::SetFileOptions(const int _nFileSteps)
{
    nFileSteps		=   countFileSteps		= _nFileSteps;
                        countIterations     = 0;
}

void OpenSMOKE_Batch::SetReactionRatesOnFile(const vector<string> _names)
{
	if (_names[0] == "ALL")
	{
		ChangeDimensions(NR, &index_reaction_rates);
		for(int k=1;k<=NR;k++)
		{
			index_reaction_rates[k] = k;
			names_reaction_rates.push_back(mix->strReaction[index_reaction_rates[k]]);
		}

	}
	else
	{
		ChangeDimensions(_names.size(), &index_reaction_rates);
		for(int k=0;k<_names.size();k++)
		{
			stringstream index_string;
			index_string << _names[k];
			index_string >> index_reaction_rates[k+1];
			if (index_reaction_rates[k+1] > NR)
				ErrorMessage("The requested reaction for the reaction rates analysis is not available: " + _names[k]);
			names_reaction_rates.push_back(mix->strReaction[index_reaction_rates[k+1]]);
		}
	}

	iVerboseReactionRates = true;
}

void OpenSMOKE_Batch::SetFormationRatesOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_formation_rates, names_formation_rates);
	iVerboseFormationRates = true;
}

void OpenSMOKE_Batch::SetROPAOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_ROPA, names_ROPA);
	iVerboseROPA = true;
}

void OpenSMOKE_Batch::SetSensitivityOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_sensitivity, names_sensitivity);
	iVerboseSensitivity = true;
}

void OpenSMOKE_Batch::lock()
{
    if (assignedKineticScheme == false)
        ErrorMessage("The kinetic scheme was not defined!");
    if (assignedGlobalKineticScheme == false)
        ErrorMessage("The global kinetic scheme was not defined!");
    if (assignedEnd == false)
		ErrorMessage("The reactor contact time was not defined!");
    if (assignedInletFlows == false)
        ErrorMessage("The inlet flow was not defined!");
    if (assignedEnergy == false)
        ErrorMessage("The kind of reactor was not defined!");
	if (assignedSoot2EModel == false)
        ErrorMessage("The 2E Model was not defined!");

	Initialize();
}

void OpenSMOKE_Batch::Initialize()
{
    T            = inletStream->T;
	P            = inletStream->P;
    omega        = inletStream->omega;
	mass         = inletStream->rho*volume;

    if (iUserDefinedTemperature == true)
        ud_temperature_profile.Check(0., T);

    if (iEnergy == false)					// Isothermal reactor
    {
        UpdateProperties_isothermal(0);
        indexProperties = 1;
    }
    else if (iEnergy == true)				// Non-isothermal reactor
        UpdateProperties(MINUSONE, NC+1);
}

void OpenSMOKE_Batch::DefineFromFile(const string inputFile)
{
    double			double_value;
    string			string_value;
    int				int_value;
	vector<string>  string_vector;

    OpenSMOKE_Dictionary_Batch dictionary;
    dictionary.ParseFile(inputFile);


    // COMPULSORY: Energy Equation solution
	AssignEnergy(dictionary.Return("#Energy"));

    // COMPULSORY: Reactor Contact Time
	if (dictionary.Return("#ContactTime", double_value, string_value))
		AssignEnd(string_value, double_value);

    // COMPULSORY: Reactor Volume
    if (dictionary.Return("#Volume", double_value, string_value))
        AssignVolume(string_value, double_value);

	// COMPULSORY: Constant Volume
	if (dictionary.Return("#ConstantVolume"))
	{}

	// COMPULSORY: Constant Volume
	if (dictionary.Return("#ConstantPressure"))
		AssignConstantPressure();


    // OPTIONAL: Output folder
    if (dictionary.Return("#Output", string_value))
        SetOutputFolder(string_value);

    // OPTIONAL: Print Options
    if (dictionary.Return("#Nvideo", int_value))
        SetVideoOptions(int_value);

    // OPTIONAL: Print Options
    if (dictionary.Return("#Nfile", int_value))
        SetFileOptions(int_value);

    // OPTIONAL: Print Options
    if (dictionary.Return("#NoVerbose"))
        UnsetVerbose();

    // OPTIONAL: Global kinetics
    if (dictionary.Return("#GlobalKinetics"))
        SetGlobalKinetics();

    // OPTIONAL: Relative tolerance
    if (dictionary.Return("#RelTol", double_value))
        SetRelativeTolerance(double_value);

    // OPTIONAL: Absolute tolerance
    if (dictionary.Return("#AbsTol", double_value))
        SetAbsoluteTolerance(double_value);

    // OPTIONAL: UserDefinedTemperature
    if (dictionary.Return("#UserDefined_T", string_value))
        SetUserDefinedTemperature(string_value);

	// OPTIONAL: Verbose Energy
    if (dictionary.Return("#VerboseEnergy"))
        SetVerboseEnergy();

	// OPTIONAL: Heat flux
    if (dictionary.Return("#Qe", double_value, string_value))
        SetConstantHeatFlux(double_value, string_value);

	// OPTIONAL: Exchange area
    if (dictionary.Return("#A", double_value, string_value))
        SetConstantExchangeArea(double_value, string_value);

	// OPTIONAL: Exchange coefficient
    if (dictionary.Return("#U", double_value, string_value))
        SetConstantHeatExchangeCoefficient(double_value, string_value);

	// OPTIONAL: Ambient temperature
    if (dictionary.Return("#Tambient", double_value, string_value))
        SetConstantAmbientTemperature(double_value, string_value);


	// OPTIONAL: UserDefined_Qe
    if (dictionary.Return("#UserDefined_Qe", string_value))
        SetUserDefinedHeatFlux(string_value);

	// OPTIONAL: UserDefined_Ae
    if (dictionary.Return("#UserDefined_A", string_value))
        SetUserDefinedExchangeArea(string_value);

	// OPTIONAL: UserDefined_U
    if (dictionary.Return("#UserDefined_U", string_value))
        SetUserDefinedHeatExchangeCoefficient(string_value);

	// OPTIONAL: UserDefined_Tambient
    if (dictionary.Return("#UserDefined_Tambient", string_value))
        SetUserDefinedAmbientTemperature(string_value);
 
	// OPTIONAL: Inlet Viscosity
    if (dictionary.Return("#Viscosity", double_value, string_value))
        SetViscosity(double_value, string_value);

	// OPTIONAL: TwoEquationModel
    if (dictionary.Return("#2E_Model"))
        SetTwoEquationModel();

	// OPTIONAL: Reaction Rates
    if (dictionary.Return("#ReactionRates", string_vector))
		SetReactionRatesOnFile(string_vector);

	// OPTIONAL: Formation Rates
    if (dictionary.Return("#FormationRates", string_vector))
		SetFormationRatesOnFile(string_vector);

	// OPTIONAL: Rate of Production Analysis
    if (dictionary.Return("#ROPA", string_vector))
		SetROPAOnFile(string_vector);

	// OPTIONAL: Sensitivity Analysis
    if (dictionary.Return("#Sensitivity", string_vector))
		SetSensitivityOnFile(string_vector);

	lock();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - ISOTHERMAL REACTOR								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenSMOKE_Batch::ODESystem_Isothermal_Batch(BzzVectorDouble &y, double t, BzzVectorDouble &dy)
{
	// Recover variables
    omega    = y;

    // User Defined Temperature profile
    if (iUserDefinedTemperature == true)
    {
        T = ud_temperature_profile.GiveMeValue(t, MINUSONE);
        indexProperties = 0;
    }

    // Properties + Reaction rates (constant temperature)
    UpdateProperties_isothermal(indexProperties);


    // Conservation equations for all the species (Batch)
    domega = R / rho;


    // Soot Two Equation Model
    if (iTwoEquationModel == true)
        UpdateTwoEquationModel(y, dy);


    // Recovering residuals (Soot 2E are already updated)
    dy   = domega;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - NON ISOTHERMAL REACTOR - CONST P					   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenSMOKE_Batch::ODESystem_NonIsothermal_ConstantPressure_Batch(BzzVectorDouble &y, double t, BzzVectorDouble &dy)
{
    int i;

    // Recovering variables;
    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    T				                = y[NC+1];

	// Updating Exchange area
	UpdateExchangeArea(t);
    UpdateHeatFlux(t);

    // Evaluation of properties
    UpdateProperties(ODE_Batch_Object.jacobianIndex, NC+1);

    // Conservation equations for all the species (Batch)
    domega	= R / rho ;											// [1/s]

    // Conservation equations for the energy (Batch)
    dT  = (volume*QReaction - A*Qe) / (mass*Cp);				// [K/s]
 
    // Soot Two Equation Model
    if (iTwoEquationModel == true)	UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]	= domega[i];
    dy[NC+1]					= dT;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//						ODE SYSTEM EQUATIONS - NON ISOTHERMAL REACTOR - CONST V					   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenSMOKE_Batch::ODESystem_NonIsothermal_ConstantVolume_Batch(BzzVectorDouble &y, double t, BzzVectorDouble &dy)
{
    int i;

    // Recovering variables;
    for (i=1;i<=NC;i++)	omega[i]    = y[i];
    T				                = y[NC+1];

	// Updating Exchange area
	UpdateExchangeArea(t);
    UpdateHeatFlux(t);

    // Evaluation of properties
    UpdateProperties(ODE_Batch_Object.jacobianIndex, NC+1);

    // Conservation equations for all the species (Batch)
    domega	= R / rho ;

	// Number of moles
	double sumRmole = 0.;
	for (i=1;i<=NC;i++)
		sumRmole += R[i]/mix->M[i];

    // Conservation equations for the energy (Batch)
    dT  = (volume*(QReaction + Constants::R_J_kmol*T*sumRmole) - A*Qe) / (mass * (Cp-Constants::R_J_kmol/MWtot) );
 
    // Soot Two Equation Model
    if (iTwoEquationModel == true)	UpdateTwoEquationModel(y, dy);

    // Recovering residuals (Soot 2E are already updated)
    for (i=1;i<=NC;i++)	dy[i]	= domega[i];
    dy[NC+1]					= dT;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									GAS MIXTURE PROPERTIES										   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::UpdateProperties_isothermal(int memoIndex)
{
    double cTot;

    mix->GetMWAndMoleFractionsFromMassFractions(MWtot, x, omega);

    // --------------------------------------------------------------------------
    // PROPERTIES FOR CONSTANT T + MEMORIZATION
    // --------------------------------------------------------------------------
    if (memoIndex==0)
    {
        // Kinetic parameters (the enthalpies are calculated in this step)
        // ----------------------------------------------------------------------
        mix->ComputeKineticParameters( T, log(T), 1./T);

        memoIndex = 1;
    }

    // --------------------------------------------------------------------------
    // PROPERTIES FOR CONSTANT T (Recover from tables)
    // --------------------------------------------------------------------------
    if (memoIndex == 1)
    {
        // a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		if (iConstantPressure == false)	// Constant volume
		{
			moles	= mass / MWtot;						// [kmol]
			rho		= mass / volume;					// [kg/m3]
			cTot	= rho  / MWtot;						// [kmol/m3]
			P		= cTot * (Constants::R_J_kmol*T);	// [Pa]
			c       = cTot * x;							// [kmol/m3]
			
		}
		else
		{
			moles	= mass  / MWtot;					// [kmol]
			volume	= moles * Constants::R_J_kmol*T/P;	// [m3]
			rho		= mass  / volume;					// [kg/m3]
			cTot	= rho   / MWtot;					// [kmol/m3]
			c       = cTot  * x;						// [kmol/m3]
		}			

        // c. Reaction Rates
        mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
        ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
    }

    // Global Kinetics
    if (iGlobalKinetics==true)
        global->GiveMeFormationRates(T,c,R);					// [kg/m3/s]
}

void OpenSMOKE_Batch::UpdateProperties(const int jacobianIndex, const int indexT)
{
	int memoIndex;

	if		(jacobianIndex<0)			memoIndex = -2;		// Calculate all without storing
	else if (jacobianIndex==0)			memoIndex = -1;		// Calculate and storing
	else if (jacobianIndex >= indexT)	memoIndex = -2;		// Calculate all without storing
	else								memoIndex = 0;		// Recover from maps

	double cTot;

	// Mole fractions
	// ----------------------------------------------------------------------
	mix->GetMWAndMoleFractionsFromMassFractions(MWtot, x, omega);

	// --------------------------------------------------------------------------
	// Update every variable which does not depend on T
	// --------------------------------------------------------------------------
	if(memoIndex == -2)
	{
        // a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		if (iConstantPressure == false)	// Constant volume
		{
			moles	= mass / MWtot;						// [kmol]
			rho		= mass / volume;					// [kg/m3]
			cTot	= rho  / MWtot;						// [kmol/m3]
			P		= cTot * (Constants::R_J_kmol*T);	// [Pa]
			c       = cTot * x;							// [kmol/m3]
			
		}
		else
		{
			moles	= mass  / MWtot;					// [kmol]
			volume	= moles * Constants::R_J_kmol*T/P;	// [m3]
			rho		= mass  / volume;					// [kg/m3]
			cTot	= rho   / MWtot;					// [kmol/m3]
			c       = cTot  * x;						// [kmol/m3]
		}

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		Cp = mix->MixCp_FromMassFractions(omega);

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters( T, log(T), 1./T);
		mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
		ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);					// [J/m3.s]
	}

	// --------------------------------------------------------------------------
	// Update only the variables depending on the temperature
	// --------------------------------------------------------------------------
	else if(memoIndex==-1)
	{
		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->SpeciesCp(T);
		CpMap = mix->Cp;

		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->ComputeKineticParameters( T, log(T), 1./T);
		k1Map			= mix->k1;
		k2Map			= mix->k2;
		uKeqMap			= mix->uKeq;
		logFcentMap		= mix->logFcent;
		reactionDHMap	= mix->reactionDH;
		reactionDSMap	= mix->reactionDS;

		memoIndex = 0;
	}

	// --------------------------------------------------------------------------
	// Update every variable which does not depend on T
	// --------------------------------------------------------------------------
	if(memoIndex == 0)
	{
        // a. Calcolo della concentrazione e della densita [kmol/m3] [kg/m3]
		if (iConstantPressure == false)	// Constant volume
		{
			moles	= mass / MWtot;						// [kmol]
			rho		= mass / volume;					// [kg/m3]
			cTot	= rho  / MWtot;						// [kmol/m3]
			P		= cTot * (Constants::R_J_kmol*T);	// [Pa]
			c       = cTot * x;							// [kmol/m3]
			
		}
		else							// Constant Pressure
		{
			moles	= mass  / MWtot;					// [kmol]
			volume	= moles * Constants::R_J_kmol*T/P;	// [m3]
			rho		= mass  / volume;					// [kg/m3]
			cTot	= rho   / MWtot;					// [kmol/m3]
			c       = cTot  * x;						// [kmol/m3]
		}

		// Specific Heat [J/kgK]
		// ----------------------------------------------------------------------
		mix->Cp = CpMap;
		Cp = mix->MixCp_FromMassFractions(omega);


		// Table Creation (Memorization)
		// ----------------------------------------------------------------------
		mix->k1			= k1Map;
		mix->k2			= k2Map;
		mix->uKeq		= uKeqMap;
		mix->logFcent	= logFcentMap;
		mix->reactionDH	= reactionDHMap;
		mix->reactionDS	= reactionDSMap;
		mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
		ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
		QReaction = - mix->ComputeQReaction(T);					// [J/m3.s]
	}

	// Global Kinetics
	if (iGlobalKinetics == true)
	{
		global->GiveMeFormationRates(T,c,R);
		global->GiveMeReactionHeat(T, R, QReaction);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//										PRINT ON FILE											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void ODEPrintExternal(BzzVectorDouble &y, double t)
{
    ptBatch->ODEPrint(y, t);
}

void OpenSMOKE_Batch::ODEPrint(BzzVectorDouble &y, double t)
{
    if (iVerbose == true)
    {
        if ( countIterations%(300*nVideoSteps) == 0)
        {
            cout    << endl
                    << "#"			<< "\t"
                    << "Tau[s]"		<< "\t\t"
                    << "T[K]"       << "\t\t"
					<< "P[atm]"		<< "\t\t"
					<< "V[l]"		<< "\t\t"
					<< "n[kmol]"	<< "\t\t";
                    
			cout	<< endl;
        }

        if (countVideoSteps == nVideoSteps)
        {
            cout	<< countIterations	<< "\t"
                    << t	  			<< "\t"
                    << T				<< "\t"
                    << P/101325.		<< "\t"
                    << volume/1000.		<< "\t"
                    << moles			<< "\t";

            // Soot Two Equation Model
            if (iTwoEquationModel == true)
            {
                cout 	<< soot2EModel->m0   << "\t"
                        << soot2EModel->fv	<< "\t";
            }

            cout	<< endl;

            countVideoSteps = 0;
        }

        if (countFileSteps == nFileSteps)
        {
            int i;

            fOutput << t					<< "\t"
                    << T					<< "\t"
                    << P/101325.			<< "\t"
                    << volume/1000.			<< "\t"
                    << QReaction			<< "\t"
                    << Qe*A/volume     		<< "\t"
                    << rho					<< "\t"
                    << MWtot				<< "\t";

            for (i=1;i<=NC;i++)
                fOutput << x[i]				<< "\t";
            for (i=1;i<=NC;i++)
                fOutput << omega[i]			<< "\t";
            fOutput << endl;

            // Soot Two Equation Model
            if (iTwoEquationModel == true)
                soot2EModel->write_on_file(fSoot2E, t, y[indexTwoEquations], y[indexTwoEquations+1]);

            if (mix->soot_manager.iSoot == 1)
            {
                // Soot Distribution data
                mix->soot_manager.large_calculateAll(omega, x, c, rho);
                mix->soot_manager.small_calculateAll(omega, x, c, rho);

                fSoot << t		<< "\t";
                fSoot << T		<< "\t";
                mix->soot_manager.print_main_data_on_file(fSoot);
                fSoot << endl;
            }

            if (mix->pah_manager.iPAH == 1)
            {
                // PAH Distribution data
                mix->pah_manager.pah_calculateAll(omega, x, c, rho);

                fPAH << t		<< "\t";
                fPAH << T		<< "\t";
                mix->pah_manager.print_main_data_on_file(fPAH);
                fPAH << endl;
            }

			if (iVerboseEnergy == true)
			{
				BzzVectorDouble h_mass(NC);
				mix->GetMixAveragedEnthalpy_Mass(h_mass, T);
	
				double H_mass  = Dot(omega, h_mass);			// [J/kg]
				double Ep_mass = P/rho;							// [J/kg]
				double Ek_mass = 0.;							// [J/kg]

				fEnergy << t									<< "\t"
						<< H_mass+Ep_mass						<< "\t"
						<< H_mass								<< "\t"
						<< Ek_mass								<< "\t"
						<< Ep_mass								<< "\t"
						<< QReaction							<< "\t"
						<< Qe*A		       						<< "\t"
						<< Cp		       						<< "\t"
						<< endl;
			}

			if (iVerboseReactionRates == true)
			{
				fReactionRates	<< t	<< "\t"
								<< 0	<< "\t"
								<< T	<< "\t";
				
				for (i=1;i<=index_reaction_rates.Size();i++)
					fReactionRates	<< mix->r[index_reaction_rates[i]]	<< "\t";
				
				fReactionRates	<< endl;
			}

			if (iVerboseFormationRates == true)
			{
				fFormationRates	<< t	<< "\t"
								<< 0	<< "\t"
								<< T	<< "\t";
				
				for (i=1;i<=index_formation_rates.Size();i++)
					fFormationRates	<< R[index_formation_rates[i]]	<< "\t";
				
				fFormationRates	<< endl;
			}

			if (iVerboseROPA == true)
			{
			}

			if (iVerboseSensitivity == true)
			{
			}

            countFileSteps = 0;
        }

    }

    countVideoSteps++;
    countFileSteps++;
    countIterations++;

    if (iHistory == true)
    {
        countGlobalIterations++;

        if(countGlobalIterations>0)
        {
            Tau_History[countGlobalIterations] = t;
            T_History[countGlobalIterations]   = T;
            mass_History.SetRow(countGlobalIterations, omega);
            mole_History.SetRow(countGlobalIterations, x);
        }
    }
}

void OpenSMOKE_Batch::LabelODEFile()
{
    int i;

    if (iVerbose == true)
    {

        {
            fOutput << "Tau[s](1)    "		<< "\t";
            fOutput << "T[K](2)      "		<< "\t";
            fOutput << "P[atm](3)    "		<< "\t";
            fOutput << "V[l](4)      "		<< "\t";
            fOutput << "QR[W/m3](5)  "		<< "\t";
            fOutput << "QE[W/m3](6)  "		<< "\t";
            fOutput << "rho[kg/m3](7)"	    << "\t";
            fOutput << "MW[kg/kml](8)"	    << "\t";

            int countFileOutput = 9;
            for (i=1;i<=NC;i++)
                fOutput << mix->names[i] << " (x" << countFileOutput++ << ")  \t";
            for (i=1;i<=NC;i++)
                fOutput << mix->names[i] << " (w" << countFileOutput++ << ")  \t";

            fOutput << endl;
            fOutput << endl;
        }

        if (iTwoEquationModel == true)
            soot2EModel->write_label_on_file(fSoot2E, 1, "Tau[s]");

        if (mix->pah_manager.iPAH == 1)
		{
            fPAH << "time[s](1) " << "\t";
			fPAH << "T[K](2)    " << "\t";
			mix->pah_manager.GnuPlotInterface(fPAH, 3);
			fPAH << endl << endl;
		}

        if (mix->soot_manager.iSoot == 1)
        {
            fSoot << "time[s](1) " << "\t";
            fSoot << "T[K](2)    " << "\t";
            mix->soot_manager.GnuPlotInterface(fSoot, 3);
            fSoot << endl << endl;
        }
		
		if (iVerboseEnergy == true)
		{
			fEnergy << "Tau[s](1)      "		<< "\t";
			fEnergy << "H+Ep[J/kg](2)  "		<< "\t";
			fEnergy << "H[J/kg](3)     "		<< "\t";
			fEnergy << "Ek[J/kg](4)    "		<< "\t";
            fEnergy << "Ep[J/kg](5)    "		<< "\t";
            fEnergy << "QR[W/m3](6)    "		<< "\t";
            fEnergy << "QE[W/m3](7)    "		<< "\t";
            fEnergy << "Cp[J/kg/K](8)  "		<< "\t";
		}

		if (iVerboseReactionRates == true)
		{
			fReactionRates << "Tau[s](1)      "		<< "\t";
            fReactionRates << "dummy(2)       "		<< "\t";
            fReactionRates << "T[K](3)        "		<< "\t";
			
			int countFileOutput = 4;
            for (i=1;i<=index_reaction_rates.Size();i++)
                fReactionRates << "R_" << index_reaction_rates[i] << "[" << names_reaction_rates[i-1] << "](" << countFileOutput++ << ")  \t";
			fReactionRates << endl;
			fReactionRates << endl;
		}

		if (iVerboseFormationRates == true)
		{
			fReactionRates << "Tau[s](1)      "		<< "\t";
            fReactionRates << "dummy(2)       "		<< "\t";
            fReactionRates << "T[K](3)        "		<< "\t";
			
			int countFileOutput = 4;
            for (i=1;i<=index_formation_rates.Size();i++)
                fFormationRates << names_formation_rates[i-1] << "(" << countFileOutput++ << ")  \t";
			fFormationRates << endl;
			fFormationRates << endl;
		}
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									SOLVE Batch       											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::Solve()
{
    ptBatch = this;

    double timeStart, timeEnd;
    BzzVectorDouble xMin;
    BzzVectorDouble xMax;
    BzzVectorDouble xInitial;

    if (iVerbose == true)
    {
        // Create output folder
        string message = "mkdir " + outputFolderName;
        system(message.c_str());

        openOutputFileAndControl(fOutput, outputBatchName.c_str());
        fOutput.setf(ios::scientific);

        if (iTwoEquationModel==true)
        {
            openOutputFileAndControl(fSoot2E, output2EName.c_str());
            fSoot2E.setf(ios::scientific);
        }

        if (mix->pah_manager.iPAH == 1)
        {
            openOutputFileAndControl(fPAH, outputPAHName.c_str());
            fPAH.setf(ios::scientific);
        }

        if (mix->soot_manager.iSoot == 1)
        {
            openOutputFileAndControl(fSoot, outputSootName.c_str());
            fSoot.setf(ios::scientific);
        }

		if (iVerboseEnergy == true)
        {
            openOutputFileAndControl(fEnergy, outputEnergyName.c_str());
            fEnergy.setf(ios::scientific);
        }

		if (iVerboseReactionRates == true)
        {
            openOutputFileAndControl(fReactionRates, outputReactionRatesName);
            fReactionRates.setf(ios::scientific);
        }

		if (iVerboseFormationRates == true)
        {
            openOutputFileAndControl(fFormationRates, outputFormationRatesName);
            fFormationRates.setf(ios::scientific);
        }

		if (iVerboseROPA == true)
        {
            openOutputFileAndControl(fROPA, outputROPAName);
            fROPA.setf(ios::scientific);
        }

		if (iVerboseSensitivity == true)
        {
            openOutputFileAndControl(fSensitivity, outputSensitivityName);
            fSensitivity.setf(ios::scientific);
        }

        LabelODEFile();
    }


    ODE_Batch_Object.assignBatch(this, iEnergy, iConstantPressure);

    // 1.A Isothermal Batch
    if (iEnergy == false && iTwoEquationModel == false)
    {
        ChangeDimensions(NC, &xMin);  xMin=ZERO;					
        ChangeDimensions(NC, &xMax);  xMax=ONE;					// Mass fractions

        xInitial = omega;
    }

    // 1.C Isothermal Batch	+ Two Equation Model
    if (iEnergy == false && iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+1;

        ChangeDimensions(NC+2, &xMin);
        xMin=ZERO;
        
		ChangeDimensions(NC+2, &xMax);
        xMax=ONE;

        xInitial = omega;
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }

    // 2.A NonIsothermal Batch
    if (iEnergy == true && iTwoEquationModel == false)
    {
        ChangeDimensions(NC+1, &xMin);
        xMin=ZERO;
        xMin[NC+1] = mix->TMIN+1.;
        
		ChangeDimensions(NC+1, &xMax);
        xMax=ONE;
        xMax[NC+1] = mix->TMAX-1.;
        
        xInitial = omega;
        xInitial.Append(T);
    }

    // 2.C NonIsothermal Batch + TwoEquations
    if (iEnergy == true && iTwoEquationModel == true)
    {
        soot2EModel->initial_values(rho);
        indexTwoEquations = NC+2;

        ChangeDimensions(NC+3, &xMin);
        xMin=ZERO;
        xMin[NC+1] = mix->TMIN+1.;
        
		ChangeDimensions(NC+3, &xMax);
        xMax=ONE;
        xMax[NC+1] = mix->TMAX-1.;

        xInitial = omega;
        xInitial.Append(T);
        xInitial.Append(soot2EModel->phiNStart);
        xInitial.Append(soot2EModel->phiMStart);
    }


    BzzOdeDoubleStiffObject o(xInitial, 0., &ODE_Batch_Object);

    o.StepPrint(ODEPrintExternal);
    o.SetMinimumConstraints(xMin);
    o.SetMaximumConstraints(xMax);

    if (iRelativeTolerance == true)	o.SetTollRel(MachEps()*relativeTolerance);
    if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

    {
        countGlobalIterations = -1;  // From -1 to avoid to store results from the first iteration

        timeStart = BzzGetCpuTime();
		
		o(TauTotal,TauTotal);
		
		timeEnd = BzzGetCpuTime();
    }


    if (iVerbose == true)
    {
        cout << endl;
        cout << "Number of function for Jacobian: " << o.GetNumFunctionForJacobian()		<< endl;
        cout << "Numerical Jacobians: "				<< o.GetNumNumericalJacobian()			<< endl;
        cout << "Time DAE solution: "				<< timeEnd - timeStart	<< " s"			<< endl << endl;
    }

    if (iVerbose == true)
    {
        fOutput.close();

        if (mix->pah_manager.iPAH == 1)		fPAH.close();
        if (mix->soot_manager.iSoot == 1)	fSoot.close();
		if (iVerboseEnergy == true)			fEnergy.close();
		if (iVerboseReactionRates == true)	fReactionRates.close();
		if (iVerboseFormationRates == true)	fFormationRates.close();
		if (iVerboseROPA == true)			fROPA.close();
		if (iVerboseSensitivity == true)	fSensitivity.close();
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									VIDEO INFORMATION											   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void OpenSMOKE_Batch::VideoGeneralInfo()
{
    cout.setf(ios::scientific);
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;

    cout << "---------------------------------------------------------------------" << endl;
    cout << "Batch Summary" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    cout << "Contact time:\t\t"			<< TauTotal							<< " [s]"		<< endl;
	cout << "Pressure:\t\t"				<< P								<< " [Pa]"		<< endl;
    cout << "Density:\t\t"				<< rho								<< " [kg/m3]"	<< endl;
    cout << "Molecular weight:\t"		<< MWtot							<< " [kg/kmol]" << endl;
    cout << "Volume:\t\t\t"				<< volume							<< " [m3]"		<< endl;
    cout << "Mass:\t\t\t"				<< mass								<< " [kg]"		<< endl;
    cout << "Moles:\t\t\t"				<< moles							<< " [kmol]"	<< endl;
	cout << "Temperature:\t\t"			<< T								<< " [K]"		<< endl;
    cout <<  endl;
}

void OpenSMOKE_Batch::videoSummary()
{
    VideoGeneralInfo();

    cout << "#\tName\tx\t\tomega" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (int i=1;i<=NC;i++)
        if (inletStream->omega[i]!=0.)
            cout << i << "\t\t" << mix->names[i] << "\t" << inletStream->x[i] << "\t" << inletStream->omega[i] << endl;
    cout << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;
}

void OpenSMOKE_Batch::MassAnalysis(OpenSMOKE_GasStream &outletStream)
{
	int j;

	// Elemental analysis
	BzzVectorDouble x_elemental(mix->NumberOfElements());
	BzzVectorDouble omega_elemental(mix->NumberOfElements());
	BzzVectorDouble x_elemental_initial(mix->NumberOfElements());
	BzzVectorDouble omega_elemental_initial(mix->NumberOfElements());
	
	mix->GetElementalMoleFractions(x, x_elemental);
	mix->GetElementalMassFractions(omega, omega_elemental);
	mix->GetElementalMoleFractions(inletStream->x, x_elemental_initial);
	mix->GetElementalMassFractions(inletStream->omega, omega_elemental_initial);

	cout << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "                       ELEMENTAL MASS ANALYSIS                       " << endl;
	cout << "---------------------------------------------------------------------" << endl;

	cout << "Name\txIN\t\txOUT" << endl;
    cout << "---------------------------------------------------------------------" << endl;
	for(j=1;j<=mix->NumberOfElements();j++)
		cout << mix->list_of_elements[j-1] << "\t" << x_elemental_initial[j]  << "\t" << x_elemental[j] << endl;

	cout << endl;
	cout << "Name\tomegaIN\t\tomegaOUT" << endl;
    cout << "---------------------------------------------------------------------" << endl;
	for(j=1;j<=mix->NumberOfElements();j++)
		cout << mix->list_of_elements[j-1] << "\t" << omega_elemental_initial[j]  << "\t" << omega_elemental[j] << endl;
}

void OpenSMOKE_Batch::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
{
	double H_mass_initial  = inletStream->massSpecificEnthalpy;			// [J/kg]
	double P_mass_initial  = inletStream->massSpecificPressureEnergy;	// [J/kg]
	double U_mass_initial  = inletStream->massSpecificInternalEnergy;	// [J/kg]
	
	double H_mass_final  = outletStream.massSpecificEnthalpy;			// [J/kg]
	double P_mass_final  = outletStream.massSpecificPressureEnergy;		// [J/kg]
	double U_mass_final  = outletStream.massSpecificInternalEnergy;		// [J/kg]

	double Work;	// Work done on the system (IN)
	if (iConstantPressure == true)	Work = -P*(volume-mass/inletStream->rho);
	else							Work = 0.;

	cout << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "                         ENERGY ANALYSIS                             " << endl;
	cout << "---------------------------------------------------------------------" << endl;
	cout << "Initial - Enthalpy:        " << H_mass_initial		<< " [J/kg]" << endl;
	cout << "Initial - Pressure Energy: " << P_mass_initial		<< " [J/kg]" << endl;
	cout << "Initial - Internal energy: " << U_mass_initial		<< " [J/kg]" << endl;
	cout << "Final   - Enthalpy:        " << H_mass_final		<< " [J/kg]" << endl;
	cout << "Final   - Pressure Energy: " << P_mass_final		<< " [J/kg]" << endl;
	cout << "Final   - Internal energy: " << U_mass_final		<< " [J/kg]" << endl;
	cout << endl;

	cout << "Initial energy:     "	<< U_mass_initial*mass								<< " [J]" << endl;
	cout << "Final energy:       "	<< U_mass_final*mass								<< " [J]" << endl;
	cout << "Work on the system: "	<< Work												<< " [J]" << endl;
	cout << "Balance (IN-OUT):   "  << Work + U_mass_initial*mass - U_mass_final*mass	<< " [J]" << endl;
}

void OpenSMOKE_Batch::OutletStream(OpenSMOKE_GasStream &outlet)
{
	outlet.AssignKineticScheme(*mix);
	outlet.AssignMassFlowRate(inletStream->massFlowRate, "kg/s");
	outlet.AssignTemperature(T, "K");
	outlet.AssignPressure(P,"Pa");
	outlet.AssignMassFractions(omega);
	outlet.lock();
}

void OpenSMOKE_Batch::videoFinalResult()
{
    int i;

    // Vectors
    BzzVectorDouble conversion(NC);

    // Print General Information about Batch
    VideoGeneralInfo();

    // Calculate Conversion of Inlet Species
    for (i=1;i<=NC;i++)
        if (inletStream->omega[i]!=0.)
            conversion[i] = (1. - omega[i]/inletStream->omega[i])*100.;

    // Final Composition
    cout << "#\tName\tx\t\tomega" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (i=1;i<=NC;i++)
        if (x[i]!=0.)
            cout << i << "\t" << mix->names[i] << "\t" << x[i] << "\t" << omega[i] << endl;
    cout << endl;

    // Conversion
    cout << "#\tName\tomegaInlet\tomegaOutlet\tConversion(%)" << endl;
    cout << "---------------------------------------------------------------------" << endl;
    for (i=1;i<=NC;i++)
        if (conversion[i]>0.)
            cout << i << "\t" << mix->names[i] << "\t"
            << inletStream->omega[i] << "\t" << omega[i] << "\t" << conversion[i] << endl;
    cout << endl;

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    cout << endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEM - CLASS DEFINITION												   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_Batch::ObjectBzzPrint(void)
{
    ::BzzPrint("\n OpenSMOKE_Batch Class\n");
}

void MyOdeSystem_Batch::GetSystemFunctions(BzzVectorDouble &x, double eta, BzzVectorDouble &f)
{
    if (iEnergy == false)
        ptBatch->ODESystem_Isothermal_Batch(x, eta, f);

    else if (iEnergy == true && iConstantPressure == false)
        ptBatch->ODESystem_NonIsothermal_ConstantVolume_Batch(x, eta, f);

    else if (iEnergy == true && iConstantPressure == true)
        ptBatch->ODESystem_NonIsothermal_ConstantPressure_Batch(x, eta, f);
}

void MyOdeSystem_Batch::assignBatch(OpenSMOKE_Batch *batch, const bool _iEnergy, const bool _iConstantPressure)
{
    ptBatch				= batch;
    iEnergy				= _iEnergy;
    iConstantPressure	= _iConstantPressure;
}

// Soot Two Equation Model
void OpenSMOKE_Batch::UpdateTwoEquationModel(BzzVectorDouble &y, BzzVectorDouble &dy)
{
    soot2EModel->update(T, P/101325., rho, 0., x, y[indexTwoEquations], y[indexTwoEquations+1]);
    soot2EModel->formation_rates();

    for (int i=1;i<=NC;i++)
        domega[i] += soot2EModel->SGas[i] / rho;

    dy[indexTwoEquations]	= soot2EModel->s / rho;
    dy[indexTwoEquations+1] = soot2EModel->S / rho;
}

double OpenSMOKE_Batch::SimplifiedViscosity()
{
	if (iUserDefinedViscosity == false)		return 1.81e-5 * sqrt(T/298.15);			// [Pa.s]
	else									return Viscosity * sqrt(T/inletStream->T);	// [Pa.s]
}	

void OpenSMOKE_Batch::UpdateHeatFlux(const double Tau)
{
	if (iUserDefinedHeatExchangeCoefficient == USERDEFINED)		U			= ud_U_profile.GiveMeValue(Tau, MINUSONE);
	if (iUserDefinedAmbientTemperature		== USERDEFINED)		Tambient	= ud_Tambient_profile.GiveMeValue(Tau, MINUSONE);
	if (iUserDefinedHeatFlux				== USERDEFINED)		Qe			= ud_Qe_profile.GiveMeValue(Tau, MINUSONE);
	
	if (iUserDefinedHeatExchangeCoefficient != NONE)			Qe = U*(T-Tambient);
}

void OpenSMOKE_Batch::UpdateExchangeArea(const double Tau)
{
	if (iUserDefinedExchangeArea == NONE)			A = pow(36.*Constants::pi*volume*volume, 1./3.);
	if (iUserDefinedExchangeArea == USERDEFINED)	A = ud_A_profile.GiveMeValue(Tau, MINUSONE);
}

void OpenSMOKE_Batch::SummaryOnFile(const string file_name)
{
	int i;

	ofstream fSummary;
	openOutputFileAndControl(fSummary, file_name);
    fSummary.setf(ios::scientific);

	double initial_volume	= mass / inletStream->rho;	// [m3]
	double initial_moles	= mass / inletStream->MW ;	// [kmol]


    fSummary << "                          \tInitial      \tEnd"								<< endl;
	fSummary << "----------------------------------------------------------------------------"	<< endl;
    fSummary << "Time[s]             \t"	<< 0.0					<< "\t" << 0.0				<< endl;
	fSummary << "Temperature[K]      \t"	<< inletStream->T		<< "\t" << T				<< endl;
	fSummary << "Pressure[Pa]        \t"	<< inletStream->P		<< "\t" << P				<< endl;
	fSummary << "Volume[m3]          \t"	<< initial_volume		<< "\t" << volume			<< endl;	
	fSummary << "Density[kg/m3]      \t"	<< inletStream->rho		<< "\t" << rho				<< endl;
	fSummary << "Mol.Weight[kg/kmol] \t"	<< inletStream->MW		<< "\t" << MWtot			<< endl;
	fSummary << "Mass[kg]            \t"	<< mass					<< "\t" << mass				<< endl;
	fSummary << "Moles[kmol]         \t"	<< initial_moles		<< "\t" << moles			<< endl;	
	fSummary << endl;
	
	
	fSummary << "Mole fractions[-]" << endl;
	for(i=1;i<=NC;i++)
	{
		fSummary.width(16);
		fSummary << mix->names[i] << "    \t"	<<	inletStream->x[i]		<< "\t" << x[i]			<< endl;
	}
	fSummary << endl;

	fSummary << endl;
	fSummary << "Mass fractions[-]" << endl;
	for(i=1;i<=NC;i++)
		fSummary << mix->names[i] << "    \t"	<<	inletStream->omega[i]	<< "\t" << x[i]			<< endl;
	fSummary << endl;
	
	fSummary << endl;
	fSummary << "Concentrations[kmol/m3]" << endl;
	for(i=1;i<=NC;i++)
		fSummary << mix->names[i] << "    \t"	<<	inletStream->c[i]		<< "\t" << x[i]			<< endl;
	fSummary << endl;

    fSummary.close();
}

void OpenSMOKE_Batch::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Batch"		<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_Batch::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Batch"		<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}


OpenSMOKE_Dictionary_Batch::OpenSMOKE_Dictionary_Batch()
{
    SetupBase();

	// Compulsory
    Add("#Energy",			'O', 'N', "The energy equation is solved");
    Add("#NoEnergy",		'O', 'N', "The energy equation is not solved");
    
	// Compulsory
    Add("#ContactTime",		'C', 'M', "Reactor contact time");

	// Compulsory
    Add("#ConstantVolume",		'O', 'N', "The volume is kept constant");
    Add("#ConstantPressure",	'O', 'N', "The pressure is kept constant");

	// Compulsory
    Add("#Volume",				'O', 'M', "Reactor volume");

	// Optional
	Add("#ReactionRates",		'O', 'V', "Reaction rates");
	Add("#FormationRates",		'O', 'V', "Formation rates");
	Add("#ROPA",				'O', 'V', "Rate of Production Analysis");
	Add("#Sensitivity",			'O', 'V', "Sensitivity Analysis");

     
	// Optional
	Add("#VerboseEnergy",			'O', 'N', "Report on energy");
	Add("#Qe",						'O', 'M', "Heat flux");
	Add("#A",						'O', 'M', "Exchange area");
	Add("#U",						'O', 'M', "Heat thermal exchange coefficient");
	Add("#Tambient",				'O', 'M', "Ambient temperature");
	Add("#Viscosity",				'O', 'M', "Inlet viscosity");
	Add("#UserDefined_T",			'O', 'S', "The temperature is assigned from external file");
	Add("#UserDefined_Qe",			'O', 'S', "The heat flux is assigned from external file");
    Add("#UserDefined_A",			'O', 'S', "The exchange area is assigned from external file");
	Add("#UserDefined_U",			'O', 'S', "The heat exchange coefficient is assigned from external file");
    Add("#UserDefined_Tambient",	'O', 'S', "The ambient temperature is assigned from external file");

	Add("#2E_Model",				'O', 'N', "The Soot 2E Model is solved");

    Conflict("#Energy",					"#NoEnergy");
    Conflict("#Energy",					"#UserDefined_T");
	Conflict("#NoEnergy",				"#Qe");
	Conflict("#NoEnergy",				"#A");
	Conflict("#UserDefined_Qe",			"#Qe");
	Conflict("#UserDefined_A",			"#A");
	Conflict("#UserDefined_U",			"#U");
	Conflict("#UserDefined_Tambient",	"#Tambient");

	Conflict("#NoEnergy",				"#U");
	Conflict("#NoEnergy",				"#Tambient");
	Conflict("#Qe",						"#U");
	Conflict("#Qe",						"#Tambient");
	Conflict("#Qe",						"#UserDefined_U");
	Conflict("#Qe",						"#UserDefined_Tambient");
	Conflict("#UserDefined_Qe",			"#U");
	Conflict("#UserDefined_Qe",			"#Tambient");
	Conflict("#UserDefined_Qe",			"#UserDefined_U");
	Conflict("#UserDefined_Qe",			"#UserDefined_Tambient");

	Compulsory("#Energy",			"#NoEnergy");
	Compulsory("#ConstantVolume",   "#ConstantPressure");

    Lock();
}
