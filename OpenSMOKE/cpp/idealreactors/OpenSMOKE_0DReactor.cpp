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
#include "idealreactors/OpenSMOKE_0DReactor.h"
#include "basic/OpenSMOKE_Conversions.h"
#include "addons/OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D.h"
#include "addons/OpenSMOKE_PolimiSoot.h"


const double    OpenSMOKE_0DReactor::ONE						=  1.0;
const double    OpenSMOKE_0DReactor::ZERO						=  0.0;
const double    OpenSMOKE_0DReactor::MAXNUMBER					=  1.e16;
const double    OpenSMOKE_0DReactor::MINPRESSURE				=  10.0;			// [Pa]
const int       OpenSMOKE_0DReactor::MAX_TIME_STEPS				=  12000;
//const int       OpenSMOKE_0DReactor::MAX_TIME_STEPS				=  5000;
const int       OpenSMOKE_0DReactor::MINUSONE					= -1;
const double    OpenSMOKE_0DReactor::MAXERRORONMASSFRACTIONS	=  1.e-5;

void OpenSMOKE_0DReactor::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Parent: OpenSMOKE_0DReactor"	<< endl;
    cout << "Class:  " << class_name		<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message			<< endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_0DReactor::WarningMessage(const string message)
{
    cout << endl;
    cout << "Parent:  OpenSMOKE_0DReactor"	<< endl;
    cout << "Class:   " << class_name		<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message			<< endl;
	cout << endl;
}

void OpenSMOKE_0DReactor::Setup()
{
    // Control Variables
    assignedKineticScheme		= false;			//
    assignedGlobalKineticScheme	= false;			//
	assignedSoot2EModel			= false;			//
    assignedInletFlows			= false;			//
    assignedEnd					= false;			//
	assignedEnergy				= false;			// 
    

    // Properties Index
    iGlobalKinetics             = false;			// Global Kinetics 0=OFF 1=ON
    iEnergy						= false;			// false = Isothermal true = Energy Equation
    iTwoEquationModel	        = false;			// 2E Model for soot predictions
    iRelativeTolerance	        = false;			// User defined relative tolerance
    iAbsoluteTolerance	        = false;			// User defined absolute tolerance
    iHistory                    = false;			// Save solution history
    iVerbose                    = true;				// Write on video and file
    iVerboseEnergy              = false;			// Report on energy
    iUserDefinedTemperature     = false;			// User defined temperature
	iUserDefinedViscosity		= false;			// User defined viscosity

	iVerboseReactionRates		= false;			// Reaction rates
	iVerboseFormationRates		= false;			// Formation rates
	iAssignedROPA				= false;			// Rate of production analysis
	iVerboseAssignedROPA		= false;			// Rate of production analysis (verbose)
	iVerboseSensitivity			= false;			// Sensitivity analysis
	iVerboseElementFluxAnalysis = false;
	iVerboseLocalElementFluxAnalysis	= false;	// Rate of production analysis (local)
	iVerboseSolid						= false;	// Solid output
	iVerboseExperiment					= false;	// Experiment
	iVerboseSelectivity					= false;	// Reaction rates

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

	// Internal variables
    timeOld = 0.;    // Old time

    // Output Folder
	SetOutputFolder("Output");

	// Name Object
	name_object = "[not assigned]";

    // Default Values
    Qe						= 0.0;
	Viscosity				= 0.0;
}

void OpenSMOKE_0DReactor::AssignKineticScheme(OpenSMOKE_ReactingGas &_mix)
{
    mix     = &_mix;

    NC = mix->NumberOfSpecies();
    NR = mix->NumberOfReactions();

    // These vectors are used for the evaluation of gas mixture properties
    ChangeDimensions(NC, &omega);
    ChangeDimensions(NC, &domega);
    ChangeDimensions(NC, &x);
    ChangeDimensions(NC, &c);
    ChangeDimensions(NC, &R);
    ChangeDimensions(NC, &Rtilde);

    ChangeDimensions(NC, &CpMap);
    ChangeDimensions(NR, &k1Map);
    ChangeDimensions(NR, &k2Map);
    ChangeDimensions(NR, &uKeqMap);
    ChangeDimensions(NR, &logFcentMap);
    ChangeDimensions(NR, &reactionDHMap);
    ChangeDimensions(NR, &reactionDSMap);

	// Output species
	ChangeDimensions(mix->NumberOfSpecies(), &iOutputSpecies);
	namesOutputSpecies.resize(mix->NumberOfSpecies()+1); 
	for(int i=1;i<=mix->NumberOfSpecies();i++)
	{
		namesOutputSpecies[i]	= mix->names[i];
		iOutputSpecies[i]		= i;
	}

    // Control Variables
    assignedKineticScheme = true;
}

void OpenSMOKE_0DReactor::AssignInletFlows(OpenSMOKE_GasStream &_inletStream)
{
    if (assignedKineticScheme == false)
        ErrorMessage("You must define the kinetic scheme before defining the inlet/initial flows!");

    inletStream = &_inletStream;

    assignedInletFlows = true;
}

void OpenSMOKE_0DReactor::AssignSoot2EModel(OpenSMOKE_2EModel &_soot2EModel)
{
    soot2EModel = &_soot2EModel;

    assignedSoot2EModel = true;
}

void OpenSMOKE_0DReactor::AssignGlobalKineticScheme(OpenSMOKE_GlobalKinetics &_global)
{
    global = &_global;

    assignedGlobalKineticScheme = true;
}

void OpenSMOKE_0DReactor::AssignEnergy(const bool _iEnergy)
{
	iEnergy = _iEnergy;

    assignedEnergy = true;
}

void OpenSMOKE_0DReactor::SetName(const string name)
{
    name_object = name;
}

void OpenSMOKE_0DReactor::SetGlobalKinetics()
{
    iGlobalKinetics = true;
}

void OpenSMOKE_0DReactor::UnsetGlobalKinetics()
{
    iGlobalKinetics = false;
}

void OpenSMOKE_0DReactor::SetTwoEquationModel()
{
    iTwoEquationModel = true;
}

void OpenSMOKE_0DReactor::UnsetTwoEquationModel()
{
    iTwoEquationModel = false;
}

void OpenSMOKE_0DReactor::SetVerbose()
{
    iVerbose = true;
}

void OpenSMOKE_0DReactor::UnsetVerbose()
{
    iVerbose = false;
}

void OpenSMOKE_0DReactor::SetRelativeTolerance(const double relative)
{
    iRelativeTolerance	= true;
    relativeTolerance	= relative;
}

void OpenSMOKE_0DReactor::UnsetRelativeTolerance()
{
    iRelativeTolerance = false;
}

void OpenSMOKE_0DReactor::SetAbsoluteTolerance(const double absolute)
{
    iAbsoluteTolerance = true;
    absoluteTolerance = absolute;
}

void OpenSMOKE_0DReactor::UnsetAbsoluteTolerance()
{
    iAbsoluteTolerance = false;
}

void OpenSMOKE_0DReactor::SetUserDefinedTemperature(const string fileName)
{
    iUserDefinedTemperature = true;
    ud_temperature_profile.AssignFromFile(fileName, "TEMPERATURE");
	ud_temperature_profile.SetName(name_object + " - Temperature Profile");
}

void OpenSMOKE_0DReactor::SetUserDefinedTemperature(OpenSMOKE_UD_Profile &udp)
{
    iUserDefinedTemperature = true;
	ud_temperature_profile = udp;
	ud_temperature_profile.SetName(name_object + " - Temperature Profile");
}

void OpenSMOKE_0DReactor::UnsetUserDefinedTemperature()
{
    iUserDefinedTemperature = false;
}

void OpenSMOKE_0DReactor::SetUserDefinedHeatFlux(const string fileName)
{
    iUserDefinedHeatFlux = USERDEFINED;
    ud_Qe_profile.AssignFromFile(fileName, "HEAT_FLUX");
	ud_Qe_profile.SetName(name_object + " - Heat Flux Profile");
}

void OpenSMOKE_0DReactor::UnsetUserDefinedHeatFlux()
{
    iUserDefinedHeatFlux = NONE;
}

void OpenSMOKE_0DReactor::SetUserDefinedHeatExchangeCoefficient(const string fileName)
{
    iUserDefinedHeatExchangeCoefficient = USERDEFINED;
    ud_U_profile.AssignFromFile(fileName, "HEAT_TRANSFER_COEFFICIENT");
	ud_U_profile.SetName(name_object + " - Heat Exchange Coefficient Profile");
}

void OpenSMOKE_0DReactor::UnsetUserDefinedHeatExchangeCoefficient()
{
    iUserDefinedHeatExchangeCoefficient = NONE;
}

void OpenSMOKE_0DReactor::SetUserDefinedAmbientTemperature(const string fileName)
{
    iUserDefinedAmbientTemperature = USERDEFINED;
    ud_Tambient_profile.AssignFromFile(fileName, "TEMPERATURE");
	ud_Tambient_profile.SetName(name_object + " - Ambient Temperature Profile");
}

void OpenSMOKE_0DReactor::UnsetUserDefinedAmbientTemperature()
{
    iUserDefinedAmbientTemperature = NONE;
}

void OpenSMOKE_0DReactor::SetVerboseEnergy()
{
    iVerboseEnergy = true;
}

void OpenSMOKE_0DReactor::UnsetVerboseEnergy()
{
    iVerboseEnergy = false;
}

void OpenSMOKE_0DReactor::SetKeySpecies(const vector<string> key_species)
{
    iKeySpecies			= true;

	for(int k=0;k<int(key_species.size());k++)
	{
		key_species_names.push_back(key_species[k]);
		key_species_index.push_back(mix->recognize_species(key_species[k]));
	}
}

void OpenSMOKE_0DReactor::UnsetKeySpecies()
{
    iKeySpecies			= false;
	
	key_species_names.clear();
	key_species_index.clear();
}

void OpenSMOKE_0DReactor::SetHistory()
{
    iHistory = true;

    ChangeDimensions(MAX_TIME_STEPS, &P_History);
    ChangeDimensions(MAX_TIME_STEPS, &Tau_History);
    ChangeDimensions(MAX_TIME_STEPS, &Csi_History);
    ChangeDimensions(MAX_TIME_STEPS, &T_History);
    ChangeDimensions(MAX_TIME_STEPS, mix->NumberOfSpecies(), &mass_History);
    ChangeDimensions(MAX_TIME_STEPS, mix->NumberOfSpecies(), &mole_History);
}

void OpenSMOKE_0DReactor::UnsetHistory()
{
    iHistory = false;
}

void OpenSMOKE_0DReactor::SetPostProcessing()
{
	SetHistory();
}

void OpenSMOKE_0DReactor::SetConstantHeatFlux(const double value, const string units)
{
	iUserDefinedHeatFlux = CONSTANT;
	Qe = OpenSMOKE_Conversions::conversion_heat_flux(value, units);
}

void OpenSMOKE_0DReactor::SetConstantHeatExchangeCoefficient(const double value, const string units)
{
	iUserDefinedHeatExchangeCoefficient = CONSTANT;
	U = OpenSMOKE_Conversions::conversion_heat_exchange_coefficient(value, units);
}

void OpenSMOKE_0DReactor::SetConstantAmbientTemperature(const double value, const string units)
{
	iUserDefinedAmbientTemperature = CONSTANT;
	Tambient = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_0DReactor::SetViscosity(const double value, const string units)
{
	iUserDefinedViscosity = true;
	Viscosity = OpenSMOKE_Conversions::conversion_dynamic_viscosity(value, units);
}

void OpenSMOKE_0DReactor::UnsetViscosity()
{
	iUserDefinedViscosity = false;
	Viscosity = 0.;
}

void OpenSMOKE_0DReactor::SetOptions(const string option)
{
    if (option == "NOTHING")
    {
    }
    else
        ErrorMessage("This option is not allowed: " + option);
}

void OpenSMOKE_0DReactor::SetOutputFolder(const string _outputFolderName)
{
    outputFolderName			= _outputFolderName;
    outputName					=  outputFolderName + "/" + out_name;
    outputSootName				=  outputFolderName + "/" + "Soot.out";
//    outputPAHName				=  outputFolderName + "/" + "PAH.out";
    output2EName				=  outputFolderName + "/" + "2E.out";
    outputEnergyName			=  outputFolderName + "/" + "Energy.out";
	outputReactionRatesName		=  outputFolderName + "/" + "ReactionRates.out";
    outputFormationRatesName	=  outputFolderName + "/" + "FormationRates.out";
    outputSolidName				=  outputFolderName + "/" + "Solid.out";
    outputConversionsName		=  outputFolderName + "/" + "Conversions.out";
    outputOSMName				=  outputFolderName + "/" + "IdealReactor.osm";	
	outputExperimentName		=  outputFolderName + "/Exp/Exp.inp";
    outputSensitivityName		=  outputFolderName + "\\" + "Sensitivity.bin";	
    outputSelectivityName		=  outputFolderName + "\\" + "Selectivity.out";	
}

void OpenSMOKE_0DReactor::SetVideoOptions(const int _nVideoSteps)
{
    nVideoSteps		=   countVideoSteps		= _nVideoSteps;
                        countIterations     = 0;
}

void OpenSMOKE_0DReactor::SetFileOptions(const int _nFileSteps)
{
    nFileSteps		=   countFileSteps		= _nFileSteps;
                        countIterations     = 0;
}

void OpenSMOKE_0DReactor::SetReactionRatesOnFile(const vector<string> _names)
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
		for(int k=0;k<int(_names.size());k++)
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

void OpenSMOKE_0DReactor::SetSelectivitiesOnFile()
{
	if (iKeySpecies == false)
		ErrorMessage("The selectivity analysis requires the #Key option!");
	iVerboseSelectivity = true;
	ChangeDimensions(NC, &selectivity_mass);
}

void OpenSMOKE_0DReactor::SetFormationRatesOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_formation_rates, names_formation_rates);
	iVerboseFormationRates = true;
}

void OpenSMOKE_0DReactor::SetVerboseROPAOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_ROPA, names_ROPA);
	iAssignedROPA = true;
	iVerboseAssignedROPA = true;

	ropa.Initialize(mix, index_ROPA);
	SetHistory();
}

void OpenSMOKE_0DReactor::SetROPAOnFile()
{
	iAssignedROPA = true;
	SetHistory();
}

void OpenSMOKE_0DReactor::SetElementFluxAnalysisOnFile(const vector<string> _names)
{
	iVerboseElementFluxAnalysis = true;

	ElementFluxAnalysis.Initialize(mix, _names);
	SetHistory();
}

void OpenSMOKE_0DReactor::SetLocalElementFluxAnalysisOnFile(const vector<string> _names)
{
	string units = _names[_names.size()-1];
	for(unsigned int j=1;j<=_names.size()-1;j++)
		index_local_ElementFluxAnalysis.Append(OpenSMOKE_Conversions::conversion_time(atof(_names[j-1].c_str()), units));
	iVerboseLocalElementFluxAnalysis = true;
}


void OpenSMOKE_0DReactor::SetExperimentOnFile(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix, _names, index_Experiment, names_Experiment);
	iVerboseExperiment = true;

	string nameFolderExperiments;
	#if LINUX_SO==1
		nameFolderExperiments = outputFolderName + "/Exp";
	#else
		nameFolderExperiments = outputFolderName + "\\Exp";
	#endif
	
	string MSDOScommand = "mkdir " + nameFolderExperiments;
	system(MSDOScommand.c_str());
}

void OpenSMOKE_0DReactor::SetSensitivityOptions(const vector<string> _names)
{
//	ErrorMessage("Sensitivity Analysis: Coming soon!!!!!!");

	sensitivityOptions = _names;
}

void OpenSMOKE_0DReactor::SetSensitivityOnFile(const vector<string> _names)
{
	iVerboseSensitivity = true;
	SetFileOptions(1);

	GiveMeIndicesAndNames(*mix, _names, index_sensitivity, names_sensitivity);

	// TODO
	sensitivity_fast = new OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D();

	kind_of_sensitivity_analysis kind;
	if (iEnergy == false)
	{
		cout << "Sensitivity Analysis: Only Species" << endl;
		kind = OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ONLY_SPECIES;
	}
	else if (iEnergy == true)
	{
		cout << "Sensitivity Analysis: Species + Temperature" << endl;
		kind = OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_SPECIES_PLUS_TEMPERATURE;
	}

	if (sensitivityOptions.size() == 0)
	{
		sensitivity_fast->Initialize(	kind, 
										OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_IMPLICIT, 
										mix, index_sensitivity, outputFolderName);
	}
	else
	{
		if (sensitivityOptions[0] == "Euler_Implicit")
			sensitivity_fast->Initialize(	kind, 
											OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_IMPLICIT, 
											mix, index_sensitivity, outputFolderName);
		else if (sensitivityOptions[0] == "Euler_Explicit")
			sensitivity_fast->Initialize(	kind, 
											OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_EXPLICIT, 
											mix, index_sensitivity, outputFolderName);
		else if (sensitivityOptions[0] == "Accurate")
			sensitivity_fast->Initialize(	kind, 
											OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE, 
											mix, index_sensitivity, outputFolderName);
		else
			ErrorMessage("Wrong option #SensitivityOptions: " + sensitivityOptions[0]);
	}
} 

void OpenSMOKE_0DReactor::SetOutputSpecies(const vector<string> string_vector)
{
	if (string_vector.size()==1 && string_vector[0] == "ALL")
	{
	}
	else
	{
		int iOutput = string_vector.size();
		ChangeDimensions(iOutput, &iOutputSpecies);
		namesOutputSpecies.resize(iOutput+1);
		for(int i=1;i<=iOutput;i++)
		{
			namesOutputSpecies[i] = string_vector[i-1];
			iOutputSpecies[i] = mix->recognize_species(string_vector[i-1]);
		}
	}
}

void OpenSMOKE_0DReactor::PrepareFiles()
{
	if (iVerbose == true)
    {
        // Create output folder
        string message = "mkdir " + outputFolderName;
        system(message.c_str());

        openOutputFileAndControl(fOutput, outputName.c_str());
        fOutput.setf(ios::scientific);

        if (iTwoEquationModel==true)
        {
            openOutputFileAndControl(fSoot2E, output2EName.c_str());
            fSoot2E.setf(ios::scientific);
        }
//
//        if (mix->pah_manager.iPAH == 1)
 //       {
 //           openOutputFileAndControl(fPAH, outputPAHName.c_str());
  //          fPAH.setf(ios::scientific);
   //     }

		mix->SetPolimiSoot();
        if (mix->polimiSoot->IsSoot() == true)
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

		if (iAssignedROPA == true)
        {
        }

		if (iVerboseElementFluxAnalysis == true)
		{
		}

		if (iVerboseSensitivity == true)
        {
        }

		if (iKeySpecies == true && key_species_names.size()>1)
        {
            openOutputFileAndControl(fConversions, outputConversionsName);
            fConversions.setf(ios::scientific);
        }

		if (iVerboseSolid == true)
        {
            openOutputFileAndControl(fSolid, outputSolidName);
            fSolid.setf(ios::scientific);
        }

		if (iVerbose == true)
		{
			// Create output folder
			string instruction;
			instruction = "copy " + mix->path_complete + "\\idealgas.bin " + outputFolderName + "\\idealgas.bin /Y";
			system(instruction.c_str());
			instruction = "copy " + mix->path_complete + "\\reactions.bin " + outputFolderName + "\\reactions.bin /Y";
			system(instruction.c_str());
		}

		if (iVerboseSelectivity)
        {
            openOutputFileAndControl(fSelectivity, outputSelectivityName);
            fSelectivity.setf(ios::scientific);
        }

        LabelODEFile();
    }
}

void OpenSMOKE_0DReactor::CloseFiles()
{
    if (iVerbose == true)
    {
        fOutput.close();

//        if (mix->pah_manager.iPAH == 1)		fPAH.close();
        if (mix->polimiSoot->IsSoot() == true)	fSoot.close();
		if (iVerboseEnergy == true)			fEnergy.close();
		if (iVerboseReactionRates == true)	fReactionRates.close();
		if (iVerboseFormationRates == true)	fFormationRates.close();
		if (iVerboseSolid == true)			fSolid.close();
		if (iVerboseSensitivity == true)	fSensitivity.close();
		if (iKeySpecies == true && 
			key_species_names.size()>1)		fConversions.close();
		if (iVerboseSelectivity == true)	fSelectivity.close();
    }
}

double OpenSMOKE_0DReactor::SimplifiedViscosity()
{
	if (iUserDefinedViscosity == false)		return 1.81e-5 * sqrt(T/298.15);			// [Pa.s]
	else									return Viscosity * sqrt(T/inletStream->T);	// [Pa.s]
}


void OpenSMOKE_0DReactor::MassAnalysis(OpenSMOKE_GasStream &outletStream)
{
	int j;

	// Elemental analysis
	BzzVector x_elemental(mix->NumberOfElements());
	BzzVector omega_elemental(mix->NumberOfElements());
	BzzVector x_elemental_initial(mix->NumberOfElements());
	BzzVector omega_elemental_initial(mix->NumberOfElements());
	
	mix->GetElementalMoleFractionsFromSpeciesMoleFractions(x_elemental, x);;
	mix->GetElementalMassFractionsFromSpeciesMassFractions(omega_elemental, omega);
	mix->GetElementalMoleFractionsFromSpeciesMoleFractions(x_elemental_initial, inletStream->x);
	mix->GetElementalMassFractionsFromSpeciesMassFractions(omega_elemental_initial, inletStream->omega);
	
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


void OpenSMOKE_0DReactor::EnergyAnalysis(OpenSMOKE_GasStream &outletStream)
{
}

void OpenSMOKE_0DReactor::OutletStream(OpenSMOKE_GasStream &outlet)
{
	outlet.AssignKineticScheme(*mix);
	outlet.AssignMassFlowRate(inletStream->massFlowRate, "kg/s");
	outlet.AssignTemperature(T, "K");
	outlet.AssignPressure(P,"Pa");

	double sum = omega.GetSumElements();
	
	if (sum > 1.+MAXERRORONMASSFRACTIONS || sum < 1.-MAXERRORONMASSFRACTIONS)
	{
		cout << "Sum of mass fractions: " << sum << endl;
		ErrorMessage("The sum of mass fractions must be equal to 1 (outlet stream)");
	}
	omega /= sum;
	outlet.AssignMassFractions(omega);

	outlet.lock();
}

