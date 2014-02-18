/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
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

#include "OpenSMOKE_KPP_DataManager.h"
#include "OpenSMOKE_KPP_Dictionary.h"
#include "OpenSMOKE_KPP_ODE_Manager.h"
#include "OpenSMOKE_KPP_NewtonMethod_Manager.h"
#include <iostream>
#include <iomanip>

void OpenSMOKE_KPP_DataManager::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_DataManager"	<< endl;
    cout << "Error: " << message_					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_DataManager::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_DataManager"	<< endl;
    cout << "Warning: " << message_					<< endl;
	cout << endl;
}

OpenSMOKE_KPP_DataManager::OpenSMOKE_KPP_DataManager(OpenSMOKE_KPP_Dictionary& dictionary) :
dictionary_(dictionary)
{
	// Default values
	nameFolderOutput_					= "Output";
	nameFolderInput_					= "Input";
	nameInputBackupFile_				= "Backup";
	nameFolderKinetics_					= "Kinetics";
	nameTopologyFile_					= "CFDNetwork.bzz";
	nameFirstGuessFile_					= "FirstGuess.bzz";
	nameInOutFile_						= "InputOutput.out";
	nameResidualFile_					= "Residuals.out";
	nameResidualFile_					= "ResidualsSpecies.out";
	nameResidualReactorStatisticsFile_	= "StatisticsReactors.out";
	nameResidualSpeciesStatisticsFile_	= "StatisticsSpecies.out";

	networkStatus_		  	= KPP_NETWORK_STATUS_START;

	iReactions_			= false;
	iKineticMaps_			= false;
	iSymbolicKinetics_		= false;
	iSaveKineticConstants_		= true;
	iOpenMP_			= false;
	iBackup_			= false;
	nThreads_			= 1;

	iTraditionalLink_ 		= false;
	

	correction_ = KPP_CORRECTION_NONE;
	
	massFlowRate_SparseLinearSolver_		= KPP_SPARSESOLVER_PARDISO;

	fillInReducingOrdering_ = PARDISO_FILLIN_METIS;

	// Predictor Corrector
	predictorCorrector_SparseLinearSolver_				= KPP_SPARSESOLVER_PARDISO;
	predictorCorrector_TimeStepIncrementFactor_			= 1.50;
	predictorCorrector_TimeStepReductionFactor_			= 0.80;
	predictorCorrector_MaxIterations_					= 10000;
	predictorCorrector_CourantCorrectionCoefficient_	= 0.20;
	predictorCorrector_MultiTimeSplitting_				= false;
	predictorCorrector_InitialTimeStep_					= 1.e-5;
	predictorCorrector_MaxTimeStep_						= 1e5*predictorCorrector_InitialTimeStep_;
	predictorCorrector_DeferredConvection_				= false;
	predictorCorrector_AlgebraicConstraint_				= false;

	// Sequence CSTR
	sequenceCSTR_MaxIterations_				= 10000;
	sequenceCSTR_VerboseStatistics_			= false;
	sequenceCSTR_UpdatingNetwork_			= false;

	// Single Reactor
	singleReactor_OdeMaxJacobian_			= 0;
	singleReactor_OdeStopResiduals_			= 0;
	singleReactor_IntegrationTime_			= 100.;
	singleReactor_OdeRelativeTolerance_		= 100.*MachEpsFloat();
	singleReactor_OdeAbsoluteTolerance_		= 1.e-10;

	temperatureMin_ = 0.;
	temperatureMax_ = 0.;


	// Global ODE
	globalODE_SparseLinearSolver_			= KPP_SPARSESOLVER_PARDISO;
	globalODE_SafetyReductionCoefficient	= OpenSMOKE_KPP_ODE_Manager::default_safetyReductionCoefficient_;
	globalODE_RelativeTolerance_			= OpenSMOKE_KPP_ODE_Manager::default_relativeTolerance_;
	globalODE_AbsoluteTolerance_			= OpenSMOKE_KPP_ODE_Manager::default_absoluteTolerance_;
	globalODE_MaxIterations_				= OpenSMOKE_KPP_ODE_Manager::default_maximumIterations_;
	globalODE_TimeStepIncrementFactor_		= OpenSMOKE_KPP_ODE_Manager::default_timeStepIncrementFactor_;
	globalODE_TimeStepReductionFactor_		= 1./OpenSMOKE_KPP_ODE_Manager::default_timeStepIncrementFactor_;
	globalODE_InitialTimeStep_				= OpenSMOKE_KPP_ODE_Manager::default_initialTimeStep_; 
	globalODE_MaxTimeStep_					= OpenSMOKE_KPP_ODE_Manager::default_maxTimeStep_; 
	globalODE_UpdatingFrequencyTimeStep_	= OpenSMOKE_KPP_ODE_Manager::default_updatingFrequencyTimeStep_;

	// Global NLS
	globalNLS_SparseLinearSolver_			= KPP_SPARSESOLVER_PARDISO;
	globalNLS_SafetyReductionCoefficient	= OpenSMOKE_KPP_NewtonMethod_Manager::default_safetyReductionCoefficient_;
	globalNLS_Method_						= OpenSMOKE_KPP_NewtonMethod_Manager::default_method_;
	globalNLS_RelativeTolerance_			= OpenSMOKE_KPP_NewtonMethod_Manager::default_relativeTolerance_;
	globalNLS_AbsoluteTolerance_			= OpenSMOKE_KPP_NewtonMethod_Manager::default_absoluteTolerance_;
	globalNLS_MaxIterations_				= OpenSMOKE_KPP_NewtonMethod_Manager::default_maxit_;
	globalNLS_MaxArmijioIterations_			= OpenSMOKE_KPP_NewtonMethod_Manager::default_maxarm_;
	
	// Check Dictionary
	CheckDictionary();
	CheckUserInput();
}

OpenSMOKE_KPP_DataManager::~OpenSMOKE_KPP_DataManager(void)
{
}

void OpenSMOKE_KPP_DataManager::CheckDictionary()
{
	int     int_value;
	double  double_value;
    string  string_value;
	vector<string> string_vector;

	if (dictionary_.Return("#Kinetics", string_value))
		AssignKinetics(string_value);

	if (dictionary_.Return("#Input", string_value))
		AssignInputFolder(string_value);

	if (dictionary_.Return("#Output", string_value))
		AssignOutputFolder(string_value);

	if (dictionary_.Return("#Backup", string_value))
		AssignBackupFile(string_value);

	if (dictionary_.Return("#Reactions", string_value))
		AssignReactions(string_value);

	if (dictionary_.Return("#Correction", string_value))
		AssignCorrection(string_value);

	if (dictionary_.Return("#OpenMP", string_value))
		AssignOpenMP(string_value);

	if (dictionary_.Return("#Threads", int_value))
		AssignThreads(int_value);

	if (dictionary_.Return("#SaveKineticConstants", string_value))
		AssignSaveKineticConstants(string_value);

	if (dictionary_.Return("#Maps", string_vector))
		AssignMaps(string_vector);

	if (dictionary_.Return("#MassFlowRate_SparseLinearSolver", string_value))
		AssignSparseLinearSolverMassFlowRate(string_value);

	if (dictionary_.Return("#FillInReducingOrdering", string_value))
		SetFillInReducingOrdering(string_value);

	if (dictionary_.Return("#KineticMaps", string_value))
		AssignKineticMaps(string_value);

	if (dictionary_.Return("#SymbolicKinetics", string_value))
		AssignSymbolicKinetics(string_value);

	if (dictionary_.Return("#MinimumTemperature", double_value, string_value))
		SetMinimumTemperature(double_value, string_value);
	if (dictionary_.Return("#MaximumTemperature", double_value, string_value))
		SetMaximumTemperature(double_value, string_value);
	
	// Single Reactor options

	if (dictionary_.Return("#SingleReactor_IntegrationTime", double_value, string_value))
		SetSingleReactor_IntegrationTime(double_value, string_value);

	if (dictionary_.Return("#SingleReactor_OdeMaxJacobian", int_value))
		SetSingleReactor_OdeMaxJacobian(int_value);

	if (dictionary_.Return("#SingleReactor_OdeStopResiduals", double_value))
		SetSingleReactor_OdeStopResiduals(double_value);

	if (dictionary_.Return("#SingleReactor_AbsoluteTolerance", double_value))
		SetSingleReactor_AbsoluteTolerance(double_value);

	if (dictionary_.Return("#SingleReactor_RelativeTolerance", double_value))
		SetSingleReactor_RelativeTolerance(double_value);


	// Sequence Options
	if (dictionary_.Return("#SequenceCSTR_UpdatingNetwork", string_value))
		AssignSequenceCSTR_UpdatingNetwork(string_value);

	if (dictionary_.Return("#SequenceCSTR_MaxIterations", int_value))
		SetSequenceCSTR_MaxIterations(int_value);

	if (dictionary_.Return("#SequenceCSTR_VerboseStatistics", string_value))
		SetSequenceCSTR_VerboseStatistics(string_value);


	// PredictorCorrector options

	if (dictionary_.Return("#PredictorCorrector_SparseLinearSolver", string_value))
		AssignPredictorCorrector_SparseLinearSolver(string_value);

	if (dictionary_.Return("#PredictorCorrector_DeferredConvection", string_value))
		AssignPredictorCorrector_DeferredConvection(string_value);

	if (dictionary_.Return("#PredictorCorrector_AlgebraicConstraint", string_value))
		AssignPredictorCorrector_AlgebraicConstraint(string_value);

	if (dictionary_.Return("#PredictorCorrector_InitialTimeStep", double_value, string_value))
		AssignPredictorCorrector_InitialTimeStep(double_value, string_value);

	if (dictionary_.Return("#PredictorCorrector_MaxTimeStep", double_value, string_value))
		SetPredictorCorrector_MaxTimeStep(double_value, string_value);

	if (dictionary_.Return("#PredictorCorrector_TimeStepIncrementFactor", double_value))
		SetPredictorCorrector_TimeStepIncrementFactor(double_value);

	if (dictionary_.Return("#PredictorCorrector_TimeStepReductionFactor", double_value))
		SetPredictorCorrector_TimeStepReductionFactor(double_value);

	if (dictionary_.Return("#PredictorCorrector_MaxIterations", int_value))
		SetPredictorCorrector_MaxIterations(int_value);

	if (dictionary_.Return("#PredictorCorrector_MultiTimeSplitting", string_value))
		SetPredictorCorrector_MultiTimeSplitting(string_value);

	if (dictionary_.Return("#PredictorCorrector_CourantCorrectionCoefficient", double_value))
		SetPredictorCorrector_CourantCorrectionCoefficient(double_value);


	// GlobalODE options

	if (dictionary_.Return("#GlobalODE_SparseLinearSolver", string_value))
		AssignGlobalODE_SparseLinearSolver(string_value);

	if (dictionary_.Return("#GlobalODE_AbsoluteTolerance", double_value))
		SetGlobalODE_AbsoluteTolerance(double_value);

	if (dictionary_.Return("#GlobalODE_RelativeTolerance", double_value))
		SetGlobalODE_RelativeTolerance(double_value);

	if (dictionary_.Return("#GlobalODE_MaxIterations", int_value))
		SetGlobalODE_MaxIterations(int_value);

	if (dictionary_.Return("#GlobalODE_TimeStepIncrementFactor", double_value))
		SetGlobalODE_TimeStepIncrementFactor(double_value);

	if (dictionary_.Return("#GlobalODE_TimeStepReductionFactor", double_value))
		SetGlobalODE_TimeStepReductionFactor(double_value);

	if (dictionary_.Return("#GlobalODE_InitialTimeStep", double_value, string_value))
		SetGlobalODE_InitialTimeStep(double_value, string_value);

	if (dictionary_.Return("#GlobalODE_MaxTimeStep", double_value, string_value))
		SetGlobalODE_MaxTimeStep(double_value, string_value);

	if (dictionary_.Return("#GlobalODE_UpdatingFrequencyTimeStep", int_value))
		SetGlobalODE_UpdatingFrequencyTimeStep(int_value);


	// GlobalNLS options

	if (dictionary_.Return("#GlobalNLS_SparseLinearSolver", string_value))
		AssignGlobalNLS_SparseLinearSolver(string_value);

	if (dictionary_.Return("#GlobalNLS_AbsoluteTolerance", double_value))
		SetGlobalNLS_AbsoluteTolerance(double_value);

	if (dictionary_.Return("#GlobalNLS_RelativeTolerance", double_value))
		SetGlobalNLS_RelativeTolerance(double_value);

	if (dictionary_.Return("#GlobalNLS_MaxIterations", int_value))
		SetGlobalNLS_MaxIterations(int_value);

	if (dictionary_.Return("#GlobalNLS_MaxArmijioIterations", int_value))
		SetGlobalNLS_MaxArmijioIterations(int_value);

	if (dictionary_.Return("#GlobalNLS_Method", string_value))
		SetGlobalNLS_Method(string_value);

	if (dictionary_.Return("#TraditionalLink"))
		SetTraditionalLink(true);

	// Final Checks
	if (iOpenMP_ == false)
		nThreads_ = 1;
}

void OpenSMOKE_KPP_DataManager::AssignKinetics(const string folderName)
{
	nameFolderKinetics_ = folderName;
}

void OpenSMOKE_KPP_DataManager::AssignInputFolder(const string folderName)
{
	nameFolderInput_ = folderName;

	nameTopologyFile_    = nameFolderInput_ + "/CFDNetwork.bzz";
	nameFirstGuessFile_  = nameFolderInput_ + "/FirstGuess.bzz";
}

void OpenSMOKE_KPP_DataManager::AssignBackupFile(const string fileName)
{
	iBackup_ = true;
	nameInputBackupFile_ = fileName;
}

void OpenSMOKE_KPP_DataManager::AssignOutputFolder(const string folderName)
{
	nameFolderOutput_ = folderName;

    // Folders
    {
		string slash = "\\";
		#if LINUX_SO==1
			slash = "/";
		#endif

		string systemCommand = "mkdir " + nameFolderOutput_;
        system(systemCommand.c_str());

		systemCommand = "mkdir " + nameFolderOutput_ + slash + "Additional";
        system(systemCommand.c_str());

		systemCommand = "mkdir " + nameFolderOutput_ + slash + "Maps";
        system(systemCommand.c_str());

		systemCommand = "mkdir " + nameFolderOutput_ + slash + "Solution";
        system(systemCommand.c_str());

		systemCommand = "mkdir " + nameFolderOutput_ + slash + "Residuals";
        system(systemCommand.c_str());

		systemCommand = "mkdir " + nameFolderOutput_ + slash + "Backup";
        system(systemCommand.c_str());
    }

	nameInOutFile_						= nameFolderOutput_ + "/Additional/InputOutput.out";
	nameResidualFile_					= nameFolderOutput_ + "/Residuals/Residuals.out";
	nameResidualSpeciesFile_			= nameFolderOutput_ + "/Residuals/ResidualsSpecies.out";
	nameResidualReactorStatisticsFile_	= nameFolderOutput_ + "/Residuals/StatisticsReactor.out";
	nameResidualSpeciesStatisticsFile_	= nameFolderOutput_ + "/Residuals/StatisticsSpecies.out";
	nameOutputBackupFile_				= nameFolderOutput_ + "/Backup/Backup";
}

void OpenSMOKE_KPP_DataManager::AssignReactions(const string option)
{
	if (option == "on")			iReactions_ = true;
	else if (option == "off")	iReactions_ = false;
	else ErrorMessage("#Reactions option: { on || off } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignPredictorCorrector_DeferredConvection(const string option)
{
	if (option == "on")			predictorCorrector_DeferredConvection_ = true;
	else if (option == "off")	predictorCorrector_DeferredConvection_ = false;
	else ErrorMessage("#PredictorCorrector_DeferredConvection option: { on || off } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignPredictorCorrector_AlgebraicConstraint(const string option)
{
	if (option == "on")			predictorCorrector_AlgebraicConstraint_ = true;
	else if (option == "off")	predictorCorrector_AlgebraicConstraint_ = false;
	else ErrorMessage("#PredictorCorrector_AlgebraicConstraint option: { on || off } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignOpenMP(const string option)
{
	if (option == "on")			iOpenMP_ = true;
	else if (option == "off")	iOpenMP_ = false;
	else ErrorMessage("#OpenMP option: { on || off } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignSaveKineticConstants(const string option)
{
	if (option == "on")			iSaveKineticConstants_ = true;
	else if (option == "off")	iSaveKineticConstants_ = false;
	else ErrorMessage("#SaveKineticConstants option: { on || off } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignSequenceCSTR_UpdatingNetwork(const string option)
{
	if (option == "continous")		sequenceCSTR_UpdatingNetwork_ = true;
	else if (option == "discrete")	sequenceCSTR_UpdatingNetwork_ = false;
	else ErrorMessage("#SequenceCSTR_UpdatingNetwork option: { continous || discrete } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignCorrection(const string correctionType)
{
	if (correctionType == "none")		correction_ = KPP_CORRECTION_NONE;
	else if (correctionType == "dirac")	correction_ = KPP_CORRECTION_DIRAC;
	else if (correctionType == "beta")	correction_ = KPP_CORRECTION_BETA;
	else if (correctionType == "sin")	correction_ = KPP_CORRECTION_SIN;
	else if (correctionType == "gauss")	correction_ = KPP_CORRECTION_GAUSS;
	else ErrorMessage("#Reactions option: { none || dirac || beta || sin || gauss } values are available!");
}

void OpenSMOKE_KPP_DataManager::SetMinimumTemperature(const double value, const string units)
{
	temperatureMin_ =  OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_KPP_DataManager::SetMaximumTemperature(const double value, const string units)
{
	temperatureMax_ = OpenSMOKE_Conversions::conversion_temperature(value, units);
}

void OpenSMOKE_KPP_DataManager::SetSingleReactor_IntegrationTime(const double value, const string units)
{
	singleReactor_IntegrationTime_ = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_KPP_DataManager::SetSequenceCSTR_MaxIterations(const int value)
{
	sequenceCSTR_MaxIterations_ = value;
}

void OpenSMOKE_KPP_DataManager::AssignPredictorCorrector_InitialTimeStep(const double value, const string units)
{
	predictorCorrector_InitialTimeStep_ = OpenSMOKE_Conversions::conversion_time(value, units);
	predictorCorrector_MaxTimeStep_ = 10.*predictorCorrector_InitialTimeStep_;
}

void OpenSMOKE_KPP_DataManager::SetPredictorCorrector_MaxTimeStep(const double value, const string units)
{
	predictorCorrector_MaxTimeStep_ = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_KPP_DataManager::SetPredictorCorrector_TimeStepIncrementFactor(const double value)
{
	predictorCorrector_TimeStepIncrementFactor_ = value;
}

void OpenSMOKE_KPP_DataManager::SetPredictorCorrector_TimeStepReductionFactor(const double value)
{
	predictorCorrector_TimeStepReductionFactor_ = value;
}

void OpenSMOKE_KPP_DataManager::SetPredictorCorrector_MaxIterations(const int value)
{
	predictorCorrector_MaxIterations_ = value;
}

void OpenSMOKE_KPP_DataManager::SetSingleReactor_OdeMaxJacobian(const int value)
{
	singleReactor_OdeMaxJacobian_ = value;
}

void OpenSMOKE_KPP_DataManager::SetSingleReactor_OdeStopResiduals(const double value)
{
	singleReactor_OdeStopResiduals_ = value;
}

void OpenSMOKE_KPP_DataManager::SetSequenceCSTR_VerboseStatistics(const string option)
{
	if (option == "on")			sequenceCSTR_VerboseStatistics_ = true;
	else if (option == "off")	sequenceCSTR_VerboseStatistics_  = false;
	else ErrorMessage("#SequenceCSTR_VerboseStatistics option: { on || off } values are available!");

}	

void OpenSMOKE_KPP_DataManager::AssignMaps(const vector<string> speciesNames)
{
	mapsSpeciesNames_ = speciesNames;
}

void OpenSMOKE_KPP_DataManager::AssignGlobalODE_SparseLinearSolver(const string solverType)
{
	if (solverType == "pardiso")			globalODE_SparseLinearSolver_ = KPP_SPARSESOLVER_PARDISO;
	else if (solverType == "mumps")			globalODE_SparseLinearSolver_ = KPP_SPARSESOLVER_MUMPS;
	else if (solverType == "lis")			globalODE_SparseLinearSolver_ = KPP_SPARSESOLVER_LIS;
	else if (solverType == "gaussSiedel")	globalODE_SparseLinearSolver_ = KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL;
	else ErrorMessage("#GlobalODE_SparseLinearSolver option: { pardiso || mumps || lis || gaussSiedel} values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignGlobalNLS_SparseLinearSolver(const string solverType)
{
		 if (solverType == "pardiso")		globalNLS_SparseLinearSolver_ = KPP_SPARSESOLVER_PARDISO;
	else if (solverType == "mumps")			globalNLS_SparseLinearSolver_ = KPP_SPARSESOLVER_MUMPS;
	else if (solverType == "lis")			globalNLS_SparseLinearSolver_ = KPP_SPARSESOLVER_LIS;
	else if (solverType == "gaussSiedel")	globalNLS_SparseLinearSolver_ = KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL;
	else ErrorMessage("#GlobalNLS_SparseLinearSolver option: { pardiso || mumps || lis || gaussSiedel } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignPredictorCorrector_SparseLinearSolver(const string solverType)
{
	     if (solverType == "pardiso")	predictorCorrector_SparseLinearSolver_ = KPP_SPARSESOLVER_PARDISO;
	else if (solverType == "mumps")		predictorCorrector_SparseLinearSolver_ = KPP_SPARSESOLVER_MUMPS;
	else if (solverType == "lis")		predictorCorrector_SparseLinearSolver_ = KPP_SPARSESOLVER_LIS;
	else ErrorMessage("#PredictorCorrector_SparseLinearSolver option: { pardiso || mumps || lis } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignSparseLinearSolverMassFlowRate(const string solverType)
{
	     if (solverType == "pardiso")	massFlowRate_SparseLinearSolver_ = KPP_SPARSESOLVER_PARDISO;
	else if (solverType == "mumps")		massFlowRate_SparseLinearSolver_ = KPP_SPARSESOLVER_MUMPS;
	else if (solverType == "bzz")		massFlowRate_SparseLinearSolver_ = KPP_SPARSESOLVER_BZZ;
	else ErrorMessage("#MassFlowRate_SparseLinearSolver_ option: { pardiso || mumps || bzz } values are available!");
}

void OpenSMOKE_KPP_DataManager::SetFillInReducingOrdering(const string fillInReducingOrdering)
{
	if (fillInReducingOrdering == "METIS")			fillInReducingOrdering_ = PARDISO_FILLIN_METIS;
	else if (fillInReducingOrdering == "MDA")		fillInReducingOrdering_ = PARDISO_FILLIN_MDA;
	else if (fillInReducingOrdering == "OpenMP")	fillInReducingOrdering_ = PARDISO_FILLIN_OPENMP;
	else ErrorMessage("#FillInReducingOrdering option: { METIS || MDA || OpenMP } values are available!");
}

void OpenSMOKE_KPP_DataManager::AssignKineticMaps(const string option)
{
	if (option == "on")			iKineticMaps_ = true;
	else if (option == "off")	iKineticMaps_ = false;
	else ErrorMessage("#KineticMaps option: { on || off } values are available!");
}
void OpenSMOKE_KPP_DataManager::AssignSymbolicKinetics(const string option)
{
	if (option == "on")			iSymbolicKinetics_ = true;
	else if (option == "off")	iSymbolicKinetics_ = false;
	else ErrorMessage("#SymbolicKinetics option: { on || off } values are available!");
}

void OpenSMOKE_KPP_DataManager::SetSingleReactor_AbsoluteTolerance(const double absTolerance)
{
	if ( (absTolerance < 1.e-12) && (absTolerance > 1.e-16) )
		singleReactor_OdeAbsoluteTolerance_ = absTolerance;
	else ErrorMessage("#SingleReactor_AbsoluteTolerance outside ranges: [1.e-16;1.e-12]");
}

void OpenSMOKE_KPP_DataManager::SetSingleReactor_RelativeTolerance(const double relTolerance)
{
	if ( (relTolerance < 1.e-5) && (relTolerance > 1.e-10) )
		singleReactor_OdeRelativeTolerance_ = relTolerance;
	else ErrorMessage("#SingleReactor_RelativeTolerance outside ranges: [1.e-10;1.e-5]");
}

void OpenSMOKE_KPP_DataManager::AssignThreads(const int value)
{
	nThreads_ = value;
}

void OpenSMOKE_KPP_DataManager::SetNetworkStatus(const KPP_Network_Status status)
{
	networkStatus_ = status;
}

void OpenSMOKE_KPP_DataManager::CheckUserInput()
{
	if (predictorCorrector_AlgebraicConstraint_ == true && predictorCorrector_DeferredConvection_ == false)
		ErrorMessage("#PredictorCorrector_AlgebraicConstraint option requires #PredictorCorrector_DeferredConvection is on...");
}

void OpenSMOKE_KPP_DataManager::SetPredictorCorrector_MultiTimeSplitting(const string option)
{
	if (option == "on")			predictorCorrector_MultiTimeSplitting_ = true;
	else if (option == "off")	predictorCorrector_MultiTimeSplitting_ = false;
	else ErrorMessage("#PredictorCorrector_MultiTimeSplitting option: { on || off } values are available!");
}
	
void OpenSMOKE_KPP_DataManager::SetPredictorCorrector_CourantCorrectionCoefficient(const double value)
{
	predictorCorrector_CourantCorrectionCoefficient_ = value;
}

void OpenSMOKE_KPP_DataManager::SetGlobalODE_AbsoluteTolerance(const double absTolerance)
{
	if ( (absTolerance <= 1.e-6) && (absTolerance >= 1.e-14) )
		globalODE_AbsoluteTolerance_ = absTolerance;
	else ErrorMessage("#GlobalODE_AbsoluteTolerance outside ranges: [1.e-14;1.e-6]");
}

void OpenSMOKE_KPP_DataManager::SetGlobalODE_RelativeTolerance(const double relTolerance)
{
	if ( (relTolerance <= 1.e-4) && (relTolerance >= 1.e-10) )
		globalODE_RelativeTolerance_ = relTolerance;
	else ErrorMessage("#GlobalODE_RelativeTolerance_ outside ranges: [1.e-10;1.e-4]");
}

void OpenSMOKE_KPP_DataManager::SetGlobalNLS_AbsoluteTolerance(const double absTolerance)
{
	if ( (absTolerance <= 1.e-8) && (absTolerance >= 1.e-16) )
		globalNLS_AbsoluteTolerance_ = absTolerance;
	else ErrorMessage("#GlobalNLS_AbsoluteTolerance outside ranges: [1.e-16;1.e-8]");
}

void OpenSMOKE_KPP_DataManager::SetGlobalNLS_RelativeTolerance(const double relTolerance)
{
	if ( (relTolerance <= 1.e-5) && (relTolerance >= 1.e-11) )
		globalNLS_RelativeTolerance_ = relTolerance;
	else ErrorMessage("#GlobalNLS_RelativeTolerance_ outside ranges: [1.e-11;1.e-5]");
}

void OpenSMOKE_KPP_DataManager::SetGlobalODE_MaxIterations(const int value)
{
	globalODE_MaxIterations_ = value;
}

void OpenSMOKE_KPP_DataManager::SetGlobalNLS_MaxIterations(const int value)
{
	globalNLS_MaxIterations_ = value;
}

void OpenSMOKE_KPP_DataManager::SetGlobalODE_TimeStepIncrementFactor(const double value)
{
	globalODE_TimeStepIncrementFactor_ = value;
}

void OpenSMOKE_KPP_DataManager::SetGlobalODE_TimeStepReductionFactor(const double value)
{
	globalODE_TimeStepReductionFactor_ = value;
}

void OpenSMOKE_KPP_DataManager::SetGlobalODE_InitialTimeStep(const double value, const string units)
{
	globalODE_InitialTimeStep_ = OpenSMOKE_Conversions::conversion_time(value, units);
	globalODE_MaxTimeStep_ = 100;
}

void OpenSMOKE_KPP_DataManager::SetGlobalODE_MaxTimeStep(const double value, const string units)
{
	globalODE_MaxTimeStep_ = OpenSMOKE_Conversions::conversion_time(value, units);
}

void OpenSMOKE_KPP_DataManager::SetGlobalODE_UpdatingFrequencyTimeStep(const int value)
{
	globalODE_UpdatingFrequencyTimeStep_ = value;
}

void OpenSMOKE_KPP_DataManager::SetGlobalNLS_MaxArmijioIterations(const int value)
{
	globalNLS_MaxArmijioIterations_ = value;
}

void OpenSMOKE_KPP_DataManager::SetGlobalNLS_Method(const string solverType)
{
	     if (solverType == "kelley")		globalNLS_Method_ = KPP_NLS_KELLEY;
	else if (solverType == "newton")		globalNLS_Method_ = KPP_NLS_NEWTON;
	else if (solverType == "chord")			globalNLS_Method_ = KPP_NLS_CHORD;
	else if (solverType == "shamanskii")	globalNLS_Method_ = KPP_NLS_SHAMANSKII;
	else ErrorMessage("#MassFlowRate_SparseLinearSolver_ option: { pardiso || mumps || bzz } values are available!");
}

void OpenSMOKE_KPP_DataManager::SetTraditionalLink(const bool value)
{
	iTraditionalLink_ = value;
}
