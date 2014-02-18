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

#ifndef OpenSMOKE_KPP_DataManager_H
#define OpenSMOKE_KPP_DataManager_H

#include "BzzMath.hpp"
#include "OpenSMOKE.hpp"
#include "OpenSMOKE_KPP_Definitions.h"
#include "linear_solvers/OpenSMOKE_PARDISO_Unsymmetric.h"

class OpenSMOKE_KPP_Dictionary;

class OpenSMOKE_KPP_DataManager
{
public:

	// Constructor-Destructor
	OpenSMOKE_KPP_DataManager(OpenSMOKE_KPP_Dictionary& dictionary);
	~OpenSMOKE_KPP_DataManager(void);

	// General
	inline bool iReactions() const { return iReactions_; };
	inline bool iKineticMaps() const { return iKineticMaps_; };
	inline bool iSymbolicKinetics() const { return iSymbolicKinetics_; };
	inline bool iSaveKineticConstants() const { return iSaveKineticConstants_; };
	inline bool iOpenMP() const { return iOpenMP_; };
	inline bool iBackup() const { return iBackup_; };
	inline int  nThreads() const { return nThreads_; };
	inline KPP_Network_Status networkStatus() { return networkStatus_; }
	inline vector<string>& mapsSpeciesNames() { return mapsSpeciesNames_; };
	inline string nameFolderOutput() const { return nameFolderOutput_; };
	inline string nameFolderInput() const { return nameFolderInput_; }
	inline string nameFolderKinetics() const { return nameFolderKinetics_; }
	inline string nameTopologyFile() const { return nameTopologyFile_; };
	inline string nameFirstGuessFile() const { return nameFirstGuessFile_; }
	inline string nameResidualFile() const { return nameResidualFile_; };
	inline string nameResidualSpeciesFile() const { return nameResidualSpeciesFile_; };
	inline string nameResidualReactorStatisticsFile() const { return nameResidualReactorStatisticsFile_; };
	inline string nameResidualSpeciesStatisticsFile() const { return nameResidualSpeciesStatisticsFile_; };
	inline string nameInOutFile() const { return nameInOutFile_; }
	inline string nameOutputBackupFile() const { return nameOutputBackupFile_; }
	inline string nameInputBackupFile() const { return nameInputBackupFile_; }
	inline PARDISO_FillInReducingOrdering FillInReducingOrdering() const { return fillInReducingOrdering_; }
	inline KPP_Correction correction() const { return correction_; }
	inline KPP_SparseLinearSolver sparseLinearSolverMassFlowRate() const {return massFlowRate_SparseLinearSolver_;}
	inline bool TraditionalLink() const { return iTraditionalLink_; }
	inline double TemperatureMin() const { return temperatureMin_; }
	inline double TemperatureMax() const { return temperatureMax_; }

	// Single reactor
	inline double SingleReactor_OdeRelativeTolerance()	const { return singleReactor_OdeRelativeTolerance_; }
	inline double SingleReactor_OdeAbsoluteTolerance()	const { return singleReactor_OdeAbsoluteTolerance_; }
	inline double SingleReactor_IntegrationTime()		const { return singleReactor_IntegrationTime_; }
	inline int	  SingleReactor_OdeMaxJacobian()		const { return singleReactor_OdeMaxJacobian_; }
	inline double SingleReactor_OdeStopResiduals()		const { return singleReactor_OdeStopResiduals_; }

	// CSTR Sequence
	inline int	  SequenceCSTR_MaxIterations()			const { return sequenceCSTR_MaxIterations_; }
	inline bool   SequenceCSTR_UpdatingNetwork()		const { return sequenceCSTR_UpdatingNetwork_; };
	inline bool	  SequenceCSTR_VerboseStatistics()		const { return sequenceCSTR_VerboseStatistics_; }

	// Predictor Corrector
	inline KPP_SparseLinearSolver PredictorCorrector_SparseLinearSolver() const {return predictorCorrector_SparseLinearSolver_;}
	inline int    PredictorCorrector_MaxIterations() const { return predictorCorrector_MaxIterations_; };
	inline double PredictorCorrector_InitialTimeStep() const { return predictorCorrector_InitialTimeStep_; };
	inline double PredictorCorrector_MaxTimeStep() const { return predictorCorrector_MaxTimeStep_; };
	inline double PredictorCorrector_TimeStepIncrementFactor() const { return predictorCorrector_TimeStepIncrementFactor_; };
	inline double PredictorCorrector_TimeStepReductionFactor() const { return predictorCorrector_TimeStepReductionFactor_; };
	inline double PredictorCorrector_CourantCorrectionCoefficient() const { return predictorCorrector_CourantCorrectionCoefficient_; };
	inline bool	  PredictorCorrector_MultiTimeSplitting() const { return predictorCorrector_MultiTimeSplitting_; };
	inline bool	  PredictorCorrector_DeferredConvection() const { return predictorCorrector_DeferredConvection_; };
	inline bool   PredictorCorrector_AlgebraicConstraint() const { return predictorCorrector_AlgebraicConstraint_; };

	// Global ODE
	inline KPP_SparseLinearSolver GlobalODE_SparseLinearSolver() const {return globalODE_SparseLinearSolver_; }
	inline double GlobalODE_SafetyReductionCoefficient() const { return globalODE_SafetyReductionCoefficient; }
	inline int    GlobalODE_MaxIterations() const { return globalODE_MaxIterations_; };
	inline int    GlobalODE_UpdatingFrequencyTimeStep() const { return globalODE_UpdatingFrequencyTimeStep_; };
	inline double GlobalODE_InitialTimeStep() const { return globalODE_InitialTimeStep_; };
	inline double GlobalODE_MaxTimeStep() const { return globalODE_MaxTimeStep_; };
	inline double GlobalODE_TimeStepIncrementFactor() const { return globalODE_TimeStepIncrementFactor_; };
	inline double GlobalODE_TimeStepReductionFactor() const { return globalODE_TimeStepReductionFactor_; };
	inline double GlobalODE_RelativeTolerance()	const { return globalODE_RelativeTolerance_; }
	inline double GlobalODE_AbsoluteTolerance()	const { return globalODE_AbsoluteTolerance_; }

	// Global NLS
	inline KPP_SparseLinearSolver GlobalNLS_SparseLinearSolver() const {return globalNLS_SparseLinearSolver_;}
	inline KPP_NonLinearSystem_Method GlobalNLS_Method() const { return globalNLS_Method_; }
	inline double GlobalNLS_SafetyReductionCoefficient() const { return globalNLS_SafetyReductionCoefficient; }
	inline double GlobalNLS_RelativeTolerance()	const { return globalNLS_RelativeTolerance_; }
	inline double GlobalNLS_AbsoluteTolerance()	const { return globalNLS_AbsoluteTolerance_; }
	inline int    GlobalNLS_MaxIterations() const { return globalNLS_MaxIterations_; };
	inline int    GlobalNLS_MaxArmijioIterations() const { return globalNLS_MaxArmijioIterations_; };

public:

	// General
	void AssignKinetics(const string folderName);
	void AssignInputFolder(const string folderName);
	void AssignOutputFolder(const string folderName);
	void AssignBackupFile(const string fileName);
	void AssignReactions(const string);
	void AssignOpenMP(const string);
	void AssignThreads(const int);
	void AssignSaveKineticConstants(const string);
	void AssignCorrection(const string correctionType);
	void AssignMaps(const vector<string> speciesNames);
	void AssignSparseLinearSolverMassFlowRate(const string solverType);
	void AssignKineticMaps(const string string_value);
	void AssignSymbolicKinetics(const string string_value);
	void SetNetworkStatus(const KPP_Network_Status status);
	void SetFillInReducingOrdering(const string fillInReducingOrdering);


	// Sequence options
	void AssignSequenceCSTR_UpdatingNetwork(const string);
	void SetSequenceCSTR_MaxIterations(const int value);
	void SetSequenceCSTR_VerboseStatistics(const string value);


	// Single Reactor
	void SetSingleReactor_IntegrationTime(const double value, const string units);
	void SetSingleReactor_AbsoluteTolerance(const double absTolerance);
	void SetSingleReactor_RelativeTolerance(const double relTolerance);
	void SetSingleReactor_OdeMaxJacobian(const int value);
	void SetSingleReactor_OdeStopResiduals(const double value);
	

	// Predictor Corrector
	void AssignPredictorCorrector_DeferredConvection(const string);
	void AssignPredictorCorrector_AlgebraicConstraint(const string);
	void AssignPredictorCorrector_InitialTimeStep(const double value, const string units);
	void AssignPredictorCorrector_SparseLinearSolver(const string solverType);
	void SetPredictorCorrector_MaxTimeStep(const double value, const string units);
	void SetPredictorCorrector_MaxIterations(const int value);
	void SetPredictorCorrector_TimeStepIncrementFactor(const double value);
	void SetPredictorCorrector_TimeStepReductionFactor(const double value);
	void SetPredictorCorrector_MultiTimeSplitting(const string option);
	void SetPredictorCorrector_CourantCorrectionCoefficient(const double value);

	// Global ODE
	void AssignGlobalODE_SparseLinearSolver(const string solverType);
	void SetGlobalODE_AbsoluteTolerance(const double absTolerance);
	void SetGlobalODE_RelativeTolerance(const double relTolerance);
	void SetGlobalODE_MaxIterations(const int value);
	void SetGlobalODE_MaxTimeStep(const double value, const string units);
	void SetGlobalODE_InitialTimeStep(const double value, const string units);
	void SetGlobalODE_TimeStepIncrementFactor(const double value);
	void SetGlobalODE_TimeStepReductionFactor(const double value);
	void SetGlobalODE_UpdatingFrequencyTimeStep(const int value);
	void SetGlobalODE_SafetyReductionCoefficient(const double value);

	// Global NLS
	void AssignGlobalNLS_SparseLinearSolver(const string solverType);
	void SetGlobalNLS_Method(const string methodType);
	void SetGlobalNLS_AbsoluteTolerance(const double absTolerance);
	void SetGlobalNLS_RelativeTolerance(const double relTolerance);
	void SetGlobalNLS_MaxIterations(const int value);
	void SetGlobalNLS_MaxArmijioIterations(const int value);
	void SetGlobalNLS_SafetyReductionCoefficient(const double value);

	void SetTraditionalLink(const bool value);
	void SetMinimumTemperature(const double value, const string units);
	void SetMaximumTemperature(const double value, const string units);


private:

	void CheckDictionary();
	void CheckUserInput();

private:

	OpenSMOKE_KPP_Dictionary& dictionary_;

	// General
	KPP_Correction correction_;
	KPP_SparseLinearSolver massFlowRate_SparseLinearSolver_;
	bool iReactions_;
	bool iKineticMaps_;
	bool iSymbolicKinetics_;
	bool iSaveKineticConstants_;
	bool iOpenMP_;
	bool iBackup_;
	int  nThreads_;
	KPP_Network_Status networkStatus_;
	PARDISO_FillInReducingOrdering fillInReducingOrdering_;

	// Single Reactor
	double	singleReactor_OdeRelativeTolerance_;
	double	singleReactor_OdeAbsoluteTolerance_;
	double	singleReactor_IntegrationTime_;
	int		singleReactor_OdeMaxJacobian_;
	double	singleReactor_OdeStopResiduals_;

	// Sequence
	bool	sequenceCSTR_VerboseStatistics_;
	int		sequenceCSTR_MaxIterations_;
	bool	sequenceCSTR_UpdatingNetwork_;

	// Predictor-Corrector
	KPP_SparseLinearSolver predictorCorrector_SparseLinearSolver_;
	int		predictorCorrector_MaxIterations_;
	double	predictorCorrector_InitialTimeStep_;
	double	predictorCorrector_MaxTimeStep_;
	double	predictorCorrector_TimeStepIncrementFactor_;
	double	predictorCorrector_TimeStepReductionFactor_;
	bool	predictorCorrector_DeferredConvection_;
	bool	predictorCorrector_AlgebraicConstraint_;
	bool	predictorCorrector_MultiTimeSplitting_;
	double  predictorCorrector_CourantCorrectionCoefficient_;

	// Global ODE
	KPP_SparseLinearSolver globalODE_SparseLinearSolver_;
	double globalODE_SafetyReductionCoefficient;
	double	globalODE_RelativeTolerance_;
	double	globalODE_AbsoluteTolerance_;
	int		globalODE_MaxIterations_;
	double	globalODE_TimeStepIncrementFactor_;
	double	globalODE_TimeStepReductionFactor_;
	double  globalODE_InitialTimeStep_;
	double  globalODE_MaxTimeStep_; 
	int		globalODE_UpdatingFrequencyTimeStep_;

	
	// Global NLS
	KPP_SparseLinearSolver globalNLS_SparseLinearSolver_;
	KPP_NonLinearSystem_Method globalNLS_Method_;
	double	globalNLS_SafetyReductionCoefficient;
	double	globalNLS_RelativeTolerance_;
	double	globalNLS_AbsoluteTolerance_;
	int		globalNLS_MaxIterations_;
	int		globalNLS_MaxArmijioIterations_;

	// Names
	string nameFolderOutput_;
	string nameFolderInput_;
	string nameFolderKinetics_;
	string nameTopologyFile_;
	string nameFirstGuessFile_;
	string nameInOutFile_;
	string nameResidualFile_;
	string nameResidualSpeciesFile_;
	string nameResidualReactorStatisticsFile_;
	string nameResidualSpeciesStatisticsFile_;
	string nameOutputBackupFile_;
	string nameInputBackupFile_;
	vector<string> mapsSpeciesNames_;

	bool iTraditionalLink_;
	double temperatureMin_;
	double temperatureMax_;

private:

	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
};

#endif	// OpenSMOKE_KPP_DataManager_H



