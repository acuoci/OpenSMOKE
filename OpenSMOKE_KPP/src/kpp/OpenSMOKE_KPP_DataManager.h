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
#include <mpi.h>

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
	inline bool iBackup() const { return iBackup_; };
	inline bool iTraditionalLink()	const { return iTraditionalLink_;};
        inline bool iEnergyClustering()	const { return iEnergyClustering_;};
	inline int  nThreads() const { return nThreads_; };
	inline int nprocs() const    { return MPI::COMM_WORLD.Get_size(); }
	inline int procrank() const    { return MPI::COMM_WORLD.Get_rank(); }
	inline KPP_Network_Status networkStatus() { return networkStatus_; }
	inline std::vector<std::string>& mapsSpeciesNames() { return mapsSpeciesNames_; };
	inline std::string nameFolderOutput() const { return nameFolderOutput_; };
	inline std::string nameFolderInput() const { return nameFolderInput_; }
	inline std::string nameFolderKinetics() const { return nameFolderKinetics_; }
	inline std::string nameTopologyFile() const { return nameTopologyFile_; };
	inline std::string nameFirstGuessFile() const { return nameFirstGuessFile_; }
	inline std::string nameResidualFile() const { return nameResidualFile_; };
	inline std::string nameResidualSpeciesFile() const { return nameResidualSpeciesFile_; };
	inline std::string nameResidualReactorStatisticsFile() const { return nameResidualReactorStatisticsFile_; };
	inline std::string nameResidualSpeciesStatisticsFile() const { return nameResidualSpeciesStatisticsFile_; };
	inline std::string nameInOutFile() const { return nameInOutFile_; }
	inline std::string nameOverallTimes() const { return nameOverallTimes_; }
	inline std::string nameOutputBackupFile() const { return nameOutputBackupFile_; }
	inline std::string nameInputBackupFile() const { return nameInputBackupFile_; }
	inline PARDISO_FillInReducingOrdering FillInReducingOrdering() const { return fillInReducingOrdering_; }
	inline KPP_Correction correction() const { return correction_; }
	inline KPP_SparseLinearSolver sparseLinearSolverMassFlowRate() const {return massFlowRate_SparseLinearSolver_;}
	inline std::string massFlowRate_LisSolvingMethod()			const {return massFlowRate_LisSolvingMethod_; }
	inline std::string massFlowRate_LisPreconditioner()			const {return massFlowRate_LisPreconditioner_; }
	inline std::string massFlowRate_LisMaxIterations()			const {return massFlowRate_LisMaxIterations_; }
	inline std::string massFlowRate_LisConvCriteria()			const {return massFlowRate_LisConvCriteria_; }
	inline std::vector<std::string>& fluctuatingSpecies() { return fluctuatingSpecies_; };

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
	inline std::string GlobalODE_LisSolvingMethod()			const {return globalODE_LisSolvingMethod_; }
	inline std::string GlobalODE_LisPreconditioner()			const {return globalODE_LisPreconditioner_; }
	inline std::string GlobalODE_LisMaxIterations()			const {return globalODE_LisMaxIterations_; }
	inline std::string GlobalODE_LisConvCriteria()			const {return globalODE_LisConvCriteria_; }
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
	inline std::string GlobalNLS_LisSolvingMethod()			const {return globalNLS_LisSolvingMethod_; }
	inline std::string GlobalNLS_LisPreconditioner()			const {return globalNLS_LisPreconditioner_; }
	inline std::string GlobalNLS_LisMaxIterations()			const {return globalNLS_LisMaxIterations_; }
	inline std::string GlobalNLS_LisConvCriteria()			const {return globalNLS_LisConvCriteria_; }
	inline KPP_NonLinearSystem_Method GlobalNLS_Method() const { return globalNLS_Method_; }
	inline double GlobalNLS_SafetyReductionCoefficient() const { return globalNLS_SafetyReductionCoefficient; }
	inline double GlobalNLS_RelativeTolerance()	const { return globalNLS_RelativeTolerance_; }
	inline double GlobalNLS_AbsoluteTolerance()	const { return globalNLS_AbsoluteTolerance_; }
	inline int    GlobalNLS_MaxIterations() const { return globalNLS_MaxIterations_; };
	inline int    GlobalNLS_MaxArmijioIterations() const { return globalNLS_MaxArmijioIterations_; };

	inline double TemperatureMin() const { return temperatureMin_;}
	inline double TemperatureMax() const { return temperatureMax_;}

	//Switch Criteria
	inline double 	NormInfConvergence()		const { return NormInfConvergence_; }
	inline double 	Norm1Convergence()		const { return Norm1Convergence_; }
	inline int	FirstLoopIterations()		const { return FirstLoopIterations_; }
	inline int	SecondLoopIterations()		const { return SecondLoopIterations_; }
	inline int	MinimumCSTRIterations()		const { return MinimumCSTRIterations_; }

public:

	// General
	void AssignKinetics(const std::string folderName);
	void AssignInputFolder(const std::string folderName);
	void AssignOutputFolder(const std::string folderName);
	void AssignBackupFile(const std::string fileName);
	void AssignClustering(const std::string options);
        void AssignEnergyClustering(const std::string options);
	void AssignReactions(const std::string);
	void AssignThreads(const int);
	void AssignSaveKineticConstants(const std::string);
	void AssignCorrection(const std::string correctionType);
	void AssignMaps(const std::vector<std::string> speciesNames);
	void AssignSparseLinearSolverMassFlowRate(const std::string solverType);
	void AssignLISMassFlowRateSolvingMethod(const std::string solver);
	void AssignLISMassFlowRateConvCriteria(const std::string residual);
	void AssignLISMassFlowRatePreconditioner(const std::string preconditioner);
	void AssignLISMassFlowRateMaxIterations(const std::string maxiter);
	void AssignKineticMaps(const std::string string_value);
	void AssignSymbolicKinetics(const std::string string_value);
	void SetNetworkStatus(const KPP_Network_Status status);
	void SetFillInReducingOrdering(const std::string fillInReducingOrdering);
	void SetFluctuatingSpecies(const std::vector<std::string> speciesNames);


	// Sequence options
	void AssignSequenceCSTR_UpdatingNetwork(const std::string);
	void SetSequenceCSTR_MaxIterations(const int value);
	void SetSequenceCSTR_VerboseStatistics(const std::string value);


	// Single Reactor
	void SetSingleReactor_IntegrationTime(const double value, const std::string units);
	void SetSingleReactor_AbsoluteTolerance(const double absTolerance);
	void SetSingleReactor_RelativeTolerance(const double relTolerance);
	void SetSingleReactor_OdeMaxJacobian(const int value);
	void SetSingleReactor_OdeStopResiduals(const double value);
	

	// Predictor Corrector
	void AssignPredictorCorrector_DeferredConvection(const std::string);
	void AssignPredictorCorrector_AlgebraicConstraint(const std::string);
	void AssignPredictorCorrector_InitialTimeStep(const double value, const std::string units);
	void AssignPredictorCorrector_SparseLinearSolver(const std::string solverType);
	void SetPredictorCorrector_MaxTimeStep(const double value, const std::string units);
	void SetPredictorCorrector_MaxIterations(const int value);
	void SetPredictorCorrector_TimeStepIncrementFactor(const double value);
	void SetPredictorCorrector_TimeStepReductionFactor(const double value);
	void SetPredictorCorrector_MultiTimeSplitting(const std::string option);
	void SetPredictorCorrector_CourantCorrectionCoefficient(const double value);

	// Global ODE
	void AssignGlobalODE_SparseLinearSolver(const std::string solverType);
	void AssignGlobalODE_LisSolvingMethod(const std::string option);
	void AssignGlobalODE_LisPreconditioner(const std::string option);
	void AssignGlobalODE_LisMaxIterations(const std::string option);
	void AssignGlobalODE_LisConvCriteria(const std::string option);
	void SetGlobalODE_AbsoluteTolerance(const double absTolerance);
	void SetGlobalODE_RelativeTolerance(const double relTolerance);
	void SetGlobalODE_MaxIterations(const int value);
	void SetGlobalODE_MaxTimeStep(const double value, const std::string units);
	void SetGlobalODE_InitialTimeStep(const double value, const std::string units);
	void SetGlobalODE_TimeStepIncrementFactor(const double value);
	void SetGlobalODE_TimeStepReductionFactor(const double value);
	void SetGlobalODE_UpdatingFrequencyTimeStep(const int value);
	void SetGlobalODE_SafetyReductionCoefficient(const double value);

	// Global NLS
	void AssignGlobalNLS_SparseLinearSolver(const std::string solverType);
	void AssignGlobalNLS_LisSolvingMethod(const std::string option);
	void AssignGlobalNLS_LisPreconditioner(const std::string option);
	void AssignGlobalNLS_LisMaxIterations(const std::string option);
	void AssignGlobalNLS_LisConvCriteria(const std::string option);
	void SetGlobalNLS_Method(const std::string methodType);
	void SetGlobalNLS_AbsoluteTolerance(const double absTolerance);
	void SetGlobalNLS_RelativeTolerance(const double relTolerance);
	void SetGlobalNLS_MaxIterations(const int value);
	void SetGlobalNLS_MaxArmijioIterations(const int value);
	void SetGlobalNLS_SafetyReductionCoefficient(const double value);

	void SetMinimumTemperature(const double value, const std::string units);
	void SetMaximumTemperature(const double value, const std::string units);

	//Switch Criteria
	void SetNormInf_Convergence(const double value);
	void SetNorm1_Convergence(const double value);
	void SetFirstLoop_Iterations(const int value);
	void SetSecondLoop_Iterations(const int value);
	void Set_CSTR_Iterations(const int value);


private:

	void CheckDictionary();
	void CheckUserInput();

private:

	OpenSMOKE_KPP_Dictionary& dictionary_;

	// General
	KPP_Correction correction_;
	KPP_SparseLinearSolver massFlowRate_SparseLinearSolver_;
	std::string	massFlowRate_LisSolvingMethod_;
	std::string	massFlowRate_LisPreconditioner_;
	std::string	massFlowRate_LisMaxIterations_;
	std::string	massFlowRate_LisConvCriteria_;
	bool iReactions_;
	bool iKineticMaps_;
	bool iSymbolicKinetics_;
	bool iSaveKineticConstants_;
	bool iBackup_;
	bool iTraditionalLink_;
        bool iEnergyClustering_;
	int  nThreads_;
	KPP_Network_Status networkStatus_;
	PARDISO_FillInReducingOrdering fillInReducingOrdering_;	
	std::vector<std::string> fluctuatingSpecies_;

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
	std::string	globalODE_LisSolvingMethod_;
	std::string	globalODE_LisPreconditioner_;
	std::string	globalODE_LisMaxIterations_;
	std::string	globalODE_LisConvCriteria_;
	double	globalODE_SafetyReductionCoefficient;
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
	std::string	globalNLS_LisSolvingMethod_;
	std::string	globalNLS_LisPreconditioner_;
	std::string	globalNLS_LisConvCriteria_;
	std::string	globalNLS_LisMaxIterations_;
	double	globalNLS_SafetyReductionCoefficient;
	double	globalNLS_RelativeTolerance_;
	double	globalNLS_AbsoluteTolerance_;
	int		globalNLS_MaxIterations_;
	int		globalNLS_MaxArmijioIterations_;

	// Switch Criteria
	double NormInfConvergence_;
	double Norm1Convergence_;
	int FirstLoopIterations_;
	int SecondLoopIterations_;
	int MinimumCSTRIterations_;

	// Names
	std::string nameFolderOutput_;
	std::string nameFolderInput_;
	std::string nameFolderKinetics_;
	std::string nameTopologyFile_;
	std::string nameFirstGuessFile_;
	std::string nameInOutFile_;
	std::string nameOverallTimes_;
	std::string nameResidualFile_;
	std::string nameResidualSpeciesFile_;
	std::string nameResidualReactorStatisticsFile_;
	std::string nameResidualSpeciesStatisticsFile_;
	std::string nameOutputBackupFile_;
	std::string nameInputBackupFile_;
	std::vector<std::string> mapsSpeciesNames_;

	double temperatureMin_;
	double temperatureMax_;

private:

	void ErrorMessage(const std::string message_);
	void WarningMessage(const std::string message_);
};

#endif	// OpenSMOKE_KPP_DataManager_H



