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

#include "OpenSMOKE_KPP_Dictionary.h"

OpenSMOKE_KPP_Dictionary::OpenSMOKE_KPP_Dictionary()
{
    SetupBase();
	SetName("OpenSMOKE_KPP_Dictionary");

	// Input data
	Add("#Kinetics",				'C', 'S', "Kinetic scheme");
	Add("#Input",					'C', 'S', "Input folder");
	Add("#Backup",					'O', 'S', "Backup file");
	Add("#Clustering",				'O', 'S', "Clustering: on || off");
        Add("#EnergyClustering",			'O', 'S', "Energy Clustering: on || off");

	// Output data
	Add("#Output",					'C', 'S', "Output folder");
	Add("#Maps",					'C', 'V', "List of species names");
	Add("#KineticMaps",				'C', 'S', "Kinetic maps: on || off");

	// Kinetics
	Add("#Reactions",				'C', 'S', "Chemical reactions: on || off");
	Add("#SymbolicKinetics",		'C', 'S', "Symbolic kinetics: on || off");
	Add("#Correction",				'C', 'S', "Kinetic correction: none || sin || dirac || gaussian || beta");

	// Numerical data
	Add("#MassFlowRate_SparseLinearSolver",	'C', 'S', "pardiso || mumps || bzz || lis");
	Add("#MassFlowRate_LisSolvingMethod",	'C', 'S', "Lis Solving Method (default = bicg)");
	Add("#MassFlowRate_LisPreconditioner",	'C', 'S', "Lis Preconditioner Method (default = ilu)");
	Add("#MassFlowRate_LisMaxIterations",	'C', 'S', "MaxIterations if solving with lis (default 1000)");
	Add("#MassFlowRate_LisConvCriteria",	'C', 'S', "Lis exit criteria (default = 1e-12)");
	Add("#FillInReducingOrdering",			'C', 'S', "METIS || MDA || OpenMP");
	Add("#Threads",							'C', 'I', "Number of threads");
	Add("#SaveKineticConstants",			'C', 'S', "Save Kinetic Constants: on || off");
	Add("#MinimumTemperature", 			'O', 'M', "Minimum allowed temperature for fluctuations");
	Add("#MaximumTemperature", 			'O', 'M', "Maximum allowed temperature for fluctuations");

	// CSTR Sequence
	Add("#SequenceCSTR_MaxIterations",		'C', 'I', "Sequence: Maximum number of steps (default 10000)");
	Add("#SequenceCSTR_UpdatingNetwork",	'C', 'S', "Sequence: updating network: discrete || continous");
	Add("#SequenceCSTR_VerboseStatistics",	'C', 'S', "Sequence: verbose statistics (default off)");

	// Single Reactor
	Add("#SingleReactor_IntegrationTime",	'C', 'M', "Single Reactor: integration time (default 100 s)");
	Add("#SingleReactor_OdeMaxJacobian",	'C', 'I', "Single Reactor: max number of Jacobians (default: not provided; suggested 20)");
	Add("#SingleReactor_OdeStopResiduals",	'C', 'D', "Single Reactor: reduction factor of ODE residuals before exit (default: not provided; suggested less than 1.e-6)");
	Add("#SingleReactor_RelativeTolerance",	'C', 'D', "Single Reactor: relative tolerance (default 1.2e-5)");
	Add("#SingleReactor_AbsoluteTolerance",	'C', 'D', "Single Reactor: absolute tolerance (default 1.e-10)");


	// GlobalODE
	Add("#GlobalODE_SparseLinearSolver",			'C', 'S', "Global ODE: pardiso || mumps || lis || gaussSiedel");
	Add("#GlobalODE_LisMaxIterations",			'C', 'S', "MaxIterations if solving with lis (default 1000)");
	Add("#GlobalODE_LisSolvingMethod",				'C', 'S', "Lis Solving Method (default = bicg)");
	Add("#GlobalODE_LisPreconditioner",			'C', 'S', "Lis Preconditioner Method (default = ilu)");
	Add("#GlobalODE_LisConvCriteria",			'C', 'S', "Lis exit criteria (default = 1e-12)");
	Add("#GlobalODE_MaxIterations",					'C', 'I', "Global ODE: Maximum number of iterations (default 10000)");

	Add("#GlobalODE_RelativeTolerance",				'C', 'D', "Global ODE: relative tolerance (default 1.2e-5)");
	Add("#GlobalODE_AbsoluteTolerance",				'C', 'D', "Global ODE: absolute tolerance (default 1.e-10)");

	Add("#GlobalODE_InitialTimeStep",					'C', 'M', "Global ODE: Initial time step");
	Add("#GlobalODE_MaxTimeStep",						'C', 'M', "Global ODE: Max time step (default 1e5*InitialTimeStep)");
	Add("#GlobalODE_TimeStepIncrementFactor",			'C', 'D', "Global ODE: Time step increment factor (default 1.50)");
	Add("#GlobalODE_TimeStepReductionFactor",			'C', 'D', "Global ODE: Time step reduction factor (default 0.80)");
	Add("#GlobalODE_UpdatingFrequencyTimeStep",			'C', 'I', "Global ODE: Updating frequency time step (default 3)");


	// GlobalNLS
	Add("#GlobalNLS_SparseLinearSolver",	'C', 'S', "Global NLS: pardiso || mumps || lis || gaussSiedel");
	Add("#GlobalNLS_LisSolvingMethod",	'C', 'S', "Lis solving method (default = bicg)");
	Add("#GlobalNLS_LisPreconditioner",	'C', 'S', "Lis Preconditioner method (default = ilu)");
	Add("#GlobalNLS_LisConvCriteria",	'C', 'S', "Lis exit criteria (default = 1e-12)");
	Add("#GlobalNLS_LisMaxIterations",	'C', 'S', "Lis Max Iterations (default = 1000)");
	Add("#GlobalNLS_Method",				'C', 'S', "Global NLS: kelley || newton || chord || shamanskii (default: kelley)");
	Add("#GlobalNLS_MaxIterations",			'C', 'I', "Global NLS: Maximum number of iterations (default 40)");
	Add("#GlobalNLS_MaxArmijioIterations",	'C', 'I', "Global NLS: Maximum number of Armijio iterations (default 20)");
	Add("#GlobalNLS_RelativeTolerance",		'C', 'D', "Global NLS: relative tolerance (default 1.2e-5)");
	Add("#GlobalNLS_AbsoluteTolerance",		'C', 'D', "Global NLS: absolute tolerance (default 1.e-10)");

	//Switch criteria
	Add("#NormInfConvergence",		'C', 'D', "NormInf Convergence Criteria to switch to global");
	Add("#Norm1Convergence",		'C', 'D', "Norm1 Convergence Criteria to switch to global");
	Add("#FirstLoopIterations",		'C', 'I', "Number of Iterations in the first solving loop");
	Add("#SecondLoopIterations",		'C', 'I', "Number of Iterations in the second solving loop");
	Add("#NumberOfSequenceCSTR",		'C', 'I', "Number of Minimum CSTR Iterations for convergence");

	Add("#FluctuatingSpecies",		'O', 'V', "List of fluctuating species");


    Lock();
}

OpenSMOKE_KPP_Dictionary::~OpenSMOKE_KPP_Dictionary(void)
{
}

std::string OpenSMOKE_KPP_Dictionary::kinetics()
{
	std::string folderName = "Kinetics";
	if (Return("#Kinetics", folderName))
		return folderName;
	return folderName;
}
