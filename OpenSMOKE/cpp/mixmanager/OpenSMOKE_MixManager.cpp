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
#include "mixmanager/OpenSMOKE_MixManager.h"
#include "basic/OpenSMOKE_Conversions.h"

void OpenSMOKE_MixManager::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_MixManager"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_MixManager::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_MixManager"		<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

void OpenSMOKE_MixManager::SetName(const string name)
{
	name_object = name;
}

void OpenSMOKE_MixManager::SetGlobalKinetics()
{
    iGlobalKinetics = true;
}

void OpenSMOKE_MixManager::UnsetGlobalKinetics()
{
    iGlobalKinetics = false;
}

void OpenSMOKE_MixManager::SetTwoEquationModel()
{
    iTwoEquationModel = true;
}

void OpenSMOKE_MixManager::UnsetTwoEquationModel()
{
    iTwoEquationModel = false;
}

OpenSMOKE_MixManager::OpenSMOKE_MixManager()
{
    // Control Variables
    assignedKineticScheme		= false;			//
    assignedGlobalKineticScheme	= false;			//
	assignedSoot2EModel			= false;			//
    assignedInletFlows			= false;			//
    
    // Properties Index
    iGlobalKinetics             = false;			// Global Kinetics 0=OFF 1=ON
    iTwoEquationModel	        = false;			// 2E Model for soot predictions

	iVerboseReactionRates		= false;			// Reaction rates
	iVerboseFormationRates		= false;			// Formation rates
	iVerboseROPA				= false;			// Rate of production analysis
	iVerboseSensitivity			= false;			// Sensitivity analysis
	
    // Output Folder
    outputFolderName			= "Output";
    outputName					=  outputFolderName + "/" + "Mix.out";
    outputSootName				=  outputFolderName + "/" + "Soot.out";
    outputPAHName				=  outputFolderName + "/" + "PAH.out";
    output2EName				=  outputFolderName + "/" + "2E.out";
	outputReactionRatesName		=  outputFolderName + "/" + "ReactionRates.out";
    outputFormationRatesName	=  outputFolderName + "/" + "FormationRates.out";
    outputROPAName				=  outputFolderName + "/" + "ROPA.out";

	// Name Object
	name_object = "[not assigned]";
}

void OpenSMOKE_MixManager::AssignKineticScheme(OpenSMOKE_ReactingGas &_mix)
{
    mix     = &_mix;

    NC = mix->NumberOfSpecies();
    NR = mix->NumberOfReactions();

    // These vectors are used for the evaluation of gas mixture properties
    ChangeDimensions(NC, &omega);
    ChangeDimensions(NC, &x);
    ChangeDimensions(NC, &c);
    ChangeDimensions(NC, &R);
    ChangeDimensions(NC, &Dmix);
    ChangeDimensions(NR, &r);

    // Control Variables
    assignedKineticScheme = true;
}

void OpenSMOKE_MixManager::AssignInletFlows(OpenSMOKE_GasStream &_inletStream)
{
    if (assignedKineticScheme == false)
        ErrorMessage("You must define the kinetic scheme before defining the inlet/initial flows!");

    inletStream = &_inletStream;

    assignedInletFlows = true;
}

void OpenSMOKE_MixManager::AssignSoot2EModel(OpenSMOKE_2EModel &_soot2EModel)
{
    soot2EModel = &_soot2EModel;

    assignedSoot2EModel = true;
}

void OpenSMOKE_MixManager::AssignGlobalKineticScheme(OpenSMOKE_GlobalKinetics &_global)
{
    global = &_global;

    assignedGlobalKineticScheme = true;
}

void OpenSMOKE_MixManager::Initialize()
{
    T            = inletStream->T;
    P            = inletStream->P;
    omega        = inletStream->omega;
    massFlowRate = inletStream->massFlowRate;
}

void OpenSMOKE_MixManager::Properties()
{
	// Concentration and density
	// ----------------------------------------------------------------------
	mix->GetMWAndMoleFractionsFromMassFractions(MWtot, x, omega);
	cTot = P  / (Constants::R_J_kmol*T);
	rho = cTot * MWtot;
	c = cTot*x;

	// Properties
	// ----------------------------------------------------------------------
	mix->SpeciesCp(T);
	Cp = mix->MixCp_FromMassFractions(omega);

	mix->SpeciesConductivityFromFitting(T);
	lambda = mix->MixConductivity_FromMolarFractions(x);

	mix->SpeciesViscosityFromFitting(T);
	mu = mix->MixViscosity_FromMolarFractions(x);

	mix->SpeciesDiffusivityFromFitting(T, P/1.e5);
	mix->MixDiffusionCoefficient_FromMolarFractionsAndPM(Dmix,x);


	// Formation rates
	// ----------------------------------------------------------------------
	mix->ComputeKineticParameters( T, log(T), 1./T, P);
	mix->ComputeFromConcentrations( T, c, cTot, &R);		// [kmol/m3/s]
	ElementByElementProduct(R, mix->M, &R);					// [kg/m3/s]
	QReaction = - mix->ComputeQReaction(T);					// [J/m3.s]

	// Global Kinetics
	// ----------------------------------------------------------------------
	if (iGlobalKinetics == true)
	{
		global->GiveMeFormationRates(T,c,R);
		global->GiveMeReactionHeat(T, R, QReaction);
	}

	H_mass  = inletStream->massSpecificEnthalpy;			// [J/kg]
	U_mass  = inletStream->massSpecificInternalEnergy;		// [J/kg]
	S_mass  = inletStream->massSpecificEntropy;				// [J/kg/K]
}

void OpenSMOKE_MixManager::VideoSummary()
{
	cout << "Density:              " << rho		<< " [kg/m3]"	<< endl;
	cout << "Total concentration:  " << cTot	<< " [kmol/m3]" << endl;
	cout << "Molecular weight:     " << MWtot	<< " [kg/kmol]" << endl;
	cout << "Specific heat:        " << Cp		<< " [J/kg/K]"	<< endl;
	cout << "Thermal conductivity: " << lambda	<< " [W/m/K]"	<< endl;
	cout << "Dynamic Viscosity:    " << mu		<< " [Pa.s]"	<< endl;
	cout << "Kinematic Viscosity:  " << mu/rho	<< " [m2/s]"	<< endl;

	cout << "Enthalpy:             " << H_mass	<< " [J/kg]"	<< endl;
	cout << "Internal Energy:      " << U_mass	<< " [J/kg]"	<< endl;
	cout << "Entropy:              " << S_mass	<< " [J/kg/K]"	<< endl;
}
