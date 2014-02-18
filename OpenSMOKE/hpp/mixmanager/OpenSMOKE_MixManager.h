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

#ifndef OPENSMOKE_MIXMANAGER_H
#define OPENSMOKE_MIXMANAGER_H

#include "OpenSMOKE.hpp"

class OpenSMOKE_Dictionary_MixManager : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_MixManager();
};

class OpenSMOKE_MixManager
{

public: // Virtual functions (compulsory)

	OpenSMOKE_MixManager();
	OpenSMOKE_GasStream stream;
	
	// 1. Compulsory
	void AssignSoot2EModel(OpenSMOKE_2EModel &_soot2EModel);
    void AssignGlobalKineticScheme(OpenSMOKE_GlobalKinetics &_global);
	void AssignKineticScheme(OpenSMOKE_ReactingGas &_mix);
    void AssignInletFlows(OpenSMOKE_GasStream &stream);

    // 2. Option settings
	void SetName(const string name);
    void SetGlobalKinetics();
	void SetTwoEquationModel();
    void UnsetGlobalKinetics();
	void UnsetTwoEquationModel();

	void Initialize();
	void Properties();
	void VideoSummary();

private:

	double T;
	double P;
	double cTot;
	double rho;
	double MWtot;
	double massFlowRate;
	double Cp;
	double lambda;
	double mu;
	double QReaction;

	double H_mass;
	double U_mass;
	double S_mass;

private:

	// Error and warning messages
	void	ErrorMessage(const string message);
	void	WarningMessage(const string message);

private: // Virtual functions (compulsory)
	
	int NC;
	int NR;

	BzzVector omega;
	BzzVector x;
	BzzVector c;
	BzzVector R;
	BzzVector Dmix;
	BzzVector r;

private: // Specific data

	// Error and warning messages
	string name_object;

	// Control Variables
	bool assignedKineticScheme;
	bool assignedGlobalKineticScheme;
	bool assignedSoot2EModel;
	bool assignedInletFlows;

	// Options
	bool iGlobalKinetics;
	bool iTwoEquationModel;

	bool iVerbose;
	bool iVerboseReactionRates;
	bool iVerboseFormationRates;
	bool iVerboseROPA;
	bool iVerboseSensitivity;

	// Print Options
	string outputFolderName;
	string outputName;
	string outputSootName;
	string outputPAHName;
	string output2EName;
	string outputReactionRatesName;
	string outputFormationRatesName;
	string outputROPAName;
	string outputSensitivityName;

private:

	OpenSMOKE_ReactingGas				*mix;						// Gas Mixture
	OpenSMOKE_GlobalKinetics			*global;					// Global Kinetics
	OpenSMOKE_2EModel					*soot2EModel;				// Two Equation Model
	OpenSMOKE_GasStream					*inletStream;				// Inlet stream
	OpenSMOKE_RateOfProductionAnalysis	 ROPA;						// Rate of production analysis
	OpenSMOKE_SensitivityAnalysis		 Sensitivity;				// Sensitivity analysis

};

#endif // OPENSMOKE_MIXMANAGER_H
