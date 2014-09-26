/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci   	                       *
 *   alberto.cuoci@polimi.it   						                       *
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

#ifndef OPENSMOKE_PREPROCESSORREACTINGGAS
#define OPENSMOKE_PREPROCESSORREACTINGGAS

#include "preprocessing/OpenSMOKE_PreProcessorIdealGas.h"
#include "preprocessing/OpenSMOKE_PreProcessorKinetics.h"
#include "engine/OpenSMOKE_ReactingGas.h"


class OpenSMOKE_PreProcessorReactingGas : public OpenSMOKE_PreProcessorIdealGas
{
public:

	OpenSMOKE_PreProcessorKinetics kinetics;

	// Main functions
	int		NumberOfReactions();				// Number of reactions
	void	PreProcessing(const std::string pathName, const std::string fileNameTransport, const std::string fileNameElements, const std::string option);
	void	PreProcessingTransportSensitivity(const std::string pathName, const std::string fileNameTransport, const std::string fileNameElements, const double eps);

	// Reaction rates
	BzzVector		r;						// Reaction rates [kmol/m3/s]	
	BzzVector		k1;						 
	BzzVector		k2;						
	BzzVector		uKeq;
	
	// Reaction names
	char **strReaction;

private:

	// Number of reactions
	int	NR;
};

enum fittingCoefficientExtraction{	OPENSMOKE_FITTING_BASE, 
									OPENSMOKE_FITTING_KE_PLUS, OPENSMOKE_FITTING_KE_MINUS,
									OPENSMOKE_FITTING_SIGMA_PLUS, OPENSMOKE_FITTING_SIGMA_MINUS,
									OPENSMOKE_FITTING_MU_PLUS, OPENSMOKE_FITTING_MU_MINUS,
									OPENSMOKE_FITTING_ALFA_PLUS, OPENSMOKE_FITTING_ALFA_MINUS,
									OPENSMOKE_FITTING_ZROT_PLUS, OPENSMOKE_FITTING_ZROT_MINUS };

enum fittingCoefficientSensitivityMode{	OPENSMOKE_FITTING_ALL, 
										OPENSMOKE_FITTING_ONLY_VISCOSITY, OPENSMOKE_FITTING_ONLY_CONDUCTIVITY,
										OPENSMOKE_FITTING_ONLY_DIFFUSIVITIES };

class OpenSMOKE_LennardJonesSensitivityCoefficientsManager
{
public:

	OpenSMOKE_LennardJonesSensitivityCoefficientsManager();
	void Setup(const std::string file_name, fittingCoefficientSensitivityMode _fittingMode);
	void SetName(const std::string name);

	void GetFittingCoefficients(const fittingCoefficientExtraction fitting, const int j, OpenSMOKE_ReactingGas &mix);

private:

	int NC;
	double epsilon;
	vector<string> names;

	fittingCoefficientSensitivityMode fittingMode;

public:

	BzzVector epsylon_over_kb;
	BzzVector sigma;
	BzzVector mu;
	BzzVector alfa;
	BzzVector zRot298;

	BzzVector delta_epsylon_over_kb;
	BzzVector delta_sigma;
	BzzVector delta_mu;
	BzzVector delta_alfa;
	BzzVector delta_zRot298;

private:

	// Fitting data
	BzzMatrix fittingEta_Base;
	BzzMatrix fittingLambda_Base;
	BzzMatrix fittingDbinary_Base;
	BzzMatrix fittingTetaBinary_Base;

	// Fitting data
	BzzMatrix fittingEta_KE_Plus;
	BzzMatrix fittingLambda_KE_Plus;
	BzzMatrix fittingDbinary_KE_Plus;
	BzzMatrix fittingTetaBinary_KE_Plus;

	// Fitting data
	BzzMatrix fittingEta_KE_Minus;
	BzzMatrix fittingLambda_KE_Minus;
	BzzMatrix fittingDbinary_KE_Minus;
	BzzMatrix fittingTetaBinary_KE_Minus;

	// Fitting data
	BzzMatrix fittingEta_Sigma_Plus;
	BzzMatrix fittingLambda_Sigma_Plus;
	BzzMatrix fittingDbinary_Sigma_Plus;
	BzzMatrix fittingTetaBinary_Sigma_Plus;

	// Fitting data
	BzzMatrix fittingEta_Sigma_Minus;
	BzzMatrix fittingLambda_Sigma_Minus;
	BzzMatrix fittingDbinary_Sigma_Minus;
	BzzMatrix fittingTetaBinary_Sigma_Minus;

	// Fitting data
	BzzMatrix fittingEta_Mu_Plus;
	BzzMatrix fittingLambda_Mu_Plus;
	BzzMatrix fittingDbinary_Mu_Plus;
	BzzMatrix fittingTetaBinary_Mu_Plus;

	// Fitting data
	BzzMatrix fittingEta_Mu_Minus;
	BzzMatrix fittingLambda_Mu_Minus;
	BzzMatrix fittingDbinary_Mu_Minus;
	BzzMatrix fittingTetaBinary_Mu_Minus;

	// Fitting data
	BzzMatrix fittingEta_Alfa_Plus;
	BzzMatrix fittingLambda_Alfa_Plus;
	BzzMatrix fittingDbinary_Alfa_Plus;
	BzzMatrix fittingTetaBinary_Alfa_Plus;

	// Fitting data
	BzzMatrix fittingEta_Alfa_Minus;
	BzzMatrix fittingLambda_Alfa_Minus;
	BzzMatrix fittingDbinary_Alfa_Minus;
	BzzMatrix fittingTetaBinary_Alfa_Minus;

	// Fitting data
	BzzMatrix fittingEta_zRot298_Plus;
	BzzMatrix fittingLambda_zRot298_Plus;
	BzzMatrix fittingDbinary_zRot298_Plus;
	BzzMatrix fittingTetaBinary_zRot298_Plus;

	// Fitting data
	BzzMatrix fittingEta_zRot298_Minus;
	BzzMatrix fittingLambda_zRot298_Minus;
	BzzMatrix fittingDbinary_zRot298_Minus;
	BzzMatrix fittingTetaBinary_zRot298_Minus;

private:

	std::string name_object;

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};

#endif	// OPENSMOKE_PREPROCESSORREACTINGGAS