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

#include "OpenSMOKE_KPP_SingleReactor_KineticsManager.h"
#include "symbolickinetics/fluent_drm22_polimi_nox/OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX.h"
#include "symbolickinetics/gri30/OpenSMOKE_SymbolicKinetics_GRI30.h"
#include "symbolickinetics/polimi_c1c3htnox_1101/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101.h"
#include "symbolickinetics/polimi_h2conox_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101.h"
#include "symbolickinetics/polimi_h2conox_nothermal_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101.h"
#include "symbolickinetics/polimi_h2conox_nonnh_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101.h"
#include "symbolickinetics/polimi_h2conox_non2o_1101/OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101.h"
#include "symbolickinetics/sandiego_avio/OpenSMOKE_SymbolicKinetics_SanDiego_AVIO.h"
#include "symbolickinetics/polimi_nc7_avio/OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO.h"
#include "symbolickinetics/polimi_c1c3htnox_avio_0702/OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702.h"

#include "OpenSMOKE_KPP_Definitions.h"
#include "distributions/OpenSMOKE_BetaDistribution.h"
#include "distributions/OpenSMOKE_ClippedGaussianAccurateDistribution.h"
#include "distributions/OpenSMOKE_DoubleDeltaDiracDistribution.h"
#include "distributions/OpenSMOKE_SinDistribution.h"
#include "distributions/OpenSMOKE_SinDistribution.h"


void WarningLargeCorrectionCoefficient(const double T, const double g, const int iReaction, const string stringReaction, const double EsuR, const double n, double &CoeffCorr, const double MaxCoeffCorr);
void WarningSmallCorrectionCoefficient(const double T, const double g, const int iReaction, const string stringReaction, const double EsuR, const double n, double &CoeffCorr);

void OpenSMOKE_KPP_SingleReactor_KineticsManager::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_SingleReactor_KineticsManager"	<< endl;
    cout << "Error: " << message_					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_SingleReactor_KineticsManager::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_SingleReactor_KineticsManager"	<< endl;
    cout << "Warning: " << message_						<< endl;
	cout << endl;
}

void CorrectionCoefficient( const KPP_Correction correction, OpenSMOKE_ReactingGas *mix, const double T, const double cTot, const double variance, const double Tmin, const double Tmax, const double MaxCoeffCorr, BzzVector &correction_k1, BzzVector &correction_uKeq, BzzVector &correction_k2)
{	
	const double R_CTOT = 0.0820578337034;
	const double UR_CTOT = 1./R_CTOT;
	double uT     = 1./T;
	double logT   = log(T);
	double uRT    = uT*UR_CTOT;
	double loguRT = log(uRT);


	double DeltaT = Tmax-Tmin;
	double csiMean = (T-Tmin)/(Tmax-Tmin);
	double g = 0.50*variance*T*T/DeltaT/DeltaT;



	OpenSMOKE_DoubleDeltaDiracDistribution 		DoubleDeltaDirac;
	OpenSMOKE_ClippedGaussianAccurateDistribution 	GaussDistribution;
	OpenSMOKE_SinDistribution 			SinDistribution;
	OpenSMOKE_BetaDistribution 			BetaDistribution;

	if       (correction == KPP_CORRECTION_DIRAC)
	{
		DoubleDeltaDirac.Set(csiMean, g+1.e-16);
	}
	else if (correction == KPP_CORRECTION_GAUSS)	
	{
		GaussDistribution.Set(csiMean, g+1.e-16);
	}
	else if (correction == KPP_CORRECTION_BETA)	
	{
		double alfa = csiMean*(csiMean*(1.-csiMean)/g-1.);
		double beta = (1.-csiMean)*(csiMean*(1.-csiMean)/g-1.);
		BetaDistribution.Set(alfa, beta);
	}

	for(int k=1;k<=mix->kinetics.numEquilibrium;k++)
	{
		int j = mix->kinetics.reactionWithEquil[k];

		double CoeffCorr = 1.;

		mix->uKeq[j] = exp(-mix->reactionDS[j] + mix->reactionDH[j] - loguRT * mix->kinetics.sumNuij[j]);
		{
			double	EsuR = -(mix->kinetics.E1[j] + T*mix->reactionDH[j]);
			double	n = mix->kinetics.beta1[j] + mix->kinetics.sumNuij[j];

			if (correction == KPP_CORRECTION_DIRAC)
				CoeffCorr = DoubleDeltaDirac.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
			else if (correction == KPP_CORRECTION_BETA)
				CoeffCorr = BetaDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
			else if (correction == KPP_CORRECTION_GAUSS)
				CoeffCorr = GaussDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
			else if (correction == KPP_CORRECTION_SIN)
				CoeffCorr = SinDistribution.CorrectionCoefficient(variance*T*T, EsuR, n, T);	

			if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(T, g, j, mix->strReaction[j], EsuR, n, CoeffCorr, MaxCoeffCorr);
			if (CoeffCorr < 0.)		WarningSmallCorrectionCoefficient(T, g, j, mix->strReaction[j], EsuR, n, CoeffCorr);
		}

		correction_uKeq[j] = CoeffCorr;

		//cout << T << " " << j << " " << CoeffCorr << " " << Tmax << " " << variance << " " << DeltaT << " " << csiMean << " " << g << " " << Tmin << endl; 
	}

	for(int j=1;j<=mix->NumberOfReactions();j++)
	{
		double CoeffCorr = 1.;
		mix->k1[j] = exp(mix->kinetics.k01[j] + mix->kinetics.beta1[j] * logT + mix->kinetics.E1[j] * uT);
		{
			double	EsuR = -mix->kinetics.E1[j];
			double	n = mix->kinetics.beta1[j];

			if (correction == KPP_CORRECTION_DIRAC)
				CoeffCorr = DoubleDeltaDirac.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
			else if (correction == KPP_CORRECTION_BETA)
				CoeffCorr = BetaDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
			else if (correction == KPP_CORRECTION_GAUSS)
				CoeffCorr = GaussDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
			else if (correction == KPP_CORRECTION_SIN)
				CoeffCorr = SinDistribution.CorrectionCoefficient(variance*T*T, EsuR, n, T);	

			if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(T, g, j, mix->strReaction[j], EsuR, n, CoeffCorr, MaxCoeffCorr);
			if (CoeffCorr < 0.)		WarningSmallCorrectionCoefficient(T, g, j, mix->strReaction[j], EsuR, n, CoeffCorr);
		}
		correction_k1[j] = CoeffCorr;

		

		if(mix->kinetics.jThirdBody[j] >= 2)
		{
			CoeffCorr = 1.;

			mix->k2[j] = exp(mix->kinetics.k02[j] + mix->kinetics.beta2[j] * logT + mix->kinetics.E2[j] * uT);
			{
				double	EsuR = -mix->kinetics.E2[j];
				double	n = mix->kinetics.beta2[j];

				if (correction == KPP_CORRECTION_DIRAC)
					CoeffCorr = DoubleDeltaDirac.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
				else if (correction == KPP_CORRECTION_BETA)
					CoeffCorr = BetaDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
				else if (correction == KPP_CORRECTION_GAUSS)
					CoeffCorr = GaussDistribution.ReactionCorrectionCoefficient(n, EsuR, T, Tmin, Tmax);
				else if (correction == KPP_CORRECTION_SIN)
					CoeffCorr = SinDistribution.CorrectionCoefficient(variance*T*T, EsuR, n, T);	

				if (CoeffCorr > MaxCoeffCorr)	WarningLargeCorrectionCoefficient(T, g, j, mix->strReaction[j], EsuR, n, CoeffCorr, MaxCoeffCorr);
				if (CoeffCorr < 0.)		WarningSmallCorrectionCoefficient(T, g, j, mix->strReaction[j], EsuR, n, CoeffCorr);
			}
			correction_k2[j] = CoeffCorr;
		}
	}
}

void WarningLargeCorrectionCoefficient(const double T, const double g, const int iReaction, const string stringReaction, const double EsuR, const double n, double &CoeffCorr, const double MaxCoeffCorr)
{
	//double csi = (T-Tmin)/DeltaT[k];
	//double Tk = KineticTemperature(CoeffCorr, T, EsuR, n);

	if (CoeffCorr >= 1.e48)
	{
		std::cout << "Correction coefficient too large..." << std::endl;
	/*	fWarning << "  Cluster:      " << k << endl;
		fWarning << "  Temperature:  " << T << " K" << endl;
		fWarning << "  g/gmax:       " << g/(1.-csi)/csi << endl;
		fWarning << "  Reaction:     " << iReaction	<<  "\t" << stringReaction << endl;
		fWarning << "  E:            " << EsuR*1.987 << " cal/mol"	<< endl;
		fWarning << "  n:            " << n			<< endl;
		fWarning << "  Correction:   " << CoeffCorr	<< endl;
		fWarning << "  Tk:           " << Tk << " K"<< endl;
		fWarning << "  Tk/T:         " << Tk/T << endl;*/
	}
	CoeffCorr = MaxCoeffCorr;
/*
	if (Tk > Tmax[k])
	{
		fWarning << "Kinetic equivalent temperature too large..." << endl;
		fWarning << "  Cluster:      " << k << endl;
		fWarning << "  Temperature:  " << T << " K" << endl;
		fWarning << "  g/gmax:       " << g/(1.-csi)/csi << endl;
		fWarning << "  Reaction:     " << iReaction	<<  "\t" << stringReaction << endl;
		fWarning << "  E:            " << EsuR*1.987 << " cal/mol"	<< endl;
		fWarning << "  n:            " << n			<< endl;
		fWarning << "  Correction:   " << CoeffCorr	<< endl;
		fWarning << "  Tk:           " << Tk << " K"<< endl;
		fWarning << "  Tmax:         " << Tmax[k] << " K"<< endl;
		fWarning << "  Tk/T:         " << Tk/T << endl;
	}*/
}

void WarningSmallCorrectionCoefficient(const double T, const double g, const int iReaction, const string stringReaction, const double EsuR, const double n, double &CoeffCorr)
{
	//double csi = (T-Tmin)/DeltaT[k];
	std::cout << "Negative correction coefficient..." << std::endl;
/*	fWarning << "  Cluster:      " << k << endl;
	fWarning << "  Temperature:  " << T << " K" << endl;
	fWarning << "  g/gmax:       " << g/(1.-csi)/csi << endl;
	fWarning << "  Reaction:     " << iReaction	<<  "\t" << stringReaction << endl;
	fWarning << "  E:            " << EsuR*1.987 << " cal/mol"	<< endl;
	fWarning << "  n:            " << n			<< endl;
	fWarning << "  Correction:   " << CoeffCorr	<< endl;
*/
	CoeffCorr = 1.;

//	ErrorMessage("Negative correction coefficient...\nSee Warning.log file...");
}

void OpenSMOKE_KPP_SingleReactor_KineticsManager::Setup(const KPP_Correction correction, OpenSMOKE_ReactingGas* mixture, const double temperature, const double pressure, const double variance)
{
	iSavedKineticConstants_ = true;
	mix_ = mixture;
	numberOfSpecies = mix_->NumberOfSpecies();
	cTot_ = pressure/Constants::R_J_kmol/temperature;
	M = mix_->M;
	uM = mix_->uM;
	mix_->ComputeKineticParameters(temperature, log(temperature), 1./temperature, pressure);


	// Corrections on reaction rates
	if (correction != KPP_CORRECTION_NONE)
	{
		BzzVector correction_k1(mix_->NumberOfReactions());
		BzzVector correction_k2(mix_->NumberOfReactions());
		BzzVector correction_uKeq(mix_->NumberOfReactions());

		double MaxCorrectionCoefficient = 1e10;

		CorrectionCoefficient( correction, mix_, temperature, cTot_, variance, Tmin_, Tmax_, MaxCorrectionCoefficient, correction_k1, correction_uKeq, correction_k2);
		mix_->CorrectKineticParameters(correction_k1, correction_uKeq, correction_k2);
	}

        if (mixture->NumberOfSpecies() ==  35 && mixture->NumberOfReactions() == 155)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX();
        else if (mixture->NumberOfSpecies() ==  53 && mixture->NumberOfReactions() == 325)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_GRI30();
        else if (mixture->NumberOfSpecies() ==  109 && mixture->NumberOfReactions() == 1764)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101();
        else if (mixture->NumberOfSpecies() ==  32 && mixture->NumberOfReactions() == 178)
        {
			if (mixture->name_kinetic_scheme == "Polimi_H2CONOX_1101")
				symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101();  
			else if (mixture->name_kinetic_scheme == "Polimi_H2CONOX_NoThermal_1101")
				symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101();   
			else if (mixture->name_kinetic_scheme == "Polimi_H2CONOX_NoNNH_1101")
				symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101();
			else if (mixture->name_kinetic_scheme == "Polimi_H2CONOX_NoN2O_1101")
				symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101();   
		}
        else if (mixture->NumberOfSpecies() ==  47 && mixture->NumberOfReactions() == 224)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_SanDiego_AVIO();
        else if (mixture->NumberOfSpecies() ==  86 && mixture->NumberOfReactions() == 1427)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO();  
 	else if (mixture->NumberOfSpecies() ==  81 && mixture->NumberOfReactions() == 1392)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702();  
        else
            ErrorMessage("No available kinetic scheme");
        
	symbolicKinetics->assignKineticConstants(mix_->k1, mix_->uKeq, mix_->logFcent, mix_->k2);

	ChangeDimensions(numberOfSpecies, &x_);
	ChangeDimensions(numberOfSpecies, &c_);
	ChangeDimensions(numberOfSpecies, &R_);
}

void OpenSMOKE_KPP_SingleReactor_KineticsManager::Setup(OpenSMOKE_ReactingGas* mixture) 
{
	iSavedKineticConstants_ = false;
	mix_ = mixture;
	numberOfSpecies = mix_->NumberOfSpecies();
	M = mix_->M;
	uM = mix_->uM;

        if (mixture->NumberOfSpecies() ==  35 && mixture->NumberOfReactions() == 155)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX();
        else if (mixture->NumberOfSpecies() ==  53 && mixture->NumberOfReactions() == 325)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_GRI30();
        else if (mixture->NumberOfSpecies() ==  109 && mixture->NumberOfReactions() == 1764)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_1101();
        else if (mixture->NumberOfSpecies() ==  32 && mixture->NumberOfReactions() == 178)
		{
			if (mixture->name_kinetic_scheme == "Polimi_H2CONOX_1101")
				symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_1101();  
			else if (mixture->name_kinetic_scheme == "Polimi_H2CONOX_NoThermal_1101")
				symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoThermal_1101();   
			else if (mixture->name_kinetic_scheme == "Polimi_H2CONOX_NoNNH_1101")
				symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoNNH_1101();
			else if (mixture->name_kinetic_scheme == "Polimi_H2CONOX_NoN2O_1101")
				symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_H2CONOX_NoN2O_1101();   
		}
        else if (mixture->NumberOfSpecies() ==  47 && mixture->NumberOfReactions() == 224)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_SanDiego_AVIO();
        else if (mixture->NumberOfSpecies() ==  86 && mixture->NumberOfReactions() == 1427)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_NC7_AVIO();  
 	else if (mixture->NumberOfSpecies() ==  81 && mixture->NumberOfReactions() == 1392)
                symbolicKinetics = new OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO_0702();  
        else
            ErrorMessage("No available kinetic scheme");
                
	ChangeDimensions(numberOfSpecies, &x_);
	ChangeDimensions(numberOfSpecies, &c_);
	ChangeDimensions(numberOfSpecies, &R_);
}

void OpenSMOKE_KPP_SingleReactor_KineticsManager::GetMoleFractionsFromMassFractions(BzzVector &x, BzzVector &y)
{
	ElementByElementProduct(y, uM, &x);
	double MWmix = 1./x.GetSumElements();
	x *= MWmix;
}

void OpenSMOKE_KPP_SingleReactor_KineticsManager::GetMWAndMoleFractionsFromMassFractions(double &MWmix, BzzVector &x, BzzVector &y)
{
	ElementByElementProduct(y, uM, &x);
	MWmix = 1./x.GetSumElements();
	x *= MWmix;
} 

double OpenSMOKE_KPP_SingleReactor_KineticsManager::GetMWFromMassFractions(BzzVector &y)
{
	BzzVector x(y.Size());
	ElementByElementProduct(y, uM, &x);
	return 1./x.GetSumElements();
}

void OpenSMOKE_KPP_SingleReactor_KineticsManager::RecoverKineticParameters(const double temperature, const double pressure)
{
	cTot_ = pressure/Constants::R_J_kmol/temperature;
	mix_->ComputeKineticParameters(temperature, log(temperature), 1./temperature, pressure);
	symbolicKinetics->assignKineticConstants(mix_->k1, mix_->uKeq, mix_->logFcent, mix_->k2);
}

void OpenSMOKE_KPP_SingleReactor_KineticsManager::UpdateProperties(BzzVector &omega, const double temperature, const double pressure, BzzVector &R)
{
	// Recover kinetic constants
	if (iSavedKineticConstants_ == false)
		RecoverKineticParameters(temperature, pressure);

	// Mole fractions from mass fractions
	GetMoleFractionsFromMassFractions(x_, omega);			// [-]

    // Concentrations
	c_ = cTot_ * x_;										// [kmol/m3]	

	// Formation Rates
	symbolicKinetics->giveReactionRates(cTot_, c_, R);		// [kmol/m3/s]
	ElementByElementProduct(R, M, &R);						// [kg/m3/s]
}

void OpenSMOKE_KPP_SingleReactor_KineticsManager::GetFormationRatesDerivatives(BzzVector& omega, const double temperature, const double pressure, BzzMatrix &dR_over_domega)
{	
	// Recover kinetic parameters
	if (iSavedKineticConstants_ == false)
		RecoverKineticParameters(temperature, pressure);

	double mw;
	GetMWAndMoleFractionsFromMassFractions(mw, x_, omega);
	Product(cTot_, x_, &c_);
	
	symbolicKinetics->giveReactionRates(cTot_, c_, R_);

	dR_over_domega = 0.;
	symbolicKinetics->giveJacobian(c_, dR_over_domega);

	double wc = mw * cTot_;
	for(int j=1;j<=numberOfSpecies;j++)
		x_[j] = wc * uM[j];
	dR_over_domega.ColumnsElementsProduct(x_);

	for(int i=1;i<=numberOfSpecies;i++)
		for(int j=1;j<=numberOfSpecies;j++)
			dR_over_domega[i][j] *= M[i];
}


