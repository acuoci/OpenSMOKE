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

#include <iomanip>
#include "engine/OpenSMOKE_ReactingGas.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "addons/OpenSMOKE_2EModel.h"

void OpenSMOKE_2EModel::setupFromFile(std::string fileName)
{
	const int SIZE = 200;
	char comments[SIZE];
	int _iNucleation, _iGrowth, _iOxidation, _iCoagulation;
	double _etaNeoh, _cCoagulation;

	ifstream fInput;
	openInputFileAndControl(fInput, fileName.c_str());

	fInput.getline(comments, SIZE);
	fInput.getline(comments, SIZE);
	fInput.getline(comments, SIZE);

	fInput >> _iNucleation;		fInput.getline(comments, SIZE);
	fInput >> _iGrowth;			fInput.getline(comments, SIZE);
	fInput >> _iOxidation;		fInput.getline(comments, SIZE);
	fInput >> _iCoagulation;	fInput.getline(comments, SIZE);

	fInput >> _etaNeoh;			fInput.getline(comments, SIZE);
	fInput >> _cCoagulation;	fInput.getline(comments, SIZE);

	fInput.close();

	setup(_iNucleation, _iGrowth, _iOxidation, _iCoagulation, _etaNeoh, _cCoagulation);

}

void OpenSMOKE_2EModel::setup(int _iNucleation, int _iGrowth, int _iOxidation, int _iCoagulation, double _etaNeoh, double _cCoagulation)
{
	if (_iNucleation>=0 && _iNucleation<=7)
		iNucleationModel	= _iNucleation;
	else MessageError("Wrong nucleation model!");

	if (_iGrowth>=0 && _iGrowth<=7)
		iGrowthModel		= _iGrowth;
	else MessageError("Wrong growth model!");
	
	if (_iOxidation>=0 && _iOxidation<=3)
		iOxidationModel		= _iOxidation;
	else MessageError("Wrong oxidation model!");

	if (_iCoagulation>=0 && _iCoagulation<=1)
		iCoagulationModel	= _iCoagulation;
	else MessageError("Wrong coagulation model!");

	if (_etaNeoh>=0 && _etaNeoh<=1)
		etaNeoh	= _etaNeoh;
	else MessageError("Wrong constant in Neoh model!");

	if ( _cCoagulation>=0)
		cCoagulation = _cCoagulation;
	else MessageError("Wrong Coagulation constant!");
}

void OpenSMOKE_2EModel::initial_values(double rho)
{
	switch(iNucleationModel)
	{
		case 1:		// Liu 2001

			dp = 0.27029949e-9*pow(90000., 1./3.);		// [m]

			break;
		
		case 2:		// Liu 2002,2006

			dp = 0.27029949e-9*pow(700., 1./3.);		// [m]

			break;

		case 3:		// Liu 2003

			dp = 0.27029949e-9*pow(700., 1./3.);		// [m]

			break;

		case 4:		// Moss 1999

			dp = 0.27029949e-9*pow(12., 1./3.);		// [m]

			break;

		case 5:		// Wen 2003

			dp = 0.27029949e-9*pow(12., 1./3.);		// [m]

			break;

		case 6:		// Lindstedt 1994

			dp = 0.27029949e-9*pow(60., 1./3.);		// [m]

			break;

		case 7:		// Leung 1991

			dp = 0.27029949e-9*pow(32., 1./3.);		// [m]

			break;

		case 0:

			dp = 0.;
			
			break;
 
	}

	double fvStart = 1.e-14;
	double m0Start = (6./Constants::pi/pow(dp,3.))*fvStart;

	
	phiNStart = m0Start/(rho*Constants::Nav_kmol);
	phiMStart = Constants::densitySoot*fvStart/rho;
}

void OpenSMOKE_2EModel::assign_mixture(OpenSMOKE_ReactingGas &_mix)
{
	mix = &_mix;

	iC2H2	= mix->recognize_species("C2H2");
	iO2		= mix->recognize_species("O2");
	iOH		= mix->recognize_species("OH");
	iH2		= mix->recognize_species("H2");
	iH		= mix->recognize_species("H");
	iCO		= mix->recognize_species("CO");

	PM_C2H2 = mix->M[iC2H2];
	PM_O2	= mix->M[iO2];
	PM_OH	= mix->M[iOH];
	PM_CO	= mix->M[iCO];
	PM_H2	= mix->M[iH2];
	PM_H	= mix->M[iH];

	ChangeDimensions(mix->NumberOfSpecies(), &SGas);
}


void OpenSMOKE_2EModel::update(double T_K, double P_atm, double rho, double DiffC, BzzVector &_x, double phiN, double phiM)
{
	T		= T_K;													// temperature [K]
	P		= P_atm;												// pressure [atm]

	double cTot	= P*101325.0/(8314.*T);								// total concentration [kmol/m3]

	double PMtot = rho/cTot;										// total molecular weight [kg/kmol]

	m0 = rho*Constants::Nav_kmol*phiN;								// soot particle number density [#/m3]

	fv = rho/Constants::densitySoot*phiM;							// volume fraction [-]
	
	ASoot = pow(36.*Constants::pi, 1./3.) * pow(m0, 1./3.) * pow(fv, 2./3.);	// soot specific area [1/m]

	dp	= pow(6.*fv/(Constants::pi*m0), 1./3.);						// soot particle diameter [m]

	MSoot	= Constants::densitySoot * fv;							// soot density [kg/m3]

	omegaSoot = MSoot / rho;										// soot mass fraction [-]
	
	xSoot	= omegaSoot * PMtot/Constants::Msoot;					// soot mole fraction [-]
	
	cSoot = xSoot * cTot;											// soot concentration [kmol/m3.s]
	
	NCarbons = BzzPow3(dp/0.27029949e-9);							// numero di carboni nella particella [-]
	
	Diff = DiffC / NCarbons;										// soot diffusivity [m2/s]

	
	// Concentrations

	p_O2	= _x[iO2]*P;						// oxygen partial pressure [atm]

	C_C2H2		= _x[iC2H2] * cTot;				// acetylene concentration [kmol/m3]
	
	C_O2		= _x[iO2] * cTot;				// oxygen concentration [kmol/m3]
	
	C_OH		= _x[iOH] * cTot;				// OH concentration [kmol/m3]

}


void OpenSMOKE_2EModel::nucleation_rate()
{
	double Mp;

	switch(iNucleationModel)
	{
		
		case 1:		// Liu 2001

			Sn = 30.0 * exp(-20643/T) * C_C2H2;			// [kmol/m3/s]

			Mp = Constants::Msoot*90000.;				// [kg/kmol]

			break;
		
		case 2:		// Liu 2002,2006

			Sn = 2.85 * exp(-16103/T) * C_C2H2;			// [kmol/m3/s]

			Mp = Constants::Msoot*700.;					// [kg/kmol]

			break;

		case 3:		// Liu 2003

			Sn = 0.004857 * exp(-7548/T) * C_C2H2;		// [kmol/m3/s]

			Mp = Constants::Msoot*700.;					// [kg/kmol]

			break;

		case 4:		// Moss 1999

			Sn = 54.0 * exp(-21100/T) * C_C2H2;			// [kmol/m3/s]

			Mp = Constants::Msoot*12.;					// [kg/kmol]

			break;

		case 5:		// Wen 2003

			Sn = 54.0 * exp(-21100/T) * C_C2H2;			// [kmol/m3/s]

			Mp = Constants::Msoot*12.;					// [kg/kmol]

			break;

		case 6:		// Lindstedt 1994

			Sn = 210.0 * exp(-21100/T) * C_C2H2;		// [kmol/m3/s]

			Mp = Constants::Msoot*60.;					// [kg/kmol]

			break;

		case 7:		// Leung 1991

			Sn = 635 * exp(-21100/T) * C_C2H2;			// [kmol/m3/s]

			Mp = Constants::Msoot*32.;					// [kg/kmol]

			break;

		case 0:

			Sn = 0.;
			
			Mp = 0.;
			
			break;
	}

	SN = Sn * Mp;										// [kg/m3/s]

	SNGas_C2H2	= - SN/(2.*Constants::Msoot) * PM_C2H2;	// [kg/m3/s] Acetylene consumption

	SNGas_H2	=   SN/(2.*Constants::Msoot) * PM_H2;	// [kg/m3/s] Hydrogen formation

}


void OpenSMOKE_2EModel::growth_rate()
{

	switch(iGrowthModel)
	{

		case 1:		// Liu 2001

			SG = 12000. * exp(-12083/T) * ASoot* C_C2H2;			// [kg/m3.s]

			break;
		
		case 2:		// Liu 2002,2006

			SG = 42000. * exp(-10064./T) * sqrt(ASoot)* C_C2H2;		// [kg/m3.s]					
		
			break;

		case 3:		// Liu 2003

			SG = 144. * exp(-6038./T) * ASoot* C_C2H2;				// [kg/m3.s]

			break;

		case 4:		// Moss 1999

			SG = 11700. * exp(-12100./T) * ASoot* C_C2H2;			// [kg/m3.s]		

			break;

		case 5:		// Wen 2003

			SG = 9000.6 * exp(-12100./T) * ASoot* C_C2H2;			// [kg/m3.s]	
		
			break;

		case 6:		// Lindstedt 1994

			SG = 18000. * exp(-12100./T) * ASoot* C_C2H2;			// [kg/m3.s]	

			break;

		case 7:		// Leung 1991

			SG = 144000. * exp(-12100./T) * sqrt(ASoot)* C_C2H2;	// [kg/m3.s]	

			break;

		case 0:

			SG = 0.;
			
			break;
	}

	SGGas_C2H2 = -SG / (2.*Constants::Msoot) * PM_C2H2;			// [kg/m3.s] Acetylene consumption

	SGGas_H2   =  SG / (2.*Constants::Msoot) * PM_H2;			// [kg/m3.s] Hydrogen formation

}

void OpenSMOKE_2EModel::oxidation_rate()
{
	double kA, kB, kT, kZ, chi;

	// WARNING: the reaction rates are expressed in [kg/m2/s] because they are multiplied by ASoot

	switch(iOxidationModel)
	{

		case 1:		// Lee 1962
		
			SO_O2 = 8903.51 * exp(-19778./T) * sqrt(T) * C_O2;		// [kg/m2.s]

			SO_OH = 0.;												// [kg/m3.s]

			break;
		
		case 2:		// Neoh 1980

			SO_OH = etaNeoh * 105.81 * sqrt(T) * C_OH;					// [kg/m2.s]

			SO_O2 = 0.;													// [kg/m3.s]
		
			break;

		case 3:		// NSC 1966

			kA = 20.	  * exp(-15098./T);									// [g/cm2.s.atm]
			
			kB = 4.46e-3  * exp(-7650./T);									// [g/cm2.s.atm]

			kT = 1.510e5  * exp(-48817./T);									// [g/cm2.s]

			kZ = 21.3	  * exp(2063./T);									// [1/atm]

			chi = 1./(1.+kT/(kB*p_O2));										// [-]

			SO_O2  = 12. * (kA*p_O2/(1.+kZ*p_O2)*chi + kB*p_O2*(1.-chi));	// [g/cm2.s]
			SO_O2 *= 10.;													// [kg/m2.s]

			SO_OH = 0.;														// [kg/m3.s]

			break;


		case 0:

			SO_O2 = 0.;			// [kg/m3.s]
			SO_OH = 0.;			// [kg/m3.s]
			
			break;

	}

	// Correction for fv --> 0
	// --------------------------------------------------------
	{
		const double fvstar = 1.e-9;
		const double ALFA  = 1.e-5;
		const double H = 1.50*log(ALFA/(1.-ALFA));
		const double K = 2.00*log((1.-ALFA)/ALFA) / fvstar;
		double delta = 1.e10;

		double m = ( tanh(K*fv + H) + 1. )/2.;
		double gammafv	=	m*pow(fv + m/delta, 2./3.) + 
							(1-m)*pow(fvstar, 2./3.-1.)*fv;

		double ASootCorrected = pow(36.*Constants::pi, 1./3.) * pow(m0, 1./3.) * gammafv;
		SO_O2 *= ASootCorrected;
		SO_OH *= ASootCorrected;
	}
 

	SOGas_O2 = - SO_O2/Constants::Msoot * 0.50 * PM_O2;						// [kg/m3.s]	O2 consumption
	
	SOGas_OH = - SO_OH/Constants::Msoot * 1.00 * PM_OH;						// [kg/m3.s]	OH consumption

	SOGas_CO =   (SO_O2+SO_OH)/Constants::Msoot * PM_CO;						// [kg/m3.s]	CO formation

	SOGas_H  =   SO_OH/Constants::Msoot * PM_H;								// [kg/m3.s]	H  formation

}

void OpenSMOKE_2EModel::coagulation_rate()
{
	double phi = 7.92e-40;

	switch(iCoagulationModel)
	{

		case 1:		// Smoluchowski
			
			{
				const double fvstar = 1.e-9;
				const double ALFA  = 1.e-5;
				const double H = 1.50*log(ALFA/(1.-ALFA));
				const double K = 2.00*log((1.-ALFA)/ALFA) / fvstar;
				double delta = 1.e10;

				double m = ( tanh(K*fv + H) + 1. )/2.;
				double gammafv	=	m*pow(fv + m/delta, 1./6.) + 
									(1-m)*pow(fvstar, 1./6.-1.)*fv;

				Sc = cCoagulation * phi * sqrt(T) * gammafv * pow(m0, 11./6.);		// [kmol/m3.s]
			}

			break;

		case 0:

			Sc = 0.;
			
			break;
	}


	// Correction for fv --> 0
	// --------------------------------------------------------
	{
	
	}

}

void OpenSMOKE_2EModel::formation_rates()
{

	// Soot formation/consumption rates
	// -----------------------------------------------------------------------
	nucleation_rate();
	growth_rate();
	oxidation_rate();
	coagulation_rate();

	// Soot source terms
	// -----------------------------------------------------------------------
	s = Sn - Sc;						// [kmol/m3/s]
	S = SN + SG - (SO_O2 + SO_OH);		// [kg/m3/s]


	// Gas species source terms
	// -----------------------------------------------------------------------

	S_C2H2	= SNGas_C2H2 + SGGas_C2H2;	// [kg/m3.s] Acetylene consumption

	S_H2	= SNGas_H2 + SGGas_H2;		// [kg/m3.s] Hydrogen formation

	S_OH	= SOGas_OH;					// [kg/m3.s] OH consumption

	S_O2	= SOGas_O2;					// [kg/m3.s] O2 consumption

	S_CO	= SOGas_CO;					// [kg/m3.s] CO formation

	S_H		= SOGas_H;					// [kg/m3.s] H formation

	SGas[iC2H2] = 	S_C2H2;				// [kg/m3.s] Acetylene consumption
	SGas[iH2]	= 	S_H2;				// [kg/m3.s] Hydrogen formation
	SGas[iOH]	= 	S_OH;				// [kg/m3.s] OH consumption
	SGas[iCO]	= 	S_CO;				// [kg/m3.s] CO formation
	SGas[iO2]	= 	S_O2;				// [kg/m3.s] O2 consumption
	SGas[iH]	= 	S_H;				// [kg/m3.s] H formation
}


void OpenSMOKE_2EModel::write_on_file(ofstream &fOutput, double phiN, double phiM)
{
	fOutput << setw(20) << left << phiN;		// phiN
	fOutput << setw(20) << left << phiM;		// phiM
	fOutput << setw(20) << left << m0;		// soot particle number density [#/m3]
	fOutput << setw(20) << left << fv;		// soot volume fraction [-]
	fOutput << setw(20) << left << omegaSoot;		// soot mass fraction [-]
	fOutput << setw(20) << left << xSoot;		// soot mole fraction [-]
	fOutput << setw(20) << left << cSoot;		// soot concentration [kmol/m3]
	fOutput << setw(20) << left << MSoot;		// soot mass density [kg/m3]
	fOutput << setw(20) << left << dp;		// soot particle diameter [m]
	fOutput << setw(20) << left << dp*1.e9;		// soot particle diameter [nm]
	fOutput << setw(20) << left << ASoot;		// soot specific area [1/m]
	fOutput << setw(20) << left << NCarbons;		// number of carbon atoms in the single particle [-]
	fOutput << setw(20) << left << Diff;		// soot diffusivity [m2/s]
	
	fOutput << setw(20) << left << s;		// soot particle number density source term [kmol/m3.s]
	fOutput << setw(20) << left << S;		// soot mass fraction source term [kg/m3.s]
	fOutput << setw(20) << left << Sn;		// soot particle number density source term: nucleation [kmol/m3.s]
	fOutput << setw(20) << left << Sc;		// soot particle number density source term: coagulation [kmol/m3.s]
	fOutput << setw(20) << left << SN;		// soot mass fraction source term: nucleation [kg/m3.s]
	fOutput << setw(20) << left << SG;		// soot mass fraction source term: growth [kg/m3.s]
	fOutput << setw(20) << left << SO_O2;		// soot mass fraction source term: O2 oxidation [kg/m3.s]
	fOutput << setw(20) << left << SO_OH;		// soot mass fraction source term: OH oxidation [kg/m3.s]
	
	fOutput << setw(20) << left << S_C2H2/PM_C2H2;	// gas source term: acetylene [kmol/m3.s]
	fOutput << setw(20) << left << S_H2/PM_H2;	// gas source term: hydrogen [kmol/m3.s]
	fOutput << setw(20) << left << S_OH/PM_OH;	// gas source term: OH [kmol/m3.s]
	fOutput << setw(20) << left << S_O2/PM_O2;	// gas source term: oxygen [kmol/m3.s]
	fOutput << setw(20) << left << S_CO/PM_CO;	// gas source term: carbon monoxide [kmol/m3.s]
	fOutput << setw(20) << left << S_H/PM_H	;	// gas source term: H [kmol/m3.s]

	fOutput << endl;
}

void OpenSMOKE_2EModel::GnuPlotInterface(ofstream &fOutput, int count)
{
	PrintTagOnGnuplotLabel(20, fOutput, "phiN",		count);
	PrintTagOnGnuplotLabel(20, fOutput, "phiM",		count);
	PrintTagOnGnuplotLabel(20, fOutput, "m0[#/m3]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "fv[-]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "wSoot[-]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "xSoot[-]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "cSoot[-]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "MSoot[-]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "dp[m]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "dp[nm]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "ASoot[1/m]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "NCarb[-]",		count);
	PrintTagOnGnuplotLabel(20, fOutput, "Diff[m2/s]",	count);

	PrintTagOnGnuplotLabel(20, fOutput, "s[kmol/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "S[kg/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "Sn[kmol/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "Sc[kmol/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "SN[kg/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "SG[kg/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "SO_O2[kg/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "SO_OH[kg/m3s]",	count);

	PrintTagOnGnuplotLabel(20, fOutput, "C2H2[kmol/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "H2[kmol/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "OH[kmol/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "O2[kmol/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "CO[kmol/m3s]",	count);
	PrintTagOnGnuplotLabel(20, fOutput, "H[kmol/m3s]",	count);
}

void OpenSMOKE_2EModel::MessageError(std::string message)
{
	std::cout << "Class: OpenSMOKE_2EModel"		<< std::endl;
	std::cout << "Error: " << message.c_str()	<< std::endl;
	std::cout << "Press enter to continue... "	<< std::endl;
	getchar();
	exit(-1);
}