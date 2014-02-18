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
#include "idealreactors/shocktube/OpenSMOKE_ShockTube_InitialConditions.h"
#include "basic/OpenSMOKE_Conversions.h"

OpenSMOKE_ShockTube_InitialConditions *pt_ics;

void OpenSMOKE_ShockTube_InitialConditions::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ShockTube_InitialConditions"		<< endl;
    cout << "Object: " << name_object							<< endl;
    cout << "Error:  " << message								<< endl;
    cout << "Press a key to continue... "						<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_ShockTube_InitialConditions::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:   OpenSMOKE_ShockTube_InitialConditions"	<< endl;
    cout << "Object:  " << name_object							<< endl;
    cout << "Warning: " << message								<< endl;
    cout << "Press a key to continue... "						<< endl;
    getchar();
}

OpenSMOKE_ShockTube_InitialConditions::OpenSMOKE_ShockTube_InitialConditions()
{
	name_object	= "Default name";	// Object Name
	ALFAMAX		= 30.;				// Maximum temperature ratio
}

void OpenSMOKE_ShockTube_InitialConditions::SetName(const string name)
{
	name_object	= name;				// Object Name
}

void OpenSMOKE_ShockTube_InitialConditions::SetPointer() 
{
	pt_ics = this;
}

double OpenSMOKE_ShockTube_InitialConditions::AlfaFirstGuessIncidentShock_CaseA()
{
	double alfa	 =	( gamma1*(M1*M1)-0.50*(gamma1-1.) ) * ( 0.50*(gamma1-1.)*(M1*M1)+1. ) / 
					( BzzPow2(0.50*(gamma1+1.)*M1) );
	return alfa;
}

double OpenSMOKE_ShockTube_InitialConditions::AlfaFirstGuessIncidentShock_CaseB()
{
	double gamma = 1.40;
	double M	 = 1.;
	double alfa	 =	( gamma*(M*M)-0.50*(gamma-1.) ) * ( 0.50*(gamma-1.)*(M*M)+1. ) / 
					( BzzPow2(0.50*(gamma+1.)*M) );
	return alfa;
}

double OpenSMOKE_ShockTube_InitialConditions::AlfaFirstGuessReflectedShock_CaseB()
{
	double eta = rho2/rho1;
	double csi = eta*(M1*M1*(eta-1.)*gamma1+eta)/(M1*M1*(eta-1.)*(gamma1-1.)+eta);
	double T5overT1 = 1.+(M1*M1)*(gamma1-1.)*(csi-1.)*(eta-1.)/(eta*(csi-eta));
	double alfa = T5overT1*T1/T2;
	return alfa;
}

double OpenSMOKE_ShockTube_InitialConditions::AlfaMaxIncidentShock_CaseA()
{
	return (group1+1.)*(group1+1.)/4./group1;
}

double OpenSMOKE_ShockTube_InitialConditions::AlfaMaxIncidentShock_CaseB()
{
	double csi = rho2*u1*u1/p2;
	if (csi >= 4.)	return ALFAMAX;
	else			return (sqrt(csi)+2.)/sqrt(csi)/(4.-csi);
}

double OpenSMOKE_ShockTube_InitialConditions::AlfaMaxReflectedShock_CaseB()
{
	if (Urs == 0.)	
		return ALFAMAX;
	else
	{
		double u2p		= Urs + u1 - u2;
		double coeff	= rho2*u2p*u2p/p2;
		return (1.+coeff)*(1.+coeff)/4./coeff; 
	}
}

void OpenSMOKE_ShockTube_InitialConditions::SetGasMixture(OpenSMOKE_ReactingGas *_mix)
{
	mix = _mix;
}

void OpenSMOKE_ShockTube_InitialConditions::Set_IncidentShock_CaseA(const double _UShock, const double _T1, const double _p1, BzzVector &_omega1)
{
	iIncidentShock = true;
	
	UShock	= _UShock;
	T1		= _T1;
	p1		= _p1;
	omega1	= _omega1;
	u1		= UShock;
	
	pm1		= mix->GetMWFromMassFractions(omega1);		
	rho1	= p1*pm1/Constants::R_J_kmol/T1;
	h1		= mix->GetMixEnthalpy_Mass(T1, omega1);
			  mix->SpeciesCp(T1);
	Cp1		= mix->MixCp_FromMassFractions(omega1);
	gamma1	= Cp1*pm1/(Cp1*pm1-Constants::R_J_kmol);
	M1		= u1*sqrt(rho1/gamma1/p1);
	mix->GetMoleFractionsFromMassFractionsAndMW(x1, omega1, pm1);

	group1  = rho1*u1*u1/p1;

	omega2	= omega1;
	pm2		= pm1;
	x2		= x1;
}

void OpenSMOKE_ShockTube_InitialConditions::Set_IncidentShock_CaseB(const double _UShock, const double _T2, const double _p2, BzzVector &_omega2)
{
	iIncidentShock = true;
	
	UShock	= _UShock;
	T2		= _T2;
	p2		= _p2;
	omega2	= _omega2;
	u1		= UShock;
	
	pm2		= mix->GetMWFromMassFractions(omega2);		
	rho2	= p2*pm2/Constants::R_J_kmol/T2;
	h2		= mix->GetMixEnthalpy_Mass(T2, omega2);
			  mix->SpeciesCp(T2);
	Cp2		= mix->MixCp_FromMassFractions(omega2);
	gamma2	= Cp2*pm2/(Cp2*pm2-Constants::R_J_kmol);
	mix->GetMoleFractionsFromMassFractionsAndMW(x2, omega2, pm2);

	omega1	= omega2;
	x1		= x2;
	pm1		= pm2;
}

void OpenSMOKE_ShockTube_InitialConditions::Set_ReflectedShock_CaseA(const double _T5, const double _p5, BzzVector &_omega5)
{
	iIncidentShock = false;

	Urs		= 0.;
	u5		= 0.;
	M5		= 0.;

	T5		= _T5;
	p5		= _p5;
	omega5	= _omega5;

	pm5		= mix->GetMWFromMassFractions(omega5);		
	rho5	= p5*pm5/Constants::R_J_kmol/T5;
	h5		= mix->GetMixEnthalpy_Mass(T5, omega5);
			  mix->SpeciesCp(T5);
	Cp5		= mix->MixCp_FromMassFractions(omega5);
	gamma5	= Cp5*pm5/(Cp5*pm5-Constants::R_J_kmol);
	mix->GetMoleFractionsFromMassFractionsAndMW(x5, omega5, pm5);

//	mix->SpeciesViscosityFromFitting(T5); 
//	mu5		= mix->MixViscosity_FromMolarFractions(x5);
}

void OpenSMOKE_ShockTube_InitialConditions::Set_ReflectedShock_CaseB(OpenSMOKE_GasStream &shockStream, const double _UShock, const double _Urs, const double _T1, const double _p1, BzzVector &_omega1)
{
	Set_IncidentShock_CaseA(_UShock, _T1, _p1, _omega1);
	Solve_IncidentShock_CaseA(shockStream);

	iIncidentShock = false;

	Urs		= _Urs;
	omega5	= omega1;
	x5		= x1;
	pm5		= pm1;
}

void OpenSMOKE_ShockTube_InitialConditions::Set_ReflectedShock_CaseC(OpenSMOKE_GasStream &shockStream, const double _UShock, const double _Urs, const double _T2, const double _p2, BzzVector &_omega2)
{
	Set_IncidentShock_CaseB(_UShock, _T2, _p2, _omega2);
	Solve_IncidentShock_CaseB(shockStream);

	iIncidentShock = false;

	Urs		= _Urs;
	omega5	= omega1;
	x5		= x1;
	pm5		= pm1;
}

double NLS_IncidentShock_CaseA(double alfa)
{
	return pt_ics->IncidentShock_CaseA(alfa);
}

double NLS_IncidentShock_CaseB(double alfa)
{
	return pt_ics->IncidentShock_CaseB(alfa);
}

double NLS_ReflectedShock_CaseB(double alfa)
{
	return pt_ics->ReflectedShock_CaseB(alfa);
}

double OpenSMOKE_ShockTube_InitialConditions::IncidentShock_CaseA(double alfa)
{
	Beta = 0.50*( (1.+group1)+sqrt(BzzPow2(1.+group1) - 4.*alfa*group1) );

	T2	= T1*alfa;
	p2	= p1*Beta;
	h2	= mix->GetMixEnthalpy_Mass(T2, omega2);

	double f = h1 + u1*u1/2.*(1.-alfa*alfa/Beta/Beta) - h2;

	return f;
}

double OpenSMOKE_ShockTube_InitialConditions::IncidentShock_CaseB(double alfa)
{
	double csi = rho2*u1*u1/p2;

	Beta = 0.50 * ( (1.+csi*alfa) + sqrt( BzzPow2(1.+csi*alfa) -4.*csi*alfa*alfa) );

	T1	= T2/alfa;
	p1	= p2/Beta;
	h1	= mix->GetMixEnthalpy_Mass(T1, omega1);
	
	double f = h1 + u1*u1/2.*(1.-alfa*alfa/Beta/Beta) - h2;

	return f;
}

void OpenSMOKE_ShockTube_InitialConditions::Solve_IncidentShock_CaseA(OpenSMOKE_GasStream &shockStream)
{
	BzzFunctionRootRobust m(AlfaFirstGuessIncidentShock_CaseA(), NLS_IncidentShock_CaseA, 1., AlfaMaxIncidentShock_CaseA());
	m.OnlyOneRoot();
	m();
	
	double alfa		= m.GetTSolution();
	double residual	= m.GetYSolution();

	T2		= alfa*T1;
	p2		= Beta*p1;
	rho2	= p2*pm2/Constants::R_J_kmol/T2;
	u2		= rho1*u1/rho2;

	h2		= mix->GetMixEnthalpy_Mass(T2, omega2);
			  mix->SpeciesCp(T2);
	Cp2		= mix->MixCp_FromMassFractions(omega2);
	gamma2	= Cp2*pm2/(Cp2*pm2-Constants::R_J_kmol);
	M2		= u2*sqrt(rho2/gamma2/p2);
	
	group2  = rho2*u2*u2/p2;

//	mix->SpeciesViscosityFromFitting(T1); 
//	mu1		= mix->MixViscosity_FromMolarFractions(x1);

//	mix->SpeciesViscosityFromFitting(T2); 
//	mu2		= mix->MixViscosity_FromMolarFractions(x2);

	shockStream.AssignPressure(p2, "Pa");
	shockStream.AssignTemperature(T2, "K");
	shockStream.AssignMassFractions(omega2);
	shockStream.AssignMassFlowRate(1., "kg/s");
	shockStream.SetVelocity(u2, "m/s");
}

void OpenSMOKE_ShockTube_InitialConditions::Solve_IncidentShock_CaseB(OpenSMOKE_GasStream &shockStream)
{
	BzzFunctionRootRobust m(AlfaFirstGuessIncidentShock_CaseB(), NLS_IncidentShock_CaseB, 1., AlfaMaxIncidentShock_CaseB());
	m.OnlyOneRoot();
	m();
	
	double alfa		= m.GetTSolution();
	double residual	= m.GetYSolution();

	T1		= T2/alfa;
	p1		= p2/Beta;
	rho1	= p1*pm1/Constants::R_J_kmol/T1;

	h1		= mix->GetMixEnthalpy_Mass(T1, omega1);
			  mix->SpeciesCp(T1);
	Cp1		= mix->MixCp_FromMassFractions(omega1);
	gamma1	= Cp1*pm1/(Cp1*pm1-Constants::R_J_kmol);

	u2		= rho1*u1/rho2;
	M1		= u1*sqrt(rho1/gamma1/p1);
	M2		= u2*sqrt(rho2/gamma2/p2);
	
	group1  = rho1*u1*u1/p1;
	group2  = rho2*u2*u2/p2;

//	mix->SpeciesViscosityFromFitting(T1); 
//	mu1		= mix->MixViscosity_FromMolarFractions(x1);

//	mix->SpeciesViscosityFromFitting(T2); 
//	mu2		= mix->MixViscosity_FromMolarFractions(x2);

	shockStream.AssignPressure(p2, "Pa");
	shockStream.AssignTemperature(T2, "K");
	shockStream.AssignMassFractions(omega2);
	shockStream.AssignMassFlowRate(1., "kg/s");
	shockStream.SetVelocity(u2, "m/s");
}

double OpenSMOKE_ShockTube_InitialConditions::ReflectedShock_CaseB(double alfa)
{
	double f;
	double u2p;

	double eta = rho2*rho1;
	double csi = eta*(M1*M1*(eta-1.)*gamma1+eta)/(M1*M1*(eta-1.)*(gamma1-1.)+eta);

	if (Urs == 0.)
	{
		double coeff	= 1.+rho2*(u1-u2)*(u1-u2)/p2 + alfa;
		Beta			= 0.50*(coeff+sqrt(coeff*coeff-4.*alfa));
		f				= h2+0.50*(u1-u2)*(u1-u2)*(1.+alfa/Beta)/(1.-alfa/Beta) - h5;
	}
	else
	{
		u2p		= Urs + u1 - u2;
		double coeff	= rho2*u2p*u2p/p2;
		Beta			= 0.50*( (1.+coeff) + sqrt(BzzPow2(1.+coeff) -4.*alfa*coeff));
	}

	T5	= alfa*T2;
	p5	= Beta*p2;
	h5	= mix->GetMixEnthalpy_Mass(T5, omega5);
	
	if (Urs == 0.)
		f				= h2+0.50*(u1-u2)*(u1-u2)*(1.+alfa/Beta)/(1.-alfa/Beta) - h5;
	else
		f				= h2+0.50*u2p*u2p*(1.-alfa*alfa/Beta/Beta)-h5;

	return f;
}

void OpenSMOKE_ShockTube_InitialConditions::Solve_ReflectedShock_CaseA(OpenSMOKE_GasStream &shockStream)
{
	shockStream.AssignPressure(p5, "Pa");
	shockStream.AssignTemperature(T5, "K");
	shockStream.AssignMassFractions(omega5);
	shockStream.AssignMassFlowRate(1., "kg/s");
	shockStream.SetVelocity(u5, "m/s");
}

void OpenSMOKE_ShockTube_InitialConditions::Solve_ReflectedShock_CaseB(OpenSMOKE_GasStream &shockStream)
{
	BzzFunctionRootRobust m(AlfaFirstGuessReflectedShock_CaseB(), NLS_ReflectedShock_CaseB, 1., AlfaMaxReflectedShock_CaseB());
	m.OnlyOneRoot();
	m();
	
	double alfa		= m.GetTSolution();
	double residual	= m.GetYSolution();

	T5		= alfa*T2;
	p5		= Beta*p2;
	rho5	= p5*pm5/Constants::R_J_kmol/T5;

	h5		= mix->GetMixEnthalpy_Mass(T5, omega5);
			  mix->SpeciesCp(T5);
	Cp5		= mix->MixCp_FromMassFractions(omega5);
	gamma5	= Cp5*pm5/(Cp5*pm5-Constants::R_J_kmol);

	if (Urs == 0.)	Urs = alfa/Beta*(u1-u2)/(1.-alfa/Beta);

	u5		= rho2/rho5*(Urs+u1-u2);
	M5		= u5*sqrt(rho5/gamma5/p5);

//	mix->SpeciesViscosityFromFitting(T5); 
//	mu5		= mix->MixViscosity_FromMolarFractions(x5);

	shockStream.AssignPressure(p5, "Pa");
	shockStream.AssignTemperature(T5, "K");
	shockStream.AssignMassFractions(omega5);
	shockStream.AssignMassFlowRate(1., "kg/s");
	shockStream.SetVelocity(u5, "m/s");
}

void OpenSMOKE_ShockTube_InitialConditions::Solve_ReflectedShock_CaseC(OpenSMOKE_GasStream &shockStream)
{
	Solve_ReflectedShock_CaseB(shockStream);
}

void OpenSMOKE_ShockTube_InitialConditions::VideoSummary()
{
	if (iIncidentShock == true)
	{
		cout << " -------------------------------------------------------------------" << endl;
		cout << "             Before\t\tAfter                                        " << endl;
		cout << " -------------------------------------------------------------------" << endl;
		cout << " T[K]        " << T1			<< "\t" << T2			<< endl;
		cout << " P[atm]      " << p1/101325.	<< "\t" << p2/101325.	<< endl;
		cout << " rho[kg/m3]  " << rho1			<< "\t" << rho2			<< endl;
		cout << " u[m/s]      " << u1			<< "\t" << u2			<< endl;
		cout << " M[-]        " << M1			<< "\t" << M2			<< endl;
		cout << " h[J/kg]     " << h1			<< "\t" << h2			<< endl;
		cout << " Cp[J/kg/K]  " << Cp1			<< "\t" << Cp2			<< endl;
		cout << " gamma[-]    " << gamma1		<< "\t" << gamma2		<< endl;
//		cout << " mu[Pa.s]    " << mu1			<< "\t" << mu2			<< endl;
	}
	else
	{
		cout << " -----------------------------------------------------------------------"		<< endl;
		cout << "             Before\t\tAfter\t\tReflected                               "		<< endl;
		cout << " -----------------------------------------------------------------------"		<< endl;
		cout << " T[K]        " << T1			<< "\t" << T2			<< "\t" << T5			<< endl;
		cout << " P[atm]      " << p1/101325.	<< "\t" << p2/101325.	<< "\t" << p5/101325.	<< endl;
		cout << " rho[kg/m3]  " << rho1			<< "\t" << rho2			<< "\t" << rho5			<< endl;
		cout << " u[m/s]      " << u1			<< "\t" << u2			<< "\t" << u5			<< endl;
		cout << " M[-]        " << M1			<< "\t" << M2			<< "\t" << M5			<< endl;
		cout << " h[J/kg]     " << h1			<< "\t" << h2			<< "\t" << h5			<< endl;
		cout << " Cp[J/kg/K]  " << Cp1			<< "\t" << Cp2			<< "\t" << Cp5			<< endl;
		cout << " gamma[-]    " << gamma1		<< "\t" << gamma2		<< "\t" << gamma5		<< endl;
//		cout << " mu[Pa.s]    " << mu1			<< "\t" << mu2			<< "\t" << mu5			<< endl;
	}
}

double OpenSMOKE_ShockTube_InitialConditions::LengthBoundaryLayerCorrection(const double d)
{
	ErrorMessage("LengthBoundaryLayerCorrection: current version requires transport properties...");
/*	double rhoWall = rho1*p2/p1;
	double muWall  = mu1;

	double C = pow( (rho2/rhoWall)*(mu2/muWall), 0.37);
	double Z = (gamma1+1.)/(gamma1-1.);
	double W = rho2/rho1;
	if (W>Z) Z=W;
	
	double B = 1.59*C*(1.+(1.796+0.802*W)/(Z*W-1.));
		
	return d*d/(16.*B*B)*BzzPow2(rho2/rhoWall)/(W-1.)*(u2/(muWall/rhoWall));*/

	return 0.;
}