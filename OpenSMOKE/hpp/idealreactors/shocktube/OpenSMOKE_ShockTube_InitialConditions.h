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

#ifndef OPENSMOKE_SHOCKTUBE_INITIALCONDITIONS_H
#define OPENSMOKE_SHOCKTUBE_INITIALCONDITIONS_H

#include "OpenSMOKE.hpp"
#include "idealreactors/OpenSMOKE_0DReactor.h"
#include "idealreactors/shocktube/OpenSMOKE_ShockTube_Geometry.h"

class   OpenSMOKE_ShockTube;

class OpenSMOKE_ShockTube_InitialConditions
{
public:
	OpenSMOKE_ShockTube_InitialConditions();
	
	void SetName(const std::string name);
	void SetGasMixture(OpenSMOKE_ReactingGas *_mix);
	void SetPointer();
	
	void Set_IncidentShock_CaseA(const double _UShock, const double _T1, const double _p1, BzzVector &_omega1);
	void Set_IncidentShock_CaseB(const double _UShock, const double _T2, const double _p2, BzzVector &_omega2);
	void Set_ReflectedShock_CaseA(const double _T5, const double _p5, BzzVector &_omega5);
	void Set_ReflectedShock_CaseB(OpenSMOKE_GasStream &shockStream, const double _UShock, const double _Urs, const double _T1, const double _p1, BzzVector &_omega1);
	void Set_ReflectedShock_CaseC(OpenSMOKE_GasStream &shockStream, const double _UShock, const double _Urs, const double _T2, const double _p2, BzzVector &_omega2);

	void Solve_IncidentShock_CaseA(OpenSMOKE_GasStream &shockStream);
	void Solve_IncidentShock_CaseB(OpenSMOKE_GasStream &shockStream);
	void Solve_ReflectedShock_CaseA(OpenSMOKE_GasStream &shockStream);
	void Solve_ReflectedShock_CaseB(OpenSMOKE_GasStream &shockStream);
	void Solve_ReflectedShock_CaseC(OpenSMOKE_GasStream &shockStream);

	void	VideoSummary();
	double	LengthBoundaryLayerCorrection(const double d);


public:

	double u1;
	double p1;
	double T1;
	double rho1;
	double h1;
	double pm1;
	double Cp1;
	double gamma1;
	double M1;
	double group1;
//	double mu1;
	BzzVector omega1;
	BzzVector x1;

	double u2;
	double p2;
	double T2;
	double rho2;
	double h2;
	double pm2;
	double Cp2;
	double gamma2;
	double M2;
	double group2;
//	double mu2;
	BzzVector omega2;
	BzzVector x2;

	double Urs;
	double u5;
	double p5;
	double T5;
	double rho5;
	double h5;
	double pm5;
	double Cp5;
	double gamma5;
	double M5;
//	double mu5;
	BzzVector omega5;
	BzzVector x5;


private:

	OpenSMOKE_ReactingGas *mix;
	double UShock;
	double Beta;
	bool iIncidentShock;
	double ALFAMAX;

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

private:

	double AlfaFirstGuessIncidentShock_CaseA();
	double AlfaFirstGuessIncidentShock_CaseB();
	double AlfaFirstGuessReflectedShock_CaseB();
	double AlfaMaxIncidentShock_CaseA();
	double AlfaMaxIncidentShock_CaseB();
	double AlfaMaxReflectedShock_CaseB();

public:

	double IncidentShock_CaseA(double alfa);
	double IncidentShock_CaseB(double alfa);
	double ReflectedShock_CaseB(double alfa);
};

#endif // OPENSMOKE_SHOCKTUBE_INITIALCONDITIONS_H
