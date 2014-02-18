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

#ifndef OPENSMOKE_PFR_H
#define OPENSMOKE_PFR_H

#include "OpenSMOKE.hpp"
#include "idealreactors/OpenSMOKE_0DReactor.h"

class   OpenSMOKE_PFR;

class MyOdeSystem_PFR : public BzzOdeSystemObject
{
private:
	bool iEnergy;
	bool iMomentum;

public:
	OpenSMOKE_PFR *ptPFR;
	void assignPFR(OpenSMOKE_PFR *pfr, bool iIsothermal, bool iMomentum);
	virtual void MyODEPrint(BzzVector &y, double eta);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Dictionary_PFR : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_PFR();
};

enum ode_solver_type { ODE_SOLVER_BZZ, ODE_SOLVER_RADAU, ODE_SOLVER_DLSODE };

class OpenSMOKE_PFR : public OpenSMOKE_0DReactor
{

public: // Virtual functions (compulsory)

	OpenSMOKE_PFR();
	
	void Lock();
	void VideoFinalResult();
	void VideoGeneralInfo();
	void VideoSummary();
	void DefineFromFile(const string inputFile);
	void Solve();
	void ReSolve();
	void AssignEnd(const string units, const double value);
	void EnergyAnalysis(OpenSMOKE_GasStream &outletStream);
	void ODEPrint(BzzVector &y, double eta);
	void SummaryOnFile();

	double GetContactTime() { return (TauTotal + initialTime); }
	double GetLenght()		{ return (L + initialLenght); }


private: // Virtual functions (compulsory)
	
	void LabelODEFile();
	void Initialize();
	void UpdateHeatFlux(const double tau, const double csi);
	void UpdateExchangeArea(const double tau, const double csi);
	void UpdateProperties_isothermal(int memoIndex);
	void UpdateProperties(int jacobianIndex, int indexT);
	void UpdateTwoEquationModel(BzzVector &y, BzzVector &dy);
	void UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr);

	void SaveOnBinaryFile(const string filename);
	void SaveOnBinaryFile(BzzSave &fOutput);

public:	// Specific Functions
	
	void AssignDiameter(const string cm_or_m, const double value);
	void AssignArea(const string units, const double value);

	void SetUserDefinedDiameter(const string fileName);
	void SetUserDefinedArea(const string fileName);

	void SetMomentumEquation();
	void UnsetMomentumEquation();

	void SetInitialTime(const string units, const double value);
	void SetInitialLenght(const string units, const double value);

	void SetConstantSpecificExchangeArea(const double value, const string units);
	void SetUserDefinedSpecificExchangeArea(const string fileName);
	void UnsetUserDefinedSpecificExchangeArea();

	void SetTimeIndependent() { iTimeIndipendent = true; }

public:	// Pseudo-Private Functions
	
	void ODESystem_Isothermal_PFR(BzzVector &x, double t, BzzVector &f);
	void ODESystem_IsothermalMomentum_PFR(BzzVector &y, double t, BzzVector &dy);
	void ODESystem_NonIsothermal_PFR(BzzVector &x, double t, BzzVector &f);
	void ODESystem_NonIsothermalMomentum_PFR(BzzVector &x, double t, BzzVector &f);

	int  status;

private: // Specific data

	bool assignedGeometry;					//
	
	bool iTimeIndipendent;					//
	bool iMomentum;							//

	double L;								// Reactor length [m]
	double Tau;								// Reactor residence time [s]
	double Csi;								// Reactor axial coordinate [m]
	double Area;							// Reactor cross section [m2]
	double D;								// Reactor internal diameter [m]
	double v;								// Gas velocity [m/s]
	double massFlowRate;					// Gas mass flow rate [kg/s]
	double Ae;								// Exchange surface per unit volume [1/m]

	double dmassFlowRate;					// 
	double dEta;							//
	double dP;								//

	// ODE System
	BzzOdeStiffObject	o;
	MyOdeSystem_PFR		ODE_PFR_Object;		//

	// Geometry
	OpenSMOKE_PFR_Geometry  geometry;		// PFR geometry
	OpenSMOKE_UD_Profile	ud_Ae_profile;	// UD Ae profile

	double initialTime;
	double initialLenght;


private: // Momentum equation data
	
	double delta;
	double FForce;
	double Beta;
	double T_coefficient;
	double T_extra_terms;

	double	FrictionFactor();
	double	FrictionForce();
	double	Delta();
	double	dv();
	void	AdditionalExtraTerms();
	void	UpdateVelocity();
	void	CheckPressure();

	int     nExperiments;
	bool	iLogExperiment;
	void	ExperimentOnFile();
	void	SetLogExperiment(const int int_value);

	BzzVector S_right;

	ode_solver_type iODE_Solver;
	void	SetODESolver(const string value);

private:

	ofstream fSootDistribution;

};

#endif // OPENSMOKE_PFR_H
