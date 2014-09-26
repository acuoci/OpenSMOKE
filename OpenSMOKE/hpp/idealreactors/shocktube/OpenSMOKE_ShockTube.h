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

#ifndef OPENSMOKE_SHOCKTUBE_H
#define OPENSMOKE_SHOCKTUBE_H

#include "OpenSMOKE.hpp"
#include "idealreactors/OpenSMOKE_0DReactor.h"
#include "idealreactors/shocktube/OpenSMOKE_ShockTube.h"
#include "idealreactors/shocktube/OpenSMOKE_ShockTube_Geometry.h"
#include "idealreactors/shocktube/OpenSMOKE_ShockTube_InitialConditions.h"

class   OpenSMOKE_ShockTube;
class	OpenSMOKE_ShockTube_InitialConditions;
void    ODEPrintExternal_ShockTube(BzzVector &omega, double t);

class MyOdeSystem_ShockTube : public BzzOdeSystemObject
{
private:
	bool iBoundaryLayerCorrection;

public:
	OpenSMOKE_ShockTube *ptShockTube;
	void assignShockTube(OpenSMOKE_ShockTube *shocktube, bool iBoundaryLayerCorrection);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Dictionary_ShockTube : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_ShockTube();
};

enum kindInitialCondition {SHOCK_NONE, SHOCK_BEFORE, SHOCK_AFTER, SHOCK_REFLECTED};


class OpenSMOKE_ShockTube : public OpenSMOKE_0DReactor
{

public: // Virtual functions (compulsory)

	OpenSMOKE_ShockTube();
	
	void Lock();
	void VideoFinalResult();
	void VideoGeneralInfo();
	void VideoSummary();
	void DefineFromFile(const std::string inputFile);
	void Solve();
	void ReSolve();
	void AssignEnd(const std::string units, const double value);
	void EnergyAnalysis(OpenSMOKE_GasStream &outletStream);
	void ODEPrint(BzzVector &y, double eta);
	void SummaryOnFile();

private: // Virtual functions (compulsory)
	
	void LabelODEFile();
	void Initialize();
	void UpdateProperties(int jacobianIndex, int indexT);
	void UpdateProperties_isothermal(int memoIndex) {};
	void UpdateHeatFlux(const double tau, const double csi) {};
	void UpdateExchangeArea(const double tau, const double csi) {};
	void UpdateTwoEquationModel(BzzVector &y, BzzVector &dy);
	double AreaBoundaryLayerCorrection();

	void SaveOnBinaryFile(const std::string filename);
	void SaveOnBinaryFile(BzzSave &fOutput);
	void UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr);


public:	// Specific Functions
	
	void AssignDiameter(const std::string cm_or_m, const double value);
	void AssignArea(const std::string units, const double value);
	void AssignIncidentShock();
	void AssignReflectedShock();
	void AssignShockVelocity(const std::string units, const double value);
	void AssignShockTemperature(const std::string units, const double value, const kindInitialCondition _iInitialCondition);
	void AssignShockPressure(const std::string units, const double value, const kindInitialCondition _iInitialCondition);
	void AssignShockDensity(const std::string units, const double value, const kindInitialCondition _iInitialCondition);
	void AssignShockMoleFractions(const vector<string> _names, const vector<double> _values);
	void AssignShockMassFractions(const vector<string> _names, const vector<double> _values);

	void AssignShockEquivalenceRatio(const double value);
	void AssignShockFuelMassFractions(const vector<string> names, const vector<double> values);
	void AssignShockFuelMoleFractions(const vector<string> names, const vector<double> values);
	void AssignShockOxidizerMassFractions(const vector<string> names, const vector<double> values);
	void AssignShockOxidizerMoleFractions(const vector<string> names, const vector<double> values);
	void AssignShockStreamComposition();
	
	void SetBoundaryLayerCorrection();
	void UnsetBoundaryLayerCorrection();
	void SetReflectedShockVelocity(const std::string units, const double value);

public:	// Pseudo-Private Functions
	
	void ODESystem_Isothermal_ShockTube(BzzVector &x, double t, BzzVector &f);
	int  status;

public:

	double	equivalence_ratio;
	bool	iShockEquivalenceRatio;
	vector<string> fuel_names;
	BzzVector moles_fuel;
	vector<string> oxidizer_names;
	BzzVector moles_oxidizer;

private: // Specific data

	bool assignedGeometry;				//
	bool assignedKindOfReactor;			//
	bool assignedShockVelocity;			//
	bool assignedShockTemperature;		//
	bool assignedShockPressure;			//
	bool assignedShockDensity;			//
	bool assignedShockComposition;		//
	
	bool iBoundaryLayerCorrection;		//
	bool iIncidentShock;

	kindInitialCondition iInitialCondition;

	double UShock;
	double UReflectedShock;
	double Tau;								// Reactor residence time [s]
	double Area;							// Reactor cross section [m2]
	double D;								// Reactor internal diameter [m]
	double lm;								// Length for viscosity correction [m]
	double Area0;							// Starting section [m2]

	double Csi;								// Reactor axial coordinate [m]
	double v;								// Gas velocity [m/s]
	double tL;								// Laboratory time [s]

	double drho;					// 
	double dtL;						//
	double dCsi;					//
	double dv;						//

	// ODE System
	MyOdeSystem_ShockTube ODE_ShockTube_Object;		//

	// Initial conditions
	OpenSMOKE_ShockTube_InitialConditions ics;

	// Streams
	OpenSMOKE_GasStream shockStream;

private:
	static const double EPSILON;
	static const double _1_PLUS_EPSILON;
	static const double _1_MINUS_EPSILON;
};

#endif // OPENSMOKE_SHOCKTUBE_H
