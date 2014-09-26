/***************************************************************************
 *   Copyright (C) 2006-2009 by Alberto Cuoci								   *
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

#ifndef OPENSMOKE_CSTR_H
#define OPENSMOKE_CSTR_H

#include "OpenSMOKE.hpp"
#include "idealreactors/OpenSMOKE_0DReactor.h"
#include "idealreactors/cstr/OpenSMOKE_CSTR_Geometry.h"

class   OpenSMOKE_CSTR;
class   OpenSMOKE_Fluctuations;
void    ODEPrintExternal_CSTR(BzzVector &omega, double t);

class MyOdeSystem_CSTR : public BzzOdeSystemObject
{
private:
	bool iEnergy;
	bool iSetGasDischargeLaw;

public:
	OpenSMOKE_CSTR *ptCSTR;
	void assignCSTR(OpenSMOKE_CSTR *_cstr, const bool _iIsothermal, const bool _iSetGasDischargeLaw);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Dictionary_CSTR : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_CSTR();
};

class OpenSMOKE_Fluctuations
{
private:

	double mean; 
	double semiAmplitude;
	double alfa;
	double frequency;
	double tWall;
	double period;

	double *pointer;
	int iKind;

public:
	
	OpenSMOKE_Fluctuations();

	void 	Setup(const double _mean, const double _frequency, const double _semiAmplitude, const double _tWall);
	void	SetPointer(double &_pointer);
	void 	Update(const double _time);	
};

class OpenSMOKE_CSTR : public OpenSMOKE_0DReactor
{

public: // Virtual functions (compulsory)

	OpenSMOKE_CSTR();
	
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
	void UpdateHeatFlux(const double tau, const double csi);
	void UpdateExchangeArea(const double tau, const double csi);
	void UpdateProperties_isothermal(int memoIndex);
	void UpdateProperties(int jacobianIndex, int indexT);
	void UpdateTwoEquationModel(BzzVector &y, BzzVector &dy);

	void SaveOnBinaryFile(const std::string filename);
	void SaveOnBinaryFile(BzzSave &fOutput);
	void UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr);

public:	// Specific Functions
	
	void AssignResidenceTime(const std::string units, const double value);
	void AssignVolume(const std::string units, const double value);
	void AssignDiameter(const std::string units, const double value);

	void SetUserDefinedVolume(const std::string fileName);
	void SetUserDefinedDiameter(const std::string fileName);
	void SetUserDefinedResidenceTime(const std::string fileName);

	void SetConstantExchangeArea(const double value, const std::string units);
	void SetUserDefinedExchangeArea(const std::string fileName);
	void UnsetUserDefinedExchangeArea();

	void SetInitialTemperature(const double value, const std::string units);
	void SetInitialMassFractions(const vector<string> _names, const vector<double> _values);
	void SetInitialMoleFractions(const vector<string> _names, const vector<double> _values);

	void SetFluctuations(const vector<string> string_vector);

	void SetValveCoefficient(const double value, const std::string units);
	void SetPressureOutlet(const double value, const std::string units);

	void SetHistoryPartial();

	void Update(const double t);

public:	// Pseudo-Private Functions
	
	void ODESystem_Isothermal_CSTR(BzzVector &x, double t, BzzVector &f);
	void ODESystem_NonIsothermal_CSTR(BzzVector &x, double t, BzzVector &f);
	void ODESystem_Isothermal_GasDischargeLaw_CSTR(BzzVector &y, double t, BzzVector &dy);
	void ODESystem_NonIsothermal_GasDischargeLaw_CSTR(BzzVector &y, double t, BzzVector &dy);

	int  status;

private: // Specific data
	
	bool assignedGeometry;				//
	bool iResidenceTime;				//

	bool iSetInitialTemperature;
	bool iSetInitialComposition;
	bool iSetFluctuations;
	bool iSetGasDischargeLaw;
	bool iHistoryPartial;

	double Tinitial;
	BzzVector omega_initial;

	double Tau;								// Reactor residence time [s]
	double Volume;							// Reactor Volume [m3]
	double Area;							// Reactor Surface [m2]
	double D;								// Reactor internal diameter [m]

	double massFlowRate;					// Gas mass flow rate [kg/s]
	double Ae;								// Exchange surface [m2]

	double specificEnthalpy;				// Specific enthalpy [J/kg]
	double specificInternalEnergy;			// Specific internal energy [J/kg]
	BzzVector h;							// Specific enthalpies [J/kg]
	BzzVector hMap;							// Specific enthalpies [J/kg]

	double tEnd;							// Integration time [s]
	double mass_flow_out;					// Outlet mass flow rate [kg/s]
	double mass_total_initial;				// Total mass (initial) [kg]
	double mass_total;						// Total mass [kg]
	double PressureOutlet;					// Outlet pressure [Pa]
	double KValve;							// Valve flow coefficient [TODO]

	// ODE System
	MyOdeSystem_CSTR ODE_CSTR_Object;		//

	// Geometry
	OpenSMOKE_CSTR_Geometry		geometry;		// CSTR geometry
	OpenSMOKE_UD_Profile		ud_Ae_profile;	// UD Ae profile


private: // Fluctuations

	OpenSMOKE_Fluctuations fluctuations;
	void ExperimentOnFile();
};

#endif // OpenSMOKE_CSTR_H
