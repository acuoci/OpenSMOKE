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

#ifndef OPENSMOKE_ICEM_H
#define OPENSMOKE_ICEM_H

#include "OpenSMOKE.hpp"
#include "idealreactors/OpenSMOKE_0DReactor.h"

class  OpenSMOKE_ICEM;

enum heat_transfer_models {ICEM_NONE, ICEM_BASE, ICEM_WOSCHINI};

class MyOdeSystem_ICEM : public BzzOdeSystemObject
{
private:

public:
	OpenSMOKE_ICEM *ptICEM;
	void assignICEM(OpenSMOKE_ICEM *icem);
	virtual void MyODEPrint(BzzVector &y, double eta);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Dictionary_ICEM : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_ICEM();
};


class OpenSMOKE_ICEM : public OpenSMOKE_0DReactor
{

public:	// Virtual functions (compulsory)

	OpenSMOKE_ICEM();

	void Lock();
	void VideoFinalResult();
	void VideoGeneralInfo();
	void VideoSummary();
	void DefineFromFile(const string inputFile);
	void Solve();
	void AssignEnd(const string units, const double value);
	void EnergyAnalysis(OpenSMOKE_GasStream &outletStream);
	void ODEPrint(BzzVector &y, double eta);
	void SummaryOnFile();

private: // Virtual functions (compulsory)
	
	void LabelODEFile();
	void Initialize();
	void UpdateProperties_isothermal(int memoIndex);
	void UpdateProperties(int jacobianIndex, int indexT);
	
	void UpdateHeatFlux(const double tau, const double csi);
	void UpdateExchangeArea(const double tau, const double csi);
	void UpdateTwoEquationModel(BzzVector &y, BzzVector &dy);


public: // Specific functions 
	
	void AssignClearanceVolume(const string units, const double value);
	void AssignRotationRate(const string units, const double value);
	void AssignCompressionRatio(const double value);
	void AssignArmRatio(const double value);
	void AssignStartAngle(const string units, const double value);
	void AssignDiameter(const string units, const double value);

	void SetConstantExchangeArea(const double value, const string units);
	void SetUserDefinedExchangeArea(const string fileName);
	void UnsetUserDefinedExchangeArea();

	void SetNumberOfCycles(const double value);
	void SetAdiabatic();


public:	// Pseudo-Private Functions

	void ODESystem_ICEM(BzzVector &x, double t, BzzVector &f);


private: // Specific data

	double mass;
	double moles;
	double A;

	double teta;
	double volume;

	double teta_start;
	double rotation_rate;
	double volume_clearance;
	double arm_ratio;
	double compression_ratio;
	double start_angle;
	double start_volume;
	double diameter;
	double area_base;
	double area_lateral;
	double number_of_cycles;
	double height;
	double dvolume_over_dt;
	double Cv;
	double UReaction;
	BzzVector u;

	// Geometry
	OpenSMOKE_UD_Profile	ud_A_profile;			// UD A profile

	// ODE System
	MyOdeSystem_ICEM ODE_ICEM_Object;

private:

	bool assignedClearanceVolume;
	bool assignedCompressionRatio;
	bool assignedArmRatio;
	bool assignedStartAngle;
	bool assignedRotationRate;
	bool assignedDiameter;

	void UpdateAngleAndVolume(const double t);

	static const double rad_to_degrees;
	static const double degrees_to_rad;

private:

	void UpdateHeatTransfer();

	heat_transfer_models heat_transfer_model;

	double Re;
	double Pr;
	double Nu;

	double a_Nu;
	double b_Nu;
	double c_Nu;

	double velocity_piston;
	double arm_La;
	double arm_Lc;

	double lambda;
	double mu;

	BzzVector lambdaMap;
	BzzVector muMap;
};

#endif // OPENSMOKE_ICEM_H
/*
		double area_effective;
		if (icem_exchange_area_model = ICEM_EXCHANGE_AREA_LATERAL)
			area_effective = area_cylinder_lateral;
		else if (icem_exchange_area_model = ICEM_EXCHANGE_AREA_TOTAL)
			area_effective = area_cylinder_lateral + 2.*area_cylinder_base;

		if (zone_count > 0)
		{
			for(int j=1;j<=N;j++)
				zone_area_effective[j] = zone_list_external_exchange[j]*area_effective;
		}
		else
			zone_area_effective[jOutmost] = area_effective;*/