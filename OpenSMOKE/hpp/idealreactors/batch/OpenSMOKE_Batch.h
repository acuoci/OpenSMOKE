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

#ifndef OPENSMOKE_Batch_H
#define OPENSMOKE_Batch_H

#include "OpenSMOKE.hpp"
#include "idealreactors/OpenSMOKE_0DReactor.h"

class   OpenSMOKE_Batch;
void    ODEPrintExternal_Batch(BzzVector &omega, double t);

class MyOdeSystem_Batch : public BzzOdeSystemObject
{
private:
	bool iEnergy;
	bool iConstantPressure;

public:
	OpenSMOKE_Batch *ptBatch;
	void assignBatch(OpenSMOKE_Batch *batch, const bool _iEnergy, const bool _iConstantPressure);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);

};

class OpenSMOKE_Dictionary_Batch : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_Batch();
};


class OpenSMOKE_Batch : public OpenSMOKE_0DReactor
{

public:	// Virtual functions (compulsory)

	OpenSMOKE_Batch();

	void Lock();
	void VideoFinalResult();
	void VideoGeneralInfo();
	void VideoSummary();
	void DefineFromFile(const std::string inputFile);
	void Solve();
	void AssignEnd(const std::string units, const double value);
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

	void SaveOnBinaryFile(const std::string filename);
	void SaveOnBinaryFile(BzzSave &fOutput);
	void UpdateReactionRates(const double TT, const double PP, BzzVector &xx, BzzVector &rr);

public: // Specific functions 
	
	void AssignConstantPressure();
	void AssignVolume(const std::string units, const double value);

	void SetConstantExchangeArea(const double value, const std::string units);
	void SetUserDefinedExchangeArea(const std::string fileName);
	void SetVolumeLaw(const double value, const std::string units);
	void UnsetUserDefinedExchangeArea();


public:	// Pseudo-Private Functions

	void ODESystem_Isothermal_Batch(BzzVector &x, double t, BzzVector &f);
	void ODESystem_NonIsothermal_ConstantPressure_Batch(BzzVector &x, double t, BzzVector &f);
	void ODESystem_NonIsothermal_ConstantVolume_Batch(BzzVector &x, double t, BzzVector &f);


private: // Specific data

	bool assignedVolume;
	bool iConstantPressure;

	double volume;
	double mass;
	double moles;
	double A;

	double specificEnthalpy;				// Specific enthalpy [J/kg]
	double specificInternalEnergy;			// Specific internal energy [J/kg]
	BzzVector h;							// Specific enthalpies [J/kg]
	BzzVector hMap;							// Specific enthalpies [J/kg]

	// Geometry
	OpenSMOKE_UD_Profile	ud_A_profile;			// UD A profile

	// ODE System
	MyOdeSystem_Batch ODE_Batch_Object;

	// Volume law
	double volumeInitial;
	double alfaVolumeLaw;
	bool iVolumeLaw;
};

#endif // OPENSMOKE_Batch_H
