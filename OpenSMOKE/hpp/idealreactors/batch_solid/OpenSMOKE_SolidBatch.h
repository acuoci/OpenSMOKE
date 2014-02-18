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

#ifndef OPENSMOKE_SOLIDBATCH_H
#define OPENSMOKE_SOLIDBATCH_H


#include "OpenSMOKE.hpp"
#include "idealreactors/OpenSMOKE_0DReactor.h"

class   OpenSMOKE_SolidBatch;
class   OpenSMOKE_SolidMixture;
void    ODEPrintExternal_SolidBatch(BzzVectorDouble &omega, double t);

class MyOdeSystem_SolidBatch : public BzzOdeSystemDoubleObject
{
private:
	bool iEnergy;
	bool iConstantPressure;

public:
	OpenSMOKE_SolidBatch *ptSolidBatch;
	void assignSolidBatch(OpenSMOKE_SolidBatch *solid_batch, const bool _iEnergy, const bool _iConstantPressure);
	virtual void GetSystemFunctions(BzzVectorDouble &y,double t,BzzVectorDouble &dy);
	virtual void ObjectBzzPrint(void);

};

class OpenSMOKE_Dictionary_SolidBatch : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_SolidBatch();
};


class OpenSMOKE_SolidBatch : public OpenSMOKE_0DReactor
{

public:	// Virtual functions (compulsory)

	OpenSMOKE_SolidBatch();

	void Lock();
	void VideoFinalResult();
	void VideoGeneralInfo();
	void VideoSummary();
	void DefineFromFile(const string inputFile);
	void Solve();
	void AssignEnd(const string units, const double value);
	void EnergyAnalysis(OpenSMOKE_GasStream &outletStream);
	void ODEPrint(BzzVectorDouble &y, double eta);
	void SummaryOnFile(const string file_name);

	
private: // Virtual functions (compulsory)
	
	void LabelODEFile();
	void Initialize();
	void UpdatePropertiesThermogravimetryIsothermal();
	void UpdateHeatFlux(const double tau, const double csi);
	void UpdateExchangeArea(const double tau, const double csi);
	void UpdateTwoEquationModel(BzzVectorDouble &y, BzzVectorDouble &dy);


public: // Specific functions 
	
	void AssignConstantPressure();
	void AssignVolume(const string units, const double value);
	void AssignSolidDiameter(const string units, const double value);
	void AssignSolidTotalMass(const string units, const double value);
	void AssignSolidKineticScheme(OpenSMOKE_SolidMixture &_solid);
	void AssignSolidMassFractions(const vector<string> _names, const vector<double> _values);

	void SetConstantExchangeArea(const double value, const string units);
	void SetUserDefinedExchangeArea(const string fileName);
	void UnsetUserDefinedExchangeArea();


public:	// Pseudo-Private Functions

	void ODESystem_Isothermal_SolidBatch_Thermogravimetry(BzzVectorDouble &x, double t, BzzVectorDouble &f);


private: // Specific data

	bool assignedVolume;
	bool assignedSolidKineticScheme;
	bool iConstantPressure;

	double volume;
	double mass;
	double moles;
	double A;

	// Solid mixture
	int NC_solid;
    int NR_solid;

	double solid_diameter_initial;
	double solid_mass_initial;
	double solid_mass;
	double solid_density;
	double solid_volume;
	double solid_diameter;
	BzzVectorDouble solid_species_mass;
	BzzVectorDouble dsolid_species_mass;
	BzzVectorDouble solid_omega;
	BzzVectorDouble solid_Rsolid;
	BzzVectorDouble solid_Rgas;
	BzzVectorDouble solid_C;

	// Solid mixture
	OpenSMOKE_SolidMixture *solid;

	// Geometry
	OpenSMOKE_UD_Profile	ud_A_profile;			// UD A profile

	// ODE System
	MyOdeSystem_SolidBatch ODE_SolidBatch_Object;

	// Virtual
	void UpdateProperties_isothermal(int memoIndex);
	void UpdateProperties(int jacobianIndex, int indexT);
};

#endif // OPENSMOKE_SOLIDBATCH_H
