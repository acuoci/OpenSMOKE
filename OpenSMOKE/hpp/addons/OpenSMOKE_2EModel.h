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

#ifndef OPENSMOKE_2EMODEL
#define OPENSMOKE_2EMODEL

#include "BzzMath.hpp"

class OpenSMOKE_ReactingGas;

class OpenSMOKE_2EModel
{
private:

	// Soot models
	void nucleation_rate();
	void growth_rate();
	void oxidation_rate();
	void coagulation_rate();

	// Gas mixture
	OpenSMOKE_ReactingGas *mix;

	// Models indices
	int iNucleationModel;
	int iGrowthModel;
	int iOxidationModel;
	int iCoagulationModel;

	// Nucleation
	double SNGas_C2H2;
	double SNGas_H2;

	// Growth
	double SGGas_C2H2;
	double SGGas_H2;

	// Oxidation
	double SO_O2;
	double SO_OH;
	double SOGas_O2;
	double SOGas_OH;
	double SOGas_CO;
	double SOGas_H;

	// Coagulation
	double Sc;

	double cCoagulation;
	double etaNeoh;

	// Gas mixxture main data
	double T;
	double P;
	double C_C2H2;
	double C_O2;
	double C_OH;
	double p_O2;

	// Molecular weights
	double PM_C2H2;
	double PM_H2;
	double PM_O2;
	double PM_OH;
	double PM_CO;
	double PM_H;

	// Species indices
	int iO2;
	int iOH;
	int iC2H2;
	int iCO;
	int iH;
	int iH2;

	// Gas species sources
	double S_C2H2;
	double S_H2;
	double S_O2;
	double S_OH;
	double S_CO;
	double S_H;

	// Soot properties
	double ASoot;
	double dp;
	double MSoot;
	double omegaSoot;
	double xSoot;
	double cSoot;
	double NCarbons;

	void MessageError(string message);

public:

	void setup(int _iNucleationModel, int _iGrowthModel, int _iOxidationModel, int _iCoagulationModel, double _etaNeoh, double _Ccoagulation);
	void assign_mixture(OpenSMOKE_ReactingGas &_mix);
	void update(double T_K, double P_atm, double rho, double DiffC, BzzVector &_x, double phiN, double phiM);

	void setupFromFile(string fileName);
	void formation_rates();
	void initial_values(double rho);

	void write_on_file(ofstream &fOutput, double phiN, double phiM);
	void GnuPlotInterface(ofstream &fOutput, int count);

	// Soot sources
	double s;
	double S;
	double SN, Sn;
	double SG;

	// Gas species sources
	BzzVector SGas;

	// Soot main variables
	double fv;
	double m0;
	double Diff;

	// Start values
	double phiNStart;
	double phiMStart;
};

#endif // OPENSMOKE_2EMODEL
