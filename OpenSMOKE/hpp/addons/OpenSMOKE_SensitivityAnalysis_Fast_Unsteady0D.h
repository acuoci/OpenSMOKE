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

#ifndef OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D
#define OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Constants.h"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_Kinetics;
class OpenSMOKE_NuManager;

enum kind_of_sensitivity_analysis {	OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_NONE, 
									OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ONLY_SPECIES, 
									OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_SPECIES_PLUS_TEMPERATURE,
									OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ICEM_MULTIZONE};

enum kind_of_integration	{	OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_IMPLICIT, 
								OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_EULER_EXPLICIT,
								OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D_ACCURATE
							};

class OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D
{

friend  class OpenSMOKE_Kinetics;

public:

	OpenSMOKE_SensitivityAnalysis_Fast_Unsteady0D();
	void SetName(const string name);
	void Initialize(const kind_of_sensitivity_analysis kind_, const kind_of_integration integration_, OpenSMOKE_ReactingGas *_mix, BzzVectorInt &_indices, const string path_name);
	void Update(const double tCurrent, BzzMatrix &J, const double T, 
				const double P_Pascal, const double rho, const double Cp, const double v, BzzVector &X);

	void Update(BzzMatrix &J, const double T, const double P_Pascal, 
				const double rho, const double Cp, const double v, BzzVector &X,
				BzzVector &x, BzzVector &r);

	void UpdateOutputFile(const double t, BzzVector &x);
	void MoveFromFiles(const string fileName, const int N, BzzSave &fOutput, const string tag);

	
	void VideoSummary();

	inline const BzzMatrix &S()					{ return _S; }
	inline const BzzMatrix &dS_over_dt()		{ return _dS_over_dt; }

	inline int NumberOfCoefficients()			{ return NV*NP; }
	inline kind_of_integration Integration()	{ return integration; }

	void PrepareBinaryFile(BzzSave &fOutput, const int N, BzzMatrix &S, const string tag);
	void SaveOnBinaryFile(BzzSave &fOutput);
	void CloseBinaryFile(BzzSave &fOutput);
	void Close();

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);

private:

	OpenSMOKE_ReactingGas	*mix;
	OpenSMOKE_NuManager		*nu;
	BzzSave					fBinary;

	int NC;		// Number of species
	int NR;		// Number of reactions
	int NP;		// Number of parameters
	int NV;		// Number of unknowns (per point)
	int NC_RED;	// Number of species (per point): reduced
	int NV_RED;	// Number of unknowns (per point): reduced

	int indexTemperature;		// index of temperature
	int indexSpecies;			// index of species
	int indexPressure;			// index of species
	int indexHeatRelease;		// index of species
	int additional_variables;	// additional variables

	double tOld;
	
	BzzVectorInt _indices_print_species;

	BzzMatrix _S;
	BzzMatrix _S_S;
	BzzMatrix _Sold;
	BzzMatrix _JAlfa;
	BzzVector _parameters;

	BzzFactorizedGauss	_A_Euler_Implicit_Factorized;
	BzzMatrix			_A_Euler_Implicit;
	BzzMatrix			_A_Euler_Explicit;
	BzzMatrix			_dS_over_dt;
	BzzMatrixDiagonal	_I;

	kind_of_sensitivity_analysis	reacting_system;
	kind_of_integration				integration;
	void BuildNuMatrix(OpenSMOKE_Kinetics *kinetics);

	double	threshold_normalization;

	ofstream fOutput;
	void PrepareOutputFile(const string file_name);
	void UpdateOutputFile(const double t);
	void Populate_S_S();

	void Euler_Explicit(const double deltat, BzzMatrix &J);
	void Euler_Implicit(const double deltat, BzzMatrix &J);

	void BuildJAlfaMatrix(const double T, const double P_Pascal, const double rho, const double denominator, const double v, BzzVector &X);

	// Old values
	BzzMatrix _JOld;
};

#endif // OPENSMOKE_SENSITIVITYANALYSIS_FAST_UNSTEADY0D

