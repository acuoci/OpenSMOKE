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

#ifndef OPENSMOKE_SENSITIVITYANALYSIS_FLAME1D_POSTPROCESSOR
#define OPENSMOKE_SENSITIVITYANALYSIS_FLAME1D_POSTPROCESSOR

#include "BzzMath.hpp"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_Kinetics;
class OpenSMOKE_NuManager;

class OpenSMOKE_SensitivityAnalysis_Flame1D_PostProcessor
{
public:
	OpenSMOKE_SensitivityAnalysis_Flame1D_PostProcessor();
	void SetName(const string name);

private:
	BzzVector x;
	BzzVector T;
	BzzMatrix omega;
	BzzVector Mtot;
	BzzVector M;
	BzzVector rho;
	BzzVector H;
	BzzVector U;
	BzzVector G;
	BzzVector parameters;
	BzzMatrix S;

	BzzVector vectorMax;
	BzzVector vectorMin;
	BzzMatrix Slocal;
	BzzVector SlocalVector;
	BzzMatrix matrixVariables;

	int	index_temperature;
	int	index_massflowrate;
	int	index_species;
	int	index_velocity;
	int	index_pressurecurvature;

	int	nExtracted;

	void GetSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, const int index_local, BzzMatrix &SMatrix, BzzVectorInt &indices);
	void GetSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, const int index_local, BzzVector &SVector, BzzVectorInt &indices);

public:
	
	int NC;
	int NP;
	int NV;
	int N;

	int		 S_NC;
	int		 S_NV;
	string	*s_names;
	BzzVectorInt indexTransfer;
	
	string	*names;
	string	*reactions;
	string	 reacting_system;

	void readFromFile(const string file_name);
	BzzVector GiveMeGrid();
	BzzVector GiveMeVelocity();
	BzzVector GiveMeTemperature();
	BzzVector GiveMeSpeciesMassFraction(const string name);
	BzzVector GiveMeSpeciesMoleFraction(const std::string name);
	BzzVector GiveMeSpeciesMassFraction(const int index);
	BzzVector GiveMeSpeciesMoleFraction(const int index);

	void GetMassFractionSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, const string name, BzzMatrix &SMatrix, BzzVectorInt &indices);
	void GetTemperatureSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, BzzMatrix &SMatrix, BzzVectorInt &indices);
	void GetFlameSpeedSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, BzzMatrix &SMatrix, BzzVectorInt &indices);
	void GetPressureCurvatureSensitivityProfiles(const bool iTotal, const bool iLocal, const double coordinate, BzzMatrix &SMatrix, BzzVectorInt &indices);

	void GetMassFractionSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, const string name, BzzVector &SVector, BzzVectorInt &indices);
	void GetTemperatureSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, BzzVector &SVector, BzzVectorInt &indices);
	void GetFlameSpeedSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, BzzVector &SVector, BzzVectorInt &indices);
	void GetPressureCurvatureSensitivityBars(const bool iTotal, const bool iLocal, const double coordinate, BzzVector &SVector, BzzVectorInt &indices);

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_SENSITIVITYANALYSIS_FLAME1D_POSTPROCESSOR

