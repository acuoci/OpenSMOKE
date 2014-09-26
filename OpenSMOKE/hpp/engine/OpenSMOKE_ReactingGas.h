/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci   	   *
 *   alberto.cuoci@polimi.it   						   *
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

#ifndef OPENSMOKE_REACTINGGAS_H
#define OPENSMOKE_REACTINGGAS_H

#include "Linux_Windows.h"
#include "engine/OpenSMOKE_Kinetics.h"
#include "engine/OpenSMOKE_IdealGas.h"

class OpenSMOKE_ReactingGas : public OpenSMOKE_IdealGas
{
public:

	OpenSMOKE_ReactingGas();
	void SetName(const std::string name);

	OpenSMOKE_Kinetics kinetics;

	// Main functions
	int		NumberOfReactions();					// Number of reactions
	void	SetupBinary(const std::string pathName);		// Setup of kinetic scheme

	BzzVector	reactionDH;
	BzzVector	reactionDS;

	BzzVector		k1;
	BzzVector		k2;
	BzzVector		uKeq;
	BzzVector		rDirC;
	BzzVector		rInvC;
	BzzVector		logFcent;	
	BzzVector		r;
	BzzVector		coeffM;
	BzzVector		coeffFallOff;

	// Kinetic scheme Path
	std::string path_complete;
	std::string folder_path;
	std::string name_kinetic_scheme;
	char **strReaction;

	// Soot mode
	bool iSootMode;
	int  indexSoot;

	// Formation rates and reaction heat
	void ComputeKineticParameters(const double T, const double logT, const double uT, double P_Pa);
	void ComputeFromConcentrations(const double T, BzzVector &c, const double cTot, BzzVector *R);
	void ComputeFromConcentrations(double T, BzzVector &c, double cTot, BzzVector &rForward, BzzVector &rBackward);
	double ComputeQReaction(const double T);

	// Matrix of stoichiometric coefficient
	void GiveMeIndexOfSpeciesInReactions(const std::string fileName, BzzVectorIntArray &indices);
	void GiveMe_Jalfa_A(BzzMatrix &JAlfa, OpenSMOKE_NuManager *nu, const double T, const int indexSpecies, const int indexTemperature, BzzVector &parameters);
	void GiveMe_Jalfa_Beta(BzzMatrix &JAlfa, OpenSMOKE_NuManager *nu, const double T, const int indexSpecies, const int indexTemperature, BzzVector &parameters);
	void GiveMe_Jalfa_Eatt(BzzMatrix &JAlfa, OpenSMOKE_NuManager *nu, const double T, const int indexSpecies, const int indexTemperature, BzzVector &parameters);

	// Kinetic data maps
	void InitializeMap(int numPoints);
	void ComputeKineticParameters_map(BzzVector &T, BzzVector &ctot, BzzVector &P_Pa);
	void ComputeFromConcentrations_map(int kReactor, BzzVector &c, BzzVector *R);
	void ComputeFromConcentrations_map(BzzVector &c, BzzMatrix &R, BzzVector &QReaction);
	void CorrectKineticParameters_map(BzzMatrix &correction_k1, BzzMatrix &correction_uKeq, BzzMatrix &correction_k2);
	void CorrectKineticParameters(BzzVector &correction_k1,BzzVector &correction_uKeq, BzzVector &correction_k2);

	// Global reactions
	void ComputeKEq(const double T, BzzMatrix &nu, BzzVector &sumNuij, BzzVector &uKeq);

	// Print On Binary File
	void SaveOnBinaryFile(BzzSave &fSave);

	void ComputeFromConcentrationsAndCorrect(double T, BzzVector &c, double cTot, BzzVector *R, BzzVector &correction);

private:

	int	 NR;

public:

	BzzMatrix	k1_map;
	BzzMatrix	k2_map;
	BzzMatrix	uKeq_map;
	BzzMatrix	logFcent_map;
	BzzMatrix	reactionDH_map;
	BzzMatrix	reactionDS_map;
	BzzVector	T_map;
	BzzVector	cTot_map;
	BzzVector	P_Pa_map;
	BzzVector	RVector;

	void VerboseDataKinetics(ofstream &fOutput, ofstream &fOutputInverseKinetics);
	void VerboseDataKinetics(ofstream &fOutput, const int k, BzzVector &T_Vector, BzzMatrix &kappa_Matrix, BzzMatrix &uKappaEquilibrium_Matrix, BzzMatrix &DH_Matrix, BzzMatrix &DS_Matrix);
	void VerboseDataSummary(ofstream &fOutput, BzzVector &DH_Tref, BzzVector &DS_Tref);
	void WriteElementsFile(const std::string fileName);

	bool IsTransportModeAvailable();

private:

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};


#endif

