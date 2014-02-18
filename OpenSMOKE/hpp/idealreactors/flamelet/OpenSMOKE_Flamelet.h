/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_FLAMELET)
#define OPENSMOKE_FLAMELET

#include "OpenSMOKE.hpp"
#include "basic/OpenSMOKE_Grid1D.h"
#include "OpenSMOKE_Flamelet_DataManager.h"
#include "OpenSMOKE_Flamelet_ScheduleClass.h"
#include "OpenSMOKE_Flamelet_ODE_Objects.h"

class OpenSMOKE_Flamelet_Solution;

class OpenSMOKE_Flamelet
{
	friend class OpenSMOKE_Flamelet_Solution;

public:

	OpenSMOKE_Flamelet();
	void SetName(const string name);

	void Assign(OpenSMOKE_Flamelet_DataManager *_data);
	void Assign(OpenSMOKE_Flamelet_ScheduleClass *_operations);
	void Assign(OpenSMOKE_ReactingGas *_mix);
	void Assign(OpenSMOKE_GlobalKinetics *_global);

	void Run();
	void ReRun();
	void RunMemorySaved();
	void PasteFromExternalSolution(OpenSMOKE_Flamelet_Solution &solution);

	OpenSMOKE_Grid1D					grid;
	OpenSMOKE_Flamelet_DataManager		*data;
	OpenSMOKE_Flamelet_ScheduleClass	*operations;
	OpenSMOKE_GlobalKinetics			*global;
	OpenSMOKE_ReactingGas				*mix;

	void FoldersAndFilesManager();
	void FoldersAndFilesManager(const string backupFolder);
	void RecoverFromBackUp(const string fileName);

	string nameOutputFolder; 
	string nameFolderBackupData;
	string nameFileBackupInputData;
	string nameFolderUnsteadyData;
	string nameFolderSteadyData;
	string nameFolderAdditionalData;
	string nameFolderProfilesData;
	string nameFileOutput;


private:	// Internal variables

	int		NC;
	int		NR;
	string	name_object;

	void Setup();
	void Solve(const double tEnd);
	void NewPoints(string kind, const char index);
	void DoubleTheGrid();
	void RefineGridPeak(const double fraction);
	void RefineLeanSide(const double fraction);

private:	// Radiation

	double Tenv4;
	int iH2O;
	int iCO2;
	int iCO;
	int iCH4;

	BzzVector asTot;
	BzzVector Qrad;

	void prepare_radiation();
	void calculate_radiation();

private:	// Derivatives

	BzzVector dT_over_dt;
	BzzMatrix dw_over_dt;
	
	BzzVector dT_over_dx;
	BzzVector dCp_over_dx;
	BzzMatrix dw_over_dx;

	BzzVector d2T_over_dx2;
	BzzMatrix d2w_over_dx2;	

	BzzVector hProfile;
	BzzVector hDefectProfile;
	BzzMatrix omega_elemental_profile;


public:		// Main variables

	BzzVector T;
	BzzMatrix w;
	BzzMatrix X;
	BzzVector c;


private:	// Ode system

	MyOdeSystem_Flamelet				odeFlamelet;
	MyOdeSystem_Flamelet_Enthalpy		odeFlamelet_Enthalpy;
	MyOdeSystem_Flamelet_EnthalpyDefect	odeFlamelet_EnthalpyDefect;
	MyOdeSystem_Flamelet_Soot			odeFlamelet_Soot;

	BzzVector initialValues;
	BzzVector xMin;
	BzzVector xMax;

	int dimTot;
	int dimBlock;

	BzzInverseErrorFunction ErfInv;


private:	// Grid

	int		Np;
	int		Ni;
	void	AddPoints(BzzVectorInt &listPoints);

	
private:	// Internal functions
	
	void Prepare();
	void InitialTemperatureProfile();
	void InitialMassFractionProfiles();
	void ErrorMessage(const string message);
	void WarningMessage(const string message);	


private:	// Properties

	BzzVector PMtot;
	BzzVector uPMtot;
	BzzVector rho;	
	BzzMatrix R;
	BzzMatrix Cpk;
	BzzVector Cp;
	BzzVector chi;
	BzzVector QReaction;
	BzzVector enthalpy;

	BzzMatrix omega_elemental;
	BzzMatrix x_elemental;

	void   ChiEvaluation(const double chi0);
	double StoichiometricChiEvaluation(const double chi0);
	void Properties();
	void Properties(const int jacobianIndex, BzzVectorInt &jacobianVariables, const int dimBlock);
	void MoleFractionsAndPMtot();
	void MassFractionsAndPMtot();


public:		// Pseudo-public functions

	void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	void GetSystemFunctions_Enthalpy(BzzVector &x,double t,BzzVector &f);
	void GetSystemFunctions_EnthalpyDefect(BzzVector &x,double t,BzzVector &f);
	void GetSystemFunctions_Soot(BzzVector &y,double t,BzzVector &dy);
	void MyODE_Print(BzzVector &x, double t);


private:	// Auxiliary variables

	BzzVector xVector;
	BzzVector wVector;
	BzzVector RVector;

	BzzMatrix CpMap;
	BzzMatrix k1Map;
	BzzMatrix k2Map;
	BzzMatrix uKeqMap;
	BzzMatrix logFcentMap;
	BzzMatrix reactionDSMap;
	BzzMatrix reactionDHMap;


private:	// Print on file

	ofstream fGnuPlot;
	ofstream fGnuPlotODE;
	ofstream fGnuPlotSoot;

	void Video_label();
	void Video_label_soot();
	void GnuPlot_label(const string name, const int N);
	void GnuPlotODE_label(const int N);
	void PrintGnuPlotODE(const double t);
	void PrintGnuPlot();

	void ElementalAnalysis();

private:	// Memory allocation

	void Allocate();
	void Allocate_Np_Variables();
	void Allocate_Master_Variables();


private:	// Soot model
	
	OpenSMOKE_2EModel	sootModel;
	BzzVector		phiN;
	BzzVector		phiM;
	BzzVector		dphiN_over_dt;
	BzzVector		dphiM_over_dt;
	BzzVector		source_phiN;
	BzzVector		source_phiM;
	BzzVector		d2phiN_over_dx2;
	BzzVector		d2phiM_over_dx2;
	BzzMatrix		SootGasCorrection;
	BzzVector		DiffusionSoot;

	void Initial_conditions_soot_module();	
};

class OpenSMOKE_Flamelet_Solution
{
public:
	int Np;
	BzzVector x;
	BzzVector U;
	BzzVector G;
	BzzVector H;
	BzzVector T;
	BzzMatrix W;
	BzzVector phiN;
	BzzVector phiM;
	BzzMatrix moments;

	void PasteFromExternalSolution(OpenSMOKE_Flamelet &flame);
};

#endif // OPENSMOKE_FLAMELET
