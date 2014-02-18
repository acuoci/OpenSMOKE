/***************************************************************************
 *   Copyright (C) 2009 by Tiziano Maffei and Alberto Cuoci  	           *
 *   tiziano.maffei@chem.polimi.it   				                       *
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

#ifndef OPENSMOKE_MYSOLIDREGRESSION
#define OPENSMOKE_MYSOLIDREGRESSION

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Dictionary.h"
#include "engine/OpenSMOKE_ReactingGas.h"

enum solid_regression_type {OPENSMOKE_REGRESSION_CHAR, OPENSMOKE_REGRESSION_DSMOKE_BIO};

class MyOpenSMOKE_SolidRegression;
class OpenSMOKE_SolidExperiment;
class OpenSMOKE_SolidExperiment_DSmoke_Bio;
class OpenSMOKE_CharKineticScheme;

class OpenSMOKE_Dictionary_SolidRegression : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_SolidRegression();
};

class MyOdeSystem : public BzzOdeSystemObject
{
public:
	MyOpenSMOKE_SolidRegression *ptMyRegression;
	void assign(MyOpenSMOKE_SolidRegression *myregression);

public:
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
};


class MyOpenSMOKE_SolidRegression
{
private:

	OpenSMOKE_CharKineticScheme *kinetics;

	MyOdeSystem ode;
	BzzOdeStiffObject o;
	
	void ErrorMessage(const string message);
	void WarningMessage(const string message);

	void PrepareOptimizerFile();
	void CheckDictionary(OpenSMOKE_Dictionary_SolidRegression &dictionary);
	void WriteOptimizedParametersChar(const string file_name);
	void WriteOptimizedParametersBio(const string file_name);
	
	BzzMatrix YY;
	BzzVector yy;

	ofstream fOutput;
	string name_object;

	int nCases;
	int numTotal;

	BzzVectorInt    indices;
	BzzVector initialconditions_x;
	BzzMatrix initialconditions_y;
	BzzVector yMin;
	BzzVector yMax;

	void Prepare();
	BzzVector temperatures_global;
	BzzVector utemperatures_global;
	BzzVector pressures_Gas_global;

	double TMeanGlobal;
	double PGasMeanGlobal;

	int iModel;
	solid_regression_parameters parameters;
	int nParameters;
	BzzMatrixInt bUserDefined;

public:
	
	double C1;
	double C2;
	double C3;

	// Constructors
	MyOpenSMOKE_SolidRegression();
	void AssignKinetics(OpenSMOKE_CharKineticScheme *_kinetics);
	void AssignKinetics(OpenSMOKE_ReactingGas *_kinetics);
	void AssignListOfFiles(const string _file_name_list_of_files);

	void SetName(const string name);
	void SetModel(const int _iModel);
	void SetRegressionAnalysis(const bool index);
	void SetSurfaceEquation(const bool index);
	void SetRobustAnalysis(const bool index);
	void SetMinimumRatio_A(const double _value);
	void SetMaximumRatio_A(const double _value);
	void SetMinimumRatio_E(const double _value);
	void SetMaximumRatio_E(const double _value);
	void SetMinimumDelta_BETA(const double _value);
	void SetMaximumDelta_BETA(const double _value);
	void SetVerbose(const bool value);
	void SetRemoveWarning(const bool value);
	void SetAbsoluteTolerance(const double _absoluteTolerance);
	void SetRelativeTolerance(const double _relativeTolerance);
	void SetParameters(const vector<string> names);
	void SetSolidSpecies(const vector<string> _names);
	void SetSolidDensity(const string units, const double value);

	// Experiments
	vector<string> list_of_names_of_files;
	OpenSMOKE_SolidExperiment			 *experiments_char;
	OpenSMOKE_SolidExperiment_DSmoke_Bio *experiments_bio;

	// Preparation
	void Setup(const string file_name);
	void Run_Char();
	void Run_DSmoke_Bio();
	
	// Regression functions
	void ModelOdeRegression_Char(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);
	void ModelOdeRegression_DSmoke_Bio(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);
	
	void MyODE_Print_Char(BzzVector &x, double t);
	void MyODE_Print_DSmoke_Bio(BzzVector &x, double t);
	
	// Ode models
	void GetModel_01(BzzVector &x,double t,BzzVector &f);
	void GetModel_02(BzzVector &x,double t,BzzVector &f);

	BzzVector teta;
	BzzVector Cadsorbed;
	BzzVector Cbulk;
	double xchar;
	double Sg;

	int ptIndex;
	double ptT;
	double ptPressure;

	double volume_total_solid;
	double mass_total_solid;
	double rho_solid;

	BzzVector R;
	BzzVector C;
	BzzVector mass;
	BzzVector dmass_over_dt;
	BzzVectorInt indices_solid;
	vector<string> names_solid;

	OpenSMOKE_ReactingGas *mix_chemkin;

	void MemoryAllocation();

	double dxchar_over_dt;
	double dSg_over_dt;
	BzzVector dteta_over_dt;
	BzzVector teta_f;
	BzzVector C_f;
	BzzVector xcarbon;
	BzzVector dxcarbon_over_dt;

	int count_global;
	bool iWriteOnFile;
	void LabelODE_File();
	ofstream fOptimizer;

	solid_regression_type solid_regression;

private:

	static const double MWchar;

	bool iRegressionAnalysis;
	bool iSurfaceEquation;
	bool iRobustAnalysis;
	bool iRelativeTolerance;
	bool iAbsoluteTolerance;
	bool iVerbose;

	string file_name_list_of_files;
	double relativeTolerance;
	double absoluteTolerance;

	double minimum_ratio_A;
	double maximum_ratio_A;
	double minimum_ratio_E;
	double maximum_ratio_E;
	double minimum_delta_BETA;
	double maximum_delta_BETA;
};


#endif // OPENSMOKE_SOLIDREGRESSION


