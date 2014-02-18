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

class MyOpenSMOKE_SolidRegression;
class OpenSMOKE_SolidExperiment;
class OpenSMOKE_CharKineticScheme;

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

public:
	
	double C1;
	double C2;
	double C3;

	// Constructors
	MyOpenSMOKE_SolidRegression();
	void SetName(const string name);
	void SetKinetics(OpenSMOKE_CharKineticScheme *_kinetics);

	// Experiments
	vector<string> list_of_names_of_files;
	OpenSMOKE_SolidExperiment *experiments;

	// Preparation
	void Setup(const string filename);
	void Run(const int _iModel, BzzVector &bFirstGuess);
	
	// Regression functions
	void ModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);
	void MyODE_Print(BzzVector &x, double t);
	
	// Ode models
	void GetModel_01(BzzVector &x,double t,BzzVector &f);

	// TODO
	int NR;
	int Nadsorbed;
	int Nbulk;
	int Ngas;

	BzzVector teta;
	BzzVector Cadsorbed;
	BzzVector Cbulk;
	double xchar;
	double Sg;

	BzzVector kappa0;
	BzzVector Ea;
	BzzVector kappa;
	BzzVector eta;

	BzzVector rr;
	BzzMatrix nuCarbon;
	BzzMatrix nuBulk;
	BzzMatrix nuGas;
	BzzMatrix nuAdsorbed;
	BzzVector nuChar;

	BzzVector lambdaCarbon;
	BzzMatrix lambdaBulk;
	BzzMatrix lambdaGas;
	BzzMatrix lambdaAdsorbed;

	BzzVector R;
	double Rchar;


	double ptT;
	double ptSg0;	
	double ptSigma;	
	BzzVector ptGas;	

	static const double MWchar;
	static const double psi;

	void PrepareKineticScheme();

	double dxchar_over_dt;
	double dSg_over_dt;
	BzzVector dteta_over_dt;

	int count_global;
	bool iWriteOnFile;
	void LabelODE_File();
	ofstream fLog;

};


#endif // OPENSMOKE_SOLIDREGRESSION