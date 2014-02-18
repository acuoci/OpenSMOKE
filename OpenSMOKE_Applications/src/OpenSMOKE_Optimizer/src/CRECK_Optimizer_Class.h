/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci							   *
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

#ifndef CRECK_OPTIMIZER_CLASS_H
#define CRECK_OPTIMIZER_CLASS_H

#include "OpenSMOKE.hpp"
#include "idealreactors/pfr/OpenSMOKE_PFR.h"
#include "idealreactors/cstr/OpenSMOKE_CSTR.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_DataManager.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_ScheduleClass.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D_FlameSpeedManager.h"
#include "idealreactors/flamelet/OpenSMOKE_Flamelet.h"
#include "idealreactors/flamelet/OpenSMOKE_Flamelet_DataManager.h"
#include "idealreactors/flamelet/OpenSMOKE_Flamelet_ScheduleClass.h"
#include "ExperimentClass.h"

class ReactionMirror
{
public:
	double A;
	double Beta;
	double Tatt;
	BzzVector lambda;

	BzzVector lambda_original;
	double Tatt_original;
};

class CRECK_Optimizer_Class
{
public:

	CRECK_Optimizer_Class();
	void SetName(const string name);
	void SetModel(const int imodel);
	void SetOnlinePostProcessing();
	void SetCountOnlinePostProcessing(const int count);

    ExperimentClass				*experiments;
    OpenSMOKE_ReactingGas		mix;				// Gas mixture (reactive)
    OpenSMOKE_GlobalKinetics	global;
    OpenSMOKE_GasStream			*inletStream;

	void setup(const string nameDetailed, const string globalName, const int index_parameter);
    double Minimization(BzzVector &newParameter);
    void PostProcessing(string flag, BzzVector &parameters);
    void Gnuplot_plots();
    void Gnuplot_plots_online();
    void Latex_report();
	void Gnuplot_plots_one_per_time();
	void Latex_report_one_per_time();

    int GiveMeNumberOfParameters();
    BzzVector GiveMeStartingParameters();
	BzzVector GiveMeMinParameters();
	BzzVector GiveMeMaxParameters();

	BzzVector Norm_C;

	void Regression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);
	void GiveMeRegressionExperimentalData(BzzVector &x_experimental, BzzVector &y_experimental);

private:

    OpenSMOKE_PFR							*pfr;
    OpenSMOKE_CSTR							*cstr;

	OpenSMOKE_Flame1D						*cfdf;
	OpenSMOKE_Flame1D_DataManager			*data_flame1D;
	OpenSMOKE_Flame1D_ScheduleClass			*operations_flame1D;
	OpenSMOKE_Flame1D_Solution				*original_flame1d_solution;
//	OpenSMOKE_Flame1D_Solution				*best_flame1d_solution;

	OpenSMOKE_Flame1D						*premix;

	OpenSMOKE_Flamelet						*flamelet;
	OpenSMOKE_Flamelet_DataManager			*data_flamelet;
	OpenSMOKE_Flamelet_ScheduleClass		*operations_flamelet;
	OpenSMOKE_Flamelet_Solution				*original_flamelet_solution;
	OpenSMOKE_Flamelet_Solution				*best_flamelet_solution;

	OpenSMOKE_Flame1D_FlameSpeedManager     *flame_speed_manager;

    string *names;
    int global_index;
    int nExperiments;
    BzzVector *fObjective;

    void RecognizeExperiments();
	
	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);

	ofstream fOutputParameters;
	ofstream fBestParameters;
	ofstream fObjectiveFunction;
	ofstream fLog;

	void SetupKineticParametersOptimization(const string fileName);
	void SetupKineticParametersOptimization(const string fileName, const int index);
	void CheckKineticParametersOptimization();
	void AssignTemperatures(BzzVector &temperatures);
	void FromKineticParametersToOptimizerParameters(BzzVector &optimizer_parameter, BzzVector &kinetic_parameter);
	void FromOptimizerParametersToKineticParameters(BzzVector &kinetic_parameter, BzzVector &optimizer_parameter);
	double	C1;
	double	C2;
	double  Tmean;
	int		iModel;
	int nParameters;

	ReactionMirror *reactionMirror;
	BzzVector startingValues;
	BzzVector minValues;
	BzzVector maxValues;
	BzzVectorInt 	iKind;
	BzzVectorInt 	iSpecies;
	BzzVectorInt 	iReaction;
	void Gnuplot_plots_residual_1D();
	void Latex_plots_residual_1D(OpenSMOKE_LatexInterface &latex);

	BzzVector dummy_v;

	BzzVectorInt number_of_experimental_points;

	int CheckFeasibilityLimits(BzzVector &optimizerParameters);

	bool iOnlinePostProcessing;
	int countOnlinePostProcessing;
	BzzVector bestOptimizerParameters;
	BzzVector bestKineticParameters;
	double bestObjectiveFunction;
	double originalObjectiveFunction;

	void WriteOptimizedKineticScheme();
	void WriteOptimizedKineticSchemeNEW();
	double Tmin;
	double Tmax;
};

#endif // CRECK_OPTIMIZER_CLASS_H
