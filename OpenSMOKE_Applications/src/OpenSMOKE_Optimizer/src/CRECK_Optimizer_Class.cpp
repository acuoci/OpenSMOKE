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

#include <cmath>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include "addons/OpenSMOKE_InverseKinetics.h"
#include "CRECK_Optimizer_Class.h"
#include "version.h"

inline bool my_isnan(double x)
{
	return x != x;
} 

void CRECK_Optimizer_Class::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Optimizer"		<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
    exit(-1);
}

void CRECK_Optimizer_Class::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_Optimizer"	<< endl;
    cout << "Object: "		<< name_object	<< endl;
    cout << "Warning:  "	<< message      << endl;
    cout << "Press a key to continue... "   << endl;
    getchar();
}

CRECK_Optimizer_Class::CRECK_Optimizer_Class()
{
	name_object					= "[Name not assigned]";
	bestObjectiveFunction		= 1.e16;
	countOnlinePostProcessing	= 100;
	iOnlinePostProcessing		= false;
}

void CRECK_Optimizer_Class::SetName(const string name)
{
	name_object = name;
}

void CRECK_Optimizer_Class::SetModel(const int imodel)
{
	iModel = imodel;
}

void CRECK_Optimizer_Class::SetOnlinePostProcessing()
{
	iOnlinePostProcessing = true;
}

void CRECK_Optimizer_Class::SetCountOnlinePostProcessing(const int count)
{
	countOnlinePostProcessing = count;
}

void CRECK_Optimizer_Class::setup(const string nameDetailed, const string globalName, const int index_parameter)
{
    int k;

    cout.setf(ios::scientific);

    cout << "" << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "                      OpenSMOKE_Optimizer                      " << endl;
	cout << "                  Version 0.3 - November 2009                  " << endl;
	cout << "             Alberto Cuoci - Politecnico di Milano             " << endl;
	cout << "                     alberto.cuoci@polimi.it                   " << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "" << endl;

    string nameOptimizationTable = "OptimizationTable.inp";

    // Kinetic scheme setup - Global kinetic scheme
    mix.SetupBinary(nameDetailed);
    global.assign_mix(&mix);
    global.read_from_file(globalName);

	if (index_parameter <=0)	SetupKineticParametersOptimization(nameOptimizationTable);
	else						SetupKineticParametersOptimization(nameOptimizationTable, index_parameter);

	global.SetupOptimization(iKind, iReaction, iSpecies);

    // Species names
    names = new string[mix.NumberOfSpecies()+1];
    for(k=1;k<=mix.NumberOfSpecies(); k++)
        names[k] = mix.names[k];

	// !!!!!! TODO !!!!!!
	{
		ChangeDimensions(mix.NumberOfSpecies(), &Norm_C);

		double cReference = 101325./8314/1000.;

		if (mix.recognize_species_without_exit("CH4") > 0)
			Norm_C[mix.recognize_species("CH4")]	= cReference * 0.10;	// [kmol/m3]
		if (mix.recognize_species_without_exit("O2") > 0)
			Norm_C[mix.recognize_species("O2")]		= cReference * 0.20;	// [kmol/m3]
		if (mix.recognize_species_without_exit("H2") > 0)
			Norm_C[mix.recognize_species("H2")]		= cReference * 0.10;	// [kmol/m3]
		if (mix.recognize_species_without_exit("H2O") > 0)
			Norm_C[mix.recognize_species("H2O")]	= cReference * 0.15;	// [kmol/m3]
		if (mix.recognize_species_without_exit("CO2") > 0)
			Norm_C[mix.recognize_species("CO2")]	= cReference * 0.15;	// [kmol/m3]
		if (mix.recognize_species_without_exit("CO") > 0)
			Norm_C[mix.recognize_species("CO")]		= cReference * 0.10;	// [kmol/m3]
		if (mix.recognize_species_without_exit("N2") > 0)
			Norm_C[mix.recognize_species("N2")]		= cReference * 0.70;	// [kmol/m3]
		if (mix.recognize_species_without_exit("OH") > 0)
			Norm_C[mix.recognize_species("OH")]		= cReference * 0.01;	// [kmol/m3]
		if (mix.recognize_species_without_exit("O") > 0)
			Norm_C[mix.recognize_species("O")]		= cReference * 0.01;	// [kmol/m3]
		if (mix.recognize_species_without_exit("H") > 0)
			Norm_C[mix.recognize_species("H")]		= cReference * 0.01;	// [kmol/m3]
		if (mix.recognize_species_without_exit("C2H4") > 0)
			Norm_C[mix.recognize_species("C2H4")]	= cReference * 0.01;	// [kmol/m3]
		if (mix.recognize_species_without_exit("C3H8") > 0)
			Norm_C[mix.recognize_species("C3H8")]	= cReference * 0.01;	// [kmol/m3]
		if (mix.recognize_species_without_exit("JET-A(G)") > 0)
			Norm_C[mix.recognize_species("JET-A(G)")]	= cReference * 0.01;	// [kmol/m3]
	}

    // Experiments setup
    RecognizeExperiments();
    inletStream = new OpenSMOKE_GasStream[nExperiments+1];

    // Setup of reactors
    pfr						= new OpenSMOKE_PFR[nExperiments+1];
    cstr					= new OpenSMOKE_CSTR[nExperiments+1];
	cfdf					= new OpenSMOKE_Flame1D[nExperiments+1];
	premix					= new OpenSMOKE_Flame1D[nExperiments+1];
	flamelet				= new OpenSMOKE_Flamelet[nExperiments+1];
	
	data_flame1D			= new OpenSMOKE_Flame1D_DataManager[nExperiments+1];
	operations_flame1D		= new OpenSMOKE_Flame1D_ScheduleClass[nExperiments+1];
	flame_speed_manager		= new OpenSMOKE_Flame1D_FlameSpeedManager[nExperiments+1];
	
	data_flamelet			= new OpenSMOKE_Flamelet_DataManager[nExperiments+1];
	operations_flamelet		= new OpenSMOKE_Flamelet_ScheduleClass[nExperiments+1];
	
	
	original_flame1d_solution	= new OpenSMOKE_Flame1D_Solution[nExperiments+1];
	original_flamelet_solution	= new OpenSMOKE_Flamelet_Solution[nExperiments+1];
	best_flamelet_solution		= new OpenSMOKE_Flamelet_Solution[nExperiments+1];

	OpenSMOKE_2EModel soot2EModel;

    for(k=1;k<=nExperiments; k++)
    { 
		stringstream reactor_number; reactor_number << k; 
		if (experiments[k].kindOfReactor == "PFR")
        {
            BzzVector solution;			// Solution

            string stream_name  = experiments[k].folderName + "/Stream.inp";
            string pfr_name     = experiments[k].folderName + "/PFR.inp";

            // Inlet stream setup
            inletStream[k].AssignKineticScheme(mix);
            inletStream[k].DefineFromFile(stream_name);

            // Plug flow reactor setup
			pfr[k].SetName("R" + reactor_number.str());
			pfr[k].UnsetVerbose();
			pfr[k].AssignSoot2EModel(soot2EModel);
            pfr[k].AssignKineticScheme(mix);
            pfr[k].AssignGlobalKineticScheme(global);
            pfr[k].AssignInletFlows(inletStream[k]);
            pfr[k].DefineFromFile(pfr_name);
            pfr[k].SetHistory();
            pfr[k].SetGlobalKinetics();
			pfr[k].Solve();
        }

		if (experiments[k].kindOfReactor == "FLAMESPEEDCURVE")
        {
			string command;
			string tab;
			#if LINUX_SO==1
				tab = "/";
			#else
				tab = "\\";
			#endif

			command     = "mkdir " + experiments[k].folderName + tab + "Steady";
			system(command.c_str());

			command     = "mkdir " + experiments[k].folderName + tab + "Additional";
			system(command.c_str());

			// Prepare flame
			premix[k].Assign(&mix);
			premix[k].Assign(&global);
			premix[k].Assign(&data_flame1D[k]);
			premix[k].Assign(&operations_flame1D[k]);
			premix[k].data->iGlobalKinetics = true;

			// Prepare data
			data_flame1D[k].Assign(&mix);
			data_flame1D[k].Assign(&premix[k]);
			
			premix[k].data->kind_of_flame = FLAME1D_PHYSICS_PREMIXED;
			
			// Prepare operations
			operations_flame1D[k].readOperations(experiments[k].folderName + "/Operations.inp");

			// Initialize flame
			premix[k].data->iBackUp = true;
			premix[k].FoldersAndFilesManager(experiments[k].folderName);
			premix[k].recoverFromBackUp(premix[k].nameFileBackupInputData);
			original_flame1d_solution[k].PasteFromExternalSolution(premix[k], *premix[k].data);
			premix[k].Run();

			// Flame speed curve
			premix[k].data->iFlameSpeedAnalysis = true;
			premix[k].data->flameSpeedAnalysisFileName = experiments[k].folderName + "/FlameSpeed.inp";
			flame_speed_manager[k].SetupFromFile(&premix[k], premix[k].data);
			flame_speed_manager[k].Run();
        }
		
		if (experiments[k].kindOfReactor == "CSTR")
        {
            BzzVector solution;			// Solution

            string stream_name  = experiments[k].folderName + "/Stream.inp";
            string cstr_name    = experiments[k].folderName + "/CSTR.inp";

            // Inlet stream setup
            inletStream[k].AssignKineticScheme(mix);
            inletStream[k].DefineFromFile(stream_name);

            // CSTR flow reactor setup
			cstr[k].SetName("CSTR_" + reactor_number.str());
			cstr[k].UnsetVerbose();
			cstr[k].AssignSoot2EModel(soot2EModel);
            cstr[k].AssignKineticScheme(mix);
            cstr[k].AssignGlobalKineticScheme(global);
            cstr[k].AssignInletFlows(inletStream[k]);
            cstr[k].DefineFromFile(cstr_name);
			cstr[k].SetHistoryPartial();
			cstr[k].SetGlobalKinetics();
			cstr[k].Solve();
        }

        if (experiments[k].kindOfReactor == "CFDF" || experiments[k].kindOfReactor == "CFDF-TWIN")
        {
			string command;
			string tab;
			#if LINUX_SO==1
				tab = "/";
			#else
				tab = "\\";
			#endif


			command     = "mkdir " + experiments[k].folderName + tab + "Steady";
			system(command.c_str());

			command     = "mkdir " + experiments[k].folderName + tab + "Additional";
			system(command.c_str());

			
			// Prepare flame
			cfdf[k].Assign(&mix);
			cfdf[k].Assign(&global);
			cfdf[k].Assign(&data_flame1D[k]);
			cfdf[k].Assign(&operations_flame1D[k]);
			cfdf[k].data->iGlobalKinetics = true;

			if (experiments[k].kindOfReactor == "CFDF-TWIN")
			{
				cfdf[k].data->kind_of_subphysics	= FLAME1D_SUBPHYSICS_NONE;
				cfdf[k].data->kind_of_flame			= FLAME1D_PHYSICS_TWIN;
			}

			// Prepare data
			data_flame1D[k].Assign(&mix);
			data_flame1D[k].Assign(&cfdf[k]);
			
			// Prepare operations
			operations_flame1D[k].readOperations(experiments[k].folderName + "/Operations.inp");

			// Initialize flame
			cfdf[k].data->iBackUp = true;
			cfdf[k].FoldersAndFilesManager(experiments[k].folderName);
			cfdf[k].recoverFromBackUp(cfdf[k].nameFileBackupInputData);
			cfdf[k].Run();

			// Run Flame
			//original_flame1d_solution[k].PasteFromExternalSolution(cfdf[k], *cfdf[k].data);
        }

		if (experiments[k].kindOfReactor == "FLAMELET")
        {
			string command;
			string tab;
			#if LINUX_SO==1
				tab = "/";
			#else
				tab = "\\";
			#endif


		//	command     = "mkdir " + experiments[k].folderName + tab + "Steady";
		//	system(command.c_str());

		//	command     = "mkdir " + experiments[k].folderName + tab + "Additional";
		//	system(command.c_str());

			
			// Prepare flame
			flamelet[k].Assign(&mix);
			flamelet[k].Assign(&global);
			flamelet[k].Assign(&data_flamelet[k]);
			flamelet[k].Assign(&operations_flamelet[k]);
			flamelet[k].data->iGlobalKinetics = true;

			// Prepare data
			data_flamelet[k].Assign(&mix);
			data_flamelet[k].Assign(&flamelet[k]);
			
			// Prepare operations
			data_flamelet[k].ReadFromFile(experiments[k].folderName + "/BackUp/Flamelet.inp");
			operations_flamelet[k].ReadOperations(experiments[k].folderName + "/Operations.inp");

			// Initialize flame
			flamelet[k].data->iBackUp = true;
			flamelet[k].FoldersAndFilesManager(experiments[k].folderName);
			flamelet[k].RecoverFromBackUp(flamelet[k].nameFileBackupInputData);
	
			// Run Flame
			original_flamelet_solution[k].PasteFromExternalSolution(flamelet[k]);
			flamelet[k].Run();
			best_flamelet_solution[k].PasteFromExternalSolution(flamelet[k]);
        }
	}

	BzzVector temperatures;
	for(k=1;k<=nExperiments; k++)
    { 
		BzzVector temp = experiments[k].ExtractTemperatures();
		temperatures.Append(temp);
		cout << temperatures.Size() << endl;
		cout << temperatures.Max() << endl;
		cout << temperatures.Min() << endl;

		Tmax = temperatures.Max();
		Tmin = temperatures.Min();
	}

	openOutputFileAndControl(fOutputParameters, "Output_Optimization/ActualParameters.out");
	fOutputParameters.setf(ios::scientific);
	fOutputParameters.precision(16);

	openOutputFileAndControl(fBestParameters,	"Output_Optimization/BestParameters.out");
	fBestParameters.setf(ios::scientific);
	fBestParameters.precision(16);

	openOutputFileAndControl(fObjectiveFunction,	"Output_Optimization/ObjectiveFunction.out");
	fObjectiveFunction.setf(ios::scientific);
	fObjectiveFunction.precision(16);

	fLog.open("Log.out", ios::out);
	fLog.setf(ios::scientific);


	global_index = 0;

	AssignTemperatures(temperatures);
}

double CRECK_Optimizer_Class::Minimization(BzzVector &optimizerParameters)
{
    double residual;
	int check_feasibility;
    BzzVector Tau_History;
    BzzVector Eta_History;
    BzzVector T_History;

	global_index++;

	BzzVector kineticParameters(nParameters);
	FromOptimizerParametersToKineticParameters(optimizerParameters, kineticParameters);
	
	// Check feasibility
	check_feasibility = CheckFeasibilityLimits(optimizerParameters);
	if(check_feasibility == 1)
	{
		cout << "Unfeasible parameters..." << endl;
		bzzUnfeasible = 1;
		return 0.;
	}

	// Update kinetic parameters
	global.ChangeKineticParameters(kineticParameters);
	
	// Objective function
	double residual_pfr			= 0.;
	double residual_cfdf		= 0.;
	double residual_twin		= 0.;
	double residual_cstr		= 0.;
	double residual_flamelet	= 0.;
	double residual_premix		= 0.;

	fLog << "*****************************************************************" << endl;
	fLog << " CALL " << global_index << endl;
	fLog << "*****************************************************************" << endl;

	for(int k=1;k<=nExperiments; k++)
	{
	    if (experiments[k].kindOfReactor == "PFR")
		{
			pfr[k].ReSolve();

			if (pfr[k].status >=0)
			{
				Tau_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].Tau_History);
				Eta_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].Csi_History);
				T_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].T_History);
				BzzMatrix mass_History(pfr[k].countGlobalIterations, mix.NumberOfSpecies(), pfr[k].mass_History); // submatrix of B										*
				BzzMatrix mole_History(pfr[k].countGlobalIterations, mix.NumberOfSpecies(), pfr[k].mole_History); // submatrix of B										*

				fObjective[k] = experiments[k].GiveMeObjectiveFunction(fLog, names, Tau_History, Eta_History, T_History, dummy_v, mass_History, mole_History);
			}
			else
				fObjective[k] = 1.e3;

			residual_pfr += fObjective[k].GetSumElements();
		}

	    if (experiments[k].kindOfReactor == "CSTR")
		{
			cstr[k].ReSolve();

			if (cstr[k].status >= 0)
			{
				double dummy_v_scalar;
				BzzVector aux1=cstr[k].mass_History.GetRow(1);
				BzzVector aux2=cstr[k].mole_History.GetRow(1);
				fObjective[k] = experiments[k].GiveMeObjectiveFunction(fLog, names, cstr[k].T_History[1], dummy_v_scalar, aux1, aux2);
			}
			else
				fObjective[k] = 1.e3;

			residual_cstr += fObjective[k].GetSumElements();
		}

		if (experiments[k].kindOfReactor == "CFDF" || experiments[k].kindOfReactor == "CFDF-TWIN")
        {
			int attempt = 1;
			global.default_set();
			for(;;)
			{
				cout << "Recover From backup..." << endl;
				cfdf[k].recoverFromBackUp(cfdf[k].nameFileBackupInputData);
				cfdf[k].Run();
				cout << "End Flame: " << k << endl;

				BzzVector cfdf_v(cfdf[k].U.Size());
				for (int m=1;m<=cfdf[k].U.Size();m++)
					cfdf_v[m] = 200.*cfdf[k].U[m]/cfdf[k].rho[m];	// 2.=nGeometry-1.

				BzzVector grid_m = cfdf[k].grid.x;
				fObjective[k] = experiments[k].GiveMeObjectiveFunction(fLog, names, grid_m, grid_m, cfdf[k].T, cfdf_v, cfdf[k].W, cfdf[k].X);

				
//				break;
				bool iMustChanged = false;
				{
					size_t found;
					stringstream fObjectiveStringStream;
					string       fObjectiveString;
					fObjectiveStringStream << fObjective[k].GetSumElements();
					fObjectiveString = fObjectiveStringStream.str();

					found = fObjectiveString.find("N");
					if (found != string::npos)
						iMustChanged = true;

					found = fObjectiveString.find("n");
					if (found != string::npos)
						iMustChanged = true;

					found = fObjectiveString.find("#");
					if (found != string::npos)
						iMustChanged = true;

					found = fObjectiveString.find("D");
					if (found != string::npos)
						iMustChanged = true;
				}

				if (iMustChanged == true)
					global.change_set();
				else
					break;

				attempt++;
			}

			if (experiments[k].kindOfReactor == "CFDF")			residual_cfdf += fObjective[k].GetSumElements();
			if (experiments[k].kindOfReactor == "CFDF-TWIN")	residual_twin += fObjective[k].GetSumElements();

		}

		if (experiments[k].kindOfReactor == "FLAMELET")
        {
			flamelet[k].PasteFromExternalSolution(original_flamelet_solution[k]);
			fLog << " Call start flame: " << k << endl;
            flamelet[k].ReRun();
			fLog << " Call end Flame:   " << k << endl;

			BzzVector grid_m = flamelet[k].grid.x;
			fObjective[k] = experiments[k].GiveMeObjectiveFunction(fLog, names, grid_m, grid_m, flamelet[k].T, dummy_v, flamelet[k].w, flamelet[k].X);

			residual_flamelet += fObjective[k].GetSumElements();
		}

		if (experiments[k].kindOfReactor == "FLAMESPEEDCURVE")
        {
			premix[k].PasteFromExternalSolution(original_flame1d_solution[k]);
			fLog << " Call start flame: " << k << endl;
            premix[k].Run();
			fLog << " Call end Flame:   " << k << endl;

			// Flame speed curve
			fLog << " Call start flame curve: " << k << endl;
			flame_speed_manager[k].Run();
			fLog << " Call end flame curve: " << k << endl;
				
			fLog << " Flame speed curve..." << endl;
			for (int i=1;i<=flame_speed_manager[k].list_of_phi.Size();i++)
				fLog << "  " << flame_speed_manager[k].list_of_phi[i] << 
						"  " << flame_speed_manager[k].list_of_flame_speeds[i] << endl;
			
			fLog << " Objective function calculations..." << endl;
			fObjective[k] = experiments[k].GiveMeObjectiveFunction(fLog, flame_speed_manager[k].list_of_phi, flame_speed_manager[k].list_of_flame_speeds);
			residual_premix += fObjective[k].GetSumElements();
		}

		fLog << endl;
	}

	residual = residual_cstr + residual_pfr + residual_twin + residual_cfdf + residual_flamelet + residual_premix;

	fLog << " Residual:          " << residual << endl;
	fLog << " Residual cstr:     " << residual_cstr << endl;
	fLog << " Residual pfr:      " << residual_pfr << endl;
	fLog << " Residual twin:     " << residual_twin << endl;
	fLog << " Residual cfdf:     " << residual_cfdf << endl;
	fLog << " Residual cfdf:     " << residual_premix << endl;
	fLog << " Residual flamelet: " << residual_flamelet << endl;
	for(int k=1;k<=nExperiments; k++)
		fLog << " fObj_" << k << " : " << fObjective[k].GetSumElements() << endl;
	fLog << endl;

	// Actual parameters
	fOutputParameters << global_index << "\t" << residual << "\t";
	for (int p=1;p<=nParameters;p++)
		fOutputParameters << kineticParameters[p] << "\t";
	fOutputParameters << endl;;

	// Best Parameters
	if (global_index == 1)
		originalObjectiveFunction = residual;
	if (residual < bestObjectiveFunction)
	{
		fLog << " Best solution found..." << endl;

		int p;
		bestOptimizerParameters = optimizerParameters;
		bestKineticParameters	= kineticParameters;
		bestObjectiveFunction	= residual;
		fBestParameters << global_index << "\t" << bestObjectiveFunction << "\t" << bestObjectiveFunction/originalObjectiveFunction << "\t";
		for (p=1;p<=nParameters;p++)
			fBestParameters << kineticParameters[p] << "\t";
		fBestParameters << endl;
		if (global_index == 1)
		{
			fObjectiveFunction << "Index(1)" << "\t" << "OF(2)" << "\t" << "OF_CFDF(3)" << "\t" << "OF_TWIN(4)" << "\t" << "OF_PFR(5)" << "\t" << "OF_FLAMELETS(6)" << "\t" << "OF_CSTR(7)" << "\t" << "OF_FS(8)" << "\t";
			
			int count = 9;
			for (p=1;p<=nExperiments;p++)
				fObjectiveFunction << "OF_Exp_" << p << "(" << count++ << ")\t";
			for (p=1;p<=nExperiments;p++)
				for (int i=1;i<=fObjective[p].Size();i++)
					fObjectiveFunction << "OF_Exp_" << p << "_" << i << "(" << count++ << ")\t";
			fObjectiveFunction << endl << endl;
		}

		fObjectiveFunction	<< global_index << "\t"	<< bestObjectiveFunction/originalObjectiveFunction                << "\t"
							<< residual_cfdf/bestObjectiveFunction << "\t" << residual_twin/bestObjectiveFunction     << "\t"
							<< residual_pfr/bestObjectiveFunction  << "\t" << residual_flamelet/bestObjectiveFunction << "\t"
							<< residual_cstr/bestObjectiveFunction << "\t" << residual_premix/bestObjectiveFunction   << "\t";
		for (p=1;p<=nExperiments;p++)
			fObjectiveFunction << fObjective[p].GetSumElements()/bestObjectiveFunction << "\t";
		for (p=1;p<=nExperiments;p++)
			for (int i=1;i<=fObjective[p].Size();i++)
				fObjectiveFunction << fObjective[p][i]/bestObjectiveFunction << "\t";

		fObjectiveFunction << endl;

		// WriteOptimizedKineticScheme();
	}

	// Online post-processing
	if (iOnlinePostProcessing == true && global_index % countOnlinePostProcessing == 1)
	{
		fLog << " Start post-processing..." << endl;

		// Post processing
		if (global_index == 1)	PostProcessing("Start", optimizerParameters);
		PostProcessing("online-post-processing", bestOptimizerParameters);
		WriteOptimizedKineticSchemeNEW();
		Gnuplot_plots_online();

		fLog << " End post-processing..." << endl;

	}

//	if (global_index == 34)
//	{
//		cout << "Pause" << endl;
//		getchar();
//	}

	return residual;//residual;
}

void CRECK_Optimizer_Class::Regression(int model, int ex, BzzVector &optimizerParameters, BzzVector &x, BzzVector &y)
{
	global_index++;

	BzzVector kineticParameters(nParameters);
	FromOptimizerParametersToKineticParameters(optimizerParameters, kineticParameters);

	global.ChangeKineticParameters(kineticParameters);

	// Write on file
	if (ex == 1)
	{
		cout << global_index << " - " << global_index/(number_of_experimental_points[nExperiments+1]-1) << endl;

		fOutputParameters << global_index << "\t";
		for (int p=1;p<=nParameters;p++)
			fOutputParameters << kineticParameters[p] << "\t";
		fOutputParameters << endl;
	}

	ChangeDimensions(0, &y);
	for(int k=1;k<=nExperiments; k++)
	{
		if (experiments[k].kindOfReactor == "PFR")
		{
			if (ex == number_of_experimental_points[k])
			{
				pfr[k].ReSolve();

				experiments[k].Tau_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].Tau_History);
				experiments[k].Eta_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].Csi_History);
				experiments[k].T_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].T_History);
				
				// TODO
				BzzMatrix B(pfr[k].countGlobalIterations, mix.NumberOfSpecies(), pfr[k].mole_History); // submatrix of B
				experiments[k].mole_History = B;
			}

			if (ex >= number_of_experimental_points[k] && ex < number_of_experimental_points[k+1])
			{
				BzzVector expected_y;
				experiments[k].GiveMeRegressionFunction(names, experiments[k].Tau_History, experiments[k].Eta_History, experiments[k].T_History, dummy_v, experiments[k].mole_History, x, expected_y);
				y.Append(expected_y);
			}
		}

		if (experiments[k].kindOfReactor == "CSTR")
		{
			if (ex == number_of_experimental_points[k])
			{
				cstr[k].ReSolve();

				experiments[k].Tau_History.GetBzzVector(cstr[k].countGlobalIterations, 1, cstr[k].Tau_History);
				experiments[k].Eta_History.GetBzzVector(cstr[k].countGlobalIterations, 1, cstr[k].Csi_History);
				experiments[k].T_History.GetBzzVector(cstr[k].countGlobalIterations, 1, cstr[k].T_History);
				
				// TODO
				BzzMatrix B(cstr[k].countGlobalIterations, mix.NumberOfSpecies(), cstr[k].mole_History); // submatrix of B
				experiments[k].mole_History = B;
			}

			if (ex >= number_of_experimental_points[k] && ex < number_of_experimental_points[k+1])
			{
				BzzVector expected_y;
				experiments[k].GiveMeRegressionFunction(names, experiments[k].Tau_History, experiments[k].Eta_History, experiments[k].T_History, dummy_v, experiments[k].mole_History, x, expected_y);
				y.Append(expected_y);
			}
		}

		if (experiments[k].kindOfReactor == "CFDF" || experiments[k].kindOfReactor == "CFDF-TWIN")
        {
			if (ex == number_of_experimental_points[k])
			{
				cfdf[k].PasteFromExternalSolution(original_flame1d_solution[k]);
				cfdf[k].ReRun();
			}

			if (ex >= number_of_experimental_points[k] && ex < number_of_experimental_points[k+1])
			{
				BzzVector cfdf_v(cfdf[k].U.Size());
				for (int m=1;m<=cfdf[k].U.Size();m++)
					cfdf_v[m] = 200.*cfdf[k].U[m]/cfdf[k].rho[m];	// 2.=nGeometry-1.


				BzzVector grid_m = cfdf[k].grid.x;
				BzzVector expected_y;
				experiments[k].GiveMeRegressionFunction(names, grid_m, grid_m, cfdf[k].T, cfdf_v, cfdf[k].X, x, expected_y);
				y.Append(expected_y);
			}
        }

		if (experiments[k].kindOfReactor == "FLAMELET")
        {
			if (ex == number_of_experimental_points[k])
			{
				flamelet[k].PasteFromExternalSolution(original_flamelet_solution[k]);
				flamelet[k].ReRun();
			}

			if (ex >= number_of_experimental_points[k] && ex < number_of_experimental_points[k+1])
			{
				BzzVector grid_m = flamelet[k].grid.x;
				BzzVector expected_y;
				experiments[k].GiveMeRegressionFunction(names, grid_m, grid_m, flamelet[k].T, dummy_v, flamelet[k].X, x, expected_y);
				y.Append(expected_y);
			}
        }
	}

}

int CRECK_Optimizer_Class::GiveMeNumberOfParameters()
{
    return nParameters;
}

BzzVector CRECK_Optimizer_Class::GiveMeStartingParameters()
{
	BzzVector startingParameters(nParameters);
	FromKineticParametersToOptimizerParameters(startingValues, startingParameters);
    return startingParameters;
}

void CRECK_Optimizer_Class::PostProcessing(string flag, BzzVector &optimizerParameters)
{
    BzzVector Tau_History;
    BzzVector Eta_History;
    BzzVector T_History;

	BzzVector kineticParameters(nParameters);
	FromOptimizerParametersToKineticParameters(optimizerParameters, kineticParameters);
	global.ChangeKineticParameters(kineticParameters);

	for(int k=1;k<=nExperiments; k++)
	{
	    if (experiments[k].kindOfReactor == "PFR")
        {
			pfr[k].ReSolve();

			Tau_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].Tau_History);
			Eta_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].Csi_History);
			T_History.GetBzzVector(pfr[k].countGlobalIterations, 1, pfr[k].T_History);
			BzzMatrix mass_History(pfr[k].countGlobalIterations, mix.NumberOfSpecies(), pfr[k].mass_History); // submatrix of B										*
			BzzMatrix mole_History(pfr[k].countGlobalIterations, mix.NumberOfSpecies(), pfr[k].mole_History); // submatrix of B										*

			experiments[k].PostProcessing(flag, names, Tau_History, Eta_History, T_History, dummy_v, mass_History, mole_History);
		}

	    if (experiments[k].kindOfReactor == "CSTR")
        {
			cstr[k].ReSolve();
			BzzVector v1(2);
			BzzVector v2(2);
			BzzMatrix masses(2, mix.NumberOfSpecies());
			BzzMatrix moles(2, mix.NumberOfSpecies());
			
			v1[1]=0.;
			v1[2]=1.;
			v2[1]=cstr[k].T_History[1];
			v2[2]=cstr[k].T_History[1];

			masses.SetRow(1, cstr[k].mass_History.GetRow(1)); 
			masses.SetRow(2, cstr[k].mass_History.GetRow(1)); 
			moles.SetRow(1, cstr[k].mole_History.GetRow(1)); 
			moles.SetRow(2, cstr[k].mole_History.GetRow(1)); 

			experiments[k].PostProcessing(flag, names, v1, v1, v2, dummy_v, masses, moles);
		}

		// ALBERTO
		if (experiments[k].kindOfReactor == "CFDF" || experiments[k].kindOfReactor == "CFDF-TWIN")
        {
		//	cfdf[k].PasteFromExternalSolution(original_flame1d_solution[k]);
        //  cfdf[k].ReRun();

			cfdf[k].recoverFromBackUp(cfdf[k].nameFileBackupInputData);
			cfdf[k].Run();
			
			BzzVector cfdf_v(cfdf[k].U.Size());
			for (int m=1;m<=cfdf[k].U.Size();m++)
				cfdf_v[m] = 200.*cfdf[k].U[m]/cfdf[k].rho[m];	// 2.=nGeometry-1.

			BzzVector grid_m = cfdf[k].grid.x;
	        experiments[k].PostProcessing(flag, names, grid_m, grid_m, cfdf[k].T, cfdf_v, cfdf[k].W, cfdf[k].X);
        }

		if (experiments[k].kindOfReactor == "FLAMESPEEDCURVE")
        {	
			premix[k].PasteFromExternalSolution(original_flame1d_solution[k]);
            premix[k].Run();
			flame_speed_manager[k].Run();
	        experiments[k].PostProcessing(flag, flame_speed_manager[k].list_of_phi, flame_speed_manager[k].list_of_flame_T, flame_speed_manager[k].list_of_flame_speeds);
        }

		if (experiments[k].kindOfReactor == "FLAMELET")
        {
	
			flamelet[k].PasteFromExternalSolution(original_flamelet_solution[k]);
            flamelet[k].ReRun();
			
			BzzVector grid_m = flamelet[k].grid.x;
	        experiments[k].PostProcessing(flag, names, grid_m, grid_m, flamelet[k].T, dummy_v, flamelet[k].w, flamelet[k].X);
        }
	}

	if (flag != "online-post-processing")
	{
		{
									cout << endl;
									cout << "------------------------------------" << endl;
			if (flag == "Start")	cout << "  Starting values " << endl;
			if (flag == "Final")	cout << "  Final values    " << endl;
									cout << "------------------------------------" << endl;

			for (int p=1;p<=nParameters;p++)
				cout << "   P[" << p << "] = "  << kineticParameters[p] << endl;
		}

		{
									fOutputParameters << endl;
									fOutputParameters << "------------------------------------" << endl;
			if (flag == "Start")	fOutputParameters << "  Starting values " << endl;
			if (flag == "Final")	fOutputParameters << "  Final values    " << endl;
									fOutputParameters << "------------------------------------" << endl;

			for (int p=1;p<=nParameters;p++)
				fOutputParameters << "   P[" << p << "] = "  << kineticParameters[p] << endl;

		}
	}
}

void CRECK_Optimizer_Class::Gnuplot_plots()
{
	for(int k=1;k<=nExperiments; k++)
        experiments[k].Gnuplot_plots();

	if(nParameters == 1)
		Gnuplot_plots_residual_1D();
}

void CRECK_Optimizer_Class::Gnuplot_plots_online()
{
	for(int k=1;k<=nExperiments; k++)
        experiments[k].Gnuplot_plots_online();
}

void CRECK_Optimizer_Class::Latex_report()
{
    OpenSMOKE_LatexInterface latex;

    latex.setup();
    latex.new_section("Kinetic scheme");
    latex.add("HERE: Description of kinetic scheme");
    latex.new_section("Optimization summary");
    latex.add("HERE: Description of optimization procedure");
    latex.new_section("Experiments");
    latex.add("HERE: Experiments description");
    latex.new_section("Results");
    latex.new_page();

    BzzMatrix A(3,3); string titlesUp[4];string titlesLeft[4];
    titlesUp[1]     = "C1";
    titlesUp[2]     = "C2";
    titlesUp[3]     = "C3";
    titlesLeft[1]   = "R1";
    titlesLeft[2]   = "R2";
    titlesLeft[3]   = "R3";

    latex.include_table(titlesUp, titlesLeft, A);

    for(int k=1;k<=nExperiments; k++)
    {
        experiments[k].Latex_figure(latex);
    }

	if (nParameters == 1)	Latex_plots_residual_1D(latex);

    latex.close();

    latex.create_pdf();
}

bool FileExists(string strFilename)
{
    struct stat   stFileInfo;
    bool          blnReturn;
    int           intStat;

    intStat = stat(strFilename.c_str(),&stFileInfo);

    if(intStat == 0)
        blnReturn = true;
    else
        blnReturn = false;

    return(blnReturn);
}

void CRECK_Optimizer_Class::RecognizeExperiments()
{

    string nameFolder;
    int count;

    count = 0;
    for(;;)
    {
        stringstream number;

        count++;
        number << count;
        if (count<10)       nameFolder = "Exp0" + number.str();
        else                nameFolder = "Exp"  + number.str();

        if (!FileExists(nameFolder))
        {
            nExperiments = count-1;
            break;
        }
    }

    // Experiments setup
    experiments = new ExperimentClass[nExperiments+1];
	fObjective  = new BzzVector[nExperiments+1];

    for(int i=1;i<=nExperiments;i++)
    {
        stringstream number;
        number << i;

        if (i<10)   nameFolder = "Exp0" + number.str();
        else        nameFolder = "Exp"  + number.str();
        experiments[i].Setup(nameFolder, i);
    }
}

void CRECK_Optimizer_Class::GiveMeRegressionExperimentalData(BzzVector &x_experimental, BzzVector &y_experimental)
{
	ChangeDimensions(nExperiments+1, &number_of_experimental_points);
	number_of_experimental_points[1] = 1;

	for(int k=1;k<=nExperiments; k++)
	{
		int sum = 0;
		for (int i=1;i<=experiments[k].nData;i++)
		{
			BzzVector x = experiments[k].experimental_data[i].GetColumn(1);
			BzzVector y = experiments[k].experimental_data[i].GetColumn(2);
			x_experimental.Append(x);
			y_experimental.Append(y);

			sum += x.Size();
		}

		number_of_experimental_points[k+1] = number_of_experimental_points[k]+sum;
	}
}

void CRECK_Optimizer_Class::SetupKineticParametersOptimization(const string fileName)
{	
	int    idummy;	
	string  dummy;
	double  minValue;
	double  maxValue;
	double  startingValue;

	ifstream fInput;
	openInputFileAndControl(fInput, fileName.c_str());

	nParameters = 0;
	for(;;)
	{
		fInput >> dummy;
		if(dummy == "END")
			break;

		if(dummy != "REACTION")
			ErrorMessage("Expected REACTION key word!");

		fInput >> idummy;
		iReaction.Append(idummy);

		fInput >> dummy;
		if(dummy == "A")
		{
			iKind.Append(1);
			iSpecies.Append(0);
		
		}
		else if(dummy == "Tatt")
		{
			iKind.Append(2);
			iSpecies.Append(0);
		}
		else if(dummy == "Beta")
		{
			iKind.Append(3);
			iSpecies.Append(0);
		}
		else if(dummy == "lambda")
		{
			iKind.Append(4);
			
			fInput >> dummy;
			iSpecies.Append(mix.recognize_species(dummy));
		}
		else
			ErrorMessage("The available parameters are: A || Tatt || Beta || lambda");

		fInput >> startingValue;
		fInput >> minValue;
		fInput >> maxValue;
		if(minValue >= maxValue) ErrorMessage("The min value must be lower than the max value");
		if(startingValue >= maxValue || startingValue <= minValue) ErrorMessage("The starting value must be between the min and the max values");
		
		minValues.Append(minValue);
		maxValues.Append(maxValue);
		startingValues.Append(startingValue);

		nParameters++;
	}	

	fInput.close();	

	CheckKineticParametersOptimization();
}

void CRECK_Optimizer_Class::SetupKineticParametersOptimization(const string fileName, const int index)
{	
	int    idummy;	
	string dummy;
	double minValue;
	double maxValue;
	double startingValue;
	const int SIZE = 200;
	char comment[SIZE];

	ifstream fInput;
	openInputFileAndControl(fInput, fileName.c_str());

	nParameters = 0;
	for(;;)
	{
		if (nParameters+1 == index)
		{
			fInput >> dummy;
			cout << dummy << endl;
			if(dummy == "END")
				break;

			if(dummy != "REACTION")
				ErrorMessage("Expected REACTION key word!");

			fInput >> idummy;
			cout << idummy << endl;
			iReaction.Append(idummy);

			fInput >> dummy;
			cout << dummy << endl;
			if(dummy == "A")
			{
				iKind.Append(1);
				iSpecies.Append(0);
			
			}
			else if(dummy == "Tatt")
			{
				iKind.Append(2);
				iSpecies.Append(0);
			}
			else if(dummy == "Beta")
			{
				iKind.Append(3);
				iSpecies.Append(0);
			}
			else if(dummy == "lambda")
			{
				iKind.Append(4);
				
				fInput >> dummy;
				iSpecies.Append(mix.recognize_species(dummy));
			}
			else
				ErrorMessage("The available parameters are: A || Tatt || Beta || lambda");

			fInput >> startingValue;
			fInput >> minValue;
			fInput >> maxValue;
			if(minValue >= maxValue) ErrorMessage("The min value must be lower than the max value");
			if(startingValue >= maxValue || startingValue <= minValue) ErrorMessage("The starting value must be between the min and the max values");
			
			minValues.Append(minValue);
			maxValues.Append(maxValue);
			startingValues.Append(startingValue);

			break;
		}
		else
		{
			nParameters++;
			fInput.getline(comment, SIZE);
		}
	}	

	nParameters = 1;

	fInput.close();	

	CheckKineticParametersOptimization();
}

void CRECK_Optimizer_Class::CheckKineticParametersOptimization()
{
	int i;
	
	for(i=1;i<=nParameters;i++)
	{
		stringstream i_string;
		stringstream iReaction_string;

		i_string 		<< i;
		iReaction_string 	<< iReaction[i];


		if (iReaction[i] > global.nReactions)
		{
			string message 	 = "The Reaction " + iReaction_string.str() + " is not included in the Global Kinetics!\n";
			message 	+= "Check the Optimization Table file!";
			ErrorMessage(message);
		}
		
		if (iKind[i] == 1)	// A
		{
			string message = "Please check A min and max values in Reaction " 	+ iReaction_string.str() + 
				             " in Optimization Table file at request #" 		+ i_string.str();

			if (minValues[i] < 0.)		ErrorMessage(message);
			if (maxValues[i] > 1.e32)	ErrorMessage(message);
		}
		
		if (iKind[i] == 2)	// Tatt
		{
			string message = "Please check Tatt min and max values in Reaction " 	+ iReaction_string.str() + 
				         " in Optimization Table file at request #" 		+ i_string.str();

			if (minValues[i] < -1e8)	ErrorMessage(message);
			if (maxValues[i] >  1e8)	ErrorMessage(message);
		}

		if (iKind[i] == 3)	// Beta
		{
			string message = "Please check Beta min and max values in Reaction " 	+ iReaction_string.str() + 
				         " in Optimization Table file at request #" 		+ i_string.str();

			if (minValues[i] < -1e1)	ErrorMessage(message);
			if (maxValues[i] >  1e1)	ErrorMessage(message);
		}
		
		if (iKind[i] == 4)	// lambda
		{
			string message = "Please check lambda min and max values in Reaction " 	+ iReaction_string.str() + 
				         " in Optimization Table file at request #" 		+ i_string.str();

			if (minValues[i] < -1e1)	ErrorMessage(message);
			if (maxValues[i] >  1e1)	ErrorMessage(message);
		}		
	}
	
	// Reaction pointers
	reactionMirror = new ReactionMirror[global.nReactions+1];
	for(i=1;i<=global.nReactions;i++)
	{
		reactionMirror[i].A			= global.A[i];
		reactionMirror[i].Beta		= global.Beta[i];
		reactionMirror[i].Tatt		= global.Tatt[i];
		reactionMirror[i].lambda	= global.lambda.GetRow(i);
		reactionMirror[i].lambda_original	= global.lambda.GetRow(i);
		reactionMirror[i].Tatt_original		= global.Tatt[i];
	}
}

void CRECK_Optimizer_Class::FromOptimizerParametersToKineticParameters(BzzVector &optimizer_parameter,
																	   BzzVector &kinetic_parameter)
{
	int i;

	if (iModel == 1)
	{
		for(int i=1;i<=nParameters;i++)
		{	
			if (iKind[i] == 1)	// A
			{
				kinetic_parameter[i]			= optimizer_parameter[i];
				reactionMirror[iReaction[i]].A	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 2)	// Tatt
			{
				kinetic_parameter[i]				= optimizer_parameter[i];
				reactionMirror[iReaction[i]].Tatt	= kinetic_parameter[i]; 
			}
			
			else if (iKind[i] == 3)	// Beta
			{
				kinetic_parameter[i]				= optimizer_parameter[i];
				reactionMirror[iReaction[i]].Beta	= kinetic_parameter[i]; 
			}
			
			else if (iKind[i] == 4)	// lambda
			{
				kinetic_parameter[i]								= optimizer_parameter[i];
				reactionMirror[iReaction[i]].lambda[iSpecies[i]]	= kinetic_parameter[i]; 
			}
		}
	}

	else if (iModel == 2)
	{
		for(int i=1;i<=nParameters;i++)
		{	
			if (iKind[i] == 1)	// A
			{
				kinetic_parameter[i]			= exp(optimizer_parameter[i]);
				reactionMirror[iReaction[i]].A	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 2)	// Tatt
			{
				kinetic_parameter[i]				= Tmean * optimizer_parameter[i];
				reactionMirror[iReaction[i]].Tatt	= kinetic_parameter[i]; 
			}
			
			else if (iKind[i] == 3)	// Beta
			{
				kinetic_parameter[i]				= optimizer_parameter[i];
				reactionMirror[iReaction[i]].Beta	= kinetic_parameter[i]; 
			}
			
			else if (iKind[i] == 4)	// lambda
			{
				kinetic_parameter[i]								= optimizer_parameter[i];
				reactionMirror[iReaction[i]].lambda[iSpecies[i]]	= kinetic_parameter[i]; 
			}
		}
	}

	else if (iModel == 4)
	{
		for(i=1;i<=nParameters;i++)
		{	
			if (iKind[i] == 2)	// Tatt
			{
				kinetic_parameter[i]				= C1*optimizer_parameter[i];
				reactionMirror[iReaction[i]].Tatt	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 3)	// Beta
			{
				kinetic_parameter[i]				= optimizer_parameter[i];
				reactionMirror[iReaction[i]].Beta	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 4)	// lambda
			{
				kinetic_parameter[i]								= optimizer_parameter[i];
				reactionMirror[iReaction[i]].lambda[iSpecies[i]]	= kinetic_parameter[i]; 
			}
		}
		
		for(i=1;i<=nParameters;i++)
			if (iKind[i] == 1)	// A
			{
				double sum = 0.;

				for(int j=1;j<=mix.NumberOfSpecies();j++)
					sum += reactionMirror[iReaction[i]].lambda[j]*log(Norm_C[j]);

				kinetic_parameter[i]			= exp(optimizer_parameter[i] - C2*reactionMirror[iReaction[i]].Tatt -sum);
				reactionMirror[iReaction[i]].A	= kinetic_parameter[i]; 
			}
	}
	
	// Video information 
	#ifdef VERBOSE_OUT
		cout << "Parameter\tKinetic\t\t\tOptimizer (O->K)" << endl;
		for(i=1;i<=nParameters;i++)
			cout << i << "\t\t" << kinetic_parameter[i] << "\t\t" << optimizer_parameter[i] << endl;
	#endif
}


void CRECK_Optimizer_Class::FromKineticParametersToOptimizerParameters(BzzVector &kinetic_parameter,
																	   BzzVector &optimizer_parameter)
{
	int i;

	if (iModel == 1)
	{
		for(i=1;i<=nParameters;i++)
		{	
			if (iKind[i] == 1)	// A
			{
				optimizer_parameter[i]			= kinetic_parameter[i];
				reactionMirror[iReaction[i]].A	= kinetic_parameter[i]; 
			}

			if (iKind[i] == 2)	// Tatt
			{
				optimizer_parameter[i]				= kinetic_parameter[i];
				reactionMirror[iReaction[i]].Tatt	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 3)	// Beta
			{
				optimizer_parameter[i]				= kinetic_parameter[i];
				reactionMirror[iReaction[i]].Beta	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 4)	// lambda
			{
				optimizer_parameter[i]								= kinetic_parameter[i];
				reactionMirror[iReaction[i]].lambda[iSpecies[i]]	= kinetic_parameter[i]; 
			}
		}
	}

	else if (iModel == 2)
	{
		for(i=1;i<=nParameters;i++)
		{	
			if (iKind[i] == 1)	// A
			{
				optimizer_parameter[i]			= log(kinetic_parameter[i]);
				reactionMirror[iReaction[i]].A	= kinetic_parameter[i]; 
			}

			if (iKind[i] == 2)	// Tatt
			{
				optimizer_parameter[i]				= kinetic_parameter[i]/Tmean;
				reactionMirror[iReaction[i]].Tatt	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 3)	// Beta
			{
				optimizer_parameter[i]				= kinetic_parameter[i];
				reactionMirror[iReaction[i]].Beta	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 4)	// lambda
			{
				optimizer_parameter[i]								= kinetic_parameter[i];
				reactionMirror[iReaction[i]].lambda[iSpecies[i]]	= kinetic_parameter[i]; 
			}
		}
	}

	else if (iModel == 4)
	{
		for(i=1;i<=nParameters;i++)
		{	
			if (iKind[i] == 2)	// Tatt
			{
				optimizer_parameter[i]				= kinetic_parameter[i]/C1;
				reactionMirror[iReaction[i]].Tatt	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 3)	// Beta
			{
				optimizer_parameter[i]				= kinetic_parameter[i];
				reactionMirror[iReaction[i]].Beta	= kinetic_parameter[i]; 
			}

			else if (iKind[i] == 4)	// lambda
			{
				optimizer_parameter[i]								= kinetic_parameter[i];
				reactionMirror[iReaction[i]].lambda[iSpecies[i]]	= kinetic_parameter[i]; 
			}
		}

		for(i=1;i<=nParameters;i++)
			if (iKind[i] == 1)	// A
			{
				double sum = 0.;
				
				for(int j=1;j<=mix.NumberOfSpecies();j++)
					sum += reactionMirror[iReaction[i]].lambda[j]*log(Norm_C[j]);

				optimizer_parameter[i]	= log(kinetic_parameter[i]) + C2*reactionMirror[iReaction[i]].Tatt + sum;
				reactionMirror[iReaction[i]].A = kinetic_parameter[i]; 

			}
	}

	cout << "Parameter\tKinetic\t\t\tOptimizer (K->O)" << endl;
	for(i=1;i<=nParameters;i++)
		cout << i << "\t\t" << kinetic_parameter[i] << "\t\t" << optimizer_parameter[i] << endl;
}

void CRECK_Optimizer_Class::AssignTemperatures(BzzVector &temperatures)
{
	BzzVector uTemperatures(temperatures.Size());

	C2	= 1./Mean(temperatures);	
	
	for(int i=1;i<=temperatures.Size();i++)
		uTemperatures[i] = 1./temperatures[i] - C2;

	C1 = 1./uTemperatures.Norm2();	

	Tmean = 1500.;

	// Output on video
	cout << "Mean temperature: " << 1./C2	<<	" [K]"		<< endl;
	cout << "Constant C1:      " << C1		<<	" [K]"		<< endl;
	cout << "Constant C2:      " << C2		<<	" [1/K]"	<< endl;
	cout << "Constant C1*C2:   " << C1*C2   <<	" [-]"		<< endl;

	// Output on file
	fOutputParameters << "Model:   " << iModel                  << endl;
	fOutputParameters << "C1:      " << C1		<<	" [K]"		<< endl;
	fOutputParameters << "C2:      " << C2		<<	" [1/K]"	<< endl;
	fOutputParameters << "C1*C2:   " << C1*C2   <<	" [-]"		<< endl;
	fOutputParameters << endl;
	fOutputParameters << "Tmin:    " << temperatures.Min()	<< " K"		<< endl;
	fOutputParameters << "Tmax:    " << temperatures.Max()	<< " K"		<< endl;
	fOutputParameters << "TmeanS:  " << 1./C2				<< " K"		<< endl;
	fOutputParameters << "Tmean:   " << Tmean				<< " K"		<< endl;
	fOutputParameters << endl;
	fOutputParameters << endl;
}

void CRECK_Optimizer_Class::Gnuplot_plots_residual_1D()
{
    string fileFigure = "Residuals";
    string fileResidual = "Output_Optimization/Output.out";

    OpenSMOKE_GnuPlotInterface gplot("Output_Optimization");
    gplot.setPlot("Residual-1D", "parameter", "RESIDUAL");
    gplot.setKey("Residual");
    gplot.setKind('p');
	if(iKind[1] == 1) gplot.setLogScaleX();
    gplot.setLogScaleY();
    gplot.plot(fileFigure, fileResidual, 2, 3);
}

void CRECK_Optimizer_Class::Latex_plots_residual_1D(OpenSMOKE_LatexInterface &latex)
{
    string fileName = "Output_Optimization/Residuals.eps";
    string caption  = "Residuals";
    latex.include_figure(fileName, caption);
}

void CRECK_Optimizer_Class::Gnuplot_plots_one_per_time()
{
	for(int k=1;k<=nExperiments; k++)
        experiments[k].Gnuplot_plots();

	if(nParameters == 1)	Gnuplot_plots_residual_1D();
}

void CRECK_Optimizer_Class::Latex_report_one_per_time()
{
    OpenSMOKE_LatexInterface latex;

    latex.setup();

    for(int k=1;k<=nExperiments; k++)
        experiments[k].Latex_figure(latex);

	if (nParameters == 1)	Latex_plots_residual_1D(latex);

    latex.close();
    latex.create_pdf();
}

int CRECK_Optimizer_Class::CheckFeasibilityLimits(BzzVector &optimizerParameters)
{
	if (iModel == 1 || iModel == 4 || iModel == 2)
	{
		for(int i=1;i<=nParameters;i++)
		{	
			if (iKind[i] == 1)	// A
			{
				if (reactionMirror[iReaction[i]].A	> maxValues[i])	return 1;
				if (reactionMirror[iReaction[i]].A	< minValues[i])	return 1;
			}

			else if (iKind[i] == 2)	// Tatt
			{
				if (reactionMirror[iReaction[i]].Tatt > maxValues[i]) return 1; 
				if (reactionMirror[iReaction[i]].Tatt < minValues[i]) return 1; 
			}
			
			else if (iKind[i] == 3)	// Beta
			{
				if (reactionMirror[iReaction[i]].Beta > maxValues[i]) return 1; 
				if (reactionMirror[iReaction[i]].Beta < minValues[i]) return 1; 
			}
			
			else if (iKind[i] == 4)	// lambda
			{
				if (reactionMirror[iReaction[i]].lambda[iSpecies[i]] > maxValues[i]) return 1; 
				if (reactionMirror[iReaction[i]].lambda[iSpecies[i]] < minValues[i]) return 1; 
			}
		}
	}

	return 0;
}


BzzVector CRECK_Optimizer_Class::GiveMeMinParameters()
{
	BzzVector minParameters(nParameters);
	if (iModel == 1)
		minParameters = minValues;
	else if (iModel == 2)
	{
		for(int i=1;i<=nParameters;i++)
		{	
			if (iKind[i] == 1)	// A
				minParameters[i] = log(minValues[i]);
			else if (iKind[i] == 2)	// Tatt
				minParameters[i] = minValues[i]/Tmean;
			else if (iKind[i] == 3)	// Beta
				minParameters[i] = minValues[i];
			else if (iKind[i] == 4)	// lambda
				minParameters[i] = minValues[i];
		}
	}

    return minParameters;
}


BzzVector CRECK_Optimizer_Class::GiveMeMaxParameters()
{
	BzzVector maxParameters(nParameters);
	if (iModel == 1)
		maxParameters = maxValues;
	else if (iModel == 2)
	{
		for(int i=1;i<=nParameters;i++)
		{	
			     if (iKind[i] == 1)	// A
				maxParameters[i] = log(maxValues[i]);
			else if (iKind[i] == 2)	// Tatt
				maxParameters[i] = maxValues[i]/Tmean;
			else if (iKind[i] == 3)	// Beta
				maxParameters[i] = maxValues[i];
			else if (iKind[i] == 4)	// lambda
				maxParameters[i] = maxValues[i];
		}
	}

    return maxParameters;
}

void CRECK_Optimizer_Class::WriteOptimizedKineticScheme()
{
	string dummy_string;

	ofstream fOut;
	openOutputFileAndControl(fOut, "Output_Optimization/OptmizedGlobalKineticScheme.out");
	fOut.setf(ios::scientific);

	fOut << "NAME Opt-Scheme" << endl;
	fOut << endl;

	int i;
	for(i=1;i<=global.nReactions;i++)
	{
		int j;

		if (global.iEquilibrium[i] == 1) dummy_string = "Y";
		else							 dummy_string = "N";

		fOut << "   Reaction " << i << endl;
		fOut << endl;
		fOut << "   A           " << setprecision(16) << reactionMirror[i].A		<< endl;
		fOut << "   Beta        " << setprecision(16) << reactionMirror[i].Beta	<< endl;
		fOut << "   Tatt        " << setprecision(16) << reactionMirror[i].Tatt	<< endl;
		fOut << "   Equilibrium " << dummy_string				<< endl;
		fOut << endl;
		
		fOut << "   Stoichiometry" << endl;
		for(j=1;j<=global.NC;j++)
			if (global.nu[i][j] != 0.)
				fOut << "   "<< setw(12) << left << mix.names[j] << " " << global.nu[i][j] << endl;;
		fOut << endl;
		
		fOut << "   Kinetics" << endl;
		for(j=1;j<=global.NC;j++)
			if (reactionMirror[i].lambda[j] != 0.)
				fOut << "   " << setw(12) << left << mix.names[j] << " " << setprecision(16) << reactionMirror[i].lambda[j] << endl;;
		fOut << endl << endl;
	}

	fOut << "END" << endl;

	// Optimized Kinetic scheme
	OpenSMOKE_GlobalKinetics	global_optimized;
	BzzVector AInverse;
	BzzVector BetaInverse;
	BzzVector TattInverse;
	BzzMatrix lambdaInverse;
	BzzVectorInt indicesInverse;
	{
		global_optimized.assign_mix(&mix);
		global_optimized.read_from_file("Output_Optimization/OptmizedGlobalKineticScheme.out");
		
		for(int i=1;i<=global_optimized.nReactions;i++)
		{
			if (global_optimized.iEquilibrium[i] == 1)
			{
				int j;
				double A, Beta, Tatt;
				BzzVector lambda;
				double Patm = 1.;	// TODO

				OpenSMOKE_InverseKinetics inverse;
				inverse.AssignKineticScheme(mix);
				inverse.AssignGlobalKineticScheme(global_optimized);
				inverse.Setup(i, Patm, Tmin, Tmax+300, (Tmax+300-Tmin)/60.,
							  A, Beta, Tatt, lambda);

				indicesInverse.Append(i);
				AInverse.Append(A);
				BetaInverse.Append(Beta);
				TattInverse.Append(Tatt);
				lambdaInverse.AppendRow(lambda);

				fOut << "   Reaction " << i << endl;
				fOut << endl;
				fOut << "   A           " << A		<< endl;
				fOut << "   Beta        " << Beta	<< endl;
				fOut << "   Tatt        " << Tatt	<< endl;
				fOut << "   Equilibrium " << "N"	<< endl;
				fOut << endl;
				
				fOut << "   Stoichiometry" << endl;
				for(j=1;j<=global.NC;j++)
					if (global.nu[i][j] != 0.)
						fOut << "   "<< setw(12) << left << mix.names[j] << " " << -global.nu[i][j] << endl;;
				fOut << endl;
				
				fOut << "   Kinetics" << endl;
				for(j=1;j<=global.NC;j++)
					if (lambda[j] != 0.)
						fOut << "   " << setw(12) << left << mix.names[j] << " " << lambda[j] << endl;;
				fOut << endl << endl;
			}
		}
	}

	fOut.close();

	{
		int j;

		openOutputFileAndControl(fOut, "Output_Optimization/FluentGlobalKineticScheme.out");

		fOut << "/ ************************************************************************ /" << endl;
		fOut << "/                            Kinetic scheme                                /" << endl;
		fOut << "/ ************************************************************************ /" << endl;
		fOut << endl;
		fOut << "nChi       1" << endl;
		fOut << "nReactions "  << global_optimized.nReactions + indicesInverse.Size() << endl;
		fOut << "nSpecies   "  << global_optimized.NC << endl;
		fOut << endl;

		fOut << "Species" << endl;
		for(j=1;j<=global_optimized.NC;j++)
			fOut << setw(12) << mix.names[j] <<  left << mix.M[j] << endl;
		fOut << endl;

		fOut << "Chi" << endl;
		fOut << 0. << endl;
		fOut << endl;

		fOut << "Stoichiometry" << endl;
		for(j=1;j<=global_optimized.NC;j++)
		{
			fOut << setw(12) << left << mix.names[j];
			for(i=1;i<=global_optimized.nReactions;i++)
				fOut << setw(8) << global_optimized.nu[i][j];
			for(i=1;i<=indicesInverse.Size();i++)
				fOut << setw(8) << -global_optimized.nu[indicesInverse[i]][j];
			fOut << endl;
		}
		fOut << endl;
		
		fOut << "A" << endl;
		for(i=1;i<=global_optimized.nReactions;i++)
			fOut << setw(17) << global_optimized.A[i];
		for(i=1;i<=indicesInverse.Size();i++)
			fOut << setw(17) << AInverse[i];
		fOut << endl << endl;

		fOut << "Beta" << endl;
		for(i=1;i<=global_optimized.nReactions;i++)
			fOut << setw(17) << global_optimized.Beta[i];
		for(i=1;i<=indicesInverse.Size();i++)
			fOut << setw(17) << BetaInverse[i];
		fOut << endl << endl;

		fOut << "Tatt" << endl;
		for(i=1;i<=global_optimized.nReactions;i++)
			fOut << setw(17) << global_optimized.Tatt[i];
		for(i=1;i<=indicesInverse.Size();i++)
			fOut << setw(17) << TattInverse[i];
		fOut << endl << endl;

		for(j=1;j<=global_optimized.NC;j++)
		{
			fOut << "lambda " << mix.names[j] << endl;
			for(i=1;i<=global_optimized.nReactions;i++)
				fOut << setw(10) << global_optimized.lambda[i][j];
			for(i=1;i<=indicesInverse.Size();i++)
				fOut << setw(10) << lambdaInverse[i][j];
			fOut << endl << endl;
		}

		fOut.close();
	}
}


void CRECK_Optimizer_Class::WriteOptimizedKineticSchemeNEW()
{
	string dummy_string;

	ofstream fOut;
	openOutputFileAndControl(fOut, "Output_Optimization/OptmizedGlobalKineticScheme.out");
	fOut.setf(ios::scientific);

	fOut << "NAME Opt-Scheme" << endl;
	fOut << endl;

	int i;
	for(i=1;i<=global.nReactions;i++)
	{
		int j;

		if (global.iEquilibrium[i] == 1) dummy_string = "Y";
		else							 dummy_string = "N";

		fOut << "   Reaction " << i << endl;
		fOut << endl;
		fOut << "   A           " << setprecision(16) << global.A[i]	<< endl;
		fOut << "   Beta        " << setprecision(16) << global.Beta[i]	<< endl;
		fOut << "   Tatt        " << setprecision(16) << global.Tatt[i]	<< endl;
		fOut << "   Equilibrium " << dummy_string				<< endl;
		fOut << endl;
		
		fOut << "   Stoichiometry" << endl;
		for(j=1;j<=global.NC;j++)
			if (global.nu[i][j] != 0.)
				fOut << "   "<< setw(12) << left << mix.names[j] << " " << global.nu[i][j] << endl;;
		fOut << endl;
		
		fOut << "   Kinetics" << endl;
		for(j=1;j<=global.NC;j++)
			if (global.lambda[i][j] != 0.)
				fOut << "   " << setw(12) << left << mix.names[j] << " " << setprecision(16) << global.lambda[i][j] << endl;;
		fOut << endl << endl;
	}

	fOut << "END" << endl;


	// Optimized Kinetic scheme
	OpenSMOKE_GlobalKinetics	global_optimized;
	BzzVector AInverse;
	BzzVector BetaInverse;
	BzzVector TattInverse;
	BzzMatrix lambdaInverse;
	BzzVectorInt indicesInverse;
	{
		{
			global_optimized.assign_mix(&mix);
			global_optimized.read_from_file("Output_Optimization/OptmizedGlobalKineticScheme.out");
			
			for(int i=1;i<=global_optimized.nReactions;i++)
			{
				if (global_optimized.iEquilibrium[i] == 1)
				{
					int j;
					double A, Beta, Tatt;
					BzzVector lambda;
					double Patm = 1.;	// TODO

					OpenSMOKE_InverseKinetics inverse;
					inverse.AssignKineticScheme(mix);
					inverse.AssignGlobalKineticScheme(global_optimized);
					inverse.Setup(i, Patm, Tmin, Tmax+300, (Tmax+300-Tmin)/60.,
								  A, Beta, Tatt, lambda);

					indicesInverse.Append(i);
					AInverse.Append(A);
					BetaInverse.Append(Beta);
					TattInverse.Append(Tatt);
					lambdaInverse.AppendRow(lambda);

					fOut << "   Reaction " << i << endl;
					fOut << endl;
					fOut << "   A           " << A		<< endl;
					fOut << "   Beta        " << Beta	<< endl;
					fOut << "   Tatt        " << Tatt	<< endl;
					fOut << "   Equilibrium " << "N"	<< endl;
					fOut << endl;
					
					fOut << "   Stoichiometry" << endl;
					for(j=1;j<=global.NC;j++)
						if (global.nu[i][j] != 0.)
							fOut << "   "<< setw(12) << left << mix.names[j] << " " << -global.nu[i][j] << endl;;
					fOut << endl;
					
					fOut << "   Kinetics" << endl;
					for(j=1;j<=global.NC;j++)
						if (lambda[j] != 0.)
							fOut << "   " << setw(12) << left << mix.names[j] << " " << lambda[j] << endl;;
					fOut << endl << endl;
				}
			}
		}
	}

	fOut.close();

	{
		int j;
	
		ofstream fOutFLUENT;
		openOutputFileAndControl(fOutFLUENT, "Output_Optimization/FluentGlobalKineticScheme.out");

		fOutFLUENT << "/ ************************************************************************ /" << endl;
		fOutFLUENT << "/                            Kinetic scheme                                /" << endl;
		fOutFLUENT << "/ ************************************************************************ /" << endl;
		fOutFLUENT << endl;
		fOutFLUENT << "nChi       1" << endl;
		fOutFLUENT << "nReactions "  << global_optimized.nReactions + indicesInverse.Size() << endl;
		fOutFLUENT << "nSpecies   "  << global_optimized.NC << endl;
		fOutFLUENT << endl;

		fOutFLUENT << "Species" << endl;
		for(j=1;j<=global_optimized.NC;j++)
			fOutFLUENT << setw(12) << left << mix.names[j] <<  left << mix.M[j] << endl;
		fOutFLUENT << endl;

		fOutFLUENT << "Chi" << endl;
		fOutFLUENT << 0. << endl;
		fOutFLUENT << endl;

		fOutFLUENT << "Stoichiometry" << endl;
		for(j=1;j<=global_optimized.NC;j++)
		{
			fOutFLUENT << setw(12) << left << mix.names[j];
			for(i=1;i<=global_optimized.nReactions;i++)
				fOutFLUENT << setw(8) << global_optimized.nu[i][j];
			for(i=1;i<=indicesInverse.Size();i++)
				fOutFLUENT << setw(8) << -global_optimized.nu[indicesInverse[i]][j];
			fOutFLUENT << endl;
		}
		fOutFLUENT << endl;
		
		fOutFLUENT << "A" << endl;
		for(i=1;i<=global_optimized.nReactions;i++)
			fOutFLUENT << setw(17) << global_optimized.A[i];
		for(i=1;i<=indicesInverse.Size();i++)
			fOutFLUENT << setw(17) << AInverse[i];
		fOutFLUENT << endl << endl;

		fOutFLUENT << "Beta" << endl;
		for(i=1;i<=global_optimized.nReactions;i++)
			fOutFLUENT << setw(17) << global_optimized.Beta[i];
		for(i=1;i<=indicesInverse.Size();i++)
			fOutFLUENT << setw(17) << BetaInverse[i];
		fOutFLUENT << endl << endl;

		fOutFLUENT << "Tatt" << endl;
		for(i=1;i<=global_optimized.nReactions;i++)
			fOutFLUENT << setw(17) << global_optimized.Tatt[i];
		for(i=1;i<=indicesInverse.Size();i++)
			fOutFLUENT << setw(17) << TattInverse[i];
		fOutFLUENT << endl << endl;

		for(j=1;j<=global_optimized.NC;j++)
		{
			fOutFLUENT << "lambda " << mix.names[j] << endl;
			for(i=1;i<=global_optimized.nReactions;i++)
				fOutFLUENT << setw(10) << global_optimized.lambda[i][j];
			for(i=1;i<=indicesInverse.Size();i++)
				fOutFLUENT << setw(10) << lambdaInverse[i][j];
			fOutFLUENT << endl << endl;
		}

		fOutFLUENT.close();
	}

}