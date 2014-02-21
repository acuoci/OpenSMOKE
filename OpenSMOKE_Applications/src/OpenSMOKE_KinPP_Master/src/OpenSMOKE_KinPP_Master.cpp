/***************************************************************************
 *   Copyright (C) 2006-2008 by                                            *
 *   Guido Buzzi-Ferraris, Alessio Frassoldati and Alberto Cuoci           *
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

#include <sstream>
#include <ctime>
#include "BzzMath.hpp"
#include "interfaces/SimpleOpt.h"
#include "basic/OpenSMOKE_Utilities.h"

void ErrorMessage(const string message)
{ 
    cout << endl;
    cout << "Executable:  OpenSMOKE_KPP_Master"		<< endl;
    cout << "Error:       "  << message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_FLAG, _T("-test"),	         SO_NONE    }, // "-test"

    { OPT_ARG,  _T("--input"),           SO_REQ_SEP }, // "--input			ARG"
    { OPT_ARG,  _T("--kinetics"),        SO_REQ_SEP }, // "--kinetics		ARG"
	{ OPT_ARG,  _T("--correction"),      SO_REQ_SEP }, // "--correction		ARG"
	{ OPT_ARG,  _T("--clustering-exe"),  SO_REQ_SEP }, // "--clustering		ARG"
	{ OPT_ARG,  _T("--kpp-exe"),		 SO_REQ_SEP }, // "--clustering		ARG"
	{ OPT_ARG,  _T("--sequence"),        SO_REQ_SEP }, // "--sequence		ARG"
	{ OPT_ARG,  _T("--backup"),			 SO_REQ_SEP }, // "-backup			ARG"
	
	


	
	{ OPT_ARG,  _T("--relaxation"),      SO_REQ_SEP }, // "--relaxation		ARG"
	{ OPT_ARG,  _T("--max-correction"),  SO_REQ_SEP }, // "--maxcorrection	ARG"
	{ OPT_ARG,  _T("-jacobian"),         SO_NONE    }, // "--jacobian		ARG"
	{ OPT_ARG,  _T("-only-ode"),		 SO_NONE    }, // "-only-ode		ARG"
	{ OPT_ARG,  _T("--t-max-delta"),     SO_REQ_SEP }, // "--t-max-delta	ARG"
	{ OPT_ARG,  _T("--t-max-ud"),		 SO_REQ_SEP }, // "--t-max-ud   	ARG"
	{ OPT_ARG,  _T("--t-max-local"),     SO_REQ_SEP }, // "--t-max-local	ARG"
	{ OPT_ARG,  _T("--species-fluctuations"),    SO_REQ_SEP }, // "--species-fluctuations	ARG"
	{ OPT_ARG,  _T("--tol-rel"),		SO_REQ_SEP }, // "--tol-rel	ARG"
	{ OPT_ARG,  _T("--tol-abs"),		SO_REQ_SEP }, // "--tol-abs	ARG"
	{ OPT_ARG,  _T("--max-tol-rel"),	SO_REQ_SEP }, // "--max-tol-rel	ARG"
	{ OPT_ARG,  _T("--max-tol-abs"),	SO_REQ_SEP }, // "--max-tol-abs	ARG"
	{ OPT_ARG,  _T("--max-newtons"),	SO_REQ_SEP }, // "--max-newtons	ARG"
	{ OPT_ARG,  _T("--f1stop"),			SO_REQ_SEP }, // "--f1stop	ARG"

    { OPT_HELP, _T("-?"),    SO_NONE    }, // "-?"
    { OPT_HELP, _T("-help"), SO_NONE    }, // "--help"

    SO_END_OF_OPTIONS                       // END
};

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_KinPP [-option] [--option STRING] [-?] [--help]"								<< endl;

    cout << "Options: "																						<< endl;
    cout << "   -?                      this help"															<< endl;
    cout << "   --input FOLDERNAME      folder containing the case"											<< endl;
    cout << "   --correction STRING		correction coefficient: none || sin || dirac || beta || gaussian"	<< endl;
    cout << "   --clustering DOUBLE		clustering ratio"													<< endl;
    cout << "   --relaxation INTEGER	relaxation factor (default 0)"										<< endl;
	cout << "   --max-correction DOUBLE	max reaction rate increment (default 1.e10)"						<< endl;

	cout << "   --t-max-delta DOUBLE	max delta T"														<< endl;
	cout << "   --t-max-ud    DOUBLE	max T"																<< endl;
	cout << "   --t-max-local DOUBLE	max correction coefficient (local)"									<< endl;
	
	cout << "   --tol-rel DOUBLE		starting relative tolerance (default 1e-2)"							<< endl;
	cout << "   --tol-abs DOUBLE		starting absolute tolerance (default 1.e-8)"						<< endl;
	cout << "   --max-tol-rel DOUBLE	max relative tolerance (default 1e-5)"								<< endl;
	cout << "   --max-tol-abs DOUBLE	max absolute tolerance (default 1.e-9)"								<< endl;
	cout << "   --max-newtons INT	    max number of newton iterations (default 10)"						<< endl;
	cout << "   --f1stop DOUBLE			single reactor stop condition (default 1e-6)"						<< endl;

	cout << "   --species-fluctuations STRING	fluctuations: name of file"																<< endl;

	cout << "   -only-ode        		only ODE system"									<< endl;
	cout << "   -jacobian        		analytical jacobian (if available)"									<< endl;
    cout << "   -backup					restart from backup"												<< endl;
    cout << "   -test					test"																<< endl;
    
	cout << "   -help                   this help"															<< endl;
    cout << "   -noverbose              no file nor video output"											<< endl;
    cout << endl;
}

int main(int argc, char* argv[])
{	
	RemoveWarningWindow();

    string				argument;
    string				folderCase;
    string				clusteringExe;
    string				kppExe;

	int					iAnalyticalJacobian;
	int					iKindOfCorrection;
	double				cicloCluster;
	int					relaxation;
	int					fromBackUp;
	bool				iSequence;
	bool				iSpeciesFluctuations;
	string				speciesFluctuations;
	string				sequenceFile;
	string				kineticSchemePath;
	double				max_correction;
	double				t_max_delta;
	double				t_max_ud;
	double				t_max_local;

	double				max_rel_tol;
	double				max_abs_tol;
	double				rel_tol;
	double				abs_tol;
	int					max_newtons;
	double				f1stop;
	bool				only_ode;


	ParserClass				parser;

	#if LINUX_SO==1
		string copystring    = "cp ";
		string removefiles   = "rm ";
		string removefolders = "rm -r ";
	#else
		string copystring = "copy ";
		string removefiles = "del ";
		string removefolders = "rmdir /S /Q";
	#endif

    // 0. Parser setup
    parser.setup(argc, argv, parser_options);

	// 1. Detailed kinetic scheme setup
    if (parser.parse("--kinetics", argument))
		kineticSchemePath = argument;
    else    ErrorMessage("The --kinetics option is compulsory");

	// 2. Clsustering Software
    if (parser.parse("--clustering-exe", argument))
		clusteringExe = argument;
    else    ErrorMessage("The --clustering-exe option is compulsory");

	// 3. KPP Software
    if (parser.parse("--kpp-exe", argument))
		kppExe = argument;
    else    ErrorMessage("The --kpp-exe option is compulsory");

	// 4. Automated sequence
    if (parser.parse("--sequence", argument))
		sequenceFile	= argument;
    else    ErrorMessage("The --sequence option is compulsory");

	// 5. Input Folder
    if (parser.parse("--input", argument))
		folderCase = argument;
    else    ErrorMessage("The --input option is compulsory");

	// 6. Correction
    if (parser.parse("--correction", argument))
	{
		if		(argument == "none")		iKindOfCorrection = 0;
		else if (argument == "sin")			iKindOfCorrection = 1;
		else if (argument == "dirac")		iKindOfCorrection = 2;
		else if (argument == "beta")		iKindOfCorrection = 3;
		else if (argument == "gaussian")	iKindOfCorrection = 4;
		else ErrorMessage("The --correction option is wrong");
	}
    else    ErrorMessage("The --correction option is compulsory");

	// 7. Recover from backup
    if (parser.parse("-backup"))
		fromBackUp = 1;
	else fromBackUp = 0;

	// 8. Species fluctuations
    if (parser.parse("--species-fluctuations", argument))
	{
		iSpeciesFluctuations	= true;
		speciesFluctuations		= argument;
	}
	else iSpeciesFluctuations = false;


	// Algorithm
	{

		int				dummy_int;
		string			command_line;

		BzzVector	    iClusterList;
		BzzVectorInt    iCyclesList;

		// Reading Sequence File
		{
			ifstream fSequence;
			openInputFileAndControl(fSequence, sequenceFile);
			for(;;)
			{
				int dummy_int;
				double dummy_double;
				string dummy_string;
				fSequence >> dummy_string;
				if	(dummy_string == "#RUN")
				{
					fSequence >> dummy_double;
					iClusterList.Append(dummy_double);
					fSequence >> dummy_int;
					iCyclesList.Append(dummy_int);
				}
				else if (dummy_string == "#END")		break;
				else ErrorMessage("Expected: #END || #RUN");
			}
			fSequence.close();
		}

		string originalFolderCase = folderCase + "\\Original"; 

		for(int i=1;i<=iClusterList.Size();i++)
		{
			stringstream number;
			number << iClusterList[i];

			{
				string command1 = "cd " + originalFolderCase;
				system(command1.c_str());

				string command2 = clusteringExe + " --input Original --kinetics " + kineticSchemePath + "  --clustering " + number.str() + " -only-clustering" +
												  " --correction none";
				cout << command2 << endl;
				system(command2.c_str());

				string command3 = "mkdir Clustering" + number.str();
				system(command3.c_str());

				string command5 = "mkdir Clustering" + number.str() + "\\Input";
				system(command5.c_str());

				string command51 = "mkdir Clustering" + number.str() + "\\Output";
				system(command51.c_str());

				string command52 = "mkdir Clustering" + number.str() + "\\Output\\Backup";
				system(command52.c_str());

				string command6 = copystring + " ClusteringTopology.out " + "Clustering" + number.str() + "\\Input\\ClusteringTopology.out";
				system(command6.c_str());

				string command7 = copystring + " Clustering.CFDNetwork.bzz " + "Clustering" + number.str() + "\\Input\\CFDNetwork.bzz";
				system(command7.c_str());

				string command8 = copystring + " Clustering.FirstGuess.bzz " + "Clustering" + number.str() + "\\Input\\FirstGuess.bzz";
				system(command8.c_str());

				string command9 = removefiles + " Clustering.FirstGuess.bzz Clustering.CFDNetwork.bzz ClusteringTopology.out";
				system(command9.c_str());

				string command10 = removefolders + " Temp";
				system(command10.c_str());

				string command11 = removefolders + " Output";
				system(command11.c_str());

			}

			{
				// Copy Backup Files
				if (i>1)
				{
					stringstream number_previous;
					number_previous << iClusterList[i-1];

					string command = copystring + "Clustering" + number_previous.str() + "\\Output\\Backup\\Backup.end " + "Clustering" + number.str() + "\\Output\\Backup\\Backup.start";
					system(command.c_str());

					string command1 = copystring + " Input.withbackup " + "Clustering" + number.str() + "\\Input.inp";
					system(command1.c_str());
				}
				else
				{
					string command1 = copystring + " Input.withoutbackup " + "Clustering" + number.str() + "\\Input.inp";
					system(command1.c_str());
				}

				string command2 = " cd Clustering" + number.str();
				system(command2.c_str());

				string command3 = kppExe;
				system(command3.c_str());

				string command4 = " cd .." + number.str();
				system(command4.c_str());
			}

		}

		getchar();

		return(0);
	}
}
/*
		// File Preparation
		ofstream fOutput;
		openOutputFileAndControl(fOutput, "LogFile.log");
		fOutput.setf(ios::scientific);
		fOutput << "-------------------------------------------------" << endl;
		fOutput << "                *** KinPP 2009 ***               " << endl;
		fOutput << "-------------------------------------------------" << endl;

		// File Preparation
		ofstream fOutlet;
		openOutputFileAndControl(fOutlet, "Outlet.out");
		fOutlet.setf(ios::scientific);
		
		
		// Starting point
		int count = 1;
		{
			OpenSMOKE_CSTRNetwork	cstr;
			cstr.AssignKineticScheme(mix);
			cstr.SetMemoTemperature();
			
			// Options
			{
				if (iSpeciesFluctuations == true)	cstr.SetFluctuationsList(speciesFluctuations);

				if (max_correction>=0)	cstr.SetMaxCorrectionCoefficient(max_correction);

				 if (t_max_delta>=0)			cstr.SetDeltaTFluctuationsMaxDelta(t_max_delta);
				else if (t_max_ud>=0)			cstr.SetDeltaTFluctuationsMaxUserDefined(t_max_ud);
				else if (t_max_local>=0)		cstr.SetDeltaTFluctuationsMaxLocal(t_max_local);

				if (max_rel_tol>0.)		cstr.SetMaxTolRel(max_rel_tol);
				if (max_abs_tol>0.)		cstr.SetMaxTolAbs(max_abs_tol);
				if (rel_tol>0.)			cstr.SetTolRel(rel_tol);
				if (abs_tol>0.)			cstr.SetTolAbs(abs_tol);
				if (max_newtons>0)		cstr.SetMaxCountNewtonIterations(max_newtons);
				if (f1stop>0.)			cstr.SetF1Stop(f1stop);
				if (only_ode==true)		cstr.SetOdeOnly(true);
			}

			time_t curr=time(0);	string current_time = ctime(&curr);
			fOutput << "#" << count << "\tClustering: " << iClusterList[1] << "\tCycle: " << "1/"<<  iCycleList[1] << " " << current_time << endl;		
			cstr(CFDNetworkFile, FirstGuessFile, TolerancesFile,  iClusterList[1], 1, fromBackUp, 0, iAnalyticalJacobian, analyticalJacobian, iKindOfCorrection);
			cstr.OutputPrint(0);
			cstr.Save('*', MassFile);
			cstr.Maps();

			cstr.GiveMeOutputLabel(fOutlet);
			cstr.GiveMeOutput(fOutlet, count);
			
			count++;
		}

		// Next points
		for(int k=1;k<=iClusterList.Size();k++)
		{
			for(int j=1;j<=iCycleList[k];j++)
			{	
				OpenSMOKE_CSTRNetwork	cstr;
				cstr.AssignKineticScheme(mix);
				cstr.SetMemoTemperature();
				
				// Options
				{
					if (iSpeciesFluctuations == true)	cstr.SetFluctuationsList(speciesFluctuations);

					if (max_correction>=0)		cstr.SetMaxCorrectionCoefficient(max_correction);

						 if (t_max_delta>=0)	cstr.SetDeltaTFluctuationsMaxDelta(t_max_delta);
					else if (t_max_ud>=0)		cstr.SetDeltaTFluctuationsMaxUserDefined(t_max_ud);
					else if (t_max_local>=0)	cstr.SetDeltaTFluctuationsMaxLocal(t_max_local);

					if (max_rel_tol>0.)		cstr.SetMaxTolRel(max_rel_tol);
					if (max_abs_tol>0.)		cstr.SetMaxTolAbs(max_abs_tol);
					if (rel_tol>0.)			cstr.SetTolRel(rel_tol);
					if (abs_tol>0.)			cstr.SetTolAbs(abs_tol);
					if (max_newtons>0)		cstr.SetMaxCountNewtonIterations(max_newtons);
					if (f1stop>0.)			cstr.SetF1Stop(f1stop);
					if (only_ode==true)		cstr.SetOdeOnly(true);
				}


				time_t curr=time(0);	string current_time = ctime(&curr);
				fOutput << "#" << count << "\tClustering: " << iClusterList[k] << "\tCycle: " << j << "/"<<  iCycleList[k]  << " " << current_time << endl;
				cstr(CFDNetworkFile, FirstGuessFile, TolerancesFile, iClusterList[k], 1, 1, 0, iAnalyticalJacobian, analyticalJacobian, iKindOfCorrection);
				cstr.OutputPrint(0);
				cstr.Save('*', MassFile);
				cstr.Maps();
				cstr.GiveMeOutput(fOutlet, count);
				count++;

				stringstream number; number << count;
				command_line = copystring + "FinalSummary.out Output\\FinalSummary_" + number.str() + ".out";
				system(command_line.c_str());
			}

			stringstream number; number << iClusterList[k];
			string folderName = "Output\\clustering_" + number.str();
			command_line = "mkdir " + folderName;
			system(command_line.c_str());
			#if LINUX_SO == 1
				command_line = copystring + "Output/maps/*.* " + folderName;
				system(command_line.c_str());
				command_line = copystring + "Temp/mass.tmp "   + folderName;
				system(command_line.c_str());
				command_line = copystring + "FinalSummary.out " + folderName;
				system(command_line.c_str());
			#else
				command_line = copystring + "Output\\maps\\*.* " + folderName;
				system(command_line.c_str());
				command_line = copystring + "Temp\\mass.tmp "   + folderName;
				system(command_line.c_str());
				command_line = copystring + "FinalSummary.out " + folderName;
				system(command_line.c_str());
			#endif
		}

		fOutput.close();
		fOutlet.close();
	}
	


















	max_correction = -1;
	if (parser.parse("--max-correction", argument))
		max_correction = atof(argument.c_str());


	t_max_delta = -1;
	if (parser.parse("--t-max-delta", argument))
		t_max_delta = atof(argument.c_str());

	t_max_ud = -1;
	if (parser.parse("--t-max-ud", argument))
		t_max_ud = atof(argument.c_str());

	t_max_local = -1;
	if (parser.parse("--t-max-local", argument))
		t_max_local = atof(argument.c_str());


	max_rel_tol = -1;
	if (parser.parse("--max-tol-rel", argument))
		max_rel_tol = atof(argument.c_str());

	max_abs_tol = -1;
	if (parser.parse("--max-tol-abs", argument))
		max_abs_tol = atof(argument.c_str());

	rel_tol = -1;
	if (parser.parse("--tol-rel", argument))
		rel_tol = atof(argument.c_str());

	abs_tol = -1;
	if (parser.parse("--tol-abs", argument))
		abs_tol = atof(argument.c_str());

	max_newtons = -1;
	if (parser.parse("--max-newtons", argument))
		max_newtons = atoi(argument.c_str());

	f1stop = -1;
	if (parser.parse("--f1stop", argument))
		f1stop = atof(argument.c_str());

	only_ode = false;
	if (parser.parse("-only-ode"))
		only_ode = true;


	// -------------------------------------------------------------- //
	//					Environment preparation						  //
	// -------------------------------------------------------------- //
	string MassFile        = "Temp/mass.tmp";
	string CFDNetworkFile  = folderCase + "/CFDNetwork.bzz";
	string FirstGuessFile  = folderCase + "/FirstGuess.bzz";
	string TolerancesFile  = folderCase + "/Tolerances.bzz";
	
	string mkdir = "mkdir ";
	string buildOutputFile		= "mkdir Output";
	string buildTempFile		= "mkdir Temp";
	string buildEnthalpyFile	= mkdir + "\"Output/enthalpy\"";
	string buildMassFile		= mkdir + "\"Output/mass\"";
	string buildMapsFile		= mkdir + "\"Output/maps\"";
	string buildResidualsFile	= mkdir + "\"Output/residuals\"";
	system(buildOutputFile.c_str());
	system(buildTempFile.c_str());
	system(buildMassFile.c_str());
	system(buildMapsFile.c_str());
	system(buildEnthalpyFile.c_str());
	system(buildResidualsFile.c_str());

	// -------------------------------------------------------------- //
	//					Network solution    						  //
	// -------------------------------------------------------------- //

	if (iSequence == false)
	{	
		OpenSMOKE_CSTRNetwork	cstr;
		cstr.AssignKineticScheme(mix);
		cstr.SetMemoTemperature();
		
		// Options
		{
			if (iSpeciesFluctuations == true)	cstr.SetFluctuationsList(speciesFluctuations);

			if (max_correction>=0)				cstr.SetMaxCorrectionCoefficient(max_correction);

				 if (t_max_delta>=0)		cstr.SetDeltaTFluctuationsMaxDelta(t_max_delta);
			else if (t_max_ud>=0)			cstr.SetDeltaTFluctuationsMaxUserDefined(t_max_ud);
			else if (t_max_local>=0)		cstr.SetDeltaTFluctuationsMaxLocal(t_max_local);
			
			if (max_rel_tol>0.)		cstr.SetMaxTolRel(max_rel_tol);
			if (max_abs_tol>0.)		cstr.SetMaxTolAbs(max_abs_tol);
			if (rel_tol>0.)			cstr.SetTolRel(rel_tol);
			if (abs_tol>0.)			cstr.SetTolAbs(abs_tol);
			if (max_newtons>0)		cstr.SetMaxCountNewtonIterations(max_newtons);
			if (f1stop>0.)			cstr.SetF1Stop(f1stop);
			if (only_ode==true)		cstr.SetOdeOnly(true);
		}

		// No Relaxation
		if(relaxation == 0)
		{
			int cicloDiffusion = 1;
			
			startUser	= BzzGetUserTime();
			startKernel = BzzGetKernelTime();
			start		= BzzGetCpuTime();

			cstr(CFDNetworkFile, FirstGuessFile, TolerancesFile, cicloCluster, cicloDiffusion, fromBackUp, relaxation, iAnalyticalJacobian, analyticalJacobian, iKindOfCorrection);
			cstr.OutputPrint(0);
			cstr.Save('*', MassFile);
		}

		// Start from previous solution
		else if(relaxation == 2)
		{
			int cicloDiffusion	= 1;
			fromBackUp			= 1;

			startUser	= BzzGetUserTime();
			startKernel = BzzGetKernelTime();
			start		= BzzGetCpuTime();

			cstr(CFDNetworkFile, FirstGuessFile, TolerancesFile, cicloCluster, cicloDiffusion, fromBackUp, relaxation, iAnalyticalJacobian, analyticalJacobian, iKindOfCorrection);
			cstr.OutputPrint(0);
			cstr.Save('*', MassFile);
		}

		// Relaxation
		else if(relaxation == 1)
		{
			int cicloDiffusion = 4;

			startUser	= BzzGetUserTime();
			startKernel = BzzGetKernelTime();
			start		= BzzGetCpuTime();

			cstr(CFDNetworkFile, FirstGuessFile, TolerancesFile, cicloCluster, cicloDiffusion, fromBackUp, relaxation, iAnalyticalJacobian, analyticalJacobian, iKindOfCorrection);
			cstr.OutputPrint(0);
			cstr.Save('*',MassFile);
			fromBackUp = 1;
			for(cicloDiffusion = 3;cicloDiffusion > 0;cicloDiffusion--)
			{
				printf("\ncicloDiffusion %d",cicloDiffusion);
				cstr(CFDNetworkFile,FirstGuessFile, TolerancesFile, cicloCluster, cicloDiffusion, fromBackUp, relaxation, iAnalyticalJacobian, analyticalJacobian, iKindOfCorrection);
				cstr.OutputPrint(0);
				cstr.Save('*',MassFile);
			}
		}

		cstr.Maps();

		cout << "Cpu Seconds for complete solution:    " << BzzGetCpuTime() - start				<< endl;
		cout << "User Seconds for complete solution:   " << BzzGetUserTime() - startUser		<< endl;
		cout << "Kernel Seconds for complete solution: " << BzzGetKernelTime() - startKernel << endl;
	}

	if (iSequence == true)
	{
	return 0;
}
/*
#if SYMBOLIC_KINETICS == 1
void test(const string fileName)
{
	cout.setf(ios::scientific);

	const double	R_CTOT = 0.0820578337034;
	const double	UR_CTOT = 1./R_CTOT;
	
	const int		numCSTRReactors = 1;
	const double	Temperature = 1500.;	// Temperature [K]
	const double	Pressure	= 1.;		// Pressure [atm]

	SymbolicKinetics analyticalJacobian;
	OpenSMOKE_ReactingGas gas;
	gas.SetupBinary(fileName);

	
	OpenSMOKE_SymbolicKinetics* reactor[numCSTRReactors+1];

	string kineticSchemeName = GiveMeFileNameFromFullPath(fileName);
	if		(kineticSchemeName == "PolimiC1C3HTNOX_0810")	analyticalJacobian = POLIMI_C1C3HTNOX_0810;
	else if (kineticSchemeName == "PolimiC1C3HTNOX_AVIO")	analyticalJacobian = POLIMI_C1C3HTNOX_AVIO;
	else if (kineticSchemeName == "Fluent_Glarborg152")		analyticalJacobian = FLUENT_GLARBORG152;
	else if (kineticSchemeName == "GRI30")					analyticalJacobian = GRI30;
	else if (kineticSchemeName == "SanDiego_AVIO")			analyticalJacobian = SANDIEGO_AVIO;
	else if (kineticSchemeName == "Fluent_DRM22_Polimi")			analyticalJacobian = FLUENT_DRM22_POLIMI;
	else if (kineticSchemeName == "Fluent_DRM22_Polimi_NOX")		analyticalJacobian = FLUENT_DRM22_POLIMI_NOX;
	else if (kineticSchemeName == "Fluent_DRM22_Polimi_ThermalNOX")	analyticalJacobian = FLUENT_DRM22_POLIMI_THERMALNOX;
	else ErrorMessage("Symbolic Jacobian not available for this kinetic scheme");

	if (analyticalJacobian == POLIMI_C1C3HTNOX_0810)
		for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
			reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_0810();

	if (analyticalJacobian == POLIMI_C1C3HTNOX_AVIO)
		for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
			reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Polimi_C1C3HTNOX_AVIO();

	if (analyticalJacobian == FLUENT_GLARBORG152)
		for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
			reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Fluent_Glarborg152();

	if (analyticalJacobian == GRI30)
		for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
			reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_GRI30();

	if (analyticalJacobian == SANDIEGO_AVIO)
		for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
			reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_SanDiego_AVIO();

	if (analyticalJacobian == FLUENT_DRM22_POLIMI)
		for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
			reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi();

	if (analyticalJacobian == FLUENT_DRM22_POLIMI_NOX)
		for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
			reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_NOX();

	if (analyticalJacobian == FLUENT_DRM22_POLIMI_THERMALNOX)
		for (int kCSTR=1;kCSTR<=numCSTRReactors;kCSTR++)
			reactor[kCSTR] = new OpenSMOKE_SymbolicKinetics_Fluent_DRM22_Polimi_ThermalNOX();

	if (reactor[1]->NC != gas.NumberOfSpecies())	ErrorMessage("The number of species in the Symbolic Kinetics is wrong");
	if (reactor[1]->NR != gas.NumberOfReactions())	ErrorMessage("The number of reactions in the Symbolic Kinetics is wrong");


	BzzVector T(numCSTRReactors);
	BzzVector cTot(numCSTRReactors);
	BzzVector P_Pa(numCSTRReactors);

	T		= Temperature;
	cTot	= Pressure / Temperature * UR_CTOT;
	P_Pa	= Pressure*101325.;

	gas.InitializeMap(numCSTRReactors);
	gas.ComputeKineticParameters_map(T, cTot, P_Pa);

	BzzVector k1_map        = gas.k1_map.GetRow(1);
	BzzVector uKeq_map      = gas.uKeq_map.GetRow(1);
	BzzVector logFCent_map  = gas.logFcent_map.GetRow(1);
	BzzVector k2_map        = gas.k2_map.GetRow(1);
	

	reactor[1]->assignKineticConstants(k1_map, uKeq_map, logFCent_map, k2_map);
	
	BzzVector omega(gas.NumberOfSpecies());
	BzzVector x(gas.NumberOfSpecies());
	BzzVector c(gas.NumberOfSpecies());
	BzzVector cOld(gas.NumberOfSpecies());
	BzzVector R1(gas.NumberOfSpecies());
	BzzVector R2(gas.NumberOfSpecies());
	BzzMatrix dRC1(gas.NumberOfSpecies(), gas.NumberOfSpecies());
	BzzMatrix dRC2(gas.NumberOfSpecies(), gas.NumberOfSpecies());
	BzzVector r1, r2;


	// Calculate Concentrations
	// -------------------------------------------------------------------
	omega = double(1./gas.NumberOfSpecies());
	double wM = gas.GetMWFromMassFractions(omega);
	gas.GetMoleFractionsFromMassFractionsAndMW(x,omega,wM);
	c	 = cTot[1]*x;
	cOld = c;

	// -------------------------------------------------------------------
	// -------------------------------------------------------------------
	//						Calculate Reaction Rates
	// -------------------------------------------------------------------
	// -------------------------------------------------------------------
	int nCycles = 50000;

	// Calculate Reaction Rates - SemiAnalytical
	// -------------------------------------------------------------------
	int kk;
	double start = BzzGetCpuTime();
	for (kk=1;kk<=nCycles;kk++)
		gas.ComputeFromConcentrations_map(1, c, &R1);
	double tEnd = BzzGetCpuTime() - start;
	r1 = gas.r;
	cout << "Reaction Rates (Semi-Analytical) - Time: " << tEnd << endl;

	// Calculate Reaction Rates - Analytical
	// -------------------------------------------------------------------
	start = BzzGetCpuTime();
	for (kk=1;kk<=nCycles;kk++)
		reactor[1]->giveReactionRates(cTot[1], c, R2);
	tEnd = BzzGetCpuTime() - start;
	r2 = gas.r;
	cout << "Reaction Rates (Analytical) - Time: " << tEnd << endl;
	cout << "Press enter..." << endl;
	getchar();
	
	int i;
	cout << endl << endl;
	for(i=1; i<=gas.NumberOfReactions();i++)
	{
		double maximum = Max(fabs(r1[i]), fabs(r2[i]));
		double rel_err = 0.;
		if (maximum>0.)	rel_err = fabs(r1[i]-r2[i])/maximum;
		cout << i << "\t" << r1[i] << "\t" << r2[i] << "\t" << rel_err << endl;
		if (rel_err > 1.e-8)
			cout << " ATTENTION: Reaction " << i << " - Relative Error " << rel_err << endl;
	}
	cout << endl << endl;

	// Write Formation rates
	// -------------------------------------------------------------------
	for(i=1; i<=gas.NumberOfSpecies();i++)
	{
		double maximum = Max(fabs(R1[i]), fabs(R2[i]));
		double rel_err = 0.;
		if (maximum>0.)	rel_err = fabs(R1[i]-R2[i])/maximum;
		cout << i << "\t" << R1[i] << "\t" << R2[i] << "\t" << rel_err << endl;
		if (rel_err > 1.e-8)
			cout << " ATTENTION: Species " << i << " - Relative Error " << rel_err << endl;
	}
	cout << endl << endl;



	// -------------------------------------------------------------------
	// -------------------------------------------------------------------
	//						Calculate Jacobian
	// -------------------------------------------------------------------
	// -------------------------------------------------------------------

	int nCyclesJ = 5000;

	// Calculate Jacobian - SemiAnalytical
	// -------------------------------------------------------------------
	start = BzzGetCpuTime();
	for (kk=1;kk<=nCyclesJ;kk++)
	{
		// Reaction Rates
		gas.ComputeFromConcentrations_map(1, c, &R1);

		// Calculation of the Jacobian
		gas.kinetics.mc		= c;
		gas.kinetics.mR		= R1;
		gas.kinetics.mr		= gas.r;
		gas.kinetics.mrDirT	= gas.coeffFallOff;
		gas.kinetics.mrDirC	= gas.rDirC;
		gas.kinetics.mrInvC	= gas.rInvC;

		// 1. Contributi dai termini di reazione
		gas.kinetics.GetDerivativesC(T[1], cTot[1], &dRC1, c, R1);
		c = cOld;
	}
	tEnd = BzzGetCpuTime() - start;
	cout << "Time: " << tEnd << endl;

	// Calculate Jacobian - Analytical
	// -------------------------------------------------------------------
	start = BzzGetCpuTime();
	for (kk=1;kk<=nCyclesJ;kk++)
	{
		reactor[1]->giveReactionRates(cTot[1], c, R2);
		reactor[1]->giveJacobian(c, dRC2);
		c = cOld;
	}
	tEnd = BzzGetCpuTime() - start;
	cout << "Time: " << tEnd << endl;


	// Print Jacobian On File
	// -------------------------------------------------------------------
	ofstream fOutput("Jacobian.out", ios::out);
	fOutput.setf(ios::scientific);

	for(i=1; i<=gas.NumberOfSpecies();i++)
		for(int j=1; j<=gas.NumberOfSpecies();j++)
			fOutput << i << " " << j << "\t\t"
			<< dRC1[i][j] << "\t"
			<< dRC2[i][j] << "\t"
			<< dRC1[i][j]-dRC2[i][j] << "\t"
			<< (dRC1[i][j]-dRC2[i][j])/(dRC2[i][j]+1.e-10) << endl;
	fOutput.close();
}
#else*/
void test(const string fileName) {}
//#endif
