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

#include <sstream>
#include "interfaces/SimpleOpt.h"
#include "CRECK_Optimizer_Class.h"

void ErrorMessage(const string message)
{
    cout << endl;
    cout << "Executable:  OpenSMOKE_Optimizer"	<< endl;
    cout << "Error:       "  << message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

// define the ID values to indentify the option
enum { OPT_HELP, OPT_FLAG, OPT_ARG };

CSimpleOpt::SOption parser_options[] =
{
    { OPT_HELP, _T("-?"),              SO_NONE    }, // "-?"
    { OPT_HELP, _T("-help"),           SO_NONE    }, // "-help"
	{ OPT_ARG,  _T("-regression"),     SO_NONE    }, // "-regression"
	{ OPT_ARG,  _T("--global"),        SO_REQ_SEP }, // "--global       ARG"
	{ OPT_ARG,  _T("--kinetics"),      SO_REQ_SEP }, // "--kinetics     ARG"
	{ OPT_ARG,  _T("--one_per_time"),  SO_REQ_SEP }, // "--one_per_time ARG"
	{ OPT_ARG,  _T("--map2d"),         SO_REQ_SEP }, // "--map2d        ARG"
	{ OPT_ARG,  _T("--model"),         SO_REQ_SEP }, // "--model        ARG"
	{ OPT_ARG,  _T("--online-pp"),     SO_REQ_SEP }, // "--online-pp    ARG"
    SO_END_OF_OPTIONS								 // END
};

CRECK_Optimizer_Class *ptOptimizer;

double Minimization(BzzVector &t);
void Regression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);


int main(int argc, char* argv[])
{
	RemoveWarningWindow();

    OpenSMOKE_logo("OpenSMOKE_Optimizer", "0.3", "November 2009");

    // ---------------------------------------------------------------------
    // Data declaration
    // ---------------------------------------------------------------------
    int nParameters;
    int nTotalParameters;
	string kineticsFolder;
	string globalName;
    BzzVector startingParameters;
    BzzVector minParameters;
    BzzVector maxParameters;
    BzzVector xSolution;
    BzzVector ySolution;

	string          argument;
    ParserClass     parser;
	bool			one_per_time;
	bool			iOptimization;
	bool			iMap2D;
	int				nTotalPoints;
	int				iModel;
	int				online_pp = 0;

	//	0. Parser setup
	parser.setup(argc, argv, parser_options);

	system("mkdir Output_Optimization");

    // ---------------------------------------------------------------------
    // Optimizer Setup
    // ---------------------------------------------------------------------
    if (parser.parse("--one_per_time", argument))
    {
		one_per_time = true;
		nTotalParameters = atoi(argument.c_str());
	}
	else
		one_per_time = false;

	if (parser.parse("--model", argument))
    {
		iModel = atoi(argument.c_str());
		if (iModel != 1 && iModel != 4 && iModel != 2)
			ErrorMessage("The --model option must be 1 || 2 || 4");
	}
	else
		ErrorMessage("The --model option is compulsory");

	if (parser.parse("--map2d", argument))
    {
		iMap2D = true;
		nTotalPoints = atoi(argument.c_str());
		iModel = 1;
	}
	else
		iMap2D = false;

	if (parser.parse("--online-pp", argument))
		online_pp = atoi(argument.c_str());

	if (parser.parse("-regression"))	iOptimization = false;
	else								iOptimization = true;

    if (parser.parse("--kinetics", argument))	kineticsFolder = argument.c_str();
	else										ErrorMessage("The --kinetics option is compulsory");

	if (parser.parse("--global", argument))		globalName = argument.c_str();
	else										ErrorMessage("The --global option is compulsory");


	// Total optimization
	if (one_per_time == false)
	{
		BzzMinimizationRobust minimization;
		CRECK_Optimizer_Class optimizer;
		ptOptimizer = &optimizer;

		// Options
		optimizer.SetModel(iModel);
		if (online_pp != 0 )
		{
			optimizer.SetOnlinePostProcessing();
			optimizer.SetCountOnlinePostProcessing(online_pp);
		}

		// Folder and kinetic scheme
		optimizer.setup(kineticsFolder, globalName, 0);

		// ---------------------------------------------------------------------
		// Memory allocation
		// ---------------------------------------------------------------------
		nParameters = optimizer.GiveMeNumberOfParameters();
		ChangeDimensions(nParameters, &xSolution);
		ChangeDimensions(nParameters, &ySolution);


		// ---------------------------------------------------------------------
		// BzzMinimizationDoubleRobust preparation
		// ---------------------------------------------------------------------
		startingParameters  = optimizer.GiveMeStartingParameters();
		minParameters		= optimizer.GiveMeMinParameters();
		maxParameters		= optimizer.GiveMeMaxParameters();

		double startTime;
		double endTime;

		// ---------------------------------------------------------------------
		// Solution of problem
		// ---------------------------------------------------------------------
		if (iMap2D == true)
		{
			minimization(startingParameters, Minimization, &minParameters, &maxParameters);

			startTime    = BzzGetCpuTime();

			ofstream fMap;
			openOutputFileAndControl(fMap, "Map2D.map");
			fMap.setf(ios::scientific);

			BzzVector parameters(2);
			for (int i=1;i<=nTotalPoints;i++)
			{
				parameters[1] = minParameters[1] + (maxParameters[1]-minParameters[1])/double(nTotalPoints-1.)*double(i-1.);
				for (int j=1;j<=nTotalPoints;j++)
				{
					parameters[2] = minParameters[2] + (maxParameters[2]-minParameters[2])/double(nTotalPoints-1.)*double(j-1.);
					
					cout << parameters[1] << "\t" << parameters[2] << endl;
					double f = optimizer.Minimization(parameters);
					
					fMap << parameters[1] << "\t" << parameters[2] << "\t" << f << endl;
				}
				fMap << endl;
			}
			fMap.close();
			endTime      = BzzGetCpuTime();
		}
		else
		{
			if (iOptimization == true)
			{
				if (iModel == 1 || iModel == 2)
					minimization(startingParameters, Minimization, &minParameters, &maxParameters);
				else if (iModel == 4)
					minimization(startingParameters, Minimization);

				startTime    = BzzGetCpuTime();
				minimization();
				endTime      = BzzGetCpuTime();
				minimization.GetXSolution(&xSolution);
			}
			else
			{
				startTime    = BzzGetCpuTime();
				
				BzzVector x_experimental;
				BzzVector y_experimental;

				optimizer.GiveMeRegressionExperimentalData(x_experimental, y_experimental);

				int numModels		= 1;
				int numX			= 1;
				int numY			= 1;
				int numExperiments	= x_experimental.Size();
			

				BzzMatrix X(numExperiments,numX);
				BzzMatrix Y(numExperiments,numY);
				X.SetColumn(1, x_experimental);
				Y.SetColumn(1, y_experimental);

				BzzNonLinearRegression nonLinReg(numModels, X,Y, Regression);
				BzzVector s2(numY,Variance(y_experimental));
				nonLinReg.SavePartial("RegressionBest.out");
				nonLinReg.SaveComplete("RegressionComplete.out");
				nonLinReg.SetVariance(numExperiments,s2);
				nonLinReg.InitializeModel(1, startingParameters, minParameters, maxParameters);

				nonLinReg.LeastSquaresAnalysis();
				nonLinReg.GetSolution(xSolution);

				endTime      = BzzGetCpuTime();
			}
		}


		// ---------------------------------------------------------------------
		// Finalizing
		// ---------------------------------------------------------------------
		cout << endl;
		for (int i=1;i<=nParameters;i++)
			cout << "Parameter " << i << ": " << xSolution[i]   << endl << endl;
		cout << "Total CPU time for optimization: " << endTime - startTime << " s" << endl;
		cout << "Press enter to continue... " << endl;

		// ---------------------------------------------------------------------
		// Post processing
		// ---------------------------------------------------------------------
		optimizer.PostProcessing("Start", startingParameters);
		optimizer.PostProcessing("Final", xSolution);
		
		optimizer.Gnuplot_plots();
		optimizer.Latex_report();
		system("dvips report.dvi");
		system("ps2pdf report.ps");
	}

	// Optimization one per time
	if (one_per_time == true)
	{
		for(int index=1;index<=nTotalParameters;index++)
		{
			BzzMinimizationRobust minimization;
			CRECK_Optimizer_Class optimizer;
			ptOptimizer = &optimizer;
			optimizer.SetModel(iModel);
			optimizer.setup(kineticsFolder, globalName, index);

			// ---------------------------------------------------------------------
			// Memory allocation
			// ---------------------------------------------------------------------
			nParameters = optimizer.GiveMeNumberOfParameters();
			ChangeDimensions(nParameters, &xSolution);
			ChangeDimensions(nParameters, &ySolution);


			// ---------------------------------------------------------------------
			// BzzMinimizationDoubleRobust preparation
			// ---------------------------------------------------------------------
			startingParameters  = optimizer.GiveMeStartingParameters();
			if (iModel == 1 || iModel == 2)
				minimization(startingParameters, Minimization, &minParameters, &maxParameters);
			else if (iModel == 4)
				minimization(startingParameters, Minimization);

			// ---------------------------------------------------------------------
			// Solution of problem
			// ---------------------------------------------------------------------
			double startTime    = BzzGetCpuTime();
			minimization();
			double endTime      = BzzGetCpuTime();
			minimization.GetXSolution(&xSolution);


			// ---------------------------------------------------------------------
			// Finalizing
			// ---------------------------------------------------------------------
			cout << endl;
			for (int i=1;i<=nParameters;i++)
				cout << "Parameter " << i << ": " << xSolution[i]   << endl << endl;
			cout << "Total CPU time for optimization: " << endTime - startTime << " s" << endl;
			cout << "Press enter to continue... " << endl;

			// ---------------------------------------------------------------------
			// Post processing
			// ---------------------------------------------------------------------
			optimizer.PostProcessing("Start", startingParameters);
			optimizer.PostProcessing("Final", xSolution);

			optimizer.Gnuplot_plots_one_per_time();
			optimizer.Latex_report_one_per_time();

			system("dvips report.dvi");
			system("ps2pdf report.ps");

			stringstream index_string;
			index_string << index; 
			string command = "rename report.pdf report_" + index_string.str() + ".pdf";
			system(command.c_str());
		}
	}


    OpenSMOKE_logo("OpenSMOKE_Optimizer", "0.2", "15th June 2009");

    cout << "Press enter to continue... " << endl;
	getchar();
    
	return EXIT_SUCCESS;
} 

void ParserClass::ShowUsage()
{
    cout << "Usage: OpenSMOKE_Optimizer [-option] [--option FILE] [-?] [--help]" << endl;

    cout << "Options: "                                                                          << endl;
    cout << "   -?                     this help"                                                << endl;
    cout << "   --global NAME          global kinetic scheme"                                    << endl;
    cout << "   --kinetics FOLDER      folder (kinetic scheme)"                                  << endl;
    cout << "   --one_per_time NUMBER  monodimensional optimization"                             << endl;
    cout << "   --map2d NUMBER         map 2d"                                                   << endl;
    cout << "   --model NUMBER         map 2d"                                                   << endl;
    cout << "   -help                  this help"                                                << endl;
    cout << "   -regression            regression instead optimization"                          << endl;
    cout << endl;
}




double Minimization(BzzVector &t)
{
    return ptOptimizer->Minimization(t);
}

void Regression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
    ptOptimizer->Regression(model, ex, b, x, y);
}