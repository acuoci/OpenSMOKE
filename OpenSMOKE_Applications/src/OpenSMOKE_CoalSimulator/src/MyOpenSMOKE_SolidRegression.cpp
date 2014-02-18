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

#include <sstream>
#include <iomanip>
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "basic/OpenSMOKE_Conversions.h"

#include "MyOpenSMOKE_CharKineticScheme.h"
#include "MyOpenSMOKE_SolidRegression.h"
#include "OpenSMOKE_SolidExperiment.h"
#include "OpenSMOKE_SolidExperiment_DSmoke_Bio.h"

MyOpenSMOKE_SolidRegression* ptRegression;
void MyBzzModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);

const double MyOpenSMOKE_SolidRegression::MWchar = 12.;	// [kg/kmol]

void ODE_Print(BzzVector &x, double t)
{
	if      (ptRegression->solid_regression == OPENSMOKE_REGRESSION_CHAR)		ptRegression->MyODE_Print_Char(x, t);
	else if (ptRegression->solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)	ptRegression->MyODE_Print_DSmoke_Bio(x, t);
}

void MyOpenSMOKE_SolidRegression::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SolidRegression"		<< endl;
    cout << "Object: " << name_object				<< endl;
    cout << "Error:  " << message					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void MyOpenSMOKE_SolidRegression::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SolidRegression"	<< endl;
    cout << "Object:  " << name_object				<< endl;
    cout << "Warning: " << message					<< endl;
	cout << endl;
}

MyOpenSMOKE_SolidRegression::MyOpenSMOKE_SolidRegression()
{
	name_object = "[Name not assigned]";
	count_global = 1;
	iWriteOnFile = false;

	// Default values
	// -----------------------------------------------------------
	iModel				= 1;						// 
	parameters			= SOLID_REGRESSION_ALL;		// 
	iSurfaceEquation	= true;
	iVerbose			= false;
	iRelativeTolerance	= false;
	iAbsoluteTolerance	= false;
	iRobustAnalysis		= false;

	minimum_ratio_A		= 0.01;
	maximum_ratio_A		= 100.;
	minimum_ratio_E		= 0.50;
	maximum_ratio_E		= 1.5;
	minimum_delta_BETA	= -1.;
	maximum_delta_BETA	=  1.;
	minimum_ratio_PHI	= 0.5;
	maximum_ratio_PHI	= 2.;
}

void MyOpenSMOKE_SolidRegression::AssignKinetics(OpenSMOKE_CharKineticScheme *_kinetics)
{
	kinetics = _kinetics;
	solid_regression = OPENSMOKE_REGRESSION_CHAR;
}

void MyOpenSMOKE_SolidRegression::AssignKinetics(OpenSMOKE_ReactingGas *_kinetics)
{
	mix_chemkin = _kinetics;
	solid_regression = OPENSMOKE_REGRESSION_DSMOKE_BIO;
}

void MyOpenSMOKE_SolidRegression::AssignListOfFiles(const string _file_name_list_of_files)
{
	file_name_list_of_files = _file_name_list_of_files;

	if (solid_regression == OPENSMOKE_REGRESSION_CHAR)
	{
		string dummy;
		vector<string> list_provisional;
		ifstream fInput;
		openInputFileAndControl(fInput, file_name_list_of_files);

		for(;;)
		{
			fInput >> dummy;
			if (dummy != "#END")
				list_provisional.push_back(dummy);
			else
				break;
		}
		fInput.close();

		nCases = list_provisional.size();
		fInput.close();
		//cout << nCases << endl;
	//	getchar();
	}
}

void MyOpenSMOKE_SolidRegression::SetName(const string name)
{
	name_object = name;
}

void MyOpenSMOKE_SolidRegression::SetModel(const int _iModel)
{
	iModel = _iModel;
}

void MyOpenSMOKE_SolidRegression::SetRegressionAnalysis(const bool index)
{
	iRegressionAnalysis = index;
}

void MyOpenSMOKE_SolidRegression::SetSurfaceEquation(const bool index)
{
	iSurfaceEquation = index;
}

void MyOpenSMOKE_SolidRegression::SetRobustAnalysis(const bool index)
{
	iRobustAnalysis = index;
}

void MyOpenSMOKE_SolidRegression::SetMinimumRatio_A(const double _value)
{
	minimum_ratio_A = _value;
}

void MyOpenSMOKE_SolidRegression::SetMaximumRatio_A(const double _value)
{
	maximum_ratio_A = _value;
}

void MyOpenSMOKE_SolidRegression::SetMinimumRatio_E(const double _value)
{
	minimum_ratio_E = _value;
}

void MyOpenSMOKE_SolidRegression::SetMaximumRatio_E(const double _value)
{
	maximum_ratio_E = _value;
}

void MyOpenSMOKE_SolidRegression::SetMinimumDelta_BETA(const double _value)
{
	minimum_delta_BETA = _value;
}

void MyOpenSMOKE_SolidRegression::SetMaximumDelta_BETA(const double _value)
{
	maximum_delta_BETA = _value;
}

void MyOpenSMOKE_SolidRegression::SetMinimumDelta_PHI(const double _value)
{
	minimum_ratio_PHI = _value;
}

void MyOpenSMOKE_SolidRegression::SetMaximumDelta_PHI(const double _value)
{
	maximum_ratio_PHI = _value;
}

void MyOpenSMOKE_SolidRegression::SetVerbose(const bool value)
{
	iVerbose = value;
}

void MyOpenSMOKE_SolidRegression::SetRemoveWarning(const bool value)
{
	if (value == true)	RemoveWarningWindow();
}

void MyOpenSMOKE_SolidRegression::SetAbsoluteTolerance(const double _absoluteTolerance)
{
	iAbsoluteTolerance = true;
	absoluteTolerance = _absoluteTolerance;
}

void MyOpenSMOKE_SolidRegression::SetRelativeTolerance(const double _relativeTolerance)
{
	iRelativeTolerance = true;
	relativeTolerance = _relativeTolerance;
}

void MyOpenSMOKE_SolidRegression::SetParameters(const vector<string> names)
{
		 if (names[0] == "ALL")				parameters = SOLID_REGRESSION_ALL;
	else if (names[0] == "A")				parameters = SOLID_REGRESSION_A;
	else if (names[0] == "BETA")			parameters = SOLID_REGRESSION_BETA;
	else if (names[0] == "E")				parameters = SOLID_REGRESSION_E;
	else if (names[0] == "A_E")				parameters = SOLID_REGRESSION_A_E;
	else if (names[0] == "BETA_E")			parameters = SOLID_REGRESSION_BETA_E;
	else if (names[0] == "A_BETA")			parameters = SOLID_REGRESSION_A_BETA;
	else if (names[0] == "USERDEFINED")		parameters = SOLID_REGRESSION_USERDEFINED;
	else if (names[0] == "ALL_PHI")			parameters = SOLID_REGRESSION_ALL_PHI;
	else if (names[0] == "A_E_PHI")			parameters = SOLID_REGRESSION_A_E_PHI;
	else 
		ErrorMessage("Available options: ALL ||ALL_PHI || A || BETA || E || A_E || A_E_PHI || BETA_E || A_BETA || USERDEFINED");
	int nr;
	if (solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)
		nr = mix_chemkin->NumberOfReactions();
	if (solid_regression == OPENSMOKE_REGRESSION_CHAR)
		nr = kinetics->number_of_reactions;

	if (parameters == SOLID_REGRESSION_ALL)
	{
		int j;
		ChangeDimensions(3*nr, 2, &bUserDefined);
		for(j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 1;
		}
		for(j=1;j<=nr;j++)
		{
			bUserDefined[nr+j][1] = j;
			bUserDefined[nr+j][2] = 2;
		}
		for(j=1;j<=nr;j++)
		{
			bUserDefined[2*nr+j][1] = j;
			bUserDefined[2*nr+j][2] = 3;
		}
  }

    //Tiziano
	if (parameters == SOLID_REGRESSION_ALL_PHI)
	{  
		int j;
		ChangeDimensions(3*nr+nCases, 2, &bUserDefined);
		for(j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 1;
		}
		for(j=1;j<=nr;j++)
		{
			bUserDefined[nr+j][1] = j;
			bUserDefined[nr+j][2] = 2;
		}
		for(j=1;j<=nr;j++)
		{
			bUserDefined[2*nr+j][1] = j;
			bUserDefined[2*nr+j][2] = 3;
		}
        
		for(j=1;j<=nCases;j++)
		{    
			bUserDefined[3*nr+j][1] = j;
 			bUserDefined[3*nr+j][2] = 4;
		}
	}

	if (parameters == SOLID_REGRESSION_A)
	{
		ChangeDimensions(nr, 2, &bUserDefined);
		for(int j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 1;
		}

	}
	
	else if (parameters == SOLID_REGRESSION_BETA)
	{
		ChangeDimensions(nr, 2, &bUserDefined);
		for(int j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 2;
		}
	}

	else if (parameters == SOLID_REGRESSION_E)
	{
		ChangeDimensions(nr, 2, &bUserDefined);
		for(int j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 3;
		}
	}

	else if (parameters == SOLID_REGRESSION_A_E)
	{
		int j;
		ChangeDimensions(2*nr, 2, &bUserDefined);
		for(j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 1;
		}
		for(j=1;j<=nr;j++)
		{
			bUserDefined[nr+j][1] = j;
			bUserDefined[nr+j][2] = 3;
		}

	}

	else if (parameters == SOLID_REGRESSION_A_E_PHI)
	{
		int j;
		ChangeDimensions(2*nr+nCases, 2, &bUserDefined);
		for(j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 1;
		}
		for(j=1;j<=nr;j++)
		{
			bUserDefined[nr+j][1] = j;
			bUserDefined[nr+j][2] = 3;
		}
		for(j=1;j<=nCases;j++)
		{
			bUserDefined[2*nr+j][1] = j; 
 			bUserDefined[2*nr+j][2] = 4;
		}
	}

	else if (parameters == SOLID_REGRESSION_A_BETA)
	{
		int j;
		ChangeDimensions(2*nr, 2, &bUserDefined);
		for(j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 1;
		}
		for(j=1;j<=nr;j++)
		{
			bUserDefined[nr+j][1] = j;
			bUserDefined[nr+j][2] = 2;
		}
	}

	else if (parameters == SOLID_REGRESSION_BETA_E)
	{
		int j;
		ChangeDimensions(2*nr, 2, &bUserDefined);
		for(j=1;j<=nr;j++)
		{
			bUserDefined[j][1] = j;
			bUserDefined[j][2] = 2;
		}
		for(j=1;j<=nr;j++)
		{
			bUserDefined[nr+j][1] = j;
			bUserDefined[nr+j][2] = 3;
		}
	}

	else if (parameters == SOLID_REGRESSION_USERDEFINED)
	{
		int i,j;
		int count = 0;
		//cout << names.size()-1 << endl;
		for(unsigned int ii=1;ii<=names.size()-1;ii++)
		{	
			cout << names[ii] << endl;
			if (names[ii] == "REACTION")
				count++;
			else if (names[ii] == "PHI")
				count+=nCases;
		}
		cout << count << endl;
		ChangeDimensions(count, 2, &bUserDefined);

		j=1;	
		i=1; 
		for(;;)
		{
			if (names[i] != "REACTION" && names[i] != "PHI")	
				ErrorMessage("Expected REACTION key word - Found: " + names[i]);
			cout << i << " " << names[i] << endl;

			if (names[i] == "REACTION")
			{
				bUserDefined[j][1] = atoi(names[i+1].c_str());
				if (bUserDefined[j][1] <=0 || bUserDefined[j][1]>nr)
					ErrorMessage("Please check the list of regression parameters: REACTION is out of boundaries...");
				
				if(names[i+2] == "A")			bUserDefined[j][2] = 1;
				else if(names[i+2] == "BETA")	bUserDefined[j][2] = 2;
				else if(names[i+2] == "E")		bUserDefined[j][2] = 3;

				else ErrorMessage("Expected A || BETA || E || PHI - Found: " + names[i+2]);
				
				j++;
				i+=3;
				if (i>int(names.size()-1))	break;
			}
			else if (names[i] == "PHI")
			{
				for(int k=1;k<=nCases;k++)
				{
					bUserDefined[j+k-1][1] = k;
					bUserDefined[j+k-1][2] = 4;
				}	

				j+=nCases;
				i+=1;
				if (i>int(names.size()-1))	break;
			}
		}	

		for(j=1;j<=bUserDefined.Rows();j++)
			for(int i=j+1;i<=bUserDefined.Rows();i++)
				if (bUserDefined[j][1] == bUserDefined[i][1])
					if (bUserDefined[j][2] == bUserDefined[i][2])
						ErrorMessage("The same regression parameter is specified more than once...");
	}

}

void MyOpenSMOKE_SolidRegression::SetSolidSpecies(const vector<string> _names)
{
	GiveMeIndicesAndNames(*mix_chemkin, _names, indices_solid, names_solid);
}

void MyOpenSMOKE_SolidRegression::SetSolidDensity(const string units, const double value)
{
    rho_solid = OpenSMOKE_Conversions::conversion_density(value, units);
}

void MyOpenSMOKE_SolidRegression::Setup(const string input_file_name)
{
	ptRegression = this;
	ode.assign(this);
	OpenSMOKE_Dictionary_SolidRegression dictionary;
	dictionary.ParseFile(input_file_name);
	CheckDictionary(dictionary);
	MemoryAllocation();

	if (solid_regression == OPENSMOKE_REGRESSION_CHAR)
	{
		int i;
		string dummy;
		ifstream fInput;
		openInputFileAndControl(fInput, file_name_list_of_files);

		for(;;)
		{
			fInput >> dummy;
			if (dummy != "#END")
				list_of_names_of_files.push_back(dummy);
			else
				break;
		}
		fInput.close();

		nCases = list_of_names_of_files.size();
		ChangeDimensions(nCases,											&indices);
		ChangeDimensions(nCases,											&initialconditions_x);
		ChangeDimensions(nCases,2+kinetics->nCarbon+kinetics->nAdsorbed,	&initialconditions_y);
		ChangeDimensions(2+kinetics->nCarbon+kinetics->nAdsorbed,			&yMin);
		ChangeDimensions(2+kinetics->nCarbon+kinetics->nAdsorbed,			&yMax);

		// Set each experiment
		experiments_char = new OpenSMOKE_SolidExperiment[nCases+1];
		for(i=1;i<=nCases;i++)
		{
			experiments_char[i].ReadFromFile(list_of_names_of_files[i-1]);
			experiments_char[i].PrepareMasses(*kinetics);
		}

		// Assembling experiments
		numTotal = 0;
		for(i=1;i<=nCases;i++)
		{
			indices[i]				= numTotal+1;

			initialconditions_x[i]	= experiments_char[i].x[1]; 
			
			int k=1;
			initialconditions_y[i][k++]	= 0.; 
			initialconditions_y[i][k++]	= 1.0;
			for(int j=1;j<=kinetics->nCarbon;j++)
				initialconditions_y[i][k++]	= experiments_char[i].initial_carbon_mass_fractions[j];
			for(int j=1;j<=kinetics->nAdsorbed;j++)
				initialconditions_y[i][k++]	= 0.0; 
			
			numTotal  += experiments_char[i].nPoints;
		}

		int k;
		k=1;
		yMax[k++] = 1.e0;
		yMax[k++] = 1.e16;
		for(int j=1;j<=kinetics->nCarbon;j++)
			yMax[k++] = 1.e0;
		for(int j=1;j<=kinetics->nAdsorbed;j++)
			yMax[k++] = 1.e0;

		k=1;
		yMin[k++] = 0.0;
		yMin[k++] = 0.0;
		for(int j=1;j<=kinetics->nCarbon;j++)
			yMin[k++] = 0.e0;
		for(int j=1;j<=kinetics->nAdsorbed;j++)
			yMin[k++] = 0.0;
		ChangeDimensions(numTotal, 1, &YY);

		for(i=1;i<=nCases;i++)
		{
			cout << "-----------------------------------------"		<< endl;
			cout << "Experiment #" << i << endl;
			cout << "-----------------------------------------"		<< endl;
			cout << "  * Start:     " << indices[i] << endl;
			cout << "  * x0:        " << initialconditions_x[i] << endl;

			int k=1;
			cout << "  * xChar:     " << initialconditions_y[i][k++] << endl;
			cout << "  * Sg:        " << initialconditions_y[i][k++] << endl;
			for(int j=1;j<=kinetics->nCarbon;j++)
				cout << "  * xc" << j << ":       " << initialconditions_y[i][k++] << endl;
			for(int j=1;j<=kinetics->nAdsorbed;j++)
				cout << "  * teta" << j << ":     " << initialconditions_y[i][k++] << endl;
			cout << endl;
		}

		Prepare();

		for(i=1;i<=nCases;i++)
			experiments_char[i].PrepareConcentrations(kinetics->speciesGas);
	}

	else if (solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)
	{
		int i;
		string dummy;
		ifstream fInput;
		openInputFileAndControl(fInput, file_name_list_of_files);

		for(;;)
		{
			fInput >> dummy;
			if (dummy != "#END")
				list_of_names_of_files.push_back(dummy);
			else
				break;
		}
		fInput.close();

		nCases = list_of_names_of_files.size();
		ChangeDimensions(nCases,									&indices);
		ChangeDimensions(nCases,									&initialconditions_x);
		ChangeDimensions(nCases, mix_chemkin->NumberOfSpecies(),	&initialconditions_y);
		ChangeDimensions(mix_chemkin->NumberOfSpecies(),			&yMin);
		ChangeDimensions(mix_chemkin->NumberOfSpecies(),			&yMax);
		
		// Set each experiment
		experiments_bio = new OpenSMOKE_SolidExperiment_DSmoke_Bio[nCases+1];
		for(i=1;i<=nCases;i++)
		{
			experiments_bio[i].ReadFromFile(list_of_names_of_files[i-1]);
			experiments_bio[i].PrepareMasses(*mix_chemkin);
		}

		// Assembling experiments
		numTotal = 0;
		for(i=1;i<=nCases;i++)
		{
			indices[i]				= numTotal+1;
			initialconditions_x[i]	= experiments_bio[i].x[1]; 
			
			for(int j=1;j<=mix_chemkin->NumberOfSpecies();j++)
				initialconditions_y[i][j] = experiments_bio[i].initial_masses[j]; 

			numTotal  += experiments_bio[i].nPoints;
		}

		for(int j=1;j<=mix_chemkin->NumberOfSpecies();j++)	yMax[j] = 1.;
		for(int j=1;j<=mix_chemkin->NumberOfSpecies();j++)	yMin[j] = 0.;
		
		ChangeDimensions(numTotal, 1, &YY);

		for(i=1;i<=nCases;i++)
		{
			cout << "-----------------------------------------"		<< endl;
			cout << "Experiment #" << i << endl;
			cout << "-----------------------------------------"		<< endl;
			cout << "  * Index:  " << indices[i] << endl;
			cout << "  * x0:     " << initialconditions_x[i] << endl;

			for(int j=1;j<=mix_chemkin->NumberOfSpecies();j++)
				cout << "  * y" << j << ":     " << experiments_bio[i].initial_masses[j] << endl;			
			cout << endl;
		}

		Prepare();
	}
}

void MyOpenSMOKE_SolidRegression::Prepare()
{
	ChangeDimensions(numTotal, &temperatures_global);
	ChangeDimensions(numTotal, &utemperatures_global);
	ChangeDimensions(numTotal, &pressures_Gas_global);

	int k=1;
	
	if (solid_regression == OPENSMOKE_REGRESSION_CHAR)
	{
		for(int i=1;i<=nCases;i++)
		{
			for(int j=1;j<=experiments_char[i].nPoints;j++)
			{
				temperatures_global[k]   = experiments_char[i].temperature_mean;
				utemperatures_global[k]  = 1./temperatures_global[k];
				pressures_Gas_global[k]  = experiments_char[i].pressure_atm;
				k++;
			}
		}
	}

	else if (solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)
	{
		for(int i=1;i<=nCases;i++)
		{
			for(int j=1;j<=experiments_bio[i].nPoints;j++)
			{
				temperatures_global[k]   = experiments_bio[i].temperature_mean;
				utemperatures_global[k]  = 1./temperatures_global[k];
				pressures_Gas_global[k]  = experiments_bio[i].pressure_atm;
				k++;
			}
		}
	}
	
	TMeanGlobal		= Mean(temperatures_global);
	PGasMeanGlobal   = Mean(pressures_Gas_global);

	C2 = -1./TMeanGlobal;
	C3 = pressures_Gas_global.Norm2();

	BzzVector auxiliar_vector = utemperatures_global;
	for(k=1;k<=numTotal;k++) auxiliar_vector[k] += C2;

	C1 = 1./auxiliar_vector.Norm2();

	cout << "Tmin:    " << temperatures_global.Min()	<< " K"		<< endl;
	cout << "Tmax:    " << temperatures_global.Max()	<< " K"		<< endl;
	cout << "Tmean:   " << TMeanGlobal					<< " K"		<< endl;
	cout << "PGasmin:  " << pressures_Gas_global.Min()	<< " atm"		<< endl;
	cout << "PGasmax:  " << pressures_Gas_global.Max()	<< " atm"		<< endl;
	cout << "PGasmean: " << PGasMeanGlobal				<< " atm"		<< endl;

	cout << "C1:    " << C1							<< " K"		<< endl;
	cout << "C2:    " << C2							<< " 1/K"	<< endl;
	cout << "C3:    " << C3							<< " atm"	<< endl;
}

void MyOpenSMOKE_SolidRegression::MemoryAllocation()
{	
	if (solid_regression == OPENSMOKE_REGRESSION_CHAR)
	{
		ChangeDimensions(kinetics->nAdsorbed,	&teta);
		ChangeDimensions(kinetics->nAdsorbed,	&Cadsorbed);
		ChangeDimensions(kinetics->nBulk,		&Cbulk);
		ChangeDimensions(kinetics->nAdsorbed,	&dteta_over_dt);
		
		ChangeDimensions(kinetics->nCarbon,		&xcarbon);
		ChangeDimensions(kinetics->nCarbon,		&dxcarbon_over_dt);
		ChangeDimensions(kinetics->nCarbon,		&teta_f);
		ChangeDimensions(kinetics->nCarbon,		&C_f);
	}

	else if (solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)
	{
		ChangeDimensions(mix_chemkin->NumberOfSpecies(),	&mass);
		ChangeDimensions(mix_chemkin->NumberOfSpecies(),	&dmass_over_dt);
		ChangeDimensions(mix_chemkin->NumberOfSpecies(),	&C);
		ChangeDimensions(mix_chemkin->NumberOfSpecies(),	&R);
	}
}

void MyOpenSMOKE_SolidRegression::PrepareOptimizerFile()
{
	openOutputFileAndControl(fOptimizer, "Output.out");
	fOptimizer.setf(ios::scientific);
	fOptimizer << "------------------------------------------------------------------------------" << endl;
	fOptimizer << "                           OpenSMOKE_SolidRegression                          " << endl;
	fOptimizer << "------------------------------------------------------------------------------" << endl;
	fOptimizer << endl;
	fOptimizer << "Tmin:     " << temperatures_global.Min()		<< " K"		<< endl;
	fOptimizer << "Tmax:     " << temperatures_global.Max()		<< " K"		<< endl;
	fOptimizer << "Tmean:    " << TMeanGlobal					<< " K"		<< endl;
	fOptimizer << "PGasmin:  " << pressures_Gas_global.Min()	<< " atm"		<< endl;
	fOptimizer << "PGasmax:  " << pressures_Gas_global.Max()	<< " atm"		<< endl;
	fOptimizer << "PGasmean: " << PGasMeanGlobal				<< " atm"		<< endl;
	fOptimizer << endl;
	fOptimizer << "Model: " << iModel						<< endl;
	fOptimizer << "C1:    " << C1							<< " K"		<< endl;
	fOptimizer << "C2:    " << C2							<< " 1/K"	<< endl;
	fOptimizer << "C3:    " << C3							<< " atm"	<< endl;
	fOptimizer << endl << endl;

	//openOutputFileAndControl(fOptimizer100, "Output100.out");
	//fOptimizer100.setf(ios::scientific);
	//for(j=1;j<=kinetics->number_of_reactions;j++)	fOptimizer100 << kinetics->A[j]	<< "\t";
}

void MyOpenSMOKE_SolidRegression::Run_Char()
{
	double start_time;
	double end_time;

	int numModels	= 1;
	int numX		= 1;
	int numY		= 1;

	BzzMatrix X(numTotal, numX);
	BzzMatrix Y(numTotal, numY);

	int i;
	int k=1;
	for(i=1;i<=nCases;i++)
	{
		for(int j=1;j<=experiments_char[i].nPoints;j++)
		{
			X[k][1]		= experiments_char[i].x[j]; 
			Y[k++][1]	= experiments_char[i].y[j]; 
		}
	}

	// First guess values
	     if (parameters == SOLID_REGRESSION_ALL)			nParameters = 3*kinetics->number_of_reactions;
	else if (parameters == SOLID_REGRESSION_A)				nParameters =   kinetics->number_of_reactions;
	else if (parameters == SOLID_REGRESSION_BETA)			nParameters =   kinetics->number_of_reactions;
	else if (parameters == SOLID_REGRESSION_E)				nParameters =   kinetics->number_of_reactions;
	else if (parameters == SOLID_REGRESSION_A_E)			nParameters = 2*kinetics->number_of_reactions;
	else if (parameters == SOLID_REGRESSION_A_BETA)			nParameters = 2*kinetics->number_of_reactions;
	else if (parameters == SOLID_REGRESSION_BETA_E)			nParameters = 2*kinetics->number_of_reactions;
	else if (parameters == SOLID_REGRESSION_USERDEFINED)	nParameters = bUserDefined.Rows();
    else if (parameters == SOLID_REGRESSION_ALL_PHI)        nParameters = nCases + 3*kinetics->number_of_reactions;   
    else if (parameters == SOLID_REGRESSION_A_E_PHI)        nParameters = nCases + 2*kinetics->number_of_reactions;   
	 

	BzzVector bFirstGuess(nParameters);
	kinetics->UpdateOptimizerParameters(iModel, parameters, bFirstGuess, bUserDefined, nCases, experiments_char);
	 

	BzzVector b = bFirstGuess;
	BzzVector bSolution = bFirstGuess;
	BzzVector bMin = bFirstGuess;
	BzzVector bMax = bFirstGuess;
	BzzVector bmin = bFirstGuess;
	BzzVector bmax = bFirstGuess;

	// Set Min-Max values
	{
		if (bUserDefined.Rows()!=nParameters)
			ErrorMessage("Wrong bUserDefined");
		
		for(int j=1;j<=bUserDefined.Rows();j++)
		{
			if (bUserDefined[j][2] == 1)	//A
			{
				bmin[j]= exp(bFirstGuess[j])* minimum_ratio_A;
				bmax[j]= exp(bFirstGuess[j])* maximum_ratio_A;
				bMin[j] = log(bmin[j]);
				bMax[j] = log(bmax[j]);
				
				//bMin[j] = bFirstGuess[j]* minimum_ratio_A;
				//bMax[j] = bFirstGuess[j]* maximum_ratio_A;

			}
			else if (bUserDefined[j][2] == 2)	//Beta
			{
				bMin[j] = bFirstGuess[j] + minimum_delta_BETA;
				bMax[j] = bFirstGuess[j] + maximum_delta_BETA;			
			}
			else if (bUserDefined[j][2] == 3)	//E
			{
				bMin[j] = bFirstGuess[j]*minimum_ratio_E;
				bMax[j] = bFirstGuess[j]*maximum_ratio_E;
			}
			else if (bUserDefined[j][2] == 4)	// PHI
			{
			    bMin[j] = bFirstGuess[j]*minimum_ratio_PHI;
				bMax[j] = bFirstGuess[j]*maximum_ratio_PHI;
			}
		}

	}

	if (iRegressionAnalysis == true)
	{
		PrepareOptimizerFile();

		// Non linear regression
		BzzNonLinearRegression nonLinReg(numModels, X, Y, MyBzzModelOdeRegression);
		
		// Options
		if (iVerbose == true)			nonLinReg.SavePartial("RegressionBest.out");
		if (iVerbose == true)			nonLinReg.SaveComplete("RegressionComplete.out");
		nonLinReg.InitializeModel(1,b, bMin, bMax);

		// Run Regression
		start_time = BzzGetCpuTime();
			if (iRobustAnalysis == false)	nonLinReg.LeastSquaresAnalysis();
			if (iRobustAnalysis == true)	nonLinReg.RobustAnalysis();
		end_time = BzzGetCpuTime();
		
		fOptimizer.close();
		//fOptimizer100.close();

		// Get Solution
		nonLinReg.GetSolution(b);

		cout << "NLR cpu time: " << end_time - start_time << " s" << endl;
	}
	
	// Update kinetic parameters
	kinetics->UpdateKineticParameters(iModel, parameters, b, bUserDefined, nCases, experiments_char);

	// Write on file
	WriteOptimizedParametersChar("Kinetics.out");

	// --------------------------------------------------------------------------------------------
	// Post-Processing
	// --------------------------------------------------------------------------------------------
	for(i=1;i<=nCases;i++)
	{
		iWriteOnFile = true;

		// Fitted data
		{
			kinetics->UpdateKineticParameters(iModel, parameters, b, bUserDefined, nCases, experiments_char);
			openOutputFileAndControl(fOutput, experiments_char[i].name_of_file + ".fitted");
			fOutput.setf(ios::scientific);
			LabelODE_File();

			//
			ptIndex		= i;									// [K]
			ptPressure	= experiments_char[i].pressure_atm;		// [atm]

			o.Deinitialize();
			o.SetInitialConditions(initialconditions_y.GetRow(i), initialconditions_x[i], &ode);
			o.StepPrint(ODE_Print);
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			if (iRelativeTolerance == true)	o.SetTollRel(relativeTolerance);
			if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

			o(experiments_char[i].x[experiments_char[i].nPoints], experiments_char[i].x[experiments_char[i].nPoints]);
			
			fOutput.close();
		}

		// Original data
		{
			kinetics->UpdateKineticParameters(iModel, parameters, bFirstGuess, bUserDefined, nCases, experiments_char);
			openOutputFileAndControl(fOutput, experiments_char[i].name_of_file + ".original");
			fOutput.setf(ios::scientific);
			LabelODE_File();

			//
			ptIndex		= i;									// [K]
			ptPressure	= experiments_char[i].pressure_atm;		// [atm]

			o.Deinitialize();
			o.SetInitialConditions(initialconditions_y.GetRow(i), initialconditions_x[i], &ode);
			o.StepPrint(ODE_Print);
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			if (iRelativeTolerance == true)	o.SetTollRel(relativeTolerance);
			if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

			o(experiments_char[i].x[experiments_char[i].nPoints], experiments_char[i].x[experiments_char[i].nPoints]);
			
			fOutput.close();
		}

		// Experimental data
		{
			ofstream fExp;
			openOutputFileAndControl(fExp, experiments_char[i].name_of_file + ".exp");
			fExp.setf(ios::scientific);
			for(int j=1;j<=experiments_char[i].nPoints;j++)
				fExp << experiments_char[i].x[j] << "\t" << experiments_char[i].y[j] << endl;
			fExp.close();
		}
	}	
}

void MyOpenSMOKE_SolidRegression::Run_DSmoke_Bio()
{
	double start_time;
	double end_time;

	int numModels	= 1;
	int numX		= 1;
	int numY		= 1;

	BzzMatrix X(numTotal, numX);
	BzzMatrix Y(numTotal, numY);

	int i;
	int k=1;
	for(i=1;i<=nCases;i++)
	{
		for(int j=1;j<=experiments_bio[i].nPoints;j++)
		{
			X[k][1]		= experiments_bio[i].x[j]; 
			Y[k++][1]	= experiments_bio[i].y[j]; 
		}
	}

	// First guess values
	if (parameters == SOLID_REGRESSION_ALL)					nParameters = 3*mix_chemkin->NumberOfReactions();
	else if (parameters == SOLID_REGRESSION_A)				nParameters =   mix_chemkin->NumberOfReactions();
	else if (parameters == SOLID_REGRESSION_BETA)			nParameters =   mix_chemkin->NumberOfReactions();
	else if (parameters == SOLID_REGRESSION_E)				nParameters =   mix_chemkin->NumberOfReactions();
	else if (parameters == SOLID_REGRESSION_A_E)			nParameters = 2*mix_chemkin->NumberOfReactions();
	else if (parameters == SOLID_REGRESSION_A_BETA)			nParameters = 2*mix_chemkin->NumberOfReactions();
	else if (parameters == SOLID_REGRESSION_BETA_E)			nParameters = 2*mix_chemkin->NumberOfReactions();
	else if (parameters == SOLID_REGRESSION_USERDEFINED)	nParameters = bUserDefined.Rows();
	

	BzzVector bFirstGuess(nParameters);
	mix_chemkin->kinetics.UpdateOptimizerParameters(iModel, parameters, bFirstGuess, bUserDefined);


	BzzVector b = bFirstGuess;
	BzzVector bSolution = bFirstGuess;
	BzzVector bMin = bFirstGuess;
	BzzVector bMax = bFirstGuess;
	BzzVector bmin = bFirstGuess;
	BzzVector bmax = bFirstGuess;


	// Set Min-Max values
	{
		if (bUserDefined.Rows()!=nParameters)
			ErrorMessage("Wrong bUserDefined");
		
		for(int j=1;j<=bUserDefined.Rows();j++)
		{
			if (bUserDefined[j][2] == 1)	//A
			{
				bmin[j]= exp(bFirstGuess[j])* minimum_ratio_A;
				bmax[j]= exp(bFirstGuess[j])* maximum_ratio_A;
				bMin[j] = log(bmin[j]);
				bMax[j] = log(bmax[j]);
				//cout<<j<<"\t"<<bMin[j]<<"\t"<< bFirstGuess[j]<<"\t"<<bMax[j]<<endl;
				//cout<<j<<"\t"<<bmin[j]<<"\t"<<exp(bFirstGuess[j])<<"\t"<<bmax[j]<<endl;

			}
			else if (bUserDefined[j][2] == 2)	//Beta
			{
				bMin[j] = bFirstGuess[j] + minimum_delta_BETA;
				bMax[j] = bFirstGuess[j] + maximum_delta_BETA;			
			}
			else if (bUserDefined[j][2] == 3)	//E
			{
				bMin[j] = bFirstGuess[j]*minimum_ratio_E;
				bMax[j] = bFirstGuess[j]*maximum_ratio_E;
			
			}
		}
     //getchar();
	}

	if (iRegressionAnalysis == true)
	{
		PrepareOptimizerFile();

		// Non linear regression
		BzzNonLinearRegression nonLinReg(numModels, X, Y, MyBzzModelOdeRegression);

		// Options
		if (iVerbose == true)			nonLinReg.SavePartial("RegressionBest.out");
		if (iVerbose == true)			nonLinReg.SaveComplete("RegressionComplete.out");
		nonLinReg.InitializeModel(1,b, bMin, bMax);
		

		// Run Regression
		start_time = BzzGetCpuTime();
			if (iRobustAnalysis == false)	nonLinReg.LeastSquaresAnalysis();
			if (iRobustAnalysis == true)	nonLinReg.RobustAnalysis();
		end_time = BzzGetCpuTime();

		fOptimizer.close();
		//fOptimizer100.close(); 

		// Get Solution
		nonLinReg.GetSolution(bSolution);

		cout << "NLR cpu time: " << end_time - start_time << " s" << endl;
	}
	
	// Update kinetic parameters
	mix_chemkin->kinetics.UpdateKineticParameters(iModel, parameters, bSolution, bUserDefined);

	WriteOptimizedParametersBio("Kinetics.out");

	// --------------------------------------------------------------------------------------------
	// Post-Processing
	// --------------------------------------------------------------------------------------------
	for(i=1;i<=nCases;i++)
	{
		iWriteOnFile = true;

		// Fitted data
		{
			mix_chemkin->kinetics.UpdateKineticParameters(iModel, parameters, bSolution, bUserDefined);
			openOutputFileAndControl(fOutput, experiments_bio[i].name_of_file + ".fitted");
			fOutput.setf(ios::scientific);
			LabelODE_File();

			// Pointers
			ptIndex		= i;									// [-]
			ptPressure	= experiments_bio[i].pressure_atm;		// [Pa]

			o.Deinitialize();
			o.SetInitialConditions(initialconditions_y.GetRow(i), initialconditions_x[i], &ode);
			o.StepPrint(ODE_Print);
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			if (iRelativeTolerance == true)	o.SetTollRel(relativeTolerance);
			if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

			o(experiments_bio[i].x[experiments_bio[i].nPoints], experiments_bio[i].x[experiments_bio[i].nPoints]);
			
			fOutput.close();
		}

		// Original data
		{
			mix_chemkin->kinetics.UpdateKineticParameters(iModel, parameters, bFirstGuess, bUserDefined);
			openOutputFileAndControl(fOutput, experiments_bio[i].name_of_file + ".original");
			fOutput.setf(ios::scientific);
			LabelODE_File();

			// Pointers
			ptIndex		= i;									// [-]
			ptPressure	= experiments_bio[i].pressure_atm;		// [Pa]

			o.Deinitialize();
			o.SetInitialConditions(initialconditions_y.GetRow(i), initialconditions_x[i], &ode);
			o.StepPrint(ODE_Print);
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			if (iRelativeTolerance == true)	o.SetTollRel(relativeTolerance);
			if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

			o(experiments_bio[i].x[experiments_bio[i].nPoints], experiments_bio[i].x[experiments_bio[i].nPoints]);
			
			fOutput.close();
		}

		// Experimental data
		{
			ofstream fExp;
			openOutputFileAndControl(fExp, experiments_bio[i].name_of_file + ".exp");
			fExp.setf(ios::scientific);
			for(int j=1;j<=experiments_bio[i].nPoints;j++)
				fExp << experiments_bio[i].x[j] << "\t" << experiments_bio[i].y[j] << endl;
			fExp.close();
		}
	}	
}

void MyOpenSMOKE_SolidRegression::MyODE_Print_Char(BzzVector &x, double t)
{
	if (iWriteOnFile == true)
	{
		int i;

		// stampo il tempo in minuti
		fOutput << setw(18) << left << t;
		
		fOutput << setw(18) << left << xchar;
		fOutput << setw(18) << left << experiments_char[ptIndex].dxchar_over_dt;
		fOutput << setw(18) << left << Sg;
		fOutput << setw(18) << left << dSg_over_dt;
		fOutput << setw(18) << left << ptT;
		fOutput << setw(18) << left << ptPressure;
		
		for(i=1;i<=kinetics->nAdsorbed;i++)
			fOutput << setw(18) << left << teta[i];

		for(i=1;i<=kinetics->nAdsorbed;i++)
			fOutput << setw(18) << left << dteta_over_dt[i];

		for(i=1;i<=kinetics->number_of_reactions;i++)
			fOutput << setw(18) << left << kinetics->rr[i];

		fOutput << setw(18) << left << kinetics->Rchar;

		for(i=1;i<=kinetics->nCarbon;i++)
			fOutput << setw(18) << left << kinetics->Rcarbon[i];

		for(i=1;i<=kinetics->nAdsorbed;i++)
			fOutput << setw(18) << left << kinetics->R[i];

		fOutput << endl;
	}
}

void MyOpenSMOKE_SolidRegression::MyODE_Print_DSmoke_Bio(BzzVector &x, double t)
{
	if (iWriteOnFile == true)
	{
		int i;

		fOutput << setw(18) << left << t;
		
		fOutput << setw(18) << left << mass_total_solid;
		fOutput << setw(18) << left << volume_total_solid;
		fOutput << setw(18) << left << ptT;
		fOutput << setw(18) << left << ptPressure;
		
		for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)
			fOutput << setw(18) << left << mass[i];

		for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)
			fOutput << setw(18) << left << dmass_over_dt[i];

		for(i=1;i<=mix_chemkin->NumberOfReactions();i++)
			fOutput << setw(18) << left << mix_chemkin->r[i];

		for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)
			fOutput << setw(18) << left << R[i];

		fOutput << endl;
	}
}

void MyOpenSMOKE_SolidRegression::LabelODE_File()
{
	if (iWriteOnFile == true)
	{

	if (solid_regression == OPENSMOKE_REGRESSION_CHAR)
	{
		int i;
		int count = 1;

		// Time
		fOutput << setw(18) << left << "t[s]" + GetLabelIndex(count++);
		//fOutput << setw(18) << left << "t[min]" + GetLabelIndex(count++);
		
		// Char
		fOutput << setw(18) << left << "xc[-]" + GetLabelIndex(count++);
		fOutput << setw(18) << left << "dxc_dt[1/s]" + GetLabelIndex(count++);
		
		// Surface
		fOutput << setw(18) << left << "Sg[-]" + GetLabelIndex(count++);
		fOutput << setw(18) << left << "dSg_dt[1/s]" + GetLabelIndex(count++);

		// Temperature and pressure
		fOutput << setw(18) << left << "T[K]" + GetLabelIndex(count++);
		fOutput << setw(18) << left << "P[atm]" + GetLabelIndex(count++);

		// Adsorbed species
		for(i=1;i<=kinetics->nAdsorbed;i++)
			fOutput << setw(18) << left << "Teta" + GetNumber(i) + "[-]" + GetLabelIndex(count++);
		for(i=1;i<=kinetics->nAdsorbed;i++)
			fOutput << setw(18) << left << "dTetadt" + GetNumber(i) + "[-]" + GetLabelIndex(count++);

		//for(i=1;i<=kinetics->nGas;i++)
		//	fOutput << setw(18) << left << "wGas[g /g char]" + GetNumber(i) + "[-]" + GetLabelIndex(count++);

		// Reaction rates
		for(i=1;i<=kinetics->number_of_reactions;i++)
			fOutput << setw(18) << left << "RR" + GetNumber(i) + "[kmol/m2/s]" + GetLabelIndex(count++);

		// Formation rates: char
		fOutput << setw(16) << left << "FRChar[kmol/m2/s]" + GetLabelIndex(count++);

		// Formation rates: carbon
		for(i=1;i<=kinetics->nCarbon;i++)
			fOutput << setw(18) << left << "FR" + GetNumber(i) + "Carb[kmol/m2/s]" + GetLabelIndex(count++);

		// Formation rates
		for(i=1;i<=kinetics->nAdsorbed;i++)
			fOutput << setw(18) << left << "FR" + GetNumber(i) + "[kmol/m2/s]" + GetLabelIndex(count++);

		fOutput << endl;
	}

	else if (solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)
	{
		int i;
		int count = 1;

		// Time
		fOutput << setw(18) << left << "t[s]" + GetLabelIndex(count++);
		
		// Particle
		fOutput << setw(18) << left << "msolid[kg]" + GetLabelIndex(count++);
		fOutput << setw(18) << left << "Vsolid[m3]" + GetLabelIndex(count++);
		fOutput << setw(18) << left << "T[K]" + GetLabelIndex(count++);
		fOutput << setw(18) << left << "P[atm]" + GetLabelIndex(count++);
		
		// Adsorbed species
		for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)
			fOutput << setw(18) << left << "m_" + mix_chemkin->names[i] + "[kg]" + GetLabelIndex(count++);

		// Adsorbed species
		for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)
			fOutput << setw(18) << left << "dm_" + mix_chemkin->names[i] + "[kg/s]" + GetLabelIndex(count++);

		// Reaction rates
		for(i=1;i<=mix_chemkin->NumberOfReactions();i++)
			fOutput << setw(18) << left << "RR" + GetNumber(i) + "[kmol/m3/s]" + GetLabelIndex(count++);

		// Formation rates
		for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)
			fOutput << setw(18) << left << "FR_" + mix_chemkin->names[i] + "[kg/m3/s]" + GetLabelIndex(count++);

		fOutput << endl;
	}

	}

}

void MyBzzModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	if (ptRegression->solid_regression == OPENSMOKE_REGRESSION_CHAR)			ptRegression->ModelOdeRegression_Char(model, ex, b, x, y);
	else if (ptRegression->solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)	ptRegression->ModelOdeRegression_DSmoke_Bio(model, ex, b, x, y);
}

void MyOpenSMOKE_SolidRegression::ModelOdeRegression_Char(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	int i;

	if(ex == 1)
	{
		int j;

		// From Regression parameters to Kinetic parameters
		kinetics->UpdateKineticParameters(iModel, parameters, b, bUserDefined, nCases, experiments_char);

		// Write on file
		{
			fOptimizer << setw(8) << left << count_global	<< "\t";
			for(j=1;j<=kinetics->number_of_reactions;j++)	fOptimizer << kinetics->A[j]	<< "\t";
			for(j=1;j<=kinetics->number_of_reactions;j++)	fOptimizer << kinetics->Beta[j]	<< "\t";
			for(j=1;j<=kinetics->number_of_reactions;j++)	fOptimizer << kinetics->E[j]	<< "\t";
			fOptimizer << endl;
		}

      /* if ((count_global-1)%100 == 0)
		{
		    fOptimizer100 << setw(8) << left << count_global	<< "\t";
			for(j=1;j<=kinetics->number_of_reactions;j++)	fOptimizer100 << kinetics->A[j]	<< "\t";
			for(j=1;j<=kinetics->number_of_reactions;j++)	fOptimizer100 << kinetics->Beta[j]	<< "\t";
			for(j=1;j<=kinetics->number_of_reactions;j++)	fOptimizer100 << kinetics->E[j]	<< "\t";
			fOptimizer100 << endl;
		}
*/
		// Write on video
		if ((count_global-1)%100 == 0)	cout << count_global << endl;

		// Update global counter
		count_global++;

		for(i=1;i<=nCases;i++)
		{
			int j;

			ptIndex		= i;
			ptPressure	= experiments_char[i].pressure_atm;

			// Solve ODE System
			o.Deinitialize();
			o.SetInitialConditions(initialconditions_y.GetRow(i), initialconditions_x[i], &ode);
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			if (iRelativeTolerance == true)	o.SetTollRel(relativeTolerance);
			if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

		    for(j=1;j<=experiments_char[i].nPoints;j++)
			{
				if (experiments_char[i].kind_of_experimental_data == Char_CONVERSION_CHAR)
				{
                    yy = o(experiments_char[i].x[j]);
				    YY.SetRow(indices[i]+j-1, yy[1]); 
				}

				else if (experiments_char[i].kind_of_experimental_data == Char_D_CONVERSION_CHAR)	
				{
					yy = o(experiments_char[i].x[j]);
					YY.SetRow(indices[i]+j-1, experiments_char[i].dxchar_over_dt);

				}
			}
		}
	}

	// Recover regression data
	YY.GetRow(ex, &y);
}

void MyOpenSMOKE_SolidRegression::ModelOdeRegression_DSmoke_Bio(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y)
{
	int i;

	if(ex == 1)
	{
		int j;

		// From Regression parameters to Kinetic parameters
		mix_chemkin->kinetics.UpdateKineticParameters(iModel, parameters, b, bUserDefined);

		// Write on file
		{
			fOptimizer << setw(8) << left << count_global	<< "\t";
			for(j=1;j<=mix_chemkin->NumberOfReactions();j++)	fOptimizer << mix_chemkin->kinetics.GiveMe_A(j)	<< "\t";
			for(j=1;j<=mix_chemkin->NumberOfReactions();j++)	fOptimizer << mix_chemkin->kinetics.GiveMe_Beta(j)	<< "\t";
			for(j=1;j<=mix_chemkin->NumberOfReactions();j++)	fOptimizer << mix_chemkin->kinetics.GiveMe_E(j)	<< "\t";
			fOptimizer << endl;
		}

		// Write on video
		if ((count_global-1)%100 == 0)	cout << count_global << endl;

		// Update global counter
		count_global++;

		for(i=1;i<=nCases;i++)
		{
			int j;

			ptIndex		= i;										// [-]
			ptPressure	= experiments_bio[i].pressure_atm;			// [Pa]

			// Solve ODE System
			o.Deinitialize();
			o.SetInitialConditions(initialconditions_y.GetRow(i), initialconditions_x[i], &ode);
			o.SetMinimumConstraints(yMin);
			o.SetMaximumConstraints(yMax);
			if (iRelativeTolerance == true)	o.SetTollRel(relativeTolerance);
			if (iAbsoluteTolerance == true)	o.SetTollAbs(absoluteTolerance);

			for(j=1;j<=experiments_bio[i].nPoints;j++)
			{
				if (experiments_bio[i].kind_of_experimental_data == DSMOKE_BIO_MASS_RATIO)
				{
					yy = o(experiments_bio[i].x[j]);
					double sum = 0.;
					for(int ii=1;ii<=experiments_bio[i].mass_ratio_indices.Size();ii++)
						sum+= yy[experiments_bio[i].mass_ratio_indices[ii]];
					YY.SetRow(indices[i]+j-1, sum); 
				}
				else if (experiments_bio[i].kind_of_experimental_data == DSMOKE_BIO_D_MASS_RATIO)
				{
					yy = o(experiments_bio[i].x[j]);
					double sum = 0.;
					for(int ii=1;ii<=experiments_bio[i].mass_ratio_indices.Size();ii++)
						sum+= dmass_over_dt[experiments_bio[i].mass_ratio_indices[ii]];
					YY.SetRow(indices[i]+j-1, sum); 
				}
			}
		}
	}

	// Recover regression data
	YY.GetRow(ex, &y);
}

void MyOpenSMOKE_SolidRegression::GetModel_01(BzzVector &y, double t, BzzVector &f)
{
	int k, i;

	// Recover unknowns
	k=1;
	xchar	= y[k++];								// [-]
	Sg		= y[k++];								// [m2/kg]
	for(i=1;i<=kinetics->nCarbon;i++)
		xcarbon[i] = y[k++];						// [-]
	for(i=1;i<=kinetics->nAdsorbed;i++)
		teta[i] = y[k++];							// [-]
	
	// Total carbon conversion
	double xcarbon_tot = xcarbon.GetSumElements();		// [-]
	
	// Free sites
	double tetaf_tot = 1. - teta.GetSumElements();		// [-]
	for(i=1;i<=kinetics->nCarbon;i++)
		teta_f[i] = tetaf_tot*xcarbon[i]/xcarbon_tot;

	// Surface concentrations
	for(i=1;i<=kinetics->nAdsorbed;i++)
//		Cadsorbed[i] = teta[i]*experiments_char[ptIndex].Sg0/5.0e4*Constants::Nav_kmol/kinetics->eta_adsorbed[i];	// [kmol/m2]
		Cadsorbed[i] = teta[i]*experiments_char[ptIndex].Sigma/Constants::Nav_kmol/kinetics->eta_adsorbed[i];	// [kmol/m2]
	
	// Free sites concentration
	for(i=1;i<=kinetics->nCarbon;i++)
		C_f[i] = teta_f[i]*experiments_char[ptIndex].Sigma/Constants::Nav_kmol/kinetics->eta_f[i];						// [kmol/m2]

	// Bulk concentrations
	Cbulk = 1.;															// [kmol/m2]



	// Update formation rates
	ptT = experiments_char[ptIndex].GetTemperature(t);
	experiments_char[ptIndex].UpdateGasConcentrations(ptT);
	kinetics->UpdateKineticConstants(ptT);



	kinetics->UpdateFormationRates(C_f, Cbulk, experiments_char[ptIndex].gas_c, Cadsorbed);				// [kmol/m2/s]


	if (iSurfaceEquation == true)
	{
		// Equation 1: Char conversion
		experiments_char[ptIndex].dxchar_over_dt = experiments_char[ptIndex].Sg0*Sg * MWchar * kinetics->Rchar;				// [1/s]

		// Equation 2: Surface
		dSg_over_dt = ( experiments_char[ptIndex].PHI*(1.-xchar)/2./Sg - Sg/(1.-xchar+1.e-12) ) * experiments_char[ptIndex].dxchar_over_dt;		// [m2/kg/s]

		// Equation 3: Carbon
		for(i=1;i<=kinetics->nCarbon;i++)
			dxcarbon_over_dt[i] = experiments_char[ptIndex].Sg0*Sg * MWchar * kinetics->Rcarbon[i];	// [1/s]

	//     Equation Tiziano: Gas
	//	for(i=1;i<=nGas;i++)
	//		dwgas_over_dt[i] = experiments_char[ptIndex].Sg0*Sg * MWgas[i] * kinetics->Rgas[i];	// [1/s]



	//	getchar();

		// Equation 4: Adsorbed 
		for(i=1;i<=kinetics->nAdsorbed;i++)
			dteta_over_dt[i] = - teta[i]/Sg*dSg_over_dt + Constants::Nav_kmol*kinetics->eta_adsorbed[i]/experiments_char[ptIndex].Sigma*kinetics->R[i] ;	
	}
	else
	{
		// Equation 1: Char conversion
		experiments_char[ptIndex].dxchar_over_dt = experiments_char[ptIndex].Sg0*Sg * MWchar * kinetics->Rchar;				// [1/s]

		// Equation 2: Surface
		dSg_over_dt = 0;													// [m2/kg/s]

		// Equation 3: Carbon
		for(i=1;i<=kinetics->nCarbon;i++)
			dxcarbon_over_dt[i] = experiments_char[ptIndex].Sg0*Sg * MWchar * kinetics->Rcarbon[i];						// [1/s]

		// Equation 3...Nadsorbed+2: 
		for(i=1;i<=kinetics->nAdsorbed;i++)
			dteta_over_dt[i] = - teta[i]/Sg*dSg_over_dt + Constants::Nav_kmol*kinetics->eta_adsorbed[i]/experiments_char[ptIndex].Sigma*kinetics->R[i] ;	
	}


	// Recover residuals
	k=1;
	f[k++] = experiments_char[ptIndex].dxchar_over_dt;
	f[k++] = dSg_over_dt;
	for(i=1;i<=kinetics->nCarbon;i++)
		f[k++] = dxcarbon_over_dt[i];
	for(i=1;i<=kinetics->nAdsorbed;i++)
		f[k++] = dteta_over_dt[i];

	if (xchar > 0.999)	f = 0.;
}


void MyOpenSMOKE_SolidRegression::GetModel_02(BzzVector &y, double t, BzzVector &f)
{
	int i;
	
	// Recover unknowns
	for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)		mass[i] = y[i];			// [-]
	
	// Total mass
	mass_total_solid = 0.;
	for(i=1;i<=indices_solid.Size();i++)
		mass_total_solid += mass[indices_solid[i]]; // la massa è la frazione massiva sul totale

	// Total volume
	volume_total_solid = mass_total_solid / rho_solid;

	// Concentrations
	for(i=1;i<=indices_solid.Size();i++)
		//C[indices_solid[i]] = mass[indices_solid[i]]/mix_chemkin->M[indices_solid[i]]/volume_total_solid;   
	    C[indices_solid[i]] = mass[indices_solid[i]]/mix_chemkin->M[indices_solid[i]];
	
	//C[indices_solid[i]]= non è esattamente una concentrazione Kg_i/kg_tot*mol_i/kg_i=mol_i/Kg_tot
	// la densita' compare nel calcolo della velocità di reazione: (per esempio reazione j)
	// dmass_over_dt[i]_j= nu_i_j * kj *rho^(sumNu_j - 1)*prod (C[indices_solid[i]]^nu_i)_j
	// R_j = nu_i * kj *rho^(sumNu_j - 1)*prod (C[indices_solid[i]]^nu_i)

	double cTot = C.GetSumElements();

	// Corrections
	BzzVector sumNu = mix_chemkin->kinetics.GiveMeSumNuDirect();
	BzzVector correction(mix_chemkin->NumberOfReactions());
	for(i=1;i<=mix_chemkin->NumberOfReactions();i++)
		correction[i] = pow(rho_solid, sumNu[i]-1.);

	// Reaction rates and formation rates
	ptT = experiments_bio[ptIndex].GetTemperature(t);
	mix_chemkin->kinetics.UpdateOmegaC(experiments_bio[ptIndex].GetOmegaC());
	mix_chemkin->ComputeKineticParameters(ptT, ptT, 1./ptT, ptPressure*101325.);	
	//mix_chemkin->ComputeFromConcentrations( ptT, C, cTot, &R);					// [kmol/m3/s]

	mix_chemkin->ComputeFromConcentrationsAndCorrect( ptT, C, cTot, &R, correction);	// [kmol/m3/s]

	ElementByElementProduct(R, mix_chemkin->M, &R);										// [kg/m3/s]
	
	
	//sumNu.BzzPrint("sumNu");
	//exit(-1);


	// Equations
	for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)		
		//dmass_over_dt[i] = R[i]*volume_total_solid;			// [-]  
       dmass_over_dt[i] = R[i];

	// Recover residuals
	for(i=1;i<=mix_chemkin->NumberOfSpecies();i++)
		f[i] = dmass_over_dt[i];	
}


void MyOdeSystem::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	if		(ptMyRegression->solid_regression == OPENSMOKE_REGRESSION_CHAR)			ptMyRegression->GetModel_01(x, t, f);
	else if (ptMyRegression->solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)	ptMyRegression->GetModel_02(x, t, f);
}

void MyOdeSystem::assign(MyOpenSMOKE_SolidRegression *myregression)
{
	ptMyRegression = myregression;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//									DICTIONARY													   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

OpenSMOKE_Dictionary_SolidRegression::OpenSMOKE_Dictionary_SolidRegression()
{
    SetupBase();
	SetName("OpenSMOKE_SolidRegression_Dictionary");

	// Compulsory options
	Add("#List",			        'C', 'S', "List of files containing experimental data");

	Add("#Regression",				'O', 'N', "Regression Analysis");
	Add("#SingleCase",				'O', 'N', "Single case running");

	// Options
	Add("#NoSurfaceEquation",		'O', 'N', "No surface equation");
	Add("#RobustAnalysis",			'O', 'N', "Robust Analysis");

	Add("#Parameters",		        'O', 'V', "Parameters (default ALL)");
	Add("#Model",					'O', 'I', "Model (default 1)");
 
	Add("#AbsoluteTolerance",		'O', 'D', "Absolute tolerance");
	Add("#RelativeTolerance",		'O', 'D', "Relative tolerance");

	Add("#Verbose",					'O', 'N', "Verbose");

	Add("#SolidSpecies",			'O', 'V', "List of solid species");
	Add("#SolidDensity",			'O', 'M', "Solid density");

	Add("#RemoveWarning",			'O', 'N', "Remove Warning message");


	Add("#MinimumRatio_A",			'O', 'D', "Minimum ratio for frequency factors");
	Add("#MaximumRatio_A",			'O', 'D', "Maximum ratio for frequency factors");
	Add("#MinimumRatio_E",			'O', 'D', "Minimum ratio for activation energy");
	Add("#MaximumRatio_E",			'O', 'D', "Maximum ratio for activation energy");
	Add("#MinimumDelta_BETA",		'O', 'D', "Minimum delta for beta");
	Add("#MaximumDelta_BETA",		'O', 'D', "Maximum delta for beta");


	// Compulsory and Conflicts
	Compulsory("#Regression",   "#SingleCase");
    Conflict("#SingleCase",		"#Regression");
    Conflict("#SingleCase",		"#RobustAnalysis");

    Lock();
}

void MyOpenSMOKE_SolidRegression::CheckDictionary(OpenSMOKE_Dictionary_SolidRegression &dictionary)
{
    int     int_value;
	double  double_value;
    string  string_value;
	vector<string> string_vector;
	vector<double> double_vector;

	if (dictionary.Return("#List", string_value))
        AssignListOfFiles(string_value);

	if (dictionary.Return("#Regression"))
        SetRegressionAnalysis(true);

	if (dictionary.Return("#SingleCase"))
        SetRegressionAnalysis(false);

	if (dictionary.Return("#NoSurfaceEquation"))
        SetSurfaceEquation(false);

	if (dictionary.Return("#RobustAnalysis"))
        SetRobustAnalysis(true);

	if (dictionary.Return("#Model", int_value))
		SetModel(int_value);

	if (dictionary.Return("#Parameters", string_vector))
		SetParameters(string_vector);

	if (dictionary.Return("#AbsoluteTolerance", double_value))
        SetAbsoluteTolerance(double_value);

	if (dictionary.Return("#RelativeTolerance", double_value))
        SetRelativeTolerance(double_value);

	if (dictionary.Return("#Verbose"))
		SetVerbose(true);

	if (dictionary.Return("#RemoveWarning"))
		SetRemoveWarning(true);

	if (dictionary.Return("#MinimumRatio_A", double_value))
        SetMinimumRatio_A(double_value);

	if (dictionary.Return("#MaximumRatio_A", double_value))
        SetMaximumRatio_A(double_value);

	if (dictionary.Return("#MinimumRatio_E", double_value))
        SetMinimumRatio_E(double_value);

	if (dictionary.Return("#MaximumRatio_A", double_value))
        SetMaximumRatio_E(double_value);

	if (dictionary.Return("#MinimumDelta_BETA", double_value))
        SetMinimumDelta_BETA(double_value);

	if (dictionary.Return("#MaximumDelta_BETA", double_value))
        SetMaximumDelta_BETA(double_value);

	if (solid_regression == OPENSMOKE_REGRESSION_CHAR)
	{
	}

	if (solid_regression == OPENSMOKE_REGRESSION_DSMOKE_BIO)
	{
		if (dictionary.Return("#SolidSpecies", string_vector))
			SetSolidSpecies(string_vector);
		else ErrorMessage("#SolidSpecies option is mandatory...");

		if (dictionary.Return("#SolidDensity", double_value, string_value))
			SetSolidDensity(string_value, double_value);
		else ErrorMessage("#SolidDensity option is mandatory...");
	}
}

void MyOpenSMOKE_SolidRegression::WriteOptimizedParametersChar(const string file_name)
{
	int j;

	cout << "--------------------------" << endl;
	cout << " Solution                 " << endl;
	cout << "  A: [kmol,m3,s]          " << endl;
	cout << "  E: [cal/mol]            " << endl;
	cout << "--------------------------" << endl;

	for(j=1;j<=kinetics->number_of_reactions;j++)
		cout << "A[" << j << "] = " << kinetics->A[j] << endl;

	for(j=1;j<=kinetics->number_of_reactions;j++)
		cout << "Beta[" << j << "] = " << kinetics->Beta[j] << endl;

	for(j=1;j<=kinetics->number_of_reactions;j++)
		cout << "E[" << j << "] = " << kinetics->E[j]/Constants::Conversion_cal_J << endl;

	cout << "--------------------------------------------------------------------" << endl;
	cout << "  PHI    " << endl;

	for(j=1;j<=nCases;j++)
		cout << "PHI[" << j << "] = " << experiments_char[j].PHI << endl;

	cout << endl;

	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);
	fOutput.setf(ios::scientific);

	fOutput << "-----------------------------------------------------------" << endl;
	fOutput << "                    DSMOKE Reactions                       " << endl;
	fOutput << "                     A: [kmol,m3,s]                        " << endl;
	fOutput << "                      E: [cal/mol]                         " << endl;
	fOutput << "-----------------------------------------------------------" << endl;

	for(j=1;j<=kinetics->number_of_reactions;j++)
	{
		fOutput << "R" << setw(6) << left << j;

		fOutput << setw(16) << left << kinetics->A[j];
		fOutput << setw(16) << left << kinetics->Beta[j];
		fOutput << setw(16) << left << kinetics->E[j]/Constants::Conversion_cal_J;
		fOutput << endl;
	}
	fOutput << endl << endl;

	fOutput << "-----------------------------------------------------------" << endl;
	fOutput << "                    Alberto Reactions                      " << endl;
	fOutput << "                     A: [kmol,m3,s]                        " << endl;
	fOutput << "                      E: [J/Kmol]                          " << endl;
	fOutput << "-----------------------------------------------------------" << endl;

	for(j=1;j<=kinetics->number_of_reactions;j++)
	{
		//fOutput.setf(ios::fixed); 
		
		fOutput << "R" << setw(6) << left << j;
		fOutput << setw(16) << setprecision(5) << left << kinetics->A[j];
		fOutput << setw(16) << setprecision(2) << left << kinetics->Beta[j];
		fOutput << setw(16) << setprecision(5) << left << kinetics->E[j];
		fOutput << endl;
	}
	fOutput << endl << endl;


	for(j=1;j<=nCases;j++)
		fOutput << "PHI[" << j << "] = " << experiments_char[j].PHI << endl;

	fOutput.close();
}

void MyOpenSMOKE_SolidRegression::WriteOptimizedParametersBio(const string file_name)
{
	int j;

	cout << "--------------------------" << endl;
	cout << " Solution" << endl;
	cout << "--------------------------" << endl;

	for(j=1;j<=mix_chemkin->NumberOfReactions();j++)
		cout << " A[" << j << "] = " << mix_chemkin->kinetics.GiveMe_A(j) << endl;

	for(j=1;j<=mix_chemkin->NumberOfReactions();j++)
		cout << " Beta[" << j << "] = " << mix_chemkin->kinetics.GiveMe_Beta(j) << endl;

	for(j=1;j<=mix_chemkin->NumberOfReactions();j++)
		cout << " E[" << j << "] = " << mix_chemkin->kinetics.GiveMe_E(j) << endl;


	ofstream fOutput;
	openOutputFileAndControl(fOutput, file_name);
	fOutput.setf(ios::scientific);

	fOutput << "-----------------------------------------------------------" << endl;
	fOutput << "                    CHEMKIN Reactions                      " << endl;
	fOutput << "                      A:[mol,cm3,s]                        " << endl;
	fOutput << "                      E:[cal/mol]                          " << endl;
	fOutput << "-----------------------------------------------------------" << endl;

	for(j=1;j<=mix_chemkin->NumberOfReactions();j++)
	{
		fOutput << "R" << setw(6) << left << j;
		fOutput << setw(16) << left << mix_chemkin->kinetics.GiveMe_A(j)*pow(1000., (mix_chemkin->kinetics.GiveMe_ForwardOrder(j)-1.));
		fOutput << setw(16) << left << mix_chemkin->kinetics.GiveMe_Beta(j);
		fOutput << setw(16) << left << mix_chemkin->kinetics.GiveMe_E(j)/Constants::Conversion_cal_J;
		fOutput << endl;
	}
	fOutput << endl << endl;


	fOutput << "-----------------------------------------------------------" << endl;
	fOutput << "                    DSMOKE Reactions                       " << endl;
	fOutput << "                      A:[kmol,m3,s]                        " << endl;
	fOutput << "                      E:[cal/mol]                          " << endl;
	fOutput << "-----------------------------------------------------------" << endl;

	for(j=1;j<=mix_chemkin->NumberOfReactions();j++)
	{
		int exponent = int(log10(mix_chemkin->kinetics.GiveMe_A(j)));
		fOutput.setf(ios::fixed); 
		fOutput << "R" << setw(6) << left << j;
		fOutput << "/";
		fOutput << setprecision(1) << left << left << fixed << mix_chemkin->kinetics.GiveMe_A(j)/pow(10., double(exponent));
		fOutput << "/";
		fOutput << setw(1) << left << fixed << exponent;
		fOutput << "/";
		fOutput << setw(1) << left << fixed << mix_chemkin->kinetics.GiveMe_Beta(j);
		fOutput <<"/";
		fOutput << setprecision(2) << left << fixed << mix_chemkin->kinetics.GiveMe_E(j)/Constants::Conversion_cal_J;
		fOutput << "/";
		fOutput << endl;
	}

	fOutput.close();
}