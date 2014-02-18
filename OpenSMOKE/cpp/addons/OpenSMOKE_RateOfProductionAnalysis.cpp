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

#include <vector>
#include <iomanip>
#include <sstream>
#include "basic/OpenSMOKE_Utilities.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "addons/OpenSMOKE_RateOfProductionAnalysis.h"

void OpenSMOKE_RateOfProductionAnalysis::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_RateOfProductionAnalysis"	<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_RateOfProductionAnalysis::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_RateOfProductionAnalysis"	<< endl;
    cout << "Object: "		<< name_object		<< endl;
    cout << "Warning:  "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

void OpenSMOKE_RateOfProductionAnalysis::SetName(const string name)
{
	name_object = name;
}

OpenSMOKE_RateOfProductionAnalysis::OpenSMOKE_RateOfProductionAnalysis()
{
	name_object	= "name not assigned";
}

void OpenSMOKE_RateOfProductionAnalysis::Initialize(OpenSMOKE_ReactingGas *_mix, BzzVectorInt &_indices)
{
	mix	= _mix;
	indices = _indices;

	NC = mix->NumberOfSpecies();
	NR = mix->NumberOfReactions();

	stoichiometry = new OpenSMOKE_Nu(mix);

	ChangeDimensions(NC, &sumIntegralProduction);
	ChangeDimensions(NC, &sumIntegralDestruction);
	
	I.Initialize(stoichiometry->number_of_reactions);

	stoichiometry->SparsityPattern(Matrix_Ip);
	stoichiometry->SparsityPattern(Matrix_Id);
}

void OpenSMOKE_RateOfProductionAnalysis::SetNumberOfPoints(const int _N)
{
	cout << "Number of points: " << _N << endl;
	cout << "Number of species: " << NC << endl;
	cout << "Number of non zero coefficients per point: " << stoichiometry->number_of_reactions.GetSumElements() << endl;
	cout << "Total number of non zero coefficients:     " << _N*stoichiometry->number_of_reactions.GetSumElements() << endl;
	cout << "Total memory:                              " << _N*stoichiometry->number_of_reactions.GetSumElements()*8/1e6 << " MB" << endl;

	N  = _N;

	cout << "Allocating memory for Rate of Production coefficients..." << endl;
	C = new OpenSMOKE_RateOfProductionCoefficient[N+1];

	cout << "Allocating additional memory..." << endl;
	ChangeDimensions(N, NC, &sumProduction);
	ChangeDimensions(N, NC, &sumDestruction);

	cout << "Initialize Rate of Production coefficients..." << endl;
	for(int i=1;i<=N;i++)
		C[i].Initialize(stoichiometry->number_of_reactions);

	cout << "Memory allocation correctly performed..." << endl;
}

void OpenSMOKE_RateOfProductionAnalysis::Run(BzzMatrix &r, BzzVector &x)
{
	Run(r, x, -1, -1);
}

void OpenSMOKE_RateOfProductionAnalysis::Run(BzzMatrix &r, BzzVector &x, const double local_coordinate)
{
	Run(r, x, local_coordinate, -1);
}

void OpenSMOKE_RateOfProductionAnalysis::Run(BzzMatrix &r, BzzVector &x, const double xA, const double xB)
{
	int i;
	int pStart = 1;
	int pEnd = N;

	I.Clean();
	sumProduction  = 0.;
	sumDestruction = 0.;
	sumIntegralProduction  = 0.;
	sumIntegralDestruction = 0.;

	// Integral Analysis
	if ( (xA == -1) && (xB == -1) )
	{
		pStart	= 1;
		pEnd	= N;
	}

	// Region Analysis
	else if ( (xA >= 0.) && (xB >= 0.) )
	{
		for(int p=1;p<=N;p++)
			if (x[p]>=xA)
			{
				pStart = p;
				break;
			}

		for(int p=N;p>=pStart;p--)
			if (x[p]<=xB)
			{
				pEnd = p;
				break;
			}
	}

	// Local Analysis
	else if ( (xA >= 0.) && (xB == -1) )
	{
		for(int p=1;p<=N;p++)
			if (x[p]>=xA)
			{
				pStart = p;
				pEnd = p;
				break;
			}
	}

	else ErrorMessage("No options are available for rate of production analysis...");

	for(int p=pStart;p<=pEnd;p++)
	{
		for(int j=1;j<=NC;j++)
		{			
			for(i=1;i<=stoichiometry->nu[j].nReactions;i++)
			{
				double element = stoichiometry->nu[j].nuReactions[i]*r[p][stoichiometry->nu[j].iReactions[i]];
				if (element > 0.)
					C[p].p[j][i] =  element;
				else
					C[p].d[j][i] = -element;
			}

			sumProduction[p][j]		=  C[p].p[j].GetSumElements();
			sumDestruction[p][j]	=  C[p].d[j].GetSumElements();
		}
	}

	if (N>1)	// 1-D systems
	{
		if (pEnd > pStart)	// Global
		{
			for(int p=pStart;p<=pEnd-1;p++)
			{
				for(int j=1;j<=NC;j++)
				{
					sumIntegralProduction[j]  +=  (sumProduction[p][j]  + sumProduction[p+1][j])  *(x[p+1]-x[p])/2.;
					sumIntegralDestruction[j] +=  (sumDestruction[p][j] + sumDestruction[p+1][j]) *(x[p+1]-x[p])/2.;

					for(i=1;i<=stoichiometry->nu[j].nReactions;i++)
					{
						I.p[j][i] += (C[p].p[j][i]+C[p+1].p[j][i])*(x[p+1]-x[p])/2.;
						I.d[j][i] += (C[p].d[j][i]+C[p+1].d[j][i])*(x[p+1]-x[p])/2.;	
					}
				}
			}
		}
		
		else	// Local
		{
				int p = pStart;
				for(int j=1;j<=NC;j++)
				{
					sumIntegralProduction[j]  =  sumProduction[p][j];
					sumIntegralDestruction[j] =  sumDestruction[p][j];

					for(i=1;i<=stoichiometry->nu[j].nReactions;i++)
					{
						I.p[j][i] = C[p].p[j][i];
						I.d[j][i] = C[p].d[j][i];	
					}
				}
		}
	}

	else		// 0-D Systems
	{
		sumIntegralProduction  =  sumProduction.GetRow(1);
		sumIntegralDestruction =  sumDestruction.GetRow(1);
		I.p = C[1].p;
		I.d = C[1].d;	
	}
/*
	for(int p=pStart;p<=pEnd;p++)
	{
		for(int j=1;j<=NC;j++)		
			for(i=1;i<=stoichiometry->nu[j].nReactions;i++)
			{
				if (sumProduction[p][j] != 0.)	C[p].p[j][i] /= sumProduction[p][j];
				else							C[p].d[j][i]  = 0.;
	
				if (sumDestruction[p][j] != 0.)	C[p].d[j][i] /= sumDestruction[p][j];
				else							C[p].d[j][i]  = 0.;
			}
	}

	for(int j=1;j<=NC;j++)
	{
		for(i=1;i<=stoichiometry->nu[j].nReactions;i++)
		{
			if (sumIntegralProduction[j] > 0.)	I.p[j][i] /= sumIntegralProduction[j];
			else								I.p[j][i]  = 0.;
			if (sumIntegralDestruction[j] > 0.)	I.d[j][i] /= sumIntegralDestruction[j];	
			else								I.d[j][i]  = 0.;
		}
	}
*/
	IntegralRateOfProductionAnalyses();
}

void OpenSMOKE_RateOfProductionAnalysis::IntegralRateOfProductionAnalyses()
{
	int i,j;
	double val;
	double *ptrVal;
	int count;

	// Build sparse matrix (production)
	count=1;
	Matrix_Ip.BeginScanning();
	while(ptrVal = Matrix_Ip.Scanning(&i,&j,&val))
	{
		*ptrVal = I.p[stoichiometry->sparse_rows[count]][stoichiometry->sparse_columns[count]];
		count++;
	}	

	// Build sparse matrix (consumption)
	count=1;
	Matrix_Id.BeginScanning();
	while(ptrVal = Matrix_Id.Scanning(&i,&j,&val))
	{
		*ptrVal = I.d[stoichiometry->sparse_rows[count]][stoichiometry->sparse_columns[count]];
		count++;
	}	
}

void OpenSMOKE_RateOfProductionAnalysis::PrintRateOfProductionAnalyses(const string fileName, BzzVector  &x, BzzVector &T)
{
	ofstream fOutput;
	openOutputFileAndControl(fOutput, fileName);
	fOutput.setf(ios::scientific);

	fOutput << setw(20) << left << "x[cm](1)";
	fOutput << setw(20) << left << "T[K](2)";
	PrintRateOfProductionAnalyses_Label(fOutput, 3);
	
	for(int p=1;p<=N;p++)
	{
		fOutput << setw(20) << left << x[p];
		fOutput << setw(20) << left << T[p];
	
		for(int j=1;j<=indices.Size();j++)
		{	
			int k = indices[j];
			
			fOutput << setw(20) << left << sumProduction[p][k];
			fOutput << setw(20) << left << sumDestruction[p][k];
			fOutput << setw(20) << left << sumProduction[p][k]-sumDestruction[p][k];

			for(int i=1;i<=stoichiometry->nu[k].nReactions;i++)
			{
				fOutput << setw(20) << left << C[p].p[k][i];
				fOutput << setw(20) << left << C[p].d[k][i];
			}
			
			fOutput << setw(20) << left << C[p].p[k].GetSumElements();
			fOutput << setw(20) << left << C[p].d[k].GetSumElements();
		}
	}

	fOutput.close();
}

void OpenSMOKE_RateOfProductionAnalysis::PrintRateOfProductionAnalyses_Label(ofstream &fOutput, const int firstColumn)
{
	int count = firstColumn;
	for(int j=1;j<=indices.Size();j++)
	{	
		int k = indices[j];

		PrintTagOnGnuplotLabel(20, fOutput, mix->names[k] + "_Pro", count);
		PrintTagOnGnuplotLabel(20, fOutput, mix->names[k] + "_Des", count);
		PrintTagOnGnuplotLabel(20, fOutput, mix->names[k] + "_Net", count);
		
		for(int i=1;i<=stoichiometry->nu[k].nReactions;i++)
		{
			stringstream reaction_index; reaction_index << stoichiometry->nu[k].iReactions[i];
			PrintTagOnGnuplotLabel(20, fOutput, mix->names[k] + "_P_" + reaction_index.str(), count);
			PrintTagOnGnuplotLabel(20, fOutput, mix->names[k] + "_D_" + reaction_index.str(), count);
		}

		PrintTagOnGnuplotLabel(20, fOutput, mix->names[k] + "_Psum", count);
		PrintTagOnGnuplotLabel(20, fOutput, mix->names[k] + "_Dsum", count);
	}
	fOutput << endl;
}

void OpenSMOKE_RateOfProductionAnalysis::PrintIntegralRateOfProductionAnalyses(const string fileName)
{
	int index;
	
	BzzVector v(NC);

	ofstream fOutput;
	openOutputFileAndControl(fOutput, fileName);
	fOutput.setf(ios::scientific);

	// Print label (production)
	fOutput << setw(49) << left << "Prod.";
	for(int j=1;j<=indices.Size();j++)
		fOutput << setw(8) << right << indices[j];
	fOutput << endl;

	// Print Coefficient Matrix (production)
	for(index=1;index<=NR;index++)
	{
		fOutput << setw(8) << left << index;
		fOutput << setw(40) << left << mix->kinetics.reactionRates->strReaction[index];
		fOutput << " ";

		v = Matrix_Ip.GetRow(index);

		for(int j=1;j<=indices.Size();j++)
			fOutput << setw(8) << fixed << setprecision(3) << right << v[indices[j]]*100.;

		fOutput << endl;
	}
	fOutput << endl << endl;

	
	fOutput << setw(49) << left << "Cons.";
	for(int j=1;j<=indices.Size();j++)
		fOutput << setw(8) << right << indices[j];
	fOutput << endl;
	for(index=1;index<=NR;index++)
	{
		fOutput << setw(8) << left << index;
		fOutput << setw(40) << left << mix->kinetics.reactionRates->strReaction[index];
		fOutput << " ";

		v = Matrix_Id.GetRow(index);

		for(int j=1;j<=indices.Size();j++)
			fOutput << setw(8) << fixed << setprecision(3) << right << v[indices[j]]*100.;

		fOutput << endl;
	}

	// Unimportant reactions
	{
		stringstream string_out_1pc;
		stringstream string_out_01pc;
		UnimportantReactions(string_out_1pc,  0.01);
		UnimportantReactions(string_out_01pc, 0.001);

		fOutput << string_out_1pc.str();
		fOutput << string_out_01pc.str();
	}

	fOutput << endl << endl;
	fOutput << "Most Important reactions" << endl;
	for(int index=1;index<=NC;index++)
	{
		BzzVector		coefficients_p = I.p[index];
		BzzVector		coefficients_d = I.d[index];
		BzzVectorInt	reactions_p(coefficients_p.Size());
		BzzVectorInt	reactions_d(coefficients_d.Size());
	
		Sort(&coefficients_p, &reactions_p);
		Sort(&coefficients_d, &reactions_d);
	
		fOutput << setw(5)  << left << index;
		fOutput << setw(15) << left << mix->names[index];
		fOutput << " (P)  ";
		for(int i=1;i<=min(coefficients_p.Size(),5);i++)
		{
			stringstream number; number << "("; number << setw(5) << fixed << setprecision(3) << coefficients_p[coefficients_p.Size()+1-i]*100.; number << ")";
			fOutput << setw(6) << right << stoichiometry->nu[index].iReactions[reactions_p[coefficients_p.Size()+1-i]];
			fOutput << setw(8) << left  << number.str();
		}
		fOutput << endl;

		fOutput << setw(5)  << left << index;
		fOutput << setw(15) << left << mix->names[index];
		fOutput << " (C)  ";
		for(int i=1;i<=min(coefficients_d.Size(),5);i++)
		{
			stringstream number; number << "("; number << setw(5) << fixed << setprecision(3) << coefficients_d[coefficients_d.Size()+1-i]*100.; number << ")";
			fOutput << setw(6) << right << stoichiometry->nu[index].iReactions[reactions_d[coefficients_d.Size()+1-i]];
			fOutput << setw(8) << left  << number.str();
		}
		fOutput << endl;
	}

	fOutput.close();
}

void OpenSMOKE_RateOfProductionAnalysis::MostImportantReactions(const int index, vector<double> &p, vector<double> &d, vector<int> &ip, vector<int> &id, 
																vector<string> &names_p, vector<string> &names_d)
{
	BzzVector		coefficients_p = I.p[index];
	BzzVector		coefficients_d = I.d[index];
	BzzVectorInt	reactions_p(coefficients_p.Size());
	BzzVectorInt	reactions_d(coefficients_d.Size());

	Sort(&coefficients_p, &reactions_p);
	Sort(&coefficients_d, &reactions_d);

	// Production
	p.resize(0);
	ip.resize(0);
	names_p.resize(0);
	for(int i=1;i<=min(coefficients_p.Size(),20);i++)
	{
//		p.push_back(coefficients_p[coefficients_p.Size()+1-i]*100.);
		p.push_back(coefficients_p[coefficients_p.Size()+1-i]*1.);		// no percentage
		ip.push_back(stoichiometry->nu[index].iReactions[reactions_p[coefficients_p.Size()+1-i]]);
//		names_p.push_back(post_processor->reactions[ip[ip.size()-1]]);
		names_p.push_back(mix->kinetics.reactionRates->strReaction[ip[ip.size()-1]]);
	}

	// Consumption
	d.resize(0);
	id.resize(0);
	names_d.resize(0);
	for(int i=1;i<=min(coefficients_d.Size(),20);i++)
	{
//		d.push_back(coefficients_d[coefficients_d.Size()+1-i]*100.);		
		d.push_back(coefficients_d[coefficients_d.Size()+1-i]*1.);		// no percentage
		id.push_back(stoichiometry->nu[index].iReactions[reactions_d[coefficients_d.Size()+1-i]]);
//		names_d.push_back(post_processor->reactions[id[id.size()-1]]);
		names_d.push_back(mix->kinetics.reactionRates->strReaction[id[id.size()-1]]);
	}
}

void OpenSMOKE_RateOfProductionAnalysis::MostImportantReactions(const int index, vector<double> &t, vector<int> &it, vector<string> &names_t)
{
	BzzVector		coefficients_p = I.p[index];
	BzzVector		coefficients_d = I.d[index];
	BzzVector		coefficients_tot(coefficients_p.Size() + coefficients_d.Size());
	BzzVectorInt	reactions_tot(coefficients_p.Size() + coefficients_d.Size());
	
	int k=1;
	for(int i=1;i<=coefficients_p.Size();i++)
		coefficients_tot[k++] = coefficients_p[i];

	for(int i=1;i<=coefficients_d.Size();i++)
		coefficients_tot[k++] = coefficients_d[i];
	
	Sort(&coefficients_tot, &reactions_tot);

	t.resize(0);
	it.resize(0);
	names_t.resize(0);
	for(int i=1;i<=min(coefficients_tot.Size(),20);i++)
	{
		if (reactions_tot[coefficients_tot.Size()+1-i] <= coefficients_p.Size())
		{
//			t.push_back(coefficients_tot[coefficients_tot.Size()+1-i]*100.);
			t.push_back(coefficients_tot[coefficients_tot.Size()+1-i]*1.);		// no percentage
			it.push_back(stoichiometry->nu[index].iReactions[reactions_tot[coefficients_tot.Size()+1-i]]);
//			names_t.push_back(post_processor->reactions[it[it.size()-1]]);
			names_t.push_back(mix->kinetics.reactionRates->strReaction[it[it.size()-1]]);
		}
		else
		{
//			t.push_back(-coefficients_tot[coefficients_tot.Size()+1-i]*100.);
			t.push_back(-coefficients_tot[coefficients_tot.Size()+1-i]*1.);		// no percentage
			it.push_back(stoichiometry->nu[index].iReactions[reactions_tot[coefficients_tot.Size()+1-i]-coefficients_p.Size()]);
//			names_t.push_back(post_processor->reactions[it[it.size()-1]]);
			names_t.push_back(mix->kinetics.reactionRates->strReaction[it[it.size()-1]]);
		}

	}
}

void OpenSMOKE_RateOfProductionAnalysis::UnimportantReactions(stringstream &string_out, const double eps_threshold)
{
	BzzVector v(NC);

	string_out << "-------------------------------------------------------" << endl;
	string_out << " Unimportant reactions (threshold " << setprecision(3) << eps_threshold*100. << "%)" << endl;
	string_out << "-------------------------------------------------------" << endl;
	string_out << setw(8)  << left << "Index";
	string_out << setw(17) << left << "Prod.";
	string_out << setw(17) << left << "Cons.";
	string_out << setw(2)  << left << "";
	string_out << setw(8)  << left << "Name";
	string_out << endl;
	string_out << "-------------------------------------------------------" << endl;

	for(int index=1;index<=NR;index++)
	{
		v = Matrix_Ip.GetRow(index); double maxp = v.Max();
		v = Matrix_Id.GetRow(index); double maxd = v.Max();
		if ( maxp<eps_threshold && maxd<eps_threshold )
		{	
			string_out << setw(8) << left << index;
			string_out << setw(17) << left << maxp*100.;
			string_out << setw(17) << left << maxd*100.;
			string_out << setw(2) << left << "";
			string_out            << left << mix->kinetics.reactionRates->strReaction[index];
			string_out << endl;
		}
	}

	string_out << endl << endl;
}

void OpenSMOKE_RateOfProductionAnalysis::SaveOnBinaryFile(BzzSave &fOutput)
{
	string dummy;
	char name[Constants::NAME_SIZE];

	dummy = "ROPA";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	dummy = "V20100417";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));

	// Indices of species for local analysis
	dummy = "INDICES";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	fOutput << indices;

	// Stoichiometry
	dummy = "NU";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	stoichiometry->SaveToBinaryFile(fOutput);

	// Integral coefficients
	dummy = "I";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));
	I.SaveToBinaryFile(fOutput);

	// Local coefficients
	dummy = "C";
	strcpy(name, dummy.c_str());
	fOutput.fileSave.write((char*) name, sizeof(name));	
	fOutput << N;
	for(int p=1;p<=N;p++)
		C[p].SaveToBinaryFile(fOutput, indices);
}















void OpenSMOKE_FluxAnalysis::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_FluxAnalysis"	<< endl;
    cout << "Object: " << name_object						<< endl;
    cout << "Error:  " << message							<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_FluxAnalysis::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_FluxAnalysis"	<< endl;
    cout << "Object: "		<< name_object					<< endl;
    cout << "Warning:  "	<< message						<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
}


OpenSMOKE_FluxAnalysis::OpenSMOKE_FluxAnalysis()
{
	name_object = "[Not assigned]";
}

void OpenSMOKE_FluxAnalysis::Initialize(OpenSMOKE_ReactingGas *mix)
{
	int i,k;

	cout << "FLUX: Start" << endl;

	NC = mix->NumberOfSpecies();
	NR = mix->NumberOfReactions();
	M  = mix->M;
	
	nu = new OpenSMOKE_NuManager[NC+1];
	ChangeDimensions(NC, &sumProduction);
	ChangeDimensions(NC, &sumDestruction);

	// Recognize C composition
	int jC = mix->recognize_element("c");
	int jO = mix->recognize_element("o");
	int jH = mix->recognize_element("h");
	int jN = mix->recognize_element("n");
	element_c = mix->elements.GetRow(jC);
	element_o = mix->elements.GetRow(jO);
	element_h = mix->elements.GetRow(jH);
	element_n = mix->elements.GetRow(jN);

	// Dimensioning
	nu_reactants = new BzzVector[NR+1];
	nu_products = new BzzVector[NR+1];
	index_reactants = new BzzVectorInt[NR+1];
	index_products = new BzzVectorInt[NR+1];

	nC = new BzzMatrix[NR+1];
	for(k=1;k<=NR;k++)
		ChangeDimensions(NC, NC, &nC[k]);

	for(i =1;i<=NC;i++)	
		nu[i].Set(i, mix->names[i]);

	BuildNuMatrix(&mix->kinetics);
	
	for(k=1;k<=NC;k++)
		nu[k].Clean();
}

void OpenSMOKE_FluxAnalysis::BuildNuMatrix(OpenSMOKE_Kinetics *kinetics)
{
	int i,k;

	int*	jD1 = kinetics->jDir1.GetHandle();
	int*	jD2 = kinetics->jDir2.GetHandle();
	int*	jD3 = kinetics->jDir3.GetHandle();
	int*	jD4 = kinetics->jDir4.GetHandle();
	int*	jD5 = kinetics->jDir5.GetHandle();
	double* vD5 = kinetics->valDir5.GetHandle();

	int*	jIT1 = kinetics->jInvTot1.GetHandle();
	int*	jIT2 = kinetics->jInvTot2.GetHandle();
	int*	jIT3 = kinetics->jInvTot3.GetHandle();
	int*	jIT4 = kinetics->jInvTot4.GetHandle();
	int*	jIT5 = kinetics->jInvTot5.GetHandle();
	double* vIT5 = kinetics->valInvTot5.GetHandle();
	
	for(i = 1;i <= NC;i++)
	{
		for(k = 1;k <= kinetics->numDir1[i];k++)
			nu[i].Set(-1., *jD1++);
		for(k = 1;k <= kinetics->numDir2[i];k++)
			nu[i].Set(-2., *jD2++);
		for(k = 1;k <= kinetics->numDir3[i];k++)
			nu[i].Set(-4., *jD3++);
		for(k = 1;k <= kinetics->numDir4[i];k++)
			nu[i].Set(-0.50, *jD4++);
		for(k = 1;k <= kinetics->numDir5[i];k++)
			nu[i].Set(-(*vD5++), *jD5++);

		for(k = 1;k <= kinetics->numInvTot1[i];k++)
			nu[i].Set(1., *jIT1++);
		for(k = 1;k <= kinetics->numInvTot2[i];k++)
			nu[i].Set(2., *jIT2++);
		for(k = 1;k <= kinetics->numInvTot3[i];k++)
			nu[i].Set(3., *jIT3++);
		for(k = 1;k <= kinetics->numInvTot4[i];k++)
			nu[i].Set(0.50, *jIT4++);
		for(k = 1;k <= kinetics->numInvTot5[i];k++)
			nu[i].Set(*vIT5++, *jIT5++);
	}


	cout << "FLUX: F" << endl;

	//
	for(k=1;k<=NC;k++)
	{
		for(i=1;i<=nu[k].nReactions;i++)
		{
			if (nu[k].nuReactions[i]<=0.)	
			{
			//	cout << "Reactants: " << k << " " << nu[k].iReactions[i] << " " << nu[k].nuReactions[i] << endl;
				nu_reactants[nu[k].iReactions[i]].Append(-nu[k].nuReactions[i]);
				index_reactants[nu[k].iReactions[i]].Append(k);
			}
			else
			{
			//	cout << "Products: " <<  k << " " << nu[k].iReactions[i] << " " << nu[k].nuReactions[i] << endl;
				nu_products[nu[k].iReactions[i]].Append(nu[k].nuReactions[i]);
				index_products[nu[k].iReactions[i]].Append(k);
			}
		}
	//	getchar();
	}

	//
	for(i=1;i<=NR;i++)
	{
		int j;
		int nReactants = nu_reactants[i].Size();
		int nProducts = nu_products[i].Size();

		cout << "R" << i << "   " << nReactants << " " << nProducts << endl;
		//getchar();
		double sum_c_reactants = 0.;
		for(j=1;j<=nReactants;j++)
			sum_c_reactants += element_c[index_reactants[i][j]]*nu_reactants[i][j];
		if (sum_c_reactants == 0.)	continue;

		int n_species_c_reactants = 0;
		int n_species_c_products = 0;
		BzzVector c_reactants(nReactants);
		BzzVector c_products(nProducts);
		for(j=1;j<=nReactants;j++)	{c_reactants[j] = element_c[index_reactants[i][j]]; if (c_reactants[j]!=0)	n_species_c_reactants++;}
		for(j=1;j<=nProducts;j++)	{c_products[j] = element_c[index_products[i][j]]; if (c_products[j]!=0)		n_species_c_products++;}

		if (n_species_c_reactants == 1 && n_species_c_products == 1)
		{
			int jReactant = 0;
			int jProduct = 0;
			for(j=1;j<=nReactants;j++)
				if (c_reactants[j]!=0)	{	jReactant = j; break; }
			for(j=1;j<=nProducts;j++)
				if (c_products[j]!=0)	{	jProduct = j; break; }
				
			nC[i][index_reactants[i][jReactant]][index_products[i][jProduct]] =  c_reactants[jReactant]*nu_reactants[i][jReactant];
			nC[i][index_products[i][jProduct]][index_reactants[i][jReactant]] = -c_reactants[jReactant]*nu_reactants[i][jReactant];

			cout << index_reactants[i][jReactant] << " " << index_products[i][jProduct] << " " << c_reactants[jReactant]*nu_reactants[i][jReactant] << endl;
		//	getchar();
		}

		else if (n_species_c_reactants == 1 && n_species_c_products == 2)
		{
			int jReactant = 0;
			BzzVectorInt jProduct;
			for(j=1;j<=nReactants;j++)
				if (c_reactants[j]!=0)	{	jReactant = j; break; }
			for(j=1;j<=nProducts;j++)
				if (c_products[j]!=0)	{	jProduct.Append(j); }
				
			nC[i][index_reactants[i][jReactant]][index_products[i][jProduct[1]]] =  c_products[jProduct[1]]*nu_products[i][jProduct[1]];
			nC[i][index_reactants[i][jReactant]][index_products[i][jProduct[2]]] =  c_products[jProduct[2]]*nu_products[i][jProduct[2]];

			nC[i][index_products[i][jProduct[1]]][index_reactants[i][jReactant]] = -c_products[jProduct[1]]*nu_products[i][jProduct[1]];
			nC[i][index_products[i][jProduct[2]]][index_reactants[i][jReactant]] = -c_products[jProduct[2]]*nu_products[i][jProduct[2]];

			cout << index_reactants[i][jReactant] << " " << index_products[i][jProduct[1]] << " " << c_products[jProduct[1]]*nu_products[i][jProduct[1]] << endl;
			cout << index_reactants[i][jReactant] << " " << index_products[i][jProduct[2]] << " " << c_products[jProduct[2]]*nu_products[i][jProduct[2]] << endl;

		//	getchar();
		}

		else if (n_species_c_reactants == 2 && n_species_c_products == 1)
		{
			BzzVectorInt jReactant;
			int jProduct = 0;
			for(j=1;j<=nReactants;j++)
				if (c_reactants[j]!=0)	{	jReactant.Append(j); }
			for(j=1;j<=nProducts;j++)
				if (c_products[j]!=0)	{	jProduct=j; break; }
				
			nC[i][index_reactants[i][jReactant[1]]][index_products[i][jProduct]] =  c_reactants[jReactant[1]]*nu_reactants[i][jReactant[1]];
			nC[i][index_reactants[i][jReactant[2]]][index_products[i][jProduct]] =  c_reactants[jReactant[2]]*nu_reactants[i][jReactant[2]];

			nC[i][index_products[i][jProduct]][index_reactants[i][jReactant[1]]] = -c_reactants[jReactant[1]]*nu_reactants[i][jReactant[1]];
			nC[i][index_products[i][jProduct]][index_reactants[i][jReactant[2]]] = -c_reactants[jReactant[2]]*nu_reactants[i][jReactant[2]];

			cout << index_reactants[i][jReactant[1]] << " " << index_products[i][jProduct] << " " << c_reactants[jReactant[1]]*nu_reactants[i][jReactant[1]] << endl;
			cout << index_reactants[i][jReactant[2]] << " " << index_products[i][jProduct] << " " << c_reactants[jReactant[2]]*nu_reactants[i][jReactant[2]] << endl;

		//	getchar();
		}


		else if (n_species_c_reactants == 2 && n_species_c_products == 2)
		{
			bool iAmbigous = false;

			BzzVectorInt jReactant;
			BzzVectorInt jProduct;
			BzzVector mwReactant;
			BzzVector mwProduct;
			for(j=1;j<=nReactants;j++)
				if (c_reactants[j]!=0)	{	jReactant.Append(j); mwReactant.Append(M(index_reactants[i][j]));	}
			for(j=1;j<=nProducts;j++)
				if (c_products[j]!=0)	{	jProduct.Append(j); mwProduct.Append(M(index_products[i][j]));		}

			BzzVectorInt index(2);
			BzzVector flux(2);

			// Case A: 1-1 = 1-1
			if (c_reactants[1] == 1. && c_reactants[2] == 1.)
			{
				if	(c_products[1] == 1. && c_products[2]==1.)	
				{
					BzzVector o_reactants(nReactants);
					BzzVector o_products(nProducts);
					BzzVector h_reactants(nReactants);
					BzzVector h_products(nProducts);
					BzzVector n_reactants(nReactants);
					BzzVector n_products(nProducts);

					for(j=1;j<=nReactants;j++)	o_reactants[j] = element_o[index_reactants[i][j]];
					for(j=1;j<=nProducts;j++)	o_products[j] = element_o[index_products[i][j]];
					for(j=1;j<=nReactants;j++)	h_reactants[j] = element_h[index_reactants[i][j]];
					for(j=1;j<=nProducts;j++)	h_products[j] = element_h[index_products[i][j]];
					for(j=1;j<=nReactants;j++)	n_reactants[j] = element_n[index_reactants[i][j]];
					for(j=1;j<=nProducts;j++)	n_products[j] = element_n[index_products[i][j]];

					BzzVector bonds(2);
					bonds(1) =	fabs(c_reactants[1]-c_products[1]) + fabs(o_reactants[1]-o_products[1]) +
								fabs(h_reactants[1]-h_products[1]) + fabs(n_reactants[1]-n_products[1]) +
								fabs(c_reactants[2]-c_products[2]) + fabs(o_reactants[2]-o_products[2]) +
								fabs(h_reactants[2]-h_products[2]) + fabs(n_reactants[2]-n_products[2]);

					bonds(2) =	fabs(c_reactants[1]-c_products[2]) + fabs(o_reactants[1]-o_products[2]) +
								fabs(h_reactants[1]-h_products[2]) + fabs(n_reactants[1]-n_products[2]) +
								fabs(c_reactants[2]-c_products[1]) + fabs(o_reactants[2]-o_products[1]) +
								fabs(h_reactants[2]-h_products[1]) + fabs(n_reactants[2]-n_products[1]);

					if (bonds[1]<=bonds[2])
					{
						index[1]=1; 
						index[2]=2; 
						flux[1]=1.; 
						flux[2]=1.;
					}
					else
					{
						index[1]=2; 
						index[2]=1; 
						flux[1]=1.; 
						flux[2]=1.;
					}

					cout << "bonds " << bonds(1) << " " << bonds(2) << endl;
				
				}
				else	
					iAmbigous = true;
			}
			
			// Case A: 1-X = 1-X
			else if (c_reactants[1] == 1. && c_reactants[2] > 1.)
			{
				if		(c_products[1] == 1. && c_products[2]==c_reactants[2])	{index[1]=1; index[2]=2; flux[1]=1.; flux[2]=c_reactants[2];}
				else if (c_products[1] == c_reactants[2] && c_products[2]==1.)	{index[1]=2; index[2]=1; flux[1]=1.; flux[2]=c_reactants[2];}
				else	iAmbigous = true;			
			}
			else if (c_reactants[1] > 1. && c_reactants[2] == 1.)
			{
				if		(c_products[1] == c_reactants[1] && c_products[2]==1.)	{index[1]=1; index[2]=2; flux[1]=c_reactants[1]; flux[2]=1.;}
				else if (c_products[1] == 1. && c_products[2]==c_reactants[1])	{index[1]=2; index[2]=1; flux[1]=c_reactants[1]; flux[2]=1.;}
				else	iAmbigous = true;
			}						
			else iAmbigous = true;

			if (iAmbigous == false)
			{	
				nC[i][index_reactants[i][jReactant[1]]][index_products[i][jProduct[index[1]]]] =  flux[1]*nu_reactants[i][jReactant[1]];
				nC[i][index_reactants[i][jReactant[2]]][index_products[i][jProduct[index[2]]]] =  flux[2]*nu_reactants[i][jReactant[2]];

				nC[i][index_products[i][jProduct[index[1]]]][index_reactants[i][jReactant[1]]] = -flux[1]*nu_reactants[i][jReactant[1]];
				nC[i][index_products[i][jProduct[index[2]]]][index_reactants[i][jReactant[2]]] = -flux[2]*nu_reactants[i][jReactant[2]];
			}

			if (iAmbigous == false)
			{
				cout << index_reactants[i][jReactant[1]] << " " << index_products[i][jProduct[index[1]]] << " " << flux[1]*nu_reactants[i][jReactant[1]] << endl;
				cout << index_reactants[i][jReactant[2]] << " " << index_products[i][jProduct[index[2]]] << " " << flux[2]*nu_reactants[i][jReactant[2]] << endl;
				getchar();
			}
			else
			{
				cout << "Ambigous" << endl; getchar();
			}
		}

		else
		{
			cout << "What to do?" << endl;
			getchar();
		}
	}

}

//void OpenSMOKE_FluxAnalysis::Run(BzzVector &r)
//{
//
//}