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
#include "basic/OpenSMOKE_Conversions.h"
#include "basic/OpenSMOKE_Utilities.h"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "addons/OpenSMOKE_ElementFluxAnalysis.h"


void OpenSMOKE_ElementFluxAnalysis::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ElementFluxAnalysis"	<< endl;
    cout << "Object: " << name_object						<< endl;
    cout << "Error:  " << message							<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_ElementFluxAnalysis::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ElementFluxAnalysis"	<< endl;
    cout << "Object: "		<< name_object					<< endl;
    cout << "Warning:  "	<< message						<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
}


OpenSMOKE_ElementFluxAnalysis::OpenSMOKE_ElementFluxAnalysis()
{
	name_object = "[Not assigned]";
}

void OpenSMOKE_ElementFluxAnalysis::Initialize(OpenSMOKE_ReactingGas *_mix, const vector<string> names)
{
	cout << "FLUX: Start" << endl;

	mix	= _mix;

	NE = mix->NumberOfElements();
	NC = mix->NumberOfSpecies();
	NR = mix->NumberOfReactions();

	stoichiometry = new OpenSMOKE_Nu(mix);

	// Recognize element list
	cout << "Recognize elements" << endl;
	for(unsigned int k=1;k<=names.size();k++)
		list_of_element_indices.Append(mix->recognize_element(names[k-1]));


	ChangeDimensions(list_of_element_indices.Size(), NR, &N);

	for(int k=1;k<=list_of_element_indices.Size();k++)
	{	
		for(int j=1;j<=NC;j++)
		{
			double n_elements = mix->elements[list_of_element_indices[k]][j];
			for(int i=1;i<=stoichiometry->nu[j].nReactions;i++)
			{
					double coefficient = stoichiometry->nu[j].nuReactions[i];
					if (coefficient>0)
						N[k][stoichiometry->nu[j].iReactions[i]] += n_elements*coefficient;
			}
		}
	}

	// Create list of reactions for each species couple
	{
		cout << "Create list of reactions for each species couple" << endl;
		Areactions = new BzzVectorInt[NC*NC+1];
		for(int k=1;k<=NC;k++)
		{
			for(int j=1;j<=NC;j++)
				for(int kk=1;kk<=stoichiometry->nu[k].nReactions;kk++)
					for(int jj=1;jj<=stoichiometry->nu[j].nReactions;jj++)
					{
						if (stoichiometry->nu[k].iReactions[kk] == stoichiometry->nu[j].iReactions[jj])
						{
							if (stoichiometry->nu[k].nuReactions[kk]*stoichiometry->nu[j].nuReactions[jj] < 0.)
							{
								if (stoichiometry->nu[k].nuReactions[kk] < 0.)
									Areactions[(k-1)*NC+j].Append(stoichiometry->nu[k].iReactions[kk]);
								else
									Areactions[(k-1)*NC+j].Append(-stoichiometry->nu[k].iReactions[kk]);
								break;
							}
						}
					}
		}
	}

	// Create matrix of species element products
	{
		cout << "Create matrix of species element products" << endl;

		n_x_n = new BzzMatrix[list_of_element_indices.Size()+1];
		for(int e=1;e<=list_of_element_indices.Size();e++)
			ChangeDimensions(NC,NC, &n_x_n[e]);

		for(int e=1;e<=list_of_element_indices.Size();e++)
			for(int k=1;k<=NC;k++)
				for(int j=1;j<=NC;j++)
					n_x_n[e][k][j] = mix->elements[list_of_element_indices[e]][k] * 
									 mix->elements[list_of_element_indices[e]][j] ;
	}

	// Create list of species indices
	{
		cout << "Create list of species indices" << endl;
		list_of_k = new BzzVectorInt[list_of_element_indices.Size()+1];
		list_of_j = new BzzVectorInt[list_of_element_indices.Size()+1];
		IsANode   = new BzzVectorInt[list_of_element_indices.Size()+1];
		for (int e=1;e<=list_of_element_indices.Size();e++)
		{
			int count=0;
			for(int k=1;k<=NC;k++)
				for(int j=1;j<=NC;j++)
					if (n_x_n[e][k][j] > 0.)
						if (Areactions[(k-1)*NC+j].Size()>0)
							count++;

			ChangeDimensions(count, &list_of_k[e]);
			ChangeDimensions(count, &list_of_j[e]);

			count = 1;
			for(int k=1;k<=NC;k++)
				for(int j=1;j<=NC;j++)
					if (n_x_n[e][k][j] > 0.)
						if (Areactions[(k-1)*NC+j].Size()>0)
						{
							list_of_k[e][count]   = k;
							list_of_j[e][count++] = j;
						}
			
			for(int k=1;k<=NC;k++)
				for(int j=1;j<=list_of_k[e].Size();j++)
					if (list_of_k[e][j]==k || list_of_j[e][j]==k)
					{	
						IsANode[e].Append(k);
						break;
					}

		}
	}

	// Allocate memory
	ChangeDimensions(NC,NC, &Alocal);

	bool iDebug = true;
	if (iDebug==true)
	{
		list_of_element_indices.BzzPrint("list_of_element_indices");

		for (int e=1;e<=list_of_element_indices.Size();e++)
		{
			list_of_k[e].BzzPrint("list_of_k-%d",e);
			list_of_j[e].BzzPrint("list_of_j-%d",e);
			IsANode[e].BzzPrint("IsANode-%d",e);
		}

		for(int k=1;k<=NC;k++)
			for(int j=1;j<=NC;j++)
				if (Areactions[(k-1)*NC+j].Size()>0) 
					Areactions[(k-1)*NC+j].BzzPrint("%d-%d", k, j);

		N.BzzPrint("N");
		
		for(int k=1;k<=list_of_element_indices.Size();k++)
			n_x_n[k].BzzPrint("%d", k);

		

	//	exit(-1);
	}
}

void OpenSMOKE_ElementFluxAnalysis::Calculate(const int e, BzzMatrix &A, BzzVector &qForward, BzzVector &qBackward)
{
	for(int k=1;k<=NC;k++)
		for(int j=1;j<=NC;j++)
		{
			if (n_x_n[e][k][j] > 0.)
			{
				int index = (k-1)*NC+j;
				if (Areactions[index].Size()>0)
				{
					for(int i=1;i<=Areactions[index].Size();i++)
					{
							if (Areactions[index][i]>0)
								A[k][j]  += qForward[Areactions[index][i]]  / N[e][Areactions[index][i]];
							else
								A[k][j]  += qBackward[-Areactions[index][i]] / N[e][-Areactions[index][i]];
					}
					A[k][j]  *= n_x_n[e][k][j];
				}
			}
		}
}

void OpenSMOKE_ElementFluxAnalysis::PrintOnFile(const int e, BzzMatrix &A, ofstream &fFluxAnalysis)
{
	double small_tag = 1e-10;

	for(int k=1;k<=IsANode[e].Size()-1;k++)
		fFluxAnalysis << mix->names[IsANode[e][k]] << " ";
	fFluxAnalysis << mix->names[IsANode[e][IsANode[e].Size()]] << endl;

	for(int k=1;k<=NC;k++)
		for(int j=k+1;j<=NC;j++)
			if (n_x_n[e][k][j] > 0.)
				if (Areactions[(k-1)*NC+j].Size()>0)
				{
					fFluxAnalysis << mix->names[k]	<< " " << mix->names[j] << " ";
					
					if(A[k][j]>small_tag) 
						fFluxAnalysis <<  A[k][j]	<< " ";
					else
						fFluxAnalysis <<  0.		<< " ";
					
					if(A[j][k]>small_tag) 
						fFluxAnalysis << -A[j][k]	<< endl;
					else
						fFluxAnalysis <<  0.		<< endl;
				}
}


void OpenSMOKE_ElementFluxAnalysis::Run(const string file_name, BzzVector &qForward, BzzVector &qBackward)
{
	for (int e=1;e<=list_of_element_indices.Size();e++)
	{
		Alocal = 0.;
		
		Calculate(e, Alocal, qForward, qBackward);

		string name;
		name = file_name + "_" + mix->list_of_elements[list_of_element_indices[e]-1] + ".ced";
			
		ofstream fFluxAnalysis;
		openOutputFileAndControl(fFluxAnalysis, name);

		fFluxAnalysis << "OpenSMOKE Flux Analysis - " << mix->list_of_elements[list_of_element_indices[e]-1] << " Element" << endl;
	
		PrintOnFile(e, Alocal, fFluxAnalysis);
	}
}

void OpenSMOKE_ElementFluxAnalysis::RunGlobal(const string file_name, BzzVector &x, const double &xA, const double &xB, BzzMatrix &qForward, BzzMatrix &qBackward)
{
	int jA=0;
	int jB=0;

	for(int j=2;j<=x.Size();j++)
		if (x[j]>=xA)	
		{
			jA = j-1;
			break;
		}

	for(int j=x.Size()-1;j>=1;j--)
		if (x[j]<=xB)	
		{
			jB = j+1;
			break;
		}

	if (jA*jB == 0)
		ErrorMessage("Wrong user defined interval for Global Element Flux Analysis...");

	BzzVector qForwardTotal(NR);
	BzzVector qBackwardTotal(NR);
	if (jA != jB)
	{
		BzzVector dx(x.Size()-1);
		for(int j=1;j<=x.Size()-1;j++)
			dx[j] = x[j+1]-x[j];

		for(int j=jA;j<=jB-1;j++)
			for(int k=1;k<=NR;k++)
			{
				qForwardTotal[k]  = 0.50*dx[j]*(qForward[j][k]  + qForward[j+1][k]);
				qBackwardTotal[k] = 0.50*dx[j]*(qBackward[j][k] + qBackward[j+1][k]);
			}
	}
	else
	{
		qForwardTotal = qForward.GetRow(jA);
		qBackwardTotal = qBackward.GetRow(jA);
	}

	Run(file_name, qForwardTotal, qBackwardTotal);
}

void OpenSMOKE_ElementFluxAnalysisManager::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ElementFluxAnalysis"	<< endl;
    cout << "Object: " << name_object						<< endl;
    cout << "Error:  " << message							<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_ElementFluxAnalysisManager::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_ElementFluxAnalysisManager"	<< endl;
    cout << "Object: "		<< name_object					<< endl;
    cout << "Warning:  "	<< message						<< endl;
    cout << "Press a key to continue... "					<< endl;
    getchar();
}


OpenSMOKE_ElementFluxAnalysisManager::OpenSMOKE_ElementFluxAnalysisManager()
{
	name_object = "[Not assigned]";
}

void OpenSMOKE_ElementFluxAnalysisManager::Initialize(const vector<string> _names)
{	
	n = _names.size()/4;
	for(int j=1;j<=n;j++)
	{
		int index = (j-1)*4;
		tags.push_back(_names[index]);
		xA.push_back(OpenSMOKE_Conversions::conversion_length(atof(_names[index+1].c_str()), _names[index+3]));
		xB.push_back(OpenSMOKE_Conversions::conversion_length(atof(_names[index+2].c_str()), _names[index+3]));
		file_names.push_back("FluxAnalysis_" + tags[j-1]);
	}
}