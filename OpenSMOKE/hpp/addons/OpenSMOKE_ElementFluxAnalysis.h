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

#ifndef OPENSMOKE_ELEMENTFLUXANALYSIS
#define OPENSMOKE_ELEMENTFLUXANALYSIS

#include "BzzMath.hpp"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_RateOfProductionAnalysis;
class OpenSMOKE_Kinetics;
class OpenSMOKE_RateOfProductionCoefficient;

class OpenSMOKE_ElementFluxAnalysis
{

friend  class OpenSMOKE_Kinetics;

public:

	OpenSMOKE_ElementFluxAnalysis();

	void Initialize(OpenSMOKE_ReactingGas *_mix, const vector<string> names);


private:

	OpenSMOKE_ReactingGas	*mix;

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
	
	void BuildNuMatrix(OpenSMOKE_Kinetics *kinetics);
	
	int NC;
	int NR;
	int NE;
	OpenSMOKE_Nu	*stoichiometry;

	BzzMatrix N;
	BzzVectorInt *Areactions;
	BzzVectorInt *list_of_k;
	BzzVectorInt *list_of_j;
	BzzMatrix    *n_x_n;
	BzzMatrix     Alocal;
	BzzVectorInt  list_of_element_indices;
	BzzVectorInt *IsANode;

	void Calculate(const int e, BzzMatrix &A, BzzVector &qForward, BzzVector &qBackward);
	void PrintOnFile(const int e, BzzMatrix &A, ofstream &fFluxAnalysis);


public:

	void Run(const std::string file_name, BzzVector &qForward, BzzVector &qBackward);
	void RunGlobal(const std::string file_name, BzzVector &x, const double &xA, const double &xB, BzzMatrix &qForward, BzzMatrix &qBackward);

};

class OpenSMOKE_ElementFluxAnalysisManager
{
public:
	
	OpenSMOKE_ElementFluxAnalysisManager();
	void Initialize(const vector<string> _names);

	int n;
	vector<double> xA;
	vector<double> xB;
	vector<string> tags;
	vector<string> file_names;

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};

#endif // OPENSMOKE_ELEMENTFLUXANALYSIS
