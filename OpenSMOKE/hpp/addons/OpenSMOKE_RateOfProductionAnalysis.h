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

#ifndef OPENSMOKE_RATEOFPRODUCTIONANALYSIS
#define OPENSMOKE_RATEOFPRODUCTIONANALYSIS

#include "BzzMath.hpp"

class OpenSMOKE_ReactingGas;
class OpenSMOKE_RateOfProductionAnalysis;
class OpenSMOKE_Kinetics;
class OpenSMOKE_RateOfProductionCoefficient;

class OpenSMOKE_RateOfProductionAnalysis
{

friend  class OpenSMOKE_Kinetics;

public:

	OpenSMOKE_RateOfProductionAnalysis();
	void SetName(const std::string name);

	void Initialize(OpenSMOKE_ReactingGas *mix, BzzVectorInt &_indices);
	void SetNumberOfPoints(const int N_);
	void Run(BzzMatrix &r, BzzVector &x);
	void Run(BzzMatrix &r, BzzVector &x, const double local_coordinate);
	void Run(BzzMatrix &r, BzzVector &x, const double xA, const double xB);
	
	// Print ASCII Files
	void PrintRateOfProductionAnalyses(const std::string fileName, BzzVector  &x, BzzVector &T);
	void PrintIntegralRateOfProductionAnalyses(const std::string fileName);

	// Print Binary Files
	void SaveOnBinaryFile(BzzSave &fOutput);

	void MostImportantReactions(const int index, vector<double> &p, vector<double> &d, vector<int> &ip, vector<int> &id, 
								vector<string> &names_p, vector<string> &names_d);
	void MostImportantReactions(const int index, vector<double> &t, vector<int> &it, vector<string> &names_t);

	void UnimportantReactions(stringstream &stringout, const double eps_threshold);

private:

	int N;			// Number of points
	int NC;			// Number of species
	int NR;			// Number of reactions

	OpenSMOKE_RateOfProductionCoefficient	*C;		// Local coefficients
	OpenSMOKE_RateOfProductionCoefficient	 I;		// Integral coefficients

	BzzMatrix		sumProduction;
	BzzMatrix		sumDestruction;
	BzzVector		sumIntegralProduction;
	BzzVector		sumIntegralDestruction;

	BzzVectorInt	indices;

	BzzMatrixSparse Matrix_Ip;			// Sparse Matrix (production)
	BzzMatrixSparse Matrix_Id;			// Sparse Matrix (consumption)

	OpenSMOKE_ReactingGas	*mix;				// Pointer to gas mixture
	OpenSMOKE_Nu			*stoichiometry;		// Pointer to kinetic stoichiometry

	void IntegralRateOfProductionAnalyses();			// Integral Rate Of Production Analysis	
	void PrintRateOfProductionAnalyses_Label(ofstream &fOutput, const int firstColumn);

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

};





class OpenSMOKE_FluxAnalysis
{

friend  class OpenSMOKE_Kinetics;

public:

	OpenSMOKE_FluxAnalysis();

	void Initialize(OpenSMOKE_ReactingGas *mix);


private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
	
	void BuildNuMatrix(OpenSMOKE_Kinetics *kinetics);
	
	int NC;
	int NR;
	BzzVector M;
	OpenSMOKE_NuManager	*nu;
	BzzVector			 sumProduction;
	BzzVector			 sumDestruction;
	BzzVector			 element_c;
	BzzVector			 element_h;
	BzzVector			 element_o;
	BzzVector			 element_n;

public:

	BzzMatrix *nC;
	BzzVector *nu_reactants;
	BzzVector *nu_products;
	BzzVectorInt *index_reactants;
	BzzVectorInt *index_products;
};

#endif // OPENSMOKE_RATEOFPRODUCTIONANALYSIS
