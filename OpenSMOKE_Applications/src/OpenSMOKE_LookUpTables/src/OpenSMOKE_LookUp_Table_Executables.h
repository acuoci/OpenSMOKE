/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_LOOKUP_TABLE_EXECUTABLES)
#define OPENSMOKE_LOOKUP_TABLE_EXECUTABLES

#include "engine/OpenSMOKE_ReactingGas.h"
#include "basic/OpenSMOKE_Dictionary.h"
#include "sources_in_flamelet_library.h"

class OpenSMOKE_Dictionary_LookUp_Table_Executables: public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_LookUp_Table_Executables();
};



class OpenSMOKE_LookUp_Table_Executables
{
public:

	OpenSMOKE_LookUp_Table_Executables();

private:
	
	bool iSootSourceTerms_Adiabatic;
	bool iSootSourceTerms_NonAdiabatic;

	int iChiPDF;								// 0=delta dirac  1=log-normal distribution
	interactions_modes interactions_mode;		

public:

	void DefineFromFile(const string file_name);
	void Run();
	void Exe_SootSourceTerms_Adiabatic();
	void Exe_SootSourceTerms_NonAdiabatic();
	void Exe_SootProfiles_NonAdiabatic();

public:

	void SetSootClosureModel(const string _iSootMode);
	void SetSootSourcesModel(const string _iSourcesModel);
	void SetLogPDF();

private:

	void SootSourceTerms_Adiabatic();
	void SootSourceTerms_NonAdiabatic();
	void SootProfiles_NonAdiabatic();

	void CalculateSourceTerms(	BzzVector &Source, sourceField_flamelet_library *field, 
								BzzMatrix &infoMatrix, int numberOfCells, int nEnthalpyDefects, double defectStep);
	void CalculateSourceTerms(	BzzVector &Source, sourceField_flamelet_library &field, BzzMatrix &infoMatrix, int numberOfCells);


private:

	vector<string> construct_name_list(const string root, const string phenomenon, const int iModel);
	double linear_interpolation_for_enthalpy_defect(double &x, double &sourceA, double &sourceB, double &xA, double &xB);
	double linear_extrapolation_for_enthalpy_defect(double &x, double &sourceA, double &sourceB, double &xA, double &xB);

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_LOOKUP_TABLE_EXECUTABLES)
