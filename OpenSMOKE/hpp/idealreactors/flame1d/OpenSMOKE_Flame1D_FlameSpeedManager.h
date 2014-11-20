/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci								   *
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

#if !defined(OpenSMOKE_Flame1D_FlameSpeedManager_H)
#define OpenSMOKE_Flame1D_FlameSpeedManager_H

#include "OpenSMOKE.hpp"

class OpenSMOKE_Flame1D;

class OpenSMOKE_Flame1D_FlameSpeedManager
{
public:

	OpenSMOKE_Flame1D_FlameSpeedManager();
	void SetName(const std::string name);
	void SetupFromFile(OpenSMOKE_Flame1D *_flame, OpenSMOKE_Flame1D_DataManager *_data);
	void Run();
	
	BzzVector GiveMeMoleFractions();
	BzzVector GiveMeMassFractions();
	void NextEquivalenceRatio();

	BzzVector list_of_phi;
	BzzVector list_of_flame_speeds;
	BzzVector list_of_flame_T;

private:
	
	OpenSMOKE_Flame1D				*flame;
	OpenSMOKE_Flame1D_DataManager	*data;

	OpenSMOKE_GasStream	inlet_stream;
	BzzVector equivalence_ratios;
	vector< vector<string> > list_names;
	vector< vector<double> > list_values;
	vector<string> tag_flame;

	bool iEquivalenceRatios;
	bool iMoleFractions;
	bool iMassFractions;
	bool iMoles;
	bool iMasses;

	int tag_unit_equivalence_ratio;
	
	int N;
	ofstream fFlameSpeed;
	ofstream fLog;

	void Print(const double equivalence_ratio);

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

	void PrintComposition(ofstream &fOut, const double equivalence);
};

#endif // !defined(OpenSMOKE_Flame1D_FlameSpeedManager_H)
