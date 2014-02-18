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

#if !defined(OpenSMOKE_Flame1D_OpposedFlameManager_H)
#define OpenSMOKE_Flame1D_OpposedFlameManager_H

#include "OpenSMOKE.hpp"

class OpenSMOKE_Flame1D;

class OpenSMOKE_Dictionary_Flame1D_OpposedFlameManager : public OpenSMOKE_Dictionary
{
public:
    OpenSMOKE_Dictionary_Flame1D_OpposedFlameManager();
};

class OpenSMOKE_Flame1D_OpposedFlameManager
{
public:

	OpenSMOKE_Flame1D_OpposedFlameManager();
	void SetName(const string name);
	void SetupFromFile(	OpenSMOKE_Flame1D *_flame, OpenSMOKE_Flame1D_DataManager *_data);
	void Run();
	
private:
	
	OpenSMOKE_Flame1D				*flame;
	OpenSMOKE_Flame1D_DataManager	*data;
	ofstream fLog;
	ofstream fOpposedFlameManager;

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);

	void PrintOnFile(ofstream &fLog, ofstream &fFlame,const int index);

	void Lock();
	void SetFuelVelocity(const vector<string> string_vector);
	void SetFuelTemperature(const vector<string> string_vector);
	void SetOxidizerVelocity(const vector<string> string_vector);
	void SetOxidizerTemperature(const vector<string> string_vector);
	void SetPressure(const vector<string> string_vector);

	double FuelDensity();
	double OxidizerDensity();

	void Update();
	void Solution(const int count);

	BzzVector vFuel;
	BzzVector vOxidizer;
	BzzVector TFuel;
	BzzVector TOxidizer;
	BzzVector Pressures;

	bool iAssignedFuelVelocity;
	bool iAssignedOxidizerVelocity;
	bool iAssignedFuelTemperature;
	bool iAssignedOxidizerTemperature;
	bool iAssignedPressure;

	int N;
};

#endif // !defined(OpenSMOKE_Flame1D_OpposedFlameManager_H)
