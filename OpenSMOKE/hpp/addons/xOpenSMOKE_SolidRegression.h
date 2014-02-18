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

#ifndef OPENSMOKE_SOLIDREGRESSION
#define OPENSMOKE_SOLIDREGRESSION

#include "BzzMath.hpp"

class OpenSMOKE_SolidExperiment
{
public:

	string name_of_file;

	double temperature;
	double pressureO2;

	int nPoints;
	BzzVector x;
	BzzVector y;

public:

	void ReadFromFile(const string filename);
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

class OpenSMOKE_SolidRegression
{
private:

	BzzOdeStiff o;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);

	BzzMatrix YY;
	BzzVector yy;

public:
	
	int nCases;
	int numTotal;
	BzzVector initialconditions_x;
	BzzVector initialconditions_y;
	BzzVectorInt    indices;
	BzzVector yMin;
	vector<string> list_of_names_of_files;
	OpenSMOKE_SolidExperiment *experiments;

	void Setup(const string filename);
	void Run();
	void ModelOdeRegression(int model, int ex, BzzVector &b, BzzVector &x, BzzVector &y);

};





#endif // OPENSMOKE_SOLIDREGRESSION