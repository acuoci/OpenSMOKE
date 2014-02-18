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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_TRANSPORTDATA_H)
#define OPENSMOKE_CHEMKININTERPRETER_TRANSPORTDATA_H

#include "BzzMath.hpp"
#include <vector>

class OpenSMOKE_CHEMKINInterpreter_TransportSpecies;

class OpenSMOKE_CHEMKINInterpreter_TransportData  
{
public:
	OpenSMOKE_CHEMKINInterpreter_TransportData();
	virtual ~OpenSMOKE_CHEMKINInterpreter_TransportData();

	void ReadTransportData(const string file_name, ofstream *_fLog);
	void CheckForDuplicates();
	double ReducedDiameter(const int i1, const int i2);
	bool IsActivated();



	void ProcessTransportData(BzzVectorInt &indices);
	void SummaryOnFile(ofstream &fOutput, BzzVectorInt &indices);
	BzzVectorInt GiveMeSpeciesIndices(vector<string> &list);

	vector<string> lines;
	int number_of_lines;
	int total_number_of_species;
	BzzVectorInt indexBlankLines;
	BzzVectorInt indexCommentLines;
	BzzVectorInt indexLines;

	ofstream *fLog;

	OpenSMOKE_CHEMKINInterpreter_TransportSpecies *species;

	bool isActivated;

public:

	BzzVectorInt	shape_factor;
	BzzVector		epsylon_over_kb;
	BzzVector		sigma;
	BzzVector		mu;
	BzzVector		alfa;
	BzzVector		zRot298;

private:
	
	void ErrorMessage(const string message);
	void ErrorMessage(const int iLine, const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_TRANSPORTDATA_H)
