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

#ifndef OPENSMOKE_POSTPROCESSOR_GENERAL
#define OPENSMOKE_POSTPROCESSOR_GENERAL

#include "BzzMath.hpp"

class OpenSMOKE_PostProcessor; 

class OpenSMOKE_PostProcessor_General
{
public:

	OpenSMOKE_PostProcessor_General(OpenSMOKE_PostProcessor *post_processor_);
	void SetName(const std::string name);
	
	virtual void ReadFromBinaryFile(BzzLoad &fLoad);
	virtual void ExportAvailableXAxis(vector<string> &x_available);
	virtual void ExportAvailableYAxis(vector<string> &y_available);
	virtual void ImportSelectedAxis(int x_axis, vector<int> y_axis, BzzMatrix &xAxis, BzzMatrix &yAxis, std::string &name_x, std::string &name_y, vector<string> &names_lines);
	virtual void Prepare();

	std::string ExportMainXAxis();

protected:

	OpenSMOKE_PostProcessor *post_processor;

	int NC;
	int NR;

	BzzVectorInt  indices;

	vector<string> list_of_x_available;
	vector<string> list_of_y_available;
	vector<string> list_of_x_labels;
	vector<string> list_of_y_labels;
	int y_start_species;


public:

	BzzVector     x;
	BzzVector     csi;
	BzzVector     MW;
	BzzVector     T;
	BzzVector     P_Pa;
	BzzVector    *X;
	BzzVector	 *omega;
	BzzVector	 *C;

	BzzVector     Ctot;
	BzzVector     rho;

protected:

	std::string name_object;
	virtual void ErrorMessage(const std::string message);
	virtual void WarningMessage(const std::string message);
};

#endif // OPENSMOKE_POSTPROCESSOR_GENERAL

