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

#ifndef OPENSMOKE_POSTPROCESSOR_PFR
#define OPENSMOKE_POSTPROCESSOR_PFR

#include "addons/OpenSMOKE_PostProcessor_General.h"

class OpenSMOKE_PostProcessor; 

class OpenSMOKE_PostProcessor_PFR : public OpenSMOKE_PostProcessor_General
{
public:

	OpenSMOKE_PostProcessor_PFR(OpenSMOKE_PostProcessor *post_processor_);
	void ReadFromBinaryFile(BzzLoad &fLoad);
	void ExportAvailableXAxis(vector<string> &x_axis);
	void ExportAvailableYAxis(vector<string> &y_axis);
	void ImportSelectedAxis(int x_axis, vector<int> y_axis, BzzMatrix &xAxis, BzzMatrix &yAxis, string &name_x, string &name_y, vector<string> &names_lines);
	void Prepare();

private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_POSTPROCESSOR_PFR

