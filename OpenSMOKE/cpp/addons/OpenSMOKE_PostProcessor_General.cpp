/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci         	                       *
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
#include "basic/OpenSMOKE_Utilities.h"
#include "basic/OpenSMOKE_Constants.h"
#include "addons/OpenSMOKE_PostProcessor.h"
#include "addons/OpenSMOKE_PostProcessor_General.h"

void OpenSMOKE_PostProcessor_General::ErrorMessage(const string message)
{
}

void OpenSMOKE_PostProcessor_General::WarningMessage(const string message)
{
}

void OpenSMOKE_PostProcessor_General::SetName(const string name)
{
	name_object = name;
}

OpenSMOKE_PostProcessor_General::OpenSMOKE_PostProcessor_General(OpenSMOKE_PostProcessor *post_processor_)
{
	name_object = "name not assigned";
	post_processor = post_processor_;
	NC = post_processor->NC;
	NR = post_processor->NR;
}

void OpenSMOKE_PostProcessor_General::ReadFromBinaryFile(BzzLoad &fLoad)
{
}

void OpenSMOKE_PostProcessor_General::Prepare()
{
}

void OpenSMOKE_PostProcessor_General::ExportAvailableXAxis(vector<string> &x_axis)
{
}	

void OpenSMOKE_PostProcessor_General::ExportAvailableYAxis(vector<string> &y_axis)
{
}	

void OpenSMOKE_PostProcessor_General::ImportSelectedAxis(int x_axis, vector<int> y_axis, BzzMatrix &xAxis, BzzMatrix &yAxis, string &name_x, string &name_y, vector<string> &names_lines)
{
}

string OpenSMOKE_PostProcessor_General::ExportMainXAxis()
{
	return list_of_x_labels.at(0);
}