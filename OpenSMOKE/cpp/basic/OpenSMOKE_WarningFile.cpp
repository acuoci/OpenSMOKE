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

#include "OpenSMOKE.hpp"
#include "basic/OpenSMOKE_WarningFile.h"

void OpenSMOKE_WarningFile::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_WarningFile"		<< endl;
    cout << "File:   " << file_name				<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_WarningFile::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:   OpenSMOKE_WarningFile"		<< endl;
    cout << "File:    " << file_name				<< endl;
    cout << "Warning: "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_WarningFile::OpenSMOKE_WarningFile()
{
	file_name = "[not assigned]";
}

void OpenSMOKE_WarningFile::Setup(const std::string _file_name)
{
	iWarning  = 0;
	file_name = _file_name;
	openOutputFileAndControl(fWarning, file_name);
}

void OpenSMOKE_WarningFile::WriteWarningMessage(const std::string warning_message)
{ 
	iWarning++;

	fWarning << "=====================================================================" << endl;
	fWarning << "                        WARNING #" << iWarning							<< endl;
	fWarning << "=====================================================================" << endl;
	fWarning << warning_message << endl;
	fWarning << "---------------------------------------------------------------------" << endl;
	fWarning << endl;
}
