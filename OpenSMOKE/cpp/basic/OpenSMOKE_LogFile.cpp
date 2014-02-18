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
#include "basic/OpenSMOKE_LogFile.h"

void OpenSMOKE_LogFile::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_LogFile"		<< endl;
    cout << "File:   " << file_name				<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LogFile::WarningMessage(const string message)
{
    cout << endl;
    cout << "Class:   OpenSMOKE_LogFile"		<< endl;
    cout << "File:    " << file_name				<< endl;
    cout << "Warning: "	<< message			<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
}

OpenSMOKE_LogFile::OpenSMOKE_LogFile()
{
	file_name = "[not assigned]";
} 

void OpenSMOKE_LogFile::Setup(const string _file_name)
{
	file_name = _file_name;
	openOutputFileAndControl(fLog, file_name);
	OpenSMOKE_logo(fLog, "OpenSMOKE Library");
}

void OpenSMOKE_LogFile::WriteLogMessage(const string log_message)
{
	fLog << log_message << endl;
	fLog << endl;
}
