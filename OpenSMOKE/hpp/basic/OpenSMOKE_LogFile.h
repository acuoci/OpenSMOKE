/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci   	                       *
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

#ifndef OPENSMOKE_LOGFILE
#define OPENSMOKE_LOGFILE

#include <vector>
#include "BzzMath.hpp"

class OpenSMOKE_LogFile
{
public:

	OpenSMOKE_LogFile();

	void Setup(const string file_name);
	void WriteLogMessage(const string log_message);

private:

	string		file_name;
	ofstream	fLog;

	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif	// OPENSMOKE_LOGFILE
