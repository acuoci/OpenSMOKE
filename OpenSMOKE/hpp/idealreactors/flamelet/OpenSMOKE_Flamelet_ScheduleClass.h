/***************************************************************************
 *   Copyright (C) 2006 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_FLAMELET_SCHEDULECLASS)
#define OPENSMOKE_FLAMELET_SCHEDULECLASS

#include "BzzMath.hpp"

class OpenSMOKE_Flamelet_ScheduleClass  
{

friend class OpenSMOKE_Flamelet;

public:

	OpenSMOKE_Flamelet_ScheduleClass();
	void ReadOperations(const string fileName);
	void SetName(const string name);

protected:

	int				nOperations;
	BzzVectorInt	iOperation;
	BzzVectorInt	iOptionA;
	BzzVector iOptionB;

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // !defined(OPENSMOKE_FLAMELET_SCHEDULECLASS)