/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
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

#if !defined(OpenSMOKE_Flame1D_SCHEDULECLASS)
#define OpenSMOKE_Flame1D_SCHEDULECLASS

#include "OpenSMOKE.hpp"

class OpenSMOKE_Flame1D_ScheduleClass  
{
public:

	OpenSMOKE_Flame1D_ScheduleClass();
	void SetName(const std::string name);

	int nOperations;
	BzzVectorInt iOperation;
	BzzVector iOptionA;
	BzzVector iOptionB;
	std::string *species;

	void readOperations(const std::string fileName);

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};

#endif // !defined(SCHEDULECLASS)
