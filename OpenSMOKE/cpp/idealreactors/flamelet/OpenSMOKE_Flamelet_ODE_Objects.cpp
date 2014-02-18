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

#include "idealreactors/flamelet/OpenSMOKE_Flamelet.h"

void MyOdeSystem_Flamelet::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlamelet->GetSystemFunctions(x, t, f);
}

void MyOdeSystem_Flamelet::assignFlamelet(OpenSMOKE_Flamelet *flamelet)
{
	ptFlamelet = flamelet;
}



void MyOdeSystem_Flamelet_Enthalpy::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlamelet->GetSystemFunctions_Enthalpy(x, t, f);
}

void MyOdeSystem_Flamelet_Enthalpy::assignFlamelet(OpenSMOKE_Flamelet *flamelet)
{
	ptFlamelet = flamelet;
}

void MyOdeSystem_Flamelet_EnthalpyDefect::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlamelet->GetSystemFunctions_EnthalpyDefect(x, t, f);
}

void MyOdeSystem_Flamelet_EnthalpyDefect::assignFlamelet(OpenSMOKE_Flamelet *flamelet)
{
	ptFlamelet = flamelet;
}

void MyOdeSystem_Flamelet_Soot::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlamelet->GetSystemFunctions_Soot(x, t, f);
}

void MyOdeSystem_Flamelet_Soot::assignFlamelet(OpenSMOKE_Flamelet *flamelet)
{
	ptFlamelet = flamelet;
}

