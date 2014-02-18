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

#if !defined(ODEOBJECTS)
#define ODEOBJECTS

#include "BzzMath.hpp"

class OpenSMOKE_Flamelet;

class MyOdeSystem_Flamelet : public BzzOdeSystemObject
{
public:
	void assignFlamelet(OpenSMOKE_Flamelet *flamelet);

	OpenSMOKE_Flamelet *ptFlamelet;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
};

class MyOdeSystem_Flamelet_Soot : public BzzOdeSystemObject
{
public:
	void assignFlamelet(OpenSMOKE_Flamelet *flamelet);

	OpenSMOKE_Flamelet *ptFlamelet;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
};

class MyOdeSystem_Flamelet_Enthalpy : public BzzOdeSystemObject
{
public:
	void assignFlamelet(OpenSMOKE_Flamelet *flamelet);

	OpenSMOKE_Flamelet *ptFlamelet;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
};

class MyOdeSystem_Flamelet_EnthalpyDefect : public BzzOdeSystemObject
{
public:
	void assignFlamelet(OpenSMOKE_Flamelet *flamelet);

	OpenSMOKE_Flamelet *ptFlamelet;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
};

#endif // !defined(ODEOBJECTS)
