/***************************************************************************
 *   Copyright (C) 2003-2008 by                                            *
 *   Guido Buzzi-Ferraris, Alessio Frassoldati and Alberto Cuoci		   *						   *
 *                                                                         *
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

#ifndef OPENSMOKE_CSTRNETWORK_MYODESYSTEMOBJECTONECSTR_H
#define OPENSMOKE_CSTRNETWORK_MYODESYSTEMOBJECTONECSTR_H

#include "BzzMath.hpp"

class OpenSMOKE_CSTRNetwork;

class OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR : public BzzOdeSystemObject
{
private:
	
	int numComponents;
	OpenSMOKE_CSTRNetwork *cstrNewtwork;

public:
	
	OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR(OpenSMOKE_CSTRNetwork *cstrN);
	
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &f);
	virtual void ObjectBzzPrint(void);
	virtual void GetJacobian(BzzVector &y,double t,BzzMatrix &JJ);

};

#endif // OPENSMOKE_CSTRNETWORK_MYODESYSTEMOBJECTONECSTR_H
