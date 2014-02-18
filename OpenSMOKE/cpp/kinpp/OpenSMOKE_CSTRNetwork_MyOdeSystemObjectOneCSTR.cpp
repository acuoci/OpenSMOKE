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

#include "kinpp/OpenSMOKE_CSTRNetwork.h"
#include "kinpp/OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR.h"

OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR::OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR(OpenSMOKE_CSTRNetwork *cstrN) 
{
	cstrNewtwork = cstrN;
	numComponents = cstrNewtwork->Reactions->NumberOfSpecies();
}

void OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR::ObjectBzzPrint(void)
{
	::BzzPrint("\nObject Print for Numerical Jacobian");
	::BzzPrint("\nReactor %d", cstrNewtwork->kReactor);
}

void OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR::GetSystemFunctions(BzzVector &x,double t, BzzVector &f)
{
	cstrNewtwork->GetResiduals(x,f);
	Minus(&f);
}

void OpenSMOKE_CSTRNetwork_MyOdeSystemObjectOneCSTR:: GetJacobian(BzzVector &x,double t,BzzMatrix &JJ)
{
	cstrNewtwork->GetJacobian(x,JJ);
}
