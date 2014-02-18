/***************************************************************************
 *   Copyright (C) 2010 by Alberto Cuoci								   *
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

#include "droplet/OpenSMOKE_DropletMicrogravity_DAE_Objects.h"
#include "droplet/OpenSMOKE_Droplet.h"
#include "droplet/OpenSMOKE_DropletMicrogravity.h"


void OpenSMOKE_DropletMicrogravity_MyDaeSystemEigenValue::ObjectBzzPrint(void)
{
}

void OpenSMOKE_DropletMicrogravity_MyDaeSystemEigenValue::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptDroplet->DAESystemEigenValue(x, t, f);
}

void OpenSMOKE_DropletMicrogravity_MyDaeSystemEigenValue::assignDroplet(OpenSMOKE_DropletMicrogravity *droplet)
{
	ptDroplet = droplet;
}


void OpenSMOKE_DropletMicrogravity_MyDaeSystemUnsteadyBatch::ObjectBzzPrint(void)
{
}

void OpenSMOKE_DropletMicrogravity_MyDaeSystemUnsteadyBatch::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptDroplet->DAESystemUnsteadyBatch(x, t, f);
}

void OpenSMOKE_DropletMicrogravity_MyDaeSystemUnsteadyBatch::assignDroplet(OpenSMOKE_DropletMicrogravity *droplet)
{
	ptDroplet = droplet;
}


void OpenSMOKE_Droplet_MyDaeSystemNoMomentum::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Droplet_MyDaeSystemNoMomentum::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptDroplet->DAESystemNoMomentum(x, t, f);
}

void OpenSMOKE_Droplet_MyDaeSystemNoMomentum::assignDroplet(OpenSMOKE_Droplet *droplet)
{
	ptDroplet = droplet;
}

void OpenSMOKE_Droplet_MyDaeSystemOnlyTemperature::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Droplet_MyDaeSystemOnlyTemperature::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptDroplet->DAESystemOnlyTemperature(x, t, f);
}

void OpenSMOKE_Droplet_MyDaeSystemOnlyTemperature::assignDroplet(OpenSMOKE_Droplet *droplet)
{
	ptDroplet = droplet;
}