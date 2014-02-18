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

#include "idealreactors/flame1d/OpenSMOKE_Flame1D_NLS_DAE_ODE_Classes.h"
#include "idealreactors/flame1d/OpenSMOKE_Flame1D.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DAE SYSTEMS													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

// 1. Opposed_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Opposed_All(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 2. Opposed_NOMOMENTUM
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOMOMENTUM::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOMOMENTUM::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Opposed_NoMomentum(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOMOMENTUM::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 2. Opposed_ONLYMOMENTUM
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ONLYMOMENTUM::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ONLYMOMENTUM::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Opposed_OnlyMomentum(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ONLYMOMENTUM::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 2. Opposed_NOMOMENTUM
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOENERGY::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Opposed_NoEnergy(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 3. Premixed_SOOT_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_SOOT_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_SOOT_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Opposed_SOOT_ALL(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_SOOT_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 3. Opposed_QMOM_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_QMOM_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_QMOM_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Opposed_QMOM_ALL(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Opposed_QMOM_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}



// 3. Premixed_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Premixed_All(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 3. Premixed_FLAMESPEED
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_FLAMESPEED::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_FLAMESPEED::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Premixed_FlameSpeed(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_FLAMESPEED::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 4. Premixed_NOENERGY
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_NOENERGY::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Premixed_NoEnergy(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 5. Premixed_QMOM_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Premixed_QMOM_All(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 5. Premixed_QMOM_NOENERGY
// ---------------------------------------------------------------------------------
void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_NOENERGY::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Premixed_QMOM_NoEnergy(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 6. Premixed_SOOT_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Premixed_SOOT_ALL(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 7. Premixed_SOOT_NOENERGY
// ---------------------------------------------------------------------------------
void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_NOENERGY::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Premixed_SOOT_NoEnergy(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEMS													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

// 1. Premixed_NOENERGY
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyOdeSystem_Premixed_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyOdeSystem_Premixed_NOENERGY::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->ODESystem_Premixed_NoEnergy(x, t, f);
}

void OpenSMOKE_Flame1D_MyOdeSystem_Premixed_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 2. Premixed_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyOdeSystem_Premixed_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyOdeSystem_Premixed_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->ODESystem_Premixed_All(x, t, f);
}

void OpenSMOKE_Flame1D_MyOdeSystem_Premixed_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									NLS SYSTEMS													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

// 1. Premixed_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_ALL::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Premixed_All(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 2. Premixed_NOENERGY
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_NOENERGY::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Premixed_NoEnergy(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 2. Premixed_FLAMESPEED
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_FLAMESPEED::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_FLAMESPEED::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Premixed_FlameSpeed(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_FLAMESPEED::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 3. Opposed_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ALL::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Opposed_All(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 3. Reduced Opposed_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Opposed_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Opposed_ALL::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Reduced_Opposed_All(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Opposed_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 4. Opposed_ONLYUGH
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYUGH::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYUGH::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Opposed_OnlyUGH(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYUGH::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 5. Opposed_ONLYMASSFRACTIONS
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYMASSFRACTIONS::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYMASSFRACTIONS::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Opposed_OnlyMassFractions(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYMASSFRACTIONS::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 6. Opposed_COLDREDUCED
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_COLDREDUCED::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_COLDREDUCED::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Opposed_ColdReduced(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_COLDREDUCED::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 7. Opposed_ONLYT
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYT::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYT::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Opposed_OnlyT(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYT::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 7. Opposed_NOENERGY
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_NOENERGY::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Opposed_NoEnergy(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}










// 1. Twin_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Twin_All(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 2. Twin_NOMOMENTUM
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOMOMENTUM::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOMOMENTUM::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Twin_NoMomentum(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOMOMENTUM::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 2. Twin_NOMOMENTUM
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOENERGY::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Twin_NoEnergy(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 3. Premixed_SOOT_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_SOOT_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_SOOT_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Twin_SOOT_ALL(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_SOOT_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 3. Twin_QMOM_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_QMOM_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_QMOM_ALL::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->DAESystem_Twin_QMOM_ALL(x, t, f);
}

void OpenSMOKE_Flame1D_MyDaeSystem_Twin_QMOM_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 3. Twin_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ALL::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Twin_All(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 3bis. Twin_Reduced_ALL
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Twin_ALL::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Twin_ALL::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Reduced_Twin_All(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Twin_ALL::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 4. Twin_ONLYUGH
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYUGH::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYUGH::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Twin_OnlyUGH(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYUGH::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 5. Twin_ONLYMASSFRACTIONS
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYMASSFRACTIONS::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYMASSFRACTIONS::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Twin_OnlyMassFractions(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYMASSFRACTIONS::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}


// 6. Twin_COLDREDUCED
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_COLDREDUCED::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_COLDREDUCED::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Twin_ColdReduced(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_COLDREDUCED::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 7. Twin_ONLYT
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYT::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYT::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Twin_OnlyT(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYT::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 7. Twin_NOENERGY
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_NOENERGY::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_NOENERGY::GetResiduals(BzzVector &x,BzzVector &f)
{
	ptFlame->nonLinearSystem_Twin_NoEnergy(x, f);
}

void OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_NOENERGY::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

// 1. SingleReactor_Isothermal
// ---------------------------------------------------------------------------------

void OpenSMOKE_Flame1D_MyOdeSystem_SingleReactor_Isothermal::ObjectBzzPrint(void)
{
}

void OpenSMOKE_Flame1D_MyOdeSystem_SingleReactor_Isothermal::GetSystemFunctions(BzzVector &x,double t,BzzVector &f)
{
	ptFlame->ODESystem_SingleReactor_Isothermal(x, t, f);
}

void OpenSMOKE_Flame1D_MyOdeSystem_SingleReactor_Isothermal::assignFlame(OpenSMOKE_Flame1D *flame)
{
	ptFlame = flame;
}

