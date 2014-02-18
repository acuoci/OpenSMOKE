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

#ifndef OpenSMOKE_Flame1D_DAE_CLASSES
#define OpenSMOKE_Flame1D_DAE_CLASSES

#include "OpenSMOKE.hpp"

class OpenSMOKE_Flame1D;

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									DAE SYSTEMS													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

class OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOMOMENTUM : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Opposed_ONLYMOMENTUM : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Opposed_NOENERGY : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Opposed_QMOM_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Opposed_SOOT_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Premixed_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Premixed_FLAMESPEED : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Premixed_NOENERGY : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Premixed_QMOM_NOENERGY : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Premixed_SOOT_NOENERGY : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEMS													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

class OpenSMOKE_Flame1D_MyOdeSystem_Premixed_NOENERGY : public BzzOdeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyOdeSystem_Premixed_ALL : public BzzOdeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									NLS SYSTEMS													   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////


class OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_ALL : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_NOENERGY : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Premixed_FLAMESPEED : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ALL : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Opposed_ALL : public BzzMyNonLinearSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYUGH : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYMASSFRACTIONS : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_COLDREDUCED : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_ONLYT : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Opposed_NOENERGY : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};





class OpenSMOKE_Flame1D_MyDaeSystem_Twin_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOMOMENTUM : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Twin_NOENERGY : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Twin_QMOM_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyDaeSystem_Twin_SOOT_ALL : public BzzDaeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};


class OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ALL : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Reduced_Twin_ALL : public BzzMyNonLinearSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYUGH : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYMASSFRACTIONS : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_COLDREDUCED : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_ONLYT : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_Flame1D_MyNonLinearSystem_Twin_NOENERGY : public BzzMyNonLinearSystemSparseObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetResiduals(BzzVector &y, BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};


class OpenSMOKE_Flame1D_MyOdeSystem_SingleReactor_Isothermal : public BzzOdeSystemObject
{
public:
	void assignFlame(OpenSMOKE_Flame1D *flame);

	OpenSMOKE_Flame1D *ptFlame;
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void ObjectBzzPrint(void);
};


#endif
