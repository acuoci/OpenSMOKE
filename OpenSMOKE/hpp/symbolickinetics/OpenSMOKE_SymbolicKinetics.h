/***************************************************************************
 *   Copyright (C) 2003-2008 by                                            *
 *   Alberto Cuoci		                                                   *
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

#if !defined(OPENSMOKE_SYMBOLICKINETICS_H)
#define OPENSMOKE_SYMBOLICKINETICS_H

#include "BzzMath.hpp"

class OpenSMOKE_SymbolicKinetics
{
public:

	OpenSMOKE_SymbolicKinetics() {};
	virtual ~OpenSMOKE_SymbolicKinetics() {};

	virtual void assignKineticConstants(BzzVector &k_, BzzVector &uK_, BzzVector &logFcent_, BzzVector &kFallOff_) = 0;
	virtual void giveReactionRates(double cTot, BzzVector &c, BzzVector &R) = 0;
	virtual void giveJacobian(BzzVector &c, BzzMatrix &J) = 0;
	        void giveNumericalJacobian(const char kind, const double eps, BzzVector &c, BzzMatrix &J);
	virtual void SetAccurateJacobian() = 0 ;
	virtual void UnsetAccurateJacobian() = 0;

public:

	int NC;
	int NR;

protected:

	int iAccurateJacobian;
};

#endif // !defined(OPENSMOKE_SYMBOLICKINETICS_H)

/*
 void OpenSMOKE_SymbolicKinetics::giveNumericalJacobian(const char kind, const double eps, BzzVector &c, BzzMatrix &J)
 {
	 // TODO

	 double cTot = c.GetSumElements();
	 double err=0.;
	 err = (pow(max(eps,err),1./3.));
	 double h = 1./double(NC)*err;

	 if (kind == 'C')
	 {
		 BzzVector cPlus(NC);
		 BzzVector cMinus(NC);
		 BzzVector RPlus(NC);
		 BzzVector RMinus(NC);

		 for(int j=1;j<=NC;j++)
		 {
			 cPlus  = c;
			 cMinus = c;

			cPlus[j]  += h;
			cMinus[j] -= h;

			giveReactionRates(cTot, cPlus, RPlus);		// [kmol/m3/s]
			giveReactionRates(cTot, cMinus, RMinus);	// [kmol/m3/s]

			double uC = 1./(cPlus[j] - cMinus[j]);
			for(int k=1;k<=NC;k++)
				J[k][j] = (RPlus[k]-RMinus[k]) * uC;
		 }
	 }
	 else if (kind == 'F')
	 {
		 BzzVector cPlus(NC);
		 BzzVector RPlus(NC);
		 BzzVector R(NC);

		 giveReactionRates(cTot, c, R);		// [kmol/m3/s]

		 for(int j=1;j<=NC;j++)
		 {
			 cPlus  = c;

			cPlus[j]  += h;

			giveReactionRates(cTot, cPlus, RPlus);		// [kmol/m3/s]

			double uC = 1./(cPlus[j] - c[j]);
			for(int k=1;k<=NC;k++)
				J[k][j] = (RPlus[k]-R[k]) * uC;
		 }
	 }

	 else if (kind == 'B')
	 {
		 BzzVector cMinus(NC);
		 BzzVector RMinus(NC);
		 BzzVector R(NC);

		 giveReactionRates(cTot, c, R);		// [kmol/m3/s]

		 for(int j=1;j<=NC;j++)
		 {
			 cMinus  = c;

			cMinus[j]  -= h;

			giveReactionRates(cTot, cMinus, RMinus);		// [kmol/m3/s]

			double uC = 1./(c[j] - cMinus[j]);
			for(int k=1;k<=NC;k++)
				J[k][j] = (R[k]-RMinus[k]) * uC;
		 }
	 }
 }*/