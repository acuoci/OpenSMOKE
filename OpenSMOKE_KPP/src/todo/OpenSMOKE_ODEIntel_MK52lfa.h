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

#ifndef OpenSMOKE_ODEIntel_MK52lfa_H
#define OpenSMOKE_ODEIntel_MK52lfa_H

#include "BzzMath.hpp"
#include "intel_ode.h"

void ResidualsQ(int *n_, double *t_, double *y_, double *f_);
void JacobianQ(int *n_, double *t_, double *y_, double *a_);

class OpenSMOKE_ODEIntel_MK52lfa
{

public:
	OpenSMOKE_ODEIntel_MK52lfa(void)
	{
		ipar_ = new int[128];
		for(int j=0;j<128;j++)
			ipar_[j] = 0;
		h_ = 1.e-7;
		hmin_ = 1.e-12;
		ep_ = 1.e-6;
		tr_ = 1.e-3;
		ierr_ = 0;
	}

	~OpenSMOKE_ODEIntel_MK52lfa(void) {};

	void SetResiduals(void (*residuals)(int *n_, double *t_, double *y_, double *f_))
	{
		Residuals = residuals;
	}
	
	void SetJacobian(void (*jacobian)(int *n_, double *t_, double *y_, double *a_))
	{
		Jacobian = jacobian; 
	}
	
	void SetInitialConditions(BzzVector& y0, const double t_start)
	{
		t_start_ = t_start;
		n_ = y0.Size();
		y_ = new double[n_];
		dpar_ = new double[(7+2*n_)*n_];
		kd_ = new int[n_];

		for(int j=0;j<n_;j++)
			y_[j] = y0[j+1];
	}



	void Solve(const double t_end, BzzVector& y_end_)
	{
		double t_ = t_start_;
		double t_end_ = t_end;

		dodesol_mk52lfa(ipar_, &n_, &t_, &t_end_, y_, Residuals, Jacobian, &h_, &hmin_, &ep_, &tr_, dpar_, kd_, &ierr_);
		
		for(int j=0;j<n_;j++)
			y_end_[j+1] = y_[j];
	}

private:

	int *ipar_;
	int n_;
	double h_;
	double hmin_;
	double ep_;
	double tr_;
	int *kd_;
	int ierr_;

	double t_start_;
	double *y_;
	double *dpar_;
	void (*Residuals)(int *n_, double *t_, double *y_, double *f_);
	void (*Jacobian)(int *n_, double *t_, double *y_, double *a_);

};


#endif // OpenSMOKE_ODEIntel_MK52lfa_H