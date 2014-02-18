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

#if !defined(OPENSMOKE_DROPLETMICROGRAVITY_GASRADIATION_H)
#define OPENSMOKE_DROPLETMICROGRAVITY_GASRADIATION_H

#include "OpenSMOKE.hpp"


class OpenSMOKE_DropletMicrogravity_GasRadiation
{
public:

	OpenSMOKE_DropletMicrogravity_GasRadiation();
	void SetName(const string name);

	void Initialize(const int N_);
	void Calculate(const BzzVector &x_, const BzzVector &T_, const BzzVector &Kp_);
	void SetOptions(const vector<string> options);
	
	BzzVector divq;
	double qSurface;

	BzzVector G0;
	BzzVector G1;
	BzzVector G2;
	double G3;

private:

	double div_q( const int j );
	double g2(const int j);
	double g3();

	double Phi(const double r, const double rho, const double z);
	void FillKernel();

	double SubIntegral_K(const double r, const double rho, const double z);
	double K(const double r, const double rho);
	
	double SubIntegral_Psi(const double rho, const double z);
	double Psi(const double rho);

	double SubIntegral_K(const double r, const double rho, const double z, const double kptilde);
	double SubIntegral_Psi(const double rho, const double z, const double kptilde);

	int N;
	int N_K;
	int N_subK;
	int N_Psi;
	int N_subPsi;

	BzzVector rTilde;
	BzzVector TTilde;
	BzzVector KpTilde;
	BzzVector rhoTilde;

	BzzVector T;
	BzzVector Kp;
	BzzVector y_aux;

	BzzMatrix Kernel;

	LinearInterpolation KpTilde_interpolation;

	enum radiation_mode { VISKANTA_MODE_NUMERICAL, VISKANTA_MODE_ANALYTICAL } mode;

private:

	string name_object;
	void ErrorMessage(const string message);
	void WarningMessage(const string message);
};

#endif // OPENSMOKE_DROPLETMICROGRAVITY_GASRADIATION_H
