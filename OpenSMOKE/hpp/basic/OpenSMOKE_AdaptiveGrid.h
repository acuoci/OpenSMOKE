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

#if !defined(OPENSMOKE_ADAPTIVEGRID)
#define OPENSMOKE_ADAPTIVEGRID

#include "BzzMath.hpp"
#include "basic/OpenSMOKE_Grid1D.h"
#include "basic/OpenSMOKE_Utilities.h"

class OpenSMOKE_AdaptiveGrid
{
public:

	OpenSMOKE_AdaptiveGrid();
	void SetName(const std::string name);

	void SetModel(adaptive_grid_model _model);
	void SetConstants(const double _alfa, const double _beta, const double _gamma);
	void SetAlfa(const double _alfa);
	void SetBeta(const double _beta);
	void SetGamma(const double _gamma);
	
	void Setup(BzzVector &_xOld, BzzVector &_yOld);
	void Ode(BzzVector &x, double t, BzzVector &f);
	void PrintOnFile(BzzVector &x);

private:

	int N;
	double alfa;
	double beta;
	double gamma;

	OpenSMOKE_Grid1D grid;
	adaptive_grid_model	model;
	void UpdateWeights();
	void UpdateWeightDerivatives();
	void UpdateProfileDerivatives(BzzVector &x);

	BzzVector xOld;
	BzzVector yOld;
	BzzVector dummy;
	BzzVector d2x_over_dcsi2;
	BzzVector dx_over_dcsi;
	BzzVector d2y_over_dx2;
	BzzVector dy_over_dx;
	BzzVector W;
	BzzVector dW_over_dcsi;
	BzzVector y;

private:

	std::string name_object;
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
};

#endif // !defined(OPENSMOKE_ADAPTIVEGRID)
