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

#if !defined(OPENSMOKE_GRID1D)
#define OPENSMOKE_GRID1D

#include "BzzMath.hpp"

class OpenSMOKE_Grid1D  
{
private:

public:

	OpenSMOKE_Grid1D();
	void SetName(const std::string name);
	void SetSpherical();

	int Np;
	int Ni;
	double xA;
	double xB;
	double L;

	BzzVector x, dxe, dxw;
	BzzVector x2, ux, udxe, udxw;
	BzzVector dxc, udxc, dxc_over_2, udxc_over_2;
	BzzVector cenp, cenc, cenm;
	BzzVector spherical, one_minus_spherical;
	BzzVector A;

	void Construct(const int Np, const double L, const double xAA);
	void Construct(const int Np, const double L, const double alfa, const double xAA);
	void ConstructStretchedStagnation(const int NP, const double LL, const double alfa, const double xAA);
	void Construct(const int N, const double LL, const double xcen, const double wmix, const double xAA);
	void Construct(const int NP, const double xStart, const double Length, 
				   const double alfaA, const double alfaB, const double fraction_of_points_A, double fraction_of_distance_A);
	void Construct(const std::string fileName);
	void Construct(BzzVector &coordinates);
	void ConstructExponential(const int NP, const double LL, const double Beta);
	
	void RecoverFromBackUp(const std::string fileName);
	void RecoverFromBackUp(BzzVector &_x);
	
	void RefineDouble();
	void Refine(const int j);
	int AddPoint(const double xPoint);
	
	BzzVectorInt QueryNewPointsDifference(const int nDiff, const double delta, BzzVector &phi);
	BzzVectorInt QueryNewPointsGradient(const int nGrad, const double delta, BzzVector &phi);

	void FirstDerivative(const char kind, BzzVector &v, BzzVector &dv);
	void FirstDerivative(const char kind, BzzVector &u, BzzVector &v, BzzVector &dv);
	
	void FirstDerivative(const char kind, BzzVector &u, BzzMatrix &v, BzzMatrix &dv);
	void FirstDerivative(const char kind, BzzMatrix &v, BzzMatrix &dv);

	void SecondDerivative(BzzVector &v, BzzVector &d2v);
	void SecondDerivative(BzzMatrix &v, BzzMatrix &d2v);

	void DoubleField(BzzVector &phi);
	void DoubleField(BzzMatrix &phi);
	void AddPointsField(BzzVector &phi, BzzVectorInt &listPoints);
	void AddPointsField(BzzVector &phi, BzzVectorInt &listPoints, const double ratio);
	void AddPointsField(BzzMatrix &phi, BzzVectorInt &listPoints);
	void AddPointsField(BzzMatrix &phi, BzzVectorInt &listPoints, const double ratio);
	void AdaptField(BzzVector &phi, BzzVector &xNew);
	void AdaptField(BzzMatrix &phi, BzzVector &xNew);

	void Rescale(const double Lnew);

private:

	bool iSpherical;

	std::string name_object;

	void Allocate();
	void Build();

	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);

};

#endif // !defined(OPENSMOKE_GRID1D)
