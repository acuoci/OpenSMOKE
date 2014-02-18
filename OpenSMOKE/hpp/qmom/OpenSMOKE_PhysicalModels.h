/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci   	                               *
 *   alberto.cuoci@polimi.it   						                       *
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

#ifndef		QMOM_PHYSICALMODELS
#define		QMOM_PHYSICALMODELS

#include "BzzMath.hpp"

class OpenSMOKE_PhysicalModels  
{
private:
	int _N;
	int _2N;
	int _2N_minus_1;

	BzzVector csi;
	BzzVector w;
	BzzVector localSource;

	BzzMatrix auxiliary_2NxN;

	BzzMatrix betaKernel;
	BzzVector aKernel;
	BzzMatrix bKernel;

	BzzVector coeffA1;
	BzzVector coeffA2;
	BzzVector coeffA3;
	BzzVector Correction;


	void prepareCoefficients();
	void prepareA3();
	void computePowerOfCsi();
	void assemblingA3();



public:
	void setup(int N);
	void update(BzzVector &_w, BzzVector &_csi);
	void sources();

	void powerLawGrowth(double alfa, double rGrowth);
	void homogeneousDispersion(BzzVector &D);
	void homogeneousAggregation();
	void homogeneousBreakage();
	void nucleation(double epsilon, double J, BzzVector &L);

	void aggregationKernel(int kind, double parameter);
	void breakageKernel(int kind, double parameterA, double parameterB, double L0, double coefficient);
	void fragmentationKernel(int kind, double parameter);

	void diffusiveTermCorrectionForDQMOM(BzzVector &dcsi_over_dx, double Dmix);
	void diffusiveTermCorrectionForDQMOM(BzzVector &Cvector);

	void aggregationKernel_Soot(BzzMatrix &Kernel);


	BzzMatrix powerOfCsi;
	BzzMatrix A3;
	BzzVector Source;
};

#endif // !defined(AFX_PHYSICALMODELS_H__3770D047_07C4_4A48_B1F3_005FBF53B929__INCLUDED_)
