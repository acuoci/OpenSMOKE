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

#ifndef QMOM_DQMOM
#define QMOM_DQMOM

#include "BzzMath.hpp"

class OpenSMOKE_PhysicalModels;

class OpenSMOKE_DQMOM  
{
private:
	int _N;
	int _2N;
	int _N_plus_1;
	int _2N_minus_1;

	BzzFactorizedGauss A;
	BzzVector coeffA1;
	BzzVector coeffA2;

public:


	BzzVector Source;
	BzzVector a;
	BzzVector b;
	BzzVector alfa;

	BzzMatrix powerOfCsi;


	int  GiveN() {return _N;};
	void setup(int N);
	void prepareCoefficients();
	void prepareA();
	void update(OpenSMOKE_PhysicalModels &_models);
	void assemblingA();

	void solveLinearSystem();
};

#endif // !defined(QMOM_DQMOM)
