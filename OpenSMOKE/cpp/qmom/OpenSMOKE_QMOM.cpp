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

#include "qmom/OpenSMOKE_QMOM.h"
#include "qmom/OpenSMOKE_PDgordon.h"

int  OpenSMOKE_QMOM::GiveN() 
{
	return _N;
}

void OpenSMOKE_QMOM::setup(int N)
{
	_N = N;
	_2N = 2*_N;
	_2N_minus_1 = _2N - 1;

	ChangeDimensions(_2N, &Source);

	PD.initialize();
}

void OpenSMOKE_QMOM::giveWeigthsAndAbscissasFromMoments( BzzVector &moments, 
											   BzzVector &w,
										       BzzVector &csi)
{
	PD.setup(moments);
	PD.calculateAbscissasAndWeigths();
		
	w = PD.w;
	csi = PD.csi;
}