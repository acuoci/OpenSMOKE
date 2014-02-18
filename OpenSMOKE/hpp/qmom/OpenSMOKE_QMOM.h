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

#ifndef		QMOM_QMOM
#define		QMOM_QMOM

#include "qmom/OpenSMOKE_PDgordon.h"

class OpenSMOKE_PhysicalModels;

class OpenSMOKE_QMOM  
{
private:
	int _N;
	int _2N;
	int _2N_minus_1;

	OpenSMOKE_PDgordon PD;

public:

	BzzVector Source;

	int  GiveN();
	void setup(int N);
	void giveWeigthsAndAbscissasFromMoments( BzzVector &moments, 
		                                     BzzVector &w,
											 BzzVector &csi);
};

#endif // !defined(AFX_QMOM_H__1A82A5C3_2586_4F21_BBD1_1235F5D16DD0__INCLUDED_)
