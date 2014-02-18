/***************************************************************************
 *   Copyright (C) 2006-2009 by Alberto Cuoci						       *
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

#if !defined(OPENSMOKE_LOOKUP_TABLE_FLAME)
#define OPENSMOKE_LOOKUP_TABLE_FLAME

#include "BzzMath.hpp"

class OpenSMOKE_LookUp_Table_Flame
{
public:

	BzzVector csi;
	BzzVector csiV;
	BzzMatrix temperature;
	BzzMatrix density;
	BzzMatrix cp;
	BzzMatrix as;

	int nCsi;
	int nCsiV;

	void read_from_file(const string fileName);
	
private:

	void ErrorMessage(const string message);
	void WarningMessage(const string message);
	string name_object;
};

#endif // !defined(OPENSMOKE_LOOKUP_TABLE_FLAME)

