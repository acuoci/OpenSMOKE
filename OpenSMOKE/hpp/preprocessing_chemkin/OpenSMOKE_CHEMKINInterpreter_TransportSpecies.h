/***************************************************************************
 *   Copyright (C) 2009 by Alberto Cuoci								   *
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

#if !defined(OPENSMOKE_CHEMKININTERPRETER_TRANSPORTSPECIES_H)
#define OPENSMOKE_CHEMKININTERPRETER_TRANSPORTSPECIES_H

#include "BzzMath.hpp"
#include <vector>

class OpenSMOKE_CHEMKINInterpreter_TransportSpecies  
{
public:
	OpenSMOKE_CHEMKINInterpreter_TransportSpecies();
	virtual ~OpenSMOKE_CHEMKINInterpreter_TransportSpecies();

	int index_line;
	
	std::string	name;
	int		shape_factor;
	double	epsylon_over_kb;
	double	sigma;
	double	mu;
	double	alfa;
	double	zRot298;

	void ReadMainData(const std::string line, const int iLine);
	void Analyze();

private:
	
	void ErrorMessage(const std::string message);
	void WarningMessage(const std::string message);
	std::string name_object;
};

#endif // !defined(OPENSMOKE_CHEMKININTERPRETER_TRANSPORTSPECIES_H)
