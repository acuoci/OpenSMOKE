/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci   	   *
 *   alberto.cuoci@polimi.it   						   *
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

#include <string>
#include "surfaceChemistry/OpenSMOKE_SurfaceMaterial.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceBulk.h"
using namespace std;


void OpenSMOKE_SurfaceBulk::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SurfaceBulk"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SurfaceBulk::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SurfaceBulk"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

OpenSMOKE_SurfaceBulk::OpenSMOKE_SurfaceBulk()
{
	name_object = "[not assigned]";
}

void OpenSMOKE_SurfaceBulk::Setup(OpenSMOKE_SurfaceMaterial* material, const std::string name_object_)
{
	ptMaterial  = material;
	name_object = name_object_;
}


