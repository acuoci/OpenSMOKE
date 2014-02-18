/***************************************************************************
 *   Copyright (C) 2006-2008 by Alberto Cuoci							   *
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

#include <sstream>
#include "basic/OpenSMOKE_Constants.h"
#include "idealreactors/cstr/OpenSMOKE_CSTR_Geometry.h"
#include "basic/OpenSMOKE_Conversions.h"

const double OpenSMOKE_CSTR_Geometry::pi_over_6  = Constants::pi / 6.;
const double OpenSMOKE_CSTR_Geometry::_6_over_pi = 6. / Constants::pi;

void OpenSMOKE_CSTR_Geometry::ErrorMessage(const string message)
{
    cout << endl;
    cout << "Class:  CSTR Geometry"		<< endl;
    cout << "Object: " << name_object					<< endl;
    cout << "Error:  " << message						<< endl;
    cout << "Press a key to continue... "				<< endl;
    getchar();
    exit(-1);
}

OpenSMOKE_CSTR_Geometry::OpenSMOKE_CSTR_Geometry()
{
    iKind		= 0;				// The diameter remains constant
	name_object	= "Default name";	// Object Name
}

void OpenSMOKE_CSTR_Geometry::SetName(const string name)
{
	name_object	= name;				// Object Name
}

void OpenSMOKE_CSTR_Geometry::Setup(const string fileName)
{
    int i;
    string label;
    string variable;
    string conversion_x;
    string conversion_y;
    BzzVector x_vector;
    BzzVector y_vector;

    ifstream fInput;
    openInputFileAndControl(fInput, fileName.c_str());

    fInput >> label;
	if (label != "#X")
        ErrorMessage("Only the following options are available: X");

	fInput >> label;
	if (label != "TIME")
        ErrorMessage("Only the following options are available: TIME");

	fInput >> conversion_x;

	fInput >> label;
	if (label != "#Y")
        ErrorMessage("Only the following options are available: Y");

	fInput >> variable;
    if (variable != "TAU" && variable != "VOLUME" && variable != "DIAMETER")
        ErrorMessage("Only the following options are available: DIAMETER || VOLUME || TAU");

	fInput >> conversion_y;

	fInput >> label;
	if (label != "#List")
		ErrorMessage("Expected #List key word...");

    for(;;)
    {
        fInput >> label;
		if (label == "#END")  break;

        x_vector.Append( atof(label.c_str()) );
        
        fInput >> label;
        y_vector.Append( atof(label.c_str()) );
    }

	// Checking values
	if (x_vector[1] != 0.)	ErrorMessage("The abscissas must start from 0.0");
	if (y_vector[1] <= 0.)	ErrorMessage("The initial y value must be larger than 0.0");

	// In case of linear profile
	if (x_vector.Size() == 2)
	{
		if (variable == "DIAMETER")
		{
			iKind = 1;

			x_vector *= OpenSMOKE_Conversions::conversion_time(1.0, conversion_x);
			y_vector *= OpenSMOKE_Conversions::conversion_length(1.0, conversion_y);

			D0 = y_vector[1];
			slope_diameter = (y_vector[2]-y_vector[1]) / x_vector[2];
		}
		else if (variable == "VOLUME")
		{
			iKind = 2;

			x_vector *= OpenSMOKE_Conversions::conversion_time(1.0, conversion_x);
			y_vector *= OpenSMOKE_Conversions::conversion_volume(1.0, conversion_y);

			Volume0 = y_vector[1];
			slope_volume = (y_vector[2]-y_vector[1]) / x_vector[2];
		}
		else if (variable == "TAU")
		{
			iKind = 3;

			x_vector *= OpenSMOKE_Conversions::conversion_time(1.0, conversion_x);
			y_vector *= OpenSMOKE_Conversions::conversion_time(1.0, conversion_y);

			Tau0 = y_vector[1];
			slope_tau = (y_vector[2]-y_vector[1]) / x_vector[2];
		}
	}
	
    else 
	{
		if (variable == "DIAMETER")
		{
			iKind = 4;

			x_vector *= OpenSMOKE_Conversions::conversion_time(1.0, conversion_x);
			y_vector *= OpenSMOKE_Conversions::conversion_length(1.0, conversion_y);

			BzzVector z_vector(x_vector.Size());
			for(i=1;i<=x_vector.Size()-1;i++)
				z_vector[i] = (y_vector[i+1]-y_vector[i]) / (x_vector[i+1]-x_vector[i]);
			i = x_vector.Size();
			z_vector[i] = (y_vector[i]-y_vector[i-1]) /	(x_vector[i]-x_vector[i-1]);

			interpolation_diameter(x_vector, y_vector);
			interpolation_dDiameter(x_vector, z_vector);

		}
		else if (variable == "VOLUME")
		{
			iKind = 5;

			x_vector *= OpenSMOKE_Conversions::conversion_time(1.0, conversion_x);
			y_vector *= OpenSMOKE_Conversions::conversion_volume(1.0, conversion_y);

			BzzVector z_vector(x_vector.Size());
			for(i=1;i<=x_vector.Size()-1;i++)
				z_vector[i] = (y_vector[i+1]-y_vector[i])/(x_vector[i+1]-x_vector[i]);
			i = x_vector.Size();
			z_vector[i] = (y_vector[i]-y_vector[i-1])/(x_vector[i]-x_vector[i-1]);
			
			interpolation_volume(x_vector, y_vector);
			interpolation_dVolume(x_vector, z_vector);
		}
		else if (variable == "TAU")
		{
			iKind = 6;

			x_vector *= OpenSMOKE_Conversions::conversion_time(1.0, conversion_x);
			y_vector *= OpenSMOKE_Conversions::conversion_time(1.0, conversion_y);

			BzzVector z_vector(x_vector.Size());
			for(i=1;i<=x_vector.Size()-1;i++)
				z_vector[i] = (y_vector[i+1]-y_vector[i])/(x_vector[i+1]-x_vector[i]);
			i = x_vector.Size();
			z_vector[i] = (y_vector[i]-y_vector[i-1])/(x_vector[i]-x_vector[i-1]);
			
			interpolation_tau(x_vector, y_vector);
			interpolation_dTau(x_vector, z_vector);
		}
	}
}


void OpenSMOKE_CSTR_Geometry::Update(const double x, double &D, double &Volume, double &Area)
{
    switch(iKind)
    {
        case 0:

        break;

        case 1:     // Linear diameter
            D		= D0 + slope_diameter*x;
			Area	= Constants::pi * D*D;
            Volume	= pi_over_6 * D*D*D;

        break;

        case 2:     // Linear volume
            Volume	= Volume0 + slope_volume*x;
            D		= pow(_6_over_pi * Volume, 1./3.);
			Area	= Constants::pi * D*D;

        break;

        case 3:    // Tau Points
            ErrorMessage("Indipendent variable cannot be the residence time");
        break;

        case 4:    // Diameter points
            D		= interpolation_diameter(x);
			Area	= Constants::pi * D*D;
            Volume	= pi_over_6 * D*D*D;
		
		break;

        case 5:    // Volume points
            Volume	= interpolation_volume(x);
            D		= pow(_6_over_pi * Volume, 1./3.);
			Area	= Constants::pi * D*D;

        case 6:    // Tau Points
            ErrorMessage("Indipendent variable cannot be the residence time");
        break;
    }
}

void OpenSMOKE_CSTR_Geometry::Update(const double x, double &Tau)
{
    switch(iKind)
    {
        case 0:
        break;

        case 1:    // Diameter Points
            ErrorMessage("Indipendent variable cannot be the reactor diameter");
        break;

        case 2:    // Volume Points
            ErrorMessage("Indipendent variable cannot be the reactor volume");
        break;

        case 3:    // Linear diameter
            Tau	= Tau0 + slope_tau*x;

        break;

        case 4:    // Diameter Points
            ErrorMessage("Indipendent variable cannot be the reactor diameter");
        break;

        case 5:    // Volume Points
            ErrorMessage("Indipendent variable cannot be the reactor volume");
        break;

        case 6:    // Tau Points
             Tau = interpolation_tau(x);
        break;
    }
}