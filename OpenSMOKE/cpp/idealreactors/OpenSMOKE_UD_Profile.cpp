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
#include "idealreactors/OpenSMOKE_UD_Profile.h"
#include "basic/OpenSMOKE_Conversions.h"

void OpenSMOKE_UD_Profile::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  User Defined Profile"		<< endl;
    cout << "Object: " << name_object			<< endl;
    cout << "Error:  " << message				<< endl;
    cout << "Press a key to continue... "		<< endl;
    getchar();
    exit(-1);
}

OpenSMOKE_UD_Profile::OpenSMOKE_UD_Profile()
{
	name_object		= "[Name not assigned]";	// Object Name
	iLinear			= false;
}

void OpenSMOKE_UD_Profile::SetName(const std::string name)
{
	name_object	= name;							// Object Name
}

bool OpenSMOKE_UD_Profile::CheckForComment(const std::string read_word)
{
    if(read_word.at(0) == '/')
        if(read_word.at(1) == '/')
            return true;
    return false;
}

bool OpenSMOKE_UD_Profile::CheckForKeyWord(const std::string read_word)
{
    if(read_word.at(0) == '#')
		return true;
    return false;
}

void OpenSMOKE_UD_Profile::AssignFromFile(const std::string fileName, const std::string name)
{
    const int SIZE = 300;
    char comment[SIZE];

	std::string kind_x;
    std::string kind_y;
    std::string unit_x;
    std::string unit_y;
    std::string dummy;
	bool    iGoBack;

    ifstream fInput;
    openInputFileAndControl(fInput, fileName.c_str());
	
	iGoBack = false;
    do
    {
		if (iGoBack == false)
			fInput >> dummy;
		iGoBack = false;

        if (dummy == "#END")
            break;
        else
        {
            if (!CheckForComment(dummy))
            {
				if (dummy == "#X")
				{
					fInput >> kind_x;
					fInput >> unit_x;
				}
				else if (dummy == "#Y")
				{
					fInput >> kind_y;
					fInput >> unit_y;
					if (kind_y != name)	ErrorMessage("Expected " + name + " key word. Found: " + kind_y);
				}
				else if (dummy == "#List")
				{
					for(;;)
					{
						fInput >> dummy;

						if (CheckForComment(dummy) == true)
						{
							fInput.getline(comment, SIZE);
							break;
						}

						if (CheckForKeyWord(dummy) == true)
						{
							iGoBack = true;
							break;
						}

						x.Append( atof(dummy.c_str()) );
						fInput >> dummy;
						y.Append( atof(dummy.c_str()) );
					}
				}
				else
					ErrorMessage("The following option is not recognized: " + dummy);
			}
			else
                fInput.getline(comment, SIZE);
		}
	}
	while(!fInput.eof());
    fInput.close();

	if		(kind_x == "TIME")		timeSupport = true;
    else if (kind_x == "LENGTH")	timeSupport = false;
    else							ErrorMessage("Expected TIME || LENGTH key word. Found: " + kind_x);

	if (timeSupport == true)		x *= OpenSMOKE_Conversions::conversion_time( 1.0, unit_x );
    else							x *= OpenSMOKE_Conversions::conversion_length( 1.0, unit_x );

	if		(name == "LENGTH")						y *= OpenSMOKE_Conversions::conversion_length( 1.0,	unit_y );
	else if (name == "AREA")						y *= OpenSMOKE_Conversions::conversion_area( 1.0, unit_y );
	else if (name == "TIME")						y *= OpenSMOKE_Conversions::conversion_time( 1.0, unit_y );
	else if	(name == "TEMPERATURE")					y *= OpenSMOKE_Conversions::conversion_temperature( 1.0, unit_y );
	else if (name == "HEAT_FLUX")					y *= OpenSMOKE_Conversions::conversion_heat_flux( 1.0, unit_y );
	else if (name == "HEAT_TRANSFER_COEFFICIENT")	y *= OpenSMOKE_Conversions::conversion_heat_exchange_coefficient( 1.0,		unit_y );
	else if (name == "SPECIFIC_AREA")				y *= OpenSMOKE_Conversions::conversion_u_length( 1.0, unit_y );
	else	ErrorMessage("Wrong key word: " + name);

	// In case of linear profile
	if (x.Size() == 2)
	{
		iLinear = true;
		yZero	= y[1];
		xMax	= 1.01*x[2];
		slope   = (y[2]-y[1])/x[2];
	}
	else
	{
		iLinear = false; 
		interpolation(x,y);
	}
}

void OpenSMOKE_UD_Profile::Setup(const std::string units_x, const std::string units_y, const double x0, const double x1, const double y0, const double y1)
{
	iLinear = true;
	timeSupport = true;
	double x0_ = OpenSMOKE_Conversions::conversion_time(x0, units_x);
	double x1_ = OpenSMOKE_Conversions::conversion_time(x1, units_x);
	double y0_ = OpenSMOKE_Conversions::conversion_temperature(y0, units_y);
	double y1_ = OpenSMOKE_Conversions::conversion_temperature(y1, units_y);

	yZero = y0_;
	xMax = 1.001*x1_;
	slope = (y1_-y0_)/(x1_-x0_);
}


/*
void OpenSMOKE_UD_Profile::AssignFromFile(const std::string fileName, const std::string name)
{
    std::string kind_x;
    std::string kind_y;
    std::string unit_x;
    std::string unit_y;
    std::string dummy;
    BzzVectorDouble x;
    BzzVectorDouble y;

    ifstream fInput;
    openInputFileAndControl(fInput, fileName.c_str());

    fInput >> dummy;
    if (dummy != "X")
        ErrorMessage("Expected X key word. Found: " + dummy);

    fInput >> kind_x;
    if		(kind_x == "TIME")  timeSupport = true;
    else if (kind_x == "SPACE") timeSupport = false;
    else                        ErrorMessage("Expected TIME || SPACE key word. Found: " + kind_x);
    fInput >> unit_x;

    fInput >> dummy;
    if (dummy != "Y")
        ErrorMessage("Expected Y key word. Found: " + dummy);

    fInput >> kind_y;
	if (kind_y != name)
		ErrorMessage("Expected " + name + " key word. Found: " + kind_y);
	fInput >> unit_y;

    for(;;)
    {
        fInput >> dummy;
        if (dummy == "//")  break;

        x.Append( atof(dummy.c_str()) );
        
        fInput >> dummy;
        y.Append( atof(dummy.c_str()) );
    }

	if (timeSupport == true)		x *= OpenSMOKE_Conversions::conversion_time( 1.0, unit_x );
    else							x *= OpenSMOKE_Conversions::conversion_length( 1.0, unit_x );

	if		(name == "SPACE")						y *= OpenSMOKE_Conversions::conversion_length( 1.0,	unit_y );
	else if (name == "TIME")						y *= OpenSMOKE_Conversions::conversion_time( 1.0, unit_y );
	else if	(name == "TEMPERATURE")					y *= OpenSMOKE_Conversions::conversion_temperature( 1.0, unit_y );
	else if (name == "HEAT_FLUX")					y *= OpenSMOKE_Conversions::conversion_heat_flux( 1.0, unit_y );
	else if (name == "HEAT_TRANSFER_COEFFICIENT")	y *= OpenSMOKE_Conversions::conversion_heat_exchange_coefficient( 1.0,		unit_y );
	else if (name == "SPECIFIC_AREA")				y *= OpenSMOKE_Conversions::conversion_u_length( 1.0, unit_y );
	else	ErrorMessage("Wrong key word: " + name);

	// In case of linear profile
	if (x.Size() == 2)
	{
		iLinear = true;
		yZero	= y[1];
		xMax	= 1.01*x[2];
		slope   = (y[2]-y[1])/x[2];
	}
	else
	{
		iLinear = false; 
		interpolation(x,y);
	}
}
*/
void OpenSMOKE_UD_Profile::Check(const double abscissa, const double expected_value)
{
	double difference;
	if (iLinear == true)	difference = fabs(expected_value - yZero);
	else					difference = fabs(expected_value - interpolation(abscissa));
    
	if (difference > 1.e-8)
        ErrorMessage("Check Failure. The actual value does not match the expected value.");
}

double OpenSMOKE_UD_Profile::GiveMeValue(const double t, const double x)
{
    if (timeSupport == true)
	{	
		if (iLinear == true)	
		{
			if (t>=0.0 && t<=xMax)	return (yZero + slope*t);
			else 
			{
				stringstream xMaxAccepted;
				stringstream xMaxFound;
				xMaxAccepted	<< xMax;
				xMaxFound		<< t;

				ErrorMessage("Out of boundaries during integration - MaxAccepted: " + xMaxAccepted.str() + " Found: " + xMaxFound.str());
				return -1;
			}
		}
		else					return interpolation(t);
	}
    else
	{
		if (iLinear == true)	
		{
			if (x>=0.0 && x<=xMax)	return (yZero + slope*x);
			else
			{
				stringstream xMaxAccepted;
				stringstream xMaxFound;
				xMaxAccepted	<< xMax;
				xMaxFound		<< x;

				ErrorMessage("Out of boundaries during integration - MaxAccepted: " + xMaxAccepted.str() + " Found: " + xMaxFound.str());
				return -1;
			}
		}
		else					return interpolation(x);
	}
}


