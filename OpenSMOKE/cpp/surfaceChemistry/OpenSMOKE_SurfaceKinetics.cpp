/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci   	                               *
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

#include <string>
#include "surfaceChemistry/OpenSMOKE_ReactingSurface.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceReaction.h"
#include "surfaceChemistry/OpenSMOKE_SurfaceKinetics.h"

void OpenSMOKE_SurfaceKinetics::ErrorMessage(const std::string message)
{
    cout << endl;
    cout << "Class:  OpenSMOKE_SurfaceKinetics"	<< endl;
    cout << "Object: " << name_object		<< endl;
    cout << "Error:  " << message           << endl;
    cout << "Press enter to continue... "   << endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_SurfaceKinetics::WarningMessage(const std::string message)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_SurfaceKinetics"	<< endl;
    cout << "Object:  " << name_object		<< endl;
    cout << "Warning: " << message          << endl;
	cout << endl;
}

OpenSMOKE_SurfaceKinetics::OpenSMOKE_SurfaceKinetics()
{
	name_object = "[not assigned]";
}

void OpenSMOKE_SurfaceKinetics::Setup(OpenSMOKE_ReactingSurface *surface, const std::string name_)
{
	ptSurface   = surface;
	name_object = name_;
}

void OpenSMOKE_SurfaceKinetics::Allocate()
{
	reaction_ = new OpenSMOKE_SurfaceReaction[numberReactions_+1];
}

void OpenSMOKE_SurfaceKinetics::ReadFromBinaryFile(BzzLoad &binaryFile)
{
	binaryFile >> numberReactions_;

	// Memory allocation
	Allocate();

	// Read reactions
	for (int j=1;j<=numberReactions_;j++)
		reaction_[j].ReadFromBinaryFile(binaryFile);

    coupleDirectAndInverseReactions();
}

void OpenSMOKE_SurfaceKinetics::ReactionRates(	const double T, const double lnT, const BzzVector &cGas, BzzVector &cSurface, BzzVector &ZSurface, BzzVector &aBulk,
						BzzVector &RGas, BzzVector &RSurface, BzzVector &RBulk, double &QReaction)
{
	QReaction = 0.;
	for (int j=1;j<=numberReactions_;j++)
	{
		double r = reaction_[j].ReactionRate(T, lnT, cGas, cSurface, ZSurface, aBulk, RGas, RSurface, RBulk);
		QReaction += dhrs[j]*r;
	}	
}

void OpenSMOKE_SurfaceKinetics::coupleDirectAndInverseReactions()
{
	BzzVectorInt inverseIndex;
    ChangeDimensions(numberReactions_, &inverseIndex);

  //  if( (numberReactions_%2) != 0)
    //    ErrorMessage("Odd number of reactions");

    for(int j = 1; j<=numberReactions_; j++)
    {
        if(inverseIndex[j] != -1)
        {

            for(int k = 2; k<=numberReactions_; k++)
            {

                if((reaction_[j].indexDirect == reaction_[k].indexInverse)&&(reaction_[j].indexInverse == reaction_[k].indexDirect))
                {
                    if((reaction_[j].nuDirect == reaction_[k].nuInverse)&&(reaction_[j].nuInverse == reaction_[k].nuDirect))
                    {
                        inverseIndex[j] = k;
                        inverseIndex[k] = -1;
                        k = numberReactions_;
                    }
                    else if(k == numberReactions_)
                        ErrorMessage("One or more inverse reactions not found");
                }
            }
        }
    }

    ChangeDimensions(numberReactions_, &dhrs);
    for (int j=1; j<=numberReactions_; j++)
    {
        if(inverseIndex[j] != -1)
            dhrs[j] =  (reaction_[j].E_over_R-reaction_[j+1].E_over_R)*Constants::R_J_kmol;
        else
            dhrs[j] = -(reaction_[j-1].E_over_R-reaction_[j].E_over_R)*Constants::R_J_kmol;
    }
}
