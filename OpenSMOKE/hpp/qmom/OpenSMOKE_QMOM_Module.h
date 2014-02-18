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

#ifndef		OpenSMOKE_QMOMMODULE
#define		OpenSMOKE_QMOMMODULE

#include "qmom/OpenSMOKE_PhysicalModels.h"
#include "qmom/OpenSMOKE_Soot_Models.h"
#include "qmom/OpenSMOKE_QMOM.h"

class OpenSMOKE_QMOM_Module
{
public:

	// These functions must be called once
	void setupGasMixture(OpenSMOKE_ReactingGas &mix);
	void readFromFile(const string fileName);
	void initialize(int N, int iFractalDimension,
					int iNucleation, int iGrowth, int iOxidation,
					int iAggregation, double L0, int iNucleationDistribution,
					double seed_initial_distribution, int iTfluctuations);

	// These functions must be called in this order every time the updating of
	// moment equations is necessary
	void updateMoments(BzzVector &moments);
	void updateData(double &T, double &P, double &rho,
					double &mu, double &tr,
					double &omega_O2, double &omega_C2H2,
					double &omega_OH, double &SNTV);
	BzzVector calculateSources();

	int jC2H2;
	int jO2;
	int jOH;
	int N;

	BzzVector w;
	BzzVector csi;

	double seed_density;
	double epsilon;

	OpenSMOKE_Soot_Models soot;
	OpenSMOKE_PhysicalModels models;

	BzzVector L;

private:

	void MessageError(char *message);
	void finalize_initialization();


	OpenSMOKE_QMOM qmom;

	double L0;
	int iFractalDimension;
	int iNucleation;
	int iGrowth;
	int iOxidation;
	int iAggregation;
	int iNucleationDistribution;
	int iTfluctuations;

	BzzVector moments;
};

#endif // !defined(AFX_QMOM_MODULE_FOR_CSTR_H__CD17DF2B_384C_45CF_A2E6_8A0164FCE18D__INCLUDED_)
