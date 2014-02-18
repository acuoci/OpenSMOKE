/***************************************************************************
 *   Copyright (C) 2003-2008 by                                            *
 *   Guido Buzzi-Ferraris, Alessio Frassoldati and Alberto Cuoci		   *
 *                                                                         *
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

#if !defined(OPENSMOKE_CSTRNETWORK_CORRECTIONCOEFFICIENTS_H)
#define OPENSMOKE_CSTRNETWORK_CORRECTIONCOEFFICIENTS_H

#include "BzzMath.hpp"

class OpenSMOKE_CSTRNetwork_CorrectionCoefficients
{
private:

	double	Tmin;
	double	Tmax;
	double	deltaT;
	static const double	epsilon;
	static const int	NSTEPSPDF;

	BzzVector x;
	BzzVector functionBeta;
	BzzVector functionCc;
	BzzVector functionI;
	

public:

	double I;
	
	void Setup(const double _Tmin, const double _Tmax);

	void GiveMeNormalizedVariables(const double T, const double TVN, double &csiMean, double &csiV);

	double Beta(const double a, const double b);
	double GiveMeCorrectionCoefficient(const double csi, const double csiMean, const double Tatt, const double n);


	double GiveMeCorrectionDoubleDirac(const double TMean, const double sigmaFluent, const double Tatt, const double n);
	double GiveMeCorrectionBetaPDF( const double csiMean, const double a, const double b, const double beta,
									const double Tatt, const double n);
	double GiveMeCorrectionClippedGaussianPDF(const double csiMean, const double x1,  const double x2, 
											  const double alfa1, const double alfa2, const double g, const double Tatt, const double n);
	void SetClippedGaussianPDF(const double csiMean, const double x1, const double x2, const double g);
};

#endif // !defined(OPENSMOKE_CSTRNETWORK_CORRECTIONCOEFFICIENTS_H)
