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

#ifndef OpenSMOKE_QMOM_DISTRIBUTIONS
#define OpenSMOKE_QMOM_DISTRIBUTIONS

#include "BzzMath.hpp"

class OpenSMOKE_Distributions
{
private:

	char KIND;
	virtual void control() = 0;
public:
	virtual void GetMomentOrderOne(double &m) = 0;
	virtual void GetMomentOrderTwo(double &m) = 0;
	virtual void GetMomentOrderThree(double &m) = 0;
	virtual void GetMomentOrderFour(double &m) = 0;
	virtual void GetMomentRatioThreeOverTwo(double &d32) = 0;
	virtual void GetMomentRatioFourOverThree(double &d43) = 0;
	virtual void GetMoment(int numberOfMoment, double &m) = 0;
	virtual void GetMoments(int numberOfMoments) = 0;
	virtual void printOnVideo();
	virtual double variance();


public:
	BzzVector moments;
};

class OpenSMOKE_DiracDistribution : public OpenSMOKE_Distributions
{
private:

	double csiN;
	void control();
public:

	void setup();
	void update(double _csiN);
	void GetMomentOrderOne(double &m);
	void GetMomentOrderTwo(double &m);
	void GetMomentOrderThree(double &m);
	void GetMomentOrderFour(double &m);
	void GetMomentRatioThreeOverTwo(double &d32);
	void GetMomentRatioFourOverThree(double &d43);
	void GetMoment(int orderOfMoment, double &m);
	void GetMoments(int maxOrderOfMoments);

};

class OpenSMOKE_MultiDiracDistribution : public OpenSMOKE_Distributions
{
private:

	int N;


	double epsilon;
	double _1_plus_epsilon;
	double _1_minus_epsilon;
	void control();
public:
	BzzVector csiN;
	BzzVector pN;
	void setup(int _N);
	void update(BzzVector &_csiN, BzzVector &_pN);
	void updateWithoutControl(BzzVector &_csiN, BzzVector &_pN);

	void GetMomentOrderOne(double &m);
	void GetMomentOrderTwo(double &m);
	void GetMomentOrderThree(double &m);
	void GetMomentOrderFour(double &m);
	void GetMomentRatioThreeOverTwo(double &d32);
	void GetMomentRatioFourOverThree(double &d43);

	void GetMoment(int orderOfMoment, double &m);
	void GetMoments(int maxOrderOfMoments);
};

class OpenSMOKE_ConstantDistribution : public OpenSMOKE_Distributions
{
private:

	double csiA, csiB;
	double H;
	void control();
public:

	void setup();
	void update(double _csiA, double _csiB, double _H);

	void GetMomentOrderOne(double &m);
	void GetMomentOrderTwo(double &m);
	void GetMomentOrderThree(double &m);
	void GetMomentOrderFour(double &m);
	void GetMomentRatioThreeOverTwo(double &d32);
	void GetMomentRatioFourOverThree(double &d43);

	void GetMoment(int orderOfMoment, double &m);
	void GetMoments(int maxOrderOfMoments);
};

class OpenSMOKE_KM_Cloud_Distribution : public OpenSMOKE_Distributions
{
private:

	double a, b;
	void control();
public:

	void setup();
	void update(double _a, double _b);

	void GetMomentOrderOne(double &m);
	void GetMomentOrderTwo(double &m);
	void GetMomentOrderThree(double &m);
	void GetMomentOrderFour(double &m);
	void GetMomentRatioThreeOverTwo(double &d32);
	void GetMomentRatioFourOverThree(double &d43);

	void GetMoment(int orderOfMoment, double &m);
	void GetMoments(int maxOrderOfMoments);
};

#endif // QMOM_DISTRIBUTIONS
