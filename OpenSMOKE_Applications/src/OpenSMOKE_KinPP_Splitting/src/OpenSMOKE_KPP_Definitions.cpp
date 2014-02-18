/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
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

#include "OpenSMOKE_KPP_Definitions.h"


// Adapted from: C. T. Kelley, April 1, 2003
// Apply three-point safeguarded parabolic model for a line search
double parabolicModel(const double lambdac, const double lambdam, const double ff0, const double ffc, const double ffm)
{
	double lambdap;

	// Internal parameters: safeguarding bounds for the linesearch
	double sigma0 = .1;
	double sigma1 = .5;

	// Compute coefficients of interpolation polynomial

	double c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);

	// Negative curvature and default solution
	if (c2 >= 0)
	{	
		lambdap = sigma1*lambdac;
		return lambdap;
	}

	// Positive curvature
	double c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
	lambdap = -c1*.5/c2;
	
	if (lambdap < sigma0*lambdac)
		lambdap = sigma0*lambdac;
	if (lambdap > sigma1*lambdac)
		lambdap = sigma1*lambdac;

	return lambdap;
}


void CleanVectorOnTheBottom(const double target, const double bottomToCheck, BzzVector &v)
{
	for(int k=1;k<=v.Size();k++)
		if (v[k] < target)
			if (v[k] > bottomToCheck)
				v[k] = target;
}


void CleanVectorOnTheTop(const double target, const double topToCheck, BzzVector &v)
{
	for(int k=1;k<=v.Size();k++)
		if (v[k] > target)
			if (v[k] < topToCheck)
				v[k] = target;	
}

void GetSideGlobalPattern(const int index, const int nBlock, BzzMatrixSparse& C, BzzVector& patternLeft, BzzVector& patternRigth)
{
	int nLeftElements=0;
	int nRigthElements=0;

	// Count elements
	{
		ElementBzzMatrixSparse *elem;
		elem = C.GetStartingElementInRow(index);
		while(elem)
		{
			if (elem->column < index)
				nLeftElements++;
			else if (elem->column > index)
				nRigthElements++;

			elem = elem->next;
		}

		ChangeDimensions(nLeftElements,  &patternLeft);
		ChangeDimensions(nRigthElements, &patternRigth);
	}

	// Assign Pattern
	{
		int iLeft=1;
		int iRigth=1;

		ElementBzzMatrixSparse *elem;
		elem = C.GetStartingElementInRow(index);
		while(elem)
		{
			if (elem->column < index)
				patternLeft[iLeft++] = (elem->column-1)*nBlock;
			else if (elem->column > index)
				patternRigth[iRigth++] = (elem->column-1)*nBlock;

			elem = elem->next;
		}
	}
}

void GetGlobalPattern(const int dimBlock, BzzMatrixSparse& C, BzzVectorInt& rows, BzzVectorInt& columns)
{
	int numberBlocks = C.Rows();

	BzzVector numberNonZeroElementsPerBlockRow(numberBlocks);
	for(int k=1;k<=numberBlocks;k++)
	{
		BzzVector patternLeft;
		BzzVector patternRigth;
		GetSideGlobalPattern(k, dimBlock, C, patternLeft, patternRigth);
		numberNonZeroElementsPerBlockRow[k] = dimBlock + patternLeft.Size() + patternRigth.Size();
	}

	int numberNonZeroElements = numberNonZeroElementsPerBlockRow.GetSumElements()*dimBlock;

	cout << numberNonZeroElements << endl;
}