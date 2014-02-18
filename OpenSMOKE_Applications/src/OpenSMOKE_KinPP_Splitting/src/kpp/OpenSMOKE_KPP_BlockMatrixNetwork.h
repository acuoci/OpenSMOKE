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

#ifndef OpenSMOKE_KPP_BlockMatrixNetwork_H
#define OpenSMOKE_KPP_BlockMatrixNetwork_H

#include "BzzMath.hpp"
#include "OpenSMOKE.hpp"
#include "OpenSMOKE_KPP_Definitions.h"

class OpenSMOKE_KPP_Dictionary;

class OpenSMOKE_KPP_BlockMatrixNetwork
{
public:

	// Constructor-Destructor
	OpenSMOKE_KPP_BlockMatrixNetwork();
	~OpenSMOKE_KPP_BlockMatrixNetwork(void);

	void SetSparsityPattern(const int NR, const int NC);
	void SetSparsityPattern(const int iBlock, const BzzVectorInt& indices);

	void SetBlockAndDiagonals(const int iBlock, BzzMatrix &block, BzzVector &diagonals);
	void SetRHS(const int iBlock, const BzzVector &rhs_local, const BzzVector &rhs_non_local);

	int GaussSiedel(BzzVector& solution);

	void SetAbsoluteTolerance(const double absTolerance);
	void SetRelativeTolerance(const double relTolerance);

private:

	int NR_;
	int NC_;
	int NE_;

	BzzMatrix* Block_;
	BzzFactorizedGauss* J_;
	BzzVector* xOld_;
	BzzVector* x_;
	BzzVector* b_;
	BzzVector* rhs_;
	BzzVectorInt* iDiagonals_;
	BzzVector* diagonals_;

	BzzVector localRow_;

	int maxIterations_;
	double absTolerance_;
	double relTolerance_;

	void SetBlock(const int iBlock, BzzMatrix &block);
	void SetDiagonals(const int iBlock, BzzVector &diagonals);

	BzzVector maxDx_;
	BzzVector maxX_;

	bool iAlreadyInitialized;

private:

	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
};

#endif	// OpenSMOKE_KPP_BlockMatrixNetwork_H



