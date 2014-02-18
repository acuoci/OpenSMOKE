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

#ifndef OpenSMOKE_LAPACK_Dense_H
#define OpenSMOKE_LAPACK_Dense_H

#include "BzzMath.hpp"
#include "mkl_lapacke.h"

class OpenSMOKE_LAPACK_Dense 
{
public:

	// Constructors and destructors
	OpenSMOKE_LAPACK_Dense(void);
	OpenSMOKE_LAPACK_Dense(const int n);
	~OpenSMOKE_LAPACK_Dense(void);

	// 1. Set from BzzMatrix
	void Set(const int n);

	// 2. Solve
	void Solve(BzzMatrix& A, BzzVector &b, BzzVector &x);
	void Solve(BzzMatrix& A, BzzVector &x);

	// 3. Delete
	void Delete();

public:

	void SetDefaultOptions();
	void UnsetMessage();

private:
	
	int error_;
	MKL_INT n_;
	MKL_INT* ipiv;		// Internal data

	MKL_INT lda;
	MKL_INT ldb;

private:
	
	string name_solver_;
	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
	void MessageOnTheScreen(const string message_);
	void CleanMemory();
	void ErrorAnalysis();
};

#endif	// OpenSMOKE_LAPACK_Dense_H
