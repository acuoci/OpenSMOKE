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

#include "mkl_lapacke.h"
#include <iomanip>
#include "OpenSMOKE_LAPACK_Dense.h"

void OpenSMOKE_LAPACK_Dense::ErrorMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class: " << "OpenSMOKE_LAPACK_Dense"	<< std::endl;
    std::cout << "Error: " << message_					<< std::endl;
    std::cout << "Press enter to continue... "			<< std::endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_LAPACK_Dense::WarningMessage(const std::string message_)
{
    std::cout << std::endl;
    std::cout << "Class:	  " << "OpenSMOKE_LAPACK_Dense"	<< std::endl;
    std::cout << "Warning: " << message_					<< std::endl;
	std::cout << std::endl;
}

OpenSMOKE_LAPACK_Dense::OpenSMOKE_LAPACK_Dense(void)
{
	//ErrorMessage("No default constructor is available");
}

OpenSMOKE_LAPACK_Dense::OpenSMOKE_LAPACK_Dense(const int n)
{
	Set(n);
}

void OpenSMOKE_LAPACK_Dense::Set(const int n)
{
	n_ = n;
	ipiv = new MKL_INT[n_];

	lda = n_;
	ldb = 1;
}

OpenSMOKE_LAPACK_Dense::~OpenSMOKE_LAPACK_Dense(void)
{
}

void OpenSMOKE_LAPACK_Dense::Solve(BzzMatrix& A, BzzVector &x)
{	
	int nrhs = 1;

	double* rhs=x.GetHandle();
	double* ptA = A.GetHandle();
	error_ = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n_, nrhs, ptA, lda, ipiv, rhs, ldb );
}

void OpenSMOKE_LAPACK_Dense::Solve(BzzMatrix& A, BzzVector &b, BzzVector &x)
{	
	int nrhs = 1;

	x = b;
	double* rhs = x.GetHandle();
	double* ptA = A.GetHandle();

	error_ = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n_, nrhs, ptA, lda, ipiv, rhs, ldb );
}

void OpenSMOKE_LAPACK_Dense::Delete()
{
	CleanMemory();
}

void OpenSMOKE_LAPACK_Dense::SetDefaultOptions()
{
	
}

void OpenSMOKE_LAPACK_Dense::CleanMemory()
{
	delete[] ipiv;
}

void OpenSMOKE_LAPACK_Dense::ErrorAnalysis()
{
	if (error_ == 0)	return;
	
	else if (error_ < 0)
	{
		ErrorMessage("Input inconsistent");
	}
	else
	{
		ErrorMessage("Matrix is singular");
	}
}



