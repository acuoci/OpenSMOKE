/***************************************************************************
 *   Copyright (C) 2007 by Alberto Cuoci								   *
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

#ifndef OpenSMOKE_Sundials_CVODE_H
#define OpenSMOKE_Sundials_CVODE_H

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <cvode/cvode_lapack.h>      /* prototype for CVLapackDense */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., and macros */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include "BzzMath.hpp"

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

class OpenSMOKE_Sundials_CVODE
{
public:
	OpenSMOKE_Sundials_CVODE(void)
	{
		y = abstol = NULL;
		cvode_mem = NULL;
	}

	~OpenSMOKE_Sundials_CVODE(void) {};

	void SetInitialConditions(BzzVector& y0, const double t_start)
	{
		t_start_ = t_start;
		NEQ = y0.Size();
		
		y = N_VNew_Serial(NEQ);
		if (check_flag((void *)y, "N_VNew_Serial", 0))
			ErrorMessage();
		
		abstol = N_VNew_Serial(NEQ); 
		if (check_flag((void *)abstol, "N_VNew_Serial", 0))
			ErrorMessage();

		for(int j=1;j<=NEQ;j++)
			Ith(y,j) = y0[j];

	}

	void SetTolerances(const double relTolerance, BzzVector& absTolerance)
	{
		reltol = relTolerance;

		for(int j=1;j<=NEQ;j++)
			Ith(abstol,j) = absTolerance[j];
	}

	void Solve(const double t_end, BzzVector& solution)
	{
		/* Call CVodeCreate to create the solver memory and specify the 
		* Backward Differentiation Formula and the use of a Newton iteration */
		cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) 
			ErrorMessage();
  
		/* Call CVodeInit to initialize the integrator memory and specify the
		* user's right hand side function in y'=f(t,y), the inital time T0, and
		* the initial dependent variable vector y. */
		flag = CVodeInit(cvode_mem, f, t_start_, y);
		if (check_flag(&flag, "CVodeInit", 1)) 
			ErrorMessage();

		/* Call CVodeSVtolerances to specify the scalar relative tolerance
		* and vector absolute tolerances */
		flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
		if (check_flag(&flag, "CVodeSVtolerances", 1)) 
			ErrorMessage();

		/* Call CVodeRootInit to specify the root function g with 2 components */
		flag = CVodeRootInit(cvode_mem, 2, g);
		if (check_flag(&flag, "CVodeRootInit", 1)) 
			ErrorMessage();

		/* Call CVDense to specify the CVDENSE dense linear solver */
		flag = CVDense(cvode_mem, NEQ);
		if (check_flag(&flag, "CVDense", 1)) 
			ErrorMessage();

		/* Set the Jacobian routine to Jac (user-supplied) */
		flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
		if (check_flag(&flag, "CVDlsSetDenseJacFn", 1))
			ErrorMessage();

		/* In loop, call CVode, print results, and test for error.
		* Break out of loop when NOUT preset output times have been reached.  */
		printf(" \n3-species kinetics problem\n\n");

		{
			iout = 0;  
			tout = t_end;

			flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
   
			cout << Ith(y,1) << endl;
			getchar();

			//PrintOutput(t, Ith(y,1), Ith(y,2), Ith(y,3));

		//	if (flag == CV_ROOT_RETURN)
		//	{
		//		flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
		//		if (check_flag(&flagr, "CVodeGetRootInfo", 1)) 
		//			ErrorMessage();
			//	PrintRootInfo(rootsfound[0],rootsfound[1]);
		//	}
		}
	}

private:

	realtype t, tout, reltol;
	N_Vector y, abstol;
	void *cvode_mem;
	int flag, flagr, iout;
	int rootsfound[2];
	int NEQ;
	double t_start_;



	static int check_flag(void *flagvalue, char *funcname, int opt);

	static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
	{
		realtype y1, y2, y3, yd1, yd3;

		y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

		yd1 = Ith(ydot,1) = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
		yd3 = Ith(ydot,3) = RCONST(3.0e7)*y2*y2;
			Ith(ydot,2) = -yd1 - yd3;

		return(0);
	}


	static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
	{
		realtype y1, y3;

		y1 = Ith(y,1); y3 = Ith(y,3);
		gout[0] = y1 - RCONST(0.0001);
		gout[1] = y3 - RCONST(0.01);

		return(0);
	}

	static int Jac(int N, realtype t,	N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
	{
		realtype y1, y2, y3;

		y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

		IJth(J,1,1) = RCONST(-0.04);
		IJth(J,1,2) = RCONST(1.0e4)*y3;
		IJth(J,1,3) = RCONST(1.0e4)*y2;

		IJth(J,2,1) = RCONST(0.04); 
		IJth(J,2,2) = RCONST(-1.0e4)*y3-RCONST(6.0e7)*y2;
		IJth(J,2,3) = RCONST(-1.0e4)*y2;

		IJth(J,3,1) = 0;
		IJth(J,3,2) = RCONST(6.0e7)*y2;
		IJth(J,3,3) = 0;

		return(0);
	}

	void ErrorMessage()
	{
		cout << "ErrorMessage" << endl;
		getchar();
		exit(-1);
	}
};

#endif // OpenSMOKE_Sundials_CVODE_H