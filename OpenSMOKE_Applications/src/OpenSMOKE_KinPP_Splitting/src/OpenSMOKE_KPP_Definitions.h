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

#ifndef OpenSMOKE_KPP_Definitions_H
#define OpenSMOKE_KPP_Definitions_H


#include "BzzMath.hpp"
#include <iomanip>
using namespace std;

enum KPP_Correction { KPP_CORRECTION_NONE, KPP_CORRECTION_DIRAC, KPP_CORRECTION_BETA, KPP_CORRECTION_SIN, KPP_CORRECTION_GAUSS };

enum KPP_SparseLinearSolver { KPP_SPARSESOLVER_PARDISO, KPP_SPARSESOLVER_BZZ,
							  KPP_SPARSESOLVER_MUMPS, KPP_SPARSESOLVER_SUPERLU,
							  KPP_SPARSESOLVER_LIS, KPP_SPARSESOLVER_OPENSMOKE_GAUSSSIEDEL };

enum KPP_TransportIntegration {	KPP_TRANSPORTINTEGRATION_IMPLICIT_EULER, 
									KPP_TRANSPORTINTEGRATION_EXPLICIT_EULER, 
									KPP_TRANSPORTINTEGRATION_BZZODE };

enum KPP_Network_Status	{ KPP_NETWORK_STATUS_START, KPP_NETWORK_STATUS_SEQUENTIAL_CSTR, 
							  KPP_NETWORK_STATUS_SPLITTING, KPP_NETWORK_STATUS_SPLITTING_ALGEBRAIC,
							  KPP_NETWORK_STATUS_GLOBALODE, KPP_NETWORK_STATUS_GLOBALNLS };

enum KPP_NonLinearSystem_Method { KPP_NLS_KELLEY,  KPP_NLS_NEWTON, KPP_NLS_CHORD, KPP_NLS_SHAMANSKII };
enum KPP_NonLinearSystem_Exit   { KPP_NLS_SUCCESS, KPP_NLS_MAXIT,  KPP_NLS_FAILED };

enum linearEquationSolver { LINEAR_EQUATION_SOLVER_CG, LINEAR_EQUATION_SOLVER_BiCG, LINEAR_EQUATION_SOLVER_CGS,
							LINEAR_EQUATION_SOLVER_BiCGSTAB, LINEAR_EQUATION_SOLVER_BiCGSTABl, LINEAR_EQUATION_SOLVER_GPBiCG,
							LINEAR_EQUATION_SOLVER_TFQMR, LINEAR_EQUATION_SOLVER_Orthomin,
							LINEAR_EQUATION_SOLVER_GMRES, LINEAR_EQUATION_SOLVER_Jacobi, LINEAR_EQUATION_SOLVER_GaussSeidel,
							LINEAR_EQUATION_SOLVER_SOR, LINEAR_EQUATION_SOLVER_BiCGSafe, LINEAR_EQUATION_SOLVER_CR,
							LINEAR_EQUATION_SOLVER_BiCR, LINEAR_EQUATION_SOLVER_CRS, LINEAR_EQUATION_SOLVER_BiCRSTAB,
							LINEAR_EQUATION_SOLVER_GPBiCR, LINEAR_EQUATION_SOLVER_BiCRSafe, LINEAR_EQUATION_SOLVER_FGMRESm,
							LINEAR_EQUATION_SOLVER_IDRs, LINEAR_EQUATION_SOLVER_MINRES };

enum lisPreconditioner { LIS_PRECONDITIONER_NONE, LIS_PRECONDITIONER_Jacobi, LIS_PRECONDITIONER_ILUk,
						 LIS_PRECONDITIONER_SSOR, LIS_PRECONDITIONER_Hybrid, LIS_PRECONDITIONER_IS,
						 LIS_PRECONDITIONER_SAINV, LIS_PRECONDITIONER_SAAMG, LIS_PRECONDITIONER_CROUTILU,
						 LIS_PRECONDITIONER_ILUT };

enum  singleReactorConvergenceType  {	KPP_SINGLEREACTOR_CONVERGENCE_DIRECT, KPP_SINGLEREACTOR_CONVERGENCE_NEWTON,
										KPP_SINGLEREACTOR_CONVERGENCE_ODE_FIRST, KPP_SINGLEREACTOR_CONVERGENCE_ODE_SECOND,
										KPP_SINGLEREACTOR_CONVERGENCE_ODE_THIRD, KPP_SINGLEREACTOR_CONVERGENCE_ODE_FOURTH };

double parabolicModel(const double lambdac, const double lambdam, const double ff0, const double ffc, const double ffm);
void CleanVectorOnTheBottom(const double target, const double bottomToCheck, BzzVector &v);
void CleanVectorOnTheTop(const double target, const double topToCheck, BzzVector &v);


void GetSideGlobalPattern(const int index, const int nBlock, BzzMatrixSparse& C, BzzVector& patternLeft, BzzVector& patternRigth);
void GetGlobalPattern(const int dimBlock, BzzMatrixSparse& C, BzzVectorInt& rows, BzzVectorInt& columns);

template<typename T>
class RingVector
{
public:

	RingVector(const int n) 
	{
		n_ = n;
		count_ = 0;
		elements_ = new T[n_];
	}

	int Size() const { return n_;}
	
	T Element(const int i) const { return elements_[i+1]; }
	
	void Append(const T value) 
	{
		if (count_<n_)
		{
			elements_[count_]	= value;
			count_++;
		}
		else
		{
			for (int k=1;k<n_;k++)
				elements_[k-1] = elements_[k];
			elements_[n_-1] = value;
		}
	}

	bool IsFilled() const 
	{ 
		if (count_>=n_)
			return true;
		else
			return false;
	}

	double Mean() const
	{
		if (IsFilled() == true)
		{
			T sum = 0;
			for (int k=0;k<n_;k++)
				sum += elements_[k];
			return sum/double(n_);
		}
		else
			return 0.;
	}

	double MeanRatios() const
	{
		if (IsFilled() == true)
		{
			T sum = 0;
			for (int k=0;k<n_-1;k++)
				sum += elements_[k+1]/elements_[k];
			return sum/double(n_-1);
		}
		else
			return 0.;
	}

private:

	T* elements_;
	int n_;
	int count_;
};

class CPUTimer
{
public:
    
    CPUTimer(const string nameFile)
    {
        iteration_=0;
        kind_=0;
        Reset();
        OpenFile(nameFile);
    }
    
    void SetStartTimeAllGlobal()        { timeStartAllGlobal_ = BzzGetCpuTime(); }
    void SetStartTimeSequenceGlobal()   { timeStartSequenceGlobal_ = BzzGetCpuTime(); }
    void SetStartTimePredictorCorrectorGlobal() { timeStartPredictorCorrectorGlobal_ = BzzGetCpuTime(); }
    void SetStartTimeODEGlobal() { timeStartODEGlobal_ = BzzGetCpuTime(); }
    void SetStartTimeNLSGlobal() { timeStartNLSGlobal_ = BzzGetCpuTime(); }
    
    void SetEndTimeAllGlobal()        { timeEndAllGlobal_ = BzzGetCpuTime(); }
    void SetEndTimeSequenceGlobal()   { timeEndSequenceGlobal_ = BzzGetCpuTime(); }
    void SetEndTimePredictorCorrectorGlobal() { timeEndPredictorCorrectorGlobal_ = BzzGetCpuTime(); }
    void SetEndTimeODEGlobal() { timeEndODEGlobal_ = BzzGetCpuTime(); }
    void SetEndTimeNLSGlobal() { timeEndNLSGlobal_ = BzzGetCpuTime(); }
    
    void SetStartTimeSequenceLocal()   { timeStartSequenceLocal_ = BzzGetCpuTime(); }
    void SetStartTimePredictorCorrectorLocal() { timeStartPredictorCorrectorLocal_ = BzzGetCpuTime(); }
    void SetStartTimeODELocal() { timeStartODELocal_ = BzzGetCpuTime(); }
    void SetStartTimeNLSLocal() { timeStartNLSLocal_ = BzzGetCpuTime(); }
    
    void SetEndTimeSequenceLocal()   { timeEndSequenceLocal_ = BzzGetCpuTime(); }
    void SetEndTimePredictorCorrectorLocal() { timeEndPredictorCorrectorLocal_ = BzzGetCpuTime(); }
    void SetEndTimeODELocal() { timeEndODELocal_ = BzzGetCpuTime(); }
    void SetEndTimeNLSLocal() { timeEndNLSLocal_ = BzzGetCpuTime(); }    
    
    void SetSequence()
    {
        kind_ = 1;
        Reset();
    }
    
    void SetPredictorCorrector()
    {
        kind_ = 2;
        Reset();
    }
    
    void SetODEGlobal()
    {
        kind_ = 3;
        Reset();
    }    
    
    void SetNLSGlobal()
    {
        kind_ = 4;
        Reset();
    }
    
    
    void WriteOnFile()
    {
        fCPU_ << setw(12)  << left << iteration_++; 
        fCPU_ << setw(6)   << left << kind_; 
     
        fCPU_ << setw(16) << left << timeEndAllGlobal_ - timeStartAllGlobal_; 
        
        fCPU_ << setw(16) << left << timeEndSequenceLocal_  - timeStartSequenceLocal_; 
        fCPU_ << setw(16) << left << timeEndSequenceGlobal_  - timeStartSequenceGlobal_; 

        fCPU_ << setw(16) << left << timeEndPredictorCorrectorLocal_  - timeStartPredictorCorrectorLocal_; 
        fCPU_ << setw(16) << left << timeEndPredictorCorrectorGlobal_  - timeStartPredictorCorrectorGlobal_; 

        fCPU_ << setw(16) << left << timeEndODELocal_  - timeStartODELocal_; 
        fCPU_ << setw(16) << left << timeEndODEGlobal_  - timeStartODEGlobal_; 

        fCPU_ << setw(16) << left << timeEndNLSLocal_  - timeStartNLSLocal_;         
        fCPU_ << setw(16) << left << timeEndNLSGlobal_  - timeStartNLSGlobal_;         

        fCPU_ << endl;
    }
    
private:
    
    ofstream fCPU_;
    
    int kind_;
    int iteration_;
    
    double timeStartAllGlobal_;
    double timeStartSequenceGlobal_;
    double timeStartPredictorCorrectorGlobal_;
    double timeStartODEGlobal_;
    double timeStartNLSGlobal_;

    double timeEndAllGlobal_;
    double timeEndSequenceGlobal_;
    double timeEndPredictorCorrectorGlobal_;
    double timeEndODEGlobal_;
    double timeEndNLSGlobal_;
    
    double timeStartSequenceLocal_;
    double timeStartPredictorCorrectorLocal_;
    double timeStartODELocal_;
    double timeStartNLSLocal_;

    double timeEndSequenceLocal_;
    double timeEndPredictorCorrectorLocal_;
    double timeEndODELocal_;
    double timeEndNLSLocal_;    
    
    void Reset()
    {
        timeStartSequenceGlobal_ = timeStartPredictorCorrectorGlobal_ = 
                timeStartODEGlobal_ = timeStartNLSGlobal_ = 0;

        timeEndSequenceGlobal_ = timeEndPredictorCorrectorGlobal_ = 
                timeEndODEGlobal_ = timeEndNLSGlobal_ = 0;   
        
        timeStartSequenceLocal_ = timeStartPredictorCorrectorLocal_ = 
                timeStartODELocal_ = timeStartNLSLocal_ = 0;

        timeEndSequenceLocal_ = timeEndPredictorCorrectorLocal_ = 
                timeEndODELocal_ = timeEndNLSLocal_ = 0;           
    }
    
    void OpenFile(const string nameFile)
    {
        fCPU_.open(nameFile.c_str(), ios::out);
        fCPU_.setf(ios::scientific);
        
        fCPU_ << setw(12)  << left << "#It(1)";
        fCPU_ << setw(6)   << left << "K(2)"; 
     
        fCPU_ << setw(16) << left << "All(3)"; 
        
        fCPU_ << setw(16) << left << "Seq-local(4)";
        fCPU_ << setw(16) << left << "Seq-global(5)";

        fCPU_ << setw(16) << left << "PC-local(6)";
        fCPU_ << setw(16) << left << "PC-global(7)";

        fCPU_ << setw(16) << left << "ODE-local(8)";
        fCPU_ << setw(16) << left << "ODE-global(9)";

        fCPU_ << setw(16) << left << "NLS-local(10)";
        fCPU_ << setw(16) << left << "NLS-global(11)";     

        fCPU_ << endl;
        
    }
};

#endif // OpenSMOKE_KPP_Definitions_H