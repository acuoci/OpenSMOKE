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

#ifndef OPENSMOKE_KPP_SingleReactor_H
#define OPENSMOKE_KPP_SingleReactor_H

#include "BzzMath.hpp"
#include "engine/OpenSMOKE_ReactingGas.h"
#include "OpenSMOKE_KPP_Definitions.h"

class OpenSMOKE_KPP_SingleReactor;
class OpenSMOKE_KPP_ReactorNetwork;
class OpenSMOKE_KPP_DataManager;
class OpenSMOKE_KPP_SingleReactor_KineticsManager;

class MyOdeSystem_KPP_ContinousReactor : public BzzOdeSystemObject
{
private:
public:
	OpenSMOKE_KPP_SingleReactor *ptReactor;
	void assignReactor(OpenSMOKE_KPP_SingleReactor *reactor);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &dy);
	virtual void GetJacobian(BzzVector &y,double t,BzzMatrix &JJ);
	virtual void ObjectBzzPrint(void);
};

class OpenSMOKE_KPP_SingleReactor
{
public:
	OpenSMOKE_KPP_SingleReactor();
	~OpenSMOKE_KPP_SingleReactor(void);

	// Functions which are called once
	void Setup(const int i, OpenSMOKE_ReactingGas* mixture, BzzVectorInt& indexReduced);
	void SetDataManager(OpenSMOKE_KPP_DataManager* data);
	void ResetStatus();
	bool ReadMassFractions(ifstream &fInput);
	bool ReadReactorProperties(ifstream &fInput);
	bool ReadTopology(ifstream &fInput);
	bool SetKinetics(OpenSMOKE_KPP_SingleReactor_KineticsManager* kinetics_);
	void UpdateInflow(const int k, const double convectionFlow, const double diffusionFlow);
	void Assembling();
	void ConvectionDiffusionMatrixCommunication(OpenSMOKE_KPP_ReactorNetwork& network);

	// Function which are called more than once
	void ODESystemContinousReactor(BzzVector &x, double t, BzzVector &dx);
	void GetJacobian(BzzVector &omega,BzzMatrix &JJ);

	void DistributeMassFlowRates(OpenSMOKE_KPP_ReactorNetwork& network);
	void SolveCSTR_CorrectorDiscrete(const double tau, BzzOdeStiffObject& o);
	void SolveCSTR_CorrectorContinous(const double tau, OpenSMOKE_KPP_ReactorNetwork& network, BzzOdeStiffObject& o);
	void SolveCSTR_Corrector_Linearized(const double deltat, BzzMatrix &tmpMatrix);
	void SolveCSTR_CorrectorDiscrete_Smart(const double tau, BzzMatrix& tmpMatrix, BzzOdeStiffObject& o);
	void SolveCSTR_CorrectorContinous_Smart(const double tau, OpenSMOKE_KPP_ReactorNetwork& network, BzzMatrix& tmpMatrix, BzzOdeStiffObject& o);

	void ReactionResiduals(BzzVector &residuals_);
	void Residuals(const int position, BzzVector &residuals, OpenSMOKE_KPP_ReactorNetwork& network);
	
	void ReconstructBlockAndDiagonals(BzzMatrix &block, BzzVector &diagonals);

	// Access functions (read-only)
	inline double volume() const { return volume_; }
	inline double temperature() const { return temperature_; }
	inline double pressure() const { return pressure_; }
	inline double variance() const { return variance_; }
	inline double M() { return (cOut_.GetSumElements() + fOut_ + diffusion_.GetSumElements()); }
	inline double mass() const { return mass_; }
	inline double fInTot()  { return fIn_.GetSumElements(); }
	inline double fOutTot() const { return fOut_; }
	inline unsigned int index() const { return index_; }
	inline const BzzVector& omega() const { return omega_; }
	inline const BzzVectorInt& in() const { return in_; }
	inline const BzzVectorInt& out() const { return out_; }
	inline const BzzVectorInt& neighbours() const { return neighbours_; }
	inline const BzzVector& cIn() const { return cIn_; }
	inline const BzzVector& cOut() const { return cOut_; }
	inline const BzzVector& diffusion() const { return diffusion_; }
	inline const BzzVector& fIn() const { return fIn_; }
	inline double RV(const int j) const { return RV_[j]; }
	inline const BzzVector& LocalRHS() const	{ return localRHS_; }
	inline const BzzVector& NonLocalRHS() const { return nonLocalRHS_; }
	inline double omegaMax()  { return omega_.Max(); }
	inline double omegaMin()  { return omega_.Min(); }

	inline int nGlobalSparsityPattern() const { return nGlobalSparsityPattern_; }
	inline int nGlobalSingleRowSparsityPattern() const { return nGlobalSingleRowSparsityPattern_; }

	inline const BzzVectorInt& globalIndicesSparsityRows() const { return globalIndicesSparsityRows_; }
//	inline const BzzVectorInt& globalIndicesSparsityColumns() const { return globalIndicesSparsityColumns_; }
	inline const BzzVector&    globalIndicesSparsityValues() const { return globalIndicesSparsityValues_; }

	inline const BzzVector& mConvectionDiffusion() const { return mConvectionDiffusion_; }
	inline const BzzVectorInt& iConvectionDiffusion() const { return iConvectionDiffusion_; }

	
	// Setting functions (write)
	inline void setcIn(const int j, const double val)  { cIn_[j]  = val; }
	inline void setcOut(const int j, const double val) { cOut_[j] = val; }
	inline void setfOutTot(const double val) { fOut_ = val; }
	inline void setOmega(const int j, const double val) { omega_[j] = val; }
	inline void setOmega(BzzVector &massFractions) { omega_ = massFractions; }
	       void setOmegaFromNetwork(const int k, const BzzVector &massFractions); 

	// Utilities
	double MassUmbalance();
	double MassFlowIn();
	double MassFlowOut();
	void UpdateOutputFlow();
	inline bool IsExternalOutputReactor()	const { return tagExternalOutput_; }
	inline bool IsExternalFeedReactor()		const { return tagExternalFeed_; }
	
//	void SetInitialConditions()	{	o.SetInitialConditions(omega_, 0., &odeSystem);	}

	void AssemblingLocalContribution(const double deltat, const bool jacobianFlag, BzzMatrix &tmpMatrix, BzzMatrix &diagonalBlockMatrix);
	void AssemblingNonLocalRHS(const double deltat, BzzVector &omegaNetwork);

private:

	OpenSMOKE_KPP_SingleReactor_KineticsManager* kinetics;
	OpenSMOKE_KPP_DataManager* data_;
//	BzzOdeStiffObject o;

	double volume_;
	double temperature_;
	double pressure_;
	double variance_;
	double mass_;

	unsigned int index_;
	int numberOfSpecies;

	BzzVector omega_;
	BzzVector R_;
	BzzVector RV_;
	BzzVector xMin_;
	BzzVector xMax_;

	BzzVector mInTot_;
	BzzVector mOut_x_omega_old_;

	BzzVectorInt in_;
	BzzVectorInt out_;
	BzzVectorInt neighbours_;

	BzzVector	cIn_;
	BzzVector	cOut_;
	BzzVector	diffusion_;
	BzzVector	fIn_;
	double		fOut_;

	double		 M_;
	
	BzzVectorInt iConvectionDiffusion_;
	BzzVector    mConvectionDiffusion_;

	BzzVectorInt jReduced_;

	bool tagExternalFeed_;
	bool tagExternalOutput_;

	// ODE solution of continously stirred reactor
//	MyOdeSystem_KPP_ContinousReactor odeSystem;

	// Solution of linear system arising from the Newton's Method
	BzzVector rhs_;
	BzzVector residuals_;

//	BzzMatrix diagonalBlockMatrix_;
	BzzVector localRHS_;
	BzzVector nonLocalRHS_;

	int nGlobalSparsityPattern_;
	int nGlobalSingleRowSparsityPattern_;
	BzzVectorInt globalIndicesSparsityRows_;
//	BzzVectorInt globalIndicesSparsityColumns_;
	BzzVector    globalIndicesSparsityValues_;

	void AssemblingDiagonalBlockMatrix(const double deltat, BzzMatrix &tmpMatrix, BzzMatrix &diagonalBlockMatrix);
	void AssemblingLocalRHS(const double deltat);

	void PrepareJacobianNewtonMethod(BzzFactorizedGauss& JacobianFactorized, BzzVector& omega, BzzMatrix &tmpMatrix);
	void CheckMassFractions(BzzVector& omega, const double epsilon);
	double NewtonStepReductionFactor(BzzVector& omega, BzzVector& direction);
	double ODESystemIntegration(const double tau, BzzVector& omega, BzzVector& residuals, BzzOdeStiffObject& o);

	BzzFactorizedGauss JacobianFactorizedGauss_;
	BzzVector d;
	BzzVector omega1;


public:
	struct
	{
		double norm1_over_nc;
		double normInf;
		singleReactorConvergenceType convergence;
		int nJacobianEvaluations;
		int nTotalCalls;
		int nNewtonIterations;
		int failure;
	}
	status;

	double Residuals(BzzVector& omega, BzzVector& residuals);

	void PrintDeFalco(const int position, BzzMatrixSparse C_, OpenSMOKE_KPP_ReactorNetwork& network, ofstream &fDeFalco, ofstream &fDeFalco2);

private:

	void ErrorMessage(const string message_);
	void WarningMessage(const string message_);
};


#endif