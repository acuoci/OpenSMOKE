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

#include <sstream>
#include "OpenSMOKE_KPP_SingleReactor.h"
#include "OpenSMOKE_KPP_ReactorNetwork.h"
#include "OpenSMOKE_KPP_DataManager.h"
#include "OpenSMOKE_KPP_SingleReactor_KineticsManager.h"

void OpenSMOKE_KPP_SingleReactor::ErrorMessage(const string message_)
{
    cout << endl;
    cout << "Class: OpenSMOKE_KPP_SingleReactor"	<< endl;
    cout << "Index: " << index_						<< endl;
    cout << "Error: " << message_					<< endl;
    cout << "Press enter to continue... "			<< endl;
    getchar();
    exit(-1);
}

void OpenSMOKE_KPP_SingleReactor::WarningMessage(const string message_)
{
    cout << endl;
    cout << "Class:	  OpenSMOKE_KPP_SingleReactor"	<< endl;
    cout << "Index:   " << index_						<< endl;
    cout << "Warning: " << message_						<< endl;
	cout << endl;
}

OpenSMOKE_KPP_SingleReactor::OpenSMOKE_KPP_SingleReactor()
{
	tagExternalOutput_	  = false;
	tagExternalFeed_	  = false;
	status.nTotalCalls	  = 0;
}

void OpenSMOKE_KPP_SingleReactor::ResetStatus()
{
	status.nTotalCalls	  = 0;
}

void OpenSMOKE_KPP_SingleReactor::Setup(const int i, OpenSMOKE_ReactingGas* mixture, BzzVectorInt& indexReduced)
{
	index_ = i;
	numberOfSpecies = mixture->NumberOfSpecies();
	jReduced_ = indexReduced;
//	odeSystem.assignReactor(this);
}

void OpenSMOKE_KPP_SingleReactor::SetDataManager(OpenSMOKE_KPP_DataManager* data)
{
	data_ = data;
}

OpenSMOKE_KPP_SingleReactor::~OpenSMOKE_KPP_SingleReactor(void)
{
}

bool OpenSMOKE_KPP_SingleReactor::ReadMassFractions(ifstream &fInput)
{
	ChangeDimensions(numberOfSpecies, &omega_);
	ChangeDimensions(numberOfSpecies, &R_);
	ChangeDimensions(numberOfSpecies, &RV_);
	ChangeDimensions(numberOfSpecies, &xMin_);
	ChangeDimensions(numberOfSpecies, &xMax_);
	ChangeDimensions(numberOfSpecies, &mInTot_);
	ChangeDimensions(numberOfSpecies, &mOut_x_omega_old_);
	ChangeDimensions(numberOfSpecies, &rhs_);
	ChangeDimensions(numberOfSpecies, &residuals_);

	ChangeDimensions(numberOfSpecies, &localRHS_);
	ChangeDimensions(numberOfSpecies, &nonLocalRHS_);

	ChangeDimensions(numberOfSpecies, &d);
	ChangeDimensions(numberOfSpecies, &omega1);


	// Max and Min for ODE integration
	xMin_ = 0.;
	xMax_ = 1.;

	double sum = 0.;
	for (int j=1; j<=jReduced_.Size(); j++)
	{
		fInput >> omega_[jReduced_[j]];
		sum += omega_[jReduced_[j]];
	}

	double epsilon = 1.e-2;
	if ( sum > (1.+epsilon) || sum < (1.-epsilon) )
	{
		stringstream sum_string;
		sum_string << sum;
		ErrorMessage("The sum of mass fractions is not equal to 1: " + sum_string.str());
	}
	else 
	{
		for (int j=1; j<=jReduced_.Size(); j++)
			omega_[jReduced_[j]] /= sum;
	}

	return true;
}

bool OpenSMOKE_KPP_SingleReactor::ReadReactorProperties(ifstream &fInput)
{
	unsigned int dummy_int;

	fInput >> dummy_int;

	if (dummy_int != index_)
		ErrorMessage("Indices do not fit!");

	fInput >> temperature_;
	fInput >> pressure_;
	fInput >> volume_;
	fInput >> variance_;

	pressure_ *= 101325.;
	
	return true;
}

bool OpenSMOKE_KPP_SingleReactor::SetKinetics(OpenSMOKE_KPP_SingleReactor_KineticsManager* kinetics_)
{
	kinetics = kinetics_;
	return true;
}

bool OpenSMOKE_KPP_SingleReactor::ReadTopology(ifstream &fInput)
{
	// External feeds
	{
		unsigned int numberInlets;
		fInput >> numberInlets;
		ChangeDimensions(numberOfSpecies, &fIn_);
		for (int k=1; k<=int(numberInlets); k++)
		{
			tagExternalFeed_ = true;

			unsigned int numberSpecies;
			double massFlow;
			fInput >> massFlow;
			fInput >> numberSpecies;
			for (unsigned int j=1; j<=numberSpecies; j++)
			{
				unsigned int indexSpecies;
				double omegaSpecies;
				fInput >> indexSpecies >> omegaSpecies;
				fIn_[indexSpecies] += omegaSpecies*massFlow;
			}
		}
	}

	// Outputs
	{
		unsigned int numberOutputs;
		fInput >> numberOutputs;

		ChangeDimensions(numberOutputs, &cOut_);
		ChangeDimensions(numberOutputs, &out_);
		ChangeDimensions(numberOutputs, &neighbours_);
		ChangeDimensions(numberOutputs, &diffusion_);

		for (int k=1; k<=out_.Size(); k++)
		{
			fInput >> out_[k] >> cOut_[k] >> diffusion_[k];
			neighbours_[k] = out_[k];
		}
	}

	return true;
}

void OpenSMOKE_KPP_SingleReactor::UpdateInflow(const int k, const double convectionFlow, const double diffusionFlow)
{
	in_.Append(k);
	cIn_.Append(convectionFlow);

	neighbours_.Append(k);
	diffusion_.Append(diffusionFlow);
}

void OpenSMOKE_KPP_SingleReactor::Assembling()
{
	// Update external outputs
	{
		BzzVectorInt iRemove;

		fOut_ = 0.;
		for (int j=1;j<=out_.Size();j++)
			if (out_[j] == 0)
			{
				tagExternalOutput_ = true;
				fOut_ += cOut_[j];
				iRemove.Append(j);
			}

		if (tagExternalOutput_ == true)
		{
			out_.DeleteElements(iRemove);
			cOut_.DeleteElements(iRemove);
			diffusion_.DeleteElements(iRemove);
			neighbours_.DeleteElements(iRemove);
		}
	}

	// Composition
	double cTot_ = pressure_/Constants::R_J_kmol/temperature_;
	double mw_ = kinetics->GetMWFromMassFractions(omega_);
	mass_ = cTot_*mw_*volume_;
}

void OpenSMOKE_KPP_SingleReactor::UpdateOutputFlow()
{
	if (tagExternalOutput_ == true)
		fOut_ = MassFlowIn() - cOut_.GetSumElements();
}

double OpenSMOKE_KPP_SingleReactor::MassUmbalance()
{
	return (cIn_.GetSumElements()+fIn_.GetSumElements()) - (cOut_.GetSumElements()+fOut_); 
}

double OpenSMOKE_KPP_SingleReactor::MassFlowIn()
{
	return (cIn_.GetSumElements()+fIn_.GetSumElements()); 
}

double OpenSMOKE_KPP_SingleReactor::MassFlowOut()
{
	return (cOut_.GetSumElements()+fOut_); 
}

void OpenSMOKE_KPP_SingleReactor::GetJacobian(BzzVector &omega,BzzMatrix &JJ)
{
	// Get dR_over_domega
	kinetics->GetFormationRatesDerivatives(omega, temperature_, pressure_, JJ);

	// JJ = dR_over_domega * volume
	JJ *= volume_;

	// Add the transport coefficient
	// This term MUST be added if we are in the sequential procedure
	if (data_->PredictorCorrector_DeferredConvection() == true || data_->networkStatus() == KPP_NETWORK_STATUS_SEQUENTIAL_CSTR )
		if (data_->networkStatus() != KPP_NETWORK_STATUS_GLOBALODE && data_->networkStatus() != KPP_NETWORK_STATUS_GLOBALNLS)
		{
			for(int i=1;i<=numberOfSpecies;i++)	
				JJ[i][i] -= M_;
		}

	// Scaling
	JJ /= mass_;
}

void OpenSMOKE_KPP_SingleReactor::ODESystemContinousReactor(BzzVector &omega, double t, BzzVector &domega)
{
	if (data_->PredictorCorrector_DeferredConvection() == true || 
		data_->networkStatus() == KPP_NETWORK_STATUS_SEQUENTIAL_CSTR )
	{
		// Formation rates
		kinetics->UpdateProperties(omega, temperature_, pressure_, R_);
		Product(volume_, R_, &RV_);

		// mass*domegaj = -mOut*omegaj + mjIn + Rj*V
		Product(-M_, omega, &domega);
		Sum(&domega, RV_);
		Sum(&domega, mInTot_);
		domega /= mass_;
	}
	else
	{
		// Formation rates
		kinetics->UpdateProperties(omega, temperature_, pressure_, R_);
		Product(volume_, R_, &RV_);

		// mass*domegaj = -mOut*omegaj(n) + mjIn + Rj*V
		domega = RV_;
		Sum(&domega, mInTot_);
		Difference(&domega, mOut_x_omega_old_);
		domega /= mass_;
	}
}

void OpenSMOKE_KPP_SingleReactor::ReactionResiduals(BzzVector &residuals)
{
	kinetics->UpdateProperties(omega_, temperature_, pressure_, R_);
	Product(volume_, R_, &RV_);
	residuals = RV_; 
}

void OpenSMOKE_KPP_SingleReactor::Residuals(const int position, BzzVector &globalResiduals, OpenSMOKE_KPP_ReactorNetwork& network)
{
	// Reactions
	kinetics->UpdateProperties(omega_, temperature_, pressure_, R_);
	Product(volume_, R_, &residuals_);

	// External input feeds
	if (tagExternalFeed_ == true)	residuals_ += fIn_;

	// Inflow
	for(int j=1;j<=iConvectionDiffusion_.Size();j++)
		residuals_ += mConvectionDiffusion_[j] * network.reactors(iConvectionDiffusion_[j]).omega();

	// Outflow
	residuals_ -= M_*omega_;

	// From local to global residuals
	globalResiduals.SetBzzVector(position, residuals_); 
}

void OpenSMOKE_KPP_SingleReactor::ConvectionDiffusionMatrixCommunication(OpenSMOKE_KPP_ReactorNetwork& network)
{
	int k;
	ElementBzzMatrixSparse *elem;
		
	elem = network.C().GetStartingElementInRow(index_);
	k = 0;
	while(elem)
	{
		k++;
		elem = elem->next;
	}

	ChangeDimensions(k-1, &iConvectionDiffusion_);
	ChangeDimensions(k-1, &mConvectionDiffusion_);

	elem = network.C().GetStartingElementInRow(index_);
	k = 0;
	while(elem)
	{
		int j = elem->column;
		if (j == index_)
		{
			M_ = elem->value;
		}
		else
		{
			k++;
			iConvectionDiffusion_[k] = j;
			mConvectionDiffusion_[k] = -elem->value;
		}
		elem = elem->next;
	}
	
	// Global sparsity structure
	nGlobalSparsityPattern_ = numberOfSpecies * (numberOfSpecies+iConvectionDiffusion_.Size());
	nGlobalSingleRowSparsityPattern_ = numberOfSpecies+iConvectionDiffusion_.Size();
//	ChangeDimensions(nGlobalSparsityPattern_, &globalIndicesSparsityColumns_);
	ChangeDimensions(nGlobalSparsityPattern_, &globalIndicesSparsityValues_);
	
	// Global columns for row block corresponding to the current reactor)
	{
		int jReactor = (index_-1)*numberOfSpecies;

		int count=1;
		for (int i=1;i<=numberOfSpecies;i++)
		{
	//		for (int k=1;k<=mConvectionDiffusion_.Size();k++)
	//			if (iConvectionDiffusion_[k] < int(index_))
	//				globalIndicesSparsityColumns_[count++] = (iConvectionDiffusion_[k]-1)*numberOfSpecies+i;

	//		for (int k=1;k<=numberOfSpecies;k++)
	//			globalIndicesSparsityColumns_[count++] = jReactor+k;

	//		for (int k=1;k<=mConvectionDiffusion_.Size();k++)
	//			if (iConvectionDiffusion_[k] > int(index_))
	//				globalIndicesSparsityColumns_[count++] = (iConvectionDiffusion_[k]-1)*numberOfSpecies+i;

			for (int k=1;k<=mConvectionDiffusion_.Size();k++)
				if (iConvectionDiffusion_[k] == int(index_))
					ErrorMessage("Wrong sparsity structure");
		}
	}
}

void OpenSMOKE_KPP_SingleReactor::DistributeMassFlowRates(OpenSMOKE_KPP_ReactorNetwork& network)
{
	// External feeds
	if (tagExternalFeed_ == true)	mInTot_ = fIn_;
	else							mInTot_ = 0.;

	// Inflow
	for(int j=1;j<=iConvectionDiffusion_.Size();j++)
		mInTot_ += mConvectionDiffusion_[j] * network.reactors(iConvectionDiffusion_[j]).omega();
}

// This function is called when the network is updated continously
void OpenSMOKE_KPP_SingleReactor::SolveCSTR_CorrectorContinous(const double tau, OpenSMOKE_KPP_ReactorNetwork& network,  BzzOdeStiffObject& o)
{
	// External feeds
	if (tagExternalFeed_ == true)	mInTot_ = fIn_;
	else							mInTot_ = 0.;

	// Inflow
	for(int j=1;j<=iConvectionDiffusion_.Size();j++)
		mInTot_ += mConvectionDiffusion_[j] * network.reactors(iConvectionDiffusion_[j]).omega();
	
	// Solution
	SolveCSTR_CorrectorDiscrete(tau, o);
}


void OpenSMOKE_KPP_SingleReactor::SolveCSTR_CorrectorDiscrete(const double tau, BzzOdeStiffObject& o)
{
	// If deferredConvection is off && the sequential procedure is off
	if (data_->PredictorCorrector_DeferredConvection() == false && 
		data_->networkStatus() != KPP_NETWORK_STATUS_SEQUENTIAL_CSTR )
		Product(M_, omega_, &mOut_x_omega_old_);

	// Solving ODE System
	{
//		o.SetInitialConditions(omega_, 0., &odeSystem);
		o.SetInitialConditions(omega_, 0.);
		o.SetAnalyticalJacobian();									
		o.SetMinimumConstraints(xMin_);
		o.SetMaximumConstraints(xMax_);
		o.SetTolAbs(data_->SingleReactor_OdeAbsoluteTolerance());
		
		// Max number of Jacobians calls
		if (data_->SingleReactor_OdeMaxJacobian() > 0)
			o.StopIntegrationBeforeRecalcuatingJacobian(data_->SingleReactor_OdeMaxJacobian());
		
		// Convergence rule
		if (data_->SingleReactor_OdeStopResiduals() > 0.)
		{
			double maxOdeSum_ = .1 * data_->SingleReactor_OdeStopResiduals() * double(numberOfSpecies);
			o.StopIntegrationWhenSumAbsY1IsLessThan(maxOdeSum_);
		}

		if(o.GetOdeCalculationState() < 0)
			ErrorMessage("ODE System (Continous Reactor): Calculation State <0");

		omega_ = o(tau);
	}
}

void OpenSMOKE_KPP_SingleReactor::SolveCSTR_Corrector_Linearized(const double deltat, BzzMatrix &tmpMatrix)
{
	// Reactions	
	kinetics->UpdateProperties(omega_, temperature_, pressure_, R_);
	Product(volume_, R_, &RV_);
	GetJacobian(omega_, tmpMatrix);				// A = d(RV)/domega * Volume/mtot;

	// The formulation is the same, both for deferred=on and deferred=off
	// The difference is the Jacobian calculation (see GetJacobian function)
	if (data_->PredictorCorrector_DeferredConvection() == false || 
		data_->PredictorCorrector_DeferredConvection() == true)
	{
		// Rigth hand side
		Product(-M_, omega_, &rhs_);		//	outflow
		Sum(&rhs_, RV_);				//	reaction
		Sum(&rhs_, mInTot_);			//  inflow
		rhs_ *= deltat/mass_;		
	
		// Matrix
		tmpMatrix *= -deltat;			
		for(int i=1;i<=numberOfSpecies;i++)
			tmpMatrix[i][i] += 1.;
	
		// Linear System solution
		BzzFactorizedGauss AGauss(tmpMatrix);
		Solve(AGauss, &rhs_);

		// Updating
		omega_ += rhs_;
	}
}

void OpenSMOKE_KPP_SingleReactor::AssemblingLocalContribution(const double deltat, const bool jacobianFlag, BzzMatrix &tmpMatrix, BzzMatrix &diagonalBlockMatrix)
{
	// Reactions	
	kinetics->UpdateProperties(omega_, temperature_, pressure_, R_);
	Product(volume_, R_, &RV_);

	// Local RHS
	AssemblingLocalRHS(deltat);

	// Block Matrix
	if (jacobianFlag == true)
		AssemblingDiagonalBlockMatrix(deltat, tmpMatrix, diagonalBlockMatrix);
}

void OpenSMOKE_KPP_SingleReactor::AssemblingDiagonalBlockMatrix(const double deltat, BzzMatrix &tmpMatrix, BzzMatrix &diagonalBlockMatrix)
{	
	if (data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALODE)
	{
		// Jacobian	
		GetJacobian(omega_, tmpMatrix);				// A =  d(RV)/domega /mtot;
		Product(-deltat, &tmpMatrix);				// A = -d(RV)/domega /mtot*deltat;

		// DiagonalBlockMatrix
		diagonalBlockMatrix = 0.;
		diagonalBlockMatrix.SetDiagonal(0, 1.+M_*deltat/mass_);
		Sum(&diagonalBlockMatrix, tmpMatrix);

		// Global values for row block corresponding to the current reactor)
		{
			BzzVector mConvectionDiffusion_x_deltat_over_mass = mConvectionDiffusion_;
			Product(deltat/mass_, &mConvectionDiffusion_x_deltat_over_mass);

			int count=1;
			for (int i=1;i<=numberOfSpecies;i++)
			{
				for (int k=1;k<=mConvectionDiffusion_.Size();k++)
					if (iConvectionDiffusion_[k] < int(index_))
						globalIndicesSparsityValues_[count++] = -mConvectionDiffusion_x_deltat_over_mass[k];

				for (int k=1;k<=numberOfSpecies;k++)
					globalIndicesSparsityValues_[count++] = diagonalBlockMatrix[i][k];

				for (int k=1;k<=mConvectionDiffusion_.Size();k++)
					if (iConvectionDiffusion_[k] > int(index_))
						globalIndicesSparsityValues_[count++] = -mConvectionDiffusion_x_deltat_over_mass[k];
			}
		}

	}
	else if (data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALNLS)
	{
		// Jacobian	
		GetJacobian(omega_, tmpMatrix);				// A =  d(RV)/domega /mtot;
		Product(-mass_, &tmpMatrix);

		// DiagonalBlockMatrix
		diagonalBlockMatrix = 0.;
		diagonalBlockMatrix.SetDiagonal(0, M_);
		Sum(&diagonalBlockMatrix, tmpMatrix);

		// Global values for row block corresponding to the current reactor)
		{

			int count=1;
			for (int i=1;i<=numberOfSpecies;i++)
			{
				for (int k=1;k<=mConvectionDiffusion_.Size();k++)
					if (iConvectionDiffusion_[k] < int(index_))
						globalIndicesSparsityValues_[count++] = -mConvectionDiffusion_[k];

				for (int k=1;k<=numberOfSpecies;k++)
					globalIndicesSparsityValues_[count++] = diagonalBlockMatrix[i][k];

				for (int k=1;k<=mConvectionDiffusion_.Size();k++)
					if (iConvectionDiffusion_[k] > int(index_))
						globalIndicesSparsityValues_[count++] = -mConvectionDiffusion_[k];
			}
		}
	}
	else
		ErrorMessage("AssemblingDiagonalBlockMatrix is available only for KINPP_NETWORK_STATUS_GLOBALODE or KINPP_NETWORK_STATUS_GLOBALNLS");		
}

void OpenSMOKE_KPP_SingleReactor::AssemblingLocalRHS(const double deltat)
{
	if (data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALODE)
	{
		Product(-M_, omega_, &localRHS_);		
		Sum(&localRHS_, RV_);	
		if (tagExternalFeed_ == true)
			Sum(&localRHS_, fIn_);			
		localRHS_ *= deltat/mass_;		
	}
	else if (data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALNLS)
	{
		Product(-M_, omega_, &localRHS_);		
		Sum(&localRHS_, RV_);	
		if (tagExternalFeed_ == true)
			Sum(&localRHS_, fIn_);				
	}
}

void OpenSMOKE_KPP_SingleReactor::AssemblingNonLocalRHS(const double deltat, BzzVector &omegaNetwork)
{
	// External feeds
	nonLocalRHS_ = 0.;
	for(int j=1;j<=iConvectionDiffusion_.Size();j++)
	{
		int position = (iConvectionDiffusion_[j]-1)*numberOfSpecies;
		for(int i=1;i<=numberOfSpecies;i++)
			nonLocalRHS_[i] += mConvectionDiffusion_[j] * omegaNetwork[position+i];
	}

	if (data_->networkStatus() == KPP_NETWORK_STATUS_GLOBALODE)
		nonLocalRHS_ *= deltat/mass_;
}

void OpenSMOKE_KPP_SingleReactor::setOmegaFromNetwork(const int k, const BzzVector &massFractions)
{
	int position = (k-1)*numberOfSpecies;
	for(int i=1;i<=numberOfSpecies;i++)
		omega_[i] = massFractions[position+i];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//																								   //
//									ODE SYSTEM - CLASS DEFINITION								   //
//																								   //
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MyOdeSystem_KPP_ContinousReactor::ObjectBzzPrint(void)
{
}

void MyOdeSystem_KPP_ContinousReactor::GetSystemFunctions(BzzVector &x, double t, BzzVector &dx)
{
	ptReactor->ODESystemContinousReactor(x, t, dx);
}

void MyOdeSystem_KPP_ContinousReactor::assignReactor(OpenSMOKE_KPP_SingleReactor *reactor)
{
	ptReactor = reactor;
}

void MyOdeSystem_KPP_ContinousReactor::GetJacobian(BzzVector &x,double t,BzzMatrix &JJ)
{
	ptReactor->GetJacobian(x,JJ);
}

// This function is called when the network is updated continously
void OpenSMOKE_KPP_SingleReactor::SolveCSTR_CorrectorContinous_Smart(const double tau, OpenSMOKE_KPP_ReactorNetwork& network, BzzMatrix &tmpMatrix, BzzOdeStiffObject& o)
{
	// External feeds
	if (tagExternalFeed_ == true)	mInTot_ = fIn_;
	else							mInTot_ = 0.;

	// Inflow
	for(int j=1;j<=iConvectionDiffusion_.Size();j++)
		mInTot_ += mConvectionDiffusion_[j] * network.reactors(iConvectionDiffusion_[j]).omega();
	
	// Solution
	SolveCSTR_CorrectorDiscrete_Smart(tau, tmpMatrix, o);
}

void OpenSMOKE_KPP_SingleReactor::SolveCSTR_CorrectorDiscrete_Smart(const double tau, BzzMatrix& tmpMatrix, BzzOdeStiffObject& o)
{
	// Variables
	double FInitial, F0, F1;

	double eps					= 1.e-6;
	
	// Should be optimized
	double hMinRequested		= 0.01;			// Minimum reduction coefficient to apply Newton
	double requestedF0			= 1.e-16;		// A reactor is ok if its residuals is lower
	double maxSumF				= 1.e-13;		// The newton method is ok if the residual is lower
	int maxIterationsNewtons	= 15;			// Maximum number of Newton iterations
	int maxJacobianAge			= 5;			// Frequency to update the jacobian during the Newton method
	double maxRatio				= 0.75;			// Max reduction factor in the newton method
	int maxSubIterations		= 4;			// Pseudo-Armijio steps
	double safetyReductionCoefficient = 0.90;

	// Status
	status.nJacobianEvaluations = 0;
	status.failure = 0;

	// 1. Checking initial composition
	CheckMassFractions(omega_, eps);

	// 2. If the network was never solved, the ODE solver is called 
	if (status.nTotalCalls == 0)
	{
		F1 = ODESystemIntegration(tau, omega_, residuals_, o);
		if (F1 > 1.e-4)
			status.failure = 1;
		status.norm1_over_nc = F1;
		status.normInf		 = residuals_.MaxAbs();
		status.convergence   = KPP_SINGLEREACTOR_CONVERGENCE_ODE_FIRST;
		status.nTotalCalls++;
		return;
	}

	// 3. Initial residuals
	F0 = Residuals(omega_, residuals_);
	FInitial = F0;

	// 4. Check if this reactor is already OK
	if(F0 < requestedF0)
	{
		status.norm1_over_nc = F0;
		status.normInf		 = residuals_.MaxAbs();
		status.convergence   = KPP_SINGLEREACTOR_CONVERGENCE_DIRECT;
		status.nTotalCalls++;
		return;
	}

	// 4. First step of newton method
	{
		PrepareJacobianNewtonMethod(JacobianFactorizedGauss_, omega_, tmpMatrix);
		status.nJacobianEvaluations++;
	}

	// 5a. Newton's correction
	Solve(JacobianFactorizedGauss_, residuals_, &d);
	double hMin = NewtonStepReductionFactor(omega_, d);

	// 6. Solve ODE system if Newton method is not satisfactory
	if(hMin<=hMinRequested)
	{
		F1 = ODESystemIntegration(tau, omega_, residuals_, o);
		if (F1 > 1.e-4)
			status.failure = 1;
		status.norm1_over_nc = F1;
		status.normInf = residuals_.MaxAbs();
		status.convergence = KPP_SINGLEREACTOR_CONVERGENCE_ODE_FIRST;
		status.nTotalCalls++;
		return;
	}

	// 7. Again with Newton's Method
	if(hMin != 1.)
	{
		hMin *= safetyReductionCoefficient;
		Product(hMin,&d);
	}
	Sum(omega_, d, &omega1);

	// 8. Checking mass fractions
	CheckMassFractions(omega1, eps);

	// 9. Update residuals
	F1 = Residuals(omega1, residuals_);

	// 10. ODE (II)
	if (F1>=F0)
	{
		F1 = ODESystemIntegration(tau, omega_, residuals_, o);
		status.norm1_over_nc = F1;
		status.normInf = residuals_.MaxAbs();
		status.convergence   = KPP_SINGLEREACTOR_CONVERGENCE_ODE_SECOND;
		status.nTotalCalls++;
		return;
	}

	int jacobianAge = 0;
	for(int iteration=1;iteration<=maxIterationsNewtons;iteration++)
	{
		jacobianAge++;

		// Check convergence
		if(F1 < maxSumF)
		{
			omega_ = omega1;
			status.norm1_over_nc = F1;
			status.normInf = residuals_.MaxAbs();
			status.convergence   = KPP_SINGLEREACTOR_CONVERGENCE_NEWTON;
			status.nNewtonIterations = iteration;
			status.nTotalCalls++;
			return;
		}

		// Updating Solution
		if (iteration>1)
		{
			if (jacobianAge == maxJacobianAge)
			{
				PrepareJacobianNewtonMethod(JacobianFactorizedGauss_, omega_, tmpMatrix);
				jacobianAge = 0;
				status.nJacobianEvaluations++;
			}

			Solve(JacobianFactorizedGauss_, residuals_, &d);
			double hMin = NewtonStepReductionFactor(omega_, d);

			if (hMin<hMinRequested)
			{
				F1 = ODESystemIntegration(tau, omega_, residuals_, o);
				status.norm1_over_nc = F1;
				status.normInf = residuals_.MaxAbs();
				status.convergence   = KPP_SINGLEREACTOR_CONVERGENCE_ODE_THIRD;
				status.nTotalCalls++;
				return;
			}

			if(hMin != 1.)
			{
				hMin *= safetyReductionCoefficient;
				Product(hMin,&d);
			}
			Sum(omega_, d, &omega1);

			CheckMassFractions(omega1, eps);
			F1 = Residuals(omega1, residuals_);
		}

		// If reduction is not satisfactory
		int subIterations = 1;
		bool flagSubIterations = true;
		while (F1/F0 > maxRatio)
		{
			if (subIterations == maxSubIterations)
			{
				flagSubIterations = false;
				break;
			}

			Product(0.50,&d);
			Sum(omega_, d, &omega1);
			
			CheckMassFractions(omega1, eps);
			F1 = Residuals(omega1, residuals_);

			subIterations++;
		}

		if (flagSubIterations == true)
		{
			omega_ = omega1;
			F0 = F1;
		}
		else
		{
			// If the jacobian is new, then it is better to perform a ODE integration
			if (jacobianAge == 0)
			{
				F1 = ODESystemIntegration(tau, omega_, residuals_, o);
				status.norm1_over_nc = F1;
				status.normInf = residuals_.MaxAbs();
				status.convergence   = KPP_SINGLEREACTOR_CONVERGENCE_ODE_FOURTH;
				status.nTotalCalls++;
				return;
			}
			// else before performing the ODE integration is better to perform an additional attempt
			else
			{
				PrepareJacobianNewtonMethod(JacobianFactorizedGauss_, omega_, tmpMatrix);
				jacobianAge = -1;
				status.nJacobianEvaluations++;
			}
		}
	}

}

void OpenSMOKE_KPP_SingleReactor::CheckMassFractions(BzzVector& omega, const double epsilon)
{
	double one_plus_epsilon  = 1.+epsilon;
	double one_minus_epsilon = 1.-epsilon;

	double sum = omega.GetSumElements();
	
	if(sum > (1.+1.e-3) || sum < (1.-1.e-3))
	{
		cout << "Error in sum of mass fractions: " << sum-1. << " (+-" << epsilon << ")" << endl;
		ErrorMessage("Fatal Error");
	}

	if(sum > one_plus_epsilon || sum < one_minus_epsilon)
		status.failure = 1;

	sum = 1./sum;
	Product(sum, &omega);
}

void OpenSMOKE_KPP_SingleReactor::PrepareJacobianNewtonMethod(BzzFactorizedGauss& JacobianFactorized, BzzVector& omega, BzzMatrix &tmpMatrix)
{
		GetJacobian(omega, tmpMatrix);				// A = d(RV)/domega * Volume/mtot;
		Product(mass_, &tmpMatrix);
		for(int i=1;i<=numberOfSpecies;i++)
			tmpMatrix[i][i] -= M_;
		JacobianFactorized = tmpMatrix;
}

double OpenSMOKE_KPP_SingleReactor::Residuals(BzzVector& omega, BzzVector& residuals)
{
	kinetics->UpdateProperties(omega, temperature_, pressure_, R_);
	Product(volume_, R_, &RV_);

	Product(-M_, omega, &residuals);		
	Sum(&residuals, RV_);				
	Sum(&residuals, mInTot_);	

	Product(-1, &residuals);

	double F = residuals.GetSumAbsElements() / double(numberOfSpecies);

	return F;
}

double OpenSMOKE_KPP_SingleReactor::NewtonStepReductionFactor(BzzVector& omega, BzzVector& direction)
{
	double hMin = 1.;
	for(int i=1;i<=numberOfSpecies;i++)
	{
		if(direction[i] < 0.)  
		{
			if(omega[i] == 0.)
			{
				if(direction[i] < -1.e-15)
				{
					hMin = 0.;
					break;
				}
				else
					direction[i] = 0.;
			}
			else
			{
				hMin = min(hMin, -omega[i]/direction[i]);
			}
		}
	}

	return hMin;
}

double OpenSMOKE_KPP_SingleReactor::ODESystemIntegration(const double tau, BzzVector& omega, BzzVector& residuals,  BzzOdeStiffObject& o)
{
	o.SetInitialConditions(omega, 0.);
	o.SetAnalyticalJacobian();									
	o.SetMinimumConstraints(xMin_);
	o.SetMaximumConstraints(xMax_);
	o.SetTolAbs(data_->SingleReactor_OdeAbsoluteTolerance());
		
	// Max number of Jacobians calls
	if (data_->SingleReactor_OdeMaxJacobian() > 0)
		o.StopIntegrationBeforeRecalcuatingJacobian(data_->SingleReactor_OdeMaxJacobian());
		
	// Convergence rule
	if (data_->SingleReactor_OdeStopResiduals() > 0.)
	{
		double maxOdeSum_ = .1 * data_->SingleReactor_OdeStopResiduals() * double(numberOfSpecies);
		o.StopIntegrationWhenSumAbsY1IsLessThan(maxOdeSum_);
	}

	// Solve
	omega = o(tau);

	// Checking solution
	if(o.GetOdeCalculationState() < 0)
		ErrorMessage("ODE System (Continous Reactor): Calculation State < 0");

	// Residui (termini di accumulo)
	residuals = o.GetY1InMeshPoint();
	double F = residuals.GetSumAbsElements() / double(numberOfSpecies);

	return F;
}

void OpenSMOKE_KPP_SingleReactor::ReconstructBlockAndDiagonals(BzzMatrix &block, BzzVector &diagonals)
{
	ChangeDimensions(iConvectionDiffusion_.Size(), &diagonals);

	int count=1;
	for (int i=1;i<=numberOfSpecies;i++)
	{
		int iDiagonal = 1;
		for (int k=1;k<=mConvectionDiffusion_.Size();k++)
			if (iConvectionDiffusion_[k] < int(index_))
			{
				if (diagonals[iDiagonal]!=0 && diagonals[iDiagonal] != globalIndicesSparsityValues_[count])
					ErrorMessage("Something wrong!");

				diagonals[iDiagonal++] = globalIndicesSparsityValues_[count++];
			}

		for (int k=1;k<=numberOfSpecies;k++)
			block[i][k] = globalIndicesSparsityValues_[count++];
				
		for (int k=1;k<=mConvectionDiffusion_.Size();k++)
			if (iConvectionDiffusion_[k] > int(index_))
			{
				diagonals[iDiagonal++] = globalIndicesSparsityValues_[count++];
			}
	}
}

void OpenSMOKE_KPP_SingleReactor::PrintDeFalco(const int position, BzzMatrixSparse C_, OpenSMOKE_KPP_ReactorNetwork& network, ofstream &fDeFalco, ofstream &fDeFalco2)
{
	// Reactions
	kinetics->UpdateProperties(omega_, temperature_, pressure_, R_);
	
	// Get dR_over_domega
	BzzMatrix JJ(omega_.Size(),omega_.Size());
	kinetics->GetFormationRatesDerivatives(omega_, temperature_, pressure_, JJ);

	// JJ = dR_over_domega * volume
	JJ *= volume_;

	fDeFalco << "Cell           " << position << endl;
	fDeFalco << "MassFractions  " << omega_.Size() << endl;
	for(int j=1;j<=omega_.Size();j++)
		fDeFalco << setw(4) << left << j << setw(16) << left << omega_[j] << endl;
	fDeFalco << "FormationRates " << omega_.Size() << endl;
	for(int j=1;j<=omega_.Size();j++)
		fDeFalco << setw(4) << left << j << setw(16) << left << R_[j]*volume_ << endl;
	fDeFalco << "Jacobian       " << omega_.Size() << " " << omega_.Size() << endl;
	for(int j=1;j<=omega_.Size();j++)
		for(int k=1;k<=omega_.Size();k++)
			fDeFalco << setw(4) << left << j << setw(4) << left << k << setw(16) << left << JJ[j][k] << endl;

	fDeFalco2 << "Cell           " << position << endl;
	for(int j=1;j<=omega_.Size();j++)
	{
		double* ptrVal;
		int i, k;
		double val;

		double sum_total_out = 0.;
		double sum_total = 0.;
		C_.BeginScanning();
		while(ptrVal = C_.Scanning(&i,&k,&val))
		{
			if (i==position && k!=position)
				sum_total_out += -C_(position, k)*network.reactors(k).omega()[j];
//			if (i==position && k!=position)
//				sum_total += -C_(position, k);

		}

		fDeFalco2 << setw(4) << left << j << setw(16) << left << omega_[j];
//		fDeFalco2 <<                         setw(16) << left << -C_(position, position)*omega_[j];
//		fDeFalco2 <<                         setw(16) << left << sum_total_out;
		fDeFalco2 <<                         setw(16) << left << -C_(position, position)*omega_[j] + sum_total_out;
		fDeFalco2 <<                         setw(16) << left << fIn_[j];
		fDeFalco2 <<                         setw(16) << left <<  R_[j]*volume_;
		fDeFalco2 <<                         setw(16) << left <<  R_[j]*volume_ + fIn_[j] + sum_total_out -C_(position, position)*omega_[j];
//		fDeFalco2 <<                         setw(16) << left <<  -C_(position, position);
//		fDeFalco2 <<                         setw(16) << left <<  sum_total;
		fDeFalco2 << endl;
	}
}
