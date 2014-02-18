#include "MyOdeSystemObjectOneCSTR.h"
#include "myBzzCSTRNetwork.hpp"


MyOdeSystemObjectOneCSTR::MyOdeSystemObjectOneCSTR(void)
{
	cstrNewtwork = 0;
}

void MyOdeSystemObjectOneCSTR::operator()(myBzzCSTRNetwork *cstrN)
{
	cstrNewtwork = cstrN;
	numComponents = cstrNewtwork->Reactions.NumberOfSpecies();
}

void MyOdeSystemObjectOneCSTR::ObjectBzzPrint(void)
{
	::BzzPrint("\nObject Print for Numerical Jacobian");
	::BzzPrint("\nReactor %d", cstrNewtwork->kReactor);
}

void MyOdeSystemObjectOneCSTR::GetSystemFunctions(BzzVector &x,double t, BzzVector &f)
{
	cstrNewtwork->GetResiduals(x,f);
	Minus(&f);
}

void MyOdeSystemObjectOneCSTR:: GetJacobian(BzzVector &x,double t,BzzMatrix &JJ)
{
	cstrNewtwork->GetJacobian(x,JJ);
}