#ifndef BZZ_MYODESYSTEMOBJECTALLCSTR
#define BZZ_MYODESYSTEMOBJECTALLCSTR

#include "BzzMath.hpp"

class myBzzCSTRNetwork;

class MyOdeSystemObjectAllCSTR : public BzzOdeSystemObject
{
private:
	int numComponents;
	myBzzCSTRNetwork *cstrNewtwork;

public:
	MyOdeSystemObjectAllCSTR(void);
	void operator()(myBzzCSTRNetwork *cstrN);

public:
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &f);
	virtual void ObjectBzzPrint(void);
	virtual void GetJacobian(BzzVector &yy,double tt);
	virtual void GetBuildJacobian(double hr,BzzMatrixSparseLockedByRows *Sh);
};

#endif

