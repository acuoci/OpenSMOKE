#ifndef BZZ_MYODESYSTEMOBJECTONECSTR
#define BZZ_MYODESYSTEMOBJECTONECSTR

#include "BzzMath.hpp"

class myBzzCSTRNetwork;

class MyOdeSystemObjectOneCSTR : public BzzOdeSystemObject
{
private:
	
	int numComponents;
	myBzzCSTRNetwork *cstrNewtwork;

public:
	
	MyOdeSystemObjectOneCSTR(void);
	
	void operator()(myBzzCSTRNetwork *cstrN);

	virtual void GetJacobian(BzzVector &y,double t,BzzMatrix &JJ);
	virtual void GetSystemFunctions(BzzVector &y,double t,BzzVector &f);
	virtual void ObjectBzzPrint(void);
};

#endif // BZZ_MYODESYSTEMOBJECTONECSTR