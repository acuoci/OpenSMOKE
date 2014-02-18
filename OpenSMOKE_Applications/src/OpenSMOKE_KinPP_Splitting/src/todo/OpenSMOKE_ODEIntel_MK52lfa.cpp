#include "OpenSMOKE_ODEIntel_MK52lfa.h"


void ResidualsQ(int *n_, double *t_, double *y_, double *f_)
{
	cout << "ResS" << endl;
	double c;
	c=1.e0-y_[0]*y_[0];
	f_[0]=y_[1];
	f_[1]=(c*y_[1]-y_[0])*1.0e6;
	cout << t_ << " " << c << " " << f_[0] << " " << f_[1] << " " << y_[0] << " " << y_[1] << endl;
	cout << "ResE" << endl;
}

void JacobianQ(int *n_, double *t_, double *y_, double *a_)
{
	cout << "JacS" << endl;
	a_[0]=0.e0;
	a_[1]=-1.e6*(1.e0+2.e0*y_[0]*y_[1]);
	a_[2]=1.e0;
	a_[3]=1.e6*(1.e0-y_[0]*y_[0]);
	cout << t_ << " " << a_[0] << " " << a_[1] << " " << a_[2] << " " << a_[3] << endl;
	cout << "JacE" << endl;
}