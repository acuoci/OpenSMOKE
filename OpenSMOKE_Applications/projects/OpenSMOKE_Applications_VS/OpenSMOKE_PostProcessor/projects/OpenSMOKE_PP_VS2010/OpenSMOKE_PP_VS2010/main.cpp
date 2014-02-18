#include "opensmoke_pp_vs2010.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	OpenSMOKE_PP_VS2010 w;
	w.show();
	return a.exec();
}
