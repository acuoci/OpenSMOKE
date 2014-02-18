#include <QtGui/QApplication>
//#include "qt_template.h"
#include "BzzMath.hpp"

#include "qbarswidget.h"
#include "qsensitivity.h"
#include "qwtplot2dwidget.h"
#include "qw_panel_chart2d.h"

#include <vector>


int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
//	QT_Template w;
//	wshow();

/*	std::vector<double> rect;
	rect.resize(12);
	rect[0] = 10.;	rect[1] = 7.;	rect[2] = -6.;
	rect[3] = 1.;	rect[4] = 0.1;	rect[5] = -0.2;
	rect[6] = -0.009;	rect[7] = -0.003;	rect[8] = -0.000005;
	rect[9] = 0.6509;	rect[10] = -0.0703;	rect[11] = -0.0000005;

	std::vector<int> irect;
	irect.resize(12);
	irect[0] = 5;	irect[1] = 14;		irect[2] = 2;
	irect[3] = 7;	irect[4] = 4;		irect[5] = 56;
	irect[6] = 2;	irect[7] = 9;		irect[8] = 47;
	irect[9] = 2;	irect[10] = 9;		irect[11] = 470;

	std::vector<QString> srect;
	srect.resize(12);
	srect[0] = "CO+H2O=CO2+CC";	srect[1] = "CO+H2O=CO2+CC";		srect[2] = "CO+H2O=CO2+CC";
	srect[3] = "CO+H2O=CO2+CC=CO2+CC";	srect[4] = "CO+H2O=CO2+CC=CO2+CC";		srect[5] = "CO+H2O=CO2+CC";
	srect[6] = "CO+H2O=CO2+CC";	srect[7] = "CO+H2O=CO2+CC";		srect[8] = "CO+H2O=CO2+CC";
	srect[9] = "CO+H2O=CO2+CC=CO2+CC";	srect[10] = "CO+H2O=CO2+CC";		srect[11] = "CO+H2O=CO2+CC=CO2+CC";

	QBarsWidget *bars;
	bars = new QBarsWidget();
	bars->setTitle("Sensitivity - CO");
	bars->setMainComment("Sensitivity analysis - CO mass fraction");
	bars->setRectangles(rect, irect, srect);
	bars->show();
*/
	qsensitivity *sensitivity;
	sensitivity = new qsensitivity();
	sensitivity->show();

/*	BzzMatrix x(5,100);
	BzzMatrix y(5,100);
	vector<QString> names;
	names.resize(5); 
	names[0] = "Line1";
	names[1] = "Line2";
	names[2] = "Line3";
	names[3] = "Line4";
	names[4] = "Line5";

	for(int j=1;j<=5;j++)
		for(int i=1;i<=100;i++)
		{	
			x[j][i] = i;
			y[j][i] = j + 1.34/double(j)*double(i);
		}

	QwtPlot2DWidget *QwtPlot2D;
	QwtPlot2D = new QwtPlot2DWidget();
//	QwtPlot2D->setup("distance from fuel nozzle [mm]", "mole fractions");
//	QwtPlot2D->setup(names, x, y);
//	QwtPlot2D->show();
 
/*	QW_Panel_Chart2D *QwtChart2D;
	QwtChart2D = new QW_Panel_Chart2D();
	std::vector<string> list_x_axis;
	std::vector<string> list_y_axis;
	list_x_axis.push_back("Axial position 1");
	list_x_axis.push_back("BAxial position 2");
	list_x_axis.push_back("CAxial position 3");

	list_y_axis.push_back("ZAxial position -1");
	list_y_axis.push_back("FAxial position -2");
	list_y_axis.push_back("XAxial position -3");
	list_y_axis.push_back("LAxial position +1");
	list_y_axis.push_back("SAxial position +2");
	list_y_axis.push_back("Axial position +3");

	QwtChart2D->setup(list_x_axis, list_y_axis);
	QwtChart2D->show();
*/	
	//DataPlot *plot = new DataPlot(&w);
	//plot->replot();

	return a.exec();
}