#ifndef QWT_PLOT2D_H
#define QWT_PLOT2D_H

#include <QTextStream>
#include <qwt_plot.h>
#include <qwt_legend.h>
#include <qwt_plot_curve.h>
#include <qwt_magnifier.h>
#include <qwt_counter.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_panner.h>
#include "BzzMath.hpp"

class Qwt_Plot2D : public QwtPlot
{
	Q_OBJECT

public:

	Qwt_Plot2D(void);
	~Qwt_Plot2D(void);

private:

	int				nCurves;
	QwtPlotCurve	*curves;
	QPen			*color_pens;
	QPen			*bw_pens;

    QwtPlotZoomer	*d_zoomer[2];
    QwtPlotPicker	*d_picker;
    QwtPlotPanner	*d_panner;

	bool			iColor;

	bool *iEnabled;
	BzzVector *x_actual;
	BzzVector *y_actual;
	string *name_actual;

	static const int max_lines;

public:

	void setTitles(const string title, const string title_x, const string title_y);
	void enableSecondYAxis(const string title);
	void assignCurveToSecondYAxis(const int index);
	void setCurve(const int index, const QString name, BzzVector &x, BzzVector &y);
	void unsetCurve(const int index);

	void attachZoom();
	void attachPanner();
	void attachPicker();
	void enableZoomMode(bool on);

	void prepareText(QString *stringText);

public slots:

	void switchColor();
};

class Qwt_Plot2D_Zoomer: public QwtPlotZoomer
{
	public:
		Qwt_Plot2D_Zoomer(int xAxis, int yAxis, QwtPlotCanvas *canvas);
};

#endif