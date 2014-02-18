#include "Qwt_Plot2D.h"
#include <QMessageBox.h>
#include <qwt_plot_grid.h>

const int Qwt_Plot2D::max_lines = 10;

Qwt_Plot2D::Qwt_Plot2D(void)
{
	curves			= new QwtPlotCurve[max_lines];
	color_pens		= new QPen[max_lines];
	bw_pens			= new QPen[max_lines];
	iEnabled		= new bool[max_lines];
	x_actual		= new BzzVector[max_lines];
	y_actual		= new BzzVector[max_lines];
	name_actual		= new string[max_lines];

	for(int j=1;j<=max_lines;j++)
		iEnabled[j-1] = false;

	iColor = true;

	setCanvasBackground(QColor(Qt::white));
	insertLegend(new QwtLegend(), QwtPlot::RightLegend);
	
	// grid 
    QwtPlotGrid *grid = new QwtPlotGrid;
    grid->enableXMin(true);
    grid->setMajPen(QPen(Qt::black, 0, Qt::DotLine));
    grid->setMinPen(QPen(Qt::gray, 0 , Qt::DotLine));
    grid->attach(this);

	setMinimumSize(400, 250);

	color_pens[0] = QPen(Qt::red, 4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[1] = QPen(Qt::blue, 4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[2] = QPen(Qt::green, 4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[3] = QPen(Qt::yellow, 4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[4] = QPen(Qt::red, 4, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[5] = QPen(Qt::blue, 4, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[6] = QPen(Qt::green, 4, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[7] = QPen(Qt::yellow, 4, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[8] = QPen(Qt::green, 4, Qt::DotLine, Qt::RoundCap, Qt::RoundJoin);
	color_pens[9] = QPen(Qt::yellow, 4, Qt::DotLine, Qt::RoundCap, Qt::RoundJoin);

	bw_pens[0] = QPen(Qt::black, 4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	bw_pens[1] = QPen(Qt::gray, 4, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
	bw_pens[2] = QPen(Qt::black, 4, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
	bw_pens[3] = QPen(Qt::gray, 4, Qt::DashLine, Qt::RoundCap, Qt::RoundJoin);
	bw_pens[4] = QPen(Qt::black, 4, Qt::DotLine, Qt::RoundCap, Qt::RoundJoin);
	bw_pens[5] = QPen(Qt::gray, 4, Qt::DotLine, Qt::RoundCap, Qt::RoundJoin);
	bw_pens[6] = QPen(Qt::black, 4, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin);
	bw_pens[7] = QPen(Qt::gray, 4, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin);
}

Qwt_Plot2D::~Qwt_Plot2D(void)
{
}

void Qwt_Plot2D::setTitles(const string title, const string title_x, const string title_y)
{
	setTitle(title.c_str());
	setAxisTitle(QwtPlot::xBottom, title_x.c_str());
	setAxisTitle(QwtPlot::yLeft, title_y.c_str());
}

void Qwt_Plot2D::enableSecondYAxis(const string title)
{
	enableAxis(QwtPlot::yRight);
	setAxisTitle(QwtPlot::yRight, title.c_str());
}

void Qwt_Plot2D::assignCurveToSecondYAxis(const int index)
{
	curves[index-1].setYAxis(yRight);
}

void Qwt_Plot2D::attachZoom()
{
	d_zoomer[0] = new Qwt_Plot2D_Zoomer( QwtPlot::xBottom, QwtPlot::yLeft, this->canvas());
    d_zoomer[0]->setRubberBand(QwtPicker::RectRubberBand);
    d_zoomer[0]->setRubberBandPen(QColor(Qt::black));
    d_zoomer[0]->setTrackerMode(QwtPicker::ActiveOnly);
    d_zoomer[0]->setTrackerPen(QColor(Qt::black));

    d_zoomer[1] = new Qwt_Plot2D_Zoomer(QwtPlot::xTop, QwtPlot::yRight, this->canvas());
}

void Qwt_Plot2D::attachPanner()
{
	d_panner = new QwtPlotPanner(this->canvas());
    d_panner->setMouseButton(Qt::MidButton);
}

void Qwt_Plot2D::attachPicker()
{
    d_picker = new QwtPlotPicker(	QwtPlot::xBottom, QwtPlot::yLeft, 
															QwtPicker::PointSelection | QwtPicker::DragSelection, 
															QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn, 
															this->canvas());
    d_picker->setRubberBandPen(QColor(Qt::gray));
    d_picker->setRubberBand(QwtPicker::CrossRubberBand);
    d_picker->setTrackerPen(QColor(Qt::black));
}

void Qwt_Plot2D::setCurve(const int index, const QString name, BzzVector &x, BzzVector &y)
{
	iEnabled[index-1] = true;
	curves[index-1].setData(x.GetHandle(), y.GetHandle(), x.Size());
	curves[index-1].setPen(color_pens[index-1]);
	curves[index-1].setTitle(name);
	curves[index-1].attach(this);

	x_actual[index-1] = x;
	y_actual[index-1] = y;
	name_actual[index-1] = name.toStdString();
}

void Qwt_Plot2D::unsetCurve(const int index)
{
	iEnabled[index-1] = false;
	curves[index-1].detach();
}

void Qwt_Plot2D::enableZoomMode(bool on)
{
    d_panner->setEnabled(on);

    d_zoomer[0]->setEnabled(on);
    d_zoomer[0]->zoom(0);

    d_zoomer[1]->setEnabled(on);
    d_zoomer[1]->zoom(0);

    d_picker->setEnabled(!on);
}

void Qwt_Plot2D::switchColor()
{
	if (iColor == true)
	{
		iColor = false;
		for(int j=0;j<max_lines;j++)
			curves[j].setPen(bw_pens[j]);
	}
	else
	{
		iColor = true;
		for(int j=0;j<max_lines;j++)
			curves[j].setPen(color_pens[j]);
	}

	replot();
}

Qwt_Plot2D_Zoomer::Qwt_Plot2D_Zoomer(int xAxis, int yAxis, QwtPlotCanvas *canvas): QwtPlotZoomer(xAxis, yAxis, canvas)
{
    setSelectionFlags(QwtPicker::DragSelection | QwtPicker::CornerToCorner);
    setTrackerMode(QwtPicker::AlwaysOff);
    setRubberBand(QwtPicker::NoRubberBand);

    // RightButton: zoom out by 1
    // Ctrl+RightButton: zoom out to full size
    setMousePattern(QwtEventPattern::MouseSelect2, Qt::RightButton, Qt::ControlModifier);
    setMousePattern(QwtEventPattern::MouseSelect3, Qt::RightButton);
}

void Qwt_Plot2D::prepareText(QString *stringText)
{
	for(int i=1;i<max_lines;i++)
	{
		if (iEnabled[i] == true)
			if (x_actual[i].Size() != x_actual[i-1].Size())
			{
				QMessageBox msgBox;
				msgBox.setText("This option is not available for multiple curves with different number of support points.");
				msgBox.exec();
				return;
			}
	}

	QTextStream textStream(stringText);
	textStream.setRealNumberNotation(QTextStream::ScientificNotation);

	textStream << qSetFieldWidth(16) << left << endl;
	textStream << left  << qSetFieldWidth(16) << "x";
	for(int i=0;i<max_lines;i++)
		if (iEnabled[i] == true)
			textStream << qSetFieldWidth(16) << left << QString::fromStdString(name_actual[i]);
	textStream << endl;

	for(int j=1;j<=x_actual[0].Size();j++)
	{
		textStream << qSetFieldWidth(16) << left << x_actual[0][j];

		for(int i=0;i<max_lines;i++)
		{
			if (iEnabled[i] == true)
			{				
				textStream << qSetFieldWidth(16) << left << y_actual[i][j];
//				textStream << x_actual[i][j] << "\t";
//				textStream << y_actual[i][j] << "\t";

			}
		}
		textStream << endl;
	}
}