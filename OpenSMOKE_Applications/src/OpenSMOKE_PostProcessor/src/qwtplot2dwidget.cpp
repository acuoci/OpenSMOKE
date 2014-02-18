#include "qwtplot2dwidget.h"
#include <qwt_scale_engine.h>
#include <QPainter>
#include <QClipboard>
#include <QFileDialog>
#include <QPrintDialog>
#include <QPixmap>
#include <QTextDocument>

QwtPlot2DWidget::QwtPlot2DWidget(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);

	ui.horizontalLayout->addWidget(&plot2d);
	plot2d.attachZoom();
	plot2d.attachPicker();
	plot2d.attachPanner();

	iLogX		= false;
	iLogY		= false;
	iNormalX	= false;
	iNormalY	= false;

	copyImageAction = new QAction(tr("&Copy Image"), this);
	copyImageAction->setStatusTip(tr("Copy plot as image"));
	copyImageAction->setIcon(QIcon(":/Resources/images/icon_copy.png"));
	copyImageAction->setShortcut(QKeySequence::Copy);
	connect(copyImageAction, SIGNAL(triggered()), this, SLOT(copyImage()));

	printImageAction = new QAction(tr("&Print Image"), this);
	printImageAction->setStatusTip(tr("Print plot as image"));
	printImageAction->setIcon(QIcon(":/Resources/images/icon_pen.png"));
	printImageAction->setShortcut(QKeySequence::Print);
	connect(printImageAction, SIGNAL(triggered()), this, SLOT(printImage()));

	saveImageAction = new QAction(tr("&Save Image"), this);
	saveImageAction->setStatusTip(tr("Save plot as image"));
	saveImageAction->setIcon(QIcon(":/Resources/images/icon_save.png"));
	saveImageAction->setShortcut(QKeySequence::Save);
	connect(saveImageAction, SIGNAL(triggered()), this, SLOT(saveImage()));

	printTextAction = new QAction(tr("Print Text"), this);
	printTextAction->setStatusTip(tr("Print plot as text"));
	printTextAction->setIcon(QIcon(":/Resources/images/icon_pen.png"));
	connect(printTextAction, SIGNAL(triggered()), this, SLOT(printText()));

	copyTextAction = new QAction(tr("Copy Text"), this);
	copyTextAction->setStatusTip(tr("Copy plot as text"));
	copyTextAction->setIcon(QIcon(":/Resources/images/icon_pen.png"));
	connect(copyTextAction, SIGNAL(triggered()), this, SLOT(copyText()));

	addAction(copyImageAction);
	addAction(printImageAction);
	addAction(saveImageAction);
	addAction(copyTextAction);
	addAction(printTextAction);
	setContextMenuPolicy(Qt::ActionsContextMenu);
}

QwtPlot2DWidget::~QwtPlot2DWidget()
{

}

void QwtPlot2DWidget::setup(const QString x_axis, const QString y_axis)
{
	plot2d.setAxisTitle(QwtPlot::xBottom,	x_axis);
	plot2d.setAxisTitle(QwtPlot::yLeft,		y_axis);
}

void QwtPlot2DWidget::setup(const vector<QString> _names, BzzMatrix &_x, BzzMatrix &_y)
{
	number_of_curves = _x.Rows();
	
	names = _names;
	x = _x;
	y = _y;

	for(int i=1;i<=number_of_curves;i++)
	{
		min_x.Append(x.GetRow(i).Min());
		max_x.Append(x.GetRow(i).Max());
		min_y.Append(y.GetRow(i).Min());
		max_y.Append(y.GetRow(i).Max());
	}

	for(int i=1;i<=number_of_curves;i++)
		plot2d.setCurve(i, names[i-1], x.GetRow(i), y.GetRow(i));
	plot2d.replot();
}

void QwtPlot2DWidget::setAutoscale()
{
	plot2d.setAxisAutoScale(QwtPlot::yLeft);
	plot2d.setAxisAutoScale(QwtPlot::xBottom);
	plot2d.replot();
}

void QwtPlot2DWidget::zoomIn(QwtScaleDiv *ax, double &x1, double &x2)
{
	x1 = ax->lowerBound();
	x2 = ax->upperBound();
	double xmean = (x1+x2)/2.;
	double delta = (x2-x1);
	x1 = xmean - delta/4.;
	x2 = xmean + delta/4.;
}

void QwtPlot2DWidget::zoomOut(QwtScaleDiv *ax, double &x1, double &x2)
{
	x1 = ax->lowerBound();
	x2 = ax->upperBound();
	double xmean = (x1+x2)/2.;
	double delta = (x2-x1);
	x1 = xmean - delta;
	x2 = xmean + delta;
}


void QwtPlot2DWidget::setZoomIn()
{
	double x1  = 0;
	double x2  = 0;

	{
		QwtScaleDiv *ax = plot2d.axisScaleDiv(QwtPlot::yLeft);
		zoomIn(ax, x1, x2);
		plot2d.setAxisScale(QwtPlot::yLeft, x1, x2);
	}

	{
		QwtScaleDiv *ax = plot2d.axisScaleDiv(QwtPlot::xBottom);
		zoomIn(ax, x1, x2);
		plot2d.setAxisScale(QwtPlot::xBottom, x1, x2);
	}

	plot2d.replot();
}

void QwtPlot2DWidget::setZoomOut()
{
	double x1  = 0;
	double x2  = 0;

	{
		QwtScaleDiv *ax = plot2d.axisScaleDiv(QwtPlot::yLeft);
		zoomOut(ax, x1, x2);
		plot2d.setAxisScale(QwtPlot::yLeft, x1, x2);
	}

	{
		QwtScaleDiv *ax = plot2d.axisScaleDiv(QwtPlot::xBottom);
		zoomOut(ax, x1, x2);
		plot2d.setAxisScale(QwtPlot::xBottom, x1, x2);
	}

	plot2d.replot();
}

void QwtPlot2DWidget::setLogX()
{
	if (iLogX == false)
	{
		plot2d.setAxisScaleEngine(QwtPlot::xBottom, new QwtLog10ScaleEngine());
		plot2d.replot();
		iLogX = true;
	}
	else
	{
		plot2d.setAxisScaleEngine(QwtPlot::xBottom, new QwtLinearScaleEngine());
		plot2d.replot();
		iLogX = false;
	}
}

void QwtPlot2DWidget::setLogY()
{
	if (iLogY == false)
	{
		plot2d.setAxisScaleEngine(QwtPlot::yLeft, new QwtLog10ScaleEngine());
		plot2d.replot();
		iLogY = true;
	}
	else
	{
		plot2d.setAxisScaleEngine(QwtPlot::yLeft, new QwtLinearScaleEngine());
		plot2d.replot();
		iLogY = false;
	}
}

void QwtPlot2DWidget::setNormalX()
{
	setNormal('x');
}

void QwtPlot2DWidget::setNormalY()
{
	setNormal('y');
}

void QwtPlot2DWidget::setNormal(char axis)
{
	if (axis == 'x')	
	{
		     if (iNormalX == true)	iNormalX = false;
		else if (iNormalX == false)	iNormalX = true;
	}

	else if (axis == 'y')	
	{
		     if (iNormalY == true)	iNormalY = false;
		else if (iNormalY == false)	iNormalY = true;
	}


	if (iNormalX == false && iNormalY == false)
	{
		for(int i=1;i<=number_of_curves;i++)
			plot2d.setCurve(i, names[i-1], x.GetRow(i), y.GetRow(i));
	}

	else if (iNormalX == true && iNormalY == false)
	{
		for(int i=1;i<=number_of_curves;i++)
		{
			BzzVector xnormal(x.GetRow(i).Size());
			for(int j=1;j<=xnormal.Size();j++)
				xnormal[j] = (x[i][j]-min_x[i])/(max_x[i]-min_x[i]); 
			plot2d.setCurve(i, names[i-1], xnormal, y.GetRow(i));
		}
	}
	
	else if (iNormalX == false && iNormalY == true)
	{
		for(int i=1;i<=number_of_curves;i++)
		{
			BzzVector ynormal(y.GetRow(i).Size());
			for(int j=1;j<=ynormal.Size();j++)
				ynormal[j] = (y[i][j]-min_y[i])/(max_y[i]-min_y[i]); 
			plot2d.setCurve(i, names[i-1], x.GetRow(i), ynormal);
		}
	}

	else if (iNormalX == true && iNormalY == true)
	{
		for(int i=1;i<=number_of_curves;i++)
		{
			BzzVector xnormal(x.GetRow(i).Size());
			BzzVector ynormal(y.GetRow(i).Size());
			for(int j=1;j<=ynormal.Size();j++)
			{
				xnormal[j] = (x[i][j]-min_x[i])/(max_x[i]-min_x[i]); 
				ynormal[j] = (y[i][j]-min_y[i])/(max_y[i]-min_y[i]); 
			}
			plot2d.setCurve(i, names[i-1], xnormal, ynormal);
		}
	}

	plot2d.replot();
}

void QwtPlot2DWidget::setColors()
{
	plot2d.switchColor();
}

void QwtPlot2DWidget::saveImage()
{
	int options = QwtPlotPrintFilter::PrintAll;
	options &= ~QwtPlotPrintFilter::PrintBackground;
	options |= QwtPlotPrintFilter::PrintFrameWithScales;
    
	QString fileName = QFileDialog::getSaveFileName(this, tr("File name"), QString(), "Graphic files (*.svg,*png)");
	QPixmap plotPxm(plot2d.width(), plot2d.height());
	plotPxm.fill(Qt::white);
	plot2d.print(plotPxm);
	plotPxm.save(fileName, "PNG");
}

void QwtPlot2DWidget::copyImage()
{
	int options = QwtPlotPrintFilter::PrintAll;
	options &= ~QwtPlotPrintFilter::PrintBackground;
	options |= QwtPlotPrintFilter::PrintFrameWithScales;
    
	QPixmap plotPxm(plot2d.width(), plot2d.height());
	plotPxm.fill(Qt::white);
	plot2d.print(plotPxm);

	QImage image = plotPxm.toImage();
	QApplication::clipboard()->setImage(image);
}

void QwtPlot2DWidget::printImage()
{
	QPrintDialog printDialog(&printer, this);
	if (printDialog.exec())
	{
		int options = QwtPlotPrintFilter::PrintAll;
		options &= ~QwtPlotPrintFilter::PrintBackground;
		options |= QwtPlotPrintFilter::PrintFrameWithScales;
    
		QPixmap plotPxm(plot2d.width(), plot2d.height());
		plotPxm.fill(Qt::white);
		plot2d.print(plotPxm);

		QImage image = plotPxm.toImage();

		QPainter painterPrinter(&printer);

		painterPrinter.setViewport( painterPrinter.window().width()/2.-plot2d.width()/2., 
									painterPrinter.window().height()/10, 
									plot2d.width(), plot2d.height() );

		QRect rect(0,0,plot2d.width(),plot2d.height());
		painterPrinter.drawImage(rect,image);
	}
}

void QwtPlot2DWidget::prepareText(QString *stringText)
{
	QTextStream textStream(stringText);
	
}

void QwtPlot2DWidget::copyText()
{
	QString textStream;
	plot2d.prepareText(&textStream);
	QApplication::clipboard()->setText(textStream);
}


void QwtPlot2DWidget::printText()
{
	QString textStream;
	plot2d.prepareText(&textStream);

	QPrintDialog printDialog(&printer, this);
	if (printDialog.exec())
	{
		QTextDocument textDocument;
		textDocument.setPlainText(textStream);
		textDocument.print(&printer);
	}
}