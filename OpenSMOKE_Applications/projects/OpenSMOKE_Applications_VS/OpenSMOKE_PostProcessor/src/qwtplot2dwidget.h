#ifndef QWTPLOT2DWIDGET_H
#define QWTPLOT2DWIDGET_H

#include <QWidget>
#include <QPrinter>
#include <QTextStream>
#include "GeneratedFiles/ui_qwtplot2dwidget.h"
#include "Qwt_Plot2D.h"

class QwtPlot2DWidget : public QWidget
{
	Q_OBJECT

public:
	QwtPlot2DWidget(QWidget *parent = 0);
	~QwtPlot2DWidget();

private:

	Ui::QwtPlot2DWidgetClass ui;

	Qwt_Plot2D plot2d;

	int number_of_curves;
	vector<QString> names;

	bool iLogX;
	bool iLogY;
	bool iNormalX;
	bool iNormalY;

	BzzMatrix x;
	BzzMatrix y;

	BzzVector min_x;
	BzzVector max_x;
	BzzVector min_y;
	BzzVector max_y;

	void setNormal(char axis);
	void zoomIn(QwtScaleDiv *ax, double &x1, double &x2);
	void zoomOut(QwtScaleDiv *ax, double &x1, double &x2);


	QAction *copyImageAction;
	QAction *printImageAction;
	QAction *saveImageAction;
	QAction *printTextAction;
	QAction *copyTextAction;

	QPrinter printer;

	void prepareText(QString *stringText);


public:

	void setup(const QString x_axis, const QString y_axis);
	void setup(const vector<QString> names, BzzMatrix &x, BzzMatrix &y);

private slots:

	void setAutoscale();
	void setZoomIn();
	void setZoomOut();
	void setLogX();
	void setLogY();
	void setNormalX();
	void setNormalY();
	void setColors();
	
	void copyImage();
	void printImage();
	void saveImage();
	void copyText();
	void printText();
};

#endif // QWTPLOT2DWIDGET_H
