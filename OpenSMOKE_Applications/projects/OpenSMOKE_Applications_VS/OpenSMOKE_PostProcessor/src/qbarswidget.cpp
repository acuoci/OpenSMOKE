#include "qbarswidget.h"
#include <QPainter>
#include <QClipboard>
#include <QPrintDialog>
#include <QTextStream>
#include <QImage>
#include <QTextDocument>
#include <math.h>
#include <iostream>

// Maximum number of bars
const int		QBarsWidget::max_bars = 20;

// Logical coordinates
const double	QBarsWidget::h_window = 1000;
const double	QBarsWidget::w_window = 1000;


QBarsWidget::QBarsWidget(QWidget *parent) : QWidget(parent)
{
	ui.setupUi(this);

	setBackgroundRole(QPalette::Light);

	// Colors
	color_PlusBars  = Qt::red;
	color_MinusBars = Qt::blue;

	// Set context menu
	setContextMenu();
}

void QBarsWidget::setContextMenu()
{
	copyImageAction = new QAction(tr("&Copy Image"), this);
	copyImageAction->setStatusTip(tr("Copy bars as image"));
	copyImageAction->setIcon(QIcon(":/Resources/images/icon_copy.png"));
	copyImageAction->setShortcut(QKeySequence::Copy);
	connect(copyImageAction, SIGNAL(triggered()), this, SLOT(copyImage()));

	printImageAction = new QAction(tr("&Print Image"), this);
	printImageAction->setStatusTip(tr("Print bars as image"));
	printImageAction->setIcon(QIcon(":/Resources/images/icon_pen.png"));
	printImageAction->setShortcut(QKeySequence::Print);
	connect(printImageAction, SIGNAL(triggered()), this, SLOT(printImage()));

	copyTextAction = new QAction(tr("Copy Text"), this);
	copyTextAction->setStatusTip(tr("Copy bars as text"));
	copyTextAction->setIcon(QIcon(":/Resources/images/icon_copy.png"));
	connect(copyTextAction, SIGNAL(triggered()), this, SLOT(copyText()));

	printTextAction = new QAction(tr("Print Text"), this);
	printTextAction->setStatusTip(tr("Print bars as text"));
	printTextAction->setIcon(QIcon(":/Resources/images/icon_pen.png"));
	connect(printTextAction, SIGNAL(triggered()), this, SLOT(printText()));

	colorAction = new QAction(tr("Change Color"), this);
	colorAction->setStatusTip(tr("Change bar colors"));
	colorAction->setIcon(QIcon(":/Resources/images/icon_colors.png"));
	connect(colorAction, SIGNAL(triggered()), this, SLOT(setColors()));

	addAction(copyImageAction);
	addAction(printImageAction);
	addAction(copyTextAction);
	addAction(printTextAction);
	addAction(colorAction);
	setContextMenuPolicy(Qt::ActionsContextMenu);
}

void QBarsWidget::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);
	draw(&painter);
}

void QBarsWidget::draw(QPainter *painter)
{
	painter->setRenderHint(QPainter::Antialiasing, true);

	// Set logical coordinate
	painter->setWindow(-w_window/2.,0, h_window,w_window);

	// Set background color
	painter->setBrush(QBrush(Qt::white, Qt::SolidPattern));
	painter->setPen(QPen(Qt::white, 0.5, Qt::SolidLine, Qt::RoundCap, Qt::MiterJoin));
	painter->drawRect(-w_window/2., 0, h_window, w_window);

	// Set title
	painter->setPen(QPen(Qt::black, 0.5, Qt::SolidLine, Qt::RoundCap, Qt::MiterJoin));
	painter->setFont(QFont("Helvetica", h_window/10./4., QFont::Normal));
	painter->drawText(-w_window/2.,0.,w_window,h_window/10., Qt::AlignVCenter| Qt::AlignCenter, mainComment);

	// Set rectangles
	for(int i=0;i<n_rectangles;i++)
	{
		if (x_coordinates[i] >= 0.)	
		{
			painter->setPen(QPen(color_PlusBars, 0.5, Qt::SolidLine, Qt::RoundCap, Qt::MiterJoin));
			painter->setBrush(QBrush(color_PlusBars, Qt::SolidPattern));

		}
		if (x_coordinates[i] <  0.)	
		{
			painter->setPen(QPen(color_MinusBars, 0.5, Qt::SolidLine, Qt::RoundCap, Qt::MiterJoin));
			painter->setBrush(QBrush(color_MinusBars, Qt::SolidPattern));
		}
		painter->drawRect(0, y_coordinates[i], x_coordinates[i], h_rectangles);

		painter->setFont(font_index);
		painter->drawText(	-w_window*0.50, y_coordinates[i], 0.10*w_window, h_rectangles, 
							Qt::AlignVCenter | Qt::AlignCenter,	r_indices[i]);

		painter->setFont(font_value);
		painter->drawText(	0.40*w_window, y_coordinates[i], 0.10*w_window, h_rectangles, 
							Qt::AlignVCenter | Qt::AlignRight,	r_values[i]);

		if (x_coordinates[i] >= 0.)
		{
			painter->setFont(font_string);
			painter->drawText(	-w_window*0.42, y_coordinates[i], w_window*0.40, h_rectangles, 
								Qt::AlignVCenter| Qt::AlignRight,	r_strings[i]);

		}

		if (x_coordinates[i]  < 0.)
		{
			painter->setFont(font_string);
			painter->drawText(	w_window*0.02, y_coordinates[i], w_window*0.40, h_rectangles, 
								Qt::AlignVCenter| Qt::AlignLeft,	r_strings[i]);

		}
	}
}

void QBarsWidget::copyImage()
{
	QImage image(width(), height(), QImage::Format_RGB32);
	QPainter painterImage(&image);
	draw(&painterImage);

	QApplication::clipboard()->setImage(image);
}


void QBarsWidget::printImage()
{
	QPrintDialog printDialog(&printer, this);
	if (printDialog.exec())
	{
		QPainter painterPrinter(&printer);
		painterPrinter.setViewport(	painterPrinter.window().width()/2.-width()/2., 
									painterPrinter.window().height()/10, 
									width(),height() );
		draw(&painterPrinter);
	}
}

void QBarsWidget::prepareText(QString *stringText)
{
	QTextStream textStream(stringText);

	int max_size_string = 0;
	for(int i=0;i<n_rectangles;i++)
		if (r_strings[i].size()>max_size_string)	max_size_string = r_strings[i].size();
	if (max_size_string > 100) max_size_string = 100;

	for(int i=0;i<n_rectangles;i++)
	{
		textStream << left  << qSetFieldWidth(7) << i+1;
		textStream << left  << qSetFieldWidth(max_size_string+2) << r_strings[i];
		textStream << right << qSetFieldWidth(16) << r_values[i];
		textStream << "\n";
	}
}

void QBarsWidget::copyText()
{
	QString stringText;
	prepareText(&stringText);
	QApplication::clipboard()->setText(stringText);
}


void QBarsWidget::printText()
{
	QString stringText;
	prepareText(&stringText);

	QPrintDialog printDialog(&printer, this);
	if (printDialog.exec())
	{
		QTextDocument textDocument;
		textDocument.setPlainText(stringText);
		textDocument.print(&printer);
	}
}

void QBarsWidget::setRectangles(const std::vector<double>		_r_coordinates,
								const std::vector<int>		    _indices,
								const std::vector<std::string>  _names)
{
	std::vector<QString> _namesQ;
	for(int j=0;j<_names.size();j++)
		_namesQ.push_back(QString::fromStdString(_names[j]));
	
	setRectangles(_r_coordinates, _indices, _namesQ); 
}

void QBarsWidget::setRectangles(const std::vector<double>   _r_coordinates,
								const std::vector<int>      _indices,
								const std::vector<QString>  _names) 
{
	r_coordinates		= _r_coordinates;
	r_strings			= _names;
	n_total_rectangles	= r_coordinates.size();
	n_rectangles		= qMin(max_bars, n_total_rectangles);

	if (n_rectangles < max_bars)
		resize(width(), height()*double(n_rectangles)/double(max_bars));

	h_rectangles = double(0.90*h_window)/double(n_rectangles) * 0.90;

	x_coordinates.resize(n_rectangles);
	y_coordinates.resize(n_rectangles);
	r_indices.resize(n_rectangles);
	r_values.resize(n_rectangles);
	
	for(int i=0;i<n_rectangles;i++)
	{
		x_coordinates[i] = r_coordinates[i]/fabs(r_coordinates[0]) * (0.40*w_window);
		y_coordinates[i] = h_window*0.10+(1.1*h_rectangles)*i;

		r_indices[i].setNum(_indices[i]); 

		if (fabs(r_coordinates[i]) > 1e-2)	r_values[i].setNum(r_coordinates[i], 'f', 3); 
		else r_values[i].setNum(r_coordinates[i], 'e', 2); 
	}

	setFonts();
}

void QBarsWidget::setFonts()
{
	QString max_string = NULL;
	int max_string_lenght = 0;
	for(int i=0;i<n_rectangles;i++)
		if (r_strings[i].size() > max_string_lenght)
		{
			max_string = r_strings[i];
			max_string_lenght = r_strings[i].size();
		}

	int count = 0;
	int w;
	do
	{
		font_string = QFont("Helvetica", h_rectangles/(3.+double(count++)*0.1), QFont::Normal);
		QFontMetrics fm(font_string);
		w = fm.width(max_string)/double(width())*double(w_window);
	} 
	while (w>=0.40*w_window && count < 10);

	font_value = QFont("Helvetica", h_rectangles/4.5, QFont::Normal);
	font_index = QFont("Helvetica", h_rectangles/4.5, QFont::Normal);
}

void QBarsWidget::setColors()
{
	color_PlusBars = Qt::black;
	color_MinusBars = Qt::gray;
	update();
}

void QBarsWidget::setTitle(const QString _title)
{
	setWindowTitle(_title);
}

void QBarsWidget::setMainComment(const QString _comment)
{
	mainComment = _comment;
}


QBarsWidget::~QBarsWidget()
{

}
