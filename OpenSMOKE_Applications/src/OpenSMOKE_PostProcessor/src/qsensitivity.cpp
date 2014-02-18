#include "qsensitivity.h"
#include "QDebugStream.h"
#include "qw_panel_chart2d.h"
#include "qw_panel_ropa.h"
#include "qw_panel_sensitivity.h"
#include "qt_opensmoke_about.h"
#include <QFileDialog>
#include <QTextEdit>
#include <iostream>

qsensitivity::qsensitivity(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);

	dialogOpen = new QFileDialog(this);
}

qsensitivity::~qsensitivity()
{

}
void qsensitivity::openFile()
{
	QFileDialog dialog(this);
	dialog.setFileMode(QFileDialog::ExistingFile);
	dialog.setNameFilter(tr("OpenSMOKE files (*.osm)"));
	dialog.setViewMode(QFileDialog::Detail);

	if (dialog.exec())
	{
		QStringList fileNames = dialog.selectedFiles();
		QString fileName = fileNames.at(0);

		char *fileName_char;
		QByteArray fileName_QByteArray = fileName.toLocal8Bit();
		fileName_char=fileName_QByteArray.data();
		string fileName_string = fileName_char;

		post_processor.ReadFromBinaryFile(fileName_string);

		if (post_processor.iChart() == true)
			ui.pushButton_chart2d->setEnabled(true);
		if (post_processor.iROPA() == true)
			ui.pushButton_ropa->setEnabled(true);
		if (post_processor.iSensitivity() == true)
			ui.pushButton_sensitivity->setEnabled(true);
		if (post_processor.iSensitivityDiffusivity() == true)
			ui.pushButton_sensitivity_diffusivity->setEnabled(true);
	}
}

void qsensitivity::openChart()
{
	vector<string> x_axis;
	vector<string> y_axis;

	post_processor.ExportAvailableXAxis(x_axis);
	post_processor.ExportAvailableYAxis(y_axis);

	QW_Panel_Chart2D *QwtChart2D;
	QwtChart2D = new QW_Panel_Chart2D();

	QwtChart2D->setParent(&post_processor);
	QwtChart2D->setKind(SPECIES);
	QwtChart2D->setup(x_axis, y_axis);
	QwtChart2D->show();
}

void qsensitivity::openROPA()
{
	vector<string> species_names;
	post_processor.ExportAllSpeciesNames(species_names);

	QW_Panel_ROPA *QwtROPA;
	QwtROPA = new QW_Panel_ROPA();

	QwtROPA->setParent(&post_processor);
	QwtROPA->setup(species_names);
	QwtROPA->show();
}

void qsensitivity::openSensitivity()
{
	const int index_sensitivity = 0;

	vector<string> species_names;
	vector<string> additional_names;
	post_processor.ExportSensitivitySpeciesNames(species_names, index_sensitivity);
	post_processor.ExportAdditionalNames(additional_names, index_sensitivity);

	QW_Panel_Sensitivity *QwtSensitivity;
	QwtSensitivity = new QW_Panel_Sensitivity();

	QwtSensitivity->setParent(&post_processor);
	QwtSensitivity->setup(species_names, additional_names, index_sensitivity);
	QwtSensitivity->show();
}

void qsensitivity::openSensitivityDiffusivity()
{
	const int index_sensitivity = 1;

	vector<string> species_names;
	vector<string> additional_names;
	post_processor.ExportSensitivitySpeciesNames(species_names, index_sensitivity);
	post_processor.ExportAdditionalNames(additional_names, index_sensitivity);

	QW_Panel_Sensitivity *QwtSensitivity;
	QwtSensitivity = new QW_Panel_Sensitivity();

	QwtSensitivity->setParent(&post_processor);
	QwtSensitivity->setup(species_names, additional_names, index_sensitivity);
	QwtSensitivity->show();
}

void qsensitivity::openAbout()
{
	Qt_OpenSMOKE_About *about;
	about = new Qt_OpenSMOKE_About();

	about->show();
}