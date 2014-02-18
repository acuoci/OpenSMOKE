#include <sstream>
#include "qbarswidget.h"
#include "qw_panel_sensitivity.h"
#include "qwtplot2dwidget.h"
#include "qw_panel_chart2d.h"
#include "addons/OpenSMOKE_PostProcessor.h"

QW_Panel_Sensitivity::QW_Panel_Sensitivity(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);
}

QW_Panel_Sensitivity::~QW_Panel_Sensitivity()
{

}

void QW_Panel_Sensitivity::setParent(OpenSMOKE_PostProcessor *post_processor_)
{
	post_processor = post_processor_;
	
	BzzVector x;
	x = post_processor->get_x();
	min_coordinate = x[1];
	max_coordinate = x[x.Size()];
}

void QW_Panel_Sensitivity::setup(std::vector<std::string> list_names, std::vector<std::string> list_additional, const int index_sensitivity)
{
	index_sensitivity_ = index_sensitivity;

	for(int j=0;j<list_names.size();j++)
		ui.comboBox_SensitivityItem_Species->addItem(list_names[j].c_str());

	for(int j=0;j<list_additional.size();j++)
		ui.comboBox_SensitivityItem_Additional->addItem(list_additional[j].c_str());

	QString coordinate; 
	coordinate.setNum(0.50*(min_coordinate+max_coordinate));
	ui.lineEdit_LocalAnalysis->setText(coordinate);
}

void QW_Panel_Sensitivity::change_species_sensitivity()
{
	ui.comboBox_SensitivityItem_Species->setEnabled(true);
	ui.comboBox_SensitivityItem_Additional->setEnabled(false);
	ui.radioButton_LocalNormalization->setEnabled(false);
	ui.radioButton_LocalNormalization->setChecked(false);
	ui.radioButton_GlobalNormalization->setChecked(true);
}

void QW_Panel_Sensitivity::change_additional_sensitivity()
{
	ui.comboBox_SensitivityItem_Species->setEnabled(false);
	ui.comboBox_SensitivityItem_Additional->setEnabled(true);
	ui.radioButton_LocalNormalization->setEnabled(true);
}
#include <QMessageBox>
void QW_Panel_Sensitivity::display_sensitivityitem()
{
	vector<double> t;
	vector<int> it;
	vector<string> names_t;
	int index = ui.comboBox_SensitivityItem_Species->currentIndex()+1;

	// Post Processing
	bool iLocalAnalysis = ui.checkBox_LocalAnalysis->isChecked();
	bool iGlobalNormalization = ui.radioButton_GlobalNormalization->isChecked();

	double coordinate = ui.lineEdit_LocalAnalysis->text().toDouble();


	// TO ADJUST
	bool iTotal = true;
	bool iLocal = true;
	if (iLocalAnalysis == true)			iTotal = false;
	if (iGlobalNormalization == true)	iLocal = false;

	if ( (coordinate < min_coordinate || coordinate > max_coordinate) && (iTotal == false))
	{
		QMessageBox msgBox;
		stringstream coordinate_string;
		coordinate_string << coordinate;
		QString coordinate_qstring = QString::fromStdString("Coordinate outside boundaries: " + coordinate_string.str());
		msgBox.setText(coordinate_qstring);
		int ret = msgBox.exec();
	}
	else
	{
		QString label;
	
		if (ui.radioButton_ComboBox_Additional->isChecked() == true)
		{
			post_processor->ImportAdditionalSensitivityBars(iTotal, iLocal, coordinate, t, it, names_t, index_sensitivity_, ui.comboBox_SensitivityItem_Additional->currentIndex()+1);
			label = ui.comboBox_SensitivityItem_Additional->currentText();
		}
		else
		{
			QString name_specie = ui.comboBox_SensitivityItem_Species->currentText();
			post_processor->ImportMassFractionSensitivityBars(iTotal, iLocal, coordinate, name_specie.toStdString(), t, it, names_t, index_sensitivity_);

			label = ui.comboBox_SensitivityItem_Species->currentText();
		}


		QBarsWidget *bars;
		bars = new QBarsWidget();
		bars->setTitle("Sensitivity analysis");
		bars->setMainComment("Sensitivity analysis - " + label);
		bars->setRectangles(t, it, names_t);
		bars->show();

		update_sensitivity_list();
	}
} 

void QW_Panel_Sensitivity::display_sensitivitychart()
{
	BzzMatrix yMatrix;
	BzzVectorInt it;
	vector<string> names_t;
	int index = ui.comboBox_SensitivityItem_Species->currentIndex()+1;

	// Post Processing
	bool iLocalAnalysis = ui.checkBox_LocalAnalysis->isChecked();
	bool iGlobalNormalization = ui.radioButton_GlobalNormalization->isChecked();

	double coordinate = ui.lineEdit_LocalAnalysis->selectedText().toDouble();

//	if (coordinate < min_coordinate || coordinate > max_coordinate)
//		MsgBox();

	QString label;
	
	if (ui.radioButton_ComboBox_Additional->isChecked() == true)
	{
		post_processor->ImportAdditionalSensitivityProfiles(iGlobalNormalization, iLocalAnalysis, coordinate, yMatrix, it, names_t, index_sensitivity_, ui.comboBox_SensitivityItem_Additional->currentIndex()+1);
		label = ui.comboBox_SensitivityItem_Additional->currentText();
	}
	else
	{
		QString name_specie = ui.comboBox_SensitivityItem_Species->currentText();
		post_processor->ImportMassFractionSensitivityProfiles(iGlobalNormalization, iLocalAnalysis, coordinate, name_specie.toStdString(), yMatrix, it, names_t, index_sensitivity_);

		label = ui.comboBox_SensitivityItem_Species->currentText();
	}

	BzzMatrix x,y;
	vector<QString> names(5);
	ChangeDimensions(5, yMatrix.Rows(), &x);
	ChangeDimensions(5, yMatrix.Rows(), &y);
	BzzVector xgrid;
	xgrid = post_processor->get_x();
	for(int j=1;j<=5;j++)
	{
		names[j-1].fromStdString(names_t[j-1]);
		for(int k=1;k<=yMatrix.Rows();k++)
		{
			x[j][k] = xgrid[k];
			y[j][k] = yMatrix[k][j];
		}
	}		

	QwtPlot2DWidget *QwtPlot2D;
	QwtPlot2D = new QwtPlot2DWidget(); 
	QwtPlot2D->setup(QString::fromStdString("x"), QString::fromStdString("y"));
	QwtPlot2D->setup(names, x, y);
	QwtPlot2D->show();

	update_sensitivity_list();

} 

void QW_Panel_Sensitivity::update_sensitivity_list()
{
	stringstream textStream;
	ui.textEdit_SensitivityList->clear();
	ui.textEdit_SensitivityList->insertPlainText(QString::fromStdString(post_processor->sensitivity_list(index_sensitivity_)));
}

void QW_Panel_Sensitivity::open_single_profileschart()
{

	int index;
	bool additional_index;
	if (ui.radioButton_ComboBox_Additional->isChecked() == true)
	{
		additional_index = true;
		index = ui.comboBox_SensitivityItem_Additional->currentIndex()+1;
	}
	else
	{
		additional_index = false;
		index = ui.comboBox_SensitivityItem_Species->currentIndex()+1;
	}


	vector<string> x_axis;
	vector<string> y_axis;

	post_processor->SetFocusSensitivityProfiles(additional_index, index);
	post_processor->ExportAvailableXAxis(x_axis);
	post_processor->ExportAvailableYAxisSensitivityCoefficients(y_axis);


	QW_Panel_Chart2D *QwtChart2D;
	QwtChart2D = new QW_Panel_Chart2D();

	QwtChart2D->setParent(post_processor);
	QwtChart2D->setKind(SENSITIVITYPROFILES);

	QwtChart2D->setup(x_axis, y_axis);
	QwtChart2D->show();
}