#include "BzzMath.hpp"
#include "Qwt_Plot2D.h"
#include "qwtplot2dwidget.h"
#include "qw_panel_chart2d.h"
#include "addons/OpenSMOKE_PostProcessor.h"

QW_Panel_Chart2D::QW_Panel_Chart2D(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);
	kind = NONE;
}

QW_Panel_Chart2D::~QW_Panel_Chart2D()
{

}

void QW_Panel_Chart2D::setParent(OpenSMOKE_PostProcessor *post_processor_)
{
	post_processor = post_processor_;
}

void QW_Panel_Chart2D::setKind(const kind_of_chart2d kind_)
{
	kind = kind_;
}

void QW_Panel_Chart2D::setup(std::vector<std::string> list_x_axis, 
							 std::vector<std::string> list_y_axis)
{
	for(int j=0;j<list_x_axis.size();j++)
		ui.comboBox_xAxis->addItem(list_x_axis[j].c_str());

	for(int j=0;j<list_y_axis.size();j++)
	{
		ui.comboBox_yAxis_01->addItem(list_y_axis[j].c_str());
		ui.comboBox_yAxis_02->addItem(list_y_axis[j].c_str());
		ui.comboBox_yAxis_03->addItem(list_y_axis[j].c_str());
		ui.comboBox_yAxis_04->addItem(list_y_axis[j].c_str());
		ui.comboBox_yAxis_05->addItem(list_y_axis[j].c_str());
	}
}

void QW_Panel_Chart2D::enable_combo_01()
{
	ui.comboBox_yAxis_01->setEnabled(!ui.comboBox_yAxis_01->isEnabled());
}

void QW_Panel_Chart2D::enable_combo_02()
{
	ui.comboBox_yAxis_02->setEnabled(!ui.comboBox_yAxis_02->isEnabled());
}

void QW_Panel_Chart2D::enable_combo_03()
{
	ui.comboBox_yAxis_03->setEnabled(!ui.comboBox_yAxis_03->isEnabled());
}

void QW_Panel_Chart2D::enable_combo_04()
{
	ui.comboBox_yAxis_04->setEnabled(!ui.comboBox_yAxis_04->isEnabled());
}

void QW_Panel_Chart2D::enable_combo_05()
{
	ui.comboBox_yAxis_05->setEnabled(!ui.comboBox_yAxis_05->isEnabled());
}

void QW_Panel_Chart2D::plot()
{
	// Recover x axis
	int x_curves = ui.comboBox_xAxis->currentIndex();

	// Recover y axis
	std::vector<int> y_curves;
	if (ui.checkBox_01->isChecked())	y_curves.push_back(ui.comboBox_yAxis_01->currentIndex());
	if (ui.checkBox_02->isChecked())	y_curves.push_back(ui.comboBox_yAxis_02->currentIndex());
	if (ui.checkBox_03->isChecked())	y_curves.push_back(ui.comboBox_yAxis_03->currentIndex());
	if (ui.checkBox_04->isChecked())	y_curves.push_back(ui.comboBox_yAxis_04->currentIndex());
	if (ui.checkBox_05->isChecked())	y_curves.push_back(ui.comboBox_yAxis_05->currentIndex());

	BzzMatrix x;
	BzzMatrix y;
	string name_x;
	string name_y;
	vector<string> names_lines;

	if (kind == SPECIES)
		post_processor->ImportSelectedAxis(x_curves, y_curves, x, y, name_x, name_y, names_lines);
	else if (kind == REACTIONRATES)
		post_processor->ImportSelectedAxisReactionRates(x_curves, y_curves, x, y, name_x, name_y, names_lines);
	else if (kind == FORMATIONRATES)
		post_processor->ImportSelectedAxisFormationRates(x_curves, y_curves, x, y, name_x, name_y, names_lines);
	else if (kind == SENSITIVITYPROFILES)
		post_processor->ImportSelectedAxisSensitivityCoefficients(x_curves, y_curves, x, y, name_x, name_y, names_lines);
	
	vector<QString> names;
	names.resize(y_curves.size()); 
	for(int k=1;k<=y_curves.size();k++)
		names[k-1] = QString::fromStdString(names_lines[k-1]);

	QwtPlot2DWidget *QwtPlot2D;
	QwtPlot2D = new QwtPlot2DWidget(); 
	QwtPlot2D->setup(QString::fromStdString(name_x), QString::fromStdString(name_y));
	QwtPlot2D->setup(names, x, y);
	QwtPlot2D->show();
} 