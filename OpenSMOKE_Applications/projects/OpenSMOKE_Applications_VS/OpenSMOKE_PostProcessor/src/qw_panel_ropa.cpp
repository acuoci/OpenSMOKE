#include <sstream>
#include <QMessageBox>
#include "qw_panel_chart2d.h"
#include "qw_panel_ropa.h"
#include "qbarswidget.h"
#include "addons/OpenSMOKE_PostProcessor.h"

QW_Panel_ROPA::QW_Panel_ROPA(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);
}

QW_Panel_ROPA::~QW_Panel_ROPA()
{

}

void QW_Panel_ROPA::setParent(OpenSMOKE_PostProcessor *post_processor_)
{
	post_processor = post_processor_;

	stringstream textStream;
	double eps_threshold = ui.doubleSpinBox_IntegralSpecies->value();
	post_processor->ImportUnimportantReactions(eps_threshold/100., textStream);
	ui.textEdit_UnimportantReactions->insertPlainText(QString::fromStdString(textStream.str()));

	ui.label_x_coordinate->setText(QString::fromStdString(post_processor->ExportMainXAxis()));

	xmin = post_processor->get_x().Min();
	xmax = post_processor->get_x().Max();
	
	ui.doubleSpinBox_xA->setMinimum(xmin);
	ui.doubleSpinBox_xA->setMaximum(xmax);
	ui.doubleSpinBox_xB->setMinimum(xmin);
	ui.doubleSpinBox_xB->setMaximum(xmax);
	
	ui.doubleSpinBox_xA->setValue(xmin);
	ui.doubleSpinBox_xB->setValue(xmax);

	ui.doubleSpinBox_xA->setSingleStep((xmax-xmin)/200.);
	ui.doubleSpinBox_xB->setSingleStep((xmax-xmin)/200.);
}

void QW_Panel_ROPA::setup(std::vector<std::string> list_names)
{
	for(int j=0;j<list_names.size();j++)
		ui.comboBox_IntegralSpecies->addItem(list_names[j].c_str());
}

void QW_Panel_ROPA::display_integralspecies()
{
	vector<double> p;
	vector<double> d;
	vector<int> ip;
	vector<int> id;
	vector<string> names_p;
	vector<string> names_d;
	int index = ui.comboBox_IntegralSpecies->currentIndex()+1;
//	post_processor->ImportIntegralROPA(index, p, d, ip, id, names_p, names_d);

	post_processor->ImportIntegralROPA(index, p, ip, names_p);


	QBarsWidget *bars;
	bars = new QBarsWidget();
	bars->setTitle("Rate of Production Analysis");
	bars->setMainComment("Rate of Production Analysis - " + ui.comboBox_IntegralSpecies->currentText());
	bars->setRectangles(p, ip, names_p);
	bars->show();
} 

void QW_Panel_ROPA::update_unimportantreactions()
{
	stringstream textStream;
	double eps_threshold = ui.doubleSpinBox_IntegralSpecies->value();
	post_processor->ImportUnimportantReactions(eps_threshold/100., textStream);
	ui.textEdit_UnimportantReactions->clear();
	ui.textEdit_UnimportantReactions->insertPlainText(QString::fromStdString(textStream.str()));
}

void QW_Panel_ROPA::update_rateofproductionanalysis()
{
	ui.pushButton_Update->setEnabled(false);

	if (ui.ropaRadioButton_Global->isChecked() == true)
		post_processor->CallIntegralAnalysis();
	else if (ui.ropaRadioButton_Region->isChecked() == true)
	{
		double xA = ui.doubleSpinBox_xA->value();
		double xB = ui.doubleSpinBox_xB->value();
		if (xB<=xA)
		{
			QMessageBox msgBox;
			msgBox.setText("Wrong region coordinates. They were updated to min/max values.");
			msgBox.exec();
			xA = xmin;	ui.doubleSpinBox_xA->setValue(xmin);
			xB = xmax;	ui.doubleSpinBox_xB->setValue(xmax);
		}
		post_processor->CallRegionAnalysis(xA, xB);
	}
	else if (ui.ropaRadioButton_Local->isChecked() == true)
	{
		double xA = ui.doubleSpinBox_xA->value();
		post_processor->CallLocalAnalysis(xA);
	}

	// Update unimportant reactions
	update_unimportantreactions();
}

void QW_Panel_ROPA::update_updatescreen()
{
	ui.pushButton_Update->setEnabled(true);

	if (ui.ropaRadioButton_Global->isChecked() == true)
	{
		ui.doubleSpinBox_xA->setEnabled(false);
		ui.doubleSpinBox_xB->setEnabled(false);
	}
	else if (ui.ropaRadioButton_Region->isChecked() == true)
	{
		ui.doubleSpinBox_xA->setEnabled(true);
		ui.doubleSpinBox_xB->setEnabled(true);
	}
	else if (ui.ropaRadioButton_Local->isChecked() == true)
	{
		ui.doubleSpinBox_xA->setEnabled(true);
		ui.doubleSpinBox_xB->setEnabled(false);
	}
}

void QW_Panel_ROPA::openReactionRatesChart()
{
	vector<string> x_axis;
	vector<string> y_axis;

	post_processor->ExportAvailableXAxis(x_axis);
	post_processor->ExportAvailableYAxisReactionRates(y_axis);

	QW_Panel_Chart2D *QwtChart2D;
	QwtChart2D = new QW_Panel_Chart2D();

	QwtChart2D->setParent(post_processor);
	QwtChart2D->setKind(REACTIONRATES);
	QwtChart2D->setup(x_axis, y_axis);
	QwtChart2D->show();
}

void QW_Panel_ROPA::openFormationRatesChart()
{
	vector<string> x_axis;
	vector<string> y_axis;

	post_processor->ExportAvailableXAxis(x_axis);
	post_processor->ExportAvailableYAxisFormationRates(y_axis);

	QW_Panel_Chart2D *QwtChart2D;
	QwtChart2D = new QW_Panel_Chart2D();

	QwtChart2D->setParent(post_processor);
	QwtChart2D->setKind(FORMATIONRATES);
	QwtChart2D->setup(x_axis, y_axis);
	QwtChart2D->show();
}