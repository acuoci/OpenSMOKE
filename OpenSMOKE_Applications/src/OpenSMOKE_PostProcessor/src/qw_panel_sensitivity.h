#ifndef QW_PANEL_SENSITIVITY_H
#define QW_PANEL_SENSITIVITY_H

#include <QWidget>
#include "GeneratedFiles/ui_qw_panel_sensitivity.h"

class OpenSMOKE_PostProcessor;

class QW_Panel_Sensitivity : public QWidget
{
	Q_OBJECT

public:
	QW_Panel_Sensitivity(QWidget *parent = 0);
	~QW_Panel_Sensitivity();

private:
	Ui::QW_Panel_SensitivityClass ui;

public:

	void setParent(OpenSMOKE_PostProcessor *post_processor_);
	void setup(std::vector<std::string> list_names, std::vector<std::string> list_additional, const int index_sensitivity);

	OpenSMOKE_PostProcessor *post_processor;

	double min_coordinate;
	double max_coordinate;

	void update_sensitivity_list();

	int index_sensitivity_;

public slots:

	void display_sensitivityitem();
	void display_sensitivitychart();
	void change_species_sensitivity();
	void change_additional_sensitivity();
	void open_single_profileschart();
};

#endif // QW_PANEL_SENSITIVITY_H
