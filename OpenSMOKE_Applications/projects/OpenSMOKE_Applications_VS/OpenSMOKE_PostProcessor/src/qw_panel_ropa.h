#ifndef QW_PANEL_ROPA_H
#define QW_PANEL_ROPA_H

#include <QWidget>
#include "GeneratedFiles/ui_qw_panel_ropa.h"

class OpenSMOKE_PostProcessor;

class QW_Panel_ROPA : public QWidget
{
	Q_OBJECT

public:
	QW_Panel_ROPA(QWidget *parent = 0);
	~QW_Panel_ROPA();

private:
	Ui::QW_Panel_ROPAClass ui;

	double xmin;
	double xmax;

public:
	void setParent(OpenSMOKE_PostProcessor *post_processor);
	void setup(std::vector<std::string> list_names);

	OpenSMOKE_PostProcessor *post_processor;

public slots:

	void display_integralspecies();
	void update_unimportantreactions();
	void update_rateofproductionanalysis();
	void update_updatescreen();

	void openReactionRatesChart();
	void openFormationRatesChart();
};

#endif // QW_PANEL_ROPA_H
