#ifndef QW_PANEL_CHART2D_H
#define QW_PANEL_CHART2D_H

#include <vector>
#include <QWidget>
#include "GeneratedFiles\ui_qw_panel_chart2d.h"

class OpenSMOKE_PostProcessor;

enum kind_of_chart2d {NONE, SPECIES, REACTIONRATES, FORMATIONRATES, SENSITIVITYPROFILES};

class QW_Panel_Chart2D : public QWidget
{
	Q_OBJECT

public:
	QW_Panel_Chart2D(QWidget *parent = 0);
	~QW_Panel_Chart2D();

	kind_of_chart2d kind;

private:
	Ui::QW_Panel_Chart2DClass ui;

public:

	void setParent(OpenSMOKE_PostProcessor *post_processor);
	void setKind(const kind_of_chart2d _kind);
	void setup(std::vector<std::string> list_x_axis, std::vector<std::string> list_y_axis);

	OpenSMOKE_PostProcessor *post_processor;

public slots:

	void plot();	
	void enable_combo_01();
	void enable_combo_02();
	void enable_combo_03();
	void enable_combo_04();
	void enable_combo_05();
};

#endif // QW_PANEL_CHART2D_H
