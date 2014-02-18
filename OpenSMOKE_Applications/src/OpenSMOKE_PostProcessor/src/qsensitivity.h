#ifndef QSENSITIVITY_H
#define QSENSITIVITY_H

#include <QWidget>
#include "BzzMath.hpp"
//#include "addons/OpenSMOKE_SensitivityAnalysis_Flame1D_PostProcessor.h"
#include "addons/OpenSMOKE_PostProcessor.h"
#include "GeneratedFiles/ui_qsensitivity.h"

class QFileDialog;

class qsensitivity : public QWidget
{
	Q_OBJECT

public:
	qsensitivity(QWidget *parent = 0);
	~qsensitivity();

private:
	Ui::qsensitivityClass ui;

	QFileDialog *dialogOpen;
	QString openFileName;

	OpenSMOKE_PostProcessor								post_processor;
	
private slots:

	void openFile();
	void openChart();
	void openROPA();
	void openSensitivity();
	void openSensitivityDiffusivity();
	void openAbout();
};

#endif // QSENSITIVITY_H
