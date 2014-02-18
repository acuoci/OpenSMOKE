#ifndef QBARWIDGET_H
#define QBARWIDGET_H

#include <Qt\QWidget.h>
#include "GeneratedFiles\ui_qbarwidget.h"

class QBarWidget : public QWidget
{
	Q_OBJECT

public:
	QBarWidget(QWidget *parent = 0);
	~QBarWidget();

private:
	Ui::QBarWidgetClass ui;
};

#endif // QBARWIDGET_H
