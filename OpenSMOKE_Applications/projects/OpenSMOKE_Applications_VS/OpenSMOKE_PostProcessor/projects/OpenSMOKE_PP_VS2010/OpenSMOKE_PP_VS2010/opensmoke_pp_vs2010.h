#ifndef OPENSMOKE_PP_VS2010_H
#define OPENSMOKE_PP_VS2010_H

#include <QtGui/QMainWindow>
#include "ui_opensmoke_pp_vs2010.h"

class OpenSMOKE_PP_VS2010 : public QMainWindow
{
	Q_OBJECT

public:
	OpenSMOKE_PP_VS2010(QWidget *parent = 0, Qt::WFlags flags = 0);
	~OpenSMOKE_PP_VS2010();

private:
	Ui::OpenSMOKE_PP_VS2010Class ui;
};

#endif // OPENSMOKE_PP_VS2010_H
