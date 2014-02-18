/********************************************************************************
** Form generated from reading UI file 'qwtplot2dwidget.ui'
**
** Created: Fri Sep 9 12:10:53 2011
**      by: Qt User Interface Compiler version 4.7.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QWTPLOT2DWIDGET_H
#define UI_QWTPLOT2DWIDGET_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QPushButton>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_QwtPlot2DWidgetClass
{
public:
    QPushButton *pushButton_autoscale;
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButton_zoomin;
    QPushButton *pushButton_zoomout;
    QCheckBox *checkBox_xLog;
    QCheckBox *checkBox_yLog;
    QCheckBox *checkBox_xNormal;
    QCheckBox *checkBox_yNormal;
    QPushButton *pushButton_autoscale_2;

    void setupUi(QWidget *QwtPlot2DWidgetClass)
    {
        if (QwtPlot2DWidgetClass->objectName().isEmpty())
            QwtPlot2DWidgetClass->setObjectName(QString::fromUtf8("QwtPlot2DWidgetClass"));
        QwtPlot2DWidgetClass->resize(475, 423);
        pushButton_autoscale = new QPushButton(QwtPlot2DWidgetClass);
        pushButton_autoscale->setObjectName(QString::fromUtf8("pushButton_autoscale"));
        pushButton_autoscale->setGeometry(QRect(260, 10, 81, 23));
        horizontalLayoutWidget = new QWidget(QwtPlot2DWidgetClass);
        horizontalLayoutWidget->setObjectName(QString::fromUtf8("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(30, 110, 431, 291));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        pushButton_zoomin = new QPushButton(QwtPlot2DWidgetClass);
        pushButton_zoomin->setObjectName(QString::fromUtf8("pushButton_zoomin"));
        pushButton_zoomin->setGeometry(QRect(260, 40, 81, 23));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/Resources/images/icon_zoomin.png"), QSize(), QIcon::Normal, QIcon::Off);
        pushButton_zoomin->setIcon(icon);
        pushButton_zoomout = new QPushButton(QwtPlot2DWidgetClass);
        pushButton_zoomout->setObjectName(QString::fromUtf8("pushButton_zoomout"));
        pushButton_zoomout->setGeometry(QRect(260, 70, 81, 23));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/Resources/images/icon_zoomout.png"), QSize(), QIcon::Normal, QIcon::Off);
        pushButton_zoomout->setIcon(icon1);
        checkBox_xLog = new QCheckBox(QwtPlot2DWidgetClass);
        checkBox_xLog->setObjectName(QString::fromUtf8("checkBox_xLog"));
        checkBox_xLog->setGeometry(QRect(40, 20, 70, 17));
        checkBox_yLog = new QCheckBox(QwtPlot2DWidgetClass);
        checkBox_yLog->setObjectName(QString::fromUtf8("checkBox_yLog"));
        checkBox_yLog->setGeometry(QRect(40, 40, 70, 17));
        checkBox_xNormal = new QCheckBox(QwtPlot2DWidgetClass);
        checkBox_xNormal->setObjectName(QString::fromUtf8("checkBox_xNormal"));
        checkBox_xNormal->setGeometry(QRect(140, 20, 70, 17));
        checkBox_yNormal = new QCheckBox(QwtPlot2DWidgetClass);
        checkBox_yNormal->setObjectName(QString::fromUtf8("checkBox_yNormal"));
        checkBox_yNormal->setGeometry(QRect(140, 40, 70, 17));
        pushButton_autoscale_2 = new QPushButton(QwtPlot2DWidgetClass);
        pushButton_autoscale_2->setObjectName(QString::fromUtf8("pushButton_autoscale_2"));
        pushButton_autoscale_2->setGeometry(QRect(350, 10, 101, 23));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/Resources/images/icon_colors.png"), QSize(), QIcon::Normal, QIcon::Off);
        pushButton_autoscale_2->setIcon(icon2);

        retranslateUi(QwtPlot2DWidgetClass);
        QObject::connect(pushButton_autoscale, SIGNAL(clicked()), QwtPlot2DWidgetClass, SLOT(setAutoscale()));
        QObject::connect(pushButton_zoomin, SIGNAL(clicked()), QwtPlot2DWidgetClass, SLOT(setZoomIn()));
        QObject::connect(pushButton_zoomout, SIGNAL(clicked()), QwtPlot2DWidgetClass, SLOT(setZoomOut()));
        QObject::connect(checkBox_xLog, SIGNAL(clicked()), QwtPlot2DWidgetClass, SLOT(setLogX()));
        QObject::connect(checkBox_yLog, SIGNAL(clicked()), QwtPlot2DWidgetClass, SLOT(setLogY()));
        QObject::connect(checkBox_xNormal, SIGNAL(clicked()), QwtPlot2DWidgetClass, SLOT(setNormalX()));
        QObject::connect(checkBox_yNormal, SIGNAL(clicked()), QwtPlot2DWidgetClass, SLOT(setNormalY()));
        QObject::connect(pushButton_autoscale_2, SIGNAL(clicked()), QwtPlot2DWidgetClass, SLOT(setColors()));

        QMetaObject::connectSlotsByName(QwtPlot2DWidgetClass);
    } // setupUi

    void retranslateUi(QWidget *QwtPlot2DWidgetClass)
    {
        QwtPlot2DWidgetClass->setWindowTitle(QApplication::translate("QwtPlot2DWidgetClass", "Plot2D", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_autoscale->setToolTip(QApplication::translate("QwtPlot2DWidgetClass", "Autoscale x and y axis", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        pushButton_autoscale->setText(QApplication::translate("QwtPlot2DWidgetClass", "Autoscale", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_zoomin->setToolTip(QApplication::translate("QwtPlot2DWidgetClass", "Zoom in", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        pushButton_zoomin->setText(QApplication::translate("QwtPlot2DWidgetClass", "Zoom In", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_zoomout->setToolTip(QApplication::translate("QwtPlot2DWidgetClass", "Zoom out", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        pushButton_zoomout->setText(QApplication::translate("QwtPlot2DWidgetClass", "Zoom Out", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        checkBox_xLog->setToolTip(QApplication::translate("QwtPlot2DWidgetClass", "logarithmic scale x axis", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        checkBox_xLog->setText(QApplication::translate("QwtPlot2DWidgetClass", "x Log", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        checkBox_yLog->setToolTip(QApplication::translate("QwtPlot2DWidgetClass", "logarithmic scale y axis", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        checkBox_yLog->setText(QApplication::translate("QwtPlot2DWidgetClass", "y Log", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        checkBox_xNormal->setToolTip(QApplication::translate("QwtPlot2DWidgetClass", "normalized profiles x axis", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        checkBox_xNormal->setText(QApplication::translate("QwtPlot2DWidgetClass", "x Normal", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        checkBox_yNormal->setToolTip(QApplication::translate("QwtPlot2DWidgetClass", "normalized profiles y axis", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        checkBox_yNormal->setText(QApplication::translate("QwtPlot2DWidgetClass", "y Normal", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_autoscale_2->setToolTip(QApplication::translate("QwtPlot2DWidgetClass", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'MS Shell Dlg 2'; font-size:8.25pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">Change plot colors</span></p></body></html>", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        pushButton_autoscale_2->setText(QApplication::translate("QwtPlot2DWidgetClass", "Change colors", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class QwtPlot2DWidgetClass: public Ui_QwtPlot2DWidgetClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QWTPLOT2DWIDGET_H
