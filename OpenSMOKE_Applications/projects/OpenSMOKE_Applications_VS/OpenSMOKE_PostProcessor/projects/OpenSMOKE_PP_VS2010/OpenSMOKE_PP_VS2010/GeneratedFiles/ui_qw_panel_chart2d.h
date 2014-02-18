/********************************************************************************
** Form generated from reading UI file 'qw_panel_chart2d.ui'
**
** Created: Wed 11. Jul 16:06:52 2012
**      by: Qt User Interface Compiler version 4.8.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QW_PANEL_CHART2D_H
#define UI_QW_PANEL_CHART2D_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_QW_Panel_Chart2DClass
{
public:
    QComboBox *comboBox_xAxis;
    QLabel *label_xAxis;
    QPushButton *pushButton_plot;
    QComboBox *comboBox_yAxis_01;
    QLabel *label_yAxis_01;
    QComboBox *comboBox_yAxis_02;
    QComboBox *comboBox_yAxis_03;
    QComboBox *comboBox_yAxis_04;
    QComboBox *comboBox_yAxis_05;
    QLabel *label_yAxis_02;
    QLabel *label_yAxis_03;
    QLabel *label_yAxis_04;
    QLabel *label_yAxis_05;
    QCheckBox *checkBox_01;
    QCheckBox *checkBox_02;
    QCheckBox *checkBox_03;
    QCheckBox *checkBox_04;
    QCheckBox *checkBox_05;

    void setupUi(QWidget *QW_Panel_Chart2DClass)
    {
        if (QW_Panel_Chart2DClass->objectName().isEmpty())
            QW_Panel_Chart2DClass->setObjectName(QString::fromUtf8("QW_Panel_Chart2DClass"));
        QW_Panel_Chart2DClass->resize(400, 300);
        comboBox_xAxis = new QComboBox(QW_Panel_Chart2DClass);
        comboBox_xAxis->setObjectName(QString::fromUtf8("comboBox_xAxis"));
        comboBox_xAxis->setGeometry(QRect(90, 40, 191, 22));
        label_xAxis = new QLabel(QW_Panel_Chart2DClass);
        label_xAxis->setObjectName(QString::fromUtf8("label_xAxis"));
        label_xAxis->setGeometry(QRect(50, 40, 46, 13));
        pushButton_plot = new QPushButton(QW_Panel_Chart2DClass);
        pushButton_plot->setObjectName(QString::fromUtf8("pushButton_plot"));
        pushButton_plot->setGeometry(QRect(310, 40, 61, 23));
        comboBox_yAxis_01 = new QComboBox(QW_Panel_Chart2DClass);
        comboBox_yAxis_01->setObjectName(QString::fromUtf8("comboBox_yAxis_01"));
        comboBox_yAxis_01->setGeometry(QRect(90, 100, 191, 22));
        comboBox_yAxis_01->setInsertPolicy(QComboBox::InsertAlphabetically);
        label_yAxis_01 = new QLabel(QW_Panel_Chart2DClass);
        label_yAxis_01->setObjectName(QString::fromUtf8("label_yAxis_01"));
        label_yAxis_01->setGeometry(QRect(40, 100, 46, 13));
        comboBox_yAxis_02 = new QComboBox(QW_Panel_Chart2DClass);
        comboBox_yAxis_02->setObjectName(QString::fromUtf8("comboBox_yAxis_02"));
        comboBox_yAxis_02->setEnabled(false);
        comboBox_yAxis_02->setGeometry(QRect(90, 140, 191, 22));
        comboBox_yAxis_02->setInsertPolicy(QComboBox::InsertAlphabetically);
        comboBox_yAxis_03 = new QComboBox(QW_Panel_Chart2DClass);
        comboBox_yAxis_03->setObjectName(QString::fromUtf8("comboBox_yAxis_03"));
        comboBox_yAxis_03->setEnabled(false);
        comboBox_yAxis_03->setGeometry(QRect(90, 180, 191, 22));
        comboBox_yAxis_04 = new QComboBox(QW_Panel_Chart2DClass);
        comboBox_yAxis_04->setObjectName(QString::fromUtf8("comboBox_yAxis_04"));
        comboBox_yAxis_04->setEnabled(false);
        comboBox_yAxis_04->setGeometry(QRect(90, 220, 191, 22));
        comboBox_yAxis_05 = new QComboBox(QW_Panel_Chart2DClass);
        comboBox_yAxis_05->setObjectName(QString::fromUtf8("comboBox_yAxis_05"));
        comboBox_yAxis_05->setEnabled(false);
        comboBox_yAxis_05->setGeometry(QRect(90, 260, 191, 22));
        label_yAxis_02 = new QLabel(QW_Panel_Chart2DClass);
        label_yAxis_02->setObjectName(QString::fromUtf8("label_yAxis_02"));
        label_yAxis_02->setGeometry(QRect(40, 140, 46, 13));
        label_yAxis_03 = new QLabel(QW_Panel_Chart2DClass);
        label_yAxis_03->setObjectName(QString::fromUtf8("label_yAxis_03"));
        label_yAxis_03->setGeometry(QRect(40, 180, 46, 13));
        label_yAxis_04 = new QLabel(QW_Panel_Chart2DClass);
        label_yAxis_04->setObjectName(QString::fromUtf8("label_yAxis_04"));
        label_yAxis_04->setGeometry(QRect(40, 220, 46, 13));
        label_yAxis_05 = new QLabel(QW_Panel_Chart2DClass);
        label_yAxis_05->setObjectName(QString::fromUtf8("label_yAxis_05"));
        label_yAxis_05->setGeometry(QRect(40, 260, 46, 13));
        checkBox_01 = new QCheckBox(QW_Panel_Chart2DClass);
        checkBox_01->setObjectName(QString::fromUtf8("checkBox_01"));
        checkBox_01->setGeometry(QRect(300, 100, 70, 17));
        checkBox_01->setChecked(true);
        checkBox_02 = new QCheckBox(QW_Panel_Chart2DClass);
        checkBox_02->setObjectName(QString::fromUtf8("checkBox_02"));
        checkBox_02->setGeometry(QRect(300, 140, 70, 17));
        checkBox_03 = new QCheckBox(QW_Panel_Chart2DClass);
        checkBox_03->setObjectName(QString::fromUtf8("checkBox_03"));
        checkBox_03->setGeometry(QRect(300, 180, 70, 17));
        checkBox_04 = new QCheckBox(QW_Panel_Chart2DClass);
        checkBox_04->setObjectName(QString::fromUtf8("checkBox_04"));
        checkBox_04->setGeometry(QRect(300, 220, 70, 17));
        checkBox_05 = new QCheckBox(QW_Panel_Chart2DClass);
        checkBox_05->setObjectName(QString::fromUtf8("checkBox_05"));
        checkBox_05->setGeometry(QRect(300, 260, 70, 17));

        retranslateUi(QW_Panel_Chart2DClass);
        QObject::connect(checkBox_01, SIGNAL(clicked()), QW_Panel_Chart2DClass, SLOT(enable_combo_01()));
        QObject::connect(checkBox_02, SIGNAL(clicked()), QW_Panel_Chart2DClass, SLOT(enable_combo_02()));
        QObject::connect(checkBox_03, SIGNAL(clicked()), QW_Panel_Chart2DClass, SLOT(enable_combo_03()));
        QObject::connect(checkBox_04, SIGNAL(clicked()), QW_Panel_Chart2DClass, SLOT(enable_combo_04()));
        QObject::connect(checkBox_05, SIGNAL(clicked()), QW_Panel_Chart2DClass, SLOT(enable_combo_05()));
        QObject::connect(pushButton_plot, SIGNAL(clicked()), QW_Panel_Chart2DClass, SLOT(plot()));

        QMetaObject::connectSlotsByName(QW_Panel_Chart2DClass);
    } // setupUi

    void retranslateUi(QWidget *QW_Panel_Chart2DClass)
    {
        QW_Panel_Chart2DClass->setWindowTitle(QApplication::translate("QW_Panel_Chart2DClass", "Chart2D Panel", 0, QApplication::UnicodeUTF8));
        label_xAxis->setText(QApplication::translate("QW_Panel_Chart2DClass", "x Axis", 0, QApplication::UnicodeUTF8));
        pushButton_plot->setText(QApplication::translate("QW_Panel_Chart2DClass", "Plot", 0, QApplication::UnicodeUTF8));
        label_yAxis_01->setText(QApplication::translate("QW_Panel_Chart2DClass", "y Axis", 0, QApplication::UnicodeUTF8));
        label_yAxis_02->setText(QApplication::translate("QW_Panel_Chart2DClass", "y Axis", 0, QApplication::UnicodeUTF8));
        label_yAxis_03->setText(QApplication::translate("QW_Panel_Chart2DClass", "y Axis", 0, QApplication::UnicodeUTF8));
        label_yAxis_04->setText(QApplication::translate("QW_Panel_Chart2DClass", "y Axis", 0, QApplication::UnicodeUTF8));
        label_yAxis_05->setText(QApplication::translate("QW_Panel_Chart2DClass", "y Axis", 0, QApplication::UnicodeUTF8));
        checkBox_01->setText(QApplication::translate("QW_Panel_Chart2DClass", "Enabled", 0, QApplication::UnicodeUTF8));
        checkBox_02->setText(QApplication::translate("QW_Panel_Chart2DClass", "Enabled", 0, QApplication::UnicodeUTF8));
        checkBox_03->setText(QApplication::translate("QW_Panel_Chart2DClass", "Enabled", 0, QApplication::UnicodeUTF8));
        checkBox_04->setText(QApplication::translate("QW_Panel_Chart2DClass", "Enabled", 0, QApplication::UnicodeUTF8));
        checkBox_05->setText(QApplication::translate("QW_Panel_Chart2DClass", "Enabled", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class QW_Panel_Chart2DClass: public Ui_QW_Panel_Chart2DClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QW_PANEL_CHART2D_H
