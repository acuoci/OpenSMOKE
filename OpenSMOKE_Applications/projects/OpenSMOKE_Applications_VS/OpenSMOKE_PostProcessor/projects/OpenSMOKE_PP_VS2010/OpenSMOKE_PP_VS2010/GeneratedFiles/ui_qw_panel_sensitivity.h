/********************************************************************************
** Form generated from reading UI file 'qw_panel_sensitivity.ui'
**
** Created: Wed 11. Jul 16:06:52 2012
**      by: Qt User Interface Compiler version 4.8.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QW_PANEL_SENSITIVITY_H
#define UI_QW_PANEL_SENSITIVITY_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QTextEdit>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_QW_Panel_SensitivityClass
{
public:
    QGroupBox *groupBox;
    QRadioButton *radioButton_GlobalNormalization;
    QRadioButton *radioButton_LocalNormalization;
    QGroupBox *groupBox_2;
    QLineEdit *lineEdit_LocalAnalysis;
    QCheckBox *checkBox_LocalAnalysis;
    QLabel *label_SensitivityItem_2;
    QGroupBox *groupBox_3;
    QComboBox *comboBox_SensitivityItem_Species;
    QRadioButton *radioButton_ComboBox_Additional;
    QComboBox *comboBox_SensitivityItem_Additional;
    QLabel *label_SensitivityItem;
    QRadioButton *radioButton_ComboBox_Species;
    QLabel *label_SensitivityItem_3;
    QPushButton *pushButton_Display;
    QPushButton *pushButton_Plot2D;
    QPushButton *pushButton_SingleProfiles;
    QTextEdit *textEdit_SensitivityList;

    void setupUi(QWidget *QW_Panel_SensitivityClass)
    {
        if (QW_Panel_SensitivityClass->objectName().isEmpty())
            QW_Panel_SensitivityClass->setObjectName(QString::fromUtf8("QW_Panel_SensitivityClass"));
        QW_Panel_SensitivityClass->resize(589, 397);
        groupBox = new QGroupBox(QW_Panel_SensitivityClass);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(60, 140, 111, 91));
        radioButton_GlobalNormalization = new QRadioButton(groupBox);
        radioButton_GlobalNormalization->setObjectName(QString::fromUtf8("radioButton_GlobalNormalization"));
        radioButton_GlobalNormalization->setGeometry(QRect(20, 30, 82, 17));
        radioButton_GlobalNormalization->setChecked(true);
        radioButton_LocalNormalization = new QRadioButton(groupBox);
        radioButton_LocalNormalization->setObjectName(QString::fromUtf8("radioButton_LocalNormalization"));
        radioButton_LocalNormalization->setGeometry(QRect(20, 60, 82, 17));
        groupBox_2 = new QGroupBox(QW_Panel_SensitivityClass);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        groupBox_2->setGeometry(QRect(210, 140, 171, 91));
        lineEdit_LocalAnalysis = new QLineEdit(groupBox_2);
        lineEdit_LocalAnalysis->setObjectName(QString::fromUtf8("lineEdit_LocalAnalysis"));
        lineEdit_LocalAnalysis->setGeometry(QRect(80, 60, 81, 20));
        lineEdit_LocalAnalysis->setInputMethodHints(Qt::ImhDigitsOnly);
        checkBox_LocalAnalysis = new QCheckBox(groupBox_2);
        checkBox_LocalAnalysis->setObjectName(QString::fromUtf8("checkBox_LocalAnalysis"));
        checkBox_LocalAnalysis->setGeometry(QRect(40, 30, 101, 17));
        label_SensitivityItem_2 = new QLabel(groupBox_2);
        label_SensitivityItem_2->setObjectName(QString::fromUtf8("label_SensitivityItem_2"));
        label_SensitivityItem_2->setGeometry(QRect(10, 60, 71, 21));
        groupBox_3 = new QGroupBox(QW_Panel_SensitivityClass);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        groupBox_3->setGeometry(QRect(20, 10, 551, 121));
        comboBox_SensitivityItem_Species = new QComboBox(groupBox_3);
        comboBox_SensitivityItem_Species->setObjectName(QString::fromUtf8("comboBox_SensitivityItem_Species"));
        comboBox_SensitivityItem_Species->setEnabled(false);
        comboBox_SensitivityItem_Species->setGeometry(QRect(180, 60, 161, 22));
        radioButton_ComboBox_Additional = new QRadioButton(groupBox_3);
        radioButton_ComboBox_Additional->setObjectName(QString::fromUtf8("radioButton_ComboBox_Additional"));
        radioButton_ComboBox_Additional->setGeometry(QRect(30, 20, 31, 17));
        radioButton_ComboBox_Additional->setChecked(true);
        comboBox_SensitivityItem_Additional = new QComboBox(groupBox_3);
        comboBox_SensitivityItem_Additional->setObjectName(QString::fromUtf8("comboBox_SensitivityItem_Additional"));
        comboBox_SensitivityItem_Additional->setEnabled(true);
        comboBox_SensitivityItem_Additional->setGeometry(QRect(180, 20, 161, 22));
        label_SensitivityItem = new QLabel(groupBox_3);
        label_SensitivityItem->setObjectName(QString::fromUtf8("label_SensitivityItem"));
        label_SensitivityItem->setGeometry(QRect(50, 60, 121, 21));
        radioButton_ComboBox_Species = new QRadioButton(groupBox_3);
        radioButton_ComboBox_Species->setObjectName(QString::fromUtf8("radioButton_ComboBox_Species"));
        radioButton_ComboBox_Species->setGeometry(QRect(30, 60, 31, 17));
        radioButton_ComboBox_Species->setChecked(false);
        label_SensitivityItem_3 = new QLabel(groupBox_3);
        label_SensitivityItem_3->setObjectName(QString::fromUtf8("label_SensitivityItem_3"));
        label_SensitivityItem_3->setGeometry(QRect(50, 20, 121, 21));
        pushButton_Display = new QPushButton(groupBox_3);
        pushButton_Display->setObjectName(QString::fromUtf8("pushButton_Display"));
        pushButton_Display->setGeometry(QRect(390, 20, 101, 31));
        pushButton_Plot2D = new QPushButton(groupBox_3);
        pushButton_Plot2D->setObjectName(QString::fromUtf8("pushButton_Plot2D"));
        pushButton_Plot2D->setGeometry(QRect(390, 50, 101, 31));
        pushButton_SingleProfiles = new QPushButton(groupBox_3);
        pushButton_SingleProfiles->setObjectName(QString::fromUtf8("pushButton_SingleProfiles"));
        pushButton_SingleProfiles->setGeometry(QRect(390, 80, 101, 31));
        textEdit_SensitivityList = new QTextEdit(QW_Panel_SensitivityClass);
        textEdit_SensitivityList->setObjectName(QString::fromUtf8("textEdit_SensitivityList"));
        textEdit_SensitivityList->setGeometry(QRect(20, 250, 551, 141));

        retranslateUi(QW_Panel_SensitivityClass);
        QObject::connect(pushButton_Display, SIGNAL(clicked()), QW_Panel_SensitivityClass, SLOT(display_sensitivityitem()));
        QObject::connect(radioButton_ComboBox_Species, SIGNAL(clicked()), QW_Panel_SensitivityClass, SLOT(change_species_sensitivity()));
        QObject::connect(radioButton_ComboBox_Additional, SIGNAL(clicked()), QW_Panel_SensitivityClass, SLOT(change_additional_sensitivity()));
        QObject::connect(pushButton_Plot2D, SIGNAL(clicked()), QW_Panel_SensitivityClass, SLOT(display_sensitivitychart()));
        QObject::connect(pushButton_SingleProfiles, SIGNAL(clicked()), QW_Panel_SensitivityClass, SLOT(open_single_profileschart()));

        QMetaObject::connectSlotsByName(QW_Panel_SensitivityClass);
    } // setupUi

    void retranslateUi(QWidget *QW_Panel_SensitivityClass)
    {
        QW_Panel_SensitivityClass->setWindowTitle(QApplication::translate("QW_Panel_SensitivityClass", "OpenSMOKE Sensitivity Analysis", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("QW_Panel_SensitivityClass", "Normalization", 0, QApplication::UnicodeUTF8));
        radioButton_GlobalNormalization->setText(QApplication::translate("QW_Panel_SensitivityClass", "Global", 0, QApplication::UnicodeUTF8));
        radioButton_LocalNormalization->setText(QApplication::translate("QW_Panel_SensitivityClass", "Local", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("QW_Panel_SensitivityClass", "Local Analysis", 0, QApplication::UnicodeUTF8));
        checkBox_LocalAnalysis->setText(QApplication::translate("QW_Panel_SensitivityClass", "Local Analysis", 0, QApplication::UnicodeUTF8));
        label_SensitivityItem_2->setText(QApplication::translate("QW_Panel_SensitivityClass", "Coordinate", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("QW_Panel_SensitivityClass", "Sensitivity Analisys", 0, QApplication::UnicodeUTF8));
        radioButton_ComboBox_Additional->setText(QString());
        label_SensitivityItem->setText(QApplication::translate("QW_Panel_SensitivityClass", "Species (mass fractions)", 0, QApplication::UnicodeUTF8));
        radioButton_ComboBox_Species->setText(QString());
        label_SensitivityItem_3->setText(QApplication::translate("QW_Panel_SensitivityClass", "Additional variables", 0, QApplication::UnicodeUTF8));
        pushButton_Display->setText(QApplication::translate("QW_Panel_SensitivityClass", "Sensitivity Bars", 0, QApplication::UnicodeUTF8));
        pushButton_Plot2D->setText(QApplication::translate("QW_Panel_SensitivityClass", "Sensitivity Profiles", 0, QApplication::UnicodeUTF8));
        pushButton_SingleProfiles->setText(QApplication::translate("QW_Panel_SensitivityClass", "Single Profiles", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class QW_Panel_SensitivityClass: public Ui_QW_Panel_SensitivityClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QW_PANEL_SENSITIVITY_H
