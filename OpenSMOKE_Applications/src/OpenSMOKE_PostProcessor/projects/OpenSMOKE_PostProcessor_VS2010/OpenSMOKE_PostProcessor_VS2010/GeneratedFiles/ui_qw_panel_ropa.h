/********************************************************************************
** Form generated from reading UI file 'qw_panel_ropa.ui'
**
** Created: Fri Sep 9 12:10:53 2011
**      by: Qt User Interface Compiler version 4.7.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QW_PANEL_ROPA_H
#define UI_QW_PANEL_ROPA_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QTextEdit>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_QW_Panel_ROPAClass
{
public:
    QComboBox *comboBox_IntegralSpecies;
    QLabel *label_IntegralSpecies;
    QPushButton *pushButton_UnimportantReactions_2;
    QTextEdit *textEdit_UnimportantReactions;
    QDoubleSpinBox *doubleSpinBox_IntegralSpecies;
    QLabel *label_IntegralSpecies_2;
    QLabel *label_IntegralSpecies_3;
    QRadioButton *ropaRadioButton_Global;
    QRadioButton *ropaRadioButton_Local;
    QRadioButton *ropaRadioButton_Region;
    QLabel *label;
    QLabel *label_2;
    QDoubleSpinBox *doubleSpinBox_xA;
    QDoubleSpinBox *doubleSpinBox_xB;
    QLabel *label_x_coordinate;
    QPushButton *pushButton_Update;
    QPushButton *pushButton_ReactionRates;
    QPushButton *pushButton_FormationRates;

    void setupUi(QWidget *QW_Panel_ROPAClass)
    {
        if (QW_Panel_ROPAClass->objectName().isEmpty())
            QW_Panel_ROPAClass->setObjectName(QString::fromUtf8("QW_Panel_ROPAClass"));
        QW_Panel_ROPAClass->resize(469, 501);
        comboBox_IntegralSpecies = new QComboBox(QW_Panel_ROPAClass);
        comboBox_IntegralSpecies->setObjectName(QString::fromUtf8("comboBox_IntegralSpecies"));
        comboBox_IntegralSpecies->setGeometry(QRect(70, 20, 151, 22));
        label_IntegralSpecies = new QLabel(QW_Panel_ROPAClass);
        label_IntegralSpecies->setObjectName(QString::fromUtf8("label_IntegralSpecies"));
        label_IntegralSpecies->setGeometry(QRect(20, 20, 46, 13));
        pushButton_UnimportantReactions_2 = new QPushButton(QW_Panel_ROPAClass);
        pushButton_UnimportantReactions_2->setObjectName(QString::fromUtf8("pushButton_UnimportantReactions_2"));
        pushButton_UnimportantReactions_2->setGeometry(QRect(240, 20, 81, 31));
        textEdit_UnimportantReactions = new QTextEdit(QW_Panel_ROPAClass);
        textEdit_UnimportantReactions->setObjectName(QString::fromUtf8("textEdit_UnimportantReactions"));
        textEdit_UnimportantReactions->setGeometry(QRect(10, 270, 431, 191));
        QFont font;
        font.setFamily(QString::fromUtf8("Courier"));
        textEdit_UnimportantReactions->setFont(font);
        doubleSpinBox_IntegralSpecies = new QDoubleSpinBox(QW_Panel_ROPAClass);
        doubleSpinBox_IntegralSpecies->setObjectName(QString::fromUtf8("doubleSpinBox_IntegralSpecies"));
        doubleSpinBox_IntegralSpecies->setGeometry(QRect(230, 470, 62, 22));
        doubleSpinBox_IntegralSpecies->setMinimum(0.01);
        doubleSpinBox_IntegralSpecies->setMaximum(10);
        doubleSpinBox_IntegralSpecies->setSingleStep(0.01);
        doubleSpinBox_IntegralSpecies->setValue(1);
        label_IntegralSpecies_2 = new QLabel(QW_Panel_ROPAClass);
        label_IntegralSpecies_2->setObjectName(QString::fromUtf8("label_IntegralSpecies_2"));
        label_IntegralSpecies_2->setGeometry(QRect(170, 240, 121, 21));
        label_IntegralSpecies_3 = new QLabel(QW_Panel_ROPAClass);
        label_IntegralSpecies_3->setObjectName(QString::fromUtf8("label_IntegralSpecies_3"));
        label_IntegralSpecies_3->setGeometry(QRect(140, 470, 81, 21));
        ropaRadioButton_Global = new QRadioButton(QW_Panel_ROPAClass);
        ropaRadioButton_Global->setObjectName(QString::fromUtf8("ropaRadioButton_Global"));
        ropaRadioButton_Global->setGeometry(QRect(350, 20, 82, 17));
        ropaRadioButton_Global->setChecked(true);
        ropaRadioButton_Local = new QRadioButton(QW_Panel_ROPAClass);
        ropaRadioButton_Local->setObjectName(QString::fromUtf8("ropaRadioButton_Local"));
        ropaRadioButton_Local->setGeometry(QRect(350, 50, 82, 17));
        ropaRadioButton_Region = new QRadioButton(QW_Panel_ROPAClass);
        ropaRadioButton_Region->setObjectName(QString::fromUtf8("ropaRadioButton_Region"));
        ropaRadioButton_Region->setGeometry(QRect(350, 80, 82, 17));
        label = new QLabel(QW_Panel_ROPAClass);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(70, 70, 46, 13));
        label->setAlignment(Qt::AlignCenter);
        label_2 = new QLabel(QW_Panel_ROPAClass);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(160, 70, 46, 13));
        label_2->setAlignment(Qt::AlignCenter);
        doubleSpinBox_xA = new QDoubleSpinBox(QW_Panel_ROPAClass);
        doubleSpinBox_xA->setObjectName(QString::fromUtf8("doubleSpinBox_xA"));
        doubleSpinBox_xA->setEnabled(false);
        doubleSpinBox_xA->setGeometry(QRect(50, 90, 81, 22));
        doubleSpinBox_xA->setDecimals(7);
        doubleSpinBox_xA->setMinimum(0);
        doubleSpinBox_xA->setMaximum(10);
        doubleSpinBox_xA->setSingleStep(0.01);
        doubleSpinBox_xA->setValue(0);
        doubleSpinBox_xB = new QDoubleSpinBox(QW_Panel_ROPAClass);
        doubleSpinBox_xB->setObjectName(QString::fromUtf8("doubleSpinBox_xB"));
        doubleSpinBox_xB->setEnabled(false);
        doubleSpinBox_xB->setGeometry(QRect(160, 90, 81, 22));
        doubleSpinBox_xB->setDecimals(7);
        doubleSpinBox_xB->setMinimum(0);
        doubleSpinBox_xB->setMaximum(1);
        doubleSpinBox_xB->setSingleStep(0.01);
        doubleSpinBox_xB->setValue(1);
        label_x_coordinate = new QLabel(QW_Panel_ROPAClass);
        label_x_coordinate->setObjectName(QString::fromUtf8("label_x_coordinate"));
        label_x_coordinate->setGeometry(QRect(60, 120, 171, 20));
        label_x_coordinate->setAlignment(Qt::AlignCenter);
        pushButton_Update = new QPushButton(QW_Panel_ROPAClass);
        pushButton_Update->setObjectName(QString::fromUtf8("pushButton_Update"));
        pushButton_Update->setEnabled(false);
        pushButton_Update->setGeometry(QRect(340, 110, 75, 23));
        QPalette palette;
        QBrush brush(QColor(0, 0, 0, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Text, brush);
        QBrush brush1(QColor(170, 0, 0, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::ButtonText, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Text, brush);
        palette.setBrush(QPalette::Inactive, QPalette::ButtonText, brush1);
        QBrush brush2(QColor(120, 120, 120, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Disabled, QPalette::Text, brush2);
        palette.setBrush(QPalette::Disabled, QPalette::ButtonText, brush2);
        pushButton_Update->setPalette(palette);
        QFont font1;
        font1.setBold(true);
        font1.setWeight(75);
        pushButton_Update->setFont(font1);
        pushButton_ReactionRates = new QPushButton(QW_Panel_ROPAClass);
        pushButton_ReactionRates->setObjectName(QString::fromUtf8("pushButton_ReactionRates"));
        pushButton_ReactionRates->setGeometry(QRect(60, 190, 111, 31));
        pushButton_FormationRates = new QPushButton(QW_Panel_ROPAClass);
        pushButton_FormationRates->setObjectName(QString::fromUtf8("pushButton_FormationRates"));
        pushButton_FormationRates->setGeometry(QRect(250, 190, 111, 31));

        retranslateUi(QW_Panel_ROPAClass);
        QObject::connect(pushButton_UnimportantReactions_2, SIGNAL(clicked()), QW_Panel_ROPAClass, SLOT(display_integralspecies()));
        QObject::connect(doubleSpinBox_IntegralSpecies, SIGNAL(valueChanged(double)), QW_Panel_ROPAClass, SLOT(update_unimportantreactions()));
        QObject::connect(pushButton_Update, SIGNAL(clicked()), QW_Panel_ROPAClass, SLOT(update_rateofproductionanalysis()));
        QObject::connect(ropaRadioButton_Global, SIGNAL(toggled(bool)), QW_Panel_ROPAClass, SLOT(update_updatescreen()));
        QObject::connect(ropaRadioButton_Local, SIGNAL(toggled(bool)), QW_Panel_ROPAClass, SLOT(update_updatescreen()));
        QObject::connect(ropaRadioButton_Region, SIGNAL(toggled(bool)), QW_Panel_ROPAClass, SLOT(update_updatescreen()));
        QObject::connect(doubleSpinBox_xA, SIGNAL(valueChanged(double)), QW_Panel_ROPAClass, SLOT(update_updatescreen()));
        QObject::connect(doubleSpinBox_xB, SIGNAL(valueChanged(double)), QW_Panel_ROPAClass, SLOT(update_updatescreen()));
        QObject::connect(pushButton_ReactionRates, SIGNAL(clicked()), QW_Panel_ROPAClass, SLOT(openReactionRatesChart()));
        QObject::connect(pushButton_FormationRates, SIGNAL(clicked()), QW_Panel_ROPAClass, SLOT(openFormationRatesChart()));

        QMetaObject::connectSlotsByName(QW_Panel_ROPAClass);
    } // setupUi

    void retranslateUi(QWidget *QW_Panel_ROPAClass)
    {
        QW_Panel_ROPAClass->setWindowTitle(QApplication::translate("QW_Panel_ROPAClass", "Rate of Production Analysis", 0, QApplication::UnicodeUTF8));
        label_IntegralSpecies->setText(QApplication::translate("QW_Panel_ROPAClass", "Species", 0, QApplication::UnicodeUTF8));
        pushButton_UnimportantReactions_2->setText(QApplication::translate("QW_Panel_ROPAClass", "Display", 0, QApplication::UnicodeUTF8));
        label_IntegralSpecies_2->setText(QApplication::translate("QW_Panel_ROPAClass", "Unimportant reactions", 0, QApplication::UnicodeUTF8));
        label_IntegralSpecies_3->setText(QApplication::translate("QW_Panel_ROPAClass", "Threshold (%)", 0, QApplication::UnicodeUTF8));
        ropaRadioButton_Global->setText(QApplication::translate("QW_Panel_ROPAClass", "Global", 0, QApplication::UnicodeUTF8));
        ropaRadioButton_Local->setText(QApplication::translate("QW_Panel_ROPAClass", "Local", 0, QApplication::UnicodeUTF8));
        ropaRadioButton_Region->setText(QApplication::translate("QW_Panel_ROPAClass", "Region", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("QW_Panel_ROPAClass", "xA", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("QW_Panel_ROPAClass", "xB", 0, QApplication::UnicodeUTF8));
        label_x_coordinate->setText(QApplication::translate("QW_Panel_ROPAClass", "TextLabel", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_Update->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        pushButton_Update->setText(QApplication::translate("QW_Panel_ROPAClass", "Update!", 0, QApplication::UnicodeUTF8));
        pushButton_ReactionRates->setText(QApplication::translate("QW_Panel_ROPAClass", "Reaction Rates", 0, QApplication::UnicodeUTF8));
        pushButton_FormationRates->setText(QApplication::translate("QW_Panel_ROPAClass", "Formation Rates", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class QW_Panel_ROPAClass: public Ui_QW_Panel_ROPAClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QW_PANEL_ROPA_H
