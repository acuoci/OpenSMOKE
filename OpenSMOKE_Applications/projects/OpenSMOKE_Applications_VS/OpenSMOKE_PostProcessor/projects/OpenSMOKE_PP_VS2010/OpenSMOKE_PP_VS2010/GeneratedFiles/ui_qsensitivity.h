/********************************************************************************
** Form generated from reading UI file 'qsensitivity.ui'
**
** Created: Wed 11. Jul 16:06:52 2012
**      by: Qt User Interface Compiler version 4.8.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QSENSITIVITY_H
#define UI_QSENSITIVITY_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QPushButton>
#include <QtGui/QTextEdit>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_qsensitivityClass
{
public:
    QPushButton *pushButton_open;
    QTextEdit *textEdit_status;
    QPushButton *pushButton_chart2d;
    QPushButton *pushButton_ropa;
    QPushButton *pushButton_sensitivity;
    QPushButton *pushButton_about;
    QPushButton *pushButton_sensitivity_diffusivity;

    void setupUi(QWidget *qsensitivityClass)
    {
        if (qsensitivityClass->objectName().isEmpty())
            qsensitivityClass->setObjectName(QString::fromUtf8("qsensitivityClass"));
        qsensitivityClass->resize(450, 336);
        pushButton_open = new QPushButton(qsensitivityClass);
        pushButton_open->setObjectName(QString::fromUtf8("pushButton_open"));
        pushButton_open->setGeometry(QRect(20, 20, 101, 41));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/Resources/images/icon_open.png"), QSize(), QIcon::Normal, QIcon::Off);
        pushButton_open->setIcon(icon);
        textEdit_status = new QTextEdit(qsensitivityClass);
        textEdit_status->setObjectName(QString::fromUtf8("textEdit_status"));
        textEdit_status->setGeometry(QRect(140, 20, 281, 291));
        pushButton_chart2d = new QPushButton(qsensitivityClass);
        pushButton_chart2d->setObjectName(QString::fromUtf8("pushButton_chart2d"));
        pushButton_chart2d->setEnabled(false);
        pushButton_chart2d->setGeometry(QRect(20, 80, 101, 41));
        pushButton_ropa = new QPushButton(qsensitivityClass);
        pushButton_ropa->setObjectName(QString::fromUtf8("pushButton_ropa"));
        pushButton_ropa->setEnabled(false);
        pushButton_ropa->setGeometry(QRect(20, 120, 101, 41));
        pushButton_sensitivity = new QPushButton(qsensitivityClass);
        pushButton_sensitivity->setObjectName(QString::fromUtf8("pushButton_sensitivity"));
        pushButton_sensitivity->setEnabled(false);
        pushButton_sensitivity->setGeometry(QRect(20, 160, 101, 41));
        pushButton_about = new QPushButton(qsensitivityClass);
        pushButton_about->setObjectName(QString::fromUtf8("pushButton_about"));
        pushButton_about->setEnabled(true);
        pushButton_about->setGeometry(QRect(20, 260, 101, 41));
        pushButton_sensitivity_diffusivity = new QPushButton(qsensitivityClass);
        pushButton_sensitivity_diffusivity->setObjectName(QString::fromUtf8("pushButton_sensitivity_diffusivity"));
        pushButton_sensitivity_diffusivity->setEnabled(false);
        pushButton_sensitivity_diffusivity->setGeometry(QRect(20, 200, 101, 41));

        retranslateUi(qsensitivityClass);
        QObject::connect(pushButton_open, SIGNAL(clicked()), qsensitivityClass, SLOT(openFile()));
        QObject::connect(pushButton_chart2d, SIGNAL(clicked()), qsensitivityClass, SLOT(openChart()));
        QObject::connect(pushButton_ropa, SIGNAL(clicked()), qsensitivityClass, SLOT(openROPA()));
        QObject::connect(pushButton_sensitivity, SIGNAL(clicked()), qsensitivityClass, SLOT(openSensitivity()));
        QObject::connect(pushButton_about, SIGNAL(clicked()), qsensitivityClass, SLOT(openAbout()));
        QObject::connect(pushButton_sensitivity_diffusivity, SIGNAL(clicked()), qsensitivityClass, SLOT(openSensitivityDiffusivity()));

        QMetaObject::connectSlotsByName(qsensitivityClass);
    } // setupUi

    void retranslateUi(QWidget *qsensitivityClass)
    {
        qsensitivityClass->setWindowTitle(QApplication::translate("qsensitivityClass", "OpenSMOKE - Sensitivity Analysis", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_open->setToolTip(QApplication::translate("qsensitivityClass", "Open file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        pushButton_open->setText(QApplication::translate("qsensitivityClass", "Open", 0, QApplication::UnicodeUTF8));
        pushButton_chart2d->setText(QApplication::translate("qsensitivityClass", "Chart", 0, QApplication::UnicodeUTF8));
        pushButton_ropa->setText(QApplication::translate("qsensitivityClass", "RoPA", 0, QApplication::UnicodeUTF8));
        pushButton_sensitivity->setText(QApplication::translate("qsensitivityClass", "Sensitivity\n"
"Frequency Factor", 0, QApplication::UnicodeUTF8));
        pushButton_about->setText(QApplication::translate("qsensitivityClass", "About...", 0, QApplication::UnicodeUTF8));
        pushButton_sensitivity_diffusivity->setText(QApplication::translate("qsensitivityClass", "Sensitivity\n"
"Diffusivity", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class qsensitivityClass: public Ui_qsensitivityClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QSENSITIVITY_H
