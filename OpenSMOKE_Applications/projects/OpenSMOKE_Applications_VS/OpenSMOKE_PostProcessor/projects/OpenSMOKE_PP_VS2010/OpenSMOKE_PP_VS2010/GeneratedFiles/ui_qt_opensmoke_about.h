/********************************************************************************
** Form generated from reading UI file 'qt_opensmoke_about.ui'
**
** Created: Wed 11. Jul 16:06:52 2012
**      by: Qt User Interface Compiler version 4.8.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QT_OPENSMOKE_ABOUT_H
#define UI_QT_OPENSMOKE_ABOUT_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Qt_OpenSMOKE_AboutClass
{
public:
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QLabel *label_6;
    QLabel *label_7;
    QLabel *label_8;

    void setupUi(QWidget *Qt_OpenSMOKE_AboutClass)
    {
        if (Qt_OpenSMOKE_AboutClass->objectName().isEmpty())
            Qt_OpenSMOKE_AboutClass->setObjectName(QString::fromUtf8("Qt_OpenSMOKE_AboutClass"));
        Qt_OpenSMOKE_AboutClass->resize(399, 244);
        label = new QLabel(Qt_OpenSMOKE_AboutClass);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(20, 20, 361, 20));
        QFont font;
        font.setPointSize(12);
        label->setFont(font);
        label->setAlignment(Qt::AlignCenter);
        label_2 = new QLabel(Qt_OpenSMOKE_AboutClass);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(70, 140, 261, 20));
        QFont font1;
        font1.setPointSize(10);
        label_2->setFont(font1);
        label_2->setAlignment(Qt::AlignCenter);
        label_3 = new QLabel(Qt_OpenSMOKE_AboutClass);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(70, 160, 261, 20));
        label_3->setFont(font1);
        label_3->setAlignment(Qt::AlignCenter);
        label_4 = new QLabel(Qt_OpenSMOKE_AboutClass);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(70, 180, 261, 20));
        label_4->setFont(font1);
        label_4->setAlignment(Qt::AlignCenter);
        label_5 = new QLabel(Qt_OpenSMOKE_AboutClass);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(70, 200, 261, 20));
        label_5->setFont(font1);
        label_5->setAlignment(Qt::AlignCenter);
        label_6 = new QLabel(Qt_OpenSMOKE_AboutClass);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(10, 50, 371, 20));
        label_6->setFont(font1);
        label_6->setAlignment(Qt::AlignCenter);
        label_7 = new QLabel(Qt_OpenSMOKE_AboutClass);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(10, 90, 371, 20));
        label_7->setFont(font1);
        label_7->setAlignment(Qt::AlignCenter);
        label_8 = new QLabel(Qt_OpenSMOKE_AboutClass);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(10, 110, 371, 20));
        label_8->setFont(font1);
        label_8->setAlignment(Qt::AlignCenter);

        retranslateUi(Qt_OpenSMOKE_AboutClass);

        QMetaObject::connectSlotsByName(Qt_OpenSMOKE_AboutClass);
    } // setupUi

    void retranslateUi(QWidget *Qt_OpenSMOKE_AboutClass)
    {
        Qt_OpenSMOKE_AboutClass->setWindowTitle(QApplication::translate("Qt_OpenSMOKE_AboutClass", "OpenSMOKE Post Processor", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("Qt_OpenSMOKE_AboutClass", "OpenSMOKE Post Processor", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("Qt_OpenSMOKE_AboutClass", "alberto.cuoci@polimi.it", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("Qt_OpenSMOKE_AboutClass", "alessio.frassoldati@polimi.it", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Qt_OpenSMOKE_AboutClass", "tiziano.faravelli@polimi.it", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("Qt_OpenSMOKE_AboutClass", "eliseo.ranzi@polimi.it", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("Qt_OpenSMOKE_AboutClass", "Version 0.1 - May 2010", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("Qt_OpenSMOKE_AboutClass", "Department of Chemistry, Materials and Chemical Engineering", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("Qt_OpenSMOKE_AboutClass", "Politecnico di Milano", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Qt_OpenSMOKE_AboutClass: public Ui_Qt_OpenSMOKE_AboutClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QT_OPENSMOKE_ABOUT_H
