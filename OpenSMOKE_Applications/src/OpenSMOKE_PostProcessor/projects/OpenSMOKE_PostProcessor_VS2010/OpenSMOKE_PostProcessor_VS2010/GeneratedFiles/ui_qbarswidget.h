/********************************************************************************
** Form generated from reading UI file 'qbarswidget.ui'
**
** Created: Fri Sep 9 12:10:53 2011
**      by: Qt User Interface Compiler version 4.7.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QBARSWIDGET_H
#define UI_QBARSWIDGET_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_QBarsWidgetClass
{
public:

    void setupUi(QWidget *QBarsWidgetClass)
    {
        if (QBarsWidgetClass->objectName().isEmpty())
            QBarsWidgetClass->setObjectName(QString::fromUtf8("QBarsWidgetClass"));
        QBarsWidgetClass->resize(800, 500);
        QBarsWidgetClass->setMinimumSize(QSize(600, 500));
        QBarsWidgetClass->setMaximumSize(QSize(1200, 1000));

        retranslateUi(QBarsWidgetClass);

        QMetaObject::connectSlotsByName(QBarsWidgetClass);
    } // setupUi

    void retranslateUi(QWidget *QBarsWidgetClass)
    {
        QBarsWidgetClass->setWindowTitle(QApplication::translate("QBarsWidgetClass", "Sensitivity bars", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class QBarsWidgetClass: public Ui_QBarsWidgetClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QBARSWIDGET_H
