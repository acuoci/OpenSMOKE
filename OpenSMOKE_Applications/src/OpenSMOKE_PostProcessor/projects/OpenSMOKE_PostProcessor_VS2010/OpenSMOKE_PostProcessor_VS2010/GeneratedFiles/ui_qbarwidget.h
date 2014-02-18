/********************************************************************************
** Form generated from reading UI file 'qbarwidget.ui'
**
** Created: Fri Sep 9 12:10:53 2011
**      by: Qt User Interface Compiler version 4.7.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QBARWIDGET_H
#define UI_QBARWIDGET_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_QBarWidgetClass
{
public:

    void setupUi(QWidget *QBarWidgetClass)
    {
        if (QBarWidgetClass->objectName().isEmpty())
            QBarWidgetClass->setObjectName(QString::fromUtf8("QBarWidgetClass"));
        QBarWidgetClass->resize(400, 300);

        retranslateUi(QBarWidgetClass);

        QMetaObject::connectSlotsByName(QBarWidgetClass);
    } // setupUi

    void retranslateUi(QWidget *QBarWidgetClass)
    {
        QBarWidgetClass->setWindowTitle(QApplication::translate("QBarWidgetClass", "QBarWidget", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class QBarWidgetClass: public Ui_QBarWidgetClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QBARWIDGET_H
