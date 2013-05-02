/********************************************************************************
** Form generated from reading UI file 'patch_basedw28657.ui'
**
** Created: Tue Mar 22 14:54:42 2011
**      by: Qt User Interface Compiler version 4.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef PATCH_BASEDW28657_H
#define PATCH_BASEDW28657_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QDialog>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpinBox>

QT_BEGIN_NAMESPACE

class Ui_Dialog
{
public:
    QLabel *label;
    QLineEdit *lineEdit_target;
    QPushButton *Browse_target;
    QPushButton *Browse_example;
    QLabel *label_2;
    QLineEdit *lineEdit_example;
    QPushButton *Browse_output;
    QLineEdit *lineEdit_output;
    QLabel *label_3;
    QSpinBox *spinBox;
    QLabel *label_4;
    QCheckBox *checkBox;
    QCheckBox *checkBox_2;
    QPushButton *pushButton_do;

    void setupUi(QDialog *Dialog)
    {
        if (Dialog->objectName().isEmpty())
            Dialog->setObjectName(QString::fromUtf8("Dialog"));
        Dialog->resize(585, 300);
        label = new QLabel(Dialog);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(10, 30, 101, 17));
        lineEdit_target = new QLineEdit(Dialog);
        lineEdit_target->setObjectName(QString::fromUtf8("lineEdit_target"));
        lineEdit_target->setGeometry(QRect(140, 30, 321, 27));
        Browse_target = new QPushButton(Dialog);
        Browse_target->setObjectName(QString::fromUtf8("Browse_target"));
        Browse_target->setGeometry(QRect(470, 30, 61, 27));
        Browse_example = new QPushButton(Dialog);
        Browse_example->setObjectName(QString::fromUtf8("Browse_example"));
        Browse_example->setGeometry(QRect(470, 80, 61, 27));
        label_2 = new QLabel(Dialog);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(10, 80, 111, 17));
        lineEdit_example = new QLineEdit(Dialog);
        lineEdit_example->setObjectName(QString::fromUtf8("lineEdit_example"));
        lineEdit_example->setGeometry(QRect(140, 80, 321, 27));
        Browse_output = new QPushButton(Dialog);
        Browse_output->setObjectName(QString::fromUtf8("Browse_output"));
        Browse_output->setGeometry(QRect(470, 130, 61, 27));
        lineEdit_output = new QLineEdit(Dialog);
        lineEdit_output->setObjectName(QString::fromUtf8("lineEdit_output"));
        lineEdit_output->setGeometry(QRect(140, 130, 321, 27));
        label_3 = new QLabel(Dialog);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(10, 130, 111, 17));
        spinBox = new QSpinBox(Dialog);
        spinBox->setObjectName(QString::fromUtf8("spinBox"));
        spinBox->setGeometry(QRect(390, 180, 59, 27));
        spinBox->setAutoFillBackground(false);
        spinBox->setMaximum(500);
        spinBox->setSingleStep(5);
        spinBox->setValue(80);
        label_4 = new QLabel(Dialog);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(300, 180, 111, 17));
        checkBox = new QCheckBox(Dialog);
        checkBox->setObjectName(QString::fromUtf8("checkBox"));
        checkBox->setGeometry(QRect(180, 170, 98, 22));
        checkBox_2 = new QCheckBox(Dialog);
        checkBox_2->setObjectName(QString::fromUtf8("checkBox_2"));
        checkBox_2->setGeometry(QRect(180, 200, 98, 22));
        pushButton_do = new QPushButton(Dialog);
        pushButton_do->setObjectName(QString::fromUtf8("pushButton_do"));
        pushButton_do->setGeometry(QRect(240, 250, 98, 27));

        retranslateUi(Dialog);

        QMetaObject::connectSlotsByName(Dialog);
    } // setupUi

    void retranslateUi(QDialog *Dialog)
    {
        Dialog->setWindowTitle(QApplication::translate("Dialog", "Dialog", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("Dialog", "Target terrain", 0, QApplication::UnicodeUTF8));
        Browse_target->setText(QApplication::translate("Dialog", "Browse", 0, QApplication::UnicodeUTF8));
        Browse_example->setText(QApplication::translate("Dialog", "Browse", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("Dialog", "Example terrain", 0, QApplication::UnicodeUTF8));
        Browse_output->setText(QApplication::translate("Dialog", "Browse", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("Dialog", "Output terrain", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("Dialog", "Patch size", 0, QApplication::UnicodeUTF8));
        checkBox->setText(QApplication::translate("Dialog", "Ridges", 0, QApplication::UnicodeUTF8));
        checkBox_2->setText(QApplication::translate("Dialog", "Valleys", 0, QApplication::UnicodeUTF8));
        pushButton_do->setText(QApplication::translate("Dialog", "START", 0, QApplication::UnicodeUTF8));
    } // retranslateUi
     
private slots:
    void getPath_target();
    void getPath_example();
    void getPath_output();
    void doSomething();

};

namespace Ui {
    class Dialog: public Ui_Dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // PATCH_BASEDH28657_H
