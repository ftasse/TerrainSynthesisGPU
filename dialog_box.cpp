#include <QtGui/QtGui>
	#include "dialog_box.h"
	
	#include "dtts_patchsynthesis.h"
	#include <iostream>
	 
	// if we include <QtGui> there is no need to include every class used: <QString>, <QFileDialog>,...
	 
	myQtApp::myQtApp(QWidget *parent)
	{
	    setupUi((QDialog*) this); // this sets up GUI
	 
	    // signals/slots mechanism in action
	    connect( Browse_target, SIGNAL( clicked() ), this, SLOT( getPath_target() ) ); 
	    connect( Browse_example, SIGNAL( clicked() ), this, SLOT( getPath_example() ) ); 
	    connect( Browse_output, SIGNAL( clicked() ), this, SLOT( getPath_output() ) ); 
	    connect( pushButton_do, SIGNAL( clicked() ), this, SLOT( doSomething() ) );
	}
	
	void myQtApp::getPath_target()
	{
	    QString path;
	    
	    path = QFileDialog::getOpenFileName(
		this,
		"Choose a file to open",
		QString::null,
		QString::null);
	 
	    lineEdit_target->setText( path );
	}
	
	void myQtApp::getPath_example()
	{
	    QString path;
	    
	    path = QFileDialog::getOpenFileName(
		this,
		"Choose a file to open",
		QString::null,
		QString::null);
	 
	    lineEdit_example->setText( path );
	} 
	 
	void myQtApp::getPath_output()
	{
	    QString path;
	    
	    path = QFileDialog::getOpenFileName(
		this,
		"Choose a file to open",
		QString::null,
		QString::null);
	 
	    lineEdit_output->setText( path );
	}
	 
	 
	void myQtApp::doSomething()
	{
	    clock_t tstart=0, tstop=0;
            tstart = clock();
            
	    Qt::CheckState state1, state2;
	    QString exemplar, sketchf, output;
	    int bsize;
	    char r = 'l';
	 
	    exemplar = lineEdit_example->text();
	    sketchf = lineEdit_target->text();
	    output = lineEdit_output->text();
	    
	    state1 = checkBox->checkState();
	    state2 = checkBox_2->checkState();
	 
	    bsize = spinBox->value();
	    
	    if ( (state1 == Qt::Checked) && (state2 == Qt::Checked)  ) r = 'a';
	    else if ( (state1 == Qt::Checked)) r = 'r';
	    else if ( (state2 == Qt::Checked)) r = 'v';
	    else r = 'l';
	    
	    terrain_synthesis((char*) exemplar.toStdString().c_str(),(char*) sketchf.toStdString().c_str(),(char*) output.toStdString().c_str(),r,bsize);
	    
	    // the stop time counter
	    tstop = clock();
	    // time in (ms)
	    std::cerr<< "Total elapsed CPU time: "<<(double)(tstop-tstart)/(double)(CLOCKS_PER_SEC)<< " (s)" << std::endl;
	}
