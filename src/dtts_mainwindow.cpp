#include "dtts_mainwindow.h"
#include "dtts_image.h"
#include "dtts_patchsynthesis.h"
#include "ui_mainwindow.h"
#include "ui_synthesis_dialog.h"

#include <QFileDialog>
#include <QMessageBox>
#include <QProgressDialog>
#include <QTime>

QImage imageToQImage(Dtts::Image& image);
Dtts::Terrain loadTerrain(std::string fname);
void saveTerrain(Dtts::Terrain& terrain, std::string fname);

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    dialog_ui = new Ui::Dialog ();
    dialog = new QDialog (this);
    dialog_ui->setupUi(dialog);

    ridges = 's';
    bsize = 80;
    nw = 512;
    nh = 512;
    output = "result.ter";

    dialog_ui->result_box->setText(output.c_str());
    dialog_ui->width_box->setValue(nw);
    dialog_ui->height_box->setValue(nh);
    dialog_ui->bsize_box->setValue(bsize);

    connect(ui->actionLoad_Target, SIGNAL(triggered()), this, SLOT(load_target()));
    connect(ui->actionLoad_Exemplar, SIGNAL(triggered()), this, SLOT(load_exemplar()));
    connect(ui->actionSave_Output, SIGNAL(triggered()), this, SLOT(save_output()));
    connect(ui->actionSave_Exemplar, SIGNAL(triggered()), this, SLOT(save_exemplar()));
    connect(ui->actionRun_Synthesis, SIGNAL(triggered()), this, SLOT(run_synthesis()));

    connect(dialog_ui->target_button, SIGNAL(clicked()), this, SLOT(load_target()));
    connect(dialog_ui->exem_button, SIGNAL(clicked()), this, SLOT(load_exemplar()));
    connect(dialog_ui->result_button, SIGNAL(clicked()), this, SLOT(setup_output()));

    connect(dialog, SIGNAL(accepted()), this, SLOT(start()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::load_exemplar()
{
    QString fileName = QFileDialog::getOpenFileName(this,
         tr("Open Terrain"), "", tr("Terrain Files (*.pgm *.ter)"));
    if (fileName.size()>0)
    {
        exemplar = fileName.toStdString();
        dialog_ui->exem_box->setText(exemplar.c_str());

        Dtts::Terrain terrain = loadTerrain(exemplar);
        QImage qimg = imageToQImage(terrain);
        ui->exem_widget->setScaledContents(false);
        ui->exem_widget->setPixmap(QPixmap::fromImage(qimg));
    }
}

void MainWindow::load_target()
{
    QString fileName = QFileDialog::getOpenFileName(this,
         tr("Open Terrain"), "", tr("Terrain Files (*.pgm *.ter)"));
    if (fileName.size()>0)
    {
        sketchf = fileName.toStdString();
        Dtts::Terrain terrain = loadTerrain(sketchf);
        nw = terrain.width();
        nh = terrain.height();
        bsize = std::min(bsize, std::min(nw, nh)/5);

        dialog_ui->target_box->setText(sketchf.c_str());
        dialog_ui->width_box->setValue(nw);
        dialog_ui->height_box->setValue(nh);
        dialog_ui->bsize_box->setValue(bsize);

        QImage qimg = imageToQImage(terrain);
        ui->target_widget->setScaledContents(false);
        ui->target_widget->setPixmap(QPixmap::fromImage(qimg));
    }
}

void MainWindow::save_output()
{
    QString fileName = QFileDialog::getSaveFileName(this,
         tr("Save Terrain"), "", tr("Terrain Files (*.pgm *.ter)"));
    if (fileName.size()>0)
    {
        if (ui->result_widget->pixmap () == 0)
        {
            QMessageBox msgBox;
            msgBox.setText("No output is available yet.");
            msgBox.setInformativeText("Do you want to run a terrain synthesis?");
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
            msgBox.setDefaultButton(QMessageBox::Yes);
            int ret = msgBox.exec();

            switch (ret) {
               case QMessageBox::Yes:
                   output = fileName.toStdString();
                   dialog_ui->result_box->setText(output.c_str());
                   run_synthesis();
                   break;
               case QMessageBox::No:
                   break;
               default:
                   // should never be reached
                   break;
             }
        } else
        {
            Dtts::Terrain terrain = loadTerrain(output);
            saveTerrain(terrain, fileName.toStdString());
        }
    }
}

void MainWindow::save_exemplar()
{
    if (exemplar.size()==0)
        return;

    QString fileName = QFileDialog::getSaveFileName(this,
         tr("Save Terrain"), "", tr("Terrain Files (*.pgm *.ter)"));
    if (fileName.size()>0)
    {
        Dtts::Terrain terrain = loadTerrain(exemplar);
        saveTerrain(terrain, fileName.toStdString());
    }
}

void MainWindow::setup_output()
{
    QString fileName = QFileDialog::getSaveFileName(this,
         tr("Save Terrain"), "", tr("Terrain Files (*.pgm *.ter)"));
    if (fileName.size()>0)
    {
        output = fileName.toStdString();
        dialog_ui->result_box->setText(output.c_str());
    }
}

void MainWindow::run_synthesis()
{
    if (ridges == 's')
        dialog_ui->features_box->setCurrentIndex(0);
    else if (ridges == 'r')
        dialog_ui->features_box->setCurrentIndex(1);
    else if (ridges == 'v')
        dialog_ui->features_box->setCurrentIndex(2);
    else
        dialog_ui->features_box->setCurrentIndex(3);

    dialog->show();
}

void MainWindow::start()
{
    nw = dialog_ui->width_box->value();
    nh = dialog_ui->height_box->value();
    bsize = dialog_ui->bsize_box->value();

    exemplar = dialog_ui->exem_box->text().toStdString();
    sketchf = dialog_ui->target_box->text().toStdString();
    output = dialog_ui->result_box->text().toStdString();

    if (dialog_ui->features_box->currentIndex() == 0)
        ridges == 's';
    if (dialog_ui->features_box->currentIndex() == 1)
        ridges == 'r';
    if (dialog_ui->features_box->currentIndex() == 2)
        ridges == 'v';
    if (dialog_ui->features_box->currentIndex() == 3)
        ridges == 'l';

    if (exemplar.size() == 0 || output.size() == 0 || nw == 0 || nh == 0 || bsize == 0)
        run_synthesis();
    else
    {
        QProgressDialog d;
        d.setMaximum(50);
        d.setValue(10);
        d.show();
        QApplication::processEvents();

        QTime myTimer;
        myTimer.start();

        int osize = bsize/4;
        if (ridges == 's')
            terrain_synthesis((char *)exemplar.c_str(),(char *)output.c_str(),nw,nh,bsize,osize);
        else
            terrain_synthesis((char *)exemplar.c_str(),(char *)sketchf.c_str(),(char *)output.c_str(),ridges,bsize);

        int nMilliseconds = myTimer.elapsed();
        qDebug("Elapsed time: %.4f secs\n", nMilliseconds/1000.0);
        d.setValue(50);

        Dtts::Terrain terrain = loadTerrain(output);
        QImage qimg = imageToQImage(terrain);
        ui->result_widget->setScaledContents(false);
        ui->result_widget->setPixmap(QPixmap::fromImage(qimg));
    }
}

void MainWindow::changeEvent(QEvent *e)
{
    QMainWindow::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

QImage imageToQImage(Dtts::Image &image)
{
    Dtts::Image newImg(image.width(), image.height());
    for (int i=0; i<image.width(); ++i)
    {
        for (int j=0; j<image.height(); ++j)
        {
            newImg.setPixel(i, j, image.getPixel(i, j));
        }
    }
    normalizeImage(newImg);

    QImage qimg(image.width(), image.height(), QImage::Format_RGB888);
    for (int i=0; i<image.width(); ++i)
    {
        for (int j=0; j<image.height(); ++j)
        {
            float val = newImg.getPixel(i, j);
            int ival = (int) val;
            QColor color(ival, ival, ival);
            qimg.setPixel(i, j, color.rgb());
        }
    }

    return qimg;
}

void saveTerrain(Dtts::Terrain &terrain, std::string fname)
{
    if (strstr (fname.c_str(),".ter") || strstr (fname.c_str(),".TER") )
    {
        qDebug("Save terrain: %s", fname.c_str());
        terrain.saveTerragen(fname.c_str());
    }
    else
    {
        terrain.savePGM(fname.c_str());
    }
}

Dtts::Terrain loadTerrain(string fname)
{
    Dtts::Terrain terrain;

    if (strstr (fname.c_str(),".ter") || strstr (fname.c_str(),".TER") )
    {
        terrain.loadTerragen(fname.c_str());
    }

    else
    {
        terrain.loadPGM(fname.c_str());

        terrain.mscale = 30;
        terrain.mbaseheight=50;
        terrain.mheightscale = 100;
    }
    return terrain;
}
