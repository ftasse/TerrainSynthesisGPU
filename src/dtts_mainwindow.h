#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
    class MainWindow;
    class Dialog;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);
    virtual ~MainWindow();

public slots:
    void load_exemplar();
    void load_target();
    void save_output();
    void setup_output();
    void save_exemplar();
    void run_synthesis();
    void start();

protected:
    void changeEvent(QEvent *e);


private:
    Ui::MainWindow *ui;
    Ui::Dialog *dialog_ui;
    QDialog *dialog;

    char ridges;
    std::string exemplar;
    std::string sketchf;
    std::string output;
    int bsize;
    int nw, nh;
};

#endif // MAINWINDOW_H
