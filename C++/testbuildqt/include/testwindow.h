#include <ui_mainwindow.h>
#include <QtWidgets/QMainWindow>
#ifndef __MAINWINDOW_H__
#define __MAINWINDOW_H__


class TestWindow : public QMainWindow
{
	Q_OBJECT

public:
	explicit TestWindow(QMainWindow *parent = 0);

private:
	 Ui_MainWindow ui;
};

#endif