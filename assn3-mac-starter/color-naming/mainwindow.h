/**********************************************************
 * File: mainwindow.h
 * Author: Keith Schwarz (htiek@cs.stanford.edu)
 *
 * Defines the main window and associated helper classes for
 * the Color Naming program.  Much of the code here is Qt-
 * specific and not standard C++.  You should not need to
 * modify this file at all in your submission.
 */

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include <QColorDialog>
#include "../KDTree.h"
#include <QThread>
#include <string>
using namespace std;

class MainWindow : public QMainWindow {
  Q_OBJECT // This is not legal C++ code.  It is interpreted by Qt's Meta Object Compiler
           // (MOC) to generate auxiliary code that makes the rest of this object work
           // correctly.
  
public:
  /* Creates and initializes the main window. */
  MainWindow(QWidget* parent = 0);
  
  /* Cleans up resources. */
  ~MainWindow();
  
private:
  /* The heart of this program is this color chooser widget, which we hook into the
   * main window to allow the user to choose colors.
   */
  QColorDialog* colorChooser;
  
  /* The kd-tree used for k-NN lookup. */
  KDTree<3, string> lookup;
  
  /* A thread class responsible for loading data in parallel with the GUI.  This
   * keeps the GUI responsive even when a huge amount of data is being loaded.
   */
  class LoadingThread;
  LoadingThread* loadingThread;
  
  friend class LoadingThread;
                            
private slots:
  /* "slots" are a Qt-specific C++ extension that represent ways for objects to
   * communicate with one another via message-passing.  "Signals" generated by
   * objects can be plugged into "slots" in another object, which causes the
   * object containing the slot to react to the signal generated by the remote
   * object.  This object has slots to react to changes in the color from the
   * color chooser, as well as notification slots to display loading progress
   * messages from the LoadingThread.
   */
  void handleColorChange(const QColor& color);
  void handleDataLoaded(int amountLoaded);
  void handleDoneIndexing();
};

/* A thread that runs in parallel with the GUI and is responsible for loading data
 * from disk.  This is done in its own thread to ensure that the GUI is responsive.
 */
class MainWindow::LoadingThread: public QThread {
  Q_OBJECT
  
public:
  explicit LoadingThread(MainWindow* master);
  virtual void run();
  
private:
  bool loadDataSet(KDTree<3, string>& dataSet);
  MainWindow* const master;
  
signals:
  /* These signals are wired into the MainWindow's slots. */
  void onDataLoaded(int amountLoaded);
  void onDoneIndexing();
};

#endif // MAINWINDOW_H