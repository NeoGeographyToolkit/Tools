// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__
#include <qapplication.h>
#include <Q3PopupMenu>
#include <qlabel.h>
#include <QtGui>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qstatusbar.h>
#include <qlayout.h>
#include <QMenu>
#include <QEvent>
#include <QTextEdit>
#include <QUrl>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "appWindow.h"
#include "chooseFilesDlg.h"
#include "aspView.h"
#include "utils.h"

using namespace std;
using namespace utils;



appWindow::appWindow(QWidget* parent, std::string progName,
                     const cmdLineOptions & options):
  QMainWindow(parent, progName.c_str()),
  m_widRatio(0.3),
  m_chooseFiles(NULL), m_poly(NULL){

  installEventFilter(this);

  m_progName = progName;
  resize(options.windowWidX, options.windowWidY);

  createMenusAndMainWidget(options);
  
  return;
}

void appWindow::resizeEvent(QResizeEvent *){
  if (m_chooseFiles)
    m_chooseFiles->setMaximumSize(int(m_widRatio*size().width()), size().height());
}

void appWindow::closeEvent(QCloseEvent *){
  forceQuit();
}

void appWindow::forceQuit(){
  exit(0); // A fix for an older buggy version of Qt
}
                     
bool appWindow::eventFilter(QObject *obj, QEvent *event){

  // If the alt or control key is hit, and the focus is on the command
  // line widget, move the focus to the m_poly widget.
  if ( event->type() == QEvent::KeyPress ||
       event->type() == QEvent::ShortcutOverride){
    QKeyEvent * keyEvent = (QKeyEvent*)event;
    Qt::KeyboardModifiers modifiers = keyEvent->modifiers();
    if ( ( modifiers & Qt::ControlModifier ) ||
         ( modifiers & Qt::AltModifier ) ){
      if (m_poly != NULL) m_poly->setFocus();
    }
  }

  if (obj == m_poly) {
    // Avoid repainting on these events
    if (event->type() == QEvent::FocusIn          ||
        event->type() == QEvent::FocusOut         ||
        event->type() == QEvent::WindowDeactivate ||
        event->type() == QEvent::Leave
        ){
      return true;
    }
  }
  
  return QWidget::eventFilter(obj, event);
}

appWindow::~appWindow(){
}

void appWindow::shiftUp (){

  if (m_poly->hasFocus()){
    m_poly->shiftUp ();
  }
}

void appWindow::shiftDown (){

  if (m_poly->hasFocus()){
    m_poly->shiftDown ();
  }
}

void appWindow::createMenusAndMainWidget(const cmdLineOptions & opt){

  // There is some twisted logic here. First initialize the menus,
  // then create the main widget, then finish creating the menus.
  // This is a workaround for a bug in a certain versions of Qt. If
  // the main widget is created before the menus then it gets
  // incorrect geometry.
  
  // Create the menu bar
  QMenuBar* menu = menuBar();
  Q3PopupMenu* file = new Q3PopupMenu( menu );
  menu->insertItem("File", file);

  // Create the main panel
  if (opt.polyOptionsVec.size() > 2){

    QWidget * mainWidget;
    mainWidget = new QWidget(this);
    setCentralWidget(mainWidget);

    QSplitter *splitter = new QSplitter(mainWidget);
    m_chooseFiles = new chooseFilesDlg(this, opt.polyOptionsVec);
    m_chooseFiles->setMaximumSize(int(m_widRatio*size().width()), size().height());
    m_poly = new aspView (mainWidget, m_chooseFiles, opt);
    splitter->addWidget(m_chooseFiles);
    splitter->addWidget(m_poly);

    QGridLayout *layout = new QGridLayout(mainWidget);
    layout->addWidget (splitter, 0, 0, 0);
    mainWidget->setLayout (layout);

  }else{
    m_poly = new aspView (this, NULL, opt);
    setCentralWidget(m_poly);
  }
  
  m_poly->setFocusPolicy(Qt::StrongFocus);
  m_poly->setFocus();

  // Finish populating the menu bar
  file->insertItem("Exit", this, SLOT(forceQuit()), Qt::Key_Q);

  Q3PopupMenu* view = new Q3PopupMenu( menu );
  menu->insertItem("View", view);
  view->insertItem("Zoom out",             m_poly, SLOT(zoomOut()),      Qt::Key_Minus);
  view->insertItem("Zoom in",              m_poly, SLOT(zoomIn()),       Qt::Key_Equal);
  view->insertItem("Move left",            m_poly, SLOT(shiftLeft()),    Qt::Key_Left);
  view->insertItem("Move right",           m_poly, SLOT(shiftRight()),   Qt::Key_Right);
  view->insertItem("Move up",              this,   SLOT(shiftUp()),      Qt::Key_Up);
  view->insertItem("Move down",            this,   SLOT(shiftDown()),    Qt::Key_Down);
  view->insertItem("Reset view",           m_poly, SLOT(resetView()),    Qt::Key_R);
  view->insertItem("Change display order", m_poly, SLOT(changeOrder()),  Qt::Key_O);
  Q3PopupMenu* help = new Q3PopupMenu( menu );
  menu->insertItem("Help", help);
  //help->insertItem("Show documentation", this, SLOT(showDoc()));
  help->insertItem("About", this, SLOT(about()));

  return;
}

void appWindow::showDoc(){
  m_docWindow.setText(utils::getDocText());
  m_docWindow.setCaption(this->caption()); // Borrow the caption from the parent
  m_docWindow.show();
  return;
}

void appWindow::about(){

  string aboutStr = string("About ") + m_progName;
  std::ostringstream about_text;
  about_text << "<h3>Vision Workbench Image Viewer (vwv)</h3>"
             << "<p>Version " << VW_PACKAGE_VERSION << "</p>"
             << "<p>Copyright &copy; 2014 NASA Ames Research Center</p>";
  static QMessageBox* about
    = new QMessageBox( aboutStr.c_str(),
                       about_text.str().c_str(),
                       QMessageBox::NoIcon, 1, 0, 0, this, 0, FALSE );
  about->setButtonText( 1, "OK" );
  about->show();
  
  return;
}


//========================================================================================
// Functions for docWindow class

docWindow::docWindow(QWidget *){
  resize(900, 800);
  // Create text area then populate it
  m_textArea = new QTextEdit (this); 
  setCentralWidget (m_textArea); 
}

docWindow::~docWindow(){
  delete m_textArea;
  m_textArea = NULL;
}

void docWindow::setText(const std::string & docText){

  m_textArea->clear();
  m_textArea->setReadOnly(true);
  m_textArea->setCurrentFont(QFont("Monospace", 10)); 
  m_textArea->insertHtml(docText.c_str());
  m_textArea->moveCursor(QTextCursor::Start);
  
  return;
}
