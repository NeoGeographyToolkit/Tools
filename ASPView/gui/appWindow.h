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
#ifndef APPWINDOW_H
#define APPWINDOW_H

#include <qmainwindow.h>
#include <qlineedit.h>
#include <QEvent>
#include <string>
#include <vector>

class aspView;
class chooseFilesDlg;
class QResizeEvent;
class QCloseEvent;
class QTextEdit;

// Do this to hide include the file.
namespace utils{
  class cmdLineOptions;
}


/// Simple window class to display HTML text
class docWindow: public QMainWindow{
  Q_OBJECT

public:
  docWindow(QWidget * parent = NULL);
  virtual ~docWindow();
  void setText(const std::string & docText);
private:
  QTextEdit * m_textArea;
};




/// Class for the main GUI window
class appWindow : public QMainWindow {
  Q_OBJECT

public:
  /// Constructor
  appWindow(QWidget* parent, std::string progName,
            const utils::cmdLineOptions & options);
  ~appWindow();
  
protected:

  /// Override default QT handling of some events.
  bool eventFilter(QObject *obj, QEvent *event);

private slots:
  void createMenusAndMainWidget(const utils::cmdLineOptions & opt);
  void showDoc();    ///< Pop up documentation window?
  void about();      ///< Pop up a little message box with program information.
  void shiftUp();    ///< Redirect call to aspView object.
  void shiftDown();  ///< Redirect call to aspView object.
  void forceQuit();  ///< Ensure the program shuts down.
  
private:
  // Event handlers
  void resizeEvent(QResizeEvent *);
  void closeEvent (QCloseEvent *);

  double           m_widRatio;    ///< This is size of the window is controlled?
  chooseFilesDlg * m_chooseFiles; ///< Panel on right side for selecting files.
  aspView        * m_poly;        ///< Central image panel.
  std::string      m_progName;    ///< Name of the program.
  docWindow        m_docWindow;   ///< Documentation window instance.
};


#endif
