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
#ifndef ASPVIEW_H
#define ASPVIEW_H

#include <QPolygon>
#include <QContextMenuEvent>
#include <QEvent>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPaintEvent>
#include <QPixmap>
#include <QWheelEvent>
#include <QWidget>
#include <vector>
#include <map>
#include "utils.h"
#include "chooseFilesDlg.h"

class aspView : public QWidget{
  Q_OBJECT
public:
  aspView(QWidget *parent, chooseFilesDlg *chooseFilesDlgPtr,
          const utils::cmdLineOptions & options);
  
public slots:

  // View menu
  void zoomOut();
  void zoomIn();
  void shiftLeft();
  void shiftRight();
  void shiftUp();
  void shiftDown();
  void resetView();
  void changeOrder();
  bool getStringVectorFromGui(std::string title,
                              std::string description,
                              std::vector<std::string> & values);
  // Edit menu

  // Transform menu

  // Options menu
  void setLineWidth();
  void setBgColor();

  // Right-click menu
  void refreshPixmap();
  
protected:

  void paintEvent(QPaintEvent *);
  void resizeEvent(QResizeEvent*);
  void popUp(std::string msg);
  bool getStringFromGui(std::string title, std::string description,
                        std::string inputStr,
                        std::string & outputStr // output
                        );
  bool getRealValuesFromGui(// Inputs
                            std::string title,
                            std::string description,
                            const std::vector<double> & inputVec,
                            // Outputs
                            std::vector<double> & values
                            );
  void setBgFgColorsFromPrefs();
  bool eventFilter(QObject *obj, QEvent *E);
  void mousePressEvent( QMouseEvent *E);
  void mouseMoveEvent( QMouseEvent *E);
  void keyPressEvent( QKeyEvent *K );
  void mouseReleaseEvent ( QMouseEvent * E );
  bool isShiftLeftMouse(QMouseEvent * E);
  void wheelEvent(QWheelEvent *E);
  void contextMenuEvent(QContextMenuEvent *E);
  void chooseFilesToShow();

private slots:
 void showFilesChosenByUser();
  
private:
  void setupViewingWindow();
  void readAllPolys();
  void setupDisplayOrder(// Inputs
                         int                 numPolys,
                         // Input-output
                         bool              & changeDisplayOrder,
                         std::vector<int>  & polyVecOrder
                         );
  bool readOnePoly(// inputs
                   std::string          & filename,
                   bool                   plotPointsOnly,
                   // output
                   utils::dPoly & poly           
                   );
  void centerViewAtPoint(double x, double y);
  void updateRubberBand(QRect & R);

  void displayData( QPainter *paint );
  void plotDPoly(bool plotPoints, bool plotEdges,
                 bool plotFilled,
                 int lineWidth,
                 int drawVertIndex, // 0 is a good choice here
                 QPainter *paint,
                 utils::dPoly currPoly // Make a local copy on purpose
                 );
  void resetTransformSettings();
  void pixelToWorldCoords(int px, int py,
                          double & wx, double & wy);
  void worldToPixelCoords(double wx, double wy,
                          int & px,  int & py);
  
  
  void setStandardCursor();
  void setPolyDrawCursor();
  
  double m_zoomFactor, m_shiftX, m_shiftY;
  int m_mousePrsX,  m_mousePrsY, m_mouseRelX,  m_mouseRelY;
  int m_screenXll,  m_screenYll, m_screenWidX, m_screenWidY;
  double m_viewXll, m_viewYll,   m_viewWidX,   m_viewWidY;
  double m_screenRatio, m_pixelSize;
  
  // Polygons
  std::vector<utils::dPoly>       m_polyVec;

  std::vector<utils::polyOptions> m_polyOptionsVec; // options for exiting polygons
  utils::polyOptions m_prefs;                       // options for future polygons

  bool m_resetView;
  bool m_firstPaintEvent;
  
  // Use double buffering: draw to a pixmap first, refresh it only
  // if really necessary, and display it when paintEvent is called.
  QPixmap m_pixmap;

  QRect   m_emptyRubberBand;
  QRect   m_rubberBand;

  int m_showEdges, m_showPoints, m_showPointsEdges, m_toggleShowPointsEdges;
  bool m_changeDisplayOrder;

  bool m_zoomToMouseSelection, m_viewChanged;
  
  double m_menuX, m_menuY;

  std::vector<int> m_polyVecOrder;

  // Choose which files to hide/show in the GUI
  chooseFilesDlg  *     m_chooseFilesDlg;
  std::set<std::string> m_filesToHide;

};

#endif // ASPVIEW_H

