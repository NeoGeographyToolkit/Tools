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

/// Class for the central image panel of the GUI
class aspView : public QWidget{
  Q_OBJECT
public:

  /// Constructor
  aspView(QWidget *parent, chooseFilesDlg *chooseFilesDlgPtr,
          const utils::cmdLineOptions & options);
  
public slots:

  // View menu
  // - These commands all modify and redraw the current view
  void zoomOut();
  void zoomIn();
  void shiftLeft();
  void shiftRight();
  void shiftUp();
  void shiftDown();
  void resetView();
  void changeOrder();

  /// Fetch multiple words from the GUI
  bool getStringVectorFromGui(std::string title,
                              std::string description,
                              std::vector<std::string> & values);
  // Edit menu - TODO

  // Transform menu - TODO

  // Options menu
  void setLineWidth(); ///< Set line width and redraw
  void setBgColor();   ///< Set the background color and redraw

  // Right-click menu????
  void refreshPixmap();
  
protected:

  void paintEvent(QPaintEvent *);  ///< Redraw cached pixmap and rubber band
  void resizeEvent(QResizeEvent*); ///< Handle by calling refreshPixmap
  void popUp(std::string msg);     ///< Create a pop-up message box

  /// Fetch a string from the GUI
  bool getStringFromGui(std::string title, std::string description,
                        std::string inputStr,
                        std::string & outputStr // output
                        );
  /// Fetch multiple comma seperated doubles from the GUI
  bool getRealValuesFromGui(// Inputs
                            std::string title,
                            std::string description,
                            const std::vector<double> & inputVec,
                            // Outputs
                            std::vector<double> & values
                            );
  /// Load colors from preferences set when object was initialized.
  void setBgFgColorsFromPrefs();


  bool eventFilter       (QObject *obj, QEvent *E); ///< Just redirect to Qt
  void mousePressEvent   (QMouseEvent *E); ///< Record click location
  void mouseMoveEvent    (QMouseEvent *E); ///< Update bounding box (rubber band)
  void keyPressEvent     (QKeyEvent   *K); ///< Handle all key-only controls
  void mouseReleaseEvent (QMouseEvent *E); ///< Zoom in or zoom out depending on movement
  bool isShiftLeftMouse  (QMouseEvent *E); ///< Return true if both are true
  void wheelEvent        (QWheelEvent *E); ///< Pan if ctrl or shift pressed
  void contextMenuEvent  (QContextMenuEvent *E); ///< Pop up menu at click location?
  void chooseFilesToShow(); ///< Parses input from choose files panel.

private slots:
 /// Process user's choice from chooseFilesToShow() and redraw.
 void showFilesChosenByUser();
  
private:
  /// Init window to correct size and location. Happens at start of displayData().
  void setupViewingWindow();

  /// Create a dPoly object for each input file in m_polyOptionsVec.
  void readAllPolys();

  /// Decide the order in which polygons are displayed.
  void setupDisplayOrder(// Inputs
                         int                 numPolys,
                         // Input-output
                         bool              & changeDisplayOrder,
                         std::vector<int>  & polyVecOrder
                         );

  /// TODO: This should probably become a class function!
  bool readOnePoly(// inputs
                   std::string          & filename,
                   bool                   plotPointsOnly,
                   // output
                   utils::dPoly & poly           
                   );
  /// Redraw the view with a new center (in world coordinates)
  void centerViewAtPoint(double x, double y);

  /// Update drawing of the rubber band
  void updateRubberBand(QRect & R);

  /// Updates the window and draws all of the polygons.
  void displayData( QPainter *paint );


  /// Plot a given dPoly with given options.
  void plotDPoly(bool plotPoints, bool plotEdges,
                 bool plotFilled,
                 int lineWidth,
                 int drawVertIndex, // 0 is a good choice here
                 QPainter *paint,
                 utils::dPoly currPoly // Make a local copy on purpose
                 );
  /// Reset zoom factor and shift to 1.0 and 0,0 respectively.
  void resetTransformSettings();

  // Convert back and forth between panel pixel coordinates and image geo location.
  void pixelToWorldCoords(int      px, int      py,
                          double & wx, double & wy);
  void worldToPixelCoords(double wx, double wy,
                          int &  px,  int &  py);
  
  
  void setStandardCursor(); ///< An arrow
  void setPolyDrawCursor(); ///< A cross
  
  /// Higher zoom factor = wider field of view (zoomed out)
  double m_zoomFactor,  m_shiftX,    m_shiftY;
  int    m_mousePrsX,   m_mousePrsY, m_mouseRelX,  m_mouseRelY;
  int    m_screenXll,   m_screenYll, m_screenWidX, m_screenWidY;
  double m_viewXll,     m_viewYll,   m_viewWidX,   m_viewWidY;
  double m_screenRatio, m_pixelSize;
  
  /// Polygons - This stores info about all image coordinates
  std::vector<utils::dPoly> m_polyVec;

  /// This stores options and the associated image file for each polygon
  std::vector<utils::polyOptions> m_polyOptionsVec; // options for exiting polygons
  utils::polyOptions m_prefs;                       // options for future polygons

  bool m_resetView; ///< Flag to indicate that the setUpViewBox() should be called again.
  bool m_firstPaintEvent; ///< Flag is true if paintEvent() has not been called yet
  
  // Use double buffering: draw to a pixmap first, refresh it only
  // if really necessary, and display it when paintEvent is called.
  QPixmap m_pixmap;

  // The rubber band objects are only used for displaying while selecting.
  QRect   m_emptyRubberBand; ///< Permanent off-screen rubberBand
  QRect   m_rubberBand;      ///< Rubber band currently being drawn

  int m_showEdges, m_showPoints, m_showPointsEdges, m_toggleShowPointsEdges;
  bool m_changeDisplayOrder;

  bool m_zoomToMouseSelection, m_viewChanged;
  
  // TODO: Are these actually used for anything?
  double m_menuX, m_menuY;

  std::vector<int> m_polyVecOrder;

  // Choose which files to hide/show in the GUI
  chooseFilesDlg  *     m_chooseFilesDlg;
  std::set<std::string> m_filesToHide;

};

#endif // ASPVIEW_H

