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
#include <QPolygon>
#include <Q3PopupMenu>
#include <QContextMenuEvent>
#include <QEvent>
#include <QFileDialog>
#include <QHoverEvent>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPaintEvent>
#include <QStyleOptionFocusRect>
#include <QStylePainter>
#include <QTableWidget>
#include <QWheelEvent>
#include <cassert>
#include <cfloat>    // defines DBL_MAX
#include <cmath>
#include <cstdlib>
#include <iomanip>   // required for use of setw()
#include <iostream>
#include <qapplication.h>
#include <qcursor.h>
#include <qdir.h>
#include <qinputdialog.h>
#include <qpainter.h>
#include "aspView.h"
#include <qmessagebox.h>
#include "utils.h"
#include <vw/Image/Manipulation.h>

using namespace std;
using namespace utils;
using namespace vw;

// To do: handle colors correctly (convert dark-gray to darkGray, etc.).
// To do: In the geom directory, put everything in a namespace, say called 'pv'.
//        Here too. Clean up, modularize, and structure the code more.
// To do: Fix other "To do" mentioned in the code.
// To do: The viewer does not render correctly in fill mode overlapping polygons
//        with each polygon having holes. A fix would require a thorough analysis
//        which would identify which hole belongs to which polygon.
// To do: Replace cmdLineOptions directly with polyOptionsVec.
// To do: Make font size a preference.

aspView::aspView(QWidget *parent, chooseFilesDlg *chooseFilesDlgPtr,
                 const cmdLineOptions & options):
  QWidget(parent), m_chooseFilesDlg(chooseFilesDlgPtr){
  
  installEventFilter(this);

  // Choose which files to hide/show in the GUI
  if (m_chooseFilesDlg){
    QObject::connect(m_chooseFilesDlg->getFilesTable(),
                     SIGNAL(itemClicked(QTableWidgetItem *)),
                     this,
                     SLOT(showFilesChosenByUser())
                   );
  }
  
  setAttribute(Qt::WA_Hover); // To be able to do hovering
  
  // Preferences per polygon file. The element in the vector
  // m_polyOptionsVec below is not associated with any polygon
  // file. Set it apart, it will be used for new polygons.
  m_polyOptionsVec = options.polyOptionsVec;
  assert(m_polyOptionsVec.size() >= 1);
  m_prefs = m_polyOptionsVec.back(); m_polyOptionsVec.pop_back();
  m_prefs.plotAsPoints = false; // most likely the user wants to see edges not points
  setBgFgColorsFromPrefs(); // must be called after m_prefs is set
  
  setStandardCursor();

  // int
  m_screenXll  = 0; m_screenYll  = 0;
  m_screenWidX = 0; m_screenWidY = 0;

  // double
  m_viewXll  = 0.0; m_viewYll  = 0.0;
  m_viewWidX = 0.0; m_viewWidY = 0.0;

  m_resetView       = true;
  m_firstPaintEvent = true;
  
  m_changeDisplayOrder = false;
  
  m_emptyRubberBand = QRect(-10, -10, 0, 0); // off-screen rubberband
  m_rubberBand      = m_emptyRubberBand;

  m_showEdges               = 1;
  m_showPointsEdges         = 2;
  m_showPoints              = 3;
  m_toggleShowPointsEdges   = m_showEdges;

  m_zoomToMouseSelection = false;
  m_viewChanged          = false;

  m_zoomFactor = 1.0;
  m_mousePrsX = 0; m_mousePrsY = 0;
  m_mouseRelX = 0; m_mouseRelY = 0;
  
  resetTransformSettings();

  // This statement must be towards the end
  readAllPolys(); // To do: avoid global variables here

  chooseFilesToShow();
  
  return;
}

bool aspView::eventFilter(QObject *obj, QEvent *E){

  return QWidget::eventFilter(obj, E);
}

void aspView::setupViewingWindow(){

  // Dimensions of the plotting window in pixels excluding any window
  // frame/menu bar/status bar
  QRect v       = this->geometry();
  m_screenXll   = v.left();
  m_screenYll   = v.top();
  m_screenWidX  = v.width();
  m_screenWidY  = v.height();
  //cout << "geom is: " << m_screenXll << ' ' << m_screenYll << ' '
  //     << m_screenWidX << ' ' << m_screenWidY << endl;
  
  if (m_resetView){
    setUpViewBox(// inputs
                 m_polyVec,
                 // outputs
                 m_viewXll, m_viewYll, m_viewWidX, m_viewWidY
                 );
    m_resetView = false;
  }

  //  This is necessary when the screen is resized
  m_screenRatio = double(m_screenWidY)/double(m_screenWidX);
  expandBoxToGivenRatio(// Inputs
                        m_screenRatio,
                        // Inputs-outputs
                        m_viewXll, m_viewYll, m_viewWidX, m_viewWidY
                        );

  // Create the new view
  double xll, yll, xur, yur, widx, widy;
    
  if (m_zoomToMouseSelection){
    
    // Form a new view based on the rectangle selected with the mouse.
    // The call to pixelToWorldCoords uses the existing view internally.
    pixelToWorldCoords(m_mousePrsX, m_mousePrsY, xll, yur); // upper-left  rect corner
    pixelToWorldCoords(m_mouseRelX, m_mouseRelY, xur, yll); // lower-right rect corner
    widx = xur - xll;
    widy = yur - yll;

  }else if (m_viewChanged){

    // Modify the view for given shift or zoom
    xll  = m_viewXll + m_viewWidX*( (1 - m_zoomFactor)/2.0 + m_shiftX );
    yll  = m_viewYll + m_viewWidY*( (1 - m_zoomFactor)/2.0 + m_shiftY );
    widx = m_viewWidX*m_zoomFactor;
    widy = m_viewWidY*m_zoomFactor;
    
    resetTransformSettings(); // Wipe the zoom and shift data 
  
  }

  if (m_zoomToMouseSelection || m_viewChanged){
    
    // If the view becomes too small, don't accept it
    if (xll + widx <= xll || yll + widy <= yll){
      cerr << "Cannot zoom to requested view."  << endl;
    }else{
      // Enlarge this rectangle if necessary to keep the aspect ratio.
      expandBoxToGivenRatio(//inputs
                            m_screenRatio,  
                            // input/outputs
                            xll, yll, widx, widy
                            );
      
      // Overwrite the view
      m_viewXll = xll; m_viewWidX = widx;
      m_viewYll = yll; m_viewWidY = widy;
    }
    m_zoomToMouseSelection = false;
    m_viewChanged          = false;
  }

  // The two ratios below will always be the same. Take the maximum
  // for robustness to floating point errors.
  m_pixelSize = max(m_viewWidX/m_screenWidX, m_viewWidY/m_screenWidY);

  return;
}

void aspView::displayData( QPainter *paint ){

  setupViewingWindow(); // Must happen before anything else
  
  // Plot the polygons
  setupDisplayOrder(m_polyVec.size(),                    //inputs
                    m_changeDisplayOrder, m_polyVecOrder // inputs-outputs
                    );
  // Will draw a vertex with a shape dependent on this index
  int drawVertIndex = -1; 

  // Loop through each of the polygons
  assert( m_polyVec.size() == m_polyOptionsVec.size() );
  for (int vi  = 0; vi < (int)m_polyVec.size(); vi++){

    int vecIter = m_polyVecOrder[vi];

    // Skip the files the user does not want to see
    string fileName = m_polyOptionsVec[vecIter].polyFileName;
    if (m_filesToHide.find(fileName) != m_filesToHide.end()) continue;
      
    int lineWidth = m_polyOptionsVec[vecIter].lineWidth;

    // Note: plotFilled, plotEdges, and plotPoints are not mutually exclusive.
    bool plotFilled = m_polyOptionsVec[vecIter].isPolyFilled;
    bool plotEdges  = (!m_polyOptionsVec[vecIter].plotAsPoints) &&
      (m_toggleShowPointsEdges != m_showPoints);
    bool plotPoints = m_polyOptionsVec[vecIter].plotAsPoints  ||
      ( m_toggleShowPointsEdges == m_showPoints )             ||
      ( m_toggleShowPointsEdges == m_showPointsEdges);
    
    if (plotPoints) drawVertIndex++;

    plotDPoly(plotPoints, plotEdges, plotFilled, lineWidth, 
              drawVertIndex, paint, m_polyVec[vecIter]
              );
    
  } // End iterating over sets of polygons

  return;
}

void aspView::plotDPoly(bool plotPoints, bool plotEdges,
                         bool plotFilled,
                         int lineWidth,
                         int drawVertIndex, // 0 is a good choice here
                         QPainter *paint,
                         dPoly currPoly // Make a local copy on purpose
                         ){

  // Plot a given dPoly with given options.
  
  // Clip the polygon a bit beyond the viewing window, as to not see
  // the edges where the cut took place. It is a bit tricky to
  // decide how much the extra should be.
  double tol    = 1e-12;
  double extra  = 2*m_pixelSize*lineWidth;
  double extraX = extra + tol*max(abs(m_viewXll), abs(m_viewXll + m_viewWidX));
  double extraY = extra + tol*max(abs(m_viewYll), abs(m_viewYll + m_viewWidY));

  dPoly clippedPoly;
  currPoly.clipPoly(//inputs
                    m_viewXll - extraX,
                    m_viewYll - extraY,
                    m_viewXll + m_viewWidX + extraX,
                    m_viewYll + m_viewWidY + extraY,
                    // output
                    clippedPoly
                    );

  const double * xv               = clippedPoly.get_xv();
  const double * yv               = clippedPoly.get_yv();
  const int    * numVerts         = clippedPoly.get_numVerts();
  int numPolys                    = clippedPoly.get_numPolys();
  const vector<string> colors     = clippedPoly.get_colors();
  //int totalNumVerts             = clippedPoly.get_totalNumVerts();

  // To do: This must go somewhere in the preferences.
  setBackgroundColor(QColor("white"));
  m_prefs.fgColor = "green";
  
  BBox2i wBox, pBox; // box in world and pixel coordinates
  int x0, y0;
  int start = 0;
  for (int pIter = 0; pIter < numPolys; pIter++){ // For each polygon
      
    if (pIter > 0) start += numVerts[pIter - 1];
      
    int pSize = numVerts[pIter];
    if (pSize <= 0) continue;

    // Expand bounding boxes to contain each vertex
    for (int vIter = 0; vIter < pSize; vIter++){

      worldToPixelCoords(xv[start + vIter], yv[start + vIter], // inputs
                         x0, y0                                // outputs
                         );
      wBox.grow(Vector2(xv[start + vIter], yv[start + vIter]));
      pBox.grow(Vector2(x0, y0));
    }

    // Compensate for the fact that an image's origin is the upper
    // left corner.
    double a = wBox.min().y(), b = wBox.max().y();
    wBox.min().y() = currPoly.m_img.rows() - b;
    wBox.max().y() = currPoly.m_img.rows() - a;
    wBox.crop(bounding_box(currPoly.m_img));
      
    //TODO: This should be turned into a function!
    // Rasterize a VW image into a QT image in the requested area
    ImageView<float> img = crop(currPoly.m_img, wBox);
    QImage qimg(img.cols(), img.rows(), QImage::Format_RGB888);
    for (int x = 0; x < img.cols(); ++x) {
      for (int y = 0; y < img.rows(); ++y) {
        qimg.setPixel(x, y, qRgb(img(x, y), img(x, y), img(x, y)));
      }
    }
      
    // Have QT draw the image in the correct location
    QRect rect(pBox.min().x(), pBox.min().y(),
               pBox.width(), pBox.height());
    paint->drawImage (rect, qimg);
  }
  return;
}
  
void aspView::zoomIn(){
  m_zoomFactor  = 0.5;
  m_viewChanged = true;
  refreshPixmap();
}

void aspView::zoomOut(){
  m_zoomFactor  = 2.0;
  m_viewChanged = true;
  refreshPixmap();
}

void aspView::shiftRight(){
  m_shiftX      = 0.25;
  m_viewChanged = true;
  refreshPixmap();
}

void aspView::shiftLeft(){
  m_shiftX      = -0.25;
  m_viewChanged = true;
  refreshPixmap();
}

void aspView::shiftUp(){
  m_shiftY      = 0.25;
  m_viewChanged = true;
  refreshPixmap();
}

void aspView::shiftDown(){
  m_shiftY      = -0.25;
  m_viewChanged = true;
  refreshPixmap();
}

void aspView::centerViewAtPoint(double x, double y){
  m_viewXll     = x - m_viewWidX/2.0;
  m_viewYll     = y - m_viewWidY/2.0;
  m_viewChanged = true;
}

void aspView::resetView(){
  m_resetView   = true;
  m_viewChanged = true;
  refreshPixmap();
}


void aspView::resetTransformSettings(){
  m_zoomFactor = 1.0;
  m_shiftX     = 0.0; m_shiftY = 0.0;
}

void aspView::mousePressEvent( QMouseEvent *E){

  const QPoint Q = E->pos();
  m_mousePrsX = Q.x();
  m_mousePrsY = Q.y();

#if 0
  cout << "Mouse pressed at "
       << m_mousePrsX << ' ' << m_mousePrsY << endl;
#endif

  // Init rubber band to empty status
  m_rubberBand = m_emptyRubberBand;

  return;
}

void aspView::mouseMoveEvent( QMouseEvent *E){

  QPoint Q = E->pos();
  int x = Q.x(), y = Q.y();

  double wx, wy;
  pixelToWorldCoords(x, y, wx, wy);

  // Redraw the rubber band at the current mouse location
  refreshPixmap(); // To do: Need to update just a small region, not the whole screen
  // Standard Qt rubberband trick (kind of confusing as to how it works).
  updateRubberBand(m_rubberBand);
  m_rubberBand = QRect( min(m_mousePrsX, x), min(m_mousePrsY, y),
                        abs(x - m_mousePrsX), abs(y - m_mousePrsY) );
  updateRubberBand(m_rubberBand);
  
  return;
}

void aspView::mouseReleaseEvent ( QMouseEvent * E ){

  const QPoint Q = E->pos();
  m_mouseRelX = Q.x();
  m_mouseRelY = Q.y();
#if 0 // Debug code
  cout << "Mouse pressed at "
       << m_mousePrsX << ' ' << m_mousePrsY << endl;
  cout << "Mouse released at "
       << m_mouseRelX << ' ' << m_mouseRelY << endl;
#endif
    
  pixelToWorldCoords(m_mouseRelX, m_mouseRelY, m_menuX, m_menuY);
  
  // Wipe the rubberband
  updateRubberBand(m_rubberBand); 
  m_rubberBand = m_emptyRubberBand;
  updateRubberBand(m_rubberBand);

  if ( E->modifiers() & Qt::ControlModifier ){
    // Control Mouse is being pressed
    refreshPixmap();
    return;
  }

  int tol = 10; 
  // Any selection smaller than 'tol' number of pixels will be ignored
  // as perhaps the user moved the mouse unintentionally between press
  // and release.
  if       (m_mouseRelX > m_mousePrsX + tol &&
            m_mouseRelY > m_mousePrsY + tol){
    
    m_zoomToMouseSelection = true; // Will zoom to the region selected with the mouse
    refreshPixmap(); 
    return;
    
  }else if (m_mouseRelX + tol < m_mousePrsX &&
            m_mouseRelY + tol < m_mousePrsY ){
    
    zoomOut();
    return;
    
  }else if (abs(m_mouseRelX - m_mousePrsX) <= tol &&
            abs(m_mouseRelY - m_mousePrsY) <= tol){
    return; // Do nothing when the selection is too small.
  }
  
  return;
}

bool aspView::isShiftLeftMouse(QMouseEvent * E){
  // This does not work in mouseReleaseEvent. 
  return ( E->buttons() & Qt::LeftButton  ) && ( E->modifiers() & Qt::ShiftModifier );
}

void aspView::wheelEvent(QWheelEvent *E){

  int delta = E->delta();

  if (E->state() == Qt::ControlModifier){

    // The control button was pressed. Zoom in/out around the current point.

    int pixelx = E->x();
    int pixely = E->y();
    double x, y;
    pixelToWorldCoords(pixelx, pixely, x, y);
    centerViewAtPoint(x, y);
    
    if (delta > 0){
      zoomIn();
    }else if (delta < 0){
      zoomOut();
    }
    
  }else{ // Control not pressed

    // Shift wheel goes left and right. Without shift we go up and down.
    if (E->state() == Qt::ShiftModifier){
      if (delta > 0){
        shiftLeft();
      }else if (delta < 0){
        shiftRight();
      }
    }else{
      if (delta > 0){
        shiftUp();
      }else if (delta < 0){
        shiftDown();
      }
    } // End Shift pressed
  } // End control not pressed case
  
  E->accept();
}

void aspView::keyPressEvent( QKeyEvent *K ){

  switch ( K->key() ) {
  case Qt::Key_Minus:
    zoomOut();
    break;
  case Qt::Key_Plus:
    zoomIn();
    break;
  case Qt::Key_Equal:
    zoomIn();
    break;
  case Qt::Key_Right:
    shiftRight();
    break;
  case Qt::Key_Left:
    shiftLeft();
    break;
  case Qt::Key_Up:
    shiftUp();
    break;
  case Qt::Key_Down:
    shiftDown();
    break;
  }
  
}

void aspView::contextMenuEvent(QContextMenuEvent *E){

  int x = E->x(), y = E->y();
  pixelToWorldCoords(x, y, m_menuX, m_menuY);

  Q3PopupMenu menu(this);
  menu.exec(E->globalPos());

  return;
}


void aspView::refreshPixmap(){

  // Draw the data onto the pixmap instead of the screen
  // directly. Later we'll display the pixmap without redrawing
  // whenever possible for reasons of speed.

  m_pixmap = QPixmap(size());
  m_pixmap.fill(this, 0, 0);
  
  QPainter paint(&m_pixmap);
  paint.initFrom(this);

  QFont F;
  F.setPointSize(m_prefs.fontSize);
  //F.setStyleStrategy(QFont::NoAntialias);
  paint.setFont(F);

  displayData( &paint );
  update();

  return;
}

void aspView::paintEvent(QPaintEvent *){

  if (m_firstPaintEvent){
    // This will be called the very first time the display is
    // initialized. There must be a better way.
    m_firstPaintEvent = false;
    refreshPixmap(); // Create m_pixmap which we will use as cache
  }

  // Note that we draw from the cached pixmap, instead of redrawing
  // the image from scratch.
  QStylePainter paint(this);
  paint.drawPixmap(0, 0, m_pixmap);

  // Drawing the rubber band
  QColor fgColor = QColor(m_prefs.fgColor.c_str());
  paint.setPen(fgColor);
  paint.drawRect(m_rubberBand.normalized().adjusted(0, 0, -1, -1));

  return;
}

void aspView::resizeEvent(QResizeEvent*){
  refreshPixmap();
  return;
}

void aspView::popUp(std::string msg){
  QMessageBox msgBox;
  msgBox.setText(msg.c_str());
  msgBox.exec();
  return;
}

bool aspView::getStringFromGui(std::string title, std::string description,
                               std::string inputStr,
                               std::string & outputStr // output
                               ){

  outputStr = "";

  bool ok = false;
  QString text = QInputDialog::getText(title.c_str(), description.c_str(),
                                       QLineEdit::Normal, inputStr.c_str(),
                                       &ok, this );

  if (ok) outputStr = text.toStdString();

  return ok;
}

bool aspView::getRealValuesFromGui(// Inputs
                                    std::string title,
                                    std::string description,
                                    const std::vector<double> & inputVec,
                                    // Outputs
                                    std::vector<double> & values
                                    ){

  values.clear();
  
  string outputStr;

  ostringstream oss;
  oss.precision(16);
  int len = (int)inputVec.size();
  for (int s = 0; s < len - 1; s++) oss << inputVec[s] << " ";
  if (len > 0) oss << inputVec[len - 1];
  string inputStr = oss.str();

  bool ok = getStringFromGui(title, description, inputStr,
                             outputStr
                             );
  
  outputStr = replaceAll(outputStr, ",", " ");
  istringstream ts(outputStr);
  double val;
  while (ts >> val) values.push_back(val);

  return ok;
}

bool aspView::getStringVectorFromGui(std::string title,
                                     std::string description,
                                     std::vector<std::string> & values){

  values.clear();
  string inputStr, outputStr;
  bool ok = getStringFromGui(title, description, inputStr,
                             outputStr // output
                             );
  
  istringstream ts(outputStr);
  string val;
  while (ts >> val) values.push_back(val);

  return ok;
}

void aspView::setLineWidth(){

  vector<double> inputVec, lineWidth;
  if ( !getRealValuesFromGui("Line width", "Enter line width", inputVec,
                             lineWidth) ) return;
  
  if ( !lineWidth.empty() && lineWidth[0] >= 1.0 ){

    int lw = (int) round(lineWidth[0]);

    for (int polyIter = 0; polyIter < (int)m_polyOptionsVec.size(); polyIter++){
      m_polyOptionsVec[polyIter].lineWidth = lw;
    }
    m_prefs.lineWidth = lw;
    
    refreshPixmap();
  
  }else{
    popUp("The line width must be a positive integer.");
  }
  return;
}

void aspView::setBgColor(){

  vector<string> values;
  if (!getStringVectorFromGui("Background", "Enter background color",
                              values)) return;

  string bgColor = "";
  if (values.size() > 0) bgColor = values[0];
  if ( QColor(bgColor.c_str()) != QColor::Invalid){
    m_prefs.bgColor = bgColor;
    setBgFgColorsFromPrefs();
    refreshPixmap();
  }else{
    popUp("Invalid background color.");
  }
  return;
}

void aspView::setBgFgColorsFromPrefs(){

  // Set the background. Watch for invalid colors.
  string bgColor   = m_prefs.bgColor;
  QColor qtBgColor = QColor(bgColor.c_str());
  if (qtBgColor == QColor::Invalid){
    bgColor   = "black";
    qtBgColor = QColor(bgColor.c_str()); // fallback color
  }
  setBackgroundColor(qtBgColor);

  string fgColor = m_prefs.fgColor;
  if ( QColor(fgColor.c_str()) == QColor::Invalid ){
    fgColor = "white";
  }

  // Make sure bg and fg have different colors
  if (bgColor == fgColor){
    if (bgColor == "black") fgColor = "white";
    else                    fgColor = "black";
  }

  // Update the preferences after doing the logic above
  m_prefs.bgColor = bgColor;
  m_prefs.fgColor = fgColor;

  // While unnecessary, also update the preferences for each
  // polygon file, for consistency with other preferences per file.
  for (int s = 0; s < (int)m_polyOptionsVec.size(); s++){
    m_polyOptionsVec[s].bgColor = bgColor; 
    m_polyOptionsVec[s].fgColor = fgColor; 
  }

  return;
}

void aspView::updateRubberBand(QRect & R){
  
  QRect rect = R.normalized();
  update(rect.left(), rect.top(),    rect.width(), 1             );
  update(rect.left(), rect.top(),    1,            rect.height() );
  update(rect.left(), rect.bottom(), rect.width(), 1             );
  update(rect.right(), rect.top(),   1,            rect.height() );
  
  return;
}

void aspView::pixelToWorldCoords(int      px, int      py,
                                 double & wx, double & wy){

  // Compensate for the Qt's origin being in the upper-left corner
  // instead of the lower-left corner.
  py = m_screenWidY - py;

  wx = px*m_pixelSize + m_viewXll;
  wy = py*m_pixelSize + m_viewYll;

}

void aspView::worldToPixelCoords(double wx, double wy,
                                  int & px,  int & py){

  px = iround((wx - m_viewXll)/m_pixelSize);
  py = iround((wy - m_viewYll)/m_pixelSize);
  
  // Compensate for the Qt's origin being in the upper-left corner
  // instead of the lower-left corner.
  py = m_screenWidY - py;
  
}

void aspView::readAllPolys(){

  int numFiles = m_polyOptionsVec.size();
  m_polyVec.resize(numFiles);

  string missingFiles = "";
  int numMissing = 0;
  
  for (int fileIter = 0; fileIter < numFiles; fileIter++){

    bool success = readOnePoly(// inputs
                               m_polyOptionsVec[fileIter].polyFileName,
                               m_polyOptionsVec[fileIter].plotAsPoints,
                               // output
                               m_polyVec[fileIter]
                               );
    if (!success){
      missingFiles += " " + m_polyOptionsVec[fileIter].polyFileName;
      numMissing++;
    }
    
  }


  if (numMissing >= 1){
    string suffix = ""; if (numMissing > 1) suffix = "s";
    popUp("Warning: Could not read file" + suffix + ":" + missingFiles + ".");
  }
  
  return;
}

void aspView::chooseFilesToShow(){
  // User's choice is processed in showFilesChosenByUser().
  if (m_chooseFilesDlg)
    m_chooseFilesDlg->chooseFiles(m_polyOptionsVec);
  return;
}

void aspView::showFilesChosenByUser(){
  
  // Process user's choice from chooseFilesToShow().

  if (!m_chooseFilesDlg)
    return;

  m_filesToHide.clear();
  QTableWidget * filesTable = m_chooseFilesDlg->getFilesTable();
  int rows = filesTable->rowCount();

  for (int rowIter = 0; rowIter < rows; rowIter++){
    QTableWidgetItem *item = filesTable->item(rowIter, 0);
    if (item->checkState() != Qt::Checked){
      string fileName
        = (filesTable->item(rowIter, 1)->data(0)).toString().toStdString();
      m_filesToHide.insert(fileName);
    }
  }
  
  refreshPixmap();
  
  return;
}

bool aspView::readOnePoly(// inputs
                           std::string   & filename,
                           bool            plotPointsOnly,
                           // output
                           dPoly & poly           
                           ){

  poly.reset();
  
  string type = getFilenameExtension(filename);
  
  if ( ! poly.readPoly(filename, plotPointsOnly) ){
    return false;
  }

  return true;
}


void aspView::changeOrder(){

  m_changeDisplayOrder = true;
  refreshPixmap();
  
}

void aspView::setStandardCursor(){
  QCursor C(Qt::ArrowCursor);
  setCursor(C);
}

void aspView::setPolyDrawCursor(){
  QCursor C(Qt::CrossCursor);
  setCursor(C);
}

void aspView::setupDisplayOrder(// Inputs
                                 int                 numPolys,
                                 // Input-output
                                 bool              & changeDisplayOrder,
                                 std::vector<int>  & polyVecOrder
                                 ){


  // Decide the order in which polygons are displayed. 

  if ((int)polyVecOrder.size() != numPolys){

    // Default order
    polyVecOrder.resize(numPolys);
    for (int c = 0; c < numPolys; c++){
      polyVecOrder[c] = c;
    }

  }else if (changeDisplayOrder && numPolys >= 1){

    changeDisplayOrder = false;

    // Cycle left
    int bk = polyVecOrder[0];
    for (int c = 1; c < numPolys; c++) polyVecOrder[c-1] = polyVecOrder[c];
    polyVecOrder[numPolys - 1] = bk;
    
  }

  return;
}
