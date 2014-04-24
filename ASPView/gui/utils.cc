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
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <cstring>
#include <algorithm>
#include <cfloat>
#include "utils.h"

#include <vw/FileIO/DiskImageView.h>
using namespace std;
using namespace utils;
using namespace vw;

#define DEBUG_CUT_POLY 0 // Must be 0 in production code

void utils::expandBoxToGivenRatio(// inputs
                                  double aspectRatio, 
                                  // inputs/outputs
                                  double & xll,  double & yll,
                                  double & widx, double & widy){

  // Expand the given box to have the aspect ratio equal to the number aspectRatio.
  assert(widx > 0.0 && widy > 0.0 && aspectRatio > 0.0);
  double nwidx = widx, nwidy = widy;
  if (widy/widx <= aspectRatio) nwidy = widx*aspectRatio;
  else                          nwidx = widy/aspectRatio;

  // Sanity checks
  double tol = 1.0e-3;
  bool check = ( nwidx >= widx*(1 - tol) && nwidy >= widy*(1 - tol)
                 && abs(nwidy/nwidx - aspectRatio) < tol*aspectRatio );
  if (!check){
    cout << "ERROR!" << endl;
    cout << "widx widy are "   << widx  << ' ' << widy  << endl;
    cout << "nwidx nwidy are " << nwidx << ' ' << nwidy << endl;
    cout << "Aspect ratio is " << aspectRatio << endl;
    cout << "|nwidy/nwidx - aspectRatio| = " << abs(nwidy/nwidx - aspectRatio) << endl;
    cout << "Max allowed error is " << tol*aspectRatio << endl;
  }
  assert(check);

  // Make the new bounding box have the same center as the old one
  xll += widx/2.0 - nwidx/2.0;
  yll += widy/2.0 - nwidy/2.0;

  // Overwrite the previous box
  widx = nwidx; 
  widy = nwidy;

  return;
}

void utils::cutEdge(double x0, double y0, double x1, double y1,
                    double nx, double ny, double H,
                    double & cutx, double & cuty){

  // Find the point at which the horizontal or vertical line
  // nx*x + ny*y = H intersects the edge (x0, y0) --> (x1, y1).
  // We assume that by now we already know that the edge
  // intersects this line.

  // To do: This needs careful reading.

  double dot0 = nx*x0 + ny*y0;
  double dot1 = nx*x1 + ny*y1;
  
  assert( (dot0 <= H && dot1 >= H) || (dot0 >= H && dot1 <= H ) );
  assert( dot0 != dot1 );

  // Find t such that (1-t)*(x0, y0) + t*(x1, y1) intersect the
  // cutting line
  
  double t = (H - dot0)/(dot1 - dot0);
  t = max(t, 0.0); t = min(t, 1.0); // extra precautions
  cutx = (1-t)*x0 + t*x1;
  cuty = (1-t)*y0 + t*y1;

  // Cut symmetrically in x0 and x1. We count on this symmetry
  // in a few places.
  t = (H - dot1)/(dot0 - dot1);
  t = max(t, 0.0); t = min(t, 1.0); // extra precautions
  double cutx2 = (1-t)*x1 + t*x0;
  double cuty2 = (1-t)*y1 + t*y0;

  cutx = 0.5*cutx + 0.5*cutx2;
  cuty = 0.5*cuty + 0.5*cuty2;
  
  // The above formulas have floating point errors. Use the fact that
  // the cutting line is either vertical or horizontal to recover
  // precisely one of the two coordinates above.
  if (nx == 0){
    cuty = H/ny;
  }else if (ny == 0){
    cutx = H/nx;
  }
  
  return;
}

void utils::cutPolyLine(// inputs -- the polygonal line
                        int numVerts,
                        const double * xv, const double * yv,
                        // inputs -- the cutting window
                        double xll, double yll, double xur, double yur,
                        // outputs -- the cut polygons
                        std::vector< double> & cutX,
                        std::vector< double> & cutY,
                        std::vector< int>    & cutNumPolys){

  // Cut a polygonal line. First make it into a polygon by traveling
  // forward and then backward on the polygonal line, then cut the
  // obtained polygon, then remove the backward portion from each
  // obtained polygon.

  vector<double> lXv, lYv, lCutX, lCutY;
  vector<int> lCutNumPolys;
  
  lXv.clear(); lYv.clear();
  for (int s = 0; s < numVerts; s++){
    lXv.push_back(xv[s]);
    lYv.push_back(yv[s]);
  }
  for (int s = numVerts - 1; s >= 0; s--){
    lXv.push_back(xv[s]);
    lYv.push_back(yv[s]);
  }
  int lNumVerts = lXv.size();
  
  cutPoly(// inputs -- the polygons
          1, &lNumVerts,  
          vecPtr(lXv), vecPtr(lYv),  
          // inputs -- the cutting window
          xll, yll, xur, yur,  
          // outputs -- the cut polygons
          lCutX, lCutY, lCutNumPolys
          );

  cutX.clear(); cutY.clear(); cutNumPolys.clear();
  
  int start = 0;
  for (int pIter = 0; pIter < (int)lCutNumPolys.size(); pIter++){

    if (pIter > 0) start += lCutNumPolys[pIter - 1];

    // Keep only half of the points of the cut polygon
    int half = lCutNumPolys[pIter]/2;
    cutNumPolys.push_back(half);
    for (int vIter = 0; vIter < half; vIter++){
      cutX.push_back(lCutX[start + vIter]);
      cutY.push_back(lCutY[start + vIter]);
    }
    
  }

  return;
}

void utils::cutPoly(// inputs -- the polygons
                    int numPolys, const int * numVerts,
                    const double * xv, const double * yv,
                    // inputs -- the cutting window
                    double xll, double yll, double xur, double yur,
                    // outputs -- the cut polygons
                    std::vector< double> & cutX,
                    std::vector< double> & cutY,
                    std::vector< int>    & cutNumPolys){
  
  // Cut a given polygon with a box.
  
  // Intersect the polygon with each of the the half-planes
  // nx*x + ny*y <= (nx + ny)*H.
  // There are four values for the triplet (nx, ny, H):
  double cutParams[] = {
   -1,  0, xll, //  -- left cut
    1,  0, xur, //  -- right cut
    0, -1, yll, //  -- bottom cut
    0,  1, yur  //  -- top cut
  };

  int totalNumVerts = 0;
  for (int s = 0; s < numPolys; s++) totalNumVerts += numVerts[s];

  vector<double> Xin(xv, xv + totalNumVerts);        // A copy of xv as vector
  vector<double> Yin(yv, yv + totalNumVerts);        // A copy of yv as vector
  vector<int>    Pin(numVerts, numVerts + numPolys); // A copy of numVerts 
  
  vector<double> cutHalfX, cutHalfY, Xout, Yout;
  vector<int>    cutHalfP, Pout;

  for (int c = 0; c < 4; c++){

    Pout.clear(); Xout.clear(); Yout.clear();
    
    double nx   = cutParams[3*c + 0];
    double ny   = cutParams[3*c + 1];
    double H    = cutParams[3*c + 2];
    double dotH = (nx + ny)*H; // This formula works only for nx*ny == 0.

    int start = 0;
    for (int pIter = 0; pIter < (int)Pin.size(); pIter++){
      
      if (pIter > 0) start += Pin[pIter - 1];
      
      int numV = Pin[pIter];
      if (numV == 0) continue;
      
      cutToHalfSpace(nx, ny, dotH,
                     numV, vecPtr(Xin) + start, vecPtr(Yin) + start,
                     cutHalfX, cutHalfY, cutHalfP);
      
      for (int pIterCut = 0; pIterCut < (int)cutHalfP.size(); pIterCut++){
        if (cutHalfP[pIterCut] > 0){
          // Append only non-empty polygons
          Pout.push_back( cutHalfP[pIterCut] );
        }
      }
      
      for (int vIter = 0; vIter < (int)cutHalfX.size(); vIter++){
        Xout.push_back( cutHalfX[vIter] );
        Yout.push_back( cutHalfY[vIter] );
      }
      
    }
    
    Pin = Pout; Xin = Xout; Yin = Yout;
    
  } // End iterating over cutting lines

  cutNumPolys = Pout; cutX = Xout; cutY = Yout;
  
  return;
}

void utils::cutToHalfSpace(// inputs 
                           double nx, double ny, double dotH,
                           int numV, 
                           const double * xv, const double * yv,
                           // outputs -- the cut polygons
                           std::vector<double> & cutX,
                           std::vector<double> & cutY,
                           std::vector<int>    & cutNumPolys){


  vector<valIndex> ptsOnCutline; ptsOnCutline.clear();
  valIndex C;
  
  cutX.clear(); cutY.clear(); cutNumPolys.clear();

  int cutPtsIndex = 0;
  
  for (int v = 0; v < numV; v++){

    int vnext = (v + 1)%numV;
    
    double xcurr = xv[v],     ycurr = yv[v];
    double xnext = xv[vnext], ynext = yv[vnext];

    double dotCurr = nx*xcurr + ny*ycurr;
    double dotNext = nx*xnext + ny*ynext;
    double cutx = 0.0, cuty = 0.0;
    
    if (dotCurr < dotH){

      // The current point is inside the half-plane
      
      cutX.push_back(xcurr);
      cutY.push_back(ycurr);
      cutPtsIndex++; 

      if (dotNext <= dotH) continue;
      
      cutEdge(xcurr, ycurr, xnext, ynext, nx, ny, dotH, cutx, cuty);

      cutX.push_back(cutx);
      cutY.push_back(cuty);

      C.val       = nx*cuty - ny*cutx;
      C.index     = cutPtsIndex;
      C.isOutward = true;
      ptsOnCutline.push_back(C);

      cutPtsIndex++; 
      
    }else if (dotCurr > dotH){

      // The current point is outside the half-plane
      
      if (dotNext >= dotH) continue;

      cutEdge(xcurr, ycurr, xnext, ynext, nx, ny, dotH, cutx, cuty);
      
      cutX.push_back(cutx);
      cutY.push_back(cuty);
      
      C.val       = nx*cuty - ny*cutx;
      C.index     = cutPtsIndex;
      C.isOutward = false;
      ptsOnCutline.push_back(C);

      cutPtsIndex++; 

    }else if (dotCurr == dotH){

      // The current point is at the edge of the half-plane

      int    vprev   = (v == 0) ? (numV - 1) : (v - 1);
      double xprev   = xv[vprev], yprev = yv[vprev];
      double dotPrev = nx*xprev + ny*yprev;
      
      if (dotPrev >= dotH && dotNext >= dotH) continue;
      
      cutX.push_back(xcurr);
      cutY.push_back(ycurr);

      if (dotPrev >= dotH || dotNext >= dotH){

        C.val       = nx*ycurr - ny*xcurr;
        C.index     = cutPtsIndex;
        C.isOutward = (dotPrev < dotH);
        ptsOnCutline.push_back(C);
        
      }

      cutPtsIndex++; 
      
    }
    
  }
  
  
  int numPtsOnCutline = ptsOnCutline.size();
  if (numPtsOnCutline == 0){
    cutNumPolys.push_back( cutX.size() );
    return;
  }

  processPointsOnCutline(ptsOnCutline);
  
  // Find the connected components in the cut polygons
  // To do: Move this to its own function.
  
  vector<double> X, Y;
  vector<int> P;
  X.clear(); Y.clear(); P.clear();

#if DEBUG_CUT_POLY
  static int c = -1;
  c++;
  char file[100];
  sprintf(file, "beforeCleanup%d.xg", c);
  cout << "\nWriting to " << file << endl;
  ofstream before(file);
  for (int s = 0; s < (int)cutX.size(); s++){
    before << cutX[s] << ' ' << cutY[s] << endl;
  }
  before.close();

  for (int s = 0; s < (int)ptsOnCutline.size(); s++){
    cout.precision(20);
    cout << "point on cutline is (index outward val) "
         << ptsOnCutline[s].index     << ' '
         << ptsOnCutline[s].isOutward << ' '
         << ptsOnCutline[s].val       << endl; 
  }
#endif

  vector<int> wasVisited;
  int numCutPts = cutX.size();
  wasVisited.assign(numCutPts, 0);

  int ptIter = 0;
  while(1){

    // Stop when all points are visited
    bool success = false;
    for (int v = 0; v < numCutPts; v++){
      if (!wasVisited[v]){
        ptIter  = v;
        success = true;
        break;
      }
    }
    if (!success) break; 

    int numPtsInComp = 0;
    
    // Visit a given connected component
    while(1){
      
      if (wasVisited[ptIter]){
        P.push_back(numPtsInComp);
        break; // Arrived back to the starting point of the given component
      }
      
      X.push_back(cutX[ptIter]);
      Y.push_back(cutY[ptIter]);
      wasVisited[ptIter] = 1;
#if DEBUG_CUT_POLY
      cout << "ptIter = " << ptIter << endl;
#endif
      numPtsInComp++;

      // Decide which point we will visit next
      
      if (nx*cutX[ptIter] + ny*cutY[ptIter] != dotH){
        // The point is not at the cutline
        ptIter = (ptIter + 1)%numCutPts;
        continue;
      }

      // The point is at the cutline. Find where exactly it is in the
      // sorted cutline points. If it is not among those sorted
      // cutline points, it means that the polygon only touches the
      // cutline at that point rather than crossing over to the other
      // side.
      // To do: Use here some faster lookup, such as a map.
      bool success = false;
      int cutlineIter = 0;
      for (cutlineIter = 0; cutlineIter < numPtsOnCutline; cutlineIter++){
        if (ptsOnCutline[cutlineIter].index == ptIter){
          success = true;
          break;
        }
      }
      if (!success){
        ptIter = (ptIter + 1)%numCutPts;
        continue;
      }

      if (ptsOnCutline[cutlineIter].isOutward){

        // We are getting out of the polygon. Find the point at
        // which we come back.
        ptIter = ptsOnCutline[cutlineIter].nextIndexInward;
        continue;
        
      }else{
        
        // The point ptIter is at the cutline on the way in.
        // The next point will be inside the current half-plane.
        ptIter = (ptIter + 1)%numCutPts;
        continue;
        
      }
      
    } // End iterating over all connected components

  } // End iterating over all points
  
  cutX        = X;
  cutY        = Y;
  cutNumPolys = P;

#if DEBUG_CUT_POLY
  sprintf(file, "afterCleanup%d.xg", c);
  cout << "Writing to " << file << endl;
  ofstream after(file);
  for (int s = 0; s < (int)cutX.size(); s++){
    after << cutX[s] << ' ' << cutY[s] << endl;
  }
  after.close();
#endif
  
}

void utils::processPointsOnCutline(std::vector<valIndex> & ptsOnCutline){

  
  // Sort the cutline points along the cutline (the sort direction
  // does not matter).
  sort( ptsOnCutline.begin(), ptsOnCutline.end(), lessThan );

  // Find the position of each outward and each inward point on the cutline
  vector<int> outwardPositions, inwardPositions;
  outwardPositions.clear(); inwardPositions.clear();
  int numPtsOnCutline = ptsOnCutline.size();
  for (int s = 0; s < numPtsOnCutline; s++){
    const valIndex & C = ptsOnCutline[s];
    if (C.isOutward) outwardPositions.push_back(s);
    else             inwardPositions.push_back (s);
    //cout << "val index isOutward "
    //     << C.val << ' ' << C.index << ' ' << C.isOutward << endl; 
  }

  // There must be an even number of points on a cutline. That holds
  // true for any closed curve, with or without self-intersections.
  // Each time we cross the cutline to the other side at some point
  // there must be a corresponding point at which we come back. Match
  // the i-th outward point to the corresponding i-th inward point.
  int numOut = outwardPositions.size();
#ifndef NDEBUG
  int numIn  = inwardPositions.size();
  assert(numIn == numOut);
#endif
  for (int i = 0; i < numOut; i++){

    valIndex & C = ptsOnCutline[outwardPositions[i]]; //alias
    
    if (C.isOutward){
      C.nextIndexInward = ptsOnCutline[inwardPositions[i]].index;
#if DEBUG_CUT_POLY
      cout << "Going from " << C.index << " to " << C.nextIndexInward << endl;
#endif
    }else{
      C.nextIndexInward = -1; // To not leave it uninitialized
    }
    
  }
  
  return;
}


// A double precision polygon class

void dPoly::reset(){
  m_isPointCloud  = false;
  m_numPolys      = 0;
  m_totalNumVerts = 0;
  m_numVerts.clear();
  m_xv.clear();
  m_yv.clear();
  m_isPolyClosed.clear();
  m_colors.clear();
}

void dPoly::bdBox(double & xll, double & yll, double & xur, double & yur) const{

  // The bounding box of all polygons
  
  if (m_totalNumVerts <= 0){
    xll = DBL_MAX/4.0, xur = -DBL_MAX/4.0; // Use 1/4.0 to avoid overflow when ...
    yll = DBL_MAX/4.0, yur = -DBL_MAX/4.0; // ... finding width and height
    return;
  }
    
  xll = *min_element( vecPtr(m_xv), vecPtr(m_xv) + m_totalNumVerts );
  yll = *min_element( vecPtr(m_yv), vecPtr(m_yv) + m_totalNumVerts );
  xur = *max_element( vecPtr(m_xv), vecPtr(m_xv) + m_totalNumVerts );
  yur = *max_element( vecPtr(m_yv), vecPtr(m_yv) + m_totalNumVerts );

  return;
};

void dPoly::bdBoxes(std::vector<double> & xll, std::vector<double> & yll,
                    std::vector<double> & xur, std::vector<double> & yur) const{

  // Bounding boxes of individual polygons
  
  xll.clear(); yll.clear(); xur.clear(); yur.clear();
  
  int start = 0;
  for (int pIter = 0; pIter < m_numPolys; pIter++){
      
    if (pIter > 0) start += m_numVerts[pIter - 1];

    int numV = m_numVerts[pIter];

    double x0, y0, x1, y1;
    if (numV <= 0){
      x0 = DBL_MAX/4.0, x1 = -DBL_MAX/4.0; // Use 1/4.0 to avoid overflow when ...
      y0 = DBL_MAX/4.0, y1 = -DBL_MAX/4.0; // ... finding width and height
    }else{
      const double * px = vecPtr(m_xv) + start;
      const double * py = vecPtr(m_yv) + start;
      x0 = *min_element( px, px + numV ); x1 = *max_element( px, px + numV );
      y0 = *min_element( py, py + numV ); y1 = *max_element( py, py + numV );
    }
    xll.push_back(x0); xur.push_back(x1);
    yll.push_back(y0); yur.push_back(y1);
    
  }
  
  return;
};

void dPoly::bdBoxCenter(double & mx, double & my) const{

  double xll, yll, xur, yur;
  bdBox(xll, yll, xur, yur);
  
  mx = (xll + xur)/2.0;
  my = (yll + yur)/2.0;

  return;
}

void dPoly::setPolygon(int numVerts,
                       const double * xv,
                       const double * yv,
                       bool isPolyClosed,
                       const std::string & color
                       ){
  reset();
  appendPolygon(numVerts, xv, yv, isPolyClosed, color);
  return;
}

void dPoly::appendPolygon(int numVerts,
                          const double * xv,
                          const double * yv,
                          bool isPolyClosed,
                          const std::string & color
                          ){

  if (numVerts <= 0) return;
  
  m_numPolys      += 1;
  m_totalNumVerts += numVerts;
  
  m_numVerts.push_back(numVerts);
  m_isPolyClosed.push_back(isPolyClosed);
  m_colors.push_back(color);
  for (int s = 0; s < numVerts; s++){
    m_xv.push_back(xv[s]);
    m_yv.push_back(yv[s]);
  }

  return;
}

void dPoly::appendRectangle(double xl, double yl, double xh, double yh,
                            bool isPolyClosed,
                            const std::string & color
                            ){
  double xv[4], yv[4];
  xv[0] = xl; xv[1] = xh; xv[2] = xh; xv[3] = xl;
  yv[0] = yl; yv[1] = yl; yv[2] = yh; yv[3] = yh;
  appendPolygon(4, xv, yv, isPolyClosed, color);
  return;
}

void dPoly::setRectangle(double xl, double yl, double xh, double yh,
                         bool isPolyClosed,
                         const std::string & color
                         ){
  reset();
  appendRectangle(xl, yl, xh, yh, isPolyClosed, color);
  return;
}

void dPoly::clipPoly(// inputs
                     double clip_xll, double clip_yll,
                     double clip_xur, double clip_yur,
                     dPoly & clippedPoly // output
                     ){

  assert(this != &clippedPoly); // source and destination must be different
  
  clippedPoly.reset();
  clippedPoly.set_isPointCloud(m_isPointCloud);
  
  const double * xv               = get_xv();
  const double * yv               = get_yv();
  const int    * numVerts         = get_numVerts();
  int numPolys                    = get_numPolys();
  const vector<char> isPolyClosed = get_isPolyClosed();
  const vector<string> colors     = get_colors();
  
  int start = 0;
  for (int pIter = 0; pIter < numPolys; pIter++){
      
    if (pIter > 0) start += numVerts[pIter - 1];
      
    int  isClosed = isPolyClosed [pIter];
    string color  = colors       [pIter];
      
    vector<double> cutXv, cutYv;
    vector<int> cutNumVerts;
    cutXv.clear(); cutYv.clear(); cutNumVerts.clear();
    
    if (m_isPointCloud){

      // To cut a point cloud to a box all is needed is to select
      // which points are in the box
      for (int vIter = 0; vIter < numVerts[pIter]; vIter++){
        double x = xv[start + vIter];
        double y = yv[start + vIter];
        if (x >= clip_xll && x <= clip_xur &&
            y >= clip_yll && y <= clip_yur
            ){
          cutXv.push_back(x);
          cutYv.push_back(y);
        }
      }
      cutNumVerts.push_back( cutXv.size() );
      
    }else if (isClosed){

      cutPoly(1, numVerts + pIter, xv + start, yv + start,
              clip_xll, clip_yll, clip_xur, clip_yur, 
              cutXv, cutYv, cutNumVerts // outputs
              );
      
    }else{
      
      cutPolyLine(numVerts[pIter], xv + start, yv + start,
                  clip_xll, clip_yll, clip_xur, clip_yur, 
                  cutXv, cutYv, cutNumVerts // outputs
                  );
      
    }
    
    int cstart = 0;
    for (int cIter = 0; cIter < (int)cutNumVerts.size(); cIter++){
        
      if (cIter > 0) cstart += cutNumVerts[cIter - 1];
      int cSize = cutNumVerts[cIter];
      clippedPoly.appendPolygon(cSize,
                                vecPtr(cutXv) + cstart,
                                vecPtr(cutYv) + cstart,
                                isClosed, color
                                );

    }
    
  }

  return;
} 


void dPoly::appendPolygons(const dPoly & poly){

  const double * xv         = poly.get_xv();
  const double * yv         = poly.get_yv();
  const int    * numVerts   = poly.get_numVerts();
  int numPolys              = poly.get_numPolys();
  vector<char> isPolyClosed = poly.get_isPolyClosed();
  vector<string> colors     = poly.get_colors();
  
  int start = 0;
  for (int pIter = 0; pIter < numPolys; pIter++){
      
    if (pIter > 0) start += numVerts[pIter - 1];
      
    bool isClosed = isPolyClosed [pIter];
    string color  = colors       [pIter];
    int pSize     = numVerts     [pIter];
    
    appendPolygon(pSize, xv + start, yv + start, isClosed, color);
    
  }

  return;
}


void dPoly::set_isPolyClosed(bool isPolyClosed){

  m_isPolyClosed.resize(m_numPolys);
  for (int s = 0; s < (int)m_isPolyClosed.size(); s++){
    m_isPolyClosed[s] = isPolyClosed;
  }

  return;
}

bool dPoly::readPoly(std::string filename,
                     // If isPointCloud is true, treat each point as a
                     // singleton polygon
                     bool isPointCloud 
                     ){

  reset();
  
  m_isPointCloud = isPointCloud;

  // To do: The test below will succeed if filename is a directory.
  // This is probably not right.
  ifstream fh(filename.c_str());
  if( !fh ){
    cerr << "Error: Could not open " << filename << endl;
    return false;
  }

  if (boost::iends_with(filename, ".tif") ||
      boost::iends_with(filename, ".ntf") ||
      boost::iends_with(filename, ".cub")
      ){
    m_img = DiskImageView<double>(filename);
    boost::shared_ptr<DiskImageResource> rsrc(DiskImageResource::open(filename));
    ChannelTypeEnum channel_type = rsrc->channel_type();
    double nodata_val = -32768;
    if ( rsrc->has_nodata_read() ) {
      nodata_val = rsrc->nodata_read();
    }

    // Must scale the image values to uint8
    if(channel_type != VW_CHANNEL_UINT8){
      double mn = DBL_MAX, mx = -DBL_MAX;
      for (int col = 0; col < m_img.cols(); col++){
        for (int row = 0; row < m_img.rows(); row++){
          if (m_img(col, row) <= nodata_val) continue;
          if (m_img(col, row) < mn) mn = m_img(col, row);
          if (m_img(col, row) > mx) mx = m_img(col, row);
        }
      }
      if (mn >= mx){
        for (int col = 0; col < m_img.cols(); col++){
          for (int row = 0; row < m_img.rows(); row++){
            m_img(col, row) = 0.0;
          }
        }
      }else{
        for (int col = 0; col < m_img.cols(); col++){
          for (int row = 0; row < m_img.rows(); row++){
            m_img(col, row) = round(255*(std::max(double(m_img(col, row)), mn)
                                         - mn)/(mx-mn));
          }
        }
      }
    }
    
    BBox2 b = bounding_box(m_img);
    setRectangle(b.min().x(), b.min().y(), b.max().x(), b.max().y(), true,
                 "red");
    return true;
  }
  
  // The current polygon has vertices in the range [beg, end)
  int beg = 0, end = 0;
  
  string color = "yellow"; // default color for polygons
  std::string line;
  while( getline(fh, line) ) {
    
    bool isLastLine = ( fh.peek() == EOF );

    // Convert to lowercase
    transform(line.begin(), line.end(), line.begin(), ::tolower);

    char * linePtr = (char*)line.c_str(); // To do: Avoid this casting hack.

    // Replace comma with space, to be able to use comma as separator
    for (int s = 0; s < (int)strlen(linePtr); s++){
      if (linePtr[s] == ',') linePtr[s] = ' ';
    }

    // Ignore any text after the comment character, which is '#' or '!'
    for (int s = 0; s < (int)strlen(linePtr); s++){
      if (linePtr[s] == '#' || linePtr[s] == '!'){
        for (int t = s; t < (int)strlen(linePtr); t++){
          linePtr[t] = '\0';
        }
        break;
      }
    }
    
    // If this is the last line in the file, or if we encountered a
    // "next" statement, or if we treat a polygon as just a set of
    // points (point cloud) then close the current polygon and start a
    // new one.
    istringstream iss_next(line);
    string val;
    bool isLastVertexOfCurrPoly = ( isLastLine                               ||
                                    ( (iss_next >> val) && (val == "next") ) ||
                                    isPointCloud
                                    );
    bool isCurrPolyNonEmpty = (beg < end);
    if (isLastVertexOfCurrPoly && isCurrPolyNonEmpty){
      
      assert( end == (int)m_xv.size() && end == (int)m_yv.size() );
      
      if (beg < end - 1              &&
          m_xv[beg] == m_xv[end - 1] &&
          m_yv[beg] == m_yv[end - 1]){
        // The first vertex equals to the last vertex in the current
        // polygon. That means that this is a true polygon rather
        // than a polygonal line. Don't store the last
        // vertex.
        end--;
        m_xv.resize(end);
        m_yv.resize(end);
        m_isPolyClosed.push_back(true);
      }else{
        m_isPolyClosed.push_back(false);
      }
      
      m_colors.push_back(color);
      
      m_numPolys++;
      m_numVerts.push_back(end - beg);
      m_totalNumVerts = end;
      
      // Start a new polygon
      beg = end;
      
    } // End processing the current polygon in the list of polygons

  } // End reading the file and processing all polygons
  
  return true; // success
  
}

void utils::bdBox(const std::vector<dPoly> & polyVec,
                  // outputs
                  double & xll, double & yll,
                  double & xur, double & yur
                  ){

  double big = DBL_MAX;
  xll = big; yll = big; xur = -big; yur = -big;

  // Iterate through input polynomial list
  for (int p = 0; p < (int)polyVec.size(); p++){

    // Skip empty objects
    if (polyVec[p].get_totalNumVerts() == 0)
      continue;

    // Determine bounding box of this object
    double xll0, yll0, xur0, yur0;
    polyVec[p].bdBox(xll0, yll0, xur0, yur0);

    // Update group bounding box
    xll = min(xll, xll0); xur = max(xur, xur0);
    yll = min(yll, yll0); yur = max(yur, yur0);
  }
  
  return;
}
  
void utils::setUpViewBox(// inputs
                         const std::vector<dPoly> & polyVec,
                         // outputs
                         double & xll,  double & yll,
                         double & widx, double & widy
                         ){
  
  // Given a set of polygons, set up a box containing these polygons.

  double xur, yur; // local variables
  
  bdBox(polyVec,           // inputs 
        xll, yll, xur, yur // outputs
        );
  
  // Treat the case of empty polygons
  if (xur < xll || yur < yll){
    xll = 0.0; yll = 0.0; xur = 1000.0; yur = 1000.0;
  }

  // Treat the case when the polygons are degenerate
  if (xur == xll){ xll -= 0.5; xur += 0.5; }
  if (yur == yll){ yll -= 0.5; yur += 0.5; }
    
  widx = xur - xll; assert(widx > 0.0);
  widy = yur - yll; assert(widy > 0.0);

  // Expand the box slightly for plotting purposes
  double factor = 0.05;
  xll -= widx*factor; xur += widx*factor; widx *= 1.0 + 2*factor;
  yll -= widy*factor; yur += widy*factor; widy *= 1.0 + 2*factor;
  
  return;
  
}

polyOptions::polyOptions(){
  plotAsPoints     = false;
  isPolyFilled     = false;
  isPolyClosed     = readClosedPolyInfoFromFile;
  fontSize         = 10; 
  lineWidth        = 1;
  useCmdLineColor  = false;
  isGridOn         = false;
  gridWidth        = 1;
  gridSize         = -1;
  readPolyFromDisk = true;
  bgColor          = "black";
  fgColor          = "white";
  cmdLineColor     = "green";
  gridColor        = "white";
  polyFileName     = "unnamed.xg";
}

cmdLineOptions::cmdLineOptions(){
  windowWidX       = 1200;
  windowWidY       = 800;
}

void utils::printUsage(std::string progName){
  cout << "Usage: " << progName << " "
       << "[ -geo[metry] <width>x<height> ] [-bg | -backgroundColor <color> ] "
       << "file1 ... fileN" << endl;
}

void utils::extractWindowDims(// inputs
                              int numArgs, char ** args,
                              // outputs
                              int & windowWidX, int & windowWidY
                              ){

  // Parse the command line arguments '-geo[metry] 500x600'

  for (int s = 1; s < numArgs; s++){

    if ( !strstr(args[s-1], "-geo") ) continue;

    string lineStr = args[s];
    char * line    = (char*) lineStr.c_str();
    
    // Blank the geometry settings once located
    // to not confuse other parsers.
    args[s-1][0] = '\0';
    args[s  ][0] = '\0';
    
    char * pch;
    char delimiter[] = "x";
    
    pch = strtok (line, delimiter);
    if (pch == NULL) continue;
    int windowWidX_tmp = atoi(pch);

    pch = strtok (NULL, delimiter);
    if (pch == NULL) continue;
    int windowWidY_tmp = atoi(pch);
    
    if (windowWidX_tmp > 0 && windowWidY_tmp > 0){
      windowWidX = windowWidX_tmp;
      windowWidY = windowWidY_tmp;
    }
    
  }
  
}

void utils::parseCmdOptions(//inputs
                            int argc, char** argv, std::string exeName,
                            // outputs
                            cmdLineOptions & cmdOpt
                            ){

  cmdOpt.polyOptionsVec.clear();
  
  polyOptions opt; // Each polygon file will have one such entry
  
  // Skip argv[0] as that's the program name
  extractWindowDims(argc - 1, argv + 1, cmdOpt.windowWidX, cmdOpt.windowWidY);

  for (int argIter = 1; argIter < argc; argIter++){

    char * currArg = argv[argIter];

    if (currArg == NULL || strlen(currArg) == 0) continue;

    if (currArg[0] == '-'){
      // Transform -P into -p, etc.
      transform(currArg, currArg + strlen(currArg), currArg, ::tolower);
    }

    if (strcmp( currArg, "-h"     ) == 0 || strcmp( currArg, "--h"    ) == 0 ||
        strcmp( currArg, "-help"  ) == 0 || strcmp( currArg, "--help" ) == 0 ||
        strcmp( currArg, "-?"     ) == 0 || strcmp( currArg, "--?"    ) == 0 ){
      printUsage(exeName);
      exit(0);
    }
    
    if ( (strcmp(currArg, "-bg") == 0 || strcmp(currArg, "-backgroundcolor") == 0 )
         &&
         argIter < argc - 1
         ){
      opt.bgColor = argv[argIter + 1];
      argIter++;
      continue;
    }

    if ( (strcmp(currArg, "-fs"      ) == 0 ||
          strcmp(currArg, "-fontsize") == 0 )
         &&
         argIter < argc - 1
         ){
      int fs = (int)round(atof(argv[argIter + 1]));
      if (fs > 0) opt.fontSize = fs;
      argIter++;
      continue;
    }

    if ( (strcmp(currArg, "-c"    ) == 0 ||
          strcmp(currArg, "-color") == 0 )
         && argIter < argc - 1){
      opt.useCmdLineColor = true;
      opt.cmdLineColor    = argv[argIter + 1];
      argIter++;
      continue;
    }

    // Other command line options are ignored
    if (currArg[0] == '-') continue;
    
    opt.polyFileName = currArg;

    cmdOpt.polyOptionsVec.push_back(opt);
  }

  // Push one more time, to guarantee that the options vector is
  // non-empty even if no polygons were provided as input, and to make
  // sure we also parsed the options after the last polygon filename.
  cmdOpt.polyOptionsVec.push_back(opt);
  
  return;
}

std::string utils::getFilenameExtension(std::string filename){

  std::string::size_type idx;
  idx = filename.rfind('.');

  if(idx != std::string::npos) return filename.substr(idx+1);
  else                         return "";
}

std::string utils::replaceAll(std::string result, 
                              const std::string & replaceWhat, 
                              const std::string & replaceWithWhat){
  
  while(1){
    const int pos = result.find(replaceWhat);
    if (pos == -1) break;
    result.replace(pos,replaceWhat.size(),replaceWithWhat);
  }
  return result;
}

