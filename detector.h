#include "CLHEP/Vector/ThreeVector.h"
#include <TPolyLine3D.h>
#include <TApplication.h>
#include <TPad.h>
#include <TView.h>
#include <TView3D.h>
#include <TCanvas.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TGeoSphere.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TRandom3.h>
#include <TAxis3D.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <math.h>

const double cellWidthXY=39.74;
const double cellWidthZ=66.66; 
const double cellLengthH=15443.00;
const double cellLengthV=15443.00;
const double gapXY=0.01;
const double gapZ=0.01;
const double xOffset=170.00;
const double yOffset=170.00;

class Trajectory
{
  CLHEP::Hep3Vector _start, _direction;
  double            _velocity;
  double            _startTime;
  double            _dE_dx;
  double            _energy;
  bool              _isMonopole;
  double            _cutoff;

  public:
  Trajectory(const CLHEP::Hep3Vector &start, const CLHEP::Hep3Vector &direction, 
             double velocity, double startTime,double dE_dx,double energy, bool isMonopole,double cutoff);

  const CLHEP::Hep3Vector &getStart() const     {return _start;}
  const CLHEP::Hep3Vector &getDirection() const {return _direction;}
  double getVelocity() const                    {return _velocity;}
  double getStartTime() const                   {return _startTime;}
  double getdE_dx() const                       {return _dE_dx;}
  double getenergy() const                      {return _energy;}
  bool   isMonopole() const                     {return _isMonopole;}
  double cutoff() const                         {return _cutoff;}

  void Draw(TCanvas *canvas) const;
};

class Cell
{
  CLHEP::Hep3Vector _midpoint;
  CLHEP::Hep3Vector _halflengths;
  CLHEP::Hep3Vector _corner1, _corner2; //perhaps use arrays for shorter access time, 
                                        //can't use default copy c'tor in this case
  int               _id;
  int               _cellnumber, _module, _layer, _diblock;
  bool              _isVertical;
  mutable bool      _drawn;

  public:
  Cell(int id, int cellnumber, int module, int layer, int diblock, bool isVertical);
  void setGeometry(const CLHEP::Hep3Vector &midpoint, const CLHEP::Hep3Vector &halflengths);

  int  getID() const         {return _id;} 
  int  getCellNumber() const {return _cellnumber;} 
  int  getModule() const     {return _module;} 
  int  getLayer() const      {return _layer;} 
  int  getDiblock() const    {return _diblock;} 
  bool isVertical() const    {return _isVertical;} 
  CLHEP::Hep3Vector getMidPoint() const {return _midpoint;}

  bool checkTrajectory(const Trajectory &trajectory, CLHEP::Hep3Vector &entrancePoint, CLHEP::Hep3Vector &exitPoint,
                       double &entranceTime, double &exitTime);

  void Draw(TCanvas *canvas, double time, double saturation) const;

  static double depositedEnergy(const CLHEP::Hep3Vector &entrancePoint, const CLHEP::Hep3Vector &exitPoint, 
                                double dE_dx, TRandom3 &r);
};

class HorizontalCell : public Cell
{
  public:
  HorizontalCell(int id, int cellnumber, int module, int layer, int diblock);
};

class VerticalCell : public Cell
{
  public:
  VerticalCell(int id, int cellnumber, int module, int layer, int diblock);
};

class Detector
{
  std::vector<Cell*> _cells;

  public:
  Detector();

  const std::vector<Cell*> &getCells() const;
  const Cell* getCell(int id) const;
  static void setColors();
};
