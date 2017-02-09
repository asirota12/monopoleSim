#include "detector.h"

#include "CLHEP/Vector/ThreeVector.h"
#include <TPolyLine.h>
#include <TBox.h>
#include <TROOT.h>
#include <TColor.h>
#include <TApplication.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <math.h>

  Trajectory::Trajectory(const CLHEP::Hep3Vector &start, const CLHEP::Hep3Vector &direction, 
			 double velocity, double startTime,double dE_dx,double energy, bool isMonopole,double cutoff) :
    _start(start), _direction(direction), _velocity(velocity), _startTime(startTime),
    _dE_dx(dE_dx) ,_energy(energy), _isMonopole(isMonopole),_cutoff(cutoff)
  {
  }

  void Trajectory::Draw(TCanvas *canvas) const
  {
    for(int i=0; i<2; i++)
    {
      canvas->cd(i+1);
      TPolyLine *line=new TPolyLine();
      line->SetLineWidth(1);
      line->SetLineColor(_isMonopole?4:46);
      line->SetPoint(0,_start[2],_start[i]);
      CLHEP::Hep3Vector end=_start+_direction*1000000.0;
      line->SetPoint(1,end[2], end[i]);
      line->Draw("same");
    }
  }


  Cell::Cell(int id, int cellnumber, int module, int layer, int diblock, bool isVertical) : 
             _id(id), _cellnumber(cellnumber), _module(module), _layer(layer), _diblock(diblock), 
             _isVertical(isVertical), _drawn(false)
  {
  }

  void Cell::setGeometry(const CLHEP::Hep3Vector &midpoint, const CLHEP::Hep3Vector &halflengths)
  {
    _midpoint=midpoint;
    _halflengths=halflengths;
    for(int i=0; i<3; i++)
    {
      _corner1[i]=_midpoint[i]-halflengths[i];
      _corner2[i]=_midpoint[i]+halflengths[i];
    }
  }

  bool Cell::checkTrajectory(const Trajectory &trajectory, 
                             CLHEP::Hep3Vector &entrancePoint, CLHEP::Hep3Vector &exitPoint,
                             double &entranceTime, double &exitTime)
  {
    double entranceTimeSide[3], exitTimeSide[3]; //these are the entrance/exit times of the trajectory 
                                                 //when crossing one of the sides of the cell
                                                 //(assuming the start time of the track is 0)
    for(int i=0; i<3; i++) //check all directions
    {
      double velocity=trajectory.getDirection()[i]*trajectory.getVelocity();
      double start=trajectory.getStart()[i];
      if(velocity!=0)
      {
        //entrance and exit in "i" direction
        entranceTimeSide[i]=(_corner1[i]-start)/velocity;  // <--- from  entranceTimeSide*velocity+start=corner1
        exitTimeSide[i]    =(_corner2[i]-start)/velocity;  // <--- from  exitTimeSide*velocity+start=corner2

        //keep all times in ascending order
        if(exitTimeSide[i]<entranceTimeSide[i]) std::swap(entranceTimeSide[i],exitTimeSide[i]); 
      }
      else //trajectories parallel to the cell borderes in "i direction",
           //which means that it will never cross these borders
      {
        if(start<_corner1[i] || start>_corner2[i]) return(false);  //outside of the borders in the "i-dimension"
                                                                   //--> the trajectory doesn't go through this cell
        //otherwise, the the track goes parallel between the two sides - inside of the cell
        //entrance time needs to be determined by other components
        entranceTimeSide[i]=-INFINITY;
        exitTimeSide[i]=INFINITY;
      }
    }

    //find overlaps of the passage times for all three dimensions
    //if there is an overlap, the trajectory went through the cell
    //and the entrance and exit time can be determined
    //example:
    //T1:     +++++++++++
    //T2:          +++++++++++++
    //T3:       +++++++
    //overlap:     ====
    //----------------------------------------------->t
    double entranceTimeOverlap=*std::max_element(entranceTimeSide,entranceTimeSide+3);
    double exitTimeOverlap=*std::min_element(exitTimeSide,exitTimeSide+3);

    if(entranceTimeOverlap<exitTimeOverlap) //trajectory passed through counter
    {
      entrancePoint=trajectory.getStart()+trajectory.getDirection()*trajectory.getVelocity()*entranceTimeOverlap;
      exitPoint=trajectory.getStart()+trajectory.getDirection()*trajectory.getVelocity()*exitTimeOverlap;
      entranceTime=entranceTimeOverlap+trajectory.getStartTime();
      exitTime=exitTimeOverlap+trajectory.getStartTime();
      return(true);
    }
    return(false);
  }

  void Cell::Draw(TCanvas *canvas, double time, double saturation) const
  {
    if(!_drawn)
    {
      int i=0;
      if(!_isVertical) i=1;

      canvas->cd(i+1);
      TBox *box = new TBox(_corner1[2],_corner1[i],_corner2[2],_corner2[i]);
      int timeint=static_cast<int>(time*100.0);
      int saturationint=static_cast<int>(saturation*10.0);
      if(timeint==100) timeint=99; 
      if(saturationint==10) saturationint=9; 
      int color=2000+timeint*10+saturationint;
      box->SetFillColor(color);
      box->SetLineWidth(0);
      box->Draw("same");
      _drawn=true;
    }
  }

  double Cell::depositedEnergy(const CLHEP::Hep3Vector &entrancePoint, const CLHEP::Hep3Vector &exitPoint, 
                                double dE_dx, TRandom3 &r)
  {
    CLHEP::Hep3Vector diff=exitPoint-entrancePoint;
    double trackLength=diff.mag();
    //double error=r.Rndm()*.2+.9; //since the photodiodes have 10:1 noise.
    double ans=dE_dx*trackLength/10; //tracklength is in mm
    if (ans>273){
      ans=273;
    }  
    return ans; //took out error for now
  }

  HorizontalCell::HorizontalCell(int id, int cellnumber, int module, int layer, int diblock) : 
          Cell(id, cellnumber, module, layer, diblock, false)
  {
    double x=cellLengthH/2.0;
    double y=((cellnumber-1)+32*(module-1))*(cellWidthXY+gapXY) + cellWidthXY/2.0;
    double z=2.0*((layer-1)+32*(diblock-1))*(cellWidthZ+gapZ) + cellWidthZ/2.0;
    y+=yOffset;
    CLHEP::Hep3Vector midpoint(x,y,z);
    CLHEP::Hep3Vector halflengths(cellLengthH/2.0, cellWidthXY/2.0, cellWidthZ/2.0);
    setGeometry(midpoint, halflengths);
  }

  VerticalCell::VerticalCell(int id, int cellnumber, int module, int layer, int diblock) : 
          Cell(id, cellnumber, module, layer, diblock, true)
  {
    double x=((cellnumber-1)+32*(module-1))*(cellWidthXY+gapXY) + cellWidthXY/2.0;
    double y=cellLengthV/2.0;
    double z=2.0*((layer-1)+32*(diblock-1))*(cellWidthZ+gapZ) + cellWidthZ/2.0;
    x+=xOffset;
    z+=cellWidthZ+gapZ;
    CLHEP::Hep3Vector midpoint(x,y,z);
    CLHEP::Hep3Vector halflengths(cellWidthXY/2.0 ,cellLengthV/2.0, cellWidthZ/2.0);
    setGeometry(midpoint, halflengths);
  }

  Detector::Detector()
  {
    int id=0;
    for(int diblock=1; diblock<=14; diblock++)
    for(int layer=1; layer<=32; layer++)
    {
      for(int module=1; module<=12; module++)
      for(int cellnumber=1; cellnumber<=32; cellnumber++)
      {
        HorizontalCell *cell = new HorizontalCell(id, cellnumber, module, layer, diblock);
        _cells.push_back(cell);
        id++;
      }
      for(int module=1; module<=12; module++)
      for(int cellnumber=1; cellnumber<=32; cellnumber++)
      {
        VerticalCell *cell = new VerticalCell(id, cellnumber, module, layer, diblock);
        _cells.push_back(cell);
        id++;
      }
    }
  }

  const std::vector<Cell*> &Detector::getCells() const
  {
    return _cells;
  }

  const Cell* Detector::getCell(int id) const
  {
    if(id>=0 && id<_cells.size()) return _cells[id];
    return(NULL);
  }

  void Detector::setColors()
  {
    for(int i=0; i<100; i++)
    for(int j=0; j<10; j++)
    {
      float r,g,b;
//      TColor::HLS2RGB(360.0/100.0*i, 0.1*(j+1), 1.0, r,g,b);
      TColor::HLS2RGB(360.0/100.0*i, 0.5, 1.0, r,g,b);
      if(!gROOT->GetColor(2000+10*i+j)) new TColor(2000+10*i+j, r,g,b);
    }
  }
