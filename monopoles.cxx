#include "detector.h"
#include <TStyle.h>
#include <TH2F.h>
#include <iostream>
#include <string>

const double twopi=TMath::TwoPi();
const double piover2=TMath::PiOver2();
using namespace std;

extern "C" {
  void hrndg2_( double*,
                const long*   const,
                const double* const,
                const double* const,
                const long*   const,
                const double* const,
                const double* const,
                double*,
                double*,
                double*,
                const float* const);
}
class Hit{
  Cell              *_cell;
  CLHEP::Hep3Vector  _entrancePoint, _exitPoint;
  double             _entranceTime, _exitTime;
  double             _depositedEnergy;
  
public:
  Hit(Cell *cell, CLHEP::Hep3Vector entrancePoint, CLHEP::Hep3Vector exitPoint,
      double entranceTime, double exitTime, double depositedEnergy) : 
                    _cell(cell), _entrancePoint(entrancePoint), _exitPoint(exitPoint),
                    _entranceTime(entranceTime), _exitTime(exitTime),
                    _depositedEnergy(depositedEnergy){}
  
  const Cell *getCell() const                       {return _cell;}
  const CLHEP::Hep3Vector &getEntrancePoint() const {return _entrancePoint;}
  const CLHEP::Hep3Vector &getExitPoint() const     {return _exitPoint;}
  double getEntranceTime() const                    {return _entranceTime;}
  double getDepositedEnergy() const                 {return _depositedEnergy;}
  double getExitTime() const                        {return _exitTime;}
  bool getSaturation() const {
    if(_depositedEnergy>=273){ //1 arbtitrary
      return true;
    } else {
      return false;
    }
  }
};

class Event
{
  Trajectory        _trajectory;
  std::vector<Hit*> _hits;

  public:
  Event(Trajectory trajectory) : _trajectory(trajectory)  { }

  std::vector<Hit*> getHits(){
    return _hits;
  }

  Trajectory getTrajectory(){
    return _trajectory;
  }

  void checkTrajectory(const Detector &detector, TRandom3 &r)
  {
    const std::vector<Cell*> &cells = detector.getCells();
    std::vector<Cell*>::const_iterator iter;
    for(iter=cells.begin(); iter!=cells.end(); iter++)
    {
      CLHEP::Hep3Vector entrancePoint, exitPoint;
      double entranceTime, exitTime;
      if((*iter)->checkTrajectory(_trajectory, entrancePoint, exitPoint, entranceTime, exitTime))
      {
        if(exitTime<(_trajectory.cutoff())){
	  double dE_dx=_trajectory.getdE_dx();
	  double isMonopole=_trajectory.isMonopole();
	  double depositedEnergy=(*iter)->depositedEnergy(entrancePoint,exitPoint,dE_dx,r);
	  Hit *hit=new Hit(*iter, entrancePoint, exitPoint, entranceTime, exitTime, depositedEnergy);
	  _hits.push_back(hit);
	}
      }
    }
  }

  void Draw(TCanvas *canvas)
  {
    _trajectory.Draw(canvas);
    std::vector<Hit*>::const_iterator iter;
    for(iter=_hits.begin(); iter!=_hits.end(); iter++)
    {
      bool isVerticalCell=(*iter)->getCell()->isVertical();
      double time=(*iter)->getEntranceTime();
      time/=1.0e-7;
      if(time<0) time=0;
      if(time>1) time=1;
      double saturation=(*iter)->getSaturation();
      (*iter)->getCell()->Draw(canvas, time, saturation);
    }
  }

  void Analyze(TNtuple *ntuple,TNtuple *ntuple2)
  {
    CLHEP::Hep3Vector pos1, pos2,start1,stop1;
    double time1=NAN;
    double time2=NAN;
    double TEnergy,b,c,d;
    bool triggered=false;
    bool missed=false;
    TEnergy=0;
    std::vector<Hit*>::const_iterator iter;
    for(iter=_hits.begin(); iter!=_hits.end(); iter++){
      double hitdepo=(*iter)->getDepositedEnergy();
      bool isMono=_trajectory.isMonopole();
      double tsublength;
      start1=(*iter)->getEntrancePoint();
      stop1=(*iter)->getExitPoint();
      tsublength=(start1-stop1).mag();
      double t1=(*iter)->getEntranceTime();
      double t2=(*iter)->getExitTime();
      
      ntuple2->Fill((hitdepo/tsublength),isMono,tsublength,tsublength/(t2-t1)/(2.998e11));
      TEnergy=TEnergy+((*iter)->getDepositedEnergy());
      if(t1<time1 || isnan(time1))
      {
        time1=t1;
        pos1=(*iter)->getEntrancePoint();
      }
      if(t2>time2 || isnan(time2))
      {
        time2=t2;
        pos2=(*iter)->getExitPoint();
      }
    }
    double muonenergy=0;
    if (!_trajectory.isMonopole()){
      muonenergy=_trajectory.getenergy();
    }
    double tracklength=(pos2-pos1).mag();
    if((TEnergy/tracklength)>1500){
      triggered=true;
    } 
    if((time1-time2/tracklength)<2.6e11){
      triggered=true;
    }
    if (_trajectory.isMonopole()&&!triggered&&(_hits.size()!=0)){
      missed=true;
    }
    const CLHEP::Hep3Vector dir= _trajectory.getDirection();
    double ctry=acos(-dir.x()/sin(acos(-1*dir.y())));       
    double stry=asin(-dir.z()/sin(acos(-1*dir.y())));
    double angle=0;
    if(dir.x()<0 && dir.y()<0 && dir.z()>0){
      angle=stry+twopi;//std::cout<<"phi quad 4 "<< fabs(phi-(stry+twopi)) << std::endl;
    } 
    else if(dir.x()<0 && dir.y()<0 && dir.z()<0){
      angle=ctry;//std::cout<<"phi quad 1 "<<fabs(phi-(ctry))<<std::endl;
    }
    else if(dir.x()>0 && dir.y()<0){ 
      angle=-stry+twopi/2;//std::cout<<"phi quad 2/3 "<<fabs(phi-(-stry+twopi/2))<<std::endl;
    }   
      
    ntuple->Fill(tracklength,_hits.size(),TEnergy,triggered,missed,muonenergy,_trajectory.getdE_dx(),_trajectory.isMonopole(),-1*dir.y(),angle);
  }
};

double de(double a){
  double e=log10(a);
  return (.2384*pow(e,6)-1.5049*pow(e,5)+2.5366*pow(e,4)+1.0809*pow(e,3)-4.0038*pow(e,2)+.4245*e+2.3786)*.8;
 // ok so this is now in de_dx in MeV/cm. The .8 is scintillator density
}

double dE_dxmono(double beta){
  double x=log10(beta);
  if (x<=-2){
    return pow(10.0,4.57+.99*x);// Mev/cm
  }
  if (-2<x && x<=-1){
    return pow(10,5.45+2.38*x+.478*x*x);
  }
  if (-1<x){ 
    return pow(10,3.79+.15*x+.013*x*x+.336*x*x*x+.1959*x*x*x*x); 
  } // taken from Guilliame
}

void simulateEvents(int numSlices,double cutoff,double beta,string name, TCanvas *canvas=NULL){
  std::vector<Event*>   Slices[numSlices*2]; //I believe i could use an array of slices. that is an array of vectors of event*s.
  std::time_t t1=0;
  std::time_t t2=0;
  std::time_t t3=std::time(0);
  std::time_t t4=0;
  Detector detector;
  TRandom3 r;
  int time=0;
  std::vector<Event*> events;
  double y=20000; //millimeters
  double xMin=-5000;
  double xMax=20000;
  double zMin=-5000;
  double zMax=65000;
  double tMin=0;
  double tMax=1e-7;
  double v=0;
  double c=2.998e11;
  double dE_dx=0; 
  double dE_dxmon=0;
  long int ii=1e6;//number of energy bins
  long int jj=100;//number of costheta bins 
    //these two can make this very innacurate and wavy distribution.
  std::vector<double> workingspace;
  workingspace.resize(ii*jj);
  double a1=.1; //energy min
  double b1=10000; //energy max
  double a2=0.00366518; //costh min
  double b2=1; // costh max
  double energy;
  double costh;
  float pro=1;//making it initialize
  double sum; //total integral. 
  hrndg2_(&workingspace[0],&ii,&a1,&b1,&jj,&a2,&b2,&sum,&energy,&costh,&pro);
  pro=111; // making it actually run
  double rate=sum*twopi/2.*.08*(xMax-xMin)/2.*(zMax-zMin)/2.;//its in millimeters
  //std::cout<<"Rate :"<<rate<<" per slice "<<rate*cutoff<<std::endl;
  int numMuons;//rate
 
  bool wantMono=false;
  for (int again=0;again<2;again++){
    if (again==1){
      wantMono=true;
    }
    for (int j=0;j<numSlices;j++){
      numMuons=r.Poisson(rate*cutoff);
      std::vector<Event*> sliceEvents;
      t4=std::time(0)-t3;
      std::cout<<(j+numSlices*again+1)/(double)numSlices*50<<"% completed "<<(t4)/60.<<" minutes"<<std::endl;
      for(double i=-1; i<numMuons; i++){
	
	double x=r.Rndm()*(xMax-xMin)+xMin;
	double z=r.Rndm()*(zMax-zMin)+zMin;
	double t=r.Rndm()*(tMax-tMin)+tMin;
	double phi=r.Rndm()*twopi;
        	
	hrndg2_(&workingspace[0],&ii,&a1,&b1,&jj,&a2,&b2,&sum,&energy,&costh,&pro);
	
	double v=pow(1-105.7*105.7/(energy*energy*1e6),1.0/2.0)*c; // getting from E to v
	double theta=acos(costh); 
	double p=pow(energy*energy-.1057*.1057,1./2.); 
	double dE_dx=0; 
	dE_dx=de(energy);//input in  GeV/c, it should be momentum 
        
	double xDir=-sin(theta)*cos(phi);
	double zDir=-sin(theta)*sin(phi);
	double yDir=-cos(theta);

	CLHEP::Hep3Vector start(x,y,z);
	CLHEP::Hep3Vector direction(xDir,yDir,zDir);

        double angle=acos(-direction.x()/sin(acos(-1*direction.y())));
	//make muon traj.

	Trajectory trajectory(start,direction,v,t,dE_dx,energy,false,cutoff);
	Event *event = new Event(trajectory);
	event->checkTrajectory(detector, r);


	std::vector<Hit*> count5=(event)->getHits();
	if(count5.size()>0){
	  events.push_back(event);
	  sliceEvents.push_back(event);	  
	}                                               // this is no muon 
	if(canvas) event->Draw(canvas);
	
	//make Monopole traj.
	if((i==numMuons-1)&&(wantMono)){
	  bool keeptrying=true;
	  int y=0;
	  Event *event2;
	  while (keeptrying){
	    x=r.Rndm()*(xMax-xMin)+xMin;
	    z=r.Rndm()*(zMax-zMin)+zMin;
	    t=r.Rndm()*(tMax-tMin)+tMin;
	    phi=r.Rndm()*twopi;
	    theta=acos(r.Rndm());
	    
	    xDir=-sin(theta)*cos(phi);
	    zDir=-sin(theta)*sin(phi);
	    yDir=-cos(theta);

	    CLHEP::Hep3Vector start2(x,y,z);
	    CLHEP::Hep3Vector direction2(xDir,yDir,zDir);
	    double rrr=dE_dxmono(beta);
	    Trajectory trajectory(start2,direction2,beta*c,t,rrr,0,true,cutoff);
	    event2 = new Event(trajectory);  
	    event2->checkTrajectory(detector, r);
	    
	    std::vector<Hit*> count=(event2)->getHits();
	    if(count.size()>1){
    	      std::vector<Hit*>::const_iterator iter0;
	      const Cell* first=(*count.begin())->getCell();
	      int cdiblock0=first->getDiblock();
	      int clayer0=first->getLayer();
	      int cmod0=first->getModule();
              for(iter0=count.begin(); iter0!=count.end(); iter0++){
	        const Cell* one=(*iter0)->getCell(); 
	        int clayer=one->getLayer();
	        int cmod=one->getModule();
	        int cdiblock=one->getDiblock();
		if(cmod0!=cmod && (cdiblock0!=cdiblock || clayer!=clayer0)){
	        keeptrying=false;
		}
	      }
	    } 
	    if(canvas) event2->Draw(canvas);	      
	  }
	  events.push_back(event2);
	  sliceEvents.push_back(event2);
	}
      }
      Slices[j+again*numSlices]=sliceEvents;
    }
  }
  TFile *file=new TFile(name.c_str(),"RECREATE");
  TNtuple *ntuple = new TNtuple("Ntuple","Ntuple","tracklength:numberhits:TdepoEnergy:triggered:missed:muonenergy:dE_dx:isMonopole:cosTh:phi");
  TNtuple *ntuple2 = new TNtuple("Ntuple2","Ntuple2","EnDep:isMonopole:tsublength:beta");
  TNtuple *ntuple3= new TNtuple("Ntuple3","Ntuple3","EnDep:hasMono:satCount:numMuons");
  std::vector<Event*>::const_iterator iter;

  t1=std::time(0);
  for(iter=events.begin(); iter!=events.end(); iter++){
    (*iter)->Analyze(ntuple,ntuple2);
  }

  double SliceDepoEnergy[numSlices*2];
  bool hasMono[numSlices*2];
  int SaturationCount[numSlices*2];
  double tempEn;
  for(int i=0;i<numSlices*2;i++){
    SliceDepoEnergy[i]=0;
    hasMono[i]=false;
    SaturationCount[i]=0;
    std::vector<Event*>::const_iterator iter2;
    int n=0;
    for(iter2=Slices[i].begin();iter2!=Slices[i].end();iter2++){
      // ok so now iter2 is an Event
      n++;
      std::vector<Hit*>::const_iterator iter3;
      std::vector<Hit*> SliceHits=(*iter2)->getHits();
      
      for(iter3=SliceHits.begin(); iter3!=SliceHits.end(); iter3++){
	// now iter3 is a hit
	tempEn=(*iter3)->getDepositedEnergy();
	if((*iter3)->getSaturation()){
          SaturationCount[i]++;
        }
       	SliceDepoEnergy[i]=SliceDepoEnergy[i]+tempEn;
      }
      Trajectory temp=(*iter2)->getTrajectory();
      if(temp.isMonopole()){
	hasMono[i]=true;
        n--;
      }
    }
    ntuple3->Fill(SliceDepoEnergy[i],hasMono[i],SaturationCount[i],n);
  }
  
  // each Event has access to everything i need. These events are in slices.

  t2=std::time(0);
  //std::cout<<"Seconds to Analyze: "<<t2-t1<<std::endl;
  ntuple->Write();
  ntuple2->Write();
  ntuple3->Write();
  file->Close(); 
}


int main(int argc, char *argv[]){
  std::time_t t1=std::time(0);
  if(argc<5 || argc>6){
    std::cout<<"Wrong number of parameters"<<std::endl;
    std::cout<<"(1) number of slices"<<std::endl;
    std::cout<<"(2) cutoff"<<std::endl;
    std::cout<<"(3) monopole beta"<<std::endl;
    std::cout<<"(4) file store name"<<std::endl;
    std::cout<<"(5) optional   draw"<<std::endl;
    return(0);
  }
  
  int numberSlices=atoi(argv[1]);
  double cutoff=atof(argv[2]);
  double beta=atof(argv[3]);
  string name=argv[4];
  
  if(argc==6){
    if(strcmp(argv[4],"draw")==0){
      TApplication *app = new TApplication("App",0,0);
      gStyle->SetOptTitle(kFALSE);
      gStyle->SetOptStat(kFALSE);
      TCanvas *canvas = new TCanvas("Canvas", "Canvas",1200,800);
      canvas->Divide(1,2);//
      canvas->cd(1);
      TH2F *h1=new TH2F("z-x","",10,-5000,65000,10,-5000,20000);
      h1->SetYTitle("x [mm]");
      h1->SetXTitle("z [mm]");
      h1->Draw();
      canvas->cd(2);
      TH2F *h2=new TH2F("z-y","",10,-5000,65000,10,-5000,20000);
      h2->SetYTitle("y [mm]");
      h2->SetXTitle("z [mm]");
      h2->Draw();
      Detector::setColors();
      simulateEvents(numberSlices,cutoff,beta,name,canvas);
      app->Run();
      return 0;
    }
  }

  simulateEvents(numberSlices,cutoff,beta,name);
  std::time_t t2=std::time(0);
  //std::cout<<"Time for program to run: "<<t2-t1<<std::endl;
  return(0);
}
