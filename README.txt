In set up we have a monopoles.cpp file, which simply put is a simulation
 of magnetic monopoles going through the detector. 
This creates a .root file, 
which you can access with X11 if you put root "saidfile.root" into the cmd line.


 

Then  TFile *_file0 = TFile::Open("beta.1.root") for example will access
beta.1.root, 
and you can do something like 
Ntuple-> Draw("pow(cosTh,3.5)","isMonopole==0") 
which makes a graph pulled from Ntuple 
in that file.
 
As you can see there are some additional programs pulled from the DAYABAY project, 

like hnrdg and gaisser which gives proper statistics on cosmic rays to simulate.



The structure has a detector which holds all the classes to do with the
physical detector such as the cells and their positions.

NOTE it takes a couple of days to run depending on your inputs