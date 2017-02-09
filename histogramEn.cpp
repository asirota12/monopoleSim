{
  int nbins=300;
  double xmin=1e3;
  double xmax=3e5;
  double logxmin=TMath::Log10(xmin);
  double logxmax =TMath::Log10(xmax);
  double binwidth=(logxmax-logxmin)/nbins;
  double xbins[nbins+1];

  for(int i=0;i<=nbins;i++) xbins[i]=TMath::Power(10,logxmin+i*binwidth);

  TFile f("beta.1.root");
  TCanvas *c=new TCanvas("Histogram","Total Energy Deposition for v/c=.1 monopoles, Graph 2");
  //c->SetLogy();
 // c->SetLogx();
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(kFALSE);
  c->Draw();
  string hname[2]={"hist1","hist2"};
  TH1F **hist = new TH1F*[2];
  for(int i=0; i<2; i++)
  {
   // hist[i] = new TH1F(hname[i].c_str(),"Histogram at .1 beta",nbins,xbins);// name,title,bins,xmin,xmax
    hist[i]= new TH1F(hname[i].c_str(),"Histogram at .1 beta",nbins, xmin,xmax);
    hist[i]->GetXaxis()->SetTitle("Energy Deposition (MeV)");
    hist[i]->GetYaxis()->SetTitle("Number of Microslices");
   // hist[i]->SetMinimum(.1);
    hist[i]->SetMaximum(1e2);
    hist[i]->SetLineColor(i+2);
    Ntuple3->Draw(("EnDep>>+"+hname[i]).c_str(),i==0?"hasMono==0 && EnDep>0":"hasMono==1 && EnDep>0","same");
  }
}
