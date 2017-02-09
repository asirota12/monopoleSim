{
  int nbins=300;
  double xmin=0;
  double xmax=1000;
 // double logxmin=TMath::Log10(xmin);
 // double logxmax=TMath::Log10(xmax);
 // double binwidth=(logxmax-logxmin)/nbins;
 // double xbins[nbins+1];

 // for(int i=0;i<=nbins;i++) xbins[i]=TMath::Power(10,logxmin+i*binwidth);

  TFile f("beta.01.root");
  TCanvas *c=new TCanvas("Histogram","Saturation Count for v/c=.01 monopoles, Graph 4");
  c->SetLogy();
  //c->SetLogx();
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(kFALSE);
  c->Draw();
  string hname[2]={"hist1","hist2"};
  TH1F **hist = new TH1F*[2];
  for(int i=0; i<2; i++)
  {
    hist[i] = new TH1F(hname[i].c_str(),"Histogram at .1 beta",nbins,xmin,xmax);// name,title,bins,xmin,xmax
    hist[i]->GetXaxis()->SetTitle("Saturation Count");
    hist[i]->GetYaxis()->SetTitle("Number of Microslices");
    hist[i]->SetMinimum(.1);
    hist[i]->SetMaximum(1e4);
    hist[i]->SetLineColor(i+2);
    Ntuple3->Draw(("satCount>>+"+hname[i]).c_str(),i==0?"hasMono==0":"hasMono==1","same");
  }
}
