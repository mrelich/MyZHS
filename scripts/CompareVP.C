
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// This is just a quick comparison script for two separate  //
// shower parameters to see what the vector potential looks //
// like.  Trying to understand why 100k 40 MeV showers has  //
// a larger vector potential while 40 100 GeV showers even  //
// though they have equivalent energy.                      // 
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "myHist.C"

TString savedir = "plots/";

//----------------------------------------//
// Main
//----------------------------------------//
void CompareVP()
{

  // Specify the the files
  const int nFiles = 2;
  TFile* files[nFiles];
  files[0] = new TFile("rootfiles/beam40MeV_100000Prim.root");
  files[1] = new TFile("rootfiles/beam100e3MeV_40Prim.root");
  
  // Specify the plot names for legend
  TString fnames[nFiles];
  fnames[0] = "E_{e-}^{prim}=40 MeV, x10^{4}";
  fnames[1] = "E_{e-}^{prim}=100 GeV, x40";

  // Colors
  int colors[] = {kBlue, kRed};
  int markers[] = {20, 25};

  // Specify the hist name and the labels
  TString pname  = "VP_avg";
  TString xtitle = "time [ns]";
  TString ytitle = "|RA(#theta_{C},t)| [Vs]";

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  // Make legend
  TLegend* leg = makeLegend(0.15,0.3,0.8,0.9);
  
  // Specify the minimum and maximum for x-range
  float xmin = -1; // ns
  float xmax = 1;  // ns

  // Loop over and plot profiles
  TProfile* profs[nFiles];
  float maximum = -999;
  for(int f=0; f<nFiles; ++f){
    profs[f] = getProfile(files[f],pname,xtitle,ytitle,colors[f],markers[f]);
    leg->AddEntry(profs[f], fnames[f].Data(), "lep");
    if( maximum < profs[f]->GetMaximum() )
      maximum = profs[f]->GetMaximum();
    profs[f]->GetXaxis()->SetRange( profs[f]->GetXaxis()->FindBin(xmin),
				    profs[f]->GetXaxis()->FindBin(xmax));
  }

  // Set maximum and draw
  profs[0]->SetMaximum(maximum*10);
  profs[0]->Draw();
  for(int f=1; f<nFiles; ++f)
    profs[f]->Draw("same");
  leg->Draw("same");
  
  c->SaveAs((savedir+"EnergyNParticleComp.png").Data());
  
  
}
