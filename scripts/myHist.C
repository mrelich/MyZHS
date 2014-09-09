
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TFile.h"
#include "TPad.h"
#include "TLine.h"

#include <vector>

//------------------------------------------------------------//
// Make line
//------------------------------------------------------------//
TLine* makeLine(float x0, float x1, float y0, float y1, 
		int color=kBlack, int style=1)
{

  TLine* line = new TLine(x0,y0,x1,y1);
  line->SetLineWidth(1);
  line->SetLineColor(color);
  line->SetLineStyle(style);

  return line;
}

//------------------------------------------------------------//
// Make pads, needed for ratio plotting
//------------------------------------------------------------//
void makePads(TCanvas* &c, TPad* &top, TPad* &bot, float ysep=0.3)
{

  c->cd();
  top = new TPad("top","",0,ysep,1,1);
  top->SetTopMargin(0.05);
  top->SetBottomMargin(0.02);
  top->SetRightMargin(0.05);
  top->SetLeftMargin(0.12);
  top->SetNumber(1);

  bot =new TPad("bot","",0,0,1,ysep);
  bot->SetTopMargin(0.015);
  bot->SetBottomMargin(0.3);
  bot->SetRightMargin(0.05);
  bot->SetLeftMargin(0.12);
  bot->SetNumber(2);

}

//------------------------------------------------------------//
// Make TLatex
//------------------------------------------------------------//
TLatex* makeLatex()
{

  TLatex* lat = new TLatex();
  lat->SetTextSize(0.04);
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextColor(kBlack);
  return lat;
}

//------------------------------------------------------------//
// Make legend
//------------------------------------------------------------//
TLegend* makeLegend(float x0, float x1, float y0, float y1)
{

  TLegend* leg = new TLegend(x0, y0, x1, y1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  return leg;

}

//------------------------------------------------------------//
// Make canvas
//------------------------------------------------------------//
TCanvas* makeCanvas(TString name)
{

  TCanvas* c = new TCanvas(name.Data(), name.Data(), 600, 500);
  c->SetTicks(1,1);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.12);
  
  return c;
}

//------------------------------------------------------------//
// Get TProfile
//------------------------------------------------------------//
TProfile* getProfile(TFile* file, TString pname, TString xtitle,
		     TString ytitle, int color, int marker)
{
  
  TProfile* prof = (TProfile*) (file->Get(pname.Data())->Clone(Form("%s_%i",pname.Data(),color)));
  prof->GetXaxis()->SetTitle(xtitle.Data());
  prof->GetYaxis()->SetTitle(ytitle.Data());
  prof->SetMarkerStyle(marker);
  prof->SetMarkerColor(color);
  prof->SetMarkerSize(0.5);
  prof->SetLineColor(color);
  prof->SetTitle("");
  prof->SetStats(0);
  prof->GetYaxis()->SetTitleOffset(1.5);
  
  return prof;

}

//------------------------------------------------------------//
// Get histogram
//------------------------------------------------------------//
TH1F* getHist(TFile* file, TString pname, TString xtitle,
		  TString ytitle, int color, int marker)
{
  
  TH1F* hist = (TH1F*) (file->Get(pname.Data())->Clone(Form("%s_%i",pname.Data(),color)));
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitle.Data());
  hist->SetMarkerStyle(marker);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize(0.5);
  hist->SetLineColor(color);
  hist->SetTitle("");
  hist->SetStats(0);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->SetLineWidth(2);
  return hist;

}

//------------------------------------------------------------//
// Get graph
//------------------------------------------------------------//
TGraph* getGraph(TFile* file, TString pname, int color, int style)
{
  
  TGraph* gr = (TGraph*) (file->Get(pname.Data()));
  gr->SetLineColor(color);
  gr->SetLineStyle(style);
  return gr;

}

//------------------------------------------------------------//
// Get 2D histogram
//------------------------------------------------------------//
TH2F* getHist2(TFile* file, TString pname, TString xtitle,
	       TString ytitle, TString ztitle)
{
  
  TH2F* hist = (TH2F*) (file->Get(pname.Data())->Clone());
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitle.Data());
  hist->GetZaxis()->SetTitle(ztitle.Data());
  hist->SetTitle("");
  hist->SetStats(0);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetZaxis()->SetTitleOffset(1.);
  hist->SetLineWidth(2);

  return hist;

}

//------------------------------------------------------------//
// Make Histogram
//------------------------------------------------------------//
TH1F* makeHist(TString name, int nbins, float xmin, float xmax,
	       TString xtitle, TString ytitle, int color, int marker)
{

  TH1F* h = new TH1F(name.Data(),"",nbins,xmin,xmax);
  setHistAtt(h,xtitle,ytitle,color,marker);
  return h;

}

//------------------------------------------------------------//
// Set histogram attributes
//------------------------------------------------------------//
void setHistAtt(TH1F* &h, TString xtitle, TString ytitle, 
		int color, int marker)
{

  h->GetXaxis()->SetTitle(xtitle.Data());
  h->GetYaxis()->SetTitle(ytitle.Data());
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetStats(0);
  h->SetTitle("");
  h->GetYaxis()->SetTitleOffset(1.5);

}

//----------------------------------------//
// Format graph
//----------------------------------------//
void formatGraph(TGraph* &gr, TString xtitle,
		 TString ytitle, int color)
{

  gr->SetTitle("");
  gr->GetYaxis()->SetTitle(ytitle.Data());
  gr->GetXaxis()->SetTitle(xtitle.Data());
  gr->SetLineColor(color);
  gr->GetYaxis()->SetTitleOffset(1.5);
  gr->SetLineWidth(2);
}

//----------------------------------------//
// Get Fourier Transform to freq space
//----------------------------------------//
void getFTF(TGraph* gr,
	    TGraph* &mag,        // placeholder for mag
	    TGraph* &phase,      // placeholder for phase
	    float min,           // freq in MHz 
	    float stepSize,      // Step size in MHz
	    const int nsteps)

{
  
  // Make place hodlers
  Double_t ReE[nsteps];
  Double_t ImE[nsteps];
  Double_t EF[nsteps];
  Double_t Phase[nsteps];
  Double_t freq[nsteps];
  Double_t pi = TMath::Pi();
  
  // Necessary for looping
  double E=0, t=0;
  int np = gr->GetN();
  int filledPoints = 0;
  
  // Get the time step info
  double dt = 0; 
  gr->GetPoint(2,dt,E);
  gr->GetPoint(1,t,E); 
  dt = dt -t;

  // Now loop and calculate
  for(int is=0; is<nsteps; ++is){
    ReE[is]  = 0;
    ImE[is]  = 0;
    freq[is] = min + is*stepSize;
    filledPoints += 1;
    
    // Now loop and calculate the real 
    // and imaginary components
    for(int ip=0; ip<np; ++ip){
      gr->GetPoint(ip,t,E);
      ReE[is] += E * cos(2*pi*freq[is]*1e6*t*1e-9) * dt * 1e-9 * 1e6;
      ImE[is] += E * sin(2*pi*freq[is]*1e6*t*1e-9) * dt * 1e-9 * 1e6;
    }// end loop over times

    EF[is]    = sqrt(ReE[is]*ReE[is] + ImE[is]*ImE[is]);
    if( ReE[is] != 0 ) Phase[is] = atan(ImE[is]/ReE[is]);
    else               Phase[is] = 0;
  }// end loop over freq points
  
  mag   = new TGraph(filledPoints, freq, EF);
  phase = new TGraph(filledPoints, freq, Phase);

}


//----------------------------------------//
// Get Fourier Transform to freq space
//----------------------------------------//
void getFTT(TGraph* gr,
	    TGraph* inPhase,
	    float min,       // time in ns 
	    float stepSize,  // Step size in ns
	    const int nsteps,
	    TGraph* &mag,
	    TGraph* &phase)
{
  
  // Make place hodlers
  Double_t ReE[nsteps];
  Double_t ImE[nsteps];
  Double_t ET[nsteps];
  Double_t Phase[nsteps];
  Double_t time[nsteps];
  Double_t pi = TMath::Pi();

  // Necessary for looping
  double E=0, f=0, inphase=0;
  int np = gr->GetN();
  int filledPoints = 0;
  
  // Get the time step info
  double df = 0; 
  gr->GetPoint(2,df,E);
  gr->GetPoint(1,f,E); 
  df = df -f;

  // Now loop and calculate
  for(int is=0; is<nsteps; ++is){
    ReE[is]  = 0;
    ImE[is]  = 0;
    time[is] = min + is*stepSize;
    filledPoints += 1;
    
    // Now loop and calculate the real 
    // and imaginary components
    for(int ip=0; ip<np; ++ip){
      gr->GetPoint(ip,f,E);
      inPhase->GetPoint(ip,f,inphase);
      ReE[is] += E * cos(2*pi*f*1e6*time[is]*1e-9 + inphase) * df; // df in MHz
      ImE[is] -= E * sin(2*pi*f*1e6*time[is]*1e-9 + inphase) * df; // df in MHz
    }// end loop over times

    ET[is] = sqrt(ReE[is]*ReE[is] + ImE[is]*ImE[is]);
    if( ReE[is] != 0 ) Phase[is] = atan(ImE[is]/ReE[is]);
    else               Phase[is] = 0;

  }// end loop over freq points

  mag = new TGraph(filledPoints, time, ET);
  phase = new TGraph(filledPoints, time, Phase);

}
