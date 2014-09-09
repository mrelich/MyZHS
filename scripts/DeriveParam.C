
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Try to derive Jaime's parameterization for our 40 MeV showers. //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "myHist.C"
#include <sstream>

bool m_fix = false;
TString savedir = "plots/";

//-----------------------------------------//
// Main
//-----------------------------------------//
void DeriveParam(TString fname,
		 TString savename,
		 double beamE, // in TeV
		 int nEvents,
		 bool fix=false)
{

  // Open up the file
  TFile* file = new TFile(fname.Data());
  
  // decide if we shoudl fix parameters
  m_fix = fix;

  // Fit function to relevant range
  findFit(file, savename, beamE);

}

//-----------------------------------------//
// Fit and draw function
//-----------------------------------------//
void findFit(TFile* file,
	     TString savename,
	     double Energy)
{

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  TPad* top = NULL;
  TPad* bot = NULL;
  makePads(c,top,bot);

  // Read in Profile
  TProfile* prof = getProfile(file,"VP_avg","time [ns]",
			      "|RA(#theta_{C},t)|", kBlue, 20);
  
  // Fit function to profile
  TF1* func = fitProfile(prof, Energy);
  //TF1* func = fitProfileM(prof, Energy);

  // Specify some maximum value
  float maximum = prof->GetMaximum();
  prof->SetMaximum( maximum*5 );
  prof->SetMinimum(maximum * 1e-3);

  // Set the range to draw
  float xmin = -1;
  float xmax = 0.8;
  prof->GetXaxis()->SetRange( prof->GetXaxis()->FindBin(xmin),
			      prof->GetXaxis()->FindBin(xmax));

  // Set some misc stats
  prof->SetLabelSize(0,"X");
  prof->SetLabelSize(0.05,"Y");
  prof->GetYaxis()->SetTitleSize(0.055);
  prof->GetYaxis()->SetTitleOffset(0.9);

  // Now draw
  c->cd();
  top->Draw();
  top->cd();
  top->SetLogy();
  prof->Draw();
  func->Draw("same");
  cout<<"Maximum: "<<(prof->ProjectionX("temp"))->GetMaximum()
      <<" Func: "<<func->Eval(1e-6)
      <<" Func2: "<<func->Eval(-1e-6)<<endl;

  // Now add Fit attributes
  TLatex* lat = makeLatex();
  lat->SetTextSize(0.05);
  TString f1 = "N*E*(e^{-|t|/a}(1+b|t|)^{c}) t > 0";
  TString f2 = "N*E*(e^{-|t|/d}(1+e|t|)^{f}) t < 0";
  string params[] = {"Energy","a","b","c","d","e","f","N"};
  int nParams = 8;
  float xi = 0.18;
  float yi = 0.8;
  float dy = 0.055;
  
  lat->DrawLatex(xi,0.9,f1.Data());
  lat->DrawLatex(xi,0.85,f2.Data());
  stringstream ss;
  for(int i=0; i<nParams; ++i){
    ss.str("");
    ss << params[i] << " = " << func->GetParameter(i);
    lat->DrawLatex(xi,yi-i*dy,ss.str().c_str());
  }

  top->Update();

  // Now set bottom ratio
  TH1D* ratio = makeRatio(prof,func);
  ratio->GetXaxis()->SetRange( ratio->FindBin(xmin), ratio->FindBin(xmax) );

  // Add some more lines
  TLine* u1  = makeLine(xmin,xmax,1.0,1.0,kBlack,2);
  TLine* u2   = makeLine(xmin,xmax,2,2,kBlack,2);
  TLine* u3   = makeLine(xmin,xmax,3,3,kBlack,2);


  // Draw ratio
  c->cd();
  bot->Draw();
  bot->cd();
  ratio->Draw();
  u1->Draw("same");
  u2->Draw("same");
  u3->Draw("same");
  ratio->Draw("same");
  bot->Update();

  c->SaveAs((savedir+savename+".png").Data());
  
}

//-----------------------------------------//
// Make ratio histogram
//-----------------------------------------//
TH1D* makeRatio(TProfile* prof,
		TF1* fit)
{

  TH1D* ratio = prof->ProjectionX("ratio");
  
  // Now divide out the fit
  int nbins = ratio->GetNbinsX();
  for(int bin=1; bin<=nbins; ++bin){
    float bc = ratio->GetBinCenter(bin);
    float co = ratio->GetBinContent(bin);
    ratio->SetBinContent(bin, co/fit->Eval(bc));
  }

  // Set some attributes
  ratio->GetYaxis()->SetTitle("ZHS/Fit");
  ratio->GetXaxis()->SetTitle("time [ns]");

  // Pretty
  ratio->SetMaximum(4);
  ratio->SetMinimum(0);
  ratio->SetStats(0);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetYaxis()->SetTitleSize(0.14);
  ratio->GetXaxis()->SetTitleOffset(1.2);
  ratio->GetYaxis()->SetTitleOffset(0.35);
  ratio->GetYaxis()->CenterTitle();
  ratio->SetLabelSize(0.12,"X");
  ratio->SetLabelSize(0.12,"Y");
  ratio->GetYaxis()->SetNdivisions(405);
  ratio->SetMarkerSize(0.5);
  ratio->SetMarkerStyle(20);

  return ratio;

}

//-----------------------------------------//
// Fit user function to input profile
//-----------------------------------------//
TF1* fitProfile(TProfile* prof, double Energy)
{

  // Define function
  TF1* func = new TF1("f",JaimeParam,-0.5,0.5,8);
  func->FixParameter(0,Energy);
  func->FixParameter(3,-3);
  func->FixParameter(6,-4);

  // If fix parameters, just set them and
  // return
  if( m_fix ){
    func->FixParameter(1,5.67715e-02);
    func->FixParameter(2,8.54433e+00);
    func->FixParameter(4,3.83746e-02);
    func->FixParameter(5,8.61173e+00);
    func->FixParameter(7,2.28959e-13);
    return func;
  }

  // Otherwise now we fit
  // Give some starting points
  func->SetParameter(1,0.057);
  func->SetParameter(2,2.87);
  func->SetParameter(4,0.03);
  func->SetParameter(5,3.05);

  // Maybe let these float...?
  //func->SetParameter(3,-3);
  //func->SetParameter(6,-3.5);

  // Fit this function
  prof->Fit("f","","",-0.15,0.2);
  
  // Return the function
  return func;

}

//-----------------------------------------//
// Fit user function to input profile
//-----------------------------------------//
TF1* fitProfileM(TProfile* prof, double Energy)
{

  // Define function
  TF1* func = new TF1("f",MattParam,-0.5,0.5,8);
  func->FixParameter(0,Energy);

  // Otherwise now we fit
  // Give some starting points
  func->SetParameter(1,2);
  func->SetParameter(2,0.02);
  func->SetParameter(3,2);
  func->SetParameter(4,0.03);

  // Fit this function
  prof->Fit("f","","",-0.15,0.2);
  
  // Return the function
  return func;

}

//-----------------------------------------//
// User fit function
//-----------------------------------------//
Double_t JaimeParam(Double_t *x, Double_t *par)
{

  // User will input the time
  Double_t t = x[0];

  // Now specify the user parameters
  Double_t E  = par[0];
  Double_t P1 = par[1];
  Double_t P2 = par[2];
  Double_t P3 = par[3];
  Double_t P4 = par[4];
  Double_t P5 = par[5];
  Double_t P6 = par[6];
  Double_t N  = par[7];

  // Now define piecwise function
  if( t > 0 )
    return N * E * (TMath::Exp(-fabs(t)/P1) + pow(1+P2*fabs(t),P3));
  else
    return N * E * (TMath::Exp(-fabs(t)/P4) + pow(1+P5*fabs(t),P6));
  
}

//-----------------------------------------//
// User fit function
//-----------------------------------------//
Double_t MattParam(Double_t *x, Double_t *par)
{
  
  // User will input the time
  Double_t t = x[0];

  // Now specify the user parameters
  Double_t E  = par[0];
  Double_t P1 = par[1];
  Double_t P2 = par[2];
  Double_t P3 = par[3];
  Double_t P4 = par[4];
  Double_t N  = par[7];

  // Now define piecwise function
  if( t > 0 )
    return N * E * 1/(P1 + pow(fabs(t),P2));
  else
    return N * E * 1/(P3 + pow(fabs(t),P4));
  
}
