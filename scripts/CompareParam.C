
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// This will compare the parameterization to what ZHS actually   //
// gives.  This is meant to be a sanity check for lower energies //
// in order to see if things are relatively accurate.            //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"

//---------------------------------//
// Main
//---------------------------------//
void CompareParam()
{

  vector<TString> fnames;  vector<double> Energies;
  //fnames.push_back("beam100e3MeV");  Energies.push_back(0.1);
  //fnames.push_back("beam1e6MeV");    Energies.push_back(1);
  //fnames.push_back("beam10e6MeV");    Energies.push_back(10);
  //fnames.push_back("beam40MeV_prof");     Energies.push_back(1000 * 40e-6);
  
  //for(unsigned int i=0; i<fnames.size(); ++i)
  //  plot(fnames.at(i), Energies.at(i));
  
  plotPeak("beam100e3MeV.root");
  
}

//---------------------------------//
// Compare peak value with
// the parameterization
//---------------------------------//
void plotPeak(TString fname)
{

  // Open the input file
  TFile* file = new TFile(("rootfiles/"+fname).Data());
  
  // make Canvas
  TCanvas* c = new TCanvas("c","",600,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  //c->SetLogy();
  
  // Make histogram to store peak 
  TH1F* h_peak = new TH1F("peak","",100,0,30);

  // The maximum
  float maximum = 4.5e-14 * 0.1;
  
  // Now loop over and get the difference
  TH1F* h = NULL;
  for(int i=0; i<50; ++i){
    h = (TH1F*) file->Get(Form("VP_evt%i",i));
    h_peak->Fill(h->GetMaximum()/maximum);
    cout<<"Max: "<<h->GetMaximum()<<" "<<maximum<<endl;
  }

  h_peak->Draw();

}

//---------------------------------//
// Plot a given event
//---------------------------------//
void plot(TString fname, double E)
{
 
  // Specify input file names
  TFile* file = new TFile((fname+".root").Data());
  
  // Specify the plot name
  TString pname = "RA";

  // Make canvas
  TCanvas* c = new TCanvas("c","",600,500);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetLogy();

  // Histogram
  //TH1F* hist = (TH1F*) file->Get(pname.Data());
  TProfile* hist = (TProfile*) file->Get(pname.Data());
  hist->SetMinimum( hist->GetMaximum() * 1e-3 );
  hist->SetStats(0);
  hist->SetTitle("");
  hist->GetYaxis()->SetTitle("|RA(t)| [V s]");
  hist->GetXaxis()->SetTitle("time [ns]");
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->Draw();

  TF1* res1TeV = new TF1("func",RA,-2,2,1);
  res1TeV->SetParameter(0,E);
  res1TeV->Draw("same");

  c->SaveAs((fname+"_VP.png").Data());
  //delete c;
}

//---------------------------------//
// Define the function
//---------------------------------//
Double_t RA(Double_t *x, 
	    Double_t *pars)
{

  // This is the time in ns
  double t = x[0];
  
  // Specify the energy in TeV
  double E  = pars[0];
  
  // Now compute result
  double result = -4.5e-14 * E;

  if( t > 0 )
    result *= (exp(-fabs(t)/0.057) + pow(1+2.87*fabs(t),-3));
  else
    result *= (exp(-fabs(t)/0.030) + pow(1+3.05*fabs(t),-3.5));
  
  return fabs(result);

}

//-------------------------------------------//
// Fit Greisen parameterization for shower
//-------------------------------------------//
TF1* fitGreisen(TProfile* prof, float E0, int color, int style)
{

  // Function 7 from this paper:
  //    * http://prd.aps.org/pdf/PRD/v65/i10/e103002
  // Function of shower energy, energy, and shower 
  // depth.

  // [0] -- A(E)
  // [1] -- a(E)

  // y term
  stringstream y;
  y << "TMath::Log("
    << E0 << "/" << E_critical
    << ")";

  // Constant Term
  stringstream C;
  C << "(0.31*[0])/(TMath::Sqrt("
    << y.str() << "))";


  // Define t1
  string t1 = "(x+[1])";  
    
  // Ln(s1)
  stringstream lns1;
  lns1 << "TMath::Log(3*"
       << t1
       << "/(" 
       << t1 << "+2*" << y.str()
       <<"))";

  // expone
  stringstream exp;
  exp << "TMath::Exp("
      << t1 <<"*(1-1.5*"
      << lns1.str()
      <<"))";
  
  //cout <<"Function: "<<C.str()+"*"+exp.str()<<endl;

  // Function
  TF1 func = TF1("greisen",(C.str()+"*"+exp.str()).c_str(),0,20);
  func.SetParameter(0,0.5);
  func.SetParameter(1,1.);
  //func.FixParameter(1,0.76);

  // Fit
  prof->Fit("greisen","RQ");
  TF1* fit = prof->GetFunction("greisen");
  fit->SetLineColor(color);
  fit->SetLineStyle(style);

  return fit;
  
}
