
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Script to plot the angular dependence of the vector potential. //
// The goal is to compare how much variation we get away from the //
// Cherenkov angle.                                               //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "myHist.C"

bool m_save = true;
//bool m_save = false;
float m_X0  = 0.38;
float m_step = m_X0 * 0.1; // step size is 1 rad length
float E_critical = 79.0255;

//--------------------------------------//
// Main
//--------------------------------------//
void CompareAng()
{

  // Specify the input file
  //TFile* file = new TFile("rootfiles/beam100e3MeV_1Prim_10psbin_Angular.root");
  TFile* file = new TFile("rootfiles/testing10e6MeV_Angular.root");
  
  // Specify the Angles
  vector<TString> angles;
  angles.push_back("10");
  angles.push_back("20");
  angles.push_back("30");
  angles.push_back("40");
  angles.push_back("50");
  angles.push_back("55.829616");
  angles.push_back("60");
  angles.push_back("70");
  angles.push_back("80");
  //angles.push_back(88);
  
  // Plot all VP on a single plot
  //plotAllAngles(file, angles);
  //plotMaxVsAngle(file, angles);
  //plotAngleVsParam(file, "60");
  plotAngleVsParam(file, "60");
  //for(int i=0; i<angles.size(); ++i)
  //plotAngleVsParam(file, angles.at(i));
  
}

//--------------------------------------//
// Compare given angle with param
//--------------------------------------//
void plotAngleVsParam(TFile* file,
		      TString angle)
{

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();
  
  // Add Legend
  TLegend* leg = makeLegend(0.15,0.3,0.85,0.9);
  TString header ="Angle = "+angle+" [deg]"; 
  leg->SetHeader(header.Data());

  // Get profile of the result
  TString pname = "VP_avg_"+angle;
  TString xtitle = "time [ns]";
  TString ytitle = "|RA(#theta,t)| [Vs]";
  TProfile* prof = getProfile(file,pname,xtitle,ytitle,kBlue,20);
  
  // Now get Histogram of a fitted result 
  // Using similar binning
  int nbins = prof->GetNbinsX();
  TH1F* h_param = getParamResult(atof(angle.Data()) + 0.02, // 0.02 due to rounding
				 nbins,
				 prof->GetXaxis()->GetBinLowEdge(0),
				 prof->GetXaxis()->GetBinLowEdge(nbins)+
				 prof->GetXaxis()->GetBinWidth(nbins));


  // Set Minimum and max
  float maximum = prof->GetMaximum();
  prof->SetMaximum(maximum * 1.8);
  prof->SetMinimum(maximum * 1e-3);

  prof->Draw();
  h_param->Draw("same");
  leg->Draw("same");

  if(m_save)
    //c->SaveAs(("plots/ParamVsActual_GreissenShower10TeV_Angle"+angle+".png").Data());
    //c->SaveAs(("plots/ParamVsActual_GreissenShower10TeV_1X0Step_Angle"+angle+".png").Data());
    c->SaveAs(("plots/ParamVsActual_GreissenShower10TeV_Angle_fixed"+angle+".png").Data());
  //delete c;
}

//--------------------------------------//
// Compare fit result
//--------------------------------------//
TH1F* getParamResult(float angle,
		     int nbins,
		     float tmin,
		     float tmax)
{

  // Make the histogram for result
  TH1F* result = makeHist("result",nbins,tmin,tmax,"","",kRed,25);
  
  // Load the charge excess profile
  //TFile* f_Qz  = new TFile("rootfiles/beam100e3MeV_1Prim_10psbin_NPart.root");
  //TProfile* Qz = (TProfile*) f_Qz->Get("NPartDiff");
  //TFile* f_Qz  = new TFile("rootfiles/TrkAna_100_100000_ice_eBeam_1MeV.root");
  TFile* f_Qz  = new TFile("rootfiles/testing10e6MeV_NPart.root");
  TProfile* Qz = (TProfile*) f_Qz->Get("NPartDiff");
  TF1* fitQz  = fitGreisen(Qz, 10e6, kRed, 2);

  TCanvas* cfit = makeCanvas("cfit");
  cfit->SetLogy();
  Qz->SetTitle("");
  Qz->SetStats(0);
  Qz->Draw();
  fitQz->Draw("same");
  cfit->SaveAs("plots/GreissenFit.png");
  delete cfit;

  // Now loop over the bins and get the result
  for(int bin=1; bin<=nbins; ++bin){
    float time = result->GetBinCenter(bin);
    //float Vp   = paramResult(Qz,time,angle);
    float Vp   = paramResult(fitQz,time,angle);
    result->Fill(time, Vp);
  }

  return result;

}

//--------------------------------------//
// Parameterized result
//--------------------------------------//
float paramResult(TProfile* Qz,
		  float t,
		  float angle)
{

  // Step size is 1 rad length in ize
  float step = m_step; //0.38*0.5; // [m]

  // Calculate the LQtot
  float LQtot = getLQtot(Qz, step);
  
  // Get constant
  float C = TMath::Sin(angle) / (LQtot * TMath::Sin(55.829616));
  
  // Now loop over charge and calculate
  float integral = 0;
  float z = 0, RA = 0, qz=0, tR = 0;
  int nbins = Qz->GetNbinsX();
  for(int bin = 1; bin<=nbins; ++bin){
    qz = Qz->GetBinContent(bin);
    z = step * (bin);
    tR = getTR(t,z,angle)*1e9;
    //cout<<"tR: "<<tR<<" time: "<<t<<endl;
    RA = getRA(tR);
    integral += step * qz * RA;
  }

  return fabs(C*integral);

}

//--------------------------------------//
// Parameterized result
//--------------------------------------//
float paramResult(TF1* fitQz,
		  float t,
		  float angle)
{

  // Now the step size can be in arbitrary
  // units of radiation length
  float ds = 0.1;
  float step = m_X0 * ds;
  int nsteps = (int) 40/ds;

  // Calculate the LQtot
  float LQtot = getLQtot(fitQz, step, nsteps);
  
  // Get constant
  float C = TMath::Sin(angle*TMath::Pi()/180) / (LQtot * TMath::Sin(55.829616*TMath::Pi()/180.));
  
  // Now loop over charge and calculate
  float integral = 0;
  float z = 0, RA = 0, qz=0, tR = 0;
  for(int i = 0; i<nsteps; ++i){
    qz = fitQz->Eval(i*step/m_X0);
    z = step * (i);
    tR = getTR(t,z,angle)*1e9;
    RA = getRA(tR);
    //cout<<"\tz: "<<z<<" qz: "<<qz<<" RA: "<<RA<<" time: "<<t<<" tR: "<<tR<<endl;
    integral += step * qz * RA;
  }

  //assert(false);
  return fabs(C*integral);

}

//--------------------------------------//
// Get the retarded time
//--------------------------------------//
float getTR(float time, float z, float angle)
{

  double c = 3e8;
  double v = c;
  double n = 1.78;
  //return time * 1e-9 - z/v + z*n*TMath::Cos(angle*TMath::Pi()/180.)/c;
  return (time * 1e-9 - z/v + z*n*TMath::Cos(angle*TMath::Pi()/180.)/c);

}

//--------------------------------------//
// getRA(time)
//--------------------------------------//
float getRA(float time)
{
  
  //float C = -4.5e-14 * 100 * 1e-3; // 100 GeV shower assumed
  float C = -4.5e-14 * 10;          // 10 TeV shower assumed
  
  float tDep = 0;
  if( time > 0)
    tDep = TMath::Exp(-fabs(time)/0.057) + pow(1+2.87*fabs(time),-3);
  else
    tDep = TMath::Exp(-fabs(time)/0.030) + pow(1+3.05*fabs(time),-3.5);
  
  return C* tDep;
  
}

//--------------------------------------//
// Get LQ tot
//--------------------------------------//
float getLQtot(TProfile* Qz, float step)
{
  
  float tot = 0;
  int nbins = Qz->GetNbinsX();
  for(int bin=1; bin<=nbins; ++bin){
    tot += step * Qz->GetBinContent(bin);
  }

  return tot;

}

//--------------------------------------//
float getLQtot(TF1* Qz, float step, int nsteps)
{
  
  float tot = 0;
  for(int i=0; i<nsteps; ++i)
    tot += step * Qz->Eval(step*i/m_X0);

  return tot;
  
}

//--------------------------------------//
// Plot Maximum VP vs. angle
//--------------------------------------//
void plotMaxVsAngle(TFile* file, 
		    vector<TString> angles)
{

  // Make Canvas
  TCanvas* c = makeCanvas("c_max");
  c->SetLogy();

  // This time save things to TGraph
  int nps = (int) angles.size();
  float ang[1000];
  float VP[1000];
  float p_VP[1000];

  // Loop over angles and get maximum
  TProfile* prof = NULL;
  for(unsigned int i=0; i<angles.size(); ++i){
    TString pname = "VP_avg_"+angles.at(i);
    
    prof = getProfile(file, pname, "", "", kBlack, 20);
    
    ang[i] = atof(angles.at(i).Data());
    VP[i] = prof->GetMaximum();
    
    // Get parameterized on
    int nbins = prof->GetNbinsX();
    TH1F* h_param = getParamResult(atof(angles.at(i).Data()),
				   nbins,
				   prof->GetXaxis()->GetBinLowEdge(0),
				   prof->GetXaxis()->GetBinLowEdge(nbins)+
				   prof->GetXaxis()->GetBinWidth(nbins));
    
    p_VP[i] = h_param->GetMaximum();
    delete h_param;
  }

  // make graph
  TGraph* gr_ang = new TGraph(nps,ang,VP);
  gr_ang->GetXaxis()->SetTitle("Angle [deg]");
  gr_ang->GetYaxis()->SetTitle("|RA(#theta,t_{max})| [Vs]");
  gr_ang->SetLineColor(kBlue);
  gr_ang->SetMarkerStyle(20);
  gr_ang->SetMarkerColor(kBlue);
  gr_ang->SetTitle("");

  // Make another graph
  TGraph* gr_param = new TGraph(nps, ang, p_VP);
  gr_param->SetMarkerStyle(25);
  gr_param->SetMarkerColor(kRed);
  gr_param->SetLineColor(kRed);
  gr_param->SetMarkerSize(1);

  // Draw
  gr_ang->Draw();
  gr_param->Draw("samepl");
  
  if(m_save)
    //c->SaveAs("plots/AngDepSummary_100GeV.png");
    //c->SaveAs("plots/AngDepSummary_10TeV.png");
    c->SaveAs("plots/AngDepSummary_10TeV_fixed.png");

}

//--------------------------------------//
// Plot all the VPs 
//--------------------------------------//
void plotAllAngles(TFile* file, 
		   vector<TString> angles)
{

  // First specify canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  // Now a legend
  TLegend* leg = makeLegend(0.7,0.85,0.5,0.93);

  // Specify some histogram attributes
  TString xtitle = "time [ns]";
  TString ytitle = "|RA(#theta,t)| [Vs]";
  int colors[] = {kRed, kBlue, kAzure-2, kGreen-3, kMagenta,
		  kBlack, kCyan+2, kViolet+3, kOrange+7};
  
  // Now plot
  TProfile* profs[10];
  float maximum = -999;
  for(unsigned int i=0; i<angles.size(); ++i){

    // set angle and plot names
    TString angle = angles.at(i);
    TString pname = "VP_avg_"+angle;
  
    // Plot and add to legend
    profs[i] = getProfile(file, pname, xtitle, ytitle, colors[i], 20);   
    leg->AddEntry(profs[i], angle.Data(), "lep");
    
    // Record maximum
    if(maximum < profs[i]->GetMaximum())
      maximum = profs[i]->GetMaximum();
    
  }
  
  // set maximum
  profs[0]->SetMaximum(maximum*10);
  profs[0]->SetMinimum(maximum*1e-3);
  profs[0]->Draw();
  for(unsigned int i=1; i<angles.size(); ++i)
    profs[i]->Draw("same");
  
  leg->Draw("same");

  if(m_save)
    c->SaveAs("plots/AllAngularResult_100GeV.png");

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
  TF1 func = TF1("greisen",(C.str()+"*"+exp.str()).c_str(),0,40);
  func.SetParameter(0,0.5);
  func.SetParameter(1,1.);
  //func.FixParameter(1,0.76);

  prof->Fit("greisen","RQ");
  TF1* fit = prof->GetFunction("greisen");
  fit->SetLineColor(color);
  fit->SetLineStyle(style);

  return fit;
  
}
