#include "TLine.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TEfficiency.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
using namespace std;

const int val = 4;
//const bool save_plots = true;
void setTDRStyle();

void limit(bool only_leptonTAG=false, bool save_plots = false, float min=0.006, float max=900., TString plot_name="Limit_prova.pdf") {
  setTDRStyle();
  TCanvas* c1  = new TCanvas("c1","c1",0,0,800,600); //c1->SetLogy();
  c1->SetLogy();
  save_plots = true;
  plot_name="limit.png";

  Double_t XMass[val]    = {1000, 1500, 2000, 2500};
	
  //LIMIT - Breit-Wigner
  Double_t Limit_Obs[val]     = {0,0,0,0};
  Double_t Limit_Exp_m2s[val] = {2.1936, 0.8895, 0.5064, 0.4483};
  Double_t Limit_Exp_m1s[val] = {3.0420, 1.2656, 0.7397, 0.6608};
  Double_t Limit_Exp[val]     = {4.5288, 1.9714, 1.1627, 1.0529};
  Double_t Limit_Exp_p1s[val] = {6.8214, 3.2445, 1.9135, 1.7327};
  Double_t Limit_Exp_p2s[val] = {9.7806, 5.1998, 3.0161, 2.7311};

  //LIMIT - GaussXLandau
  Double_t Limit_Obs2[val]     = {0,0,0,0};
  Double_t Limit_Exp_m2s2[val] = {2.2500 , 8.0346 , 27.1690 , 63.4689 };
  Double_t Limit_Exp_m1s2[val] = {3.1133 , 11.0869, 38.2007 , 90.3546 };
  Double_t Limit_Exp2[val]     = {4.6265 , 16.2598, 58.2031 , 138.2812};
  Double_t Limit_Exp_p1s2[val] = {7.0792 , 24.7502, 93.0032 , 219.3069};
  Double_t Limit_Exp_p2s2[val] = {10.3827, 36.8075, 148.7738, 337.1894};

  Double_t XSEC[val] = {1.5, 0.176, 0.031, 0.003};

  for(int i=0; i<val; i++){
    Limit_Obs[i]     = Limit_Obs[i]       ;//* XSEC[i];
    Limit_Exp_m2s[i] = Limit_Exp_m2s[i]	  ;//* XSEC[i];
    Limit_Exp_m1s[i] = Limit_Exp_m1s[i]	  ;//* XSEC[i];
    Limit_Exp[i]     = Limit_Exp[i]    	  ;//* XSEC[i];
    Limit_Exp_p1s[i] = Limit_Exp_p1s[i]	  ;//* XSEC[i];
    Limit_Exp_p2s[i] = Limit_Exp_p2s[i]   ;//* XSEC[i]; 
    Limit_Obs2[i]     = Limit_Obs2[i]     ;//* XSEC[i];
    Limit_Exp_m2s2[i] = Limit_Exp_m2s2[i] ;//* XSEC[i];
    Limit_Exp_m1s2[i] = Limit_Exp_m1s2[i] ;//* XSEC[i];
    Limit_Exp2[i]     = Limit_Exp2[i]     ;//* XSEC[i];
    Limit_Exp_p1s2[i] = Limit_Exp_p1s2[i] ;//* XSEC[i];
    Limit_Exp_p2s2[i] = Limit_Exp_p2s2[i] ;//* XSEC[i]; 
  }
	
  //Draw limit
  TLegend *legendr = new TLegend(0.38,0.70,0.87,0.92);
  legendr->SetTextSize(0.030); 
  legendr->SetFillColor(0);
  legendr->SetBorderSize(1);
	
  //2 sigma 
  TGraph *band_2s = new TGraph(2*val+1);
  for(int i=0;i<val;i++){
    band_2s->SetPoint(i,XMass[i],Limit_Exp_p2s[i]);
    band_2s->SetPoint(i+val,XMass[val-1-i],Limit_Exp_m2s[val-1-i]);
  }
  band_2s->SetPoint(2*val,XMass[0],Limit_Exp_p2s[0]);
  band_2s->SetLineStyle(2);
  band_2s->SetFillColor(kYellow);
  band_2s->Draw("ALF2");
  band_2s->SetMinimum(min);
  band_2s->SetMaximum(max);
  band_2s->SetTitle(0);
  band_2s->GetXaxis()->SetTitle("M(Z,H) [GeV]");
  band_2s->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z' #rightarrow ZH) x BR(Z #rightarrow qq) X BR(H #rightarrow #tau#tau) [fb]");
  band_2s->GetXaxis()->SetTitleOffset(1.1);
  band_2s->GetYaxis()->SetTitleOffset(1.3);
  band_2s->GetXaxis()->SetRangeUser(1000,2500);
  band_2s->GetYaxis()->SetTitleSize(0.040);
  band_2s->GetXaxis()->SetNdivisions(1005);
  band_2s->Draw("LF2");
	
  //1 sigma
  TGraph *band_1s = new TGraph(2*val+1);
  for(int i=0;i<val;i++){
    band_1s->SetPoint(i,XMass[i],Limit_Exp_p1s[i]);
    band_1s->SetPoint(i+val,XMass[val-1-i],Limit_Exp_m1s[val-1-i]);
  }
  band_1s->SetPoint(2*val,XMass[0],Limit_Exp_p1s[0]);
  band_1s->SetLineStyle(2);
  band_1s->SetFillColor(kGreen);
  band_1s->Draw("LF2");
  //Expected
  TGraph *limit_exp = new TGraph(val,XMass,Limit_Exp);
  limit_exp->SetLineWidth(2);
  limit_exp->SetLineStyle(2);
  limit_exp->SetMarkerSize(1.3);
  limit_exp->Draw("LP");
  TGraph *limit_exp2 = new TGraph(val,XMass,Limit_Exp2);
  limit_exp2->SetLineWidth(2);
  limit_exp2->SetLineStyle(2);
  limit_exp2->SetLineColor(kBlue);
  limit_exp2->SetMarkerColor(kBlue);
  limit_exp2->SetMarkerSize(1.3);
  //limit_exp2->Draw("LPsame");
  //Observed
  TGraph *limit_obs = new TGraph(val,XMass,Limit_Obs);
  limit_obs->SetMarkerStyle(21);
  //limit_obs->Draw("LP");

  TGraph *xsec = new TGraph(val,XMass,XSEC);
  xsec->SetMarkerStyle(21);
  xsec->SetLineColor(2);
  xsec->SetMarkerColor(2);
  xsec->SetLineWidth(2);
  xsec->Draw("LP");

  TGraph *band_1s2 = new TGraph(2*val+1);
  TGraph *band_2s2 = new TGraph(2*val+1);
  band_1s2->SetLineStyle(2);
  band_1s2->SetFillColor(kGreen);
  band_1s2->SetLineColor(kBlue);
  band_1s2->SetMarkerColor(kBlue);
  band_2s2->SetLineStyle(2);
  band_2s2->SetFillColor(kYellow);
  band_2s2->SetLineColor(kBlue);
  band_2s2->SetMarkerColor(kBlue);

  //Legend
  legendr->AddEntry(xsec,"Theoretical cross-section","LP");
  legendr->AddEntry(band_1s, "Expected #pm 1#sigma","flP");
  legendr->AddEntry(band_2s, "Expected #pm 2#sigma","flP");
  //legendr->AddEntry(band_1s2,"Expected #pm 1#sigma - GaussianXLandau","flP");
  //legendr->AddEntry(band_2s2,"Expected #pm 2#sigma - GaussianXLandau","flP");
  legendr->Draw();
	
  //Save plot
  if(save_plots) c1->SaveAs(plot_name); 
}
/////
//   Set setTDRStyle_modified (from link https://twiki.cern.ch/twiki/pub/CMS/TRK10001/setTDRStyle_modified.C)
/////
void setTDRStyle(){
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);
  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);
  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  tdrStyle->SetErrorX(0.);
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("e");
  tdrStyle->SetStatColor(kGray);
  tdrStyle->SetStatFont(42);

  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(0);
  tdrStyle->SetStatX(1.); //Starting position on X axis
  tdrStyle->SetStatY(1.); //Starting position on Y axis
  tdrStyle->SetStatFontSize(0.025); //Vertical Size
  tdrStyle->SetStatW(0.15); //Horizontal size 
  // tdrStyle->SetStatStyle(Style_t style = 1001);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.125);
  tdrStyle->SetPadLeftMargin(0.105);
  tdrStyle->SetPadRightMargin(0.1);

  // For the Global title:

  //  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.05, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.8);
  tdrStyle->SetTitleOffset(0.7, "Y"); // Another way to set the Offset

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.045, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  tdrStyle->cd();
}



/*
//2 APRILE
//LIMIT - Breit-Wigner
Double_t Limit_Obs[val]     = {0,0,0,0};
Double_t Limit_Exp_m2s[val] = {2.5857 , 6.4391 , 27.5311 , 71.5248 };
Double_t Limit_Exp_m1s[val] = {3.5378 , 8.8715 , 38.6027 , 100.4287};
Double_t Limit_Exp[val]     = {5.1514 , 13.1348, 58.0078 , 151.9531};
Double_t Limit_Exp_p1s[val] = {7.6359 , 19.9934, 90.1476 , 236.1443};
Double_t Limit_Exp_p2s[val] = {10.8662, 30.3567, 136.6394, 354.4008};
//LIMIT - GaussXLandau
Double_t Limit_Obs2[val]     = {0,0,0,0};
Double_t Limit_Exp_m2s2[val] = {2.2500 , 8.0346 , 27.1690 , 63.4689 };
Double_t Limit_Exp_m1s2[val] = {3.1133 , 11.0869, 38.2007 , 90.3546 };
Double_t Limit_Exp2[val]     = {4.6265 , 16.2598, 58.2031 , 138.2812};
Double_t Limit_Exp_p1s2[val] = {7.0792 , 24.7502, 93.0032 , 219.3069};
Double_t Limit_Exp_p2s2[val] = {10.3827, 36.8075, 148.7738, 337.1894};
*/
