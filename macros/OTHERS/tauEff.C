{
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5); //title X location 
  gStyle->SetTitleY(0.96); //title Y location 
  gStyle->SetPaintTextFormat(".2f");

  using namespace std;

  bool QCD = false;
  bool save = false;
  char channel   [500];   sprintf(channel,   "MuoTau"); 
  char demo      [500];   sprintf(demo,      "demo/TreeEff"); 
  char demoUsual [500];   sprintf(demoUsual, "demo/TreeEffUsual"); 
  char openTree  [500];   sprintf(openTree,  "%s%s",demo,channel); 
  char openTreeU [500];   sprintf(openTreeU, "%s%s",demoUsual,channel); 
  TString CHANNEL = channel;

  TH1D* tauSelection1 = new TH1D("tauSelection1","tauSelection1",9,-0.5,8.5);
  TH1D* tauSelection2 = new TH1D("tauSelection2","tauSelection2",9,-0.5,8.5);
  TH1D* tauSelection3 = new TH1D("tauSelection3","tauSelection3",9,-0.5,8.5);
  TH1D* tauSelection4 = new TH1D("tauSelection4","tauSelection4",9,-0.5,8.5);
  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
  TFile *file1 = TFile::Open("../RISULTATI/AnalysisNotePlots_070314/ZH2500_tauMC.root");

  TTree *Tree1;  Tree1 = (TTree*)file1->Get(openTree);  
  TTree *Tree2;	 Tree2 = (TTree*)file1->Get(openTreeU);

  int N1=0;
  int N2=0;

  char input1[500]; sprintf(input1,  "uno==1");
  char input2[500]; sprintf(input2,  "due==1");
  char input3[500]; sprintf(input3,  "tre==1");
  char input4[500]; sprintf(input4,  "qua==1");
  char input5[500]; sprintf(input5,  "cin==1");
  char input6[500]; sprintf(input6,  "sei==1");
  char input7[500]; sprintf(input7,  "set==1");

  Tree1->Draw("uno>>h1(4,-0.5,3.5)", input1, "E"); 
  N1=h1->GetBinContent(2); 
  tauSelection1->SetBinContent(2,N1); 
  tauSelection1->SetBinError(2,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree1->Draw("due>>h1(4,-0.5,3.5)", input2, "E"); 
  tauSelection1->SetBinContent(3,h1->GetBinContent(2));	       
  tauSelection1->SetBinError(3,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree1->Draw("tre>>h1(4,-0.5,3.5)", input3, "E"); 
  tauSelection1->SetBinContent(4,h1->GetBinContent(2));	       
  tauSelection1->SetBinError(4,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree1->Draw("qua>>h1(4,-0.5,3.5)", input4, "E"); 
  tauSelection1->SetBinContent(5,h1->GetBinContent(2));	       
  tauSelection1->SetBinError(5,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree1->Draw("cin>>h1(4,-0.5,3.5)", input5, "E"); 
  tauSelection1->SetBinContent(6,h1->GetBinContent(2));	       
  tauSelection1->SetBinError(6,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree1->Draw("sei>>h1(4,-0.5,3.5)", input6, "E"); 
  tauSelection1->SetBinContent(7,h1->GetBinContent(2));	       
  tauSelection1->SetBinError(7,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree1->Draw("set>>h1(4,-0.5,3.5)", input7, "E"); 
  tauSelection1->SetBinContent(8,h1->GetBinContent(2));        
  tauSelection1->SetBinError(8,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));

  Tree2->Draw("uno>>h1(4,-0.5,3.5)", input1, "E"); 
  N1=h1->GetBinContent(2); 
  tauSelection2->SetBinContent(2,N1); 
  tauSelection2->SetBinError(2,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree2->Draw("due>>h1(4,-0.5,3.5)", input2, "E"); 
  tauSelection2->SetBinContent(3,h1->GetBinContent(2));	       
  tauSelection2->SetBinError(3,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree2->Draw("tre>>h1(4,-0.5,3.5)", input3, "E"); 
  tauSelection2->SetBinContent(4,h1->GetBinContent(2));	       
  tauSelection2->SetBinError(4,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree2->Draw("qua>>h1(4,-0.5,3.5)", input4, "E"); 
  tauSelection2->SetBinContent(5,h1->GetBinContent(2));	       
  tauSelection2->SetBinError(5,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree2->Draw("cin>>h1(4,-0.5,3.5)", input5, "E"); 
  tauSelection2->SetBinContent(6,h1->GetBinContent(2));	       
  tauSelection2->SetBinError(6,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree2->Draw("sei>>h1(4,-0.5,3.5)", input6, "E"); 
  tauSelection2->SetBinContent(7,h1->GetBinContent(2));	       
  tauSelection2->SetBinError(7,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));
  Tree2->Draw("set>>h1(4,-0.5,3.5)", input7, "E"); 
  tauSelection2->SetBinContent(8,h1->GetBinContent(2));        
  tauSelection2->SetBinError(8,(1./(N1*N1))*sqrt(N1*h1->GetBinContent(2)*(N1-h1->GetBinContent(2))));


  float a1 = 100./tauSelection1->GetBinContent(2);
  tauSelection1->Scale(a1);
  tauSelection1->SetLineWidth(2);
  tauSelection1->SetLineColor(1);
  tauSelection1->SetMarkerColor(1);
  tauSelection1->SetMarkerStyle(20); 
  tauSelection1->SetMarkerSize(1.3); 
  tauSelection1->SetTitle("");
  tauSelection1->GetYaxis()->SetTitle("Tau selection efficiency [%]");
  tauSelection1->SetMaximum(120);
  tauSelection1->SetMinimum(0.00);
  tauSelection1->GetYaxis()->SetTitleSize(0.045);
  tauSelection1->GetXaxis()->SetTitleSize(0.045);
  tauSelection1->GetXaxis()->SetLabelSize(0.040);
  tauSelection1->GetYaxis()->SetLabelSize(0.045);
  tauSelection1->GetYaxis()->SetTitleOffset(1.05);
  tauSelection1->GetXaxis()->SetRangeUser(0,8);
  tauSelection1->GetXaxis()->SetNdivisions(0015);
  tauSelection1->Draw("E");

  float a2 = 100./tauSelection2->GetBinContent(2);
  tauSelection2->Scale(a2);
  tauSelection2->SetLineWidth(2);
  tauSelection2->SetLineColor(2);
  tauSelection2->SetMarkerColor(2);
  tauSelection2->SetMarkerStyle(20); 
  tauSelection2->SetMarkerSize(1.3); 
  tauSelection2->Draw("E same");

  for(int i=0; i<10; i++){
    tauSelection3->SetBinContent(i,-1000);
    tauSelection4->SetBinContent(i,-1000);
  }
  tauSelection3->SetBinContent(8,tauSelection1->GetBinContent(8));        
  tauSelection3->SetBinError(8,tauSelection1->GetBinError(8));
  tauSelection4->SetBinContent(8,tauSelection2->GetBinContent(8));        
  tauSelection4->SetBinError(8,tauSelection2->GetBinError(8));
  tauSelection3->SetLineWidth(2);
  tauSelection3->SetLineColor(1);
  tauSelection3->SetMarkerColor(1);
  tauSelection3->SetMarkerStyle(1); 
  tauSelection3->SetMarkerSize(1.8); 
  tauSelection4->SetLineWidth(2);
  tauSelection4->SetLineColor(2);
  tauSelection4->SetMarkerColor(2);
  tauSelection4->SetMarkerStyle(1); 
  tauSelection4->SetMarkerSize(1.8);  
  tauSelection3->Draw("E same text1"); 
  tauSelection4->Draw("E same text1");

  TLegend *pl = new TLegend(0.44,0.77,0.89,0.89);
  pl->SetTextSize(0.040); 
  pl->SetFillColor(0);
  TLegendEntry *ple = pl->AddEntry(tauSelection1, "Cleaned PFJet collection",  "L");
  ple = pl->AddEntry(tauSelection2, "Standard PFJet collection",  "L");
  pl->Draw();

  //if(save)
}
