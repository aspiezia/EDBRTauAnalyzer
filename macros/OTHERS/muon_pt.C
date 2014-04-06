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

  bool save = false;
  char efficiency[500];      sprintf(efficiency,     "Muo_mt"); 
  char efficiencyReco[500];  sprintf(efficiencyReco, "Highpt%s", efficiency); 
  char demo      [500];      sprintf(demo,           "demo/Tree"); 
  char openTree  [500];      sprintf(openTree,       "%s",demo); 
  TString EFFICIENCY = efficiency;
  int N = 5;

  //Double_t xbins[N] = {400, 500, 600, 700, 800, 900, 1000, 1500};
  Double_t xbins[N] = {0,10,200,400,1000};

  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);

  TFile *file1 = TFile::Open("../RISULTATI/MuonSelection/ZH1000_Eff.root");
  TFile *file2 = TFile::Open("../RISULTATI/MuonSelection/ZH1500_Eff.root");
  TFile *file3 = TFile::Open("../RISULTATI/MuonSelection/ZH2000_Eff.root");
  TFile *file4 = TFile::Open("../RISULTATI/MuonSelection/ZH2500_Eff.root");

  TTree *Tree1;  
  TTree *Tree2;
  TTree *Tree3;
  TTree *Tree4;
  
  Tree1 = (TTree*)file1->Get(openTree); 
  Tree2 = (TTree*)file2->Get(openTree); 
  Tree3 = (TTree*)file3->Get(openTree); 
  Tree4 = (TTree*)file4->Get(openTree); 
  
  char input1[500]; char input2[500]; 
  sprintf(input1,  "gen%s>0.5 && matched%s",efficiency,efficiencyReco);
  sprintf(input2,  "gen%s>0.5 && matched%s && reco%s_PFIso<0.2",efficiency,efficiencyReco,efficiencyReco);
  
  char variable1[500]; sprintf(variable1, "gen%s_Pt>>h1(4000,0,2500)",efficiency); Tree1->Draw(variable1, input2, "E");
  char variable2[500]; sprintf(variable2, "gen%s_Pt>>h2(4000,0,2500)",efficiency); Tree2->Draw(variable2, input2, "Esame");
  char variable3[500]; sprintf(variable3, "gen%s_Pt>>h3(4000,0,2500)",efficiency); Tree3->Draw(variable3, input2, "Esame");
  char variable4[500]; sprintf(variable4, "gen%s_Pt>>h4(4000,0,2500)",efficiency); Tree4->Draw(variable4, input2, "Esame");
  char variable5[500]; sprintf(variable5, "gen%s_Pt>>h5(4000,0,2500)",efficiency); Tree1->Draw(variable5, input1, "Esame");
  char variable6[500]; sprintf(variable6, "gen%s_Pt>>h6(4000,0,2500)",efficiency); Tree2->Draw(variable6, input1, "Esame");
  char variable7[500]; sprintf(variable7, "gen%s_Pt>>h7(4000,0,2500)",efficiency); Tree3->Draw(variable7, input1, "Esame");
  char variable8[500]; sprintf(variable8, "gen%s_Pt>>h8(4000,0,2500)",efficiency); Tree4->Draw(variable8, input1, "Esame");
  
  hnew1 = h1->Rebin(N-1,"hnew1",xbins);
  hnew2 = h2->Rebin(N-1,"hnew2",xbins);
  hnew3 = h3->Rebin(N-1,"hnew3",xbins);
  hnew4 = h4->Rebin(N-1,"hnew4",xbins);
  hnew5 = h5->Rebin(N-1,"hnew5",xbins);
  hnew6 = h6->Rebin(N-1,"hnew6",xbins);
  hnew7 = h7->Rebin(N-1,"hnew7",xbins);
  hnew8 = h8->Rebin(N-1,"hnew8",xbins);

  hnew1->Divide(hnew5);
  hnew2->Divide(hnew6);
  hnew3->Divide(hnew7);
  hnew4->Divide(hnew8);
  
  hnew1->SetLineWidth(2);
  hnew2->SetLineWidth(2);
  hnew3->SetLineWidth(2);
  hnew4->SetLineWidth(2);
  
  hnew1->SetLineColor(1);
  hnew2->SetLineColor(2);
  hnew3->SetLineColor(kGreen+2);
  hnew4->SetLineColor(4);
  
  hnew1->SetMarkerColor(1);
  hnew2->SetMarkerColor(2);
  hnew3->SetMarkerColor(kGreen+2);
  hnew4->SetMarkerColor(4);

  hnew1->SetMarkerStyle(20); 
  hnew2->SetMarkerStyle(20); 
  hnew3->SetMarkerStyle(20); 
  hnew4->SetMarkerStyle(20); 
  
  hnew1->SetMarkerSize(1.3); 
  hnew2->SetMarkerSize(1.3); 
  hnew3->SetMarkerSize(1.3); 
  hnew4->SetMarkerSize(1.3); 
  
  hnew1->Draw("E");
  hnew2->Draw("Esame");
  hnew3->Draw("Esame");
  hnew4->Draw("Esame");

  hnew1->SetTitle("");
  hnew1->GetYaxis()->SetTitle("Efficiency");
  hnew1->SetMaximum(1.50);
  hnew1->SetMinimum(0.00);
  if(EFFICIENCY=="Jet")      hnew1->GetXaxis()->SetTitle("gen jet pt [GeV]");
  if(EFFICIENCY=="Tau_mt")   hnew1->GetXaxis()->SetTitle("gen tau pt [GeV]");
  if(EFFICIENCY=="Tau_et")   hnew1->GetXaxis()->SetTitle("gen tau pt [GeV]");
  if(EFFICIENCY=="Ele_et")   hnew1->GetXaxis()->SetTitle("gen electron pt [GeV]");
  if(EFFICIENCY=="Ele1_ee")  hnew1->GetXaxis()->SetTitle("gen electron pt [GeV]");
  if(EFFICIENCY=="Ele2_ee")  hnew1->GetXaxis()->SetTitle("gen electron pt [GeV]");
  if(EFFICIENCY=="Ele_em")   hnew1->GetXaxis()->SetTitle("gen electron pt [GeV]");
  if(EFFICIENCY=="Muo_mt")   hnew1->GetXaxis()->SetTitle("gen muon pt [GeV]");
  if(EFFICIENCY=="Muo1_mm")  hnew1->GetXaxis()->SetTitle("gen muon pt [GeV]");
  if(EFFICIENCY=="Muo2_mm")  hnew1->GetXaxis()->SetTitle("gen muon pt [GeV]");
  if(EFFICIENCY=="Muo_em")   hnew1->GetXaxis()->SetTitle("gen muon pt [GeV]");

  TLegend *pl2 = new TLegend(0.68,0.72,0.89,0.89);
  pl2->SetTextSize(0.03); 
  pl2->SetFillColor(0);
  TLegendEntry *ple2 = pl2->AddEntry(hnew1, "M(X) = 1.0 TeV",  "LP");
  ple2 = pl2->AddEntry(hnew2, "M(X) = 1.5 TeV",  "LP");
  ple2 = pl2->AddEntry(hnew3, "M(X) = 2.0 TeV",  "LP");
  ple2 = pl2->AddEntry(hnew4, "M(X) = 2.5 TeV",  "LP");
  pl2->Draw();
  
  if(save){
    c1->SaveAs("eff_"+EFFICIENCY+".png");
    cout<<"Saving eff_"<<EFFICIENCY<<".png"<<endl;
  }
}
