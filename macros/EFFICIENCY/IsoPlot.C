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
  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);

  bool save = true;
  bool MUON=true; bool QCD=true;
  char openTree[500];             sprintf(openTree,    "demo/Tree"); 
  char demoQCD[500];              sprintf(demoQCD,     "demo/TreeQCD"); 
  char openTreeQCDElectron[500];  sprintf(openTreeQCDElectron, "%sElectron", demoQCD); 
  char openTreeQCDMuon[500];      sprintf(openTreeQCDMuon,     "%sMuon",     demoQCD); 
  
  if(!QCD){
    TFile *file1 = TFile::Open("../../RISULTATI/efficiency_220514/ZH2500.root"); 
    TFile *file2 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/ZH2500.root"); 
    TTree *Tree1; Tree1 = (TTree*)file1->Get(openTree); 
    TTree *Tree2; Tree2 = (TTree*)file2->Get(openTree); 
  } else {
    TFile *file1 = TFile::Open("../../RISULTATI/efficiency_220514/QCD1000.root"); 
    TFile *file2 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/QCD1000.root"); 
    TTree *Tree1; Tree1 = (TTree*)file1->Get(openTreeQCDElectron);
    TTree *Tree2; Tree2 = (TTree*)file1->Get(openTreeQCDMuon);
    TTree *Tree3; Tree3 = (TTree*)file2->Get(openTreeQCDElectron);
    TTree *Tree4; Tree4 = (TTree*)file2->Get(openTreeQCDMuon);
  }

  if(!QCD){
    if(MUON){
      Tree1->Draw("recoHighptMuo_mt_PFIso>>h1(100,0,3)","genMuo_mt>0.5 && matchedHighptMuo_mt","");
      Tree2->Draw("recoHighptMuo_mt_CorrPFIso>>h2(100,0,3)","genMuo_mt>0.5 && matchedHighptMuo_mt","same");
    } else {
      Tree1->Draw("recoEle_et_PFIso>>h1(100,0,3)","genEle_et>0.5 && matchedEle_et","");
      Tree2->Draw("recoEle_et_CorrPFIso>>h2(100,0,3)","genEle_et>0.5 && matchedEle_et","same");
    }
  } else {
    if(MUON){
      Tree2->Draw("JetMuon_PFIso>>h1(100,0,10)","","");
      Tree4->Draw("JetMuon_CorrPFIso>>h2(100,0,10)","","same");
    } else {
      Tree1->Draw("JetElectron_PFIso>>h1(100,0,10)","","");
      Tree3->Draw("JetElectron_CorrPFIso>>h2(100,0,10)","","same");
    }
  }

  h1->SetTitle();
  h1->SetLineWidth(2);
  h1->SetLineColor(1);
  h1->SetXTitle("PFIso/pt");
  h1->SetYTitle("Events");

  h2->SetTitle();
  h2->SetLineWidth(2);
  h2->SetLineColor(2);
  h2->SetXTitle("PFIso/pt");
  h2->SetYTitle("Events");

  h2->Draw();
  h1->Draw("same");


  TLegend *pl2 = new TLegend(0.50,0.75,0.89,0.89);
  pl2->SetTextSize(0.03); 
  pl2->SetFillColor(0);
  if(!MUON){
    TLegendEntry *ple2 = pl2->AddEntry(h1, "ele-tau: Standard Isolation",  "LP");
    ple2 = pl2->AddEntry(h2, "ele-tau: Modified Isolation",  "LP");
  } else {
    TLegendEntry *ple2 = pl2->AddEntry(h1, "muo-tau: Standard Isolation",  "LP");
    ple2 = pl2->AddEntry(h2, "muo-tau: Modified Isolation",  "LP");
  }
  pl2->Draw();
  c1->SetLogy();

  if(save){
    if(!QCD){
      if(MUON) c1->SaveAs("PFIsoCorr_muon.pdf");
      else     c1->SaveAs("PFIsoCorr_electron.pdf");
    } else {
      if(MUON) c1->SaveAs("PFIsoCorr_muon_QCD.pdf");
      else     c1->SaveAs("PFIsoCorr_electron_QCD.pdf");
    }
  }
  
}
