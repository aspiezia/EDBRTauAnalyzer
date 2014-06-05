{
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5);  //title X location 
  gStyle->SetTitleY(0.96); //title Y location 
  gStyle->SetPaintTextFormat(".2f");

  using namespace std;

  bool save = true;
  bool QCD = true;
  bool MUON = false; bool DY=false; bool dR = false;
  char channel1[500]; char channel2[500];
  if(!DY) {sprintf(channel1,"_et"); sprintf(channel2,"_mt");}
  else    {sprintf(channel1,   ""); sprintf(channel2,   "");}
  char efficiency1[500];  sprintf(efficiency1, "Ele%s",channel1); 
  char efficiency2[500];  sprintf(efficiency2, "Muo%s",channel2); 
  char efficiency3[500];  sprintf(efficiency3, "Ele%s",channel1); 
  char efficiency4[500];  sprintf(efficiency4, "HighptMuo%s",channel2); 
  char demo       [500];  sprintf(demo,        "demo/Tree"); 
  char openTree   [500];  sprintf(openTree,    "%s",demo); 
  char demoQCD    [500];  sprintf(demoQCD,     "demo/TreeQCD"); 
  char openTreeQCDElectron[500];  sprintf(openTreeQCDElectron, "%sElectron", demoQCD); 
  char openTreeQCDMuon[500];      sprintf(openTreeQCDMuon,     "%sMuon",     demoQCD); 
  if(!dR){
    int N = 10;
    Double_t xbins[N] = {10,50,100,200,300,400,500,600,800,1000};
  } else {
    int N = 7;
    Double_t xbins[N] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8};
  }

  if(!QCD){
    if(!DY){
      TFile *file1 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/ZH2500.root");
      TFile *file2 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/ZH2500.root");
      TFile *file3 = TFile::Open("../../RISULTATI/efficiency_220514/ZH2500.root");
      TFile *file4 = TFile::Open("../../RISULTATI/efficiency_220514/ZH2500.root");
    } else {
      TFile *file1 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/DY100.root");
      TFile *file2 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/DY100.root");
      TFile *file3 = TFile::Open("../../RISULTATI/efficiency_220514/DY100.root");
      TFile *file4 = TFile::Open("../../RISULTATI/efficiency_220514/DY100.root");
    }
    TTree *Tree1; Tree1 = (TTree*)file1->Get(openTree);
    TTree *Tree2; Tree2 = (TTree*)file2->Get(openTree);
    TTree *Tree3; Tree3 = (TTree*)file3->Get(openTree);
    TTree *Tree4; Tree4 = (TTree*)file4->Get(openTree);
    
    if(!dR){
      char input1[500]; char input2[500]; char input3[500]; char input4[500]; char input5[500]; char input6[500];
      sprintf(input1,  "gen%s>0.5 && matched%s",efficiency1,efficiency3);
      sprintf(input2,  "gen%s>0.5 && matched%s",efficiency2,efficiency4);
      sprintf(input3,  "gen%s>0.5 && matched%s && reco%s_CorrPFIso<0.1",efficiency1,efficiency3,efficiency3);
      sprintf(input4,  "gen%s>0.5 && matched%s && reco%s_CorrPFIso<0.2",efficiency2,efficiency4,efficiency4);
      sprintf(input5,  "gen%s>0.5 && matched%s && reco%s_PFIso<0.1",efficiency1,efficiency3,efficiency3);
      sprintf(input6,  "gen%s>0.5 && matched%s && reco%s_PFIso<0.2",efficiency2,efficiency4,efficiency4);
      char variable1[500]; sprintf(variable1, "gen%s_Pt>>h1(4000,0,2500)",efficiency1); Tree1->Draw(variable1, input3, "E");     //ele-tau -> Corr 
      char variable2[500]; sprintf(variable2, "gen%s_Pt>>h2(4000,0,2500)",efficiency2); Tree2->Draw(variable2, input4, "Esame"); //muo-tau -> Corr 
      char variable3[500]; sprintf(variable3, "gen%s_Pt>>h3(4000,0,2500)",efficiency1); Tree3->Draw(variable3, input5, "Esame"); //ele-tau -> UnCorr 
      char variable4[500]; sprintf(variable4, "gen%s_Pt>>h4(4000,0,2500)",efficiency2); Tree4->Draw(variable4, input6, "Esame"); //muo-tau -> UnCorr 
      char variable5[500]; sprintf(variable5, "gen%s_Pt>>h5(4000,0,2500)",efficiency1); Tree1->Draw(variable5, input1, "Esame"); //ele-tau -> Corr 
      char variable6[500]; sprintf(variable6, "gen%s_Pt>>h6(4000,0,2500)",efficiency2); Tree2->Draw(variable6, input2, "Esame"); //muo-tau -> Corr 
      char variable7[500]; sprintf(variable7, "gen%s_Pt>>h7(4000,0,2500)",efficiency1); Tree3->Draw(variable7, input1, "Esame"); //ele-tau -> UnCorr 
      char variable8[500]; sprintf(variable8, "gen%s_Pt>>h8(4000,0,2500)",efficiency2); Tree4->Draw(variable8, input2, "Esame"); //muo-tau -> UnCorr
    } else {
      char input1[500]; char input2[500]; char input3[500]; char input4[500]; char input5[500]; char input6[500];
      sprintf(input1,  "gen%s>0.5 && matched%s",efficiency1,efficiency3);
      sprintf(input2,  "gen%s>0.5 && matched%s",efficiency2,efficiency4);
      sprintf(input3,  "gen%s>0.5 && matched%s && reco%s_CorrPFIso<0.1",efficiency1,efficiency3,efficiency3);
      sprintf(input4,  "gen%s>0.5 && matched%s && reco%s_CorrPFIso<0.2",efficiency2,efficiency4,efficiency4);
      sprintf(input5,  "gen%s>0.5 && matched%s && reco%s_PFIso<0.1",efficiency1,efficiency3,efficiency3);
      sprintf(input6,  "gen%s>0.5 && matched%s && reco%s_PFIso<0.2",efficiency2,efficiency4,efficiency4);
      char variable1[500]; sprintf(variable1, "gen%s_deltaR>>h1(4000,0,10)",channel1); Tree1->Draw(variable1, input3, "E");     //ele-tau -> Corr 
      char variable2[500]; sprintf(variable2, "gen%s_deltaR>>h2(4000,0,10)",channel2); Tree2->Draw(variable2, input4, "Esame"); //muo-tau -> Corr 
      char variable3[500]; sprintf(variable3, "gen%s_deltaR>>h3(4000,0,10)",channel1); Tree3->Draw(variable3, input5, "Esame"); //ele-tau -> UnCorr 
      char variable4[500]; sprintf(variable4, "gen%s_deltaR>>h4(4000,0,10)",channel2); Tree4->Draw(variable4, input6, "Esame"); //muo-tau -> UnCorr 
      char variable5[500]; sprintf(variable5, "gen%s_deltaR>>h5(4000,0,10)",channel1); Tree1->Draw(variable5, input1, "Esame"); //ele-tau -> Corr 
      char variable6[500]; sprintf(variable6, "gen%s_deltaR>>h6(4000,0,10)",channel2); Tree2->Draw(variable6, input2, "Esame"); //muo-tau -> Corr 
      char variable7[500]; sprintf(variable7, "gen%s_deltaR>>h7(4000,0,10)",channel1); Tree3->Draw(variable7, input1, "Esame"); //ele-tau -> UnCorr 
      char variable8[500]; sprintf(variable8, "gen%s_deltaR>>h8(4000,0,10)",channel2); Tree4->Draw(variable8, input2, "Esame"); //muo-tau -> UnCorr
    }
  } else {
    TFile *file1 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/QCD1000.root");
    TFile *file2 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/QCD1000.root");
    TFile *file3 = TFile::Open("../../RISULTATI/efficiency_220514/QCD1000.root");
    TFile *file4 = TFile::Open("../../RISULTATI/efficiency_220514/QCD1000.root");
    
    TTree *Tree1; Tree1 = (TTree*)file1->Get(openTreeQCDElectron);
    TTree *Tree2; Tree2 = (TTree*)file2->Get(openTreeQCDMuon);
    TTree *Tree3; Tree3 = (TTree*)file3->Get(openTreeQCDElectron);
    TTree *Tree4; Tree4 = (TTree*)file4->Get(openTreeQCDMuon);
    
    char input1[500]; char input2[500]; char input3[500]; char input4[500]; char input5[500]; char input6[500];
    sprintf(input1,  "JetElectron_Pt>10", efficiency1);
    sprintf(input2,  "JetMuon_Pt>10",     efficiency2);
    sprintf(input3,  "JetElectron_Pt>10 && JetElectron_CorrPFIso<0.1",efficiency1);
    sprintf(input4,  "JetMuon_Pt>10     && JetMuon_CorrPFIso<0.2",    efficiency2);
    sprintf(input5,  "JetElectron_Pt>10 && JetElectron_PFIso<0.1",    efficiency1);
    sprintf(input6,  "JetMuon_Pt>10     && JetMuon_PFIso<0.2",        efficiency2);
    
    char variable1[500]; sprintf(variable1, "JetElectron_Pt>>h1(4000,0,2500)",efficiency1); Tree1->Draw(variable1, input3, "E");     //ele-tau -> Corr 
    char variable2[500]; sprintf(variable2, "JetMuon_Pt>>h2(4000,0,2500)",    efficiency2); Tree2->Draw(variable2, input4, "Esame"); //muo-tau -> Corr 
    char variable3[500]; sprintf(variable3, "JetElectron_Pt>>h3(4000,0,2500)",efficiency1); Tree3->Draw(variable3, input5, "Esame"); //ele-tau -> UnCorr 
    char variable4[500]; sprintf(variable4, "JetMuon_Pt>>h4(4000,0,2500)",    efficiency2); Tree4->Draw(variable4, input6, "Esame"); //muo-tau -> UnCorr 
    char variable5[500]; sprintf(variable5, "JetElectron_Pt>>h5(4000,0,2500)",efficiency1); Tree1->Draw(variable5, input1, "Esame"); //ele-tau -> Corr 
    char variable6[500]; sprintf(variable6, "JetMuon_Pt>>h6(4000,0,2500)",    efficiency2); Tree2->Draw(variable6, input2, "Esame"); //muo-tau -> Corr 
    char variable7[500]; sprintf(variable7, "JetElectron_Pt>>h7(4000,0,2500)",efficiency1); Tree3->Draw(variable7, input1, "Esame"); //ele-tau -> UnCorr 
    char variable8[500]; sprintf(variable8, "JetMuon_Pt>>h8(4000,0,2500)",    efficiency2); Tree4->Draw(variable8, input2, "Esame"); //muo-tau -> UnCorr
  }

  Hnew1 = h1->Rebin(N-1,"Hnew1",xbins);
  Hnew2 = h2->Rebin(N-1,"Hnew2",xbins);
  Hnew3 = h3->Rebin(N-1,"Hnew3",xbins);
  Hnew4 = h4->Rebin(N-1,"Hnew4",xbins);
  Hnew5 = h5->Rebin(N-1,"Hnew5",xbins);
  Hnew6 = h6->Rebin(N-1,"Hnew6",xbins);
  Hnew7 = h7->Rebin(N-1,"Hnew7",xbins);
  Hnew8 = h8->Rebin(N-1,"Hnew8",xbins);

  TGraphAsymmErrors *hnew1 = new TGraphAsymmErrors(Hnew1,Hnew5); //ELE - CORR 
  TGraphAsymmErrors *hnew2 = new TGraphAsymmErrors(Hnew2,Hnew6); //MUO - CORR
  TGraphAsymmErrors *hnew3 = new TGraphAsymmErrors(Hnew3,Hnew7); //ELE - UNCO
  TGraphAsymmErrors *hnew4 = new TGraphAsymmErrors(Hnew4,Hnew8); //MUO - UNCO

  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);

  //hnew1->Divide(hnew5);
  //hnew2->Divide(hnew6);
  //hnew3->Divide(hnew7);
  //hnew4->Divide(hnew8);
  
  hnew1->SetLineWidth(2);
  hnew2->SetLineWidth(2);
  hnew3->SetLineWidth(2);
  hnew4->SetLineWidth(2);
  
  hnew1->SetLineColor(2);
  hnew2->SetLineColor(2);
  hnew3->SetLineColor(1);
  hnew4->SetLineColor(1);
  
  hnew1->SetMarkerColor(2);
  hnew2->SetMarkerColor(2);
  hnew3->SetMarkerColor(1);
  hnew4->SetMarkerColor(1);

  hnew1->SetMarkerStyle(20); 
  hnew2->SetMarkerStyle(20); 
  hnew3->SetMarkerStyle(20); 
  hnew4->SetMarkerStyle(20); 
  
  hnew1->SetMarkerSize(1.3); 
  hnew2->SetMarkerSize(1.3); 
  hnew3->SetMarkerSize(1.3); 
  hnew4->SetMarkerSize(1.3); 
  
  if(!MUON){
    hnew1->Draw("AP");
    hnew3->Draw("Psame");
  } else {
    hnew2->Draw("AP");
    hnew4->Draw("Psame");
  }

  hnew1->SetTitle("");
  hnew1->GetYaxis()->SetTitle("Efficiency");
  hnew1->GetYaxis()->SetTitleOffset(1.2);
  hnew1->SetMinimum(0.00);
  hnew2->SetTitle("");
  hnew2->GetYaxis()->SetTitle("Efficiency");
  hnew2->GetYaxis()->SetTitleOffset(1.2);
  hnew2->SetMinimum(0.00);
  if(!QCD){
    hnew1->SetMaximum(1.02);
    hnew1->GetXaxis()->SetTitle("gen lepton pt [GeV]");
    if(dR) hnew1->GetXaxis()->SetTitle("#DeltaR(gen lepton, gen tau)");
    hnew2->SetMaximum(1.02);
    hnew2->GetXaxis()->SetTitle("gen lepton pt [GeV]");
    if(dR) hnew2->GetXaxis()->SetTitle("#DeltaR(gen lepton, gen tau)");
  } else {
      hnew1->GetXaxis()->SetTitle("lepton pt [GeV]");
      hnew2->GetXaxis()->SetTitle("lepton pt [GeV]");
      if(MUON){
	hnew1->SetMaximum(0.40);
	hnew2->SetMaximum(0.40);
      } else {
	hnew1->SetMaximum(0.12);
	hnew2->SetMaximum(0.12);
      }
  }

  if(!dR) hnew1->GetXaxis()->SetRangeUser(0,1000);
  else    hnew1->GetXaxis()->SetRangeUser(0.2,0.7);
  if(!dR) hnew2->GetXaxis()->SetRangeUser(0,1000);
  else    hnew2->GetXaxis()->SetRangeUser(0.2,0.7);

  TLegend *pl2 = new TLegend(0.45,0.12,0.89,0.29);
  pl2->SetTextSize(0.03); 
  pl2->SetFillColor(0);
  if(!MUON){
    TLegendEntry *ple2 = pl2->AddEntry(hnew3, "ele-tau: Standard Isolation",  "LP");
    ple2 = pl2->AddEntry(hnew1, "ele-tau: Modified Isolation",  "LP");
  } else {
    TLegendEntry *ple2 = pl2->AddEntry(hnew4, "muo-tau: Standard Isolation",  "LP");
    ple2 = pl2->AddEntry(hnew2, "muo-tau: Modified Isolation",  "LP");
  }
  pl2->Draw();
  
  if(save){
    if(!DY){
      if(QCD  && MUON)         c1->SaveAs("eff_isolation_muon_QCD.pdf");
      if(QCD  && !MUON)        c1->SaveAs("eff_isolation_electron_QCD.pdf");
      if(!QCD && MUON  && !dR) c1->SaveAs("eff_isolation_muon_pt.pdf");
      if(!QCD && !MUON && !dR) c1->SaveAs("eff_isolation_electron_pt.pdf");
      if(!QCD && MUON  && dR)  c1->SaveAs("eff_isolation_muon_dR.pdf");
      if(!QCD && !MUON && dR)  c1->SaveAs("eff_isolation_electron_dR.pdf");
    } else {
      if(QCD  && MUON)         c1->SaveAs("eff_DY_isolation_muon_QCD.pdf");
      if(QCD  && !MUON)        c1->SaveAs("eff_DY_isolation_electron_QCD.pdf");
      if(!QCD && MUON  && !dR) c1->SaveAs("eff_DY_isolation_muon_pt.pdf");
      if(!QCD && !MUON && !dR) c1->SaveAs("eff_DY_isolation_electron_pt.pdf");
      if(!QCD && MUON  && dR)  c1->SaveAs("eff_DY_isolation_muon_dR.pdf");
      if(!QCD && !MUON && dR)  c1->SaveAs("eff_DY_isolation_electron_dR.pdf");
    }
  }
}
