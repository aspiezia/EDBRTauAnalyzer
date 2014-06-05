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
  bool save=false;
  
  int N = 7;
  Double_t xbins[N] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8};

  int i=1;
  char channel[500];
  if(i==0) sprintf(channel,"EleMuo"); 
  if(i==1) sprintf(channel,"MuoMuo"); 
  if(i==2) sprintf(channel,"EleEle"); 
  if(i==3) sprintf(channel,"MuoTau"); 
  if(i==4) sprintf(channel,"EleTau"); 
  char demo      [500];  sprintf(demo,      "demo/Tree"); 
  char openTree  [500];  sprintf(openTree,  "%s%s",demo,channel); 
  TString CHANNEL = channel;

  TFile *file1 = TFile::Open("../../RISULTATI/efficiency_280514_CorrPFIso/ZH1000.root");
  TFile *file2 = TFile::Open("../../RISULTATI/efficiency_280514_CorrPFIso/ZH1500.root");
  TFile *file3 = TFile::Open("../../RISULTATI/efficiency_280514_CorrPFIso/ZH2000.root");
  TFile *file4 = TFile::Open("../../RISULTATI/efficiency_280514_CorrPFIso/ZH2500.root");
  
  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
  TTree *Tree1;  Tree1 = (TTree*)file1->Get(openTree); 
  TTree *Tree2;  Tree2 = (TTree*)file2->Get(openTree); 
  TTree *Tree3;  Tree3 = (TTree*)file3->Get(openTree); 
  TTree *Tree4;  Tree4 = (TTree*)file4->Get(openTree);

  vector<TString> NAME;
  NAME.push_back("recoLep1Pt");
  NAME.push_back("recoLep1Eta");
  NAME.push_back("recoLep1Phi");
  NAME.push_back("recoLep2Pt");
  NAME.push_back("recoLep2Eta");
  NAME.push_back("recoLep2Phi");
  NAME.push_back("recoLep1dRJet");
  NAME.push_back("recoLep2dRJet");
  if(i==0){
    //ELE
    NAME.push_back("recoLep1HEEP");
    NAME.push_back("recoLep1PFIso");
    NAME.push_back("recoLep1CorrPFIso");
    NAME.push_back("recoLep1DETIso");
    NAME.push_back("recoLep1HEEPIso");
    NAME.push_back("recoLep1dEtaIn");
    NAME.push_back("recoLep1dPhiIn");
    NAME.push_back("recoLep1SigmaIEIE");
    NAME.push_back("recoLep1HoverE");
    NAME.push_back("recoLep1dxy");
    NAME.push_back("recoLep1dz");
    NAME.push_back("recoLep1EandP");
    NAME.push_back("recoLep1Conv1");
    NAME.push_back("recoLep1Conv2");
    //MUO
    NAME.push_back("recoLep2PFIso");
    NAME.push_back("recoLep2CorrPFIso");
    NAME.push_back("recoLep2DETIso");
    NAME.push_back("recoLep2IsGlobal");
    NAME.push_back("recoLep2ptErr");
    NAME.push_back("recoLep2MuonHits");
    NAME.push_back("recoLep2Matches");
    NAME.push_back("recoLep2dxy");
    NAME.push_back("recoLep2dz");
    NAME.push_back("recoLep2PixelHits");
    NAME.push_back("recoLep2TrackerL");
    NAME.push_back("recoLep2IsPFMuon");
    NAME.push_back("recoLep2Chi2");
  } else if(i==1){
    //MUO1
    NAME.push_back("recoLep1PFIso");
    NAME.push_back("recoLep1CorrPFIso");
    NAME.push_back("recoLep1DETIso");
    NAME.push_back("recoLep1IsGlobal");
    NAME.push_back("recoLep1ptErr");
    NAME.push_back("recoLep1MuonHits");
    NAME.push_back("recoLep1Matches");
    NAME.push_back("recoLep1dxy");
    NAME.push_back("recoLep1dz");
    NAME.push_back("recoLep1PixelHits");
    NAME.push_back("recoLep1TrackerL");
    NAME.push_back("recoLep1IsPFMuon");
    NAME.push_back("recoLep1Chi2");
    //MUO2
    NAME.push_back("recoLep2PFIso");
    NAME.push_back("recoLep2CorrPFIso");
    NAME.push_back("recoLep2DETIso");
    NAME.push_back("recoLep2IsGlobal");
    NAME.push_back("recoLep2ptErr");
    NAME.push_back("recoLep2MuonHits");
    NAME.push_back("recoLep2Matches");
    NAME.push_back("recoLep2dxy");
    NAME.push_back("recoLep2dz");
    NAME.push_back("recoLep2PixelHits");
    NAME.push_back("recoLep2TrackerL");
    NAME.push_back("recoLep2IsPFMuon");
    NAME.push_back("recoLep2Chi2");
  } else if(i==2){
    //ELE1
    NAME.push_back("recoLep1HEEP");
    NAME.push_back("recoLep1PFIso");
    NAME.push_back("recoLep1CorrPFIso");
    NAME.push_back("recoLep1DETIso");
    NAME.push_back("recoLep1HEEPIso");
    NAME.push_back("recoLep1dEtaIn");
    NAME.push_back("recoLep1dPhiIn");
    NAME.push_back("recoLep1SigmaIEIE");
    NAME.push_back("recoLep1HoverE");
    NAME.push_back("recoLep1dxy");
    NAME.push_back("recoLep1dz");
    NAME.push_back("recoLep1EandP");
    NAME.push_back("recoLep1Conv1");
    NAME.push_back("recoLep1Conv2");
    //ELE2
    NAME.push_back("recoLep2HEEP");
    NAME.push_back("recoLep2PFIso");
    NAME.push_back("recoLep2CorrPFIso");
    NAME.push_back("recoLep2DETIso");
    NAME.push_back("recoLep2HEEPIso");
    NAME.push_back("recoLep2dEtaIn");
    NAME.push_back("recoLep2dPhiIn");
    NAME.push_back("recoLep2SigmaIEIE");
    NAME.push_back("recoLep2HoverE");
    NAME.push_back("recoLep2dxy");
    NAME.push_back("recoLep2dz");
    NAME.push_back("recoLep2EandP");
    NAME.push_back("recoLep2Conv1");
    NAME.push_back("recoLep2Conv2");
  } else if(i==3){
    //TAU
    NAME.push_back("recoLep1Discr1");
    NAME.push_back("recoLep1Discr2");
    NAME.push_back("recoLep1Discr3");
    NAME.push_back("recoLep1Discr4");
    //MUO
    NAME.push_back("recoLep2PFIso");
    NAME.push_back("recoLep2CorrPFIso");
    NAME.push_back("recoLep2DETIso");
    NAME.push_back("recoLep2IsGlobal");
    NAME.push_back("recoLep2ptErr");
    NAME.push_back("recoLep2MuonHits");
    NAME.push_back("recoLep2Matches");
    NAME.push_back("recoLep2dxy");
    NAME.push_back("recoLep2dz");
    NAME.push_back("recoLep2PixelHits");
    NAME.push_back("recoLep2TrackerL");
    NAME.push_back("recoLep2IsPFMuon");
    NAME.push_back("recoLep2Chi2");
  } else if(i==4){
    //TAU
    NAME.push_back("recoLep1Discr1");
    NAME.push_back("recoLep1Discr2");
    NAME.push_back("recoLep1Discr3");
    NAME.push_back("recoLep1Discr4");
    //ELE
    NAME.push_back("recoLep2HEEP");
    NAME.push_back("recoLep2PFIso");
    NAME.push_back("recoLep2CorrPFIso");
    NAME.push_back("recoLep2DETIso");
    NAME.push_back("recoLep2HEEPIso");
    NAME.push_back("recoLep2dEtaIn");
    NAME.push_back("recoLep2dPhiIn");
    NAME.push_back("recoLep2SigmaIEIE");
    NAME.push_back("recoLep2HoverE");
    NAME.push_back("recoLep2dxy");
    NAME.push_back("recoLep2dz");
    NAME.push_back("recoLep2EandP");
    NAME.push_back("recoLep2Conv1");
    NAME.push_back("recoLep2Conv2");
  }
  NAME.push_back("AllLep1");
  NAME.push_back("AllLep2");
  NAME.push_back("AllLep1PlusIso");
  NAME.push_back("AllLep2PlusIso");
  NAME.push_back("AllLep1HEEPPlusIso");
  NAME.push_back("AllLep2HEEPPlusIso");

  for(int j=0; j<NAME.size(); j++){
    char *name; name=NAME[j].Data();
    
    char input1[500]; char input2[500];
    if(NAME[j]=="AllLep1" && (i==0 || i==2)){
      sprintf(input1,  "recoLep1PtF>0");
      sprintf(input2,  "recoLep1Pt==1 && recoLep1Eta==1 && recoLep1Phi==1 && recoLep1dEtaIn==1 && recoLep1dPhiIn==1 && recoLep1SigmaIEIE==1 && recoLep1HoverE==1 && recoLep1dxy==1 && recoLep1dz==1 && recoLep1EandP==1 && recoLep1Conv1==1 && recoLep1Conv2==1");
    } else if(NAME[j]=="AllLep1" && i==1){
      sprintf(input1,  "recoLep1PtF>0");
      sprintf(input2,  "recoLep1Pt==1 && recoLep1Eta==1 && recoLep1Phi==1 && recoLep1Matches==1 && recoLep1dxy==1 && recoLep1dz==1 && recoLep1PixelHits==1 && recoLep1TrackerL==1&& recoLep1ptErr==1 ");
    } else  if(NAME[j]=="AllLep1" && (i==3 || i==4)){
      sprintf(input1,  "recoLep1PtF>0");
      sprintf(input2,  "recoLep1Pt==1 && recoLep1Eta==1 && recoLep1Phi==1 && recoLep1Discr1==1 && recoLep1Discr2==1 && recoLep1Discr3==1 && recoLep1Discr4==1");
    } else if(NAME[j]=="AllLep2" && i==1){
      sprintf(input1,  "recoLep2PtF>0");
      sprintf(input2,  "recoLep2Pt==1 && recoLep2Eta==1 && recoLep2Phi==1 && recoLep2Matches==1 && recoLep2dxy==1 && recoLep2dz==1 && recoLep2PixelHits==1 && recoLep2TrackerL==1&& recoLep2ptErr==1 ");
    } else if(NAME[j]=="AllLep2" && (i==0 ||i==3)){
      sprintf(input1,  "recoLep2PtF>0");
      sprintf(input2,  "recoLep2Pt==1 && recoLep2Eta==1 && recoLep2Phi==1 && recoLep2IsGlobal==1 && recoLep2ptErr==1 && recoLep2MuonHits==1 && recoLep2Matches==1 && recoLep2dxy==1 && recoLep2dz==1 && recoLep2PixelHits==1 && recoLep2TrackerL==1");
    } else if(NAME[j]=="AllLep2" && (i==2 || i==4)){
      sprintf(input1,  "recoLep2PtF>0");
      sprintf(input2,  "recoLep2Pt==1 && recoLep2Eta==1 && recoLep2Phi==1 && recoLep2dEtaIn==1 && recoLep2dPhiIn==1 && recoLep2SigmaIEIE==1 && recoLep2HoverE==1 && recoLep2dxy==1 && recoLep2dz==1 && recoLep2EandP==1 && recoLep2Conv1==1 && recoLep2Conv2==1");
    } else if(NAME[j]=="AllLep1PlusIso" && (i==0 || i==2)){
      sprintf(input1,  "recoLep1PtF>0");
      sprintf(input2,  "recoLep1Pt==1 && recoLep1Eta==1 && recoLep1Phi==1 && recoLep1dEtaIn==1 && recoLep1dPhiIn==1 && recoLep1SigmaIEIE==1 && recoLep1HoverE==1 && recoLep1dxy==1 && recoLep1dz==1 && recoLep1EandP==1 && recoLep1Conv1==1 && recoLep1Conv2==1 && recoLep1CorrPFIso==1");
    } else if(NAME[j]=="AllLep1PlusIso" && i==1){
      sprintf(input1,  "recoLep1PtF>0");
      sprintf(input2,  "recoLep1Pt==1 && recoLep1Eta==1 && recoLep1Phi==1 && recoLep1Matches==1 && recoLep1dxy==1 && recoLep1dz==1 && recoLep1PixelHits==1 && recoLep1TrackerL==1&& recoLep1ptErr==1 && recoLep1DETIso==1");
    } else if(NAME[j]=="AllLep1PlusIso" && (i==3 || i==4)){
      continue;
    } else if(NAME[j]=="AllLep2PlusIso" && i==1){
      sprintf(input1,  "recoLep2PtF>0");
      sprintf(input2,  "recoLep2Pt==1 && recoLep2Eta==1 && recoLep2Phi==1 && recoLep2Matches==1 && recoLep2dxy==1 && recoLep2dz==1 && recoLep2PixelHits==1 && recoLep2TrackerL==1&& recoLep2ptErr==1 && recoLep2DETIso==1");
    } else if(NAME[j]=="AllLep2PlusIso" && (i==0 || i==3)){
      sprintf(input1,  "recoLep2PtF>0");
      sprintf(input2,  "recoLep2Pt==1 && recoLep2Eta==1 && recoLep2Phi==1 && recoLep2IsGlobal==1 && recoLep2ptErr==1 && recoLep2MuonHits==1 && recoLep2Matches==1 && recoLep2dxy==1 && recoLep2dz==1 && recoLep2PixelHits==1 && recoLep2TrackerL==1 && recoLep2CorrPFIso==1");
    } else if(NAME[j]=="AllLep2PlusIso" && (i==2 || i==4)){
      sprintf(input1,  "recoLep2PtF>0");
      sprintf(input2,  "recoLep2Pt==1 && recoLep2Eta==1 && recoLep2Phi==1 && recoLep2dEtaIn==1 && recoLep2dPhiIn==1 && recoLep2SigmaIEIE==1 && recoLep2HoverE==1 && recoLep2dxy==1 && recoLep2dz==1 && recoLep2EandP==1 && recoLep2Conv1==1 && recoLep2Conv2==1 && recoLep2CorrPFIso==1");
    } else if(NAME[j]=="AllLep1HEEPPlusIso" && (i==0 || i==2)){
      sprintf(input1,  "recoLep1PtF>0");
      sprintf(input2,  "recoLep1Pt==1 && recoLep1Eta==1 && recoLep1Phi==1 && recoLep1HEEP==1 && recoLep1HEEPIso==1");
    } else if(NAME[j]=="AllLep1HEEPPlusIso" && i==1){
      continue;
    } else if(NAME[j]=="AllLep1HEEPPlusIso" && (i==3 || i==4)){
      continue;
    } else if(NAME[j]=="AllLep2HEEPPlusIso" && (i==0 || i==1 || i==3)){
      continue;
    } else if(NAME[j]=="AllLep2HEEPPlusIso" && (i==2 || i==4)){
      sprintf(input1,  "recoLep2PtF>0");
      sprintf(input2,  "recoLep2Pt==1 && recoLep2Eta==1 && recoLep2Phi==1 && recoLep2HEEP==1 && recoLep2HEEPIso==1");
    } else {
      sprintf(input1,  "%s==0 || %s==1",name,name);
      sprintf(input2,  "%s==1",name);
    }
    
    if(NAME[j].Contains("Lep1")){
      char variable1[500]; sprintf(variable1, "recoDeltaR>>h1(600,-6,6)"); Tree1->Draw(variable1, input2, "E");
      char variable2[500]; sprintf(variable2, "recoDeltaR>>h2(600,-6,6)"); Tree2->Draw(variable2, input2, "Esame");
      char variable3[500]; sprintf(variable3, "recoDeltaR>>h3(600,-6,6)"); Tree3->Draw(variable3, input2, "Esame");
      char variable4[500]; sprintf(variable4, "recoDeltaR>>h4(600,-6,6)"); Tree4->Draw(variable4, input2, "Esame");
      char variable5[500]; sprintf(variable5, "recoDeltaR>>h5(600,-6,6)"); Tree1->Draw(variable5, input1, "Esame");
      char variable6[500]; sprintf(variable6, "recoDeltaR>>h6(600,-6,6)"); Tree2->Draw(variable6, input1, "Esame");
      char variable7[500]; sprintf(variable7, "recoDeltaR>>h7(600,-6,6)"); Tree3->Draw(variable7, input1, "Esame");
      char variable8[500]; sprintf(variable8, "recoDeltaR>>h8(600,-6,6)"); Tree4->Draw(variable8, input1, "Esame");
    } else if(NAME[j].Contains("Lep2")){
      char variable1[500]; sprintf(variable1, "recoDeltaR>>h1(600,-6,6)"); Tree1->Draw(variable1, input2, "E");
      char variable2[500]; sprintf(variable2, "recoDeltaR>>h2(600,-6,6)"); Tree2->Draw(variable2, input2, "Esame");
      char variable3[500]; sprintf(variable3, "recoDeltaR>>h3(600,-6,6)"); Tree3->Draw(variable3, input2, "Esame");
      char variable4[500]; sprintf(variable4, "recoDeltaR>>h4(600,-6,6)"); Tree4->Draw(variable4, input2, "Esame");
      char variable5[500]; sprintf(variable5, "recoDeltaR>>h5(600,-6,6)"); Tree1->Draw(variable5, input1, "Esame");
      char variable6[500]; sprintf(variable6, "recoDeltaR>>h6(600,-6,6)"); Tree2->Draw(variable6, input1, "Esame");
      char variable7[500]; sprintf(variable7, "recoDeltaR>>h7(600,-6,6)"); Tree3->Draw(variable7, input1, "Esame");
      char variable8[500]; sprintf(variable8, "recoDeltaR>>h8(600,-6,6)"); Tree4->Draw(variable8, input1, "Esame");
    }
    
    if(j>0){
      hnew1->Delete();
      hnew2->Delete();
      hnew3->Delete();
      hnew4->Delete();
      Hnew1->Delete();
      Hnew2->Delete();
      Hnew3->Delete();
      Hnew4->Delete();
      Hnew5->Delete();
      Hnew6->Delete();
      Hnew7->Delete();
      Hnew8->Delete();
    }

    Hnew1 = h1->Rebin(N-1,"Hnew1",xbins);
    Hnew2 = h2->Rebin(N-1,"Hnew2",xbins);
    Hnew3 = h3->Rebin(N-1,"Hnew3",xbins);
    Hnew4 = h4->Rebin(N-1,"Hnew4",xbins);
    Hnew5 = h5->Rebin(N-1,"Hnew5",xbins);
    Hnew6 = h6->Rebin(N-1,"Hnew6",xbins);
    Hnew7 = h7->Rebin(N-1,"Hnew7",xbins);
    Hnew8 = h8->Rebin(N-1,"Hnew8",xbins);
  
    TGraphAsymmErrors *hnew1 = new TGraphAsymmErrors(Hnew1,Hnew5);
    TGraphAsymmErrors *hnew2 = new TGraphAsymmErrors(Hnew2,Hnew6);
    TGraphAsymmErrors *hnew3 = new TGraphAsymmErrors(Hnew3,Hnew7);
    TGraphAsymmErrors *hnew4 = new TGraphAsymmErrors(Hnew4,Hnew8);
  
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
  
    hnew4->Draw("AP");
    hnew2->Draw("Psame");
    hnew3->Draw("Psame");
    hnew1->Draw("Psame");

    hnew4->SetTitle("");
    hnew4->GetYaxis()->SetTitle("Efficiency");
    hnew4->SetMaximum(1.50);
    hnew4->SetMinimum(0.00);
    hnew4->GetXaxis()->SetTitle("#DeltaR(lep1,lep2)");
    hnew4->GetXaxis()->SetRangeUser(xbins[0],xbins[N-1]);

    TLegend *pl2 = new TLegend(0.68,0.72,0.89,0.89);
    pl2->SetTextSize(0.03); 
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(hnew1, "M(X) = 1.0 TeV",  "LP");
    ple2 = pl2->AddEntry(hnew2, "M(X) = 1.5 TeV",  "LP");
    ple2 = pl2->AddEntry(hnew3, "M(X) = 2.0 TeV",  "LP");
    ple2 = pl2->AddEntry(hnew4, "M(X) = 2.5 TeV",  "LP");
    pl2->Draw();
    
    c1->SaveAs("Eff_dR_"+CHANNEL+"_"+NAME[j]+".pdf");
    cout<<"Eff_dR_"<<CHANNEL<<"_"<<NAME[j]<<".pdf"<<endl;
  }
}
