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

  bool QCD = false;
  bool save = true;
  char channel   [500];   sprintf(channel,   "EleTau"); 
  char demo      [500];   sprintf(demo,      "demo/Tree"); 
  char demoUsual [500];   sprintf(demoUsual, "demo/TreeUsual"); 
  char openTree  [500];   sprintf(openTree,  "%s%s",demo,channel); 
  char openTreeU [500];   sprintf(openTreeU, "%s%s",demoUsual,channel); 
  char openTree1 [500];   sprintf(openTree1, "demo/TreeQCDMuoTau");
  char openTree2 [500];   sprintf(openTree2, "demo/TreeQCDEleTau");
  char openTree3 [500];   sprintf(openTree3, "demo/TreeQCDUsual");
  TString CHANNEL = channel;

  int Nbins = 10;
  Double_t xbins[Nbins] = {20, 50, 80, 110, 150, 200, 250, 300, 400, 500};
  
  TFile *file1 = TFile::Open("../../RISULTATI/tauAnalyzer_120514/ZH1000.root");
  TFile *file2 = TFile::Open("../../RISULTATI/tauAnalyzer_120514/ZH2500.root");
  TFile *file3 = TFile::Open("../../RISULTATI/tauAnalyzer_120514/ZH1000.root");
  TFile *file4 = TFile::Open("../../RISULTATI/tauAnalyzer_120514/ZH2500.root");
  TFile *file5 = TFile::Open("../../RISULTATI/tauAnalyzer_120514/ZH1000.root"); 
  
  TTree *Tree1;  
  TTree *Tree2;
  TTree *Tree3;
  TTree *Tree4;
  
  if(!QCD){
    Tree1 = (TTree*)file1->Get(openTree);  
    Tree2 = (TTree*)file2->Get(openTree);  
    Tree3 = (TTree*)file3->Get(openTreeU); 
    Tree4 = (TTree*)file4->Get(openTreeU); 
  } else {
    Tree1 = (TTree*)file5->Get(openTree1); 
    Tree2 = (TTree*)file5->Get(openTree2); 
    Tree3 = (TTree*)file5->Get(openTree3); 
    Tree4 = (TTree*)file5->Get(openTree1); 
  }
  
  TString discriminant = "decay"; int N = 2; 
  for(int j=0; j<6; j++){
    if(j==0) {discriminant = "decay";    N = 2;    }
    if(j==1) {discriminant = "iso";      N = 40*2; }
    if(j==2) {discriminant = "electron"; N = 11*2; }
    if(j==3) {discriminant = "muon";     N = 12*2; }
    if(j==4) {discriminant = "final";    N = 1;    }
    if(j==5) {discriminant = "finalOLD"; N = 1;    }

    vector<TString> NAME;
    if(discriminant=="decay") {
      NAME.push_back("decayModeFinding");
      NAME.push_back("decayModeFindingNewDMs");
    } else if(discriminant=="iso") {
      NAME.push_back("byVLooseCombinedIsolationDeltaBetaCorr");
      NAME.push_back("byLooseIsolation");
      NAME.push_back("byLooseCombinedIsolationDeltaBetaCorr");
      NAME.push_back("byMediumCombinedIsolationDeltaBetaCorr");
      NAME.push_back("byTightCombinedIsolationDeltaBetaCorr");
      NAME.push_back("byCombinedIsolationDeltaBetaCorrRaw");
      NAME.push_back("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      NAME.push_back("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      NAME.push_back("byTightCombinedIsolationDeltaBetaCorr3Hits");
      NAME.push_back("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      NAME.push_back("chargedIsoPtSum");
      NAME.push_back("neutralIsoPtSum");
      NAME.push_back("byIsolationMVA3oldDMwoLTraw");
      NAME.push_back("byVLooseIsolationMVA3oldDMwoLT");
      NAME.push_back("byLooseIsolationMVA3oldDMwoLT");
      NAME.push_back("byMediumIsolationMVA3oldDMwoLT");
      NAME.push_back("byTightIsolationMVA3oldDMwoLT");
      NAME.push_back("byVTightIsolationMVA3oldDMwoLT");
      NAME.push_back("byVVTightIsolationMVA3oldDMwoLT");
      NAME.push_back("byIsolationMVA3oldDMwLTraw");
      NAME.push_back("byVLooseIsolationMVA3oldDMwLT");
      NAME.push_back("byLooseIsolationMVA3oldDMwLT");
      NAME.push_back("byMediumIsolationMVA3oldDMwLT");
      NAME.push_back("byTightIsolationMVA3oldDMwLT");
      NAME.push_back("byVTightIsolationMVA3oldDMwLT");
      NAME.push_back("byVVTightIsolationMVA3oldDMwLT");
      NAME.push_back("byIsolationMVA3newDMwoLTraw");
      NAME.push_back("byVLooseIsolationMVA3newDMwoLT");
      NAME.push_back("byLooseIsolationMVA3newDMwoLT");
      NAME.push_back("byMediumIsolationMVA3newDMwoLT");
      NAME.push_back("byTightIsolationMVA3newDMwoLT");
      NAME.push_back("byVTightIsolationMVA3newDMwoLT");
      NAME.push_back("byVVTightIsolationMVA3newDMwoLT");
      NAME.push_back("byIsolationMVA3newDMwLTraw");
      NAME.push_back("byVLooseIsolationMVA3newDMwLT");
      NAME.push_back("byLooseIsolationMVA3newDMwLT");
      NAME.push_back("byMediumIsolationMVA3newDMwLT");
      NAME.push_back("byTightIsolationMVA3newDMwLT");
      NAME.push_back("byVTightIsolationMVA3newDMwLT");
      NAME.push_back("byVVTightIsolationMVA3newDMwLT");
    } else if(discriminant=="electron") {
      NAME.push_back("againstElectronLoose");
      NAME.push_back("againstElectronMedium");
      NAME.push_back("againstElectronTight");
      NAME.push_back("againstElectronMVA5raw");
      NAME.push_back("againstElectronMVA5category");
      NAME.push_back("againstElectronVLooseMVA5");
      NAME.push_back("againstElectronLooseMVA5");
      NAME.push_back("againstElectronMediumMVA5");
      NAME.push_back("againstElectronTightMVA5");
      NAME.push_back("againstElectronVTightMVA5");
      NAME.push_back("againstElectronDeadECAL");
    } else if(discriminant=="muon") {
      NAME.push_back("againstMuonLoose");
      NAME.push_back("againstMuonMedium");
      NAME.push_back("againstMuonTight");
      NAME.push_back("againstMuonLoose2");
      NAME.push_back("againstMuonMedium2");
      NAME.push_back("againstMuonTight2");
      NAME.push_back("againstMuonLoose3");
      NAME.push_back("againstMuonTight3");
      NAME.push_back("againstMuonMVAraw");
      NAME.push_back("againstMuonLooseMVA");
      NAME.push_back("againstMuonMediumMVA");
      NAME.push_back("againstMuonTightMVA");
    }

    for(int i=0; i<N; i++){
      if(discriminant=="iso" && (i+1==6 || i+1==10 || i+1==11 || i+1==12 || i+1==13 || i+1==20 || i+1==27 || i+1==34)) continue;
      if(discriminant=="muon" && i+1==9) continue;
      if(discriminant=="electron" && (i+1==4 || i+1==5)) continue;
      if(discriminant=="iso" && (i+1==6+N/2 || i+1==10+N/2 || i+1==11+N/2 || i+1==12+N/2 || i+1==13+N/2 || i+1==20+N/2 || i+1==27+N/2 || i+1==34+N/2)) continue;
      if(discriminant=="muon" && i+1==9+N/2) continue;
      if(discriminant=="electron" && (i+1==4+N/2 || i+1==5+N/2)) continue;
      char input1[500]; 
    
      if(i>=0 && i<9){
	if(discriminant=="decay")          sprintf(input1,  "tauDecay0%d>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1);
	else if (discriminant=="iso")      sprintf(input1,  "tauIso0%d>0.5 && tauDecay01>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1);
	else if (discriminant=="electron") sprintf(input1,  "tauAgainstElectron0%d>0.5 && tauDecay01>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1);
	else if (discriminant=="muon")     sprintf(input1,  "tauAgainstMuon0%d>0.5 && tauDecay01>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1);
      } else if (i>=9 && i<N/2){
	if(discriminant=="decay")          sprintf(input1,  "tauDecay%d>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1);
	else if (discriminant=="iso")      sprintf(input1,  "tauIso%d>0.5 && tauDecay01>0.5  && tauPt>20 && abs(tauEta)<2.4", i+1);
	else if (discriminant=="electron") sprintf(input1,  "tauAgainstElectron%d>0.5 && tauDecay01>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1);
	else if (discriminant=="muon")     sprintf(input1,  "tauAgainstMuon%d>0.5 && tauDecay01>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1);
      } else if(i>=N/2 && i<N/2+9){
	if(discriminant=="decay")          sprintf(input1,  "tauDecay0%d>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1-N/2);
	else if (discriminant=="iso")      sprintf(input1,  "tauIso0%d>0.5 && tauDecay02>0.5  && tauPt>20 && abs(tauEta)<2.4", i+1-N/2);
	else if (discriminant=="electron") sprintf(input1,  "tauAgainstElectron0%d>0.5 && tauDecay02>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1-N/2);
	else if (discriminant=="muon")     sprintf(input1,  "tauAgainstMuon0%d>0.5 && tauDecay02>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1-N/2);
      } else if(i>=N/2+9) {
	if(discriminant=="decay")          sprintf(input1,  "tauDecay%d>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1-N/2);
	else if (discriminant=="iso")      sprintf(input1,  "tauIso%d>0.5 && tauDecay02>0.5  && tauPt>20 && abs(tauEta)<2.4", i+1-N/2);
	else if (discriminant=="electron") sprintf(input1,  "tauAgainstElectron%d>0.5 && tauDecay02>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1-N/2);
	else if (discriminant=="muon")     sprintf(input1,  "tauAgainstMuon%d>0.5 && tauDecay02>0.5 && tauPt>20 && abs(tauEta)<2.4", i+1-N/2);
      }
      if(discriminant=="final") sprintf(input1, "tauDecay02>0.5 && tauIso28>0.5 && tauAgainstElectron01>0.5 && tauAgainstMuon01>0.5 && tauPt>20 && abs(tauEta)<2.4");
      if(discriminant=="finalOLD") sprintf(input1, "tauDecay01>0.5 && tauIso01>0.5 && tauAgainstElectron01>0.5 && tauAgainstMuon01>0.5 && tauPt>20 && abs(tauEta)<2.4");

      Tree1->Draw("tauPt>>h1(1000,0,500)", input1, "E");
      Tree2->Draw("tauPt>>h2(1000,0,500)", input1, "Esame");
      Tree3->Draw("tauPt>>h3(1000,0,500)", input1, "Esame");
      Tree4->Draw("tauPt>>h4(1000,0,500)", input1, "Esame");
      Tree1->Draw("tauPt>>h5(1000,0,500)", "", "Esame");
      Tree2->Draw("tauPt>>h6(1000,0,500)", "", "Esame");
      Tree3->Draw("tauPt>>h7(1000,0,500)", "", "Esame");
      Tree4->Draw("tauPt>>h8(1000,0,500)", "", "Esame");
 
      if(i>0){
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
    
      Hnew1 = h1->Rebin(Nbins-1,"Hnew1",xbins);
      Hnew2 = h2->Rebin(Nbins-1,"Hnew2",xbins);
      Hnew3 = h3->Rebin(Nbins-1,"Hnew3",xbins);
      Hnew4 = h4->Rebin(Nbins-1,"Hnew4",xbins);
      Hnew5 = h5->Rebin(Nbins-1,"Hnew5",xbins);
      Hnew6 = h6->Rebin(Nbins-1,"Hnew6",xbins);
      Hnew7 = h7->Rebin(Nbins-1,"Hnew7",xbins);
      Hnew8 = h8->Rebin(Nbins-1,"Hnew8",xbins);
    
      TGraphAsymmErrors *hnew1 = new TGraphAsymmErrors(Hnew1,Hnew5);
      TGraphAsymmErrors *hnew2 = new TGraphAsymmErrors(Hnew2,Hnew6);
      TGraphAsymmErrors *hnew3 = new TGraphAsymmErrors(Hnew3,Hnew7);
      TGraphAsymmErrors *hnew4 = new TGraphAsymmErrors(Hnew4,Hnew8);

      //hnew1 = h1->Rebin(8,"hnew1",xbins);
      //hnew2 = h2->Rebin(8,"hnew2",xbins);
      //hnew3 = h3->Rebin(8,"hnew3",xbins);
      //hnew4 = h4->Rebin(8,"hnew4",xbins);
      //hnew5 = h5->Rebin(8,"hnew5",xbins);
      //hnew6 = h6->Rebin(8,"hnew6",xbins);
      //hnew7 = h7->Rebin(8,"hnew7",xbins);
      //hnew8 = h8->Rebin(8,"hnew8",xbins);
      //
      //hnew1->Divide(hnew5);
      //hnew2->Divide(hnew6);
      //hnew3->Divide(hnew7);
      //hnew4->Divide(hnew8);

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

      hnew1->Draw("AP");
      hnew2->Draw("Psame");
      hnew3->Draw("Psame");
      if(!QCD) hnew4->Draw("Psame");

      if(discriminant=="final" || discriminant=="finalOLD") hnew1->SetTitle("Tau Selection");
      else if(discriminant=="decay" || (discriminant=="iso" && i<40) || (discriminant=="electron" && i<11) || (discriminant=="muon" && i<12)) hnew1->SetTitle(NAME[i]);
      else hnew1->SetTitle(NAME[i-N/2]);
      if(!QCD){
	hnew1->GetXaxis()->SetTitle("tau pt [GeV]");
	hnew1->GetYaxis()->SetTitle("Efficiency");
	hnew1->SetMaximum(1.50);
	hnew1->SetMinimum(0.00);
      } else {
	hnew1->GetXaxis()->SetTitle("jet pt [GeV]");
	hnew1->GetYaxis()->SetTitle("Fake Rate");
	if(discriminant=="iso" || discriminant=="final" || discriminant=="finalOLD"){
	  hnew1->SetMaximum(1.0);
	  hnew1->SetMinimum(0.0001);
	  c1->SetLogy();
	} else {
	  hnew1->SetMaximum(1.50);
	  hnew1->SetMinimum(0.00);
	}
      }

      TLegend *pl2 = new TLegend(0.43,0.72,0.89,0.89);
      pl2->SetTextSize(0.03); 
      pl2->SetFillColor(0);
      if(!QCD){
	TLegendEntry *ple2 = pl2->AddEntry(hnew1, "Cleaned Collection - M(X) = 1.0 TeV",  "LP");
	ple2 = pl2->AddEntry(hnew2, "Cleaned Collection - M(X) = 2.5 TeV",  "LP");
	ple2 = pl2->AddEntry(hnew3, "Usual Collection - M(X) = 1.0 TeV",  "LP");
	ple2 = pl2->AddEntry(hnew4, "Usual Collection - M(X) = 2.5 TeV",  "LP");
	pl2->Draw();
      } else {
	TLegendEntry *ple2 = pl2->AddEntry(hnew1, "Cleaned muon-tau collection",  "LP");
	ple2 = pl2->AddEntry(hnew2, "Cleaned electron-tau collection",  "LP");
	ple2 = pl2->AddEntry(hnew3, "Usual collection",  "LP");
	pl2->Draw();
      }

      if(save){
	if (!QCD) {
	  if(discriminant=="final" || discriminant=="finalOLD"){
	    c1->SaveAs(CHANNEL+"_"+discriminant+".pdf");
	    cout<<i+1<<") Saving "<<CHANNEL<<"_"<<discriminant<<".pdf"<<endl;
	  } else if(discriminant=="decay" || (discriminant=="iso" && i<40) || (discriminant=="electron" && i<11) || (discriminant=="muon" && i<12)) {
	    c1->SaveAs(CHANNEL+"_"+discriminant+"_"+NAME[i]+".pdf");
	    cout<<i+1<<") Saving "<<CHANNEL<<"_"<<discriminant<<"_"<<NAME[i]<<".pdf"<<endl;
	  } else {
	    c1->SaveAs(CHANNEL+"_"+discriminant+"2_"+NAME[i-N/2]+".pdf");
	    cout<<i+1<<") Saving "<<CHANNEL<<"_"<<discriminant<<"2_"<<NAME[i-N/2]<<".pdf"<<endl;
	  }
	} else {
	  if(discriminant=="final" || discriminant=="finalOLD"){
	    c1->SaveAs("QCD_"+discriminant+".pdf");
	    cout<<i+1<<") Saving "<<"QCD_"<<discriminant<<".pdf"<<endl;
	  } else if(discriminant=="decay" || (discriminant=="iso" && i<40) || (discriminant=="electron" && i<11) || (discriminant=="muon" && i<12)) {
	    c1->SaveAs("QCD_"+discriminant+"_"+NAME[i]+".pdf");
	    cout<<i+1<<") Saving "<<"QCD_"<<discriminant<<"_"<<NAME[i]<<".pdf"<<endl;
	  } else {
	    c1->SaveAs("QCD_"+discriminant+"2_"+NAME[i-N/2]+".pdf");
	    cout<<i+1<<") Saving "<<"QCD_"<<discriminant<<"2_"<<NAME[i-N/2]<<".pdf"<<endl;
	  }
	}
      }
    }
  }
}
