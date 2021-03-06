{
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat("rme");
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5); //title X location 
  gStyle->SetTitleY(0.96); //title Y location 
  gStyle->SetPaintTextFormat(".2f");

  using namespace std;

  vector<string> PLOT;                vector<string> RANGE;               vector<TString> AXIS;
  PLOT.push_back("MassVis");          RANGE.push_back("(30,0,300)");      AXIS.push_back("M_{vis}(lep1, lep2) [GeV]");
  PLOT.push_back("MassEff");          RANGE.push_back("(30,0,600)");      AXIS.push_back("M_{eff}(lep1, lep2) [GeV]");
  PLOT.push_back("MassCA");           RANGE.push_back("(30,0,300)");      AXIS.push_back("M_{CA}(lep1, lep2) [GeV]");
  PLOT.push_back("MassSvfit");        RANGE.push_back("(30,0,300)");      AXIS.push_back("M_{svfit}(lep1, lep2) [GeV]");
  PLOT.push_back("XMassSVFit");       RANGE.push_back("(100,0,3000)");    AXIS.push_back("M_{svfit}(Z.H) [GeV]");
  PLOT.push_back("PtSvfit");          RANGE.push_back("(50,0,2000)");     AXIS.push_back("pt(#tau,#tau) [GeV] (SvFit)");
    

  for(int i=0; i<PLOT.size(); i++){
    char *plot  = PLOT[i].c_str();
    char *range = RANGE[i].c_str();
    bool save = true;
    char channel   [500];   sprintf(channel,   "EleMuo"); 
    char demo      [500];   sprintf(demo,      "demo/Tree"); 
    TString name = plot;
    for(int j=0; j<5; j++){
      if(j==0) sprintf(channel,   "EleMuo"); 
      if(j==1) sprintf(channel,   "MuoMuo"); 
      if(j==2) sprintf(channel,   "EleEle"); 
      if(j==3) sprintf(channel,   "MuoTau"); 
      if(j==4) sprintf(channel,   "EleTau"); 
      char openTree  [500];   sprintf(openTree,  "%s%s",demo,channel); 
      TString CHANNEL = channel;
    
      TFile *file1 = TFile::Open("../../RISULTATI/analyzer_290514/ZH1000.root");
      TFile *file2 = TFile::Open("../../RISULTATI/analyzer_290514/ZH1500.root");
      TFile *file3 = TFile::Open("../../RISULTATI/analyzer_290514/ZH2000.root");
      TFile *file4 = TFile::Open("../../RISULTATI/analyzer_290514/ZH2500.root");
      TTree *Tree1;  Tree1 = (TTree*)file1->Get(openTree); 
      TTree *Tree2;  Tree2 = (TTree*)file2->Get(openTree); 
      TTree *Tree3;  Tree3 = (TTree*)file3->Get(openTree); 
      TTree *Tree4;  Tree4 = (TTree*)file4->Get(openTree); 

      char input1[500]; sprintf(input1, "%s>>h1%s", plot, range);
      char input2[500]; sprintf(input2, "%s>>h2%s", plot, range);
      char input3[500]; sprintf(input3, "%s>>h3%s", plot, range);
      char input4[500]; sprintf(input4, "%s>>h4%s", plot, range);
      Tree1->Draw(input1, "", "E");
      Tree2->Draw(input2, "", "Esame");
      Tree3->Draw(input3, "", "Esame");
      Tree4->Draw(input4, "", "Esame");

      h1->SetLineWidth(2);
      h2->SetLineWidth(2);
      h3->SetLineWidth(2);
      h4->SetLineWidth(2);
  
      h1->SetLineColor(1);
      h2->SetLineColor(2);
      h3->SetLineColor(kGreen+2);
      h4->SetLineColor(4);
  
      h1->SetMarkerColor(1);
      h2->SetMarkerColor(2);
      h3->SetMarkerColor(kGreen+2);
      h4->SetMarkerColor(4);
  
      h1->SetMarkerStyle(20); 
      h2->SetMarkerStyle(20); 
      h3->SetMarkerStyle(20); 
      h4->SetMarkerStyle(20); 
  
      h1->SetMarkerSize(1.3); 
      h2->SetMarkerSize(1.3); 
      h3->SetMarkerSize(1.3); 
      h4->SetMarkerSize(1.3); 

      float a1 = 1./h1->Integral();
      float a2 = 1./h2->Integral();
      float a3 = 1./h3->Integral();
      float a4 = 1./h4->Integral();
  
      h1->Scale(a1); 
      h2->Scale(a2); 
      h3->Scale(a3); 
      h4->Scale(a4); 
  
      h1->Draw("histo");
      gPad->Update();
      TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
      ps1->SetY1NDC(0.75); ps1->SetY2NDC(0.89); ps1->SetTextColor(1);
      ps1->SetX1NDC(0.70); ps1->SetX2NDC(0.89); 
      h2->Draw("histo sames");
      gPad->Update();
      TPaveStats *ps2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
      ps2->SetY1NDC(0.60); ps2->SetY2NDC(0.74); ps2->SetTextColor(2);
      ps2->SetX1NDC(0.70); ps2->SetX2NDC(0.89); 
      h3->Draw("histo sames");
      gPad->Update();
      TPaveStats *ps3 = (TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
      ps3->SetY1NDC(0.45); ps3->SetY2NDC(0.59); ps3->SetTextColor(kGreen+2);
      ps3->SetX1NDC(0.70); ps3->SetX2NDC(0.89); 
      h4->Draw("histo sames");
      gPad->Update();
      TPaveStats *ps4 = (TPaveStats*)h4->GetListOfFunctions()->FindObject("stats");
      ps4->SetY1NDC(0.30); ps4->SetY2NDC(0.44); ps4->SetTextColor(4);
      ps4->SetX1NDC(0.70); ps4->SetX2NDC(0.89); 

      if(CHANNEL=="EleMuo") h1->SetTitle("Electron + Muon Category");
      if(CHANNEL=="MuoMuo") h1->SetTitle("Muon + Muon Category");
      if(CHANNEL=="EleEle") h1->SetTitle("Electron + Electron Category");
      if(CHANNEL=="MuoTau") h1->SetTitle("Muon + Tau_{h} Category");
      if(CHANNEL=="EleTau") h1->SetTitle("Electron + Tau_{h} Category");
      h1->GetYaxis()->SetTitle("Events (normalized unit)");
      h1->GetXaxis()->SetTitle(AXIS[i]);
      h1->GetYaxis()->SetTitleSize(0.045);
      h1->GetXaxis()->SetTitleSize(0.045);
      h1->GetYaxis()->SetLabelSize(0.045);
      h1->GetXaxis()->SetLabelSize(0.045);
      h1->GetYaxis()->SetTitleOffset(1.0); 

      TLegend *pl2 = new TLegend(0.49,0.72,0.69,0.89);
      pl2->SetTextSize(0.03); 
      pl2->SetFillColor(0);
      TLegendEntry *ple2 = pl2->AddEntry(h1, "M(X) = 1.0 TeV",  "LP");
      ple2 = pl2->AddEntry(h2, "M(X) = 1.5 TeV",  "LP");
      ple2 = pl2->AddEntry(h3, "M(X) = 2.0 TeV",  "LP");
      ple2 = pl2->AddEntry(h4, "M(X) = 2.5 TeV",  "LP");
      pl2->Draw();

  
      if(save){
	c1->SaveAs("ZH_"+name+"_"+CHANNEL+".pdf");
	cout<<"Saving ZH_"<<name<<"_"<<CHANNEL<<".pdf"<<endl;
      }
    }
  }
}
