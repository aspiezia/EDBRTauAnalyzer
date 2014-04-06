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
  PLOT.push_back("jetPt");            RANGE.push_back("(150,0,1500)");      AXIS.push_back("jet pt [GeV]");
  PLOT.push_back("jetMass");          RANGE.push_back("(100,0,200)");       AXIS.push_back("jet mass [GeV]");
  PLOT.push_back("jetTau21");         RANGE.push_back("(100,0,1)");         AXIS.push_back("jet subjettiness #tau_{2}/#tau_{1}");

  for(int i=0; i<PLOT.size(); i++){
    char *plot  = PLOT[i].c_str();
    char *range = RANGE[i].c_str();
    bool save = true;
    char channel   [500];   sprintf(channel,   ""); 
    char demo      [500];   sprintf(demo,      "demo/TreeNOTaggedJet"); 
    char openTree  [500];   sprintf(openTree,  "%s",demo); 
    TString CHANNEL = channel;
    TString name = plot;
    
    TFile *file1 = TFile::Open("/afs/cern.ch/work/a/aspiezia/EXOInizio/EDBRTauAnalyzer/CMSSW_5_3_11_patch6/src/ExoDiBosonResonances/EDBRTauAnalyzer/CRAB/AnalysisNoteTrigger/data.root");
    TTree *Tree1;  Tree1 = (TTree*)file1->Get(openTree); 

    char input1[500]; sprintf(input1, "%s>>h1%s", plot, range);
  
    Tree1->Draw(input1, "(trigger320 || trigger650) && jetPt>100", "E");

    h1->SetLineWidth(2);
    h1->SetLineColor(1);
    h1->SetFillColor(1);
    h1->SetFillStyle(3001);
    h1->SetMarkerColor(1);
    h1->SetMarkerStyle(20); 
    h1->SetMarkerSize(1.3); 
    //float a1 = 1./h1->Integral();
    //h1->Scale(a1); 
  
    h1->Draw("histo");
    gPad->Update();
    TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    ps1->SetY1NDC(0.75); ps1->SetY2NDC(0.89); ps1->SetTextColor(1);
    ps1->SetX1NDC(0.70); ps1->SetX2NDC(0.89); 

    h1->SetTitle("");
    h1->GetYaxis()->SetTitle("Events");
    h1->GetXaxis()->SetTitle(AXIS[i]);
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetTitleOffset(1.2); 
    //if(PLOT[i]=="jetPt")           h1->SetMaximum(0.12);
    //if(PLOT[i]=="jetMass")         h1->SetMaximum(0.12);
    //if(PLOT[i]=="jetEta")          h1->SetMaximum(0.07);
    if(PLOT[i]=="jetTau21") h1->SetMaximum(2500000);

    //TLegend *pl2 = new TLegend(0.49,0.72,0.69,0.89);
    //pl2->SetTextSize(0.03); 
    //pl2->SetFillColor(0);
    //TLegendEntry *ple2 = pl2->AddEntry(h1, "M(X) = 1.0 TeV",  "LP");
    //ple2 = pl2->AddEntry(h2, "M(X) = 1.5 TeV",  "LP");
    //ple2 = pl2->AddEntry(h3, "M(X) = 2.0 TeV",  "LP");
    //ple2 = pl2->AddEntry(h4, "M(X) = 2.5 TeV",  "LP");
    //pl2->Draw();

    TLatex latexLabel1;
    latexLabel1.SetTextSize(0.04);
    latexLabel1.SetNDC();
    latexLabel1.DrawLatex(0.18, 0.91, "CMS");	
    TLatex latexLabel2;
    latexLabel2.SetTextSize(0.04);
    latexLabel2.SetTextFont(42);
    latexLabel2.SetNDC();
    latexLabel2.DrawLatex(0.60, 0.91, "L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");

    if(save){
      c1->SaveAs("JetId_DATA_"+name+".png");
      cout<<"Saving JetId_DATA_"+name+".png"<<endl;
    }
  }

}
