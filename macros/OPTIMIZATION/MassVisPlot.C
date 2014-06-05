{
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5); //title X location 
  gStyle->SetTitleY(0.96); //title Y location 
  gStyle->SetPaintTextFormat(".2f");

  TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
  using namespace std;
  bool save=true; 
  bool log=false;
  bool SL=false;
  float lep1PtCut=10; float tauPtCut=35; int bTagCut=1; float deltaRCutFL=1; float deltaRCutSL=1; float MassSvfitCut=10000; float METCutFL = 100; float METCutSL = 50; 
  float MassVisCutSL=0;  
  float MassVisCutFL1=0; //ELE-MUO
  float MassVisCutFL2=0; //ELE-ELE e MUO-MUO
  bool beginning=true;

  vector<string> PLOT;           vector<int> BIN;   vector<float> MIN;  vector<float> MAX;   vector<float> MAXY;   vector<TString> AXIS;
  PLOT.push_back("MassVis");     BIN.push_back(30); MIN.push_back(0);   MAX.push_back(300);   MAXY.push_back(6);   AXIS.push_back("M(lep1,lep2) [GeV]");
  
  TH1F *massoptimization1  = new TH1F("massoptimization1","massoptimization1",20,-2.5,97.5);
  TH1F *massoptimization2  = new TH1F("massoptimization2","massoptimization2",20,-2.5,97.5);
  TH1F *SignalEfficiency1 = new TH1F("SignalEfficiency1","SignalEfficiency1",20,-2.5,97.5);
  TH1F *SignalEfficiency2 = new TH1F("SignalEfficiency2","SignalEfficiency2",20,-2.5,97.5);
  TH1F *BackgroundYield1  = new TH1F("BackgroundYield1","BackgroundYield1",20,-2.5,97.5);
  TH1F *BackgroundYield2  = new TH1F("BackgroundYield2","BackgroundYield2",20,-2.5,97.5);

  for(int type=0; type<2; type++){
    bool SL=true;
    if(type==0) SL=true;
    if(type==1) SL=false;
    float ZH1000Integral = 0;
    float ZH2500Integral = 0;
    float ADD = (1+ massoptimization1->GetXaxis()->GetBinCenter(massoptimization1->GetXaxis()->GetNbins())-
		 massoptimization1->GetXaxis()->GetBinCenter(1))/massoptimization1->GetXaxis()->GetNbins();
    float mass = 0. - ADD;
    for(int i=0; i<massoptimization1->GetXaxis()->GetNbins(); i++){
      if(beginning){ if(i!=0) continue;}
      mass = mass + ADD;
      char *plot = PLOT[0].c_str();
      TString name = PLOT[0];
      int bin=BIN[0]; 
      float min=MIN[0]; 
      float max=MAX[0]; 
      float maxy=MAXY[0];
      TString axis = AXIS[0];
      
      if(!SL) {maxy=3;}
      if(SL ) {maxy=6;}
      
      TH1F *ZH1000 = new TH1F("","",bin,min,max);
      TH1F *ZH1500 = new TH1F("","",bin,min,max);
      TH1F *ZH2000 = new TH1F("","",bin,min,max);
      TH1F *ZH2500 = new TH1F("","",bin,min,max);
      TH1F *DY100 = new TH1F("","",bin,min,max);
      TH1F *DY70 = new TH1F("","",bin,min,max);
      TH1F *DYM50_100 = new TH1F("","",bin,min,max);
      TH1F *DYM50_70 = new TH1F("","",bin,min,max);											     
      TH1F *QCD1000 = new TH1F("","",bin,min,max);
      TH1F *QCD250 = new TH1F("","",bin,min,max);
      TH1F *QCD500 = new TH1F("","",bin,min,max);
      TH1F *TT = new TH1F("","",bin,min,max);
      TH1F *WJets180 = new TH1F("","",bin,min,max);
      TH1F *WW = new TH1F("","",bin,min,max);
      TH1F *WZ = new TH1F("","",bin,min,max);
      TH1F *ZZ = new TH1F("","",bin,min,max);
      
      int jMIN=100; int jMAX=10;
      if(SL){
	jMIN=3;
	jMAX=5;
      } else {
	jMIN=0;
	jMAX=3;
      }
      for(int j=jMIN; j<jMAX; j++){
	char BTAG[1000]; sprintf(BTAG, "");
	if(bTagCut==1) {sprintf(BTAG, " && nbtagsL1<1");}
	if(bTagCut==2) {sprintf(BTAG, " && nbtagsL1<2");}
	if(bTagCut==3) {sprintf(BTAG, " && nbtagsM1<1");}
	if(bTagCut==4) {sprintf(BTAG, " && nbtagsM1<2");}
	if(bTagCut==5) {sprintf(BTAG, " && nbtagsT1<1");}
	if(bTagCut==6) {sprintf(BTAG, " && nbtagsT1<2");}
	char CUTPre[1000]; sprintf(CUTPre, "PUWeight*(trigger==1 && MassSvfit<%f %s",MassSvfitCut,BTAG);
	char CUT [1000]; 
	char openTree[500];
	if(j==0){ 
	  sprintf(openTree, "demo/TreeEleMuo"); 
	  sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && MassVis>%f)",CUTPre,METCutFL,MassVisCutFL1,lep1PtCut,deltaRCutFL,mass);
	}else if(j==1){ 
	  sprintf(openTree, "demo/TreeMuoMuo"); 
	  sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && MassVis>%f)",CUTPre,METCutFL,MassVisCutFL2,lep1PtCut,deltaRCutFL,mass);
	}else if(j==2){ 
	  sprintf(openTree, "demo/TreeEleEle"); 
	  sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && MassVis>%f)",CUTPre,METCutFL,MassVisCutFL2,lep1PtCut,deltaRCutFL,mass);
	}else if(j==3){ 
	  sprintf(openTree, "demo/TreeMuoTau");
	  sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && MassVis>%f)",CUTPre,METCutSL,MassVisCutSL, tauPtCut, deltaRCutSL,mass);
	}else if(j==4){ 
	  sprintf(openTree, "demo/TreeEleTau"); 
	  sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && MassVis>%f)",CUTPre,METCutSL,MassVisCutSL, tauPtCut, deltaRCutSL,mass);
	}
	
	TFile *fileAA = TFile::Open("../../RISULTATI/analyzer_290514/ZH1000.root");    TTree *TreeAA = (TTree*)fileAA->Get(openTree); 
	TFile *fileBB = TFile::Open("../../RISULTATI/analyzer_290514/ZH1500.root");    TTree *TreeBB = (TTree*)fileBB->Get(openTree); 
	TFile *fileCC = TFile::Open("../../RISULTATI/analyzer_290514/ZH2000.root");    TTree *TreeCC = (TTree*)fileCC->Get(openTree);
	TFile *fileDD = TFile::Open("../../RISULTATI/analyzer_290514/ZH2500.root");    TTree *TreeDD = (TTree*)fileDD->Get(openTree); 
	TFile *file02 = TFile::Open("../../RISULTATI/analyzer_290514/DY100.root");     TTree *Tree02 = (TTree*)file02->Get(openTree); 
	TFile *file03 = TFile::Open("../../RISULTATI/analyzer_290514/DY70.root");      TTree *Tree03 = (TTree*)file03->Get(openTree); 
	TFile *file04 = TFile::Open("../../RISULTATI/analyzer_290514/DYM50_100.root"); TTree *Tree04 = (TTree*)file04->Get(openTree); 
	TFile *file05 = TFile::Open("../../RISULTATI/analyzer_290514/DYM50_70.root");  TTree *Tree05 = (TTree*)file05->Get(openTree); 
	TFile *file06 = TFile::Open("../../RISULTATI/analyzer_290514/QCD1000.root");   TTree *Tree06 = (TTree*)file06->Get(openTree); 
	TFile *file07 = TFile::Open("../../RISULTATI/analyzer_290514/QCD250.root");    TTree *Tree07 = (TTree*)file07->Get(openTree); 
	TFile *file08 = TFile::Open("../../RISULTATI/analyzer_290514/QCD500.root");    TTree *Tree08 = (TTree*)file08->Get(openTree); 
	TFile *file09 = TFile::Open("../../RISULTATI/analyzer_290514/TT.root");        TTree *Tree09 = (TTree*)file09->Get(openTree); 
	TFile *file10 = TFile::Open("../../RISULTATI/analyzer_290514/WJets180.root");  TTree *Tree10 = (TTree*)file10->Get(openTree); 
	TFile *file11 = TFile::Open("../../RISULTATI/analyzer_290514/WW.root");        TTree *Tree11 = (TTree*)file11->Get(openTree); 
	TFile *file12 = TFile::Open("../../RISULTATI/analyzer_290514/WZ.root");        TTree *Tree12 = (TTree*)file12->Get(openTree); 
	TFile *file13 = TFile::Open("../../RISULTATI/analyzer_290514/ZZ.root");        TTree *Tree13 = (TTree*)file13->Get(openTree); 
	
	char inputAA[50]; sprintf(inputAA, "%s>>hAA(%i,%f,%f)", plot,bin,min,max);
	char inputBB[50]; sprintf(inputBB, "%s>>hBB(%i,%f,%f)", plot,bin,min,max);
	char inputCC[50]; sprintf(inputCC, "%s>>hCC(%i,%f,%f)", plot,bin,min,max);
	char inputDD[50]; sprintf(inputDD, "%s>>hDD(%i,%f,%f)", plot,bin,min,max);
	char input02[50]; sprintf(input02, "%s>>h02(%i,%f,%f)", plot,bin,min,max);
	char input03[50]; sprintf(input03, "%s>>h03(%i,%f,%f)", plot,bin,min,max);
	char input04[50]; sprintf(input04, "%s>>h04(%i,%f,%f)", plot,bin,min,max);
	char input05[50]; sprintf(input05, "%s>>h05(%i,%f,%f)", plot,bin,min,max);
	char input06[50]; sprintf(input06, "%s>>h06(%i,%f,%f)", plot,bin,min,max);
	char input07[50]; sprintf(input07, "%s>>h07(%i,%f,%f)", plot,bin,min,max);
	char input08[50]; sprintf(input08, "%s>>h08(%i,%f,%f)", plot,bin,min,max);
	char input09[50]; sprintf(input09, "%s>>h09(%i,%f,%f)", plot,bin,min,max);
	char input10[50]; sprintf(input10, "%s>>h10(%i,%f,%f)", plot,bin,min,max);
	char input11[50]; sprintf(input11, "%s>>h11(%i,%f,%f)", plot,bin,min,max);
	char input12[50]; sprintf(input12, "%s>>h12(%i,%f,%f)", plot,bin,min,max);
	char input13[50]; sprintf(input13, "%s>>h13(%i,%f,%f)", plot,bin,min,max);
	
	TreeAA->Draw(inputAA,CUT); if(TreeAA->Draw(inputAA,CUT)) {ZH1000->Add(hAA);}
	TreeBB->Draw(inputBB,CUT); if(TreeBB->Draw(inputBB,CUT)) {ZH1500->Add(hBB);}
	TreeCC->Draw(inputCC,CUT); if(TreeCC->Draw(inputCC,CUT)) {ZH2000->Add(hCC);}
	TreeDD->Draw(inputDD,CUT); if(TreeDD->Draw(inputDD,CUT)) {ZH2500->Add(hDD);}
	Tree02->Draw(input02,CUT); if(Tree02->Draw(input02,CUT)) {DY100->Add(h02);}
	Tree03->Draw(input03,CUT); if(Tree03->Draw(input03,CUT)) {DY70->Add(h03);}
	Tree04->Draw(input04,CUT); if(Tree04->Draw(input04,CUT)) {DYM50_100->Add(h04);}
	Tree05->Draw(input05,CUT); if(Tree05->Draw(input05,CUT)) {DYM50_70->Add(h05);}
	Tree06->Draw(input06,CUT); if(Tree06->Draw(input06,CUT)) {QCD1000->Add(h06);}
	Tree07->Draw(input07,CUT); if(Tree07->Draw(input07,CUT)) {QCD250->Add(h07);}
	Tree08->Draw(input08,CUT); if(Tree08->Draw(input08,CUT)) {QCD500->Add(h08);}
	Tree09->Draw(input09,CUT); if(Tree09->Draw(input09,CUT)) {TT->Add(h09);}
	Tree10->Draw(input10,CUT); if(Tree10->Draw(input10,CUT)) {WJets180->Add(h10);}
	Tree11->Draw(input11,CUT); if(Tree11->Draw(input11,CUT)) {WW->Add(h11);}
	Tree12->Draw(input12,CUT); if(Tree12->Draw(input12,CUT)) {WZ->Add(h12);}
	Tree13->Draw(input13,CUT); if(Tree13->Draw(input13,CUT)) {ZZ->Add(h13);}
	
	fileAA->Close(); fileBB->Close(); fileCC->Close(); fileDD->Close();
	file13->Close(); file02->Close(); file03->Close(); file04->Close();
	file05->Close(); file06->Close(); file07->Close(); file08->Close();
	file09->Close(); file10->Close(); file11->Close(); file12->Close();
	fileAA->Delete(); fileBB->Delete(); fileCC->Delete(); fileDD->Delete();
	file13->Delete(); file02->Delete(); file03->Delete(); file04->Delete();
	file05->Delete(); file06->Delete(); file07->Delete(); file08->Delete();
	file09->Delete(); file10->Delete(); file11->Delete(); file12->Delete();
      }
      
      DY100->Sumw2();
      DY70->Sumw2();
      DYM50_100->Sumw2();
      DYM50_70->Sumw2();
      QCD1000->Sumw2();
      QCD500->Sumw2();
      TT->Sumw2();
      WW->Sumw2();
      WZ->Sumw2();
      ZZ->Sumw2();
      WJets180->Sumw2();
      
      float w_DY100     = ( 39.100*19702./12511326.);
      float w_DY70      = ( 62.900*19702./11764538.);
      float w_DYM50_100 = (  4.220*19702./4146124.0);
      float w_DYM50_70  = ( 11.050*19702./5389313.0);
      float w_QCD1000   = (204.000*19702./13843863.);
      float w_QCD500    = (8426.00*19702./30599292.);
      float w_QCD250    = (276000.*19702./27062078.);
      float w_TT        = (225.197*19702./21675970.);
      float w_WW        = (57.1097*19702./10000431.);
      float w_WZ        = ( 33.210*19702./10000283.);
      float w_ZZ        = (  8.059*19702./9799908.0);
      float w_WJets180  = ( 23.500*19702./9739464.0);
      
      DY100->Scale(w_DY100);
      DY70->Scale(w_DY70);
      DYM50_100->Scale(w_DYM50_100);
      DYM50_70->Scale(w_DYM50_70);
      QCD1000->Scale(w_QCD1000);
      QCD500->Scale(w_QCD500);
      QCD250->Scale(w_QCD250);
      TT->Scale(w_TT);
      WW->Scale(w_WW);
      WZ->Scale(w_WZ);
      ZZ->Scale(w_ZZ);
      WJets180->Scale(w_WJets180);
      
      DY100->Add(DY70);
      DY100->Add(DYM50_70);
      DY100->Add(DYM50_100);
      QCD1000->Add(QCD250);
      QCD1000->Add(QCD500);
      WW->Add(WZ);
      WW->Add(ZZ);
      QCD1000->Scale(1.9);

      float scaleSig1=(DY100->Integral()+QCD1000->Integral()+WW->Integral()+WJets180->Integral())/ZH1000->Integral();
      float scaleSig2=(DY100->Integral()+QCD1000->Integral()+WW->Integral()+WJets180->Integral())/ZH2500->Integral();
      ZH1000->Scale(scaleSig1);
      ZH2500->Scale(scaleSig2);
      
      DY100->SetFillColor(kBlue-9);
      QCD1000->SetFillColor(kRed-2);
      WW->SetFillColor(kGreen+1);
      TT->SetFillColor(kBlue);
      WJets180->SetFillColor(kMagenta-5);
      ZH1000->SetLineColor(1);
      ZH1500->SetLineColor(2);
      ZH2000->SetLineColor(4);
      ZH2500->SetLineColor(8);
      ZH1000->SetFillColor(1);
      ZH1500->SetFillColor(2);
      ZH2000->SetFillColor(4);
      ZH2500->SetFillColor(8);
      ZH1000->SetFillStyle(3003);
      ZH1500->SetFillStyle(3003);
      ZH2000->SetFillStyle(3003);
      ZH2500->SetFillStyle(3003);
      ZH1000->SetLineWidth(2);
      ZH1500->SetLineWidth(2);
      ZH2000->SetLineWidth(2);
      ZH2500->SetLineWidth(2);
      
      THStack *hs = new THStack("hs","hs");
      hs->SetTitle();
      hs->Add(DY100);
      hs->Add(WW);
      hs->Add(TT);
      hs->Add(WJets180);
      hs->Add(QCD1000);
      
      hs->Draw("histo");
      hs->GetYaxis()->SetTitleSize(0.045);
      hs->GetXaxis()->SetTitleSize(0.045);
      hs->GetYaxis()->SetLabelSize(0.045);
      hs->GetXaxis()->SetLabelSize(0.045);
      hs->GetYaxis()->SetTitleOffset(0.9); 
      hs->GetYaxis()->SetTitle(TString("Events / ")+TString::Format("%.2f",(max-min)/bin));
      hs->GetXaxis()->SetTitle(axis); 
      hs->SetMinimum(0);
      hs->SetMaximum(maxy);
      
      ZH1000->Draw("histo same");
      //ZH1500->Draw("histo same");
      //ZH2000->Draw("histo same");
      ZH2500->Draw("histo same");
      
      TLatex latexLabel1;
      latexLabel1.SetTextSize(0.04);
      latexLabel1.SetNDC();
      latexLabel1.DrawLatex(0.10, 0.91, "CMS");	
      TLatex latexLabel2;
      latexLabel2.SetTextSize(0.04);
      latexLabel2.SetTextFont(42);
      latexLabel2.SetNDC();
      latexLabel2.DrawLatex(0.605, 0.91, "L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
      
      TLegend *pl2 = new TLegend(0.43,0.55,0.89,0.89);
      pl2->SetTextSize(0.03); 
      pl2->SetFillColor(0);
      TLegendEntry *ple2 = pl2->AddEntry(DY100, "Drell-Yan",  "F");
      ple2 = pl2->AddEntry(WW, "Diboson",  "F");
      ple2 = pl2->AddEntry(TT, "ttbar",  "F");
      ple2 = pl2->AddEntry(WJets180, "WJets (pt180)",  "F");
      ple2 = pl2->AddEntry(QCD1000, "QCD",  "F");
      ple2 = pl2->AddEntry(ZH1000, "Signal (M_{X}=1.0TeV - Normalized to bkg yield)",  "F");
      ple2->SetTextSize(0.025); 
      ple2 = pl2->AddEntry(ZH2500, "Signal (M_{X}=2.5TeV - Normalized to bkg yield)",  "F");
      ple2->SetTextSize(0.025); 
      pl2->Draw();
      
      if(log) {
	hs->SetMinimum(0.05);
	c1->SetLogy();
      }
      char Lep1pt      [2];  sprintf(Lep1pt,      "%.0f", lep1PtCut);    TString lep1pt       = Lep1pt      ;
      char Taupt       [2];  sprintf(Taupt,       "%.0f", tauPtCut);	   TString taupt        = Taupt       ;
      char Btag        [2];  sprintf(Btag,        "%.0f", bTagCut);	   TString btag         = Btag        ;
      char Drcut       [2];  sprintf(Drcut,       "%.0f", deltaRCutSL);	   TString drcut        = Drcut       ;
      char Masssvfitcut[2];  sprintf(Masssvfitcut,"%.0f", MassSvfitCut); TString masssvfitcut = Masssvfitcut;
      char Massviscut  [2];  sprintf(Massviscut,  "%.0f", MassVisCutSL); TString massviscut   = Massviscut  ;
      char Met         [2];  sprintf(Met,         "%.0f", METCutSL);	   TString met          = Met         ;
      if(save && i==0 && !SL) 
	c1->SaveAs("massVis_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_FL1000.pdf");
      if(save && i==0 && SL)  
	c1->SaveAs("massVis_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_SL1000.pdf");
    }
  }
}
