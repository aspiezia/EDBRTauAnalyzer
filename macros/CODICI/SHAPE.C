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
  float lep1PtCut=10; int bTagCut=0; float deltaRCut=10; float MassSvfitCut=10000; float METCut = 20; float MassVisCut=0; 
  char demo[500]; sprintf(demo,"demo/TreeSB3");
  TString DEMO = demo;
  char channel[500];  

  vector<string> PLOT;              vector<int> BIN;   vector<float> MIN;  vector<float> MAX;   vector<float> MAXY;   vector<TString> AXIS;
  PLOT.push_back("XMassSVFit");     BIN.push_back(30); MIN.push_back(0);   MAX.push_back(3000); MAXY.push_back(90);   AXIS.push_back("M(Z,H) [GeV] ");

  for(int type=0; type<5; type++){
    if(type==0) sprintf(channel,   "EleMuo"); 
    if(type==1) sprintf(channel,   "MuoMuo"); 
    if(type==2) sprintf(channel,   "EleEle"); 
    if(type==3) sprintf(channel,   "MuoTau"); 
    if(type==4) sprintf(channel,   "EleTau"); 
    char openTree[500]; sprintf(openTree, "%s%s",demo,channel); 

    char *plot = PLOT[0].c_str();
    TString name = PLOT[0];
    int bin=BIN[0]; 
    float min=MIN[0]; 
    float max=MAX[0]; 
    float maxy=MAXY[0];
    TString axis = AXIS[0];
      
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
    
    TFile *fileAA = TFile::Open("../../RISULTATI/analyzer_290514/ZH1000.root");
    TFile *fileBB = TFile::Open("../../RISULTATI/analyzer_290514/ZH1500.root");
    TFile *fileCC = TFile::Open("../../RISULTATI/analyzer_290514/ZH2000.root");
    TFile *fileDD = TFile::Open("../../RISULTATI/analyzer_290514/ZH2500.root");
    TFile *file02 = TFile::Open("../../RISULTATI/analyzer_290514/DY100.root");
    TFile *file03 = TFile::Open("../../RISULTATI/analyzer_290514/DY70.root");
    TFile *file04 = TFile::Open("../../RISULTATI/analyzer_290514/DYM50_100.root");
    TFile *file05 = TFile::Open("../../RISULTATI/analyzer_290514/DYM50_70.root");
    TFile *file06 = TFile::Open("../../RISULTATI/analyzer_290514/QCD1000.root");
    TFile *file07 = TFile::Open("../../RISULTATI/analyzer_290514/QCD250.root");
    TFile *file08 = TFile::Open("../../RISULTATI/analyzer_290514/QCD500.root");
    TFile *file09 = TFile::Open("../../RISULTATI/analyzer_290514/TT.root");
    TFile *file10 = TFile::Open("../../RISULTATI/analyzer_290514/WJets180.root");
    TFile *file11 = TFile::Open("../../RISULTATI/analyzer_290514/WW.root");
    TFile *file12 = TFile::Open("../../RISULTATI/analyzer_290514/WZ.root");
    TFile *file13 = TFile::Open("../../RISULTATI/analyzer_290514/ZZ.root");
    
    char BTAG[1000]; sprintf(BTAG, "");
    if(bTagCut==1) {sprintf(BTAG, " && nbtagsL1<1");}
    if(bTagCut==2) {sprintf(BTAG, " && nbtagsL1<2");}
    if(bTagCut==3) {sprintf(BTAG, " && nbtagsM1<1");}
    if(bTagCut==4) {sprintf(BTAG, " && nbtagsM1<2");}
    if(bTagCut==5) {sprintf(BTAG, " && nbtagsT1<1");}
    if(bTagCut==6) {sprintf(BTAG, " && nbtagsT1<2");}
    char CUTPre[1000]; sprintf(CUTPre, "PUWeight*(trigger==1 && MassSvfit<%f %s",MassSvfitCut,BTAG);
    char CUT [1000]; 
    sprintf(CUT, "%s && met>%f && MassVis>%f && dRLep1Lep2<%f)",CUTPre,METCut,MassVisCut,deltaRCut);
    
    TTree *TreeAA = (TTree*)fileAA->Get(openTree);  char inputAA[50]; sprintf(inputAA, "%s>>hAA(%i,%f,%f)", plot,bin,min,max);
    TTree *TreeBB = (TTree*)fileBB->Get(openTree);  char inputBB[50]; sprintf(inputBB, "%s>>hBB(%i,%f,%f)", plot,bin,min,max);
    TTree *TreeCC = (TTree*)fileCC->Get(openTree);  char inputCC[50]; sprintf(inputCC, "%s>>hCC(%i,%f,%f)", plot,bin,min,max);
    TTree *TreeDD = (TTree*)fileDD->Get(openTree);  char inputDD[50]; sprintf(inputDD, "%s>>hDD(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree02 = (TTree*)file02->Get(openTree);  char input02[50]; sprintf(input02, "%s>>h02(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree03 = (TTree*)file03->Get(openTree);  char input03[50]; sprintf(input03, "%s>>h03(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree04 = (TTree*)file04->Get(openTree);  char input04[50]; sprintf(input04, "%s>>h04(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree05 = (TTree*)file05->Get(openTree);  char input05[50]; sprintf(input05, "%s>>h05(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree06 = (TTree*)file06->Get(openTree);  char input06[50]; sprintf(input06, "%s>>h06(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree07 = (TTree*)file07->Get(openTree);  char input07[50]; sprintf(input07, "%s>>h07(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree08 = (TTree*)file08->Get(openTree);  char input08[50]; sprintf(input08, "%s>>h08(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree09 = (TTree*)file09->Get(openTree);  char input09[50]; sprintf(input09, "%s>>h09(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree10 = (TTree*)file10->Get(openTree);  char input10[50]; sprintf(input10, "%s>>h10(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree11 = (TTree*)file11->Get(openTree);  char input11[50]; sprintf(input11, "%s>>h11(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree12 = (TTree*)file12->Get(openTree);  char input12[50]; sprintf(input12, "%s>>h12(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree13 = (TTree*)file13->Get(openTree);  char input13[50]; sprintf(input13, "%s>>h13(%i,%f,%f)", plot,bin,min,max);
    
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
    
    if(DEMO=="demo/TreeSB1"){
      char CUTPreSB2[1000]; sprintf(CUTPreSB2, "PUWeight*(trigger==1 && MassSvfit<%f && jetMass>150 %s",MassSvfitCut,BTAG);
      char CUTSB2[1000];    sprintf(CUTSB2, "%s && met>%f && MassVis>%f && dRLep1Lep2<%f)",CUTPreSB2,METCut,MassVisCut,deltaRCut);
      char openTreeSB2[500]; sprintf(openTreeSB2, "demo/TreeSB2%s",channel); 
      
      TTree *TreeAASB2 = (TTree*)fileAA->Get(openTreeSB2);  char inputAASB2[50]; sprintf(inputAASB2, "%s>>hAASB2(%i,%f,%f)", plot,bin,min,max);
      TTree *TreeBBSB2 = (TTree*)fileBB->Get(openTreeSB2);  char inputBBSB2[50]; sprintf(inputBBSB2, "%s>>hBBSB2(%i,%f,%f)", plot,bin,min,max);
      TTree *TreeCCSB2 = (TTree*)fileCC->Get(openTreeSB2);  char inputCCSB2[50]; sprintf(inputCCSB2, "%s>>hCCSB2(%i,%f,%f)", plot,bin,min,max);
      TTree *TreeDDSB2 = (TTree*)fileDD->Get(openTreeSB2);  char inputDDSB2[50]; sprintf(inputDDSB2, "%s>>hDDSB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree02SB2 = (TTree*)file02->Get(openTreeSB2);  char input02SB2[50]; sprintf(input02SB2, "%s>>h02SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree03SB2 = (TTree*)file03->Get(openTreeSB2);  char input03SB2[50]; sprintf(input03SB2, "%s>>h03SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree04SB2 = (TTree*)file04->Get(openTreeSB2);  char input04SB2[50]; sprintf(input04SB2, "%s>>h04SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree05SB2 = (TTree*)file05->Get(openTreeSB2);  char input05SB2[50]; sprintf(input05SB2, "%s>>h05SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree06SB2 = (TTree*)file06->Get(openTreeSB2);  char input06SB2[50]; sprintf(input06SB2, "%s>>h06SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree07SB2 = (TTree*)file07->Get(openTreeSB2);  char input07SB2[50]; sprintf(input07SB2, "%s>>h07SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree08SB2 = (TTree*)file08->Get(openTreeSB2);  char input08SB2[50]; sprintf(input08SB2, "%s>>h08SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree09SB2 = (TTree*)file09->Get(openTreeSB2);  char input09SB2[50]; sprintf(input09SB2, "%s>>h09SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree10SB2 = (TTree*)file10->Get(openTreeSB2);  char input10SB2[50]; sprintf(input10SB2, "%s>>h10SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree11SB2 = (TTree*)file11->Get(openTreeSB2);  char input11SB2[50]; sprintf(input11SB2, "%s>>h11SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree12SB2 = (TTree*)file12->Get(openTreeSB2);  char input12SB2[50]; sprintf(input12SB2, "%s>>h12SB2(%i,%f,%f)", plot,bin,min,max);
      TTree *Tree13SB2 = (TTree*)file13->Get(openTreeSB2);  char input13SB2[50]; sprintf(input13SB2, "%s>>h13SB2(%i,%f,%f)", plot,bin,min,max);
      
      TreeAASB2->Draw(inputAASB2,CUTSB2); if(TreeAASB2->Draw(inputAASB2,CUTSB2)) {ZH1000->Add(hAASB2);}
      TreeBBSB2->Draw(inputBBSB2,CUTSB2); if(TreeBBSB2->Draw(inputBBSB2,CUTSB2)) {ZH1500->Add(hBBSB2);}
      TreeCCSB2->Draw(inputCCSB2,CUTSB2); if(TreeCCSB2->Draw(inputCCSB2,CUTSB2)) {ZH2000->Add(hCCSB2);}
      TreeDDSB2->Draw(inputDDSB2,CUTSB2); if(TreeDDSB2->Draw(inputDDSB2,CUTSB2)) {ZH2500->Add(hDDSB2);}
      Tree02SB2->Draw(input02SB2,CUTSB2); if(Tree02SB2->Draw(input02SB2,CUTSB2)) {DY100->Add(h02SB2);}
      Tree03SB2->Draw(input03SB2,CUTSB2); if(Tree03SB2->Draw(input03SB2,CUTSB2)) {DY70->Add(h03SB2);}
      Tree04SB2->Draw(input04SB2,CUTSB2); if(Tree04SB2->Draw(input04SB2,CUTSB2)) {DYM50_100->Add(h04SB2);}
      Tree05SB2->Draw(input05SB2,CUTSB2); if(Tree05SB2->Draw(input05SB2,CUTSB2)) {DYM50_70->Add(h05SB2);}
      Tree06SB2->Draw(input06SB2,CUTSB2); if(Tree06SB2->Draw(input06SB2,CUTSB2)) {QCD1000->Add(h06SB2);}
      Tree07SB2->Draw(input07SB2,CUTSB2); if(Tree07SB2->Draw(input07SB2,CUTSB2)) {QCD250->Add(h07SB2);}
      Tree08SB2->Draw(input08SB2,CUTSB2); if(Tree08SB2->Draw(input08SB2,CUTSB2)) {QCD500->Add(h08SB2);}
      Tree09SB2->Draw(input09SB2,CUTSB2); if(Tree09SB2->Draw(input09SB2,CUTSB2)) {TT->Add(h09SB2);}
      Tree10SB2->Draw(input10SB2,CUTSB2); if(Tree10SB2->Draw(input10SB2,CUTSB2)) {WJets180->Add(h10SB2);}
      Tree11SB2->Draw(input11SB2,CUTSB2); if(Tree11SB2->Draw(input11SB2,CUTSB2)) {WW->Add(h11SB2);}
      Tree12SB2->Draw(input12SB2,CUTSB2); if(Tree12SB2->Draw(input12SB2,CUTSB2)) {WZ->Add(h12SB2);}
      Tree13SB2->Draw(input13SB2,CUTSB2); if(Tree13SB2->Draw(input13SB2,CUTSB2)) {ZZ->Add(h13SB2);}
    }
    
    fileAA->Close(); fileBB->Close(); fileCC->Close(); fileDD->Close();
    file13->Close(); file02->Close(); file03->Close(); file04->Close();
    file05->Close(); file06->Close(); file07->Close(); file08->Close();
    file09->Close(); file10->Close(); file11->Close(); file12->Close();
    fileAA->Delete(); fileBB->Delete(); fileCC->Delete(); fileDD->Delete();
    file13->Delete(); file02->Delete(); file03->Delete(); file04->Delete();
    file05->Delete(); file06->Delete(); file07->Delete(); file08->Delete();
    file09->Delete(); file10->Delete(); file11->Delete(); file12->Delete();
    
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
      
    TLatex latexLabel1;
    latexLabel1.SetTextSize(0.04);
    latexLabel1.SetNDC();
    latexLabel1.DrawLatex(0.10, 0.91, "CMS");	
    TLatex latexLabel2;
    latexLabel2.SetTextSize(0.04);
    latexLabel2.SetTextFont(42);
    latexLabel2.SetNDC();
    latexLabel2.DrawLatex(0.605, 0.91, "L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
      
    TLegend *pl2 = new TLegend(0.55,0.70,0.89,0.89);
    pl2->SetTextSize(0.03); 
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(DY100, "Drell-Yan",  "F");
    ple2 = pl2->AddEntry(WW, "Diboson",  "F");
    ple2 = pl2->AddEntry(TT, "ttbar",  "F");
    ple2 = pl2->AddEntry(WJets180, "WJets (pt180)",  "F");
    ple2 = pl2->AddEntry(QCD1000, "QCD",  "F"); 
    pl2->Draw();
      
    if(log) {
      hs->SetMinimum(0.05);
      c1->SetLogy();
    }
      
    TString CHANNEL = channel;
    char Btag        [2];  sprintf(Btag,        "%.0f", bTagCut);      TString btag         = Btag        ;
    char Drcut       [2];  sprintf(Drcut,       "%.0f", deltaRCut);    TString drcut        = Drcut       ;
    char Masssvfitcut[2];  sprintf(Masssvfitcut,"%.0f", MassSvfitCut); TString masssvfitcut = Masssvfitcut;
    char Massviscut  [2];  sprintf(Massviscut,  "%.0f", MassVisCut);   TString massviscut   = Massviscut  ;
    char Met         [2];  sprintf(Met,         "%.0f", METCut);       TString met          = Met         ;
    char sb[3];
    if(DEMO=="demo/TreeSB1"){ sprintf(sb,"SB1");     TString SB = sb;}
    if(DEMO=="demo/TreeSB2"){ sprintf(sb,"SB2");     TString SB = sb;}
    if(DEMO=="demo/TreeSB3"){ sprintf(sb,"SB3");     TString SB = sb;}
    if(DEMO=="demo/Tree"   ){ sprintf(sb,"SR" );     TString SB = sb;}
    if(save)c1->SaveAs("shape_"+CHANNEL+"_"+SB+"_MET"+met+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+".pdf");
  }
}
