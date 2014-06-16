void Datacard_NOshape(){
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleY(0.96);
  gStyle->SetPaintTextFormat(".2f");
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  
  float XMassWidth=200; float XMassBin=15; float XMassMin=0; float XMassMax=XMassWidth; bool save=true;
  float lep1PtCut=10; float tauPtCut=35; int bTagCut=1; float deltaRCutSL=1; float deltaRCutFL=1; float MassSvfitCut=10000; float METCutFL = 100; float METCutSL = 50; 
  float MassVisCutSL=0;  
  float MassVisCutFL1=0; //ELE-MUO
  float MassVisCutFL2=0; //ELE-ELE e MUO-MUO

  float SigYield  = 0; 
  float BkgYield  = 0;
  float AlphaErr  = 0;
  float Nsideband = 0;
  int WIDTH=150;
  int MIN=850;
  int MAX=1000;
	
  float SignalMass=1000;
  for(int i=0; i<4; i++){
    if(SignalMass==1000){
      WIDTH=150;
      MIN=850;
      MAX=1000;
    } else if(SignalMass==1500){
      WIDTH=250;
      MIN=1250;
      MAX=1500;
    } else if(SignalMass==2000){
      WIDTH=300;
      MIN=1700;
      MAX=2000;
    } else if(SignalMass==2500){
      WIDTH=450;
      MIN=2050;
      MAX=2500;
    }
		
    TH1F *histo_PRE1 = new TH1F("histo_PRE1","histo_PRE1",2, 850,1150); 
    TH1F *ERR1       = new TH1F("ERR1",      "ERR1",      2, 850,1150); 
    THStack *hs1     = new THStack("hs1","hs1");
    TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600); 
    BackgroundEstimation("EleMuo",WIDTH,MIN,MAX,histo_PRE1,ERR1,hs1,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL1,BkgYield,Nsideband,AlphaErr);
    SignalYield("EleMuo",WIDTH,MIN,MAX,SignalMass,save,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL1,SigYield);
    SaveDatacard("EleMuo",BkgYield,Nsideband,SigYield,AlphaErr,SignalMass,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL1);
    SavePlot(c1,MAX,histo_PRE1,ERR1,hs1,"EleMuo",save,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL1,SignalMass);
		
    TH1F *histo_PRE2 = new TH1F("histo_PRE2","histo_PRE2",2, 850,1150); 
    TH1F *ERR2       = new TH1F("ERR2",      "ERR2",      2, 850,1150); 
    THStack *hs2     = new THStack("hs2","hs2");
    TCanvas* c2 = new TCanvas("c2","c2",0,0,800,600); 
    BackgroundEstimation("MuoMuo",WIDTH,MIN,MAX,histo_PRE2,ERR2,hs2,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL2,BkgYield,Nsideband,AlphaErr);
    SignalYield("MuoMuo",WIDTH,MIN,MAX,SignalMass,save,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL2,SigYield);
    SaveDatacard("MuoMuo",BkgYield,Nsideband,SigYield,AlphaErr,SignalMass,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL2);
    SavePlot(c2,MAX,histo_PRE2,ERR2,hs2,"MuoMuo",save,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL2,SignalMass);
		
    TH1F *histo_PRE3 = new TH1F("histo_PRE3","histo_PRE3",2, 850,1150); 
    TH1F *ERR3       = new TH1F("ERR3",      "ERR3",      2, 850,1150); 
    THStack *hs3     = new THStack("hs3","hs3");
    TCanvas* c3 = new TCanvas("c3","c3",0,0,800,600); 
    BackgroundEstimation("EleEle",WIDTH,MIN,MAX,histo_PRE3,ERR3,hs3,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL2,BkgYield,Nsideband,AlphaErr);
    SignalYield("EleEle",WIDTH,MIN,MAX,SignalMass,save,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL2,SigYield);
    SaveDatacard("EleEle",BkgYield,Nsideband,SigYield,AlphaErr,SignalMass,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL2);
    SavePlot(c3,MAX,histo_PRE3,ERR3,hs3,"EleEle",save,lep1PtCut,tauPtCut,bTagCut,METCutFL,deltaRCutFL,MassSvfitCut,MassVisCutFL2,SignalMass);
		
    TH1F *histo_PRE4 = new TH1F("histo_PRE4","histo_PRE4",2, 850,1150); 
    TH1F *ERR4       = new TH1F("ERR4",      "ERR4",      2, 850,1150); 
    THStack *hs4     = new THStack("hs4","hs4");
    TCanvas* c4 = new TCanvas("c4","c4",0,0,800,600); 
    BackgroundEstimation("MuoTau",WIDTH,MIN,MAX,histo_PRE4,ERR4,hs4,lep1PtCut,tauPtCut,bTagCut,METCutSL,deltaRCutSL,MassSvfitCut,MassVisCutSL,BkgYield,Nsideband,AlphaErr);
    SignalYield("MuoTau",WIDTH,MIN,MAX,SignalMass,save,lep1PtCut,tauPtCut,bTagCut,METCutSL,deltaRCutSL,MassSvfitCut,MassVisCutSL,SigYield);
    SaveDatacard("MuoTau",BkgYield,Nsideband,SigYield,AlphaErr,SignalMass,lep1PtCut,tauPtCut,bTagCut,METCutSL,deltaRCutSL,MassSvfitCut,MassVisCutSL);
    SavePlot(c4,MAX,histo_PRE4,ERR4,hs4,"MuoTau",save,lep1PtCut,tauPtCut,bTagCut,METCutSL,deltaRCutSL,MassSvfitCut,MassVisCutSL,SignalMass);
		
    TH1F *histo_PRE5 = new TH1F("histo_PRE5","histo_PRE5",2, 850,1150); 
    TH1F *ERR5       = new TH1F("ERR5",      "ERR5",      2, 850,1150); 
    THStack *hs5     = new THStack("hs5","hs5");
    TCanvas* c5 = new TCanvas("c5","c5",0,0,800,600); 
    BackgroundEstimation("EleTau",WIDTH,MIN,MAX,histo_PRE5,ERR5,hs5,lep1PtCut,tauPtCut,bTagCut,METCutSL,deltaRCutSL,MassSvfitCut,MassVisCutSL,BkgYield,Nsideband,AlphaErr);
    SignalYield("EleTau",WIDTH,MIN,MAX,SignalMass,save,lep1PtCut,tauPtCut,bTagCut,METCutSL,deltaRCutSL,MassSvfitCut,MassVisCutSL,SigYield);
    SaveDatacard("EleTau",BkgYield,Nsideband,SigYield,AlphaErr,SignalMass,lep1PtCut,tauPtCut,bTagCut,METCutSL,deltaRCutSL,MassSvfitCut,MassVisCutSL);
    SavePlot(c5,MAX,histo_PRE5,ERR5,hs5,"EleTau",save,lep1PtCut,tauPtCut,bTagCut,METCutSL,deltaRCutSL,MassSvfitCut,MassVisCutSL,SignalMass);

    SignalMass=SignalMass+500;
  }
}



void SavePlot(TCanvas* c, float MAX, TH1F *histo_PRE, TH1F *ERR, THStack *hs, 
	      char *channel, bool save, float lep1PtCut, float tauPtCut, int bTagCut, float METCut, float deltaRCut, float MassSvfitCut, float MassVisCut, float SignalMass){
  TString CHANNEL = channel; 
	
  histo_PRE->SetLineWidth(2); 
  histo_PRE->SetLineColor(2);
  histo_PRE->SetMarkerColor(2); 
  histo_PRE->SetMarkerStyle(20); 
  histo_PRE->SetMarkerSize(1.3);
  histo_PRE->SetTitle("");
  histo_PRE->SetMinimum(0);
  histo_PRE->GetYaxis()->SetTitleSize(0.045);
  histo_PRE->GetXaxis()->SetTitleSize(0.045);
  histo_PRE->GetYaxis()->SetLabelSize(0.045);
  histo_PRE->GetXaxis()->SetLabelSize(0.045); 
  histo_PRE->GetYaxis()->SetTitle("Events");
  histo_PRE->GetXaxis()->SetTitle("M(Z,H) [GeV]"); 

  ERR->SetFillStyle(3005);
  ERR->SetFillColor(12);
  ERR->SetLineColor(12);

  TH1F *DYFake       = new TH1F("DYFake",       "DYFake",       1,0,10);
  TH1F *TTbarFake    = new TH1F("TTbarFake",    "TTbarFake",    1,0,10);
  TH1F *DIBOSONFake  = new TH1F("DIBOSONFake",  "DIBOSONFake",  1,0,10); 
  TH1F *WJETSFake    = new TH1F("WJETSFake",    "WJETSFake",    1,0,10); 
  TH1F *QCDFake      = new TH1F("QCDFake",      "QCDFake",      1,0,10); 
  DYFake->SetFillColor(kBlue-9);
  TTbarFake->SetFillColor(kBlue);
  QCDFake->SetFillColor(kRed-2);
  DIBOSONFake->SetFillColor(kGreen+1);
  WJETSFake->SetFillColor(kMagenta-5);

  TLegend *pl2 = new TLegend(0.57,0.65,0.89,0.89);
  pl2->SetTextSize(0.035); 
  pl2->SetFillColor(0);
  TLegendEntry *ple2 = pl2->AddEntry(histo_PRE, "BKG prediction",  "LP");
  ple2 = pl2->AddEntry(DYFake, "Drell-Yan",  "F");
  ple2 = pl2->AddEntry(TTbarFake, "TTbar",  "F");
  ple2 = pl2->AddEntry(DIBOSONFake, "Diboson",  "F");
  ple2 = pl2->AddEntry(WJETSFake, "WJets",  "F");
  ple2 = pl2->AddEntry(QCDFake, "QCD",  "F");

  histo_PRE->Draw();	
  hs->Draw("histosame");
  ERR->Draw("E2same");
  histo_PRE->Draw("same");
  pl2->Draw();
	
  if(save){
    char lep1pt[2];      sprintf(lep1pt,      "%.0f", lep1PtCut);
    char taupt[2];       sprintf(taupt,       "%.0f", tauPtCut);
    char btag [2];       sprintf(btag,        "%.0f", bTagCut);
    char drcut [2];      sprintf(drcut,       "%.0f", deltaRCut);
    char masssvfitcut[2];sprintf(masssvfitcut,"%.0f", MassSvfitCut);
    char massviscut[2];  sprintf(massviscut,  "%.0f", MassVisCut);
    char met[2];         sprintf(met,         "%.0f", METCut);
    char SIG[2];         sprintf(SIG,         "%.0f", SignalMass);
    //c->SaveAs("datacard_"+CHANNEL+"_"+SIG+"_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+".pdf");
  }
}


void BackgroundEstimation(char *channel, float XMassWidth, float XMassMin, float XMassMax, TH1F *histo_PRE, TH1F *ERR, THStack *hs, float lep1PtCut, float tauPtCut, 
			  int bTagCut, float METCut, float deltaRCut, float MassSvfitCut, float MassVisCut, float & BkgYield, float & Nsideband, float & AlphaErr){
  char *plot = "met";
  int bin=1; 
  float min=0; 
  float max=50000; 
  char demoSB1    [1000]; sprintf(demoSB1,    "demo/TreeSB1"); 
  char openTreeSB1[1000]; sprintf(openTreeSB1,"%s%s",demoSB1,channel);
  char demoSB2    [1000]; sprintf(demoSB2,    "demo/TreeSB2"); 
  char openTreeSB2[1000]; sprintf(openTreeSB2,"%s%s",demoSB2,channel);
  char demoSB3    [1000]; sprintf(demoSB3,    "demo/TreeSB3"); 
  char openTreeSB3[1000]; sprintf(openTreeSB3,"%s%s",demoSB3,channel);
  char demo       [1000]; sprintf(demo,       "demo/Tree"); 
  char openTree   [1000]; sprintf(openTree,   "%s%s",demo,channel); 
  TString CHANNEL = channel;  
  Nsideband=0;
  AlphaErr=0;
  float Alpha=0;

  TH1F *DY       = new TH1F("DY",       "DY",       2, XMassMin, XMassMax+XMassWidth);
  TH1F *TTbar    = new TH1F("TTbar",    "TTbar",    2, XMassMin, XMassMax+XMassWidth);
  TH1F *DIBOSON  = new TH1F("DIBOSON",  "DIBOSON",  2, XMassMin, XMassMax+XMassWidth); 
  TH1F *WJETS    = new TH1F("WJETS",    "WJETS",    2, XMassMin, XMassMax+XMassWidth); 
  TH1F *QCD      = new TH1F("QCD",      "QCD",      2, XMassMin, XMassMax+XMassWidth);  
	
  TFile *file01=TFile::Open("../../RISULTATI/analyzer_290514/data.root");      
  TFile *file02=TFile::Open("../../RISULTATI/analyzer_290514/DY100.root");     
  TFile *file03=TFile::Open("../../RISULTATI/analyzer_290514/DY70.root");      
  TFile *file04=TFile::Open("../../RISULTATI/analyzer_290514/DYM50_100.root"); 
  TFile *file05=TFile::Open("../../RISULTATI/analyzer_290514/DYM50_70.root");  
  TFile *file06=TFile::Open("../../RISULTATI/analyzer_290514/QCD1000.root");   
  TFile *file07=TFile::Open("../../RISULTATI/analyzer_290514/QCD250.root");    
  TFile *file08=TFile::Open("../../RISULTATI/analyzer_290514/QCD500.root");    
  TFile *file09=TFile::Open("../../RISULTATI/analyzer_290514/TT.root");        
  TFile *file10=TFile::Open("../../RISULTATI/analyzer_290514/WJetsHT.root");  
  TFile *file11=TFile::Open("../../RISULTATI/analyzer_290514/WW.root");         
  TFile *file12=TFile::Open("../../RISULTATI/analyzer_290514/WZ.root");         
  TFile *file13=TFile::Open("../../RISULTATI/analyzer_290514/ZZ.root");        

  TTree *Tree01=(TTree*)file01->Get(openTreeSB1);  TTree *Tree14=(TTree*)file01->Get(openTree); 
  TTree *Tree02=(TTree*)file02->Get(openTreeSB1);  TTree *Tree15=(TTree*)file02->Get(openTree); 
  TTree *Tree03=(TTree*)file03->Get(openTreeSB1);  TTree *Tree16=(TTree*)file03->Get(openTree); 
  TTree *Tree04=(TTree*)file04->Get(openTreeSB1);  TTree *Tree17=(TTree*)file04->Get(openTree); 
  TTree *Tree05=(TTree*)file05->Get(openTreeSB1);  TTree *Tree18=(TTree*)file05->Get(openTree); 
  TTree *Tree06=(TTree*)file06->Get(openTreeSB1);  TTree *Tree19=(TTree*)file06->Get(openTree); 
  TTree *Tree07=(TTree*)file07->Get(openTreeSB1);  TTree *Tree20=(TTree*)file07->Get(openTree); 
  TTree *Tree08=(TTree*)file08->Get(openTreeSB1);  TTree *Tree21=(TTree*)file08->Get(openTree); 
  TTree *Tree09=(TTree*)file09->Get(openTreeSB1);  TTree *Tree22=(TTree*)file09->Get(openTree); 
  TTree *Tree10=(TTree*)file10->Get(openTreeSB1);  TTree *Tree23=(TTree*)file10->Get(openTree); 
  TTree *Tree11=(TTree*)file11->Get(openTreeSB1);  TTree *Tree24=(TTree*)file11->Get(openTree); 
  TTree *Tree12=(TTree*)file12->Get(openTreeSB1);  TTree *Tree25=(TTree*)file12->Get(openTree); 
  TTree *Tree13=(TTree*)file13->Get(openTreeSB1);  TTree *Tree26=(TTree*)file13->Get(openTree); 
  
  TTree *Tree27=(TTree*)file01->Get(openTreeSB2);  TTree *Tree40=(TTree*)file01->Get(openTreeSB3); 
  TTree *Tree28=(TTree*)file02->Get(openTreeSB2);  TTree *Tree41=(TTree*)file02->Get(openTreeSB3); 
  TTree *Tree29=(TTree*)file03->Get(openTreeSB2);  TTree *Tree42=(TTree*)file03->Get(openTreeSB3); 
  TTree *Tree30=(TTree*)file04->Get(openTreeSB2);  TTree *Tree43=(TTree*)file04->Get(openTreeSB3); 
  TTree *Tree31=(TTree*)file05->Get(openTreeSB2);  TTree *Tree44=(TTree*)file05->Get(openTreeSB3); 
  TTree *Tree32=(TTree*)file06->Get(openTreeSB2);  TTree *Tree45=(TTree*)file06->Get(openTreeSB3); 
  TTree *Tree33=(TTree*)file07->Get(openTreeSB2);  TTree *Tree46=(TTree*)file07->Get(openTreeSB3); 
  TTree *Tree34=(TTree*)file08->Get(openTreeSB2);  TTree *Tree47=(TTree*)file08->Get(openTreeSB3); 
  TTree *Tree35=(TTree*)file09->Get(openTreeSB2);  TTree *Tree48=(TTree*)file09->Get(openTreeSB3); 
  TTree *Tree36=(TTree*)file10->Get(openTreeSB2);  TTree *Tree49=(TTree*)file10->Get(openTreeSB3); 
  TTree *Tree37=(TTree*)file11->Get(openTreeSB2);  TTree *Tree50=(TTree*)file11->Get(openTreeSB3); 
  TTree *Tree38=(TTree*)file12->Get(openTreeSB2);  TTree *Tree51=(TTree*)file12->Get(openTreeSB3); 
  TTree *Tree39=(TTree*)file13->Get(openTreeSB2);  TTree *Tree52=(TTree*)file13->Get(openTreeSB3); 
	
  for(int i=0; i<2; i++){
    TH1F *data      = new TH1F("","",bin,min,max);       TH1F *data_SR      = new TH1F("","",bin,min,max);
    TH1F *DY100     = new TH1F("","",bin,min,max);	 TH1F *DY100_SR     = new TH1F("","",bin,min,max);
    TH1F *DY70      = new TH1F("","",bin,min,max);	 TH1F *DY70_SR      = new TH1F("","",bin,min,max);
    TH1F *DYM50_100 = new TH1F("","",bin,min,max);	 TH1F *DYM50_100_SR = new TH1F("","",bin,min,max);
    TH1F *DYM50_70  = new TH1F("","",bin,min,max);	 TH1F *DYM50_70_SR  = new TH1F("","",bin,min,max);
    TH1F *QCD1000   = new TH1F("","",bin,min,max);	 TH1F *TT_SR        = new TH1F("","",bin,min,max);
    TH1F *QCD250    = new TH1F("","",bin,min,max);	 TH1F *QCD1000_SR   = new TH1F("","",bin,min,max);
    TH1F *QCD500    = new TH1F("","",bin,min,max);	 TH1F *QCD250_SR    = new TH1F("","",bin,min,max);
    TH1F *TT        = new TH1F("","",bin,min,max);	 TH1F *QCD500_SR    = new TH1F("","",bin,min,max);
    TH1F *WJetsHT   = new TH1F("","",bin,min,max);	 TH1F *WJetsHT_SR   = new TH1F("","",bin,min,max);
    TH1F *WW        = new TH1F("","",bin,min,max);	 TH1F *WW_SR        = new TH1F("","",bin,min,max);
    TH1F *WZ        = new TH1F("","",bin,min,max);	 TH1F *WZ_SR        = new TH1F("","",bin,min,max);
    TH1F *ZZ        = new TH1F("","",bin,min,max);	 TH1F *ZZ_SR        = new TH1F("","",bin,min,max);
		
    char BTAG[1000]; sprintf(BTAG, "");
    if(bTagCut==1) {sprintf(BTAG, " && nbtagsL1<1");}
    if(bTagCut==2) {sprintf(BTAG, " && nbtagsL1<2");}
    if(bTagCut==3) {sprintf(BTAG, " && nbtagsM1<1");}
    if(bTagCut==4) {sprintf(BTAG, " && nbtagsM1<2");}
    if(bTagCut==5) {sprintf(BTAG, " && nbtagsT1<1");}
    if(bTagCut==6) {sprintf(BTAG, " && nbtagsT1<2");}
    
    //PLOTS IN SIDEBAND1
    char CUTPre[1000]; sprintf(CUTPre, "PUWeight*(trigger==1 && PtSvfit>100 && XMassSVFit>%f && XMassSVFit<%f && met>%f && dRLep1Lep2<%f && MassSvfit<%f && MassVis>%f %s",
			       XMassMin,XMassMax,METCut,deltaRCut,MassSvfitCut,MassVisCut,BTAG);
    char CUT[1000]; 
    if(CHANNEL=="EleMuo")     sprintf(CUT, "%s )",CUTPre);
    else if(CHANNEL=="MuoMuo")sprintf(CUT, "%s && EleMuo==0 && lep1Pt>%f)",CUTPre,lep1PtCut);
    else if(CHANNEL=="EleEle")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && lep1Pt>%f)",CUTPre,lep1PtCut);
    else if(CHANNEL=="EleTau")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && EleEle==0 && lep1Pt>%f)",CUTPre,tauPtCut);
    else if(CHANNEL=="MuoTau")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && EleEle==0 && EleTau==0 && lep1Pt>%f)",CUTPre,tauPtCut);
    char input01[50];sprintf(input01,  "%s>>h01(%i,%f,%f)", plot,bin,min,max);      Tree01->Draw(input01,CUT,"E"); if(Tree01->Draw(input01,CUT,"E")) {data      = h01; }
    char input02[50];sprintf(input02,  "%s>>h02(%i,%f,%f)", plot,bin,min,max);      Tree02->Draw(input02,CUT);     if(Tree02->Draw(input02,CUT))     {DY100     = h02; }
    char input03[50];sprintf(input03,  "%s>>h03(%i,%f,%f)", plot,bin,min,max);      Tree03->Draw(input03,CUT);     if(Tree03->Draw(input03,CUT))     {DY70      = h03; }
    char input04[50];sprintf(input04,  "%s>>h04(%i,%f,%f)", plot,bin,min,max);      Tree04->Draw(input04,CUT);     if(Tree04->Draw(input04,CUT))     {DYM50_100 = h04; }
    char input05[50];sprintf(input05,  "%s>>h05(%i,%f,%f)", plot,bin,min,max);      Tree05->Draw(input05,CUT);     if(Tree05->Draw(input05,CUT))     {DYM50_70  = h05; }
    char input06[50];sprintf(input06,  "%s>>h06(%i,%f,%f)", plot,bin,min,max);      Tree06->Draw(input06,CUT);     if(Tree06->Draw(input06,CUT))     {QCD1000   = h06; }
    char input07[50];sprintf(input07,  "%s>>h07(%i,%f,%f)", plot,bin,min,max);      Tree07->Draw(input07,CUT);     if(Tree07->Draw(input07,CUT))     {QCD250    = h07; }
    char input08[50];sprintf(input08,  "%s>>h08(%i,%f,%f)", plot,bin,min,max);      Tree08->Draw(input08,CUT);     if(Tree08->Draw(input08,CUT))     {QCD500    = h08; }
    char input09[50];sprintf(input09,  "%s>>h09(%i,%f,%f)", plot,bin,min,max);      Tree09->Draw(input09,CUT);     if(Tree09->Draw(input09,CUT))     {TT        = h09; }
    char input10[50];sprintf(input10,  "%s>>h10(%i,%f,%f)", plot,bin,min,max);      Tree10->Draw(input10,CUT);     if(Tree10->Draw(input10,CUT))     {WJetsHT   = h10; }
    char input11[50];sprintf(input11,  "%s>>h11(%i,%f,%f)", plot,bin,min,max);      Tree11->Draw(input11,CUT);     if(Tree11->Draw(input11,CUT))     {WW        = h11; }
    char input12[50];sprintf(input12,  "%s>>h12(%i,%f,%f)", plot,bin,min,max);      Tree12->Draw(input12,CUT);     if(Tree12->Draw(input12,CUT))     {WZ        = h12; }
    char input13[50];sprintf(input13,  "%s>>h13(%i,%f,%f)", plot,bin,min,max);      Tree13->Draw(input13,CUT);     if(Tree13->Draw(input13,CUT))     {WZ        = h13; }

    //PLOTS IN SIDEBAND2
    sprintf(CUTPre, "PUWeight*(trigger==1 && PtSvfit>100 && jetMass>150 && XMassSVFit>%f && XMassSVFit<%f && met>%f && dRLep1Lep2<%f && MassSvfit<%f && MassVis>%f %s",
	    XMassMin,XMassMax,METCut,deltaRCut,MassSvfitCut,MassVisCut,BTAG);
    if(CHANNEL=="EleMuo")     sprintf(CUT, "%s )",CUTPre);
    else if(CHANNEL=="MuoMuo")sprintf(CUT, "%s && EleMuo==0 && lep1Pt>%f)",CUTPre,lep1PtCut);
    else if(CHANNEL=="EleEle")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && lep1Pt>%f)",CUTPre,lep1PtCut);
    else if(CHANNEL=="EleTau")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && EleEle==0 && lep1Pt>%f)",CUTPre,tauPtCut);
    else if(CHANNEL=="MuoTau")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && EleEle==0 && EleTau==0 && lep1Pt>%f)",CUTPre,tauPtCut);
    char input27[50];sprintf(input27,  "%s>>h27(%i,%f,%f)", plot,bin,min,max);      Tree27->Draw(input27,CUT,"E"); if(Tree27->Draw(input27,CUT,"E")) {data->Add(h27); }
    char input28[50];sprintf(input28,  "%s>>h28(%i,%f,%f)", plot,bin,min,max);      Tree28->Draw(input28,CUT);     if(Tree28->Draw(input28,CUT))     {DY100->Add(h28); }
    char input29[50];sprintf(input29,  "%s>>h29(%i,%f,%f)", plot,bin,min,max);      Tree29->Draw(input29,CUT);     if(Tree29->Draw(input29,CUT))     {DY70->Add(h29); }
    char input30[50];sprintf(input30,  "%s>>h30(%i,%f,%f)", plot,bin,min,max);      Tree30->Draw(input30,CUT);     if(Tree30->Draw(input30,CUT))     {DYM50_100->Add(h30); }
    char input31[50];sprintf(input31,  "%s>>h31(%i,%f,%f)", plot,bin,min,max);      Tree31->Draw(input31,CUT);     if(Tree31->Draw(input31,CUT))     {DYM50_70->Add(h31); }
    char input32[50];sprintf(input32,  "%s>>h32(%i,%f,%f)", plot,bin,min,max);      Tree32->Draw(input32,CUT);     if(Tree32->Draw(input32,CUT))     {QCD1000->Add(h32); }
    char input33[50];sprintf(input33,  "%s>>h33(%i,%f,%f)", plot,bin,min,max);      Tree33->Draw(input33,CUT);     if(Tree33->Draw(input33,CUT))     {QCD250->Add(h33); }
    char input34[50];sprintf(input34,  "%s>>h34(%i,%f,%f)", plot,bin,min,max);      Tree34->Draw(input34,CUT);     if(Tree34->Draw(input34,CUT))     {QCD500->Add(h34); }
    char input35[50];sprintf(input35,  "%s>>h35(%i,%f,%f)", plot,bin,min,max);      Tree35->Draw(input35,CUT);     if(Tree35->Draw(input35,CUT))     {TT->Add(h35); }
    char input36[50];sprintf(input36,  "%s>>h36(%i,%f,%f)", plot,bin,min,max);      Tree36->Draw(input36,CUT);     if(Tree36->Draw(input36,CUT))     {WJetsHT->Add(h36); }
    char input37[50];sprintf(input37,  "%s>>h37(%i,%f,%f)", plot,bin,min,max);      Tree37->Draw(input37,CUT);     if(Tree37->Draw(input37,CUT))     {WW->Add(h37); }
    char input38[50];sprintf(input38,  "%s>>h38(%i,%f,%f)", plot,bin,min,max);      Tree38->Draw(input38,CUT);     if(Tree38->Draw(input38,CUT))     {WZ->Add(h38); }
    char input39[50];sprintf(input39,  "%s>>h39(%i,%f,%f)", plot,bin,min,max);      Tree39->Draw(input39,CUT);     if(Tree39->Draw(input39,CUT))     {WZ->Add(h39); }

    //PLOTS IN SIDEBAND3
    sprintf(CUTPre, "PUWeight*(trigger==1 && PtSvfit>100 && XMassSVFit>%f && XMassSVFit<%f && met>%f && dRLep1Lep2<%f && MassSvfit<%f && MassVis>%f %s",
	    XMassMin,XMassMax,METCut,deltaRCut,MassSvfitCut,MassVisCut,BTAG);
    if(CHANNEL=="EleMuo")     sprintf(CUT, "%s )",CUTPre);
    else if(CHANNEL=="MuoMuo")sprintf(CUT, "%s && EleMuo==0 && lep1Pt>%f)",CUTPre,lep1PtCut);
    else if(CHANNEL=="EleEle")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && lep1Pt>%f)",CUTPre,lep1PtCut);
    else if(CHANNEL=="EleTau")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && EleEle==0 && lep1Pt>%f)",CUTPre,tauPtCut);
    else if(CHANNEL=="MuoTau")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && EleEle==0 && EleTau==0 && lep1Pt>%f)",CUTPre,tauPtCut);
    char input40[50];sprintf(input40,  "%s>>h40(%i,%f,%f)", plot,bin,min,max);      Tree40->Draw(input40,CUT,"E"); if(Tree40->Draw(input40,CUT,"E")) {data->Add(h40); }
    char input41[50];sprintf(input41,  "%s>>h41(%i,%f,%f)", plot,bin,min,max);      Tree41->Draw(input41,CUT);     if(Tree41->Draw(input41,CUT))     {DY100->Add(h41); }
    char input42[50];sprintf(input42,  "%s>>h42(%i,%f,%f)", plot,bin,min,max);      Tree42->Draw(input42,CUT);     if(Tree42->Draw(input42,CUT))     {DY70->Add(h42); }
    char input43[50];sprintf(input43,  "%s>>h43(%i,%f,%f)", plot,bin,min,max);      Tree43->Draw(input43,CUT);     if(Tree43->Draw(input43,CUT))     {DYM50_100->Add(h43); }
    char input44[50];sprintf(input44,  "%s>>h44(%i,%f,%f)", plot,bin,min,max);      Tree44->Draw(input44,CUT);     if(Tree44->Draw(input44,CUT))     {DYM50_70->Add(h44); }
    char input45[50];sprintf(input45,  "%s>>h45(%i,%f,%f)", plot,bin,min,max);      Tree45->Draw(input45,CUT);     if(Tree45->Draw(input45,CUT))     {QCD1000->Add(h45); }
    char input46[50];sprintf(input46,  "%s>>h46(%i,%f,%f)", plot,bin,min,max);      Tree46->Draw(input46,CUT);     if(Tree46->Draw(input46,CUT))     {QCD250->Add(h46); }
    char input47[50];sprintf(input47,  "%s>>h47(%i,%f,%f)", plot,bin,min,max);      Tree47->Draw(input47,CUT);     if(Tree47->Draw(input47,CUT))     {QCD500->Add(h47); }
    char input48[50];sprintf(input48,  "%s>>h48(%i,%f,%f)", plot,bin,min,max);      Tree48->Draw(input48,CUT);     if(Tree48->Draw(input48,CUT))     {TT->Add(h48); }
    char input49[50];sprintf(input49,  "%s>>h49(%i,%f,%f)", plot,bin,min,max);      Tree49->Draw(input49,CUT);     if(Tree49->Draw(input49,CUT))     {WJetsHT->Add(h49); }
    char input50[50];sprintf(input50,  "%s>>h50(%i,%f,%f)", plot,bin,min,max);      Tree50->Draw(input50,CUT);     if(Tree50->Draw(input50,CUT))     {WW->Add(h50); }
    char input51[50];sprintf(input51,  "%s>>h51(%i,%f,%f)", plot,bin,min,max);      Tree51->Draw(input51,CUT);     if(Tree51->Draw(input51,CUT))     {WZ->Add(h51); }
    char input52[50];sprintf(input52,  "%s>>h52(%i,%f,%f)", plot,bin,min,max);      Tree52->Draw(input52,CUT);     if(Tree52->Draw(input52,CUT))     {WZ->Add(h52); }
		
    //PLOTS IN SIGNAL REGION
    char input14 [50];   sprintf(input14, "%s>>h14(%i,%f,%f)", plot,bin,min,max);
    char input15 [50];   sprintf(input15, "%s>>h15(%i,%f,%f)", plot,bin,min,max);
    char input16 [50];   sprintf(input16, "%s>>h16(%i,%f,%f)", plot,bin,min,max);
    char input17 [50];   sprintf(input17, "%s>>h17(%i,%f,%f)", plot,bin,min,max);
    char input18 [50];   sprintf(input18, "%s>>h18(%i,%f,%f)", plot,bin,min,max);
    char input19 [50];   sprintf(input19, "%s>>h19(%i,%f,%f)", plot,bin,min,max);
    char input20 [50];   sprintf(input20, "%s>>h20(%i,%f,%f)", plot,bin,min,max);
    char input21 [50];   sprintf(input21, "%s>>h21(%i,%f,%f)", plot,bin,min,max);
    char input22 [50];   sprintf(input22, "%s>>h22(%i,%f,%f)", plot,bin,min,max);
    char input23 [50];   sprintf(input23, "%s>>h23(%i,%f,%f)", plot,bin,min,max);
    char input24 [50];   sprintf(input24, "%s>>h24(%i,%f,%f)", plot,bin,min,max);
    char input25 [50];   sprintf(input25, "%s>>h25(%i,%f,%f)", plot,bin,min,max);
    char input26 [50];   sprintf(input26, "%s>>h26(%i,%f,%f)", plot,bin,min,max);
    Tree14->Draw(input14,CUT,"E"); if(Tree14->Draw(input14,CUT,"E")) {data_SR      = h14; }
    Tree15->Draw(input15,CUT);     if(Tree15->Draw(input15,CUT))     {DY100_SR     = h15; }
    Tree16->Draw(input16,CUT);     if(Tree16->Draw(input16,CUT))     {DY70_SR      = h16; }
    Tree17->Draw(input17,CUT);     if(Tree17->Draw(input17,CUT))     {DYM50_100_SR = h17; }
    Tree18->Draw(input18,CUT);     if(Tree18->Draw(input18,CUT))     {DYM50_70_SR  = h18; }
    Tree19->Draw(input19,CUT);     if(Tree19->Draw(input19,CUT))     {QCD1000_SR   = h19; }
    Tree20->Draw(input20,CUT);     if(Tree20->Draw(input20,CUT))     {QCD250_SR    = h20; }
    Tree21->Draw(input21,CUT);     if(Tree21->Draw(input21,CUT))     {QCD500_SR    = h21; }
    Tree22->Draw(input22,CUT);     if(Tree22->Draw(input22,CUT))     {TT_SR        = h22; }
    Tree23->Draw(input23,CUT);     if(Tree23->Draw(input23,CUT))     {WJetsHT_SR   = h23; }
    Tree24->Draw(input24,CUT);     if(Tree24->Draw(input24,CUT))     {WW_SR        = h24; }
    Tree25->Draw(input25,CUT);     if(Tree25->Draw(input25,CUT))     {WZ_SR        = h25; }
    Tree26->Draw(input26,CUT);     if(Tree26->Draw(input26,CUT))     {WZ_SR        = h26; }
		
    float w_DY100     = ( 39.100*19702./12511326.);
    float w_DY70      = ( 62.900*19702./11764538.);
    float w_DYM50_100 = (  4.220*19702./4146124.0);
    float w_DYM50_70  = ( 11.050*19702./5389313.0);
    float w_TT        = (225.197*19702./21675970.);
    float w_QCD1000   = (204.000*19702./13843863.);
    float w_QCD500    = (8426.00*19702./30599292.);
    float w_QCD250    = (276000.*19702./27062078.);
    float w_WW        = (57.1097*19702./10000431.);
    float w_WZ        = ( 33.210*19702./10000283.);
    float w_ZZ        = (  8.059*19702./9799908.0);
    float w_WJetsHT   = ( 25.220*19702./4971847.0);
		
    double N_sb     = data->Integral();
    double N_sb_err = sqrt(data->Integral());
    double f_sb_err = 0;  
    double f_sb = 0; 
    double alpha = 0; 
    double alpha_err = 0;

    double N_QCD_SR      = 1.9*(w_QCD1000*QCD1000_SR->Integral() + w_QCD250*QCD250_SR->Integral() + w_QCD500*QCD500_SR->Integral());
    double N_QCD_SR_err  = sqrt(1.9*1.9*(w_QCD1000*w_QCD1000*QCD1000_SR->Integral()+w_QCD250*w_QCD250*QCD250_SR->Integral()+w_QCD500*w_QCD500*QCD500_SR->Integral()));
    double N_QCD_SB      = 1.9*(w_QCD1000*QCD1000->Integral()  +  w_QCD250*QCD250->Integral()  +  w_QCD500*QCD500->Integral());
    double N_QCD_SB_err  = sqrt(1.9*1.9*(w_QCD1000*w_QCD1000*QCD1000->Integral()  +  w_QCD250*w_QCD250*QCD250->Integral()  +  w_QCD500*w_QCD500*QCD500->Integral()));

    double N_DY_SR      = w_DY100*DY100_SR->Integral() + w_DY70*DY70_SR->Integral() + w_DYM50_100*DYM50_100_SR->Integral() + w_DYM50_70*DYM50_70_SR->Integral();
    double N_DY_SR_err  = sqrt(w_DY100*w_DY100*DY100_SR->Integral()+w_DY70*w_DY70*DY70_SR->Integral()+w_DYM50_100*w_DYM50_100*DYM50_100_SR->Integral()
			       +w_DYM50_70*w_DYM50_70*DYM50_70_SR->Integral());
    double N_DY_SB      = w_DY100*DY100->Integral() + w_DY70*DY70->Integral() + w_DYM50_100*DYM50_100->Integral() + w_DYM50_70*DYM50_70->Integral();
    double N_DY_SB_err  = sqrt(w_DY100*w_DY100*DY100->Integral()+w_DY70*w_DY70*DY70->Integral()+w_DYM50_100*w_DYM50_100*DYM50_100->Integral()
			       +w_DYM50_70*w_DYM50_70*DYM50_70->Integral());

    double N_VV_SR      = w_WZ*WZ_SR->Integral() + w_WW*WW_SR->Integral() + w_ZZ*ZZ_SR->Integral();
    double N_VV_SR_err  = sqrt(w_WZ*w_WZ*WZ_SR->Integral() + w_WW*w_WW*WW_SR->Integral() + w_ZZ*w_ZZ*ZZ_SR->Integral());
    double N_VV_SB      = w_WZ*WZ->Integral() + w_WW*WW->Integral() + w_ZZ*ZZ->Integral();
    double N_VV_SB_err  = sqrt(w_WZ*w_WZ*WZ->Integral() + w_WW*w_WW*WW->Integral() + w_ZZ*w_ZZ*ZZ->Integral());

    double N_TT_SR      = w_TT*TT_SR->Integral();
    double N_TT_SR_err  = sqrt(w_TT*w_TT*TT_SR->Integral());
    double N_TT_SB      = w_TT*TT->Integral();
    double N_TT_SB_err  = sqrt(w_TT*w_TT*TT->Integral());

    double N_Wjets_SR      = w_WJetsHT*WJetsHT_SR->Integral();
    double N_Wjets_SR_err  = sqrt(w_WJetsHT*w_WJetsHT*WJetsHT_SR->Integral());
    double N_Wjets_SB      = w_WJetsHT*WJetsHT->Integral();
    double N_Wjets_SB_err  = sqrt(w_WJetsHT*w_WJetsHT*WJetsHT->Integral());
		
    double alpha_num     = N_QCD_SR + N_DY_SR + N_VV_SR + N_TT_SR + N_Wjets_SR;
    double alpha_num_err = sqrt(N_QCD_SR_err*N_QCD_SR_err + N_DY_SR_err*N_DY_SR_err + N_VV_SR_err*N_VV_SR_err + N_TT_SR_err*N_TT_SR_err + N_Wjets_SR_err*N_Wjets_SR_err);
    double alpha_den     = N_QCD_SB + N_DY_SB + N_VV_SB + N_TT_SB + N_Wjets_SB;
    double alpha_den_err = sqrt(N_QCD_SB_err*N_QCD_SB_err + N_DY_SB_err*N_DY_SB_err + N_VV_SB_err*N_VV_SB_err + N_TT_SB_err*N_TT_SB_err + N_Wjets_SB_err*N_Wjets_SB_err);
		
    if(alpha_den!=0) {
      alpha = alpha_num/alpha_den;
      alpha_err = sqrt(alpha_den*alpha_den*alpha_num_err*alpha_num_err + alpha_num*alpha_num*alpha_den_err*alpha_den_err)/(alpha_den*alpha_den);
    }
    Alpha=alpha/2;
    AlphaErr=AlphaErr+alpha_err*alpha_err/4;
		
    double N_PRE = N_sb * alpha;
    double N_PRE_err = sqrt(N_sb_err*N_sb_err*alpha*alpha + alpha_err*alpha_err*N_sb*N_sb);
    Nsideband=Nsideband+N_sb;

    histo_PRE->SetBinContent( i+1, N_PRE);
    histo_PRE->SetBinError(   i+1, N_PRE_err);
    DY->SetBinContent(        i+1, N_DY_SR);
    TTbar->SetBinContent(     i+1, N_TT_SR);
    DIBOSON->SetBinContent(   i+1, N_VV_SR);
    WJETS->SetBinContent(     i+1, N_Wjets_SR);
    QCD->SetBinContent(       i+1, N_QCD_SR);
    ERR->SetBinContent(       i+1, alpha_num);
    ERR->SetBinError(         i+1, alpha_num_err);
    
    XMassMin=XMassMin+XMassWidth;
    XMassMax=XMassMax+XMassWidth;
  }

  double N1_err = 0.; double N1    = histo_PRE->IntegralAndError(1,2,N1_err);
  double N2_err = 0.; double N2    = ERR->IntegralAndError(1,2,N2_err);
  BkgYield=N1;
  if(Alpha!=0) AlphaErr=sqrt(AlphaErr)/Alpha;
  else AlphaErr=0;
  
  cout<<endl;
  cout<<CHANNEL<<endl;
  cout<<"Number of predicted BKG events = "<<N1<<"+/-"<<N1_err<<endl;
  cout<<"Number of BKG events from MC   = "<<N2<<"+/-"<<N2_err<<endl;

  DY->SetFillColor(kBlue-9);
  TTbar->SetFillColor(kBlue);
  QCD->SetFillColor(kRed-2);
  DIBOSON->SetFillColor(kGreen+1);
  WJETS->SetFillColor(kMagenta-5);
  
  hs->Add(DY);
  hs->Add(TTbar);
  hs->Add(DIBOSON);
  hs->Add(WJETS);
  hs->Add(QCD);
}

void SignalYield(char *channel, float XMassWidth, float XMassMin, float XMassMax, int SignalMass, bool save, float lep1PtCut, float tauPtCut, 
		 int bTagCut, float METCut, float deltaRCut, float MassSvfitCut, float MassVisCut, float & SigYield){ 
  char demo       [500]; sprintf(demo,       "demo/Tree"); 
  char openTree   [500]; sprintf(openTree,   "%s%s",demo,channel);  
  TTree *TreeSig;
  double WeightSig = 1;
  if(SignalMass==1000)      {
    TFile *fileSig = TFile::Open("../../RISULTATI/analyzer_290514/ZH1000.root");  
    TreeSig=(TTree*)fileSig->Get(openTree); 
    WeightSig = 19.700/17958;
  } else if(SignalMass==1500) {
    TFile *fileSig = TFile::Open("../../RISULTATI/analyzer_290514/ZH1500.root");  
    TreeSig=(TTree*)fileSig->Get(openTree); 
    WeightSig = 19.700/14099;
  } else if(SignalMass==2000) {
    TFile *fileSig = TFile::Open("../../RISULTATI/analyzer_290514/ZH2000.root");  
    TreeSig=(TTree*)fileSig->Get(openTree); 
    WeightSig = 19.700/13080;
  } else if(SignalMass==2500) {
    TFile *fileSig = TFile::Open("../../RISULTATI/analyzer_290514/ZH2500.root");  
    TreeSig=(TTree*)fileSig->Get(openTree); 
    WeightSig = 19.700/13091;
  }  
  
  TString CHANNEL = channel; 
  char BTAG[1000]; sprintf(BTAG, "");
  if(bTagCut==1) {sprintf(BTAG, " && nbtagsL1<1");}
  if(bTagCut==2) {sprintf(BTAG, " && nbtagsL1<2");}
  if(bTagCut==3) {sprintf(BTAG, " && nbtagsM1<1");}
  if(bTagCut==4) {sprintf(BTAG, " && nbtagsM1<2");}
  if(bTagCut==5) {sprintf(BTAG, " && nbtagsT1<1");}
  if(bTagCut==6) {sprintf(BTAG, " && nbtagsT1<2");}
  char CUTPre[1000]; sprintf(CUTPre, "trigger==1 && PtSvfit>100 && met>%f && dRLep1Lep2<%f && MassSvfit<%f && MassVis>%f %s",
			     METCut,deltaRCut,MassSvfitCut,MassVisCut,BTAG);
  char CUT[1000]; 
  if(CHANNEL=="EleMuo")     sprintf(CUT, "%s ",CUTPre);
  else if(CHANNEL=="MuoMuo")sprintf(CUT, "%s && EleMuo==0 && lep1Pt>%f",CUTPre,lep1PtCut);
  else if(CHANNEL=="EleEle")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && lep1Pt>%f",CUTPre,lep1PtCut);
  else if(CHANNEL=="EleTau")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && EleEle==0 && lep1Pt>%f",CUTPre,tauPtCut);
  else if(CHANNEL=="MuoTau")sprintf(CUT, "%s && EleMuo==0 && MuoMuo==0 && EleEle==0 && EleTau==0 && lep1Pt>%f",CUTPre,tauPtCut);
 
  char CUT1[500]; sprintf(CUT1,  "XMassSVFit>>data_sig(2,%f,%f)",XMassMin,XMassMax+XMassWidth);
  char CUT2[500]; sprintf(CUT2,  "PUWeight*(%s)",CUT);
  TreeSig->SetWeight(WeightSig);
  TreeSig->Draw(CUT1,CUT2);  
  SigYield = data_sig->Integral(); 
}


void SaveDatacard(char *channel, float BkgYield, float NSideband, float SigYield, float AlphaErr, int SignalMass, float lep1PtCut, float tauPtCut, 
		  int bTagCut, float METCut, float deltaRCut, float MassSvfitCut, float MassVisCut){
  char lep1pt[2];      sprintf(lep1pt,      "%.0f", lep1PtCut);
  char taupt[2];       sprintf(taupt,       "%.0f", tauPtCut);
  char btag [2];       sprintf(btag,        "%.0f", bTagCut);
  char drcut [2];      sprintf(drcut,       "%.0f", deltaRCut);
  char masssvfitcut[2];sprintf(masssvfitcut,"%.0f", MassSvfitCut);
  char massviscut[2];  sprintf(massviscut,  "%.0f", MassVisCut);
  char met[2];         sprintf(met,         "%.0f", METCut);
  char saveName [500]; sprintf(saveName, "datacard_%s_%i_MET%s_lep1Pt%s_tauPt%s_btag%s_dR%s_massSvfit%s_massVis%s.txt",
			       channel,SignalMass,met,lep1pt,taupt,btag,drcut,masssvfitcut,massviscut);

  //SYSTEMATICS HERE
  float TauIDStandard = 0;
  if(channel=="MuoTau" && SignalMass==1000) TauIDStandard = 1 + 0.089;
  if(channel=="EleTau" && SignalMass==1000) TauIDStandard = 1 + 0.085;
  if(channel=="MuoTau" && SignalMass==1500) TauIDStandard = 1 + 0.099;
  if(channel=="EleTau" && SignalMass==1500) TauIDStandard = 1 + 0.097;
  if(channel=="MuoTau" && SignalMass==2000) TauIDStandard = 1 + 0.111;
  if(channel=="EleTau" && SignalMass==2000) TauIDStandard = 1 + 0.108;
  if(channel=="MuoTau" && SignalMass==2500) TauIDStandard = 1 + 0.124;
  if(channel=="EleTau" && SignalMass==2500) TauIDStandard = 1 + 0.119;
  float TauCleaning   = 0;
  if(channel=="MuoTau" && SignalMass==1000) TauCleaning = 1 + 0.100;
  if(channel=="EleTau" && SignalMass==1000) TauCleaning = 1 + 0.170;
  if(channel=="MuoTau" && SignalMass==1500) TauCleaning = 1 + 0.006;
  if(channel=="EleTau" && SignalMass==1500) TauCleaning = 1 + 0.005;
  if(channel=="MuoTau" && SignalMass==2000) TauCleaning = 1 + 0.007;
  if(channel=="EleTau" && SignalMass==2000) TauCleaning = 1 + 0.008;
  if(channel=="MuoTau" && SignalMass==2500) TauCleaning = 1 + 0.008;
  if(channel=="EleTau" && SignalMass==2500) TauCleaning = 1 + 0.018;

  //DATACARD
  if(BkgYield==0) BkgYield=0.00001;
  ofstream myfile;
  myfile.open(saveName); 
  myfile<<"imax 1"<<endl;
  myfile<<"jmax 1"<<endl;
  myfile<<"kmax *"<<endl;
  myfile<<"---------------------------------------------"<<endl;
  myfile<<"bin                         "<<channel<<endl;
  myfile<<"observation                 -1"<<endl;
  myfile<<"---------------------------------------------"<<endl;
  myfile<<"bin                         "<<channel<<"           "<<channel<<endl;
  myfile<<"process                     signal"<<channel<<"     background"<<channel<<endl;
  myfile<<"process                     0                1"<<endl;
  myfile<<"rate                        "<<SigYield<<"         "<<BkgYield<<endl;
  myfile<<"-----------------------------------------------"<<endl;
  myfile<<"lumi              lnN       1.026            -"<<endl;
  if(NSideband!=0)                           myfile<<"Bkg"<<channel<<"         gmN "<<NSideband<<"    -                "<<BkgYield/NSideband<<endl;
  else                                       myfile<<"Bkg"<<channel<<"         gmN "<<NSideband<<"    -                0"<<endl;
  if(channel=="MuoTau" || channel=="EleTau") myfile<<"TauIDStandard     lnN       "<<TauIDStandard<<"            -"<<endl;
  if(channel=="MuoTau" || channel=="EleTau") myfile<<"TauCleaning"<<channel<<" lnN       "<<TauCleaning<<"             -"<<endl;
  //myfile<<"alphaSyst lnN   -                "<<1+AlphaErr<<endl;
}
