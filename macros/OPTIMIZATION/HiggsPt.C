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
  float lep1PtCut=10; float tauPtCut=35; int bTagCut=1; float deltaRCutFL=1; float deltaRCutSL=1; float MassSvfitCut=10000; float METCutFL = 100; float METCutSL = 50; 
  float MassVisCutSL=0;  
  float MassVisCutFL1=0; //ELE-MUO
  float MassVisCutFL2=0; //ELE-ELE e MUO-MUO
  bool beginning=true;

  vector<string> PLOT;              vector<int> BIN;   vector<float> MIN;  vector<float> MAX;   vector<float> MAXY;   vector<TString> AXIS;
  PLOT.push_back("PtSvfit");        BIN.push_back(50); MIN.push_back(0);   MAX.push_back(1000); MAXY.push_back(90);   AXIS.push_back("pt(#tau-#tau) [GeV] (SVFit)");

  TH1F *Optimization1  = new TH1F("Optimization1", "Optimization1", 60,0,600);
  TH1F *Optimization2  = new TH1F("Optimization2", "Optimization2", 60,0,600);
  TH1F *Optimization3  = new TH1F("Optimization3", "Optimization3", 60,0,600);
  TH1F *Optimization4  = new TH1F("Optimization4", "Optimization4", 60,0,600);
  TH1F *SignalEfficiency1 = new TH1F("SignalEfficiency1","SignalEfficiency1",60,0,600);
  TH1F *SignalEfficiency2 = new TH1F("SignalEfficiency2","SignalEfficiency2",60,0,600);
  TH1F *BackgroundYield1  = new TH1F("BackgroundYield1", "BackgroundYield1", 60,0,600);
  TH1F *BackgroundYield2  = new TH1F("BackgroundYield2", "BackgroundYield2", 60,0,600);

  float max1 = 0;
  float max2 = 0;

  for(int type=0; type<2; type++){
    bool SL=true;
    if(type==0) SL=true;
    if(type==1) SL=false;
    for(int k=0; k<2; k++){
      float ZH1000Integral = 0;
      float ZH2500Integral = 0;
      int KK=1;
      float PtSvfit = -10;
      for(int i=0; i<Optimization1->GetXaxis()->GetNbins(); i++){
	PtSvfit = PtSvfit + 10;
	char *plot = PLOT[0].c_str();
	TString name = PLOT[0];
	int bin=BIN[0]; 
	float min=MIN[0]; 
	float max=MAX[0]; 
	float maxy=MAXY[0];
	TString axis = AXIS[0];

	if(!SL && k==0) {maxy=2; }
	if(!SL && k==1) {maxy=3; }
	if(SL  && k==0) {maxy=5; }
	if(SL  && k==1) {maxy=3; }

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
	  if(k==0){
	    if(j==0){ 
	      sprintf(openTree, "demo/TreeEleMuo"); 
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>837  && XMassSVFit<1163&& PtSvfit>%f)",
		      CUTPre,METCutFL,MassVisCutFL1,lep1PtCut,deltaRCutFL,PtSvfit);
	    }else if(j==1){ 
	      sprintf(openTree, "demo/TreeMuoMuo"); 
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>837  && XMassSVFit<1163&& PtSvfit>%f)",
		      CUTPre,METCutFL,MassVisCutFL2,lep1PtCut,deltaRCutFL,PtSvfit);
	    }else if(j==2){ 
	      sprintf(openTree, "demo/TreeEleEle"); 
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>837  && XMassSVFit<1163&& PtSvfit>%f)",
		      CUTPre,METCutFL,MassVisCutFL2,lep1PtCut,deltaRCutFL,PtSvfit);
	    }else if(j==3){ 
	      sprintf(openTree, "demo/TreeMuoTau");
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>837  && XMassSVFit<1163&& PtSvfit>%f)",
		      CUTPre,METCutSL,MassVisCutSL, tauPtCut, deltaRCutSL,PtSvfit);
	    }else if(j==4){ 
	      sprintf(openTree, "demo/TreeEleTau"); 
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>837  && XMassSVFit<1163&& PtSvfit>%f)",
		      CUTPre,METCutSL,MassVisCutSL, tauPtCut, deltaRCutSL,PtSvfit);
	    }
	  }else{
	    if(j==0){ 
	      sprintf(openTree, "demo/TreeEleMuo"); 
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>2125 && XMassSVFit<2875&& PtSvfit>%f)",
		      CUTPre,METCutFL,MassVisCutFL1,lep1PtCut,deltaRCutFL,PtSvfit);
	    }else if(j==1){ 
	      sprintf(openTree, "demo/TreeMuoMuo"); 
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>2125 && XMassSVFit<2875&& PtSvfit>%f)",
		      CUTPre,METCutFL,MassVisCutFL2,lep1PtCut,deltaRCutFL,PtSvfit);
	    }else if(j==2){ 
	      sprintf(openTree, "demo/TreeEleEle"); 
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>2125 && XMassSVFit<2875&& PtSvfit>%f)",
		      CUTPre,METCutFL,MassVisCutFL2,lep1PtCut,deltaRCutSL,PtSvfit);
	    }else if(j==3){ 
	      sprintf(openTree, "demo/TreeMuoTau");
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>2125 && XMassSVFit<2875&& PtSvfit>%f)",
		      CUTPre,METCutSL,MassVisCutSL, tauPtCut, deltaRCutSL,PtSvfit);
	    }else if(j==4){ 
	      sprintf(openTree, "demo/TreeEleTau"); 
	      sprintf(CUT, "%s && met>%f && MassVis>%f && lep1Pt>%f && dRLep1Lep2<%f && XMassSVFit>2125 && XMassSVFit<2875&& PtSvfit>%f)",
		      CUTPre,METCutSL,MassVisCutSL, tauPtCut, deltaRCutSL,PtSvfit);
	    }
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
	float Bkg = (w_DY100*DY100->Integral() + w_DY70*DY70->Integral() + w_DYM50_100*DYM50_100->Integral() + 
		     w_DYM50_70*DYM50_70->Integral() + 1.9*w_QCD250*QCD250->Integral() + 1.9*w_QCD500*QCD500->Integral() + 
		     1.9*w_QCD1000*QCD1000->Integral() + w_TT*TT->Integral() + w_WW*WW->Integral() + w_WZ*WZ->Integral() + 
		     w_ZZ*ZZ->Integral() + w_WJets180*WJets180->Integral()
		     );
	float BkgErr = sqrt(w_DY100*w_DY100*DY100->Integral() + w_DY70*w_DY70*DY70->Integral() + w_DYM50_100*w_DYM50_100*DYM50_100->Integral() + 
			    w_DYM50_70*w_DYM50_70*DYM50_70->Integral() + 1.9*1.9*w_QCD250*w_QCD250*QCD250->Integral() + 1.9*1.9*w_QCD500*w_QCD500*QCD500->Integral() + 
			    1.9*1.9*w_QCD1000*w_QCD1000*QCD1000->Integral() + w_TT*w_TT*TT->Integral() + w_WW*w_WW*WW->Integral() + w_WZ*w_WZ*WZ->Integral() + 
			    w_ZZ*w_ZZ*ZZ->Integral() + w_WJets180*w_WJets180*WJets180->Integral()
			    );
	cout<<i<<" "<<PtSvfit<<"                        "<<Bkg<<" +/- "<<BkgErr<<endl;

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
	DY100->Add(WW);
	DY100->Add(TT);
	DY100->Add(WJets180);
	DY100->Add(QCD1000);

	if(i==0 && k==0) max1=DY100->Integral();
	if(i==0 && k==1) max2=DY100->Integral();

	if(i==0) ZH2500Integral = ZH2500->Integral(); //SEMI LEP - 2500
	if(i==0) ZH1000Integral = ZH1000->Integral(); //SEMI LEP - 1000
	if(k==0){
	  Optimization1->SetBinContent(KK, (ZH1000->Integral()/ZH1000Integral)/(1+sqrt(DY100->Integral())));
	  Optimization3->SetBinContent(KK, (ZH1000->Integral()/ZH1000Integral)/(1+sqrt(DY100->Integral() + BkgErr*BkgErr/(Bkg*Bkg))));
	  SignalEfficiency1->SetBinContent(KK, ZH1000->Integral()/ZH1000Integral);
	  BackgroundYield1->SetBinContent(KK, DY100->Integral());
	  Optimization1->SetBinError(KK, 0.00001);
	  Optimization3->SetBinError(KK, 0.00001);
	  SignalEfficiency1->SetBinError(KK, 0.00001);
	  BackgroundYield1->SetBinError(KK, 0.00001);
	} else {
	  Optimization2->SetBinContent(KK, (ZH2500->Integral()/ZH2500Integral)/(1+sqrt(DY100->Integral())));
	  Optimization4->SetBinContent(KK, (ZH2500->Integral()/ZH2500Integral)/(1+sqrt(DY100->Integral() + BkgErr*BkgErr/(Bkg*Bkg))));
	  SignalEfficiency2->SetBinContent(KK, ZH2500->Integral()/ZH2500Integral);
	  BackgroundYield2->SetBinContent(KK, DY100->Integral());
	  Optimization2->SetBinError(KK, 0.00001);
	  Optimization4->SetBinError(KK, 0.00001);
	  SignalEfficiency2->SetBinError(KK, 0.00001);
	  BackgroundYield2->SetBinError(KK, 0.00001);
	}
	KK=KK+1;
      }
    }
    
    char Lep1pt      [2];  sprintf(Lep1pt,      "%.0f", lep1PtCut);    TString lep1pt       = Lep1pt      ;
    char Taupt       [2];  sprintf(Taupt,       "%.0f", tauPtCut);	   TString taupt        = Taupt       ;
    char Btag        [2];  sprintf(Btag,        "%.0f", bTagCut);	   TString btag         = Btag        ;
    char Drcut       [2];  sprintf(Drcut,       "%.0f", deltaRCutSL);	   TString drcut        = Drcut       ;
    char Masssvfitcut[2];  sprintf(Masssvfitcut,"%.0f", MassSvfitCut);     TString masssvfitcut = Masssvfitcut;
    char Massviscut  [2];  sprintf(Massviscut,  "%.0f", MassVisCutSL);     TString massviscut   = Massviscut  ;
    char Met         [2];  sprintf(Met,         "%.0f", METCutSL);	   TString met          = Met         ;
    
    TCanvas* c2 = new TCanvas("c2","c2",0,0,800,600);
    Optimization1->Draw("E");
    Optimization1->SetMinimum(0.0);
    Optimization1->SetMaximum(0.4);
    Optimization1->SetMarkerStyle(21);
    Optimization1->SetLineColor(1);
    Optimization1->GetYaxis()->SetTitleSize(0.045);
    Optimization1->GetXaxis()->SetTitleSize(0.045);
    Optimization1->GetYaxis()->SetLabelSize(0.045);
    Optimization1->GetXaxis()->SetLabelSize(0.045);
    Optimization1->SetTitle("");
    Optimization1->GetXaxis()->SetTitle("pt(#tau-#tau) [GeV] (SVFit)");
    Optimization1->GetYaxis()->SetTitle("#epsilon_{S}/1+#sqrt{B}");
    if(save && !SL) 
      c2->SaveAs("PtSvfit_1000_OPTIMIZATION_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_FL.pdf");
    if(save && SL)  
      c2->SaveAs("PtSvfit_1000_OPTIMIZATION_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_SL.pdf");
    
    TCanvas* c3 = new TCanvas("c3","c3",0,0,800,600);
    Optimization2->Draw("E");
    Optimization2->SetMinimum(0.0);
    Optimization2->SetMaximum(0.5);
    Optimization2->SetMarkerStyle(21);
    Optimization2->SetLineColor(1);
    Optimization2->GetYaxis()->SetTitleSize(0.045);
    Optimization2->GetXaxis()->SetTitleSize(0.045);
    Optimization2->GetYaxis()->SetLabelSize(0.045);
    Optimization2->GetXaxis()->SetLabelSize(0.045);
    Optimization2->SetTitle("");
    Optimization2->GetXaxis()->SetTitle("pt(#tau-#tau) [GeV] (SVFit)");
    Optimization2->GetYaxis()->SetTitle("#epsilon_{S}/1+#sqrt{B}");
    if(save && !SL) 
      c3->SaveAs("PtSvfit_2500_OPTIMIZATION_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_FL.pdf");
    if(save && SL)  
      c3->SaveAs("PtSvfit_2500_OPTIMIZATION_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_SL.pdf");

    
    TCanvas* c4 = new TCanvas("c4","c4",0,0,800,600);
    TPad *pad1 = new TPad("pad1","",0,0,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000); //will be transparent
    pad1->Draw();
    pad1->cd();

    SignalEfficiency1->Draw("E");
    SignalEfficiency1->SetMinimum(0.0);
    SignalEfficiency1->SetMaximum(1.4);
    SignalEfficiency1->SetMarkerStyle(21);
    SignalEfficiency1->SetLineColor(1);
    SignalEfficiency1->GetYaxis()->SetTitleSize(0.045);
    SignalEfficiency1->GetXaxis()->SetTitleSize(0.045);
    SignalEfficiency1->GetYaxis()->SetLabelSize(0.045);
    SignalEfficiency1->GetXaxis()->SetLabelSize(0.045);
    SignalEfficiency1->SetTitle("");
    SignalEfficiency1->GetXaxis()->SetTitle("pt(#tau-#tau) [GeV] (SVFit)");
    SignalEfficiency1->GetYaxis()->SetTitle("Signal Efficiency");

    c4->cd();
    Double_t ymin = 0;
    Double_t ymax = max1*1.5;
    Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
    Double_t xmin = SignalEfficiency1->GetXaxis()->GetXmin();
    Double_t xmax = SignalEfficiency1->GetXaxis()->GetXmax();
    Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
    pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad2->Draw();
    pad2->cd();

    BackgroundYield1->SetMinimum(0.0);
    BackgroundYield1->GetYaxis()->SetTitleSize(0.045);
    BackgroundYield1->GetXaxis()->SetTitleSize(0.045);
    BackgroundYield1->GetYaxis()->SetLabelSize(0.045);
    BackgroundYield1->GetXaxis()->SetLabelSize(0.045);
    BackgroundYield1->SetTitle("");
    BackgroundYield1->GetXaxis()->SetTitle("pt(#tau-#tau) [GeV] (SVFit)");
    BackgroundYield1->GetYaxis()->SetTitle("Background yield");
    BackgroundYield1->SetMarkerStyle(21);
    BackgroundYield1->SetMarkerColor(2);
    BackgroundYield1->SetLineColor(2);
    BackgroundYield1->Draw("][sames");
   
    // draw axis on the right side of the pad
    TGaxis *Axis1 = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
    Axis1->SetLabelColor(kRed);
    Axis1->SetTitle("Background Yield");
    Axis1->SetTitleColor(kRed);
    Axis1->SetTitleSize(0.045);
    Axis1->SetLabelSize(0.045);
    Axis1->Draw();
  
    TLegend *pl = new TLegend(0.57,0.75,0.89,0.89);
    pl->SetTextSize(0.03); 
    pl->SetFillColor(0);
    TLegendEntry *ple = pl->AddEntry(SignalEfficiency1, "Signal Efficiency",  "LP");
    ple               = pl->AddEntry(BackgroundYield1,  "Background Yield",  "LP");
    pl->Draw();

    if(save && !SL) 
      c4->SaveAs("PtSvfit_1000_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_FL.pdf");
    if(save && SL)  
      c4->SaveAs("PtSvfit_1000_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_SL.pdf");

    
    TCanvas* c5 = new TCanvas("c5","c5",0,0,800,600);
    TPad *pad3 = new TPad("pad3","",0,0,1,1);
    TPad *pad4 = new TPad("pad4","",0,0,1,1);
    pad4->SetFillStyle(4000); //will be transparent
    pad3->Draw();
    pad3->cd();

    SignalEfficiency2->Draw("E");
    SignalEfficiency2->SetMinimum(0.0);
    SignalEfficiency2->SetMaximum(1.4);
    SignalEfficiency2->SetMarkerStyle(21);
    SignalEfficiency2->SetLineColor(1);
    SignalEfficiency2->GetYaxis()->SetTitleSize(0.045);
    SignalEfficiency2->GetXaxis()->SetTitleSize(0.045);
    SignalEfficiency2->GetYaxis()->SetLabelSize(0.045);
    SignalEfficiency2->GetXaxis()->SetLabelSize(0.045);
    SignalEfficiency2->SetTitle("");
    SignalEfficiency2->GetXaxis()->SetTitle("pt(#tau-#tau) [GeV] (SVFit)");
    SignalEfficiency2->GetYaxis()->SetTitle("Signal Efficiency");

    c5->cd();
    Double_t ymin = 0;
    Double_t ymax = max2*1.5;
    Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
    Double_t xmin = SignalEfficiency1->GetXaxis()->GetXmin();
    Double_t xmax = SignalEfficiency1->GetXaxis()->GetXmax();
    Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
    pad4->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad4->Draw();
    pad4->cd();

    BackgroundYield2->SetMinimum(0.0);
    BackgroundYield2->GetYaxis()->SetTitleSize(0.045);
    BackgroundYield2->GetXaxis()->SetTitleSize(0.045);
    BackgroundYield2->GetYaxis()->SetLabelSize(0.045);
    BackgroundYield2->GetXaxis()->SetLabelSize(0.045);
    BackgroundYield2->SetTitle("");
    BackgroundYield2->GetXaxis()->SetTitle("pt(#tau-#tau) [GeV] (SVFit)");
    BackgroundYield2->GetYaxis()->SetTitle("Background yield");
    BackgroundYield2->SetMarkerStyle(21);
    BackgroundYield2->SetMarkerColor(2);
    BackgroundYield2->SetLineColor(2);
    BackgroundYield2->Draw("][sames");
   
    // draw axis on the right side of the pad
    TGaxis *Axis2 = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,510,"+L");
    Axis2->SetLabelColor(kRed);
    Axis2->SetTitle("Background Yield");
    Axis2->SetTitleColor(kRed);
    Axis2->SetTitleSize(0.045);
    Axis2->SetLabelSize(0.045);
    Axis2->Draw();
  
    TLegend *pl = new TLegend(0.57,0.75,0.89,0.89);
    pl->SetTextSize(0.03); 
    pl->SetFillColor(0);
    TLegendEntry *ple = pl->AddEntry(SignalEfficiency2, "Signal Efficiency",  "LP");
    ple               = pl->AddEntry(BackgroundYield2,  "Background Yield",  "LP");
    pl->Draw();

    if(save && !SL) 
      c5->SaveAs("PtSvfit_2500_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_FL.pdf");
    if(save && SL)  
      c5->SaveAs("PtSvfit_2500_MET"+met+"_lep1Pt"+lep1pt+"_tauPt"+taupt+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+"_SL.pdf");
  }
}
