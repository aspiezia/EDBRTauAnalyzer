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


  char channel  [500]; sprintf(channel,  "EleTau"); 
  char demo     [500]; sprintf(demo,     "demo/Tree"); 
  char openTree [500]; sprintf(openTree, "%s%s",demo,channel); 
  TString CHANNEL = channel;
  bool save=true; 
  bool log=false;
  bool lepTauIso = true;

  vector<string> PLOT;                vector<int> BIN;   vector<float> MIN;  vector<float> MAX;   vector<float> MAXY;   vector<TString> AXIS;
  if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo" || CHANNEL=="EleMuo"){
    /*PLOT.push_back("MassCA");         BIN.push_back(30); MIN.push_back(0);   MAX.push_back(300);  MAXY.push_back(4);  AXIS.push_back("M(lep1,lep2) [GeV] (Collinear Approximation)");    	
    PLOT.push_back("MassEff");        BIN.push_back(35); MIN.push_back(0);   MAX.push_back(350);  MAXY.push_back(4);  AXIS.push_back("M(lep1,lep2) [GeV] (effective)");
    PLOT.push_back("MassSvfit");      BIN.push_back(25); MIN.push_back(0);   MAX.push_back(250);  MAXY.push_back(6);  AXIS.push_back("M(lep1,lep2) [GeV] (SvFit)");
    PLOT.push_back("MassVis");        BIN.push_back(25); MIN.push_back(0);   MAX.push_back(250);  MAXY.push_back(8);  AXIS.push_back("M(lep1,lep2) [GeV] (visible)");
    PLOT.push_back("NVertices");      BIN.push_back(41); MIN.push_back(-0.5);MAX.push_back(40.5); MAXY.push_back(2);  AXIS.push_back("N(vertices)");  
    PLOT.push_back("XMassCA");        BIN.push_back(30); MIN.push_back(0);   MAX.push_back(3000); MAXY.push_back(5);  AXIS.push_back("M(Z,H) [GeV] (Collinear Approximation)");
    PLOT.push_back("XMassEff");       BIN.push_back(30); MIN.push_back(0);   MAX.push_back(3000); MAXY.push_back(5);  AXIS.push_back("M(Z,H) [GeV] (effective)");*/
    PLOT.push_back("XMassSVFit");     BIN.push_back(30); MIN.push_back(0);   MAX.push_back(3000); MAXY.push_back(15); AXIS.push_back("M(Z,H) [GeV] (SvFit)");/*
    PLOT.push_back("XMassVis");       BIN.push_back(30); MIN.push_back(0);   MAX.push_back(3000); MAXY.push_back(4);  AXIS.push_back("M(Z,H) [GeV] (visible)");
    PLOT.push_back("dPhiJetMet");     BIN.push_back(20); MIN.push_back(-4);  MAX.push_back(4);    MAXY.push_back(8);  AXIS.push_back("#Delta#phi(jet,met)");  
    PLOT.push_back("dPhiLepMet1");    BIN.push_back(20); MIN.push_back(-4);  MAX.push_back(4);    MAXY.push_back(6);  AXIS.push_back("#Delta#phi(lep1,met)");
    PLOT.push_back("dPhiLepMet2");    BIN.push_back(20); MIN.push_back(-4);  MAX.push_back(4);    MAXY.push_back(6);  AXIS.push_back("#Delta#phi(lep2,met)");
    PLOT.push_back("dRJetLep1");      BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(6);  AXIS.push_back("#DeltaR(jet,lep1)");
    PLOT.push_back("dRJetLep2");      BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(6);  AXIS.push_back("#DeltaR(jet,lep2)");
    PLOT.push_back("dRJetMet");       BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(8);  AXIS.push_back("#DeltaR(jet,met)");
    PLOT.push_back("dRLepLep");       BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(7);  AXIS.push_back("#DeltaR(lep1,lep2)");
    PLOT.push_back("dRLepMet1");      BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(4);  AXIS.push_back("#DeltaR(lep1,met)");
    PLOT.push_back("dRLepMet2");      BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(4);  AXIS.push_back("#DeltaR(lep2,met)");
    PLOT.push_back("dRZZCA");         BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(6);  AXIS.push_back("#DeltaR(Z,H) (Collinear Approximation)");
    PLOT.push_back("dRZZEff");        BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(8);  AXIS.push_back("#DeltaR(Z,H) (effective)");
    PLOT.push_back("dRZZSvFit");      BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(8);  AXIS.push_back("#DeltaR(Z,H) (SvFit)");
    PLOT.push_back("dRZZVis");        BIN.push_back(30); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(8);  AXIS.push_back("#DeltaR(Z,H) (visible)");
    PLOT.push_back("jetEta");         BIN.push_back(32); MIN.push_back(-4);  MAX.push_back(4);    MAXY.push_back(3);  AXIS.push_back("jet #eta");
    PLOT.push_back("jetMass");        BIN.push_back(20); MIN.push_back(70);  MAX.push_back(110);  MAXY.push_back(2);  AXIS.push_back("jet mass [GeV]");
    PLOT.push_back("jetPt");          BIN.push_back(30); MIN.push_back(400); MAX.push_back(1400); MAXY.push_back(6);  AXIS.push_back("jet pt [GeV]");
    PLOT.push_back("jetSubjettiness");BIN.push_back(20); MIN.push_back(0);   MAX.push_back(1);    MAXY.push_back(4);  AXIS.push_back("jet #tau_{21}");
    PLOT.push_back("lepEta1");        BIN.push_back(20); MIN.push_back(-4);  MAX.push_back(4);    MAXY.push_back(4);  AXIS.push_back("lep1 #eta");
    PLOT.push_back("lepEta2");        BIN.push_back(20); MIN.push_back(-4);  MAX.push_back(4);    MAXY.push_back(4);  AXIS.push_back("lep2 #eta");
    PLOT.push_back("lepPt1");         BIN.push_back(35); MIN.push_back(0);   MAX.push_back(700);  MAXY.push_back(3);  AXIS.push_back("lep1 pt [GeV]");
    PLOT.push_back("lepPt2");         BIN.push_back(25); MIN.push_back(0);   MAX.push_back(500);  MAXY.push_back(5);  AXIS.push_back("lep2 pt [GeV]");
    PLOT.push_back("met");            BIN.push_back(50); MIN.push_back(0);   MAX.push_back(1000); MAXY.push_back(10); AXIS.push_back("met [GeV]");
    PLOT.push_back("metPhi");         BIN.push_back(20); MIN.push_back(-4);  MAX.push_back(4);    MAXY.push_back(2);  AXIS.push_back("met #phi");*/
  } 
  else {/*
    PLOT.push_back("MassCA");          BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(300);  MAXY.push_back(4);  AXIS.push_back("M(tau,lep) [GeV] (Collinear Approximation)");    	
    PLOT.push_back("MassEff");         BIN.push_back(35);  MIN.push_back(0);    MAX.push_back(350);  MAXY.push_back(4);  AXIS.push_back("M(tau,lep) [GeV] (effective)");
    PLOT.push_back("MassSvfit");       BIN.push_back(25);  MIN.push_back(0);    MAX.push_back(250);  MAXY.push_back(6);  AXIS.push_back("M(tau,lep) [GeV] (SvFit)");
    PLOT.push_back("MassVis");         BIN.push_back(25);  MIN.push_back(0);    MAX.push_back(250);  MAXY.push_back(8);  AXIS.push_back("M(tau,lep) [GeV] (visible)");
    PLOT.push_back("NVertices");       BIN.push_back(41);  MIN.push_back(-0.5); MAX.push_back(40.5); MAXY.push_back(2);  AXIS.push_back("N(vertices)");  
    PLOT.push_back("XMassCA");         BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(3000); MAXY.push_back(5);  AXIS.push_back("M(Z,H) [GeV] (Collinear Approximation)");  
    PLOT.push_back("XMassEff");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(3000); MAXY.push_back(5);  AXIS.push_back("M(Z,H) [GeV] (effective)");*/
    PLOT.push_back("XMassSVFit");      BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(3000); MAXY.push_back(15); AXIS.push_back("M(Z,H) [GeV] (SvFit)");/*
    PLOT.push_back("XMassVis");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(3000); MAXY.push_back(4);  AXIS.push_back("M(Z,H) [GeV] (visible)");
    PLOT.push_back("dPhiJetMet");      BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(8);  AXIS.push_back("#Delta#phi(jet,met)");  
    PLOT.push_back("dPhiTauMet");      BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(6);  AXIS.push_back("#Delta#phi(tau,met)");
    PLOT.push_back("dPhiLepMet");      BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(6);  AXIS.push_back("#Delta#phi(lep,met)");
    PLOT.push_back("dRJetTau");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(6);  AXIS.push_back("#DeltaR(jet,tau)");
    PLOT.push_back("dRJetLep");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(6);  AXIS.push_back("#DeltaR(jet,lep)");
    PLOT.push_back("dRJetMet");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(8);  AXIS.push_back("#DeltaR(jet,met)");
    PLOT.push_back("dRTauLep");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(7);  AXIS.push_back("#DeltaR(tau,lep)");
    PLOT.push_back("dRTauMet");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(4);  AXIS.push_back("#DeltaR(tau,met)");
    PLOT.push_back("dRLepMet");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(4);  AXIS.push_back("#DeltaR(lep,met)");
    PLOT.push_back("dRZZCA");          BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(6);  AXIS.push_back("#DeltaR(Z,H) (Collinear Approximation)");
    PLOT.push_back("dRZZEff");         BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(8);  AXIS.push_back("#DeltaR(Z,H) (effective)");
    PLOT.push_back("dRZZSvFit");       BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(8);  AXIS.push_back("#DeltaR(Z,H) (SvFit)");
    PLOT.push_back("dRZZVis");         BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(8);  AXIS.push_back("#DeltaR(Z,H) (visible)");
    PLOT.push_back("jetEta");          BIN.push_back(32);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(3);  AXIS.push_back("jet #eta");
    PLOT.push_back("jetMass");         BIN.push_back(20);  MIN.push_back(70);   MAX.push_back(110);  MAXY.push_back(2);  AXIS.push_back("jet mass [GeV]");
    PLOT.push_back("jetPt");           BIN.push_back(30);  MIN.push_back(400);  MAX.push_back(1400); MAXY.push_back(6);  AXIS.push_back("jet pt [GeV]");
    PLOT.push_back("jetSubjettiness"); BIN.push_back(20);  MIN.push_back(0);    MAX.push_back(1);    MAXY.push_back(4);  AXIS.push_back("jet #tau_{21}");
    PLOT.push_back("tauEta");          BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(4);  AXIS.push_back("tau #eta");
    PLOT.push_back("lepEta");          BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(4);  AXIS.push_back("lep #eta");
    PLOT.push_back("tauPt");           BIN.push_back(35);  MIN.push_back(0);    MAX.push_back(700);  MAXY.push_back(3);  AXIS.push_back("tau pt [GeV]");
    PLOT.push_back("lepPt");           BIN.push_back(25);  MIN.push_back(0);    MAX.push_back(500);  MAXY.push_back(5);  AXIS.push_back("lep pt [GeV]");
    PLOT.push_back("met");             BIN.push_back(25);  MIN.push_back(50);   MAX.push_back(1050); MAXY.push_back(10); AXIS.push_back("met [GeV]");
    PLOT.push_back("metPhi");          BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(2);  AXIS.push_back("met #phi");*/
  }

  //PLOT.push_back("metPx");           BIN.push_back(25);  MIN.push_back(-100); MAX.push_back(100);  MAXY.push_back(2);  AXIS.push_back("met_x [GeV]");
  //PLOT.push_back("metPy");           BIN.push_back(25);  MIN.push_back(-100); MAX.push_back(100);  MAXY.push_back(2);  AXIS.push_back("met_y [GeV]");
  //PLOT.push_back("nbtagsL");         BIN.push_back(7);   MIN.push_back(-0.5); MAX.push_back(6.5);  MAXY.push_back(10); AXIS.push_back("N(loose btag)");
  //PLOT.push_back("nbtagsM");         BIN.push_back(7);   MIN.push_back(-0.5); MAX.push_back(6.5);  MAXY.push_back(15); AXIS.push_back("N(medium btag)");
  //PLOT.push_back("nbtagsT");         BIN.push_back(7);   MIN.push_back(-0.5); MAX.push_back(6.5);  MAXY.push_back(15); AXIS.push_back("N(tight btag)");
  //PLOT.push_back("uncorrmet");       BIN.push_back(25);  MIN.push_back(30);   MAX.push_back(500);  MAXY.push_back(10);   AXIS.push_back("uncorrected met [GeV]");
  //PLOT.push_back("uncorrmetPhi");    BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(10);   AXIS.push_back("uncorrected met #phi");
  //PLOT.push_back("lepNumberEle");    BIN.push_back(4);   MIN.push_back(-0.5); MAX.push_back(3.5);  MAXY.push_back(30); AXIS.push_back("Number of Electron");
  //PLOT.push_back("lepNumberMuo");    BIN.push_back(4);   MIN.push_back(-0.5); MAX.push_back(3.5);  MAXY.push_back(30); AXIS.push_back("Number of Muon");
  //PLOT.push_back("lepPFIso1");       BIN.push_back(15);  MIN.push_back(0);    MAX.push_back(0.15); MAXY.push_back(10); AXIS.push_back("lep1 PFIso");
  //PLOT.push_back("lepPFIso2");       BIN.push_back(15);  MIN.push_back(0);    MAX.push_back(0.15); MAXY.push_back(10); AXIS.push_back("lep2 PFIso");
  //PLOT.push_back("lepDetIso1");      BIN.push_back(15);  MIN.push_back(0);    MAX.push_back(0.15); MAXY.push_back(8);  AXIS.push_back("lep1 DetIso");
  //PLOT.push_back("lepDetIso2");      BIN.push_back(15);  MIN.push_back(0);    MAX.push_back(0.15); MAXY.push_back(8);  AXIS.push_back("lep2 DetIso");
  
  for(int i=0; i<PLOT.size(); i++){
    //for(int i=7; i<8; i++){
    char *plot = PLOT[i].c_str();
    TString name = PLOT[i];
    int bin=BIN[i]; 
    float min=MIN[i]; 
    float max=MAX[i]; 
    float maxy=MAXY[i];
    TString axis = AXIS[i];

    char CUT [500]; 
    if(!lepTauIso){
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo") sprintf(CUT, "trigger*PUWeight*(met>100 && uncorrmet>30 && MassVis>10)");
      else if(CHANNEL=="EleMuo") sprintf(CUT, "trigger*PUWeight*(met>100 && uncorrmet>30)");
      else sprintf(CUT, "trigger*PUWeight*(met>100 && uncorrmet>30)");
    } else {
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo") sprintf(CUT, "trigger*PUWeight*(met>100 && uncorrmet>30 && MassVis>10)");
      else if(CHANNEL=="EleMuo") sprintf(CUT, "trigger*PUWeight*(met>100 && uncorrmet>30)");
      else if(CHANNEL=="EleTau")  sprintf(CUT, "trigger*PUWeight*(met>100 && uncorrmet>30 && lepCorrPFIso<0.1)");
      else   sprintf(CUT, "trigger*PUWeight*(met>100 && uncorrmet>30 && lepCorrPFIso<0.2)");
    }

    cout<<endl;
    cout<<"Plotting plot "<<i+1<<"/"<<PLOT.size()<<", i.e. "<<name<<endl;

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
  
    TFile *fileAA = TFile::Open("../RISULTATI/LepIso_010413/ZH1000_SL.root");    TTree *TreeAA = (TTree*)fileAA->Get(openTree); 
    TFile *fileBB = TFile::Open("../RISULTATI/LepIso_010413/ZH1500_SL.root");    TTree *TreeBB = (TTree*)fileBB->Get(openTree); 
    TFile *fileCC = TFile::Open("../RISULTATI/LepIso_010413/ZH2000_SL.root");    TTree *TreeCC = (TTree*)fileCC->Get(openTree);
    TFile *fileDD = TFile::Open("../RISULTATI/LepIso_010413/ZH2500_SL.root");    TTree *TreeDD = (TTree*)fileDD->Get(openTree); 
    TFile *file02 = TFile::Open("../RISULTATI/LepIso_010413/DY100_SL.root");     TTree *Tree02 = (TTree*)file02->Get(openTree); 
    TFile *file03 = TFile::Open("../RISULTATI/LepIso_010413/DY70_SL.root");      TTree *Tree03 = (TTree*)file03->Get(openTree); 
    TFile *file04 = TFile::Open("../RISULTATI/LepIso_010413/DYM50_100_SL.root"); TTree *Tree04 = (TTree*)file04->Get(openTree); 
    TFile *file05 = TFile::Open("../RISULTATI/LepIso_010413/DYM50_70_SL.root");  TTree *Tree05 = (TTree*)file05->Get(openTree); 
    TFile *file06 = TFile::Open("../RISULTATI/LepIso_010413/QCD1000_SL.root");   TTree *Tree06 = (TTree*)file06->Get(openTree); 
    TFile *file07 = TFile::Open("../RISULTATI/LepIso_010413/QCD250_SL.root");    TTree *Tree07 = (TTree*)file07->Get(openTree); 
    TFile *file08 = TFile::Open("../RISULTATI/LepIso_010413/QCD500_SL.root");    TTree *Tree08 = (TTree*)file08->Get(openTree); 
    TFile *file09 = TFile::Open("../RISULTATI/LepIso_010413/TT_SL.root");        TTree *Tree09 = (TTree*)file09->Get(openTree); 
    TFile *file10 = TFile::Open("../RISULTATI/LepIso_010413/WJets180_SL.root");  TTree *Tree10 = (TTree*)file10->Get(openTree); 
    TFile *file11 = TFile::Open("../RISULTATI/LepIso_010413/WW_SL.root");        TTree *Tree11 = (TTree*)file11->Get(openTree); 
    TFile *file12 = TFile::Open("../RISULTATI/LepIso_010413/WZ_SL.root");        TTree *Tree12 = (TTree*)file12->Get(openTree); 
    TFile *file13 = TFile::Open("../RISULTATI/LepIso_010413/ZZ_SL.root");        TTree *Tree13 = (TTree*)file13->Get(openTree); 

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

    TreeAA->Draw(inputAA,CUT); if(TreeAA->Draw(inputAA,CUT)) {ZH1000    = hAA;}
    TreeBB->Draw(inputBB,CUT); if(TreeBB->Draw(inputBB,CUT)) {ZH1500    = hBB;}
    TreeCC->Draw(inputCC,CUT); if(TreeCC->Draw(inputCC,CUT)) {ZH2000    = hCC;}
    TreeDD->Draw(inputDD,CUT); if(TreeDD->Draw(inputDD,CUT)) {ZH2500    = hDD;}
    Tree02->Draw(input02,CUT); if(Tree02->Draw(input02,CUT)) {DY100     = h02;}
    Tree03->Draw(input03,CUT); if(Tree03->Draw(input03,CUT)) {DY70      = h03;}
    Tree04->Draw(input04,CUT); if(Tree04->Draw(input04,CUT)) {DYM50_100 = h04;}
    Tree05->Draw(input05,CUT); if(Tree05->Draw(input05,CUT)) {DYM50_70  = h05;}
    Tree06->Draw(input06,CUT); if(Tree06->Draw(input06,CUT)) {QCD1000   = h06;}
    Tree07->Draw(input07,CUT); if(Tree07->Draw(input07,CUT)) {QCD250    = h07;}
    Tree08->Draw(input08,CUT); if(Tree08->Draw(input08,CUT)) {QCD500    = h08;}
    Tree09->Draw(input09,CUT); if(Tree09->Draw(input09,CUT)) {TT        = h09;}
    Tree10->Draw(input10,CUT); if(Tree10->Draw(input10,CUT)) {WJets180  = h10;}
    Tree11->Draw(input11,CUT); if(Tree11->Draw(input11,CUT)) {WW        = h11;}
    Tree12->Draw(input12,CUT); if(Tree12->Draw(input12,CUT)) {WZ        = h12;}
    Tree13->Draw(input13,CUT); if(Tree13->Draw(input13,CUT)) {ZZ        = h13;}
  

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
    ZH1000->Scale(0.03);
    ZH1500->Scale(0.03);
    ZH2000->Scale(0.03);
    ZH2500->Scale(0.03);

    DY100->SetFillColor(kBlue-9);
    QCD1000->SetFillColor(kRed-2);
    WW->SetFillColor(kGreen+1);
    TT->SetFillColor(kBlue);
    WJets180->SetFillColor(kMagenta-5);
    ZH1000->SetLineColor(1);
    ZH1500->SetLineColor(2);
    ZH2000->SetLineColor(4);
    ZH2500->SetLineColor(8);
    ZH1000->SetLineWidth(2);
    ZH1500->SetLineWidth(2);
    ZH2000->SetLineWidth(2);
    ZH2500->SetLineWidth(2);
    ZH1000->SetFillColor(1);
    ZH1500->SetFillColor(2);
    ZH2000->SetFillColor(4);
    ZH2500->SetFillColor(8);
    ZH1000->SetFillStyle(3003);
    ZH1500->SetFillStyle(3003);
    ZH2000->SetFillStyle(3003);
    ZH2500->SetFillStyle(3003);
	
    //plots 
    TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);

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
    //hs->SetMaximum(maxy);

    TH1F *ERR = new TH1F("","",DY100->GetNbinsX(),DY100->GetXaxis()->GetXmin(),DY100->GetXaxis()->GetXmax());
    ERR->Sumw2();
    ERR->Add(DY100);
    ERR->Add(WW);
    ERR->Add(TT);
    ERR->Add(WJets180);
    ERR->Add(QCD1000);
    ERR->SetFillStyle(3005);
    ERR->SetFillColor(12);
    ERR->SetLineColor(12);
    //ERR->Draw("E2same");
    ZH1000->Draw("histo same");
    ZH1500->Draw("histo same");
    ZH2000->Draw("histo same");
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
    ple2 = pl2->AddEntry(ZH1000, "Signal (M_{X} = 1.0 TeV - arbitrary unit)",  "L");
    ple2 = pl2->AddEntry(ZH1500, "Signal (M_{X} = 1.5 TeV - arbitrary unit)",  "L");
    ple2 = pl2->AddEntry(ZH2000, "Signal (M_{X} = 2.0 TeV - arbitrary unit)",  "L");
    ple2 = pl2->AddEntry(ZH2500, "Signal (M_{X} = 2.5 TeV - arbitrary unit)",  "L");
    pl2->Draw();

    if(log) {
      hs->SetMinimum(0.05);
      c1->SetLogy();
    }
    if(save && !lepTauIso){
      c1->SaveAs(CHANNEL+"_"+name+"_SR.png");
      //c1->SaveAs(CHANNEL+"_"+name+"_SR.root");
      //c1->SaveAs(CHANNEL+"_"+name+"_SR.C");
    }
    if(save && lepTauIso){
      c1->SaveAs(CHANNEL+"_"+name+"_SR_tauIso.png");
      //c1->SaveAs(CHANNEL+"_"+name+"_SR_tauIso.root");
      //c1->SaveAs(CHANNEL+"_"+name+"_SR_tauIso.C");
    }
  }
}
