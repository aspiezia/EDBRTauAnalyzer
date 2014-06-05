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
  
  for(int j=0; j<5; j++){

    if(j==0){ char channel[500]; sprintf(channel,  "EleMuo"); }
    if(j==1){ char channel[500]; sprintf(channel,  "MuoMuo"); }
    if(j==2){ char channel[500]; sprintf(channel,  "EleEle"); }
    if(j==3){ char channel[500]; sprintf(channel,  "MuoTau"); }
    if(j==4){ char channel[500]; sprintf(channel,  "EleTau"); }

    char demo0     [500]; sprintf(demo0,     "demo/Tree"); 
    char demo1     [500]; sprintf(demo1,     "demo/TreeSB1"); 
    char demo2     [500]; sprintf(demo2,     "demo/TreeSB2"); 
    char openTree0 [500]; sprintf(openTree0, "%s%s",demo0,channel); 
    char openTree1 [500]; sprintf(openTree1, "%s%s",demo1,channel); 
    char openTree2 [500]; sprintf(openTree2, "%s%s",demo2,channel); 
    TString CHANNEL = channel;
    bool save=true; 
    bool log=true;

    vector<string> PLOT;              vector<int> BIN;   vector<float> MIN;  vector<float> MAX;   vector<float> MAXY; vector<TString> AXIS;
    PLOT.push_back("MassVis");        BIN.push_back(25); MIN.push_back(0);   MAX.push_back(250);  MAXY.push_back(500); AXIS.push_back("M(tau,lep) [GeV] (visible)");
    PLOT.push_back("NVertices");      BIN.push_back(41); MIN.push_back(-0.5);MAX.push_back(40.5); MAXY.push_back(500); AXIS.push_back("N(vertices)");
    PLOT.push_back("XMassSVFit");     BIN.push_back(60); MIN.push_back(0);   MAX.push_back(3000); MAXY.push_back(500); AXIS.push_back("M(Z,H) [GeV] (SvFit)");
    PLOT.push_back("jetMass");        BIN.push_back(20); MIN.push_back(20);  MAX.push_back(200);  MAXY.push_back(500); AXIS.push_back("jet mass [GeV]");
    PLOT.push_back("jetPt");          BIN.push_back(30); MIN.push_back(400); MAX.push_back(1000); MAXY.push_back(500); AXIS.push_back("jet pt [GeV]");
    PLOT.push_back("jetSubjettiness");BIN.push_back(20); MIN.push_back(0);   MAX.push_back(0.75); MAXY.push_back(500); AXIS.push_back("jet #tau_{21}");
    PLOT.push_back("met");            BIN.push_back(25); MIN.push_back(20);  MAX.push_back(520);  MAXY.push_back(500); AXIS.push_back("met [GeV]");
    PLOT.push_back("PtSvfit");        BIN.push_back(50); MIN.push_back(0);   MAX.push_back(700);  MAXY.push_back(500); AXIS.push_back("pt(#tau,#tau) [GeV] (SvFit)");
    PLOT.push_back("nbtagsL1");       BIN.push_back(6);  MIN.push_back(-0.5);MAX.push_back(5.5);  MAXY.push_back(500); AXIS.push_back("# loose btagged jet");
    PLOT.push_back("dRLep1Lep2");     BIN.push_back(60); MIN.push_back(0);   MAX.push_back(6);    MAXY.push_back(500); AXIS.push_back("#DeltaR (lep1,lep2)");
    PLOT.push_back("lep1Pt");         BIN.push_back(50); MIN.push_back(0);   MAX.push_back(500);  MAXY.push_back(500); AXIS.push_back("lep1 pt [GeV]");
    PLOT.push_back("lep2Pt");         BIN.push_back(50); MIN.push_back(0);   MAX.push_back(500);  MAXY.push_back(500); AXIS.push_back("lep2 pt [GeV]");

  

    for(int i=0; i<PLOT.size(); i++){
      //for(int i=3; i<4; i++){
      char *plot = PLOT[i].c_str();
      TString name = PLOT[i];
      int bin=BIN[i]; 
      float min=MIN[i]; 
      float max=MAX[i]; 
      float maxy=MAXY[i];
      TString axis = AXIS[i];
      if(PLOT[i]=="jetSubjettiness" && (CHANNEL=="MuoTau" || CHANNEL=="EleTau")) maxy=4000;

      //int N = 17; Double_t xbins[N] = {0,20,30,40,50,60,70,80,90,100,110,120,140,160,180,200,250}; //VIS MASS
      //int N = 15;  Double_t xbins[N] = {0,4,6,8,10,12,14,16,18,20,22,24,26,28,34}; //VERTICES
      //int N = 16;  Double_t xbins[N] = {0,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1800,2000,2500}; //X MASS

      char CUT [500];
      if(CHANNEL=="EleMuo")      sprintf(CUT, "trigger*PUWeight*(met>20)");
      else if(CHANNEL=="MuoMuo") sprintf(CUT, "trigger*PUWeight*(met>20)");
      else if(CHANNEL=="EleEle") sprintf(CUT, "trigger*PUWeight*(met>20)");
      else if(CHANNEL=="EleTau") sprintf(CUT, "trigger*PUWeight*(met>20)");
      else if(CHANNEL=="MuoTau") sprintf(CUT, "trigger*PUWeight*(met>20)");


      cout<<endl;
      cout<<"Plotting plot "<<i+1<<"/"<<PLOT.size()<<", i.e. "<<name<<endl;

      TH1F *data     = new TH1F("","",bin,min,max);
      TH1F *DY100    = new TH1F("","",bin,min,max);
      TH1F *DY70     = new TH1F("","",bin,min,max);
      TH1F *DYM50_100= new TH1F("","",bin,min,max);
      TH1F *DYM50_70 = new TH1F("","",bin,min,max);											     
      TH1F *QCD1000  = new TH1F("","",bin,min,max);
      TH1F *QCD250   = new TH1F("","",bin,min,max);
      TH1F *QCD500   = new TH1F("","",bin,min,max);
      TH1F *TT       = new TH1F("","",bin,min,max);
      TH1F *WJets180 = new TH1F("","",bin,min,max);
      TH1F *WW       = new TH1F("","",bin,min,max);
      TH1F *WZ       = new TH1F("","",bin,min,max);
      TH1F *ZZ       = new TH1F("","",bin,min,max);
	
      TFile *file01=TFile::Open("../../RISULTATI/analyzer_290514/data.root");      
      TFile *file02=TFile::Open("../../RISULTATI/analyzer_290514/DY100.root");      
      TFile *file03=TFile::Open("../../RISULTATI/analyzer_290514/DY70.root");       
      TFile *file04=TFile::Open("../../RISULTATI/analyzer_290514/DYM50_100.root");  
      TFile *file05=TFile::Open("../../RISULTATI/analyzer_290514/DYM50_70.root");   
      TFile *file06=TFile::Open("../../RISULTATI/analyzer_290514/QCD1000.root");    
      TFile *file07=TFile::Open("../../RISULTATI/analyzer_290514/QCD250.root");     
      TFile *file08=TFile::Open("../../RISULTATI/analyzer_290514/QCD500.root");     
      TFile *file09=TFile::Open("../../RISULTATI/analyzer_290514/TT.root");        
      TFile *file10=TFile::Open("../../RISULTATI/analyzer_290514/WJets180.root");   
      TFile *file11=TFile::Open("../../RISULTATI/analyzer_290514/WW.root");          
      TFile *file12=TFile::Open("../../RISULTATI/analyzer_290514/WZ.root");          
      TFile *file13=TFile::Open("../../RISULTATI/analyzer_290514/ZZ.root");        

      TTree *Tree01=(TTree*)file01->Get(openTree0);    TTree *Tree14=(TTree*)file01->Get(openTree1);    TTree *Tree27=(TTree*)file01->Get(openTree2);  
      TTree *Tree02=(TTree*)file02->Get(openTree0);    TTree *Tree15=(TTree*)file02->Get(openTree1);    TTree *Tree28=(TTree*)file02->Get(openTree2);  
      TTree *Tree03=(TTree*)file03->Get(openTree0);    TTree *Tree16=(TTree*)file03->Get(openTree1);    TTree *Tree29=(TTree*)file03->Get(openTree2);  
      TTree *Tree04=(TTree*)file04->Get(openTree0);    TTree *Tree17=(TTree*)file04->Get(openTree1);    TTree *Tree30=(TTree*)file04->Get(openTree2);  
      TTree *Tree05=(TTree*)file05->Get(openTree0);    TTree *Tree18=(TTree*)file05->Get(openTree1);    TTree *Tree31=(TTree*)file05->Get(openTree2);  
      TTree *Tree06=(TTree*)file06->Get(openTree0);    TTree *Tree19=(TTree*)file06->Get(openTree1);    TTree *Tree32=(TTree*)file06->Get(openTree2);  
      TTree *Tree07=(TTree*)file07->Get(openTree0);    TTree *Tree20=(TTree*)file07->Get(openTree1);    TTree *Tree33=(TTree*)file07->Get(openTree2);  
      TTree *Tree08=(TTree*)file08->Get(openTree0);    TTree *Tree21=(TTree*)file08->Get(openTree1);    TTree *Tree34=(TTree*)file08->Get(openTree2);  
      TTree *Tree09=(TTree*)file09->Get(openTree0);    TTree *Tree22=(TTree*)file09->Get(openTree1);    TTree *Tree35=(TTree*)file09->Get(openTree2);  
      TTree *Tree10=(TTree*)file10->Get(openTree0);    TTree *Tree23=(TTree*)file10->Get(openTree1);    TTree *Tree36=(TTree*)file10->Get(openTree2);  
      TTree *Tree11=(TTree*)file11->Get(openTree0);    TTree *Tree24=(TTree*)file11->Get(openTree1);    TTree *Tree37=(TTree*)file11->Get(openTree2);  
      TTree *Tree12=(TTree*)file12->Get(openTree0);    TTree *Tree25=(TTree*)file12->Get(openTree1);    TTree *Tree38=(TTree*)file12->Get(openTree2);  
      TTree *Tree13=(TTree*)file13->Get(openTree0);    TTree *Tree26=(TTree*)file13->Get(openTree1);    TTree *Tree39=(TTree*)file13->Get(openTree2);  
    
      char input01 [50]; sprintf(input01, "%s>>h01(%i,%f,%f)",plot,bin,min,max);
      char input02 [50]; sprintf(input02, "%s>>h02(%i,%f,%f)",plot,bin,min,max);
      char input03 [50]; sprintf(input03, "%s>>h03(%i,%f,%f)",plot,bin,min,max);
      char input04 [50]; sprintf(input04, "%s>>h04(%i,%f,%f)",plot,bin,min,max);
      char input05 [50]; sprintf(input05, "%s>>h05(%i,%f,%f)",plot,bin,min,max);
      char input06 [50]; sprintf(input06, "%s>>h06(%i,%f,%f)",plot,bin,min,max);
      char input07 [50]; sprintf(input07, "%s>>h07(%i,%f,%f)",plot,bin,min,max);
      char input08 [50]; sprintf(input08, "%s>>h08(%i,%f,%f)",plot,bin,min,max);
      char input09 [50]; sprintf(input09, "%s>>h09(%i,%f,%f)",plot,bin,min,max);
      char input10 [50]; sprintf(input10, "%s>>h10(%i,%f,%f)",plot,bin,min,max);
      char input11 [50]; sprintf(input11, "%s>>h11(%i,%f,%f)",plot,bin,min,max);
      char input12 [50]; sprintf(input12, "%s>>h12(%i,%f,%f)",plot,bin,min,max);
      char input13 [50]; sprintf(input13, "%s>>h13(%i,%f,%f)",plot,bin,min,max);
      char input14 [50]; sprintf(input14, "%s>>h14(%i,%f,%f)",plot,bin,min,max);
      char input15 [50]; sprintf(input15, "%s>>h15(%i,%f,%f)",plot,bin,min,max);
      char input16 [50]; sprintf(input16, "%s>>h16(%i,%f,%f)",plot,bin,min,max);
      char input17 [50]; sprintf(input17, "%s>>h17(%i,%f,%f)",plot,bin,min,max);
      char input18 [50]; sprintf(input18, "%s>>h18(%i,%f,%f)",plot,bin,min,max);
      char input19 [50]; sprintf(input19, "%s>>h19(%i,%f,%f)",plot,bin,min,max);
      char input20 [50]; sprintf(input20, "%s>>h20(%i,%f,%f)",plot,bin,min,max);
      char input21 [50]; sprintf(input21, "%s>>h21(%i,%f,%f)",plot,bin,min,max);
      char input22 [50]; sprintf(input22, "%s>>h22(%i,%f,%f)",plot,bin,min,max);
      char input23 [50]; sprintf(input23, "%s>>h23(%i,%f,%f)",plot,bin,min,max);
      char input24 [50]; sprintf(input24, "%s>>h24(%i,%f,%f)",plot,bin,min,max);
      char input25 [50]; sprintf(input25, "%s>>h25(%i,%f,%f)",plot,bin,min,max);
      char input26 [50]; sprintf(input26, "%s>>h26(%i,%f,%f)",plot,bin,min,max);
      char input27 [50]; sprintf(input27, "%s>>h27(%i,%f,%f)",plot,bin,min,max);
      char input28 [50]; sprintf(input28, "%s>>h28(%i,%f,%f)",plot,bin,min,max);
      char input29 [50]; sprintf(input29, "%s>>h29(%i,%f,%f)",plot,bin,min,max);
      char input30 [50]; sprintf(input30, "%s>>h30(%i,%f,%f)",plot,bin,min,max);
      char input31 [50]; sprintf(input31, "%s>>h31(%i,%f,%f)",plot,bin,min,max);
      char input32 [50]; sprintf(input32, "%s>>h32(%i,%f,%f)",plot,bin,min,max);
      char input33 [50]; sprintf(input33, "%s>>h33(%i,%f,%f)",plot,bin,min,max);
      char input34 [50]; sprintf(input34, "%s>>h34(%i,%f,%f)",plot,bin,min,max);
      char input35 [50]; sprintf(input35, "%s>>h35(%i,%f,%f)",plot,bin,min,max);
      char input36 [50]; sprintf(input36, "%s>>h36(%i,%f,%f)",plot,bin,min,max);
      char input37 [50]; sprintf(input37, "%s>>h37(%i,%f,%f)",plot,bin,min,max);
      char input38 [50]; sprintf(input38, "%s>>h38(%i,%f,%f)",plot,bin,min,max);
      char input39 [50]; sprintf(input39, "%s>>h39(%i,%f,%f)",plot,bin,min,max);
      
      Tree01->Draw(input01,CUT,"E"); if(Tree01->Draw(input01,CUT,"E")) {data      = h01; }
      Tree02->Draw(input02,CUT);     if(Tree02->Draw(input02,CUT))     {DY100     = h02; }
      Tree03->Draw(input03,CUT);     if(Tree03->Draw(input03,CUT))     {DY70      = h03; }
      Tree04->Draw(input04,CUT);     if(Tree04->Draw(input04,CUT))     {DYM50_100 = h04; }
      Tree05->Draw(input05,CUT);     if(Tree05->Draw(input05,CUT))     {DYM50_70  = h05; }
      Tree06->Draw(input06,CUT);     if(Tree06->Draw(input06,CUT))     {QCD1000   = h06; }
      Tree07->Draw(input07,CUT);     if(Tree07->Draw(input07,CUT))     {QCD250    = h07; }
      Tree08->Draw(input08,CUT);     if(Tree08->Draw(input08,CUT))     {QCD500    = h08; }
      Tree09->Draw(input09,CUT);     if(Tree09->Draw(input09,CUT))     {TT        = h09; }
      Tree10->Draw(input10,CUT);     if(Tree10->Draw(input10,CUT))     {WJets180  = h10; }
      Tree11->Draw(input11,CUT);     if(Tree11->Draw(input11,CUT))     {WW        = h11; }
      Tree12->Draw(input12,CUT);     if(Tree12->Draw(input12,CUT))     {WZ        = h12; }
      Tree13->Draw(input13,CUT);     if(Tree13->Draw(input13,CUT))     {WZ        = h13; }
      Tree14->Draw(input14,CUT,"E"); if(Tree14->Draw(input14,CUT,"E")) {data      ->Add(h14); }
      Tree15->Draw(input15,CUT);     if(Tree15->Draw(input15,CUT))     {DY100     ->Add(h15); }
      Tree16->Draw(input16,CUT);     if(Tree16->Draw(input16,CUT))     {DY70      ->Add(h16); }
      Tree17->Draw(input17,CUT);     if(Tree17->Draw(input17,CUT))     {DYM50_100 ->Add(h17); }
      Tree18->Draw(input18,CUT);     if(Tree18->Draw(input18,CUT))     {DYM50_70  ->Add(h18); }
      Tree19->Draw(input19,CUT);     if(Tree19->Draw(input19,CUT))     {QCD1000   ->Add(h19); }
      Tree20->Draw(input20,CUT);     if(Tree20->Draw(input20,CUT))     {QCD250    ->Add(h20); }
      Tree21->Draw(input21,CUT);     if(Tree21->Draw(input21,CUT))     {QCD500    ->Add(h21); }
      Tree22->Draw(input22,CUT);     if(Tree22->Draw(input22,CUT))     {TT        ->Add(h22); }
      Tree23->Draw(input23,CUT);     if(Tree23->Draw(input23,CUT))     {WJets180  ->Add(h23); }
      Tree24->Draw(input24,CUT);     if(Tree24->Draw(input24,CUT))     {WW        ->Add(h24); }
      Tree25->Draw(input25,CUT);     if(Tree25->Draw(input25,CUT))     {WZ        ->Add(h25); }
      Tree26->Draw(input26,CUT);     if(Tree26->Draw(input26,CUT))     {WZ        ->Add(h26); }
      Tree27->Draw(input27,CUT,"E"); if(Tree27->Draw(input27,CUT,"E")) {data      ->Add(h27); }
      Tree28->Draw(input28,CUT);     if(Tree28->Draw(input28,CUT))     {DY100     ->Add(h28); }
      Tree29->Draw(input29,CUT);     if(Tree29->Draw(input29,CUT))     {DY70      ->Add(h29); }
      Tree30->Draw(input30,CUT);     if(Tree30->Draw(input30,CUT))     {DYM50_100 ->Add(h30); }
      Tree31->Draw(input31,CUT);     if(Tree31->Draw(input31,CUT))     {DYM50_70  ->Add(h31); }
      Tree32->Draw(input32,CUT);     if(Tree32->Draw(input32,CUT))     {QCD1000   ->Add(h32); }
      Tree33->Draw(input33,CUT);     if(Tree33->Draw(input33,CUT))     {QCD250    ->Add(h33); }
      Tree34->Draw(input34,CUT);     if(Tree34->Draw(input34,CUT))     {QCD500    ->Add(h34); }
      Tree35->Draw(input35,CUT);     if(Tree35->Draw(input35,CUT))     {TT        ->Add(h35); }
      Tree36->Draw(input36,CUT);     if(Tree36->Draw(input36,CUT))     {WJets180  ->Add(h36); }
      Tree37->Draw(input37,CUT);     if(Tree37->Draw(input37,CUT))     {WW        ->Add(h37); }
      Tree38->Draw(input38,CUT);     if(Tree38->Draw(input38,CUT))     {WZ        ->Add(h38); }
      Tree39->Draw(input39,CUT);     if(Tree39->Draw(input39,CUT))     {WZ        ->Add(h39); }

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

      TH1F *ERR = new TH1F("","",data->GetNbinsX(),data->GetXaxis()->GetXmin(),data->GetXaxis()->GetXmax());
      for(int m=1; m<ERR->GetNbinsX()+1; m++){
	ERR->SetBinContent(m,w_DY100     *  DY100->GetBinContent(m)+
			     w_DY70      *  DY70->GetBinContent(m)+
			     w_DYM50_100 *  DYM50_100->GetBinContent(m)+
			     w_DYM50_70  *  DYM50_70->GetBinContent(m)+
			     w_QCD1000   *  1.9 * QCD1000->GetBinContent(m)+
			     w_QCD500    *  1.9 * QCD500->GetBinContent(m)+
			     w_QCD250    *  1.9 * QCD250->GetBinContent(m)+
			     w_TT        *  TT->GetBinContent(m)+
			     w_WW        *  WW->GetBinContent(m)+
			     w_WZ        *  WZ->GetBinContent(m)+
			     w_ZZ        *  ZZ->GetBinContent(m)+
			     w_WJets180  *  WJets180->GetBinContent(m)
			   ); 
	ERR->SetBinError(m,sqrt(
				w_DY100     *  w_DY100     *  DY100->GetBinContent(m)+
				w_DY70      *  w_DY70      *  DY70->GetBinContent(m)+
				w_DYM50_100 *  w_DYM50_100 *  DYM50_100->GetBinContent(m)+
				w_DYM50_70  *  w_DYM50_70  *  DYM50_70->GetBinContent(m)+
				w_QCD1000   *  w_QCD1000   *  QCD1000->GetBinContent(m)+
				w_QCD500    *  w_QCD500    *  QCD500->GetBinContent(m)+
				w_QCD250    *  w_QCD250    *  QCD250->GetBinContent(m)+
				w_TT        *  w_TT        *  TT->GetBinContent(m)+
				w_WW        *  w_WW        *  WW->GetBinContent(m)+
				w_WZ        *  w_WZ        *  WZ->GetBinContent(m)+
				w_ZZ        *  w_ZZ        *  ZZ->GetBinContent(m)+
				w_WJets180  *  w_WJets180  *  WJets180->GetBinContent(m)
				)); 
      }

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

      //hnew1 = DY100->Rebin(N-1,"hnew1",xbins);
      //hnew2 = WW->Rebin(N-1,"hnew2",xbins);
      //hnew3 = TT->Rebin(N-1,"hnew3",xbins);
      //hnew4 = WJets180->Rebin(N-1,"hnew4",xbins);
      //hnew5 = QCD1000->Rebin(N-1,"hnew5",xbins);
      //hnew6 = data->Rebin(N-1,"hnew6",xbins);
      //hnew7 = ERR->Rebin(N-1,"hnew7",xbins);

      TH1D *RATIO = new TH1D("","",ERR->GetNbinsX(),ERR->GetXaxis()->GetXmin(),ERR->GetXaxis()->GetXmax());
      for(int m=1; m<ERR->GetNbinsX()+1; m++){ 
	if(ERR->GetBinContent(m)!=0 && data->GetBinContent(m)!=0) {
	  RATIO->SetBinContent(m,data->GetBinContent(m)/ERR->GetBinContent(m));
	  RATIO->SetBinError(m,sqrt(ERR->GetBinContent(m)*ERR->GetBinContent(m)*data->GetBinError(m)*data->GetBinError(m)
				    +data->GetBinContent(m)*data->GetBinContent(m)*ERR->GetBinError(m)*ERR->GetBinError(m))/(ERR->GetBinContent(m)*ERR->GetBinContent(m)));
	}
      }
      //TH1D *RATIO = new TH1D("","",N-1,xbins);
      //for(int m=1; m<hnew7->GetNbinsX()+1; m++){ 
      //  if(hnew7->GetBinContent(m)!=0 && hnew6->GetBinContent(m)!=0) {
      //	RATIO->SetBinContent(m,hnew6->GetBinContent(m)/hnew7->GetBinContent(m));
      //	RATIO->SetBinError(m,sqrt(hnew7->GetBinContent(m)*hnew7->GetBinContent(m)*hnew6->GetBinError(m)*hnew6->GetBinError(m)
      //				  +hnew6->GetBinContent(m)*hnew6->GetBinContent(m)*hnew7->GetBinError(m)*hnew7->GetBinError(m))/
      //			   (hnew7->GetBinContent(m)*hnew7->GetBinContent(m)));
      //  }
      //}
      RATIO->Draw("E");
      RATIO->SetMaximum(3);
      RATIO->SetMinimum(-1);
      RATIO->SetLineColor(1);
      RATIO->SetLineWidth(2);
      RATIO->GetXaxis()->SetTitleOffset(0.9);
      RATIO->GetXaxis()->SetLabelSize(0.08);
      RATIO->GetXaxis()->SetTitleSize(0.09);
      RATIO->GetYaxis()->SetTitleOffset(0.35);
      RATIO->GetYaxis()->SetTitle("data/MC");
      RATIO->GetXaxis()->SetTitle(axis);
      RATIO->GetYaxis()->SetLabelSize(0.08);
      RATIO->GetYaxis()->SetTitleSize(0.09);
	
      TLine* line = new TLine(ERR->GetXaxis()->GetXmin(),1,ERR->GetXaxis()->GetXmax(),1);
      line->SetLineColor(2);
      line->SetLineWidth(2);
	
      //plots 
      TCanvas* c1 = new TCanvas("c1","c1",0,0,600,600);

      TLatex latexLabel1;
      latexLabel1.SetTextSize(0.04);
      latexLabel1.SetNDC();
      latexLabel1.DrawLatex(0.08, 0.96, "CMS");	
      TLatex latexLabel2;
      latexLabel2.SetTextSize(0.04);
      latexLabel2.SetTextFont(42);
      latexLabel2.SetNDC();
      latexLabel2.DrawLatex(0.57, 0.96, "L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");

      TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.98,0.32);
      c1_1->Draw();
      c1_1->cd();
      c1_1->SetTopMargin(0.05);
      c1_1->SetBottomMargin(0.2);
      c1_1->SetRightMargin(0.02);
      c1_1->SetLeftMargin(0.07);
      RATIO->Draw("E");
      line->Draw("same");	
      c1->cd();
	
      TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.98,0.96);
      c1_2->Draw();
      c1_2->cd();
      c1_2->SetTopMargin(0.01);
      c1_2->SetBottomMargin(0.1);
      c1_2->SetRightMargin(0.02);
      c1_2->SetLeftMargin(0.07);
      THStack *hs = new THStack("hs","hs");
      hs->SetTitle();
      //hs->Add(hnew1);
      //hs->Add(hnew2);
      //hs->Add(hnew3);
      //hs->Add(hnew4);
      //hs->Add(hnew5);
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
      hs->GetYaxis()->SetTitleOffset(0.7); 
      hs->GetXaxis()->SetTitle(axis); 
      hs->SetMinimum(0);
      hs->SetMaximum(maxy);
      hs->GetYaxis()->SetTitle(TString("Events / ")+TString::Format("%.2f",(max-min)/bin));
      data->SetLineWidth(2); 
      data->SetLineColor(1);
      data->SetMarkerStyle(20); 
      data->SetMarkerSize(1.3); 
      data->Draw("Esame");
      ERR->SetFillStyle(3005);
      ERR->SetFillColor(12);
      ERR->SetLineColor(12);
      ERR->Draw("E2same");
      //hs->GetYaxis()->SetTitle("Events");
      //hnew6->SetLineWidth(2); 
      //hnew6->SetLineColor(1);
      //hnew6->SetMarkerStyle(20); 
      //hnew6->SetMarkerSize(1.3); 
      //hnew6->Draw("Esame");
      //hnew7->SetFillStyle(3005);
      //hnew7->SetFillColor(12);
      //hnew7->SetLineColor(12);
      //hnew7->Draw("E2same");
  
      TLegend *pl2 = new TLegend(0.79,0.70,0.97,0.98);
      pl2->SetTextSize(0.035); 
      pl2->SetFillColor(0);
      TLegendEntry *ple2 = pl2->AddEntry(DY100, "Drell-Yan",  "F");
      ple2 = pl2->AddEntry(WW, "Diboson",  "F");
      ple2 = pl2->AddEntry(TT, "ttbar",  "F");
      ple2 = pl2->AddEntry(WJets180, "WJets (pt180)",  "F");
      ple2 = pl2->AddEntry(QCD1000, "QCD",  "F");
      ple2 = pl2->AddEntry(data, "data",  "LP");
      pl2->Draw();

      if(log) {
	hs->SetMinimum(0.5);
	c1_2->SetLogy();
      }
      if(save) c1->SaveAs("DataMC_"+CHANNEL+"_"+name+".pdf");
      if(save) c1->SaveAs("DataMC_"+CHANNEL+"_"+name+".png");
    }
  }
}
