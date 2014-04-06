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

  vector<string> PLOT;                     vector<int> BIN;    vector<float> MIN;   vector<float> MAX;   vector<float> MAXY;   vector<TString> AXIS;/*
  PLOT.push_back("MassCA");          BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(300);  MAXY.push_back(6);    AXIS.push_back("M(lep1,lep1) [GeV] (Collinear Approximation)");    	
  PLOT.push_back("MassEff");         BIN.push_back(35);  MIN.push_back(0);    MAX.push_back(350);  MAXY.push_back(6);    AXIS.push_back("M(lep1,lep1) [GeV] (effective)");
  PLOT.push_back("MassSvfit");       BIN.push_back(25);  MIN.push_back(0);    MAX.push_back(250);  MAXY.push_back(40);   AXIS.push_back("M(lep1,lep1) [GeV] (SvFit)");
  PLOT.push_back("MassVis");         BIN.push_back(25);  MIN.push_back(0);    MAX.push_back(250);  MAXY.push_back(4000);   AXIS.push_back("M(lep1,lep1) [GeV] (visible)");
  PLOT.push_back("NVertices");       BIN.push_back(41);  MIN.push_back(-0.5); MAX.push_back(40.5); MAXY.push_back(1000);   AXIS.push_back("N(vertices)");
  PLOT.push_back("XMassCA");         BIN.push_back(50);  MIN.push_back(0);    MAX.push_back(2500); MAXY.push_back(8);    AXIS.push_back("M(Z,H) [GeV] (Collinear Approximation)");  
  PLOT.push_back("XMassEff");        BIN.push_back(50);  MIN.push_back(0);    MAX.push_back(2500); MAXY.push_back(10);   AXIS.push_back("M(Z,H) [GeV] (effective)");
  PLOT.push_back("XMassSVFit");      BIN.push_back(50);  MIN.push_back(0);    MAX.push_back(2500); MAXY.push_back(20);   AXIS.push_back("M(Z,H) [GeV] (SvFit)");
  PLOT.push_back("XMassVis");        BIN.push_back(50);  MIN.push_back(0);    MAX.push_back(2500); MAXY.push_back(10);   AXIS.push_back("M(Z,H) [GeV] (visible)");
  PLOT.push_back("dPhiJetMet");      BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(20);   AXIS.push_back("#Delta#phi(jet,met)");  
  PLOT.push_back("dPhiLepMet1");     BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(20);   AXIS.push_back("#Delta#phi(lep1,met)");
  PLOT.push_back("dPhiLepMet2");     BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(15);   AXIS.push_back("#Delta#phi(lep2,met)");
  PLOT.push_back("dRJetLep1");       BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(20);   AXIS.push_back("#DeltaR(jet,lep1)");
  PLOT.push_back("dRJetLep2");       BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(15);   AXIS.push_back("#DeltaR(jet,lep2)");
  PLOT.push_back("dRJetMet");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(20);   AXIS.push_back("#DeltaR(jet,met)");
  PLOT.push_back("dRLepLep");        BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(15);   AXIS.push_back("#DeltaR(lep1,lep2)");
  PLOT.push_back("dRLepMet1");       BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(12);   AXIS.push_back("#DeltaR(lep1,met)");
  PLOT.push_back("dRLepMet2");       BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(10);   AXIS.push_back("#DeltaR(lep2,met)");
  PLOT.push_back("dRZZCA");          BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(20);   AXIS.push_back("#DeltaR(Z,H) (Collinear Approximation)");
  PLOT.push_back("dRZZEff");         BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(30);   AXIS.push_back("#DeltaR(Z,H) (effective)");
  PLOT.push_back("dRZZSvFit");       BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(30);   AXIS.push_back("#DeltaR(Z,H) (SvFit)");
  PLOT.push_back("dRZZVis");         BIN.push_back(30);  MIN.push_back(0);    MAX.push_back(6);    MAXY.push_back(30);   AXIS.push_back("#DeltaR(Z,H) (visible)");
  PLOT.push_back("jetEta");          BIN.push_back(32);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(12);   AXIS.push_back("jet #eta");
  PLOT.push_back("jetMass");         BIN.push_back(25);  MIN.push_back(20);   MAX.push_back(70);   MAXY.push_back(6000);   AXIS.push_back("jet mass [GeV]");*/
  PLOT.push_back("jetPt");           BIN.push_back(30);  MIN.push_back(400);  MAX.push_back(1000); MAXY.push_back(5000);   AXIS.push_back("jet pt [GeV]");/*
  PLOT.push_back("jetSubjettiness"); BIN.push_back(20);  MIN.push_back(0);    MAX.push_back(1);    MAXY.push_back(20);   AXIS.push_back("jet #tau_{21}");
  PLOT.push_back("lepDetIso1");      BIN.push_back(15);  MIN.push_back(0);    MAX.push_back(0.15); MAXY.push_back(25);   AXIS.push_back("lep1 DetIso");
  PLOT.push_back("lepDetIso2");      BIN.push_back(15);  MIN.push_back(0);    MAX.push_back(0.15); MAXY.push_back(15);   AXIS.push_back("lep2 DetIso");
  PLOT.push_back("lepEta1");         BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(20);   AXIS.push_back("lep1 #eta");
  PLOT.push_back("lepEta2");         BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(15);   AXIS.push_back("lep2 #eta");
  PLOT.push_back("lepNumberEle");    BIN.push_back(4);   MIN.push_back(-0.5); MAX.push_back(3.5);  MAXY.push_back(80);   AXIS.push_back("Number of Electron");
  PLOT.push_back("lepNumberMuo");    BIN.push_back(4);   MIN.push_back(-0.5); MAX.push_back(3.5);  MAXY.push_back(80);   AXIS.push_back("Number of Muon");
  PLOT.push_back("lepPFIso1");       BIN.push_back(15);  MIN.push_back(0);    MAX.push_back(0.15); MAXY.push_back(30);   AXIS.push_back("lep1 PFIso");
  PLOT.push_back("lepPFIso2");       BIN.push_back(15);  MIN.push_back(0);    MAX.push_back(0.15); MAXY.push_back(30);   AXIS.push_back("lep2 PFIso");
  PLOT.push_back("lepPt1");          BIN.push_back(35);  MIN.push_back(0);    MAX.push_back(700);  MAXY.push_back(15);   AXIS.push_back("lep1 pt [GeV]");
  PLOT.push_back("lepPt2");          BIN.push_back(25);  MIN.push_back(0);    MAX.push_back(500);  MAXY.push_back(20);   AXIS.push_back("lep2 pt [GeV]");
  PLOT.push_back("met");             BIN.push_back(25);  MIN.push_back(50);   MAX.push_back(500);  MAXY.push_back(30);   AXIS.push_back("met [GeV]");
  PLOT.push_back("metPhi");          BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(10);   AXIS.push_back("met #phi");
  PLOT.push_back("metPx");           BIN.push_back(25);  MIN.push_back(-100); MAX.push_back(100);  MAXY.push_back(8);    AXIS.push_back("met_x [GeV]");
  PLOT.push_back("metPy");           BIN.push_back(25);  MIN.push_back(-100); MAX.push_back(100);  MAXY.push_back(8);    AXIS.push_back("met_y [GeV]");
  PLOT.push_back("nbtagsL");         BIN.push_back(7);   MIN.push_back(-0.5); MAX.push_back(6.5);  MAXY.push_back(25);   AXIS.push_back("N(loose btag)");
  PLOT.push_back("nbtagsM");         BIN.push_back(7);   MIN.push_back(-0.5); MAX.push_back(6.5);  MAXY.push_back(50000);   AXIS.push_back("N(medium btag)");
  PLOT.push_back("nbtagsT");         BIN.push_back(7);   MIN.push_back(-0.5); MAX.push_back(6.5);  MAXY.push_back(50);   AXIS.push_back("N(tight btag)");
  PLOT.push_back("uncorrmet");       BIN.push_back(25);  MIN.push_back(30);   MAX.push_back(500);  MAXY.push_back(30);   AXIS.push_back("uncorrected met [GeV]");
  PLOT.push_back("uncorrmetPhi");    BIN.push_back(20);  MIN.push_back(-4);   MAX.push_back(4);    MAXY.push_back(10);   AXIS.push_back("uncorrected met #phi");*/

  for(int i=0; i<PLOT.size(); i++){
  //for(int i=3; i<4; i++){
    char *plot = PLOT[i].c_str();
    TString name = PLOT[i];
    int bin=BIN[i]; 
    float min=MIN[i]; 
    float max=MAX[i]; 
    float maxy=MAXY[i];
    TString axis = AXIS[i];
    bool save=false; 
    bool log=true;


    char channel  [500]; sprintf(channel,  "MuoTau"); 
    char demo     [500]; sprintf(demo,     "demo/Tree"); 
    char openTree [500]; sprintf(openTree, "%s%s",demo,channel); 
    TString CHANNEL = channel;
    char CUT [500]; 

    sprintf(CUT, "weight*trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)");

    cout<<endl;
    cout<<"Plotting plot "<<i+1<<"/"<<PLOT.size()<<", i.e. "<<name<<endl;

    TH1F *data1    = new TH1F("","",bin,min,max);
    TH1F *data2    = new TH1F("","",bin,min,max);
    TH1F *data3    = new TH1F("","",bin,min,max);
    TH1F *data4    = new TH1F("","",bin,min,max);
    TH1F *DY100    = new TH1F("","",bin,min,max);
    TH1F *DY70     = new TH1F("","",bin,min,max);											     
    TH1F *QCD1000  = new TH1F("","",bin,min,max);
    TH1F *QCD250   = new TH1F("","",bin,min,max);
    TH1F *QCD500   = new TH1F("","",bin,min,max);
    TH1F *TT       = new TH1F("","",bin,min,max);
    TH1F *WJets180 = new TH1F("","",bin,min,max);
    TH1F *WW       = new TH1F("","",bin,min,max);
    TH1F *WZ       = new TH1F("","",bin,min,max);
    TH1F *ZZ       = new TH1F("","",bin,min,max);
  
    TFile *file1  = TFile::Open("../RISULTATI/05_260114/RunA_QCD.root");     TTree *Tree1  = (TTree*)file1->Get(openTree); 
    TFile *file2  = TFile::Open("../RISULTATI/05_260114/RunB_QCD.root");     TTree *Tree2  = (TTree*)file2->Get(openTree); 
    TFile *file3  = TFile::Open("../RISULTATI/05_260114/RunC_QCD.root");     TTree *Tree3  = (TTree*)file3->Get(openTree); 
    TFile *file4  = TFile::Open("../RISULTATI/05_260114/RunD_QCD.root");     TTree *Tree4  = (TTree*)file4->Get(openTree); 
    TFile *file5  = TFile::Open("../RISULTATI/05_260114/DY100_QCD.root");    TTree *Tree5  = (TTree*)file5->Get(openTree); 
    TFile *file6  = TFile::Open("../RISULTATI/05_260114/DY70_QCD.root");     TTree *Tree6  = (TTree*)file6->Get(openTree); 
    TFile *file7  = TFile::Open("../RISULTATI/05_260114/QCD1000_QCD.root");  TTree *Tree7  = (TTree*)file7->Get(openTree);
    TFile *file8  = TFile::Open("../RISULTATI/05_260114/QCD250_QCD.root");   TTree *Tree8  = (TTree*)file8->Get(openTree); 
    TFile *file9  = TFile::Open("../RISULTATI/05_260114/QCD500_QCD.root");   TTree *Tree9  = (TTree*)file9->Get(openTree); 
    TFile *file10 = TFile::Open("../RISULTATI/05_260114/TT_QCD.root");       TTree *Tree10 = (TTree*)file10->Get(openTree); 
    TFile *file11 = TFile::Open("../RISULTATI/05_260114/WJets180_QCD.root"); TTree *Tree11 = (TTree*)file11->Get(openTree);
    TFile *file12 = TFile::Open("../RISULTATI/05_260114/WW_QCD.root");       TTree *Tree12 = (TTree*)file12->Get(openTree); 
    TFile *file13 = TFile::Open("../RISULTATI/05_260114/WZ_QCD.root");       TTree *Tree13 = (TTree*)file13->Get(openTree); 
    TFile *file14 = TFile::Open("../RISULTATI/05_260114/ZZ_QCD.root");       TTree *Tree14 = (TTree*)file14->Get(openTree); 

    char input1  [50]; sprintf(input1,  "%s>>h1(%i,%f,%f)", plot,bin,min,max);
    char input2  [50]; sprintf(input2,  "%s>>h2(%i,%f,%f)", plot,bin,min,max);
    char input3  [50]; sprintf(input3,  "%s>>h3(%i,%f,%f)", plot,bin,min,max);
    char input4  [50]; sprintf(input4,  "%s>>h4(%i,%f,%f)", plot,bin,min,max);
    char input5  [50]; sprintf(input5,  "%s>>h5(%i,%f,%f)", plot,bin,min,max);
    char input6  [50]; sprintf(input6,  "%s>>h6(%i,%f,%f)", plot,bin,min,max);
    char input7  [50]; sprintf(input7,  "%s>>h7(%i,%f,%f)", plot,bin,min,max);
    char input8  [50]; sprintf(input8,  "%s>>h8(%i,%f,%f)", plot,bin,min,max);
    char input9  [50]; sprintf(input9,  "%s>>h9(%i,%f,%f)", plot,bin,min,max);
    char input10 [50]; sprintf(input10, "%s>>h10(%i,%f,%f)",plot,bin,min,max);
    char input11 [50]; sprintf(input11, "%s>>h11(%i,%f,%f)",plot,bin,min,max);
    char input12 [50]; sprintf(input12, "%s>>h12(%i,%f,%f)",plot,bin,min,max);
    char input13 [50]; sprintf(input13, "%s>>h13(%i,%f,%f)",plot,bin,min,max);
    char input14 [50]; sprintf(input14, "%s>>h14(%i,%f,%f)",plot,bin,min,max);

    Tree1->Draw(input1,CUT,"E"); 
    if(Tree1->Draw(input1,CUT,"E")) data1 = h1;
    Tree2->Draw(input2,CUT,"E"); 
    if(Tree2->Draw(input2,CUT,"E")) data2 = h2;
    Tree3->Draw(input3,CUT,"E"); 
    if(Tree3->Draw(input3,CUT,"E")) data3 = h3;
    Tree4->Draw(input4,CUT,"E"); 
    if(Tree4->Draw(input4,CUT,"E")) data4 = h4;
    Tree5->Draw(input5,CUT);     
    if(Tree5->Draw(input5,CUT))     DY100 = h5;
    Tree6->Draw(input6,CUT);         
    if(Tree6->Draw(input6,CUT))     DY70 = h6; 
    Tree7->Draw(input7,CUT);         
    if(Tree7->Draw(input7,CUT))     QCD1000 = h7;
    Tree8->Draw(input8,CUT);         
    if(Tree8->Draw(input8,CUT))     QCD250 = h8;
    Tree9->Draw(input9,CUT);         
    if(Tree9->Draw(input9,CUT))     QCD500 = h9;
    Tree10->Draw(input10,CUT);      
    if(Tree10->Draw(input10,CUT))   TT = h10; 
    Tree11->Draw(input11,CUT);      
    if(Tree11->Draw(input11,CUT))   WJets180 = h11;
    Tree12->Draw(input12,CUT);      
    if(Tree12->Draw(input12,CUT))   WW = h12;
    Tree13->Draw(input13,CUT);      
    if(Tree13->Draw(input13,CUT))   WZ = h13;
    Tree14->Draw(input14,CUT);      
    if(Tree14->Draw(input14,CUT))   ZZ = h14;
  
    DY100->Sumw2();
    DY70->Sumw2();
    QCD1000->Sumw2();
    QCD500->Sumw2();
    TT->Sumw2();
    WW->Sumw2();
    WZ->Sumw2();
    ZZ->Sumw2();
    WJets180->Sumw2();

    DY100->Scale(39.1/32.9);
    DY70->Scale(62.9/53.0);
    data1->Add(data2);
    data1->Add(data3);
    data1->Add(data4);
    DY100->Add(DY70);
    QCD1000->Add(QCD250);
    QCD1000->Add(QCD500);
    WW->Add(WZ);
    WW->Add(ZZ);
    if(CHANNEL=="MuoTau") QCD1000->Scale(1.7);
    if(CHANNEL=="EleTau") QCD1000->Scale(2.0);

    DY100->SetFillColor(kBlue-9);
    QCD1000->SetFillColor(kRed-2);
    WW->SetFillColor(kGreen+1);
    TT->SetFillColor(kBlue);
    WJets180->SetFillColor(kMagenta-5);

    cout<<"data    = "<<data1->Integral()<<endl;
    cout<<"DY      = "<<DY100->Integral()<<endl;
    cout<<"diboson = "<<WW->Integral()<<endl;
    cout<<"TT      = "<<TT->Integral()<<endl;
    cout<<"WJets   = "<<WJets180->Integral()<<endl;
    cout<<"QCD     = "<<QCD1000->Integral()<<endl;

    TH1F *ERR = new TH1F("","",data1->GetNbinsX(),data1->GetXaxis()->GetXmin(),data1->GetXaxis()->GetXmax());
    ERR->Sumw2();
    ERR->Add(DY100);
    ERR->Add(WW);
    ERR->Add(TT);
    ERR->Add(WJets180);
    ERR->Add(QCD1000);
    TH1D *RATIO = new TH1D("","",ERR->GetNbinsX(),ERR->GetXaxis()->GetXmin(),ERR->GetXaxis()->GetXmax());
    for(int m=1; m<ERR->GetNbinsX()+1; m++){ 
      if(ERR->GetBinContent(m)!=0 && data1->GetBinContent(m)!=0) {
	RATIO->SetBinContent(m,data1->GetBinContent(m)/ERR->GetBinContent(m));
	RATIO->SetBinError(m,sqrt(ERR->GetBinContent(m)*ERR->GetBinContent(m)*data1->GetBinError(m)*data1->GetBinError(m)
				  +data1->GetBinContent(m)*data1->GetBinContent(m)*ERR->GetBinError(m)*ERR->GetBinError(m))/(ERR->GetBinContent(m)*ERR->GetBinContent(m)));
      }
    }
    RATIO->Draw("E");
    RATIO->SetMaximum(5);
    RATIO->SetMinimum(-3);
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
    latexLabel2.DrawLatex(0.57, 0.96, "L = 18.9 fb^{-1} at #sqrt{s} = 8 TeV");

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
    hs->GetYaxis()->SetTitle(TString("Events / ")+TString::Format("%.2f",(max-min)/bin));
    hs->GetXaxis()->SetTitle(axis); 
    hs->SetMinimum(0);
    hs->SetMaximum(maxy);
    data1->SetLineWidth(2); 
    data1->SetLineColor(1);
    data1->SetMarkerStyle(20); 
    data1->SetMarkerSize(1.3); 
    data1->Draw("Esame");
    ERR->SetFillStyle(3005);
    ERR->SetFillColor(12);
    ERR->SetLineColor(12);
    ERR->Draw("E2same");
  
    TLegend *pl2 = new TLegend(0.79,0.70,0.97,0.98);
    pl2->SetTextSize(0.035); 
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(DY100, "Drell-Yan",  "F");
    ple2 = pl2->AddEntry(WW, "Diboson",  "F");
    ple2 = pl2->AddEntry(TT, "ttbar",  "F");
    ple2 = pl2->AddEntry(WJets180, "WJets (pt180)",  "F");
    ple2 = pl2->AddEntry(QCD1000, "QCD",  "F");
    ple2 = pl2->AddEntry(data1, "data",  "LP");
    pl2->Draw();

    if(log) {
      hs->SetMinimum(0.05);
      c1_2->SetLogy();
    }
    if(save){
      c1->SaveAs(CHANNEL+"_SB_"+name+"_LepIso.png");
      //c1->SaveAs(CHANNEL+"_SB_"+name+".root");
      //c1->SaveAs(CHANNEL+"_SB_"+name+".C");
    }
  }
}
