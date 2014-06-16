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

  vector<string> PLOT;              vector<int> BIN;   vector<float> MIN;  vector<float> MAX;   vector<float> MAXY;   vector<TString> AXIS;
  PLOT.push_back("XMassSVFit");     BIN.push_back(10); MIN.push_back(0);   MAX.push_back(3000); MAXY.push_back(90);   AXIS.push_back("M(Z,H) [GeV] ");

  char channel[500];  
  for(int type=0; type<5; type++){
    if(type==0) sprintf(channel,   "EleMuo"); 
    if(type==1) sprintf(channel,   "MuoMuo"); 
    if(type==2) sprintf(channel,   "EleEle"); 
    if(type==3) sprintf(channel,   "MuoTau"); 
    if(type==4) sprintf(channel,   "EleTau"); 
    
    TFile *file01 = TFile::Open("../../RISULTATI/analyzer_290514/data.root");
    TFile *file02 = TFile::Open("../../RISULTATI/analyzer_290514/DY100.root");
    TFile *file03 = TFile::Open("../../RISULTATI/analyzer_290514/DY70.root");
    TFile *file04 = TFile::Open("../../RISULTATI/analyzer_290514/DYM50_100.root");
    TFile *file05 = TFile::Open("../../RISULTATI/analyzer_290514/DYM50_70.root");
    TFile *file06 = TFile::Open("../../RISULTATI/analyzer_290514/QCD1000.root");
    TFile *file07 = TFile::Open("../../RISULTATI/analyzer_290514/QCD250.root");
    TFile *file08 = TFile::Open("../../RISULTATI/analyzer_290514/QCD500.root");
    TFile *file09 = TFile::Open("../../RISULTATI/analyzer_290514/TT.root");
    TFile *file10 = TFile::Open("../../RISULTATI/analyzer_290514/WJetsHT.root");
    TFile *file11 = TFile::Open("../../RISULTATI/analyzer_290514/WW.root");
    TFile *file12 = TFile::Open("../../RISULTATI/analyzer_290514/WZ.root");
    TFile *file13 = TFile::Open("../../RISULTATI/analyzer_290514/ZZ.root");

    char *plot = PLOT[0].c_str();
    TString name = PLOT[0];
    int bin=BIN[0]; 
    float min=MIN[0]; 
    float max=MAX[0]; 
    float maxy=MAXY[0];
    TString axis = AXIS[0];
      
    TH1F *SB1           = new TH1F("","",bin,min,max);
    TH1F *SB1_data      = new TH1F("","",bin,min,max);
    TH1F *SB1_DY100     = new TH1F("","",bin,min,max);
    TH1F *SB1_DY70      = new TH1F("","",bin,min,max);
    TH1F *SB1_DYM50_100 = new TH1F("","",bin,min,max);
    TH1F *SB1_DYM50_70  = new TH1F("","",bin,min,max);											     
    TH1F *SB1_QCD1000   = new TH1F("","",bin,min,max);
    TH1F *SB1_QCD250    = new TH1F("","",bin,min,max);
    TH1F *SB1_QCD500    = new TH1F("","",bin,min,max);
    TH1F *SB1_TT        = new TH1F("","",bin,min,max);
    TH1F *SB1_WJetsHT   = new TH1F("","",bin,min,max);
    TH1F *SB1_WW        = new TH1F("","",bin,min,max);
    TH1F *SB1_WZ        = new TH1F("","",bin,min,max);
    TH1F *SB1_ZZ        = new TH1F("","",bin,min,max);
    TH1F *SB2           = new TH1F("","",bin,min,max);
    TH1F *SB2_data      = new TH1F("","",bin,min,max);
    TH1F *SB2_DY100     = new TH1F("","",bin,min,max);
    TH1F *SB2_DY70      = new TH1F("","",bin,min,max);
    TH1F *SB2_DYM50_100 = new TH1F("","",bin,min,max);
    TH1F *SB2_DYM50_70  = new TH1F("","",bin,min,max);											     
    TH1F *SB2_QCD1000   = new TH1F("","",bin,min,max);
    TH1F *SB2_QCD250    = new TH1F("","",bin,min,max);
    TH1F *SB2_QCD500    = new TH1F("","",bin,min,max);
    TH1F *SB2_TT        = new TH1F("","",bin,min,max);
    TH1F *SB2_WJetsHT   = new TH1F("","",bin,min,max);
    TH1F *SB2_WW        = new TH1F("","",bin,min,max);
    TH1F *SB2_WZ        = new TH1F("","",bin,min,max);
    TH1F *SB2_ZZ        = new TH1F("","",bin,min,max);
    
    char BTAG[1000]; sprintf(BTAG, "");
    if(bTagCut==1) {sprintf(BTAG, " && nbtagsL1<1");}
    if(bTagCut==2) {sprintf(BTAG, " && nbtagsL1<2");}
    if(bTagCut==3) {sprintf(BTAG, " && nbtagsM1<1");}
    if(bTagCut==4) {sprintf(BTAG, " && nbtagsM1<2");}
    if(bTagCut==5) {sprintf(BTAG, " && nbtagsT1<1");}
    if(bTagCut==6) {sprintf(BTAG, " && nbtagsT1<2");}


    //SB1
    char demo[500];     sprintf(demo,"demo/TreeSB1");
    char openTree[500]; sprintf(openTree, "%s%s",demo,channel); 
    char CUT [1000]; 
    sprintf(CUT, "PUWeight*(trigger==1 && jetMass>40 && MassSvfit<%f %s && met>%f && MassVis>%f && dRLep1Lep2<%f)",MassSvfitCut,BTAG,METCut,MassVisCut,deltaRCut);
    TTree *Tree01 = (TTree*)file01->Get(openTree);  char input01[50]; sprintf(input01, "%s>>h01(%i,%f,%f)", plot,bin,min,max);
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
    Tree01->Draw(input01,CUT); if(Tree01->Draw(input01,CUT)) {SB1_data->Add(h01);}
    Tree02->Draw(input02,CUT); if(Tree02->Draw(input02,CUT)) {SB1_DY100->Add(h02);}
    Tree03->Draw(input03,CUT); if(Tree03->Draw(input03,CUT)) {SB1_DY70->Add(h03);}
    Tree04->Draw(input04,CUT); if(Tree04->Draw(input04,CUT)) {SB1_DYM50_100->Add(h04);}
    Tree05->Draw(input05,CUT); if(Tree05->Draw(input05,CUT)) {SB1_DYM50_70->Add(h05);}
    Tree06->Draw(input06,CUT); if(Tree06->Draw(input06,CUT)) {SB1_QCD1000->Add(h06);}
    Tree07->Draw(input07,CUT); if(Tree07->Draw(input07,CUT)) {SB1_QCD250->Add(h07);}
    Tree08->Draw(input08,CUT); if(Tree08->Draw(input08,CUT)) {SB1_QCD500->Add(h08);}
    Tree09->Draw(input09,CUT); if(Tree09->Draw(input09,CUT)) {SB1_TT->Add(h09);}
    Tree10->Draw(input10,CUT); if(Tree10->Draw(input10,CUT)) {SB1_WJetsHT->Add(h10);}
    Tree11->Draw(input11,CUT); if(Tree11->Draw(input11,CUT)) {SB1_WW->Add(h11);}
    Tree12->Draw(input12,CUT); if(Tree12->Draw(input12,CUT)) {SB1_WZ->Add(h12);}
    Tree13->Draw(input13,CUT); if(Tree13->Draw(input13,CUT)) {SB1_ZZ->Add(h13);}

    //SB2
    sprintf(demo,"demo/TreeSB2");
    sprintf(openTree, "%s%s",demo,channel);
    sprintf(CUT, "PUWeight*(trigger==1 && jetMass>140 && MassSvfit<%f %s && met>%f && MassVis>%f && dRLep1Lep2<%f)",MassSvfitCut,BTAG,METCut,MassVisCut,deltaRCut);
    TTree *Tree14 = (TTree*)file01->Get(openTree);  char input14[50]; sprintf(input14, "%s>>h14(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree15 = (TTree*)file02->Get(openTree);  char input15[50]; sprintf(input15, "%s>>h15(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree16 = (TTree*)file03->Get(openTree);  char input16[50]; sprintf(input16, "%s>>h16(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree17 = (TTree*)file04->Get(openTree);  char input17[50]; sprintf(input17, "%s>>h17(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree18 = (TTree*)file05->Get(openTree);  char input18[50]; sprintf(input18, "%s>>h18(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree19 = (TTree*)file06->Get(openTree);  char input19[50]; sprintf(input19, "%s>>h19(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree20 = (TTree*)file07->Get(openTree);  char input20[50]; sprintf(input20, "%s>>h20(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree21 = (TTree*)file08->Get(openTree);  char input21[50]; sprintf(input21, "%s>>h21(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree22 = (TTree*)file09->Get(openTree);  char input22[50]; sprintf(input22, "%s>>h22(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree23 = (TTree*)file10->Get(openTree);  char input23[50]; sprintf(input23, "%s>>h23(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree24 = (TTree*)file11->Get(openTree);  char input24[50]; sprintf(input24, "%s>>h24(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree25 = (TTree*)file12->Get(openTree);  char input25[50]; sprintf(input25, "%s>>h25(%i,%f,%f)", plot,bin,min,max);
    TTree *Tree26 = (TTree*)file13->Get(openTree);  char input26[50]; sprintf(input26, "%s>>h26(%i,%f,%f)", plot,bin,min,max);
    Tree14->Draw(input14,CUT); if(Tree14->Draw(input14,CUT)) {SB2_data->Add(h14);}
    Tree15->Draw(input15,CUT); if(Tree15->Draw(input15,CUT)) {SB2_DY100->Add(h15);}
    Tree16->Draw(input16,CUT); if(Tree16->Draw(input16,CUT)) {SB2_DY70->Add(h16);}
    Tree17->Draw(input17,CUT); if(Tree17->Draw(input17,CUT)) {SB2_DYM50_100->Add(h17);}
    Tree18->Draw(input18,CUT); if(Tree18->Draw(input18,CUT)) {SB2_DYM50_70->Add(h18);}
    Tree19->Draw(input19,CUT); if(Tree19->Draw(input19,CUT)) {SB2_QCD1000->Add(h19);}
    Tree20->Draw(input20,CUT); if(Tree20->Draw(input20,CUT)) {SB2_QCD250->Add(h20);}
    Tree21->Draw(input21,CUT); if(Tree21->Draw(input21,CUT)) {SB2_QCD500->Add(h21);}
    Tree22->Draw(input22,CUT); if(Tree22->Draw(input22,CUT)) {SB2_TT->Add(h22);}
    Tree23->Draw(input23,CUT); if(Tree23->Draw(input23,CUT)) {SB2_WJetsHT->Add(h23);}
    Tree24->Draw(input24,CUT); if(Tree24->Draw(input24,CUT)) {SB2_WW->Add(h24);}
    Tree25->Draw(input25,CUT); if(Tree25->Draw(input25,CUT)) {SB2_WZ->Add(h25);}
    Tree26->Draw(input26,CUT); if(Tree26->Draw(input26,CUT)) {SB2_ZZ->Add(h26);}
      
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
    float w_WJetsHT   = ( 25.220*19702./4971847.0);
      
    SB1->GetYaxis()->SetTitleSize(0.045);
    SB1->GetXaxis()->SetTitleSize(0.045);
    SB1->GetYaxis()->SetLabelSize(0.045);
    SB1->GetXaxis()->SetLabelSize(0.045);
    SB1->GetYaxis()->SetTitleOffset(0.9); 
    SB1->GetYaxis()->SetTitle(TString("Normalized Events / ")+TString::Format("%.2f",(max-min)/bin));
    SB1->GetXaxis()->SetTitle(axis); 
    SB1->SetLineWidth(2);
    SB1->SetLineColor(1);
      
    SB2->GetYaxis()->SetTitleSize(0.045);
    SB2->GetXaxis()->SetTitleSize(0.045);
    SB2->GetYaxis()->SetLabelSize(0.045);
    SB2->GetXaxis()->SetLabelSize(0.045);
    SB2->GetYaxis()->SetTitleOffset(0.9); 
    SB2->GetYaxis()->SetTitle(TString("Normalized Events / ")+TString::Format("%.2f",(max-min)/bin));
    SB2->GetXaxis()->SetTitle(axis); 
    SB2->SetLineWidth(2);
    SB2->SetLineColor(2);

    for(int i=0; i<SB1->GetXaxis()->GetNbins(); i++){
      //float Bkg1 = (w_DY100*SB1_DY100->GetBinContent(i) + w_DY70*SB1_DY70->GetBinContent(i) + w_DYM50_100*SB1_DYM50_100->GetBinContent(i) + 
      //		   w_DYM50_70*SB1_DYM50_70->GetBinContent(i) + 1.9*w_QCD250*SB1_QCD250->GetBinContent(i) + 1.9*w_QCD500*SB1_QCD500->GetBinContent(i) + 
      //		   1.9*w_QCD1000*SB1_QCD1000->GetBinContent(i) + w_TT*SB1_TT->GetBinContent(i) + w_WW*SB1_WW->GetBinContent(i) + w_WZ*SB1_WZ->GetBinContent(i) + 
      //		   w_ZZ*SB1_ZZ->GetBinContent(i) + w_WJetsHT*SB1_WJetsHT->GetBinContent(i)
      //		   );
      //float BkgErr1 = sqrt(w_DY100*w_DY100*SB1_DY100->GetBinContent(i) + w_DY70*w_DY70*SB1_DY70->GetBinContent(i) + w_DYM50_100*w_DYM50_100*SB1_DYM50_100->GetBinContent(i) + 
      //			  w_DYM50_70*w_DYM50_70*SB1_DYM50_70->GetBinContent(i) + 1.9*1.9*w_QCD250*w_QCD250*SB1_QCD250->GetBinContent(i) + 
      //			  1.9*1.9*w_QCD500*w_QCD500*SB1_QCD500->GetBinContent(i) + 1.9*1.9*w_QCD1000*w_QCD1000*SB1_QCD1000->GetBinContent(i) + 
      //			  w_TT*w_TT*SB1_TT->GetBinContent(i) + w_WW*w_WW*SB1_WW->GetBinContent(i) + w_WZ*w_WZ*SB1_WZ->GetBinContent(i) + 
      //			  w_ZZ*w_ZZ*SB1_ZZ->GetBinContent(i) + w_WJetsHT*w_WJetsHT*SB1_WJetsHT->GetBinContent(i)
      //			  );
      //float Bkg2 = (w_DY100*SB2_DY100->GetBinContent(i) + w_DY70*SB2_DY70->GetBinContent(i) + w_DYM50_100*SB2_DYM50_100->GetBinContent(i) + 
      //		   w_DYM50_70*SB2_DYM50_70->GetBinContent(i) + 1.9*w_QCD250*SB2_QCD250->GetBinContent(i) + 1.9*w_QCD500*SB2_QCD500->GetBinContent(i) + 
      //		   1.9*w_QCD1000*SB2_QCD1000->GetBinContent(i) + w_TT*SB2_TT->GetBinContent(i) + w_WW*SB2_WW->GetBinContent(i) + w_WZ*SB2_WZ->GetBinContent(i) + 
      //		   w_ZZ*SB2_ZZ->GetBinContent(i) + w_WJetsHT*SB2_WJetsHT->GetBinContent(i)
      //		   );
      //float BkgErr2 = sqrt(w_DY100*w_DY100*SB2_DY100->GetBinContent(i) + w_DY70*w_DY70*SB2_DY70->GetBinContent(i) + w_DYM50_100*w_DYM50_100*SB2_DYM50_100->GetBinContent(i) + 
      //			  w_DYM50_70*w_DYM50_70*SB2_DYM50_70->GetBinContent(i) + 1.9*1.9*w_QCD250*w_QCD250*SB2_QCD250->GetBinContent(i) + 
      //			  1.9*1.9*w_QCD500*w_QCD500*SB2_QCD500->GetBinContent(i) + 1.9*1.9*w_QCD1000*w_QCD1000*SB2_QCD1000->GetBinContent(i) + 
      //			  w_TT*w_TT*SB2_TT->GetBinContent(i) + w_WW*w_WW*SB2_WW->GetBinContent(i) + w_WZ*w_WZ*SB2_WZ->GetBinContent(i) + 
      //			  w_ZZ*w_ZZ*SB2_ZZ->GetBinContent(i) + w_WJetsHT*w_WJetsHT*SB2_WJetsHT->GetBinContent(i)
      //			  );
      float Bkg1 = SB1_data->GetBinContent(i);
      float Bkg2 = SB2_data->GetBinContent(i);
      float BkgErr1 = sqrt(SB1_data->GetBinContent(i));
      float BkgErr2 = sqrt(SB2_data->GetBinContent(i));
      SB2->SetBinContent(i,Bkg1);
      SB2->SetBinError(i,BkgErr1);
      SB1->SetBinContent(i,Bkg2);
      SB1->SetBinError(i,BkgErr2);
    }

    float SB1norm = 1/SB1->Integral();
    float SB2norm = 1/SB2->Integral();
    SB1->Scale(SB1norm);
    SB2->Scale(SB2norm);


    float chi2=0;
    for(int i=0; i<SB1->GetXaxis()->GetNbins(); i++){
      float Bkg1    = SB1->GetBinContent(i);
      float Bkg2    = SB2->GetBinContent(i);
      float BkgErr1 = SB1->GetBinError(i);
      float BkgErr2 = SB2->GetBinError(i);
      if((BkgErr1*BkgErr1+BkgErr2*BkgErr2)!=0) chi2 = chi2 + (Bkg1-Bkg2)*(Bkg1-Bkg2)/(BkgErr1*BkgErr1+BkgErr2*BkgErr2);
    }

    TH1D *RATIO = new TH1D("","",SB1->GetNbinsX(),SB1->GetXaxis()->GetXmin(),SB1->GetXaxis()->GetXmax());
    for(int m=1; m<SB1->GetNbinsX()+1; m++){ 
      float A=SB2->GetBinContent(m);
      float B=SB1->GetBinContent(m);
      float Aerr=SB1->GetBinError(m);
      float Berr=SB2->GetBinError(m);
      if(B!=0) {
	RATIO->SetBinContent(m,(A/B));
	RATIO->SetBinError(m,sqrt(Aerr*Aerr/(B*B) + Berr*Berr*A*A/(B*B*B*B)));
      }
      //cout<<SB1->GetBinContent(m)<<" "<<SB2->GetBinContent(m)<<endl;
    }
    RATIO->Draw("E");
    RATIO->SetMaximum(3);
    RATIO->SetMinimum(-1);
    RATIO->SetLineColor(1);
    RATIO->SetLineWidth(2);
    RATIO->GetXaxis()->SetTitleOffset(0.9);
    RATIO->GetXaxis()->SetLabelSize(0.08);
    RATIO->GetXaxis()->SetTitleSize(0.09);
    RATIO->GetYaxis()->SetTitleOffset(0.35);
    RATIO->GetYaxis()->SetTitle("A/B");
    RATIO->GetXaxis()->SetTitle(axis);
    RATIO->GetYaxis()->SetLabelSize(0.08);
    RATIO->GetYaxis()->SetTitleSize(0.09);
    RATIO->SetMarkerStyle(21);
    
    TLine* line = new TLine(SB1->GetXaxis()->GetXmin(),1,SB1->GetXaxis()->GetXmax(),1);
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

    SB2->Draw("E");
    SB1->Draw("E same");
    SB1->SetMinimum(0);
    SB2->SetMinimum(0);
    SB1->SetMarkerStyle(21);
    SB2->SetMarkerStyle(21);
    SB1->SetMarkerColor(1);
    SB2->SetMarkerColor(2);

    TString CHANNEL = channel;
    Double_t res[SB1->GetXaxis()->GetNbins()];
    cout<<CHANNEL<<endl;
    cout<<"Prob="<<TMath::Prob(chi2,SB1->GetNbinsX())<<",   chi2="<<chi2<<",   nbin="<<SB1->GetNbinsX()<<endl;
    //cout<<"Chi2 Test = "<<SB1->Chi2Test(SB2,"WW P",res)<<endl;
    cout<<"Kolmogorov Test = "<<SB1->KolmogorovTest(SB2)<<endl;
      
    TLegend *pl2 = new TLegend(0.60,0.75,0.97,0.98);
    //TLegend *pl2 = new TLegend(0.55,0.70,0.89,0.89);
    pl2->SetTextSize(0.04); 
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(SB1, "20 GeV < mjet < 70 GeV",  "LP");
    ple2 = pl2->AddEntry(SB2, "mjet > 140 GeV",  "LP"); 
    pl2->Draw();
      
    if(log) {
      hs->SetMinimum(0.05);
      c1->SetLogy();
    }
      
    char Btag        [2];  sprintf(Btag,        "%.0f", bTagCut);      TString btag         = Btag        ;
    char Drcut       [2];  sprintf(Drcut,       "%.0f", deltaRCut);    TString drcut        = Drcut       ;
    char Masssvfitcut[2];  sprintf(Masssvfitcut,"%.0f", MassSvfitCut); TString masssvfitcut = Masssvfitcut;
    char Massviscut  [2];  sprintf(Massviscut,  "%.0f", MassVisCut);   TString massviscut   = Massviscut  ;
    char Met         [2];  sprintf(Met,         "%.0f", METCut);       TString met          = Met         ;
    if(save)c1->SaveAs("shapeComparison_"+CHANNEL+"_MET"+met+"_btag"+btag+"_dR"+drcut+"_massSvfit"+masssvfitcut+"_massVis"+massviscut+".pdf");
  }
}
