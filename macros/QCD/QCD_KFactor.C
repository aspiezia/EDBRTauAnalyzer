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

  vector<string> PLOT;               vector<int> BIN;    vector<float> MIN;   vector<float> MAX;    vector<TString> AXIS;/*
  PLOT.push_back("NVertices");       BIN.push_back(41);  MIN.push_back(-0.5); MAX.push_back(40.5);  AXIS.push_back("N(vertices)");
  PLOT.push_back("jetEta");          BIN.push_back(25);  MIN.push_back(-2.5); MAX.push_back(2.5);   AXIS.push_back("jet #eta");
  PLOT.push_back("jetMass");         BIN.push_back(25);  MIN.push_back(20);   MAX.push_back(70);    AXIS.push_back("jet mass [GeV]");
  PLOT.push_back("jetPt");           BIN.push_back(20);  MIN.push_back(400);  MAX.push_back(800);   AXIS.push_back("jet pt [GeV]");
  PLOT.push_back("lepEta");          BIN.push_back(10);  MIN.push_back(-2.5); MAX.push_back(2.5);   AXIS.push_back("lep #eta");*
  PLOT.push_back("lepPt");           BIN.push_back(25);  MIN.push_back(0);    MAX.push_back(500);   AXIS.push_back("lep pt [GeV]");
  PLOT.push_back("met");             BIN.push_back(20);  MIN.push_back(50);   MAX.push_back(250);   AXIS.push_back("met [GeV]");*/
  PLOT.push_back("metPhi");          BIN.push_back(14);  MIN.push_back(-3.5); MAX.push_back(3.5);   AXIS.push_back("met #phi");/*
  PLOT.push_back("tauEta");          BIN.push_back(10);  MIN.push_back(-2.5); MAX.push_back(2.5);   AXIS.push_back("tau #eta");
  PLOT.push_back("tauPt");           BIN.push_back(20);  MIN.push_back(20);   MAX.push_back(420);   AXIS.push_back("tau pt [GeV]");*/

  float medium=0;
  float error=0;
  for(int i=0; i<PLOT.size(); i++){
    //for(int i=33; i<34; i++){
    char *plot = PLOT[i].c_str();
    TString name = PLOT[i];
    int bin=BIN[i]; 
    float min=MIN[i]; 
    float max=MAX[i]; 
    TString axis = AXIS[i];
    bool save=false; 
    bool log=true;

    char CUT [500]; sprintf(CUT, "weight*trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)");

    char channel  [500]; sprintf(channel,  "EleTau"); 
    char demo     [500]; sprintf(demo,     "demo/Tree"); 
    char openTree [500]; sprintf(openTree, "%s%s",demo,channel); 

    cout<<endl;cout<<endl;cout<<endl;
    cout<<"Fitting plot "<<i+1<<"/"<<PLOT.size()<<", i.e. "<<name<<endl;

    TH1F *data1 = new TH1F("","",bin,min,max);
    TH1F *data2 = new TH1F("","",bin,min,max);
    TH1F *data3 = new TH1F("","",bin,min,max);
    TH1F *data4 = new TH1F("","",bin,min,max);
    TH1F *DY100 = new TH1F("","",bin,min,max);
    TH1F *DY70 = new TH1F("","",bin,min,max);											     
    TH1F *QCD1000 = new TH1F("","",bin,min,max);
    TH1F *QCD250 = new TH1F("","",bin,min,max);
    TH1F *QCD500 = new TH1F("","",bin,min,max);
    TH1F *TT = new TH1F("","",bin,min,max);
    TH1F *WJets180 = new TH1F("","",bin,min,max);
    TH1F *WW = new TH1F("","",bin,min,max);
    TH1F *WZ = new TH1F("","",bin,min,max);
    TH1F *ZZ = new TH1F("","",bin,min,max);

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
    data1->Sumw2();
    data2->Sumw2();
    data3->Sumw2();
    data4->Sumw2();

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
    //QCD1000->Scale(2);

    DY100->SetFillColor(kBlue-9);
    QCD1000->SetFillColor(kRed-2);
    WW->SetFillColor(kGreen+1);
    TT->SetFillColor(kBlue);
    WJets180->SetFillColor(kMagenta-5);
 
    for(int j=1; j<data1->GetNbinsX()+1; j++){
      data1->SetBinContent(j,data1->GetBinContent(j)-(DY100->GetBinContent(j)+TT->GetBinContent(j)+WW->GetBinContent(j)+WJets180->GetBinContent(j)));
    }
    
    TH1D *RATIO = new TH1D("","",QCD1000->GetNbinsX(),QCD1000->GetXaxis()->GetXmin(),QCD1000->GetXaxis()->GetXmax());
    for(int m=1; m<QCD1000->GetNbinsX()+1; m++){ 
      if(QCD1000->GetBinContent(m)!=0 && data1->GetBinContent(m)!=0) {
	RATIO->SetBinContent(m,data1->GetBinContent(m)/QCD1000->GetBinContent(m));
	RATIO->SetBinError(m,sqrt(QCD1000->GetBinContent(m)*QCD1000->GetBinContent(m)*data1->GetBinError(m)*data1->GetBinError(m)
				  +data1->GetBinContent(m)*data1->GetBinContent(m)*QCD1000->GetBinError(m)*QCD1000->GetBinError(m))/(QCD1000->GetBinContent(m)*QCD1000->GetBinContent(m)));
      }
    }
    
    RATIO->Draw("E");
    RATIO->SetMaximum(4);
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

    TF1 *func = new TF1("fit","pol0",QCD1000->GetXaxis()->GetXmin(),QCD1000->GetXaxis()->GetXmax());
    func->SetParameters(500,RATIO->GetMean(),RATIO->GetRMS());
    func->SetParNames("Constant");
    func->SetLineColor(kBlue);
    func->SetLineWidth(3);
    RATIO->Fit("fit");
    double * p = func->GetParameters();
    double * epar = func->GetParErrors();
    medium=medium+p[0];
    error=error+epar[0]*epar[0];
    //cout<<p[0]<<" +/- "<<epar[0]<<endl;
	
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

    TLatex latexLabel3;
    latexLabel3.SetTextSize(0.14);
    latexLabel3.SetNDC();
    latexLabel3.DrawLatex(0.64, 0.29, TString("#color[4]{SF = }")+TString::Format("#color[4]{%.2f +/- %.2f}",p[0],epar[0]));
    latexLabel3.SetTextColor(kBlue);

    c1->cd();

	
    TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.98,0.96);
    c1_2->Draw();
    c1_2->cd();
    c1_2->SetTopMargin(0.01);
    c1_2->SetBottomMargin(0.1);
    c1_2->SetRightMargin(0.02);
    c1_2->SetLeftMargin(0.07);
    QCD1000->SetTitle();
  
    QCD1000->Draw("histo");
    QCD1000->GetYaxis()->SetTitleOffset(0.80);
    QCD1000->GetYaxis()->SetTitleSize(0.045);
    QCD1000->GetXaxis()->SetTitleSize(0.045);
    QCD1000->GetYaxis()->SetLabelSize(0.045);
    QCD1000->GetXaxis()->SetLabelSize(0.045);
    QCD1000->GetYaxis()->SetTitleOffset(0.8); 
    QCD1000->GetYaxis()->SetTitle(TString("Events / ")+TString::Format("%.2f",(max-min)/bin));
    QCD1000->GetXaxis()->SetTitle(axis); 
    QCD1000->SetMinimum(0);
    QCD1000->SetMaximum(10000);
    data1->SetLineWidth(2); 
    data1->SetLineColor(1);
    data1->SetMarkerStyle(20); 
    data1->SetMarkerSize(1.3); 
    data1->Draw("Esame");
    //TH1F *ERR = QCD1000;
    //ERR->SetFillStyle(3005);
    //ERR->SetFillColor(12);
    //ERR->SetLineColor(12);
    //ERR->Draw("E2same");
  
    TLegend *pl2 = new TLegend(0.79,0.80,0.97,0.98);
    pl2->SetTextSize(0.035); 
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(QCD1000, "QCD",  "F");
    ple2 = pl2->AddEntry(data1, "data",  "LP");
    pl2->Draw();

    if(log){
      QCD1000->SetMinimum(0.05);
      c1_2->SetLogy();
    }

    TString CHANNEL = channel;
    if(save) c1->SaveAs("KfactorForQCD_"+CHANNEL+"_"+name+".png");
  }

  cout<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"The Final Scale Factor is "<<medium/PLOT.size()<<" +/- "<<sqrt(error/PLOT.size())<<endl;
}
