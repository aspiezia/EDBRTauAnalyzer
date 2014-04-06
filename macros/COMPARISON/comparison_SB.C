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
  
  char channel   [500]; sprintf(channel,  "MuoMuo"); 
  char demo1     [500]; sprintf(demo1,     "demo/TreeSB1"); 
  char demo2     [500]; sprintf(demo2,     "demo/TreeSB2"); 
  char openTree1 [500]; sprintf(openTree1, "%s%s",demo1,channel); 
  char openTree2 [500]; sprintf(openTree2, "%s%s",demo2,channel); 
  TString CHANNEL = channel;
  bool save=true; 
  bool log=true;
  bool lepTauIso = false;

  vector<string> PLOT;                vector<int> BIN;   vector<float> MIN;  vector<float> MAX;   vector<float> MAXY; vector<TString> AXIS;
  if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo" || CHANNEL=="EleMuo"){
    PLOT.push_back("MassVis");        BIN.push_back(20); MIN.push_back(0);   MAX.push_back(240);  MAXY.push_back(300); AXIS.push_back("M(lep1,lep2) [GeV] (visible)");
    PLOT.push_back("NVertices");      BIN.push_back(20); MIN.push_back(-0.5);MAX.push_back(39.5); MAXY.push_back(300); AXIS.push_back("N(vertices)");
    PLOT.push_back("XMassSVFit");     BIN.push_back(50); MIN.push_back(0);   MAX.push_back(2500); MAXY.push_back(300); AXIS.push_back("M(Z,H) [GeV] (SvFit)");
    PLOT.push_back("jetMass");        BIN.push_back(20); MIN.push_back(20);  MAX.push_back(220);  MAXY.push_back(300); AXIS.push_back("jet mass [GeV]");
    PLOT.push_back("jetPt");          BIN.push_back(30); MIN.push_back(400); MAX.push_back(1000); MAXY.push_back(300); AXIS.push_back("jet pt [GeV]");
    PLOT.push_back("jetSubjettiness");BIN.push_back(20); MIN.push_back(0);   MAX.push_back(1);    MAXY.push_back(300); AXIS.push_back("jet #tau_{21}");
    PLOT.push_back("met");            BIN.push_back(25); MIN.push_back(50);  MAX.push_back(500);  MAXY.push_back(300); AXIS.push_back("met [GeV]");
  } else {
    PLOT.push_back("MassVis");        BIN.push_back(25); MIN.push_back(0);   MAX.push_back(250);  MAXY.push_back(40); AXIS.push_back("M(tau,lep) [GeV] (visible)");
    PLOT.push_back("NVertices");      BIN.push_back(41); MIN.push_back(-0.5);MAX.push_back(40.5); MAXY.push_back(10); AXIS.push_back("N(vertices)");
    PLOT.push_back("XMassSVFit");     BIN.push_back(50); MIN.push_back(0);   MAX.push_back(2500); MAXY.push_back(25); AXIS.push_back("M(Z,H) [GeV] (SvFit)");
    PLOT.push_back("jetMass");        BIN.push_back(25); MIN.push_back(20);  MAX.push_back(70);   MAXY.push_back(25); AXIS.push_back("jet mass [GeV]");
    PLOT.push_back("jetPt");          BIN.push_back(30); MIN.push_back(400); MAX.push_back(1000); MAXY.push_back(300); AXIS.push_back("jet pt [GeV]");
    PLOT.push_back("jetSubjettiness");BIN.push_back(20); MIN.push_back(0);   MAX.push_back(1);    MAXY.push_back(20); AXIS.push_back("jet #tau_{21}");
    PLOT.push_back("met");            BIN.push_back(25); MIN.push_back(50);  MAX.push_back(500);  MAXY.push_back(30); AXIS.push_back("met [GeV]");
  }

  for(int i=0; i<PLOT.size(); i++){
  //for(int i=3; i<4; i++){
    char *plot = PLOT[i].c_str();
    TString name = PLOT[i];
    int bin=BIN[i]; 
    float min=MIN[i]; 
    float max=MAX[i]; 
    float maxy=MAXY[i];
    TString axis = AXIS[i];

    //int N = 17; Double_t xbins[N] = {0,20,30,40,50,60,70,80,90,100,110,120,140,160,180,200,250}; //VIS MASS
    //int N = 15;  Double_t xbins[N] = {0,4,6,8,10,12,14,16,18,20,22,24,26,28,34}; //VERTICES
    //int N = 16;  Double_t xbins[N] = {0,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1800,2000,2500}; //X MASS

    char CUT [500]; 
    if(!lepTauIso){
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo") sprintf(CUT, "trigger*PUWeight*(met>30 && uncorrmet>30 && MassVis>10)");
      else if(CHANNEL=="EleMuo") sprintf(CUT, "trigger*PUWeight*(met>30 && uncorrmet>30)");
      else sprintf(CUT, "trigger*PUWeight*(met>30 && uncorrmet>30 && charge==-1)");
    } else {
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo") sprintf(CUT, "trigger*PUWeight*(met>30 && uncorrmet>30 && MassVis>10)");
      else if(CHANNEL=="EleMuo") sprintf(CUT, "trigger*PUWeight*(met>30 && uncorrmet>30)");
      else if(CHANNEL=="EleTau")  sprintf(CUT, "trigger*PUWeight*(met>30 && uncorrmet>30 && lepCorrPFIso<0.1 && charge==-1)");
      else   sprintf(CUT, "trigger*PUWeight*(met>30 && uncorrmet>30 && lepCorrPFIso<0.2 && charge==-1)");
    }

    cout<<endl;
    cout<<"Plotting plot "<<i+1<<"/"<<PLOT.size()<<", i.e. "<<name<<endl;

    TH1F *data1    = new TH1F("","",bin,min,max);
    TH1F *data2    = new TH1F("","",bin,min,max);
    TH1F *data3    = new TH1F("","",bin,min,max);
    TH1F *data4    = new TH1F("","",bin,min,max);
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
  
    TFile *file1 =TFile::Open("../RISULTATI/SB2_050214/RunA_FL.root");       TTree *Tree1  = (TTree*)file1->Get(openTree1);  TTree *Tree1_1 = (TTree*)file1->Get(openTree2);
    TFile *file2 =TFile::Open("../RISULTATI/SB2_050214/RunB_Final_FL.root"); TTree *Tree2  = (TTree*)file2->Get(openTree1);  TTree *Tree2_1 = (TTree*)file2->Get(openTree2);
    TFile *file3 =TFile::Open("../RISULTATI/SB2_050214/RunC_Final_FL.root"); TTree *Tree3  = (TTree*)file3->Get(openTree1);  TTree *Tree3_1 = (TTree*)file3->Get(openTree2);
    TFile *file4 =TFile::Open("../RISULTATI/SB2_050214/RunD_Final_FL.root"); TTree *Tree4  = (TTree*)file4->Get(openTree1);  TTree *Tree4_1 = (TTree*)file4->Get(openTree2);
    TFile *file5 =TFile::Open("../RISULTATI/SB2_050214/DY100_FL.root");      TTree *Tree5  = (TTree*)file5->Get(openTree1);  TTree *Tree5_1 = (TTree*)file5->Get(openTree2);
    TFile *file6 =TFile::Open("../RISULTATI/SB2_050214/DY70_FL.root");       TTree *Tree6  = (TTree*)file6->Get(openTree1);  TTree *Tree6_1 = (TTree*)file6->Get(openTree2);
    TFile *file51=TFile::Open("../RISULTATI/SB2_050214/DYM50_100_FL.root");  TTree *Tree51 = (TTree*)file51->Get(openTree1); TTree *Tree51_1= (TTree*)file51->Get(openTree2);
    TFile *file61=TFile::Open("../RISULTATI/SB2_050214/DYM50_70_FL.root");   TTree *Tree61 = (TTree*)file61->Get(openTree1); TTree *Tree61_1= (TTree*)file61->Get(openTree2);
    TFile *file7 =TFile::Open("../RISULTATI/SB2_050214/QCD1000_FL.root");    TTree *Tree7  = (TTree*)file7->Get(openTree1);  TTree *Tree7_1 = (TTree*)file7->Get(openTree2);
    TFile *file8 =TFile::Open("../RISULTATI/SB2_050214/QCD250_FL.root");     TTree *Tree8  = (TTree*)file8->Get(openTree1);  TTree *Tree8_1 = (TTree*)file8->Get(openTree2);
    TFile *file9 =TFile::Open("../RISULTATI/SB2_050214/QCD500_FL.root");     TTree *Tree9  = (TTree*)file9->Get(openTree1);  TTree *Tree9_1 = (TTree*)file9->Get(openTree2);
    TFile *file10=TFile::Open("../RISULTATI/SB2_050214/TT_FL.root");         TTree *Tree10 = (TTree*)file10->Get(openTree1); TTree *Tree10_1= (TTree*)file10->Get(openTree2);
    TFile *file11=TFile::Open("../RISULTATI/SB2_050214/WJets180_FL.root");   TTree *Tree11 = (TTree*)file11->Get(openTree1); TTree *Tree11_1= (TTree*)file11->Get(openTree2);
    TFile *file12=TFile::Open("../RISULTATI/SB2_050214/WW_FL.root");         TTree *Tree12 = (TTree*)file12->Get(openTree1); TTree *Tree12_1= (TTree*)file12->Get(openTree2);
    TFile *file13=TFile::Open("../RISULTATI/SB2_050214/WZ_FL.root");         TTree *Tree13 = (TTree*)file13->Get(openTree1); TTree *Tree13_1= (TTree*)file13->Get(openTree2);
    TFile *file14=TFile::Open("../RISULTATI/SB2_050214/ZZ_FL.root");         TTree *Tree14 = (TTree*)file14->Get(openTree1); TTree *Tree14_1= (TTree*)file14->Get(openTree2);
    
    char input1  [50]; sprintf(input1,  "%s>>h1(%i,%f,%f)", plot,bin,min,max);
    char input2  [50]; sprintf(input2,  "%s>>h2(%i,%f,%f)", plot,bin,min,max);
    char input3  [50]; sprintf(input3,  "%s>>h3(%i,%f,%f)", plot,bin,min,max);
    char input4  [50]; sprintf(input4,  "%s>>h4(%i,%f,%f)", plot,bin,min,max);
    char input5  [50]; sprintf(input5,  "%s>>h5(%i,%f,%f)", plot,bin,min,max);
    char input6  [50]; sprintf(input6,  "%s>>h6(%i,%f,%f)", plot,bin,min,max);
    char input51 [50]; sprintf(input51, "%s>>h51(%i,%f,%f)",plot,bin,min,max);
    char input61 [50]; sprintf(input61, "%s>>h61(%i,%f,%f)",plot,bin,min,max);
    char input7  [50]; sprintf(input7,  "%s>>h7(%i,%f,%f)", plot,bin,min,max);
    char input8  [50]; sprintf(input8,  "%s>>h8(%i,%f,%f)", plot,bin,min,max);
    char input9  [50]; sprintf(input9,  "%s>>h9(%i,%f,%f)", plot,bin,min,max);
    char input10 [50]; sprintf(input10, "%s>>h10(%i,%f,%f)",plot,bin,min,max);
    char input11 [50]; sprintf(input11, "%s>>h11(%i,%f,%f)",plot,bin,min,max);
    char input12 [50]; sprintf(input12, "%s>>h12(%i,%f,%f)",plot,bin,min,max);
    char input13 [50]; sprintf(input13, "%s>>h13(%i,%f,%f)",plot,bin,min,max);
    char input14 [50]; sprintf(input14, "%s>>h14(%i,%f,%f)",plot,bin,min,max);

    Tree1->Draw(input1,CUT,"E"); if(Tree1->Draw(input1,CUT,"E")) {data1 = h1;      }
    Tree2->Draw(input2,CUT,"E"); if(Tree2->Draw(input2,CUT,"E")) {data2 = h2;	   }
    Tree3->Draw(input3,CUT,"E"); if(Tree3->Draw(input3,CUT,"E")) {data3 = h3;	   }
    Tree4->Draw(input4,CUT,"E"); if(Tree4->Draw(input4,CUT,"E")) {data4 = h4;	   }
    Tree5->Draw(input5,CUT);     if(Tree5->Draw(input5,CUT))     {DY100 = h5;	   }
    Tree6->Draw(input6,CUT);     if(Tree6->Draw(input6,CUT))     {DY70 = h6; 	   }
    Tree51->Draw(input51,CUT);   if(Tree51->Draw(input51,CUT))   {DYM50_100 = h51; }
    Tree61->Draw(input61,CUT);   if(Tree61->Draw(input61,CUT))   {DYM50_70 = h61;  }
    Tree7->Draw(input7,CUT);     if(Tree7->Draw(input7,CUT))     {QCD1000 = h7;	   }
    Tree8->Draw(input8,CUT);     if(Tree8->Draw(input8,CUT))     {QCD250 = h8;	   }
    Tree9->Draw(input9,CUT);     if(Tree9->Draw(input9,CUT))     {QCD500 = h9;	   }
    Tree10->Draw(input10,CUT);   if(Tree10->Draw(input10,CUT))   {TT = h10; 	   }
    Tree11->Draw(input11,CUT);   if(Tree11->Draw(input11,CUT))   {WJets180 = h11;  }
    Tree12->Draw(input12,CUT);   if(Tree12->Draw(input12,CUT))   {WW = h12;	   }
    Tree13->Draw(input13,CUT);   if(Tree13->Draw(input13,CUT))   {WZ = h13;	   }
    Tree14->Draw(input14,CUT);   if(Tree14->Draw(input14,CUT))   {ZZ = h14;        }
    
    sprintf(input1,  "%s>>h1_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input2,  "%s>>h2_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input3,  "%s>>h3_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input4,  "%s>>h4_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input5,  "%s>>h5_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input6,  "%s>>h6_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input51, "%s>>h51_1(%i,%f,%f)",plot,bin,min,max);
    sprintf(input61, "%s>>h61_1(%i,%f,%f)",plot,bin,min,max);
    sprintf(input7,  "%s>>h7_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input8,  "%s>>h8_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input9,  "%s>>h9_1(%i,%f,%f)", plot,bin,min,max);
    sprintf(input10, "%s>>h10_1(%i,%f,%f)",plot,bin,min,max);
    sprintf(input11, "%s>>h11_1(%i,%f,%f)",plot,bin,min,max);
    sprintf(input12, "%s>>h12_1(%i,%f,%f)",plot,bin,min,max);
    sprintf(input13, "%s>>h13_1(%i,%f,%f)",plot,bin,min,max);
    sprintf(input14, "%s>>h14_1(%i,%f,%f)",plot,bin,min,max);
    Tree1_1->Draw(input1,CUT,"E"); if(Tree1_1->Draw(input1,CUT,"E")) {data1->Add(h1_1);      }
    Tree2_1->Draw(input2,CUT,"E"); if(Tree2_1->Draw(input2,CUT,"E")) {data2->Add(h2_1);	     }
    Tree3_1->Draw(input3,CUT,"E"); if(Tree3_1->Draw(input3,CUT,"E")) {data3->Add(h3_1);	     }
    Tree4_1->Draw(input4,CUT,"E"); if(Tree4_1->Draw(input4,CUT,"E")) {data4->Add(h4_1);	     }
    Tree5_1->Draw(input5,CUT);     if(Tree5_1->Draw(input5,CUT))     {DY100->Add(h5_1);	     }
    Tree6_1->Draw(input6,CUT);     if(Tree6_1->Draw(input6,CUT))     {DY70->Add(h6_1); 	     }
    Tree51_1->Draw(input51,CUT);   if(Tree51_1->Draw(input51,CUT))   {DYM50_100->Add(h51_1); }
    Tree61_1->Draw(input61,CUT);   if(Tree61_1->Draw(input61,CUT))   {DYM50_70->Add(h61_1);  }
    Tree7_1->Draw(input7,CUT);     if(Tree7_1->Draw(input7,CUT))     {QCD1000->Add(h7_1);    }
    Tree8_1->Draw(input8,CUT);     if(Tree8_1->Draw(input8,CUT))     {QCD250->Add(h8_1);     }
    Tree9_1->Draw(input9,CUT);     if(Tree9_1->Draw(input9,CUT))     {QCD500->Add(h9_1);     }
    Tree10_1->Draw(input10,CUT);   if(Tree10_1->Draw(input10,CUT))   {TT->Add(h10_1); 	     }
    Tree11_1->Draw(input11,CUT);   if(Tree11_1->Draw(input11,CUT))   {WJets180->Add(h11_1);  }
    Tree12_1->Draw(input12,CUT);   if(Tree12_1->Draw(input12,CUT))   {WW->Add(h12_1);	     }
    Tree13_1->Draw(input13,CUT);   if(Tree13_1->Draw(input13,CUT))   {WZ->Add(h13_1);        }
    Tree14_1->Draw(input14,CUT);   if(Tree14_1->Draw(input14,CUT))   {ZZ->Add(h14_1);        }

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

    TH1F *ERR = new TH1F("","",data1->GetNbinsX(),data1->GetXaxis()->GetXmin(),data1->GetXaxis()->GetXmax());
    for(int m=1; m<ERR->GetNbinsX()+1; m++){
      ERR->SetBinError(m,sqrt(
			      w_DY100     *  w_DY100     *DY100->GetBinContent(m)+
			      w_DY70      *  w_DY70      *DY70->GetBinContent(m)+
			      w_DYM50_100 *  w_DYM50_100 *DYM50_100->GetBinContent(m)+
			      w_DYM50_70  *  w_DYM50_70  *DYM50_70->GetBinContent(m)+
			      w_QCD1000   *  w_QCD1000   *QCD1000->GetBinContent(m)+
			      w_QCD500    *  w_QCD500    *QCD500->GetBinContent(m)+
			      w_QCD250    *  w_QCD250    *QCD250->GetBinContent(m)+
			      w_TT        *  w_TT        *TT->GetBinContent(m)+
			      w_WW        *  w_WW        *WW->GetBinContent(m)+
			      w_WZ        *  w_WZ        *WZ->GetBinContent(m)+
			      w_ZZ        *  w_ZZ        *ZZ->GetBinContent(m)+
			      w_WJets180  *  w_WJets180  *WJets180->GetBinContent(m)
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

    data1->Add(data2);
    data1->Add(data3);
    data1->Add(data4);
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

    ERR->Sumw2();
    ERR->Add(DY100);
    ERR->Add(WW);
    ERR->Add(TT);
    ERR->Add(WJets180);
    ERR->Add(QCD1000);

    //hnew1 = DY100->Rebin(N-1,"hnew1",xbins);
    //hnew2 = WW->Rebin(N-1,"hnew2",xbins);
    //hnew3 = TT->Rebin(N-1,"hnew3",xbins);
    //hnew4 = WJets180->Rebin(N-1,"hnew4",xbins);
    //hnew5 = QCD1000->Rebin(N-1,"hnew5",xbins);
    //hnew6 = data1->Rebin(N-1,"hnew6",xbins);
    //hnew7 = ERR->Rebin(N-1,"hnew7",xbins);

    TH1D *RATIO = new TH1D("","",ERR->GetNbinsX(),ERR->GetXaxis()->GetXmin(),ERR->GetXaxis()->GetXmax());
    for(int m=1; m<ERR->GetNbinsX()+1; m++){ 
      if(ERR->GetBinContent(m)!=0 && data1->GetBinContent(m)!=0) {
    	RATIO->SetBinContent(m,data1->GetBinContent(m)/ERR->GetBinContent(m));
    	RATIO->SetBinError(m,sqrt(ERR->GetBinContent(m)*ERR->GetBinContent(m)*data1->GetBinError(m)*data1->GetBinError(m)
    				  +data1->GetBinContent(m)*data1->GetBinContent(m)*ERR->GetBinError(m)*ERR->GetBinError(m))/(ERR->GetBinContent(m)*ERR->GetBinContent(m)));
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
    data1->SetLineWidth(2); 
    data1->SetLineColor(1);
    data1->SetMarkerStyle(20); 
    data1->SetMarkerSize(1.3); 
    data1->Draw("Esame");
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
    ple2 = pl2->AddEntry(data1, "data",  "LP");
    pl2->Draw();

    if(log) {
      hs->SetMinimum(0.5);
      c1_2->SetLogy();
    }
    if(save && !lepTauIso) c1->SaveAs(CHANNEL+"_SB_"+name+".png");
    if(save && lepTauIso)  c1->SaveAs(CHANNEL+"_SB_"+name+"_tauIso.png");
  }

}
