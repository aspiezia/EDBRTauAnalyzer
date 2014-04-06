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

  TH1F *histo1 = new TH1F("","",10,0,2500);
  TH1F *histo2 = new TH1F("","",10,0,2500);
  TH1F *histo3 = new TH1F("","",10,0,2500);
  TH1F *ERR    = new TH1F("","",10,0,2500);

  char channel  [500]; sprintf(channel,   "MuoTau"); 
  char demo     [500]; sprintf(demo,     "demo/Tree");  
  char openTree [500]; sprintf(openTree, "%s%s",demo,channel);
  TString CHANNEL = channel;
  bool save = true;
  bool lepTauIso = true;
		     
  char *charge = "charge";
  TH1F *data1_enr = new TH1F("","",3,-1.5,1.5);
  TH1F *data2_enr = new TH1F("","",3,-1.5,1.5);
  TH1F *data3_enr = new TH1F("","",3,-1.5,1.5);
  TH1F *data4_enr = new TH1F("","",3,-1.5,1.5);

  TFile *file15 = TFile::Open("../RISULTATI/05_260114/RunA_QCD.root");  TTree *Tree15 = (TTree*)file15->Get(openTree); 
  TFile *file16 = TFile::Open("../RISULTATI/05_260114/RunB_QCD.root");  TTree *Tree16 = (TTree*)file16->Get(openTree); 
  TFile *file17 = TFile::Open("../RISULTATI/05_260114/RunC_QCD.root");  TTree *Tree17 = (TTree*)file17->Get(openTree); 
  TFile *file18 = TFile::Open("../RISULTATI/05_260114/RunD_QCD.root");  TTree *Tree18 = (TTree*)file18->Get(openTree);

  char input15 [50]; sprintf(input15, "%s>>h15(%i,%f,%f)",charge,3,-1.5,1.5);
  char input16 [50]; sprintf(input16, "%s>>h16(%i,%f,%f)",charge,3,-1.5,1.5);
  char input17 [50]; sprintf(input17, "%s>>h17(%i,%f,%f)",charge,3,-1.5,1.5);
  char input18 [50]; sprintf(input18, "%s>>h18(%i,%f,%f)",charge,3,-1.5,1.5);
  Tree15->Draw(input15,"trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)");     
  if(Tree15->Draw(input15,"trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)")) data1_enr = h15;
  Tree16->Draw(input16,"trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)");    
  if(Tree16->Draw(input16,"trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)")) data2_enr = h16;
  Tree17->Draw(input17,"trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)");    
  if(Tree17->Draw(input17,"trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)")) data3_enr = h17;
  Tree18->Draw(input18,"trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)");    
  if(Tree18->Draw(input18,"trigger*PUWeight*(met>50 && uncorrmet>30 && MassVis>20 && lepCorrPFIso>0.2 && tauIso<0.5 && nbtagsM==0)")) data4_enr = h18;
  data1_enr->Add(data2_enr);
  data1_enr->Add(data3_enr);
  data1_enr->Add(data4_enr);

  double a = data1_enr->GetBinContent(1);
  double b = data1_enr->GetBinContent(3);
  double NOSoverNSS     = a/b;
  double NOSoverNSS_err = (sqrt(b*b*a + a*a*b))/(b*b);

  cout<<NOSoverNSS<<" "<<NOSoverNSS_err<<endl;

  char *plot = "XMassSVFit";
  int bin=250; 
  float min=0; 
  float max=250; 
  for(int i=0; i<10; i++){
    TH1F *data1 = new TH1F("","",bin,min,max);
    TH1F *data2 = new TH1F("","",bin,min,max);
    TH1F *data3 = new TH1F("","",bin,min,max);
    TH1F *data4 = new TH1F("","",bin,min,max);
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
    
    char CUT1    [500]; 
    char CUT_QCD1[500]; 
    if(lepTauIso && CHANNEL=="MuoTau"){
      sprintf(CUT1,    "trigger*PUWeight*(sideband==0 && met>50 && uncorrmet>30 && MassVis>10 && lepCorrPFIso<0.2 && charge==+1");
      sprintf(CUT_QCD1,"trigger*PUWeight*(sideband==0 && met>50 && uncorrmet>30 && MassVis>10 && lepCorrPFIso<0.2 && charge==-1");
    }
    if(lepTauIso && CHANNEL=="EleTau"){
      sprintf(CUT1,    "trigger*PUWeight*(sideband==0 && met>50 && uncorrmet>30 && MassVis>10 && lepCorrPFIso<0.1 && charge==+1");
      sprintf(CUT_QCD1,"trigger*PUWeight*(sideband==0 && met>50 && uncorrmet>30 && MassVis>10 && lepCorrPFIso<0.1 && charge==-1");
    }
    if(!lepTauIso){
      sprintf(CUT1,    "trigger*PUWeight*(sideband==0 && met>50 && uncorrmet>30 && MassVis>10 && charge==+1");
      sprintf(CUT_QCD1,"trigger*PUWeight*(sideband==0 && met>50 && uncorrmet>30 && MassVis>10 && charge==-1");
    }
    char CUT2[500]; sprintf(CUT2,        " && XMassSVFit>%f && XMassSVFit<%f)",min,max);
    char CUT_QCD2[500]; sprintf(CUT_QCD2," && XMassSVFit>%f && XMassSVFit<%f)",min,max);
    char CUT [500]; sprintf(CUT,         "%s%s",CUT1,CUT2);    
    char CUT_QCD [500]; sprintf(CUT_QCD, "%s%s",CUT_QCD1,CUT_QCD2);

    TFile *file1  = TFile::Open("../RISULTATI/05_260114/RunA_SL.root");       TTree *Tree1  = (TTree*)file1->Get(openTree); 
    TFile *file2  = TFile::Open("../RISULTATI/05_260114/RunB_Final_SL.root"); TTree *Tree2  = (TTree*)file2->Get(openTree); 
    TFile *file3  = TFile::Open("../RISULTATI/05_260114/RunC_Final_SL.root"); TTree *Tree3  = (TTree*)file3->Get(openTree); 
    TFile *file4  = TFile::Open("../RISULTATI/05_260114/RunD_Final_SL.root"); TTree *Tree4  = (TTree*)file4->Get(openTree);
    TFile *file5  = TFile::Open("../RISULTATI/05_260114/DY100_SL.root");      TTree *Tree5  = (TTree*)file5->Get(openTree); 
    TFile *file6  = TFile::Open("../RISULTATI/05_260114/DY70_SL.root");       TTree *Tree6  = (TTree*)file6->Get(openTree); 
    TFile *file51 = TFile::Open("../RISULTATI/05_260114/DYM50_100_SL.root");  TTree *Tree51 = (TTree*)file51->Get(openTree); 
    TFile *file61 = TFile::Open("../RISULTATI/05_260114/DYM50_70_SL.root");   TTree *Tree61 = (TTree*)file61->Get(openTree); 
    TFile *file7  = TFile::Open("../RISULTATI/05_260114/QCD1000_SL.root");    TTree *Tree7  = (TTree*)file7->Get(openTree);
    TFile *file8  = TFile::Open("../RISULTATI/05_260114/QCD250_SL.root");     TTree *Tree8  = (TTree*)file8->Get(openTree); 
    TFile *file9  = TFile::Open("../RISULTATI/05_260114/QCD500_SL.root");     TTree *Tree9  = (TTree*)file9->Get(openTree); 
    TFile *file10 = TFile::Open("../RISULTATI/05_260114/TT_SL.root");         TTree *Tree10 = (TTree*)file10->Get(openTree);
    TFile *file11 = TFile::Open("../RISULTATI/05_260114/WJets180_SL.root");   TTree *Tree11 = (TTree*)file11->Get(openTree); 
    TFile *file12 = TFile::Open("../RISULTATI/05_260114/WW_SL.root");         TTree *Tree12 = (TTree*)file12->Get(openTree); 
    TFile *file13 = TFile::Open("../RISULTATI/05_260114/WZ_SL.root");         TTree *Tree13 = (TTree*)file13->Get(openTree); 
    TFile *file14 = TFile::Open("../RISULTATI/05_260114/ZZ_SL.root");         TTree *Tree14 = (TTree*)file14->Get(openTree);  

    char input1  [50]; sprintf(input1,  "%s>>h1(%i,%f,%f)",  plot,bin,min,max);
    char input2  [50]; sprintf(input2,  "%s>>h2(%i,%f,%f)",  plot,bin,min,max);
    char input3  [50]; sprintf(input3,  "%s>>h3(%i,%f,%f)",  plot,bin,min,max);
    char input4  [50]; sprintf(input4,  "%s>>h4(%i,%f,%f)",  plot,bin,min,max);
    char input5  [50]; sprintf(input5,  "%s>>h5(%i,%f,%f)",  plot,bin,min,max);
    char input6  [50]; sprintf(input6,  "%s>>h6(%i,%f,%f)",  plot,bin,min,max);
    char input51 [50]; sprintf(input51, "%s>>h51(%i,%f,%f)", plot,bin,min,max);
    char input61 [50]; sprintf(input61, "%s>>h61(%i,%f,%f)", plot,bin,min,max);
    char input7  [50]; sprintf(input7,  "%s>>h7(%i,%f,%f)",  plot,bin,min,max);
    char input8  [50]; sprintf(input8,  "%s>>h8(%i,%f,%f)",  plot,bin,min,max);
    char input9  [50]; sprintf(input9,  "%s>>h9(%i,%f,%f)",  plot,bin,min,max);
    char input10 [50]; sprintf(input10, "%s>>h10(%i,%f,%f)", plot,bin,min,max);
    char input11 [50]; sprintf(input11, "%s>>h11(%i,%f,%f)", plot,bin,min,max);
    char input12 [50]; sprintf(input12, "%s>>h12(%i,%f,%f)", plot,bin,min,max);
    char input13 [50]; sprintf(input13, "%s>>h13(%i,%f,%f)", plot,bin,min,max);
    char input14 [50]; sprintf(input14, "%s>>h14(%i,%f,%f)", plot,bin,min,max);

    Tree1->Draw(input1,CUT,"E"); 
    if(Tree1->Draw(input1,CUT,"E")) data1 = h1;
    Tree2->Draw(input2,CUT,"E"); 
    if(Tree2->Draw(input2,CUT,"E")) data2 = h2;
    Tree3->Draw(input3,CUT,"E"); 
    if(Tree3->Draw(input3,CUT,"E")) data3 = h3;
    Tree4->Draw(input4,CUT,"E"); 
    if(Tree4->Draw(input4,CUT,"E")) data4 = h4;
  
    Tree5->Draw(input5,CUT);     
    if(Tree5->Draw(input5,CUT)) DY100 = h5;
    Tree6->Draw(input6,CUT);     
    if(Tree6->Draw(input6,CUT)) DY70 = h6;
    Tree51->Draw(input51,CUT);     
    if(Tree51->Draw(input51,CUT)) DYM50_100 = h51;
    Tree61->Draw(input61,CUT);     
    if(Tree61->Draw(input61,CUT)) DYM50_70 = h61; 
    Tree7->Draw(input7,CUT_QCD);     
    if(Tree7->Draw(input7,CUT_QCD)) QCD1000 = h7;
    Tree8->Draw(input8,CUT_QCD);     
    if(Tree8->Draw(input8,CUT_QCD)) QCD250 = h8;
    Tree9->Draw(input9,CUT_QCD);     
    if(Tree9->Draw(input9,CUT_QCD)) QCD500 = h9; 
    Tree10->Draw(input10,CUT);    
    if(Tree10->Draw(input10,CUT)) TT = h10;
    Tree11->Draw(input11,CUT);    
    if(Tree11->Draw(input11,CUT)) WJets180 = h11;
    Tree12->Draw(input12,CUT);    
    if(Tree12->Draw(input12,CUT)) WW = h12;
    Tree13->Draw(input13,CUT);    
    if(Tree13->Draw(input13,CUT)) WZ = h13;
    Tree14->Draw(input14,CUT);    
    if(Tree14->Draw(input14,CUT)) ZZ = h14;
    
    data1->Add(data2);
    data1->Add(data3);
    data1->Add(data4);

    float w_DY100     = ( 32.900*19702./12511326.);
    float w_DY70      = ( 53.000*19702./11764538.);
    float w_DYM50_100 = (  4.220*19702./4146124.0);
    float w_DYM50_70  = ( 11.050*19702./5389313.0);
    float w_TT        = (225.197*19702./21675970.);
    float w_QCD1000   = (204.000*19702./13843863.);
    float w_QCD500    = (8426.00*19702./30599292.);
    float w_QCD250    = (276000.*19702./27062078.);
    float w_WW        = (57.1097*19702./10000431.);
    float w_WZ        = ( 33.210*19702./10000283.);
    float w_ZZ        = (  8.059*19702./9799908.0);
    float w_WJets180  = ( 23.500*19702./9739464.0);
    
    double N_SS = (data1->Integral() - (w_DY100*DY100->Integral()+w_DY70*DY70->Integral()+
					w_DYM50_100*DYM50_100->Integral()+w_DYM50_70*DYM50_70->Integral()+w_WZ*WZ->Integral()+
					w_WW*WW->Integral()+w_ZZ*ZZ->Integral()+
					w_WJets180*WJets180->Integral()+w_TT*TT->Integral())) * NOSoverNSS;
    double N_SS_err = sqrt(data1->Integral()+w_DY100*w_DY100*DY100->Integral()+w_DY70*w_DY70*DY70->Integral()+
			   w_DYM50_100*w_DYM50_100*DYM50_100->Integral()+w_DYM50_70*w_DYM50_70*DYM50_70->Integral()+w_WZ*w_WZ*WZ->Integral()+
			   w_WW*w_WW*WW->Integral()+w_ZZ*w_ZZ*ZZ->Integral()+w_WJets180*w_WJets180*WJets180->Integral()+w_TT*w_TT*TT->Integral() + 
			   NOSoverNSS_err*NOSoverNSS_err); 

    double N_QCD = 1.9*(w_QCD1000*QCD1000->Integral()+w_QCD250*QCD250->Integral()+w_QCD500*QCD500->Integral());
    double N_QCD_err = sqrt(1.9*1.9*(w_QCD1000*w_QCD1000*QCD1000->Integral()+w_QCD250*w_QCD250*QCD250->Integral()+w_QCD500*w_QCD500*QCD500->Integral()));
     
    histo1->SetBinContent(i+1,N_SS);
    histo1->SetBinError(i+1,N_SS_err);
    histo2->SetBinContent(i+1,N_QCD);
    histo2->SetBinError(i+1,N_QCD_err);
    ERR->SetBinContent(i+1,N_QCD);
    ERR->SetBinError(i+1,N_QCD_err);

    min=min+bin;
    max=max+bin;
  }


  histo1->SetLineWidth(2); 
  histo1->SetLineColor(kRed);
  histo1->SetMarkerColor(kRed);
  histo1->SetFillColor(kRed); 
  histo1->SetMarkerStyle(20); 
  histo1->SetMarkerSize(1.3);
  histo2->SetFillColor(kRed-2);
  histo2->SetLineColor(kRed-2);
  histo2->SetMarkerColor(kRed-2); 

  histo2->Draw("histo");
  histo2->GetYaxis()->SetTitleSize(0.045);
  histo2->GetXaxis()->SetTitleSize(0.045);
  histo2->GetYaxis()->SetLabelSize(0.045);
  histo2->GetXaxis()->SetLabelSize(0.045);
  histo2->GetYaxis()->SetTitleOffset(1.0); 
  histo2->GetYaxis()->SetTitle("Events");
  histo2->GetXaxis()->SetTitle("M(Z,H) [GeV] (SvFit)"); 
  histo2->GetXaxis()->SetNdivisions(905);
  histo2->SetMinimum(0);
  if(CHANNEL=="MuoTau" && !lepTauIso) histo2->SetMaximum(200);
  if(CHANNEL=="EleTau" && !lepTauIso) histo2->SetMaximum(50);
  if(CHANNEL=="MuoTau" &&  lepTauIso) histo2->SetMaximum(30);
  if(CHANNEL=="EleTau" &&  lepTauIso) histo2->SetMaximum(20);

  histo1->Draw("Esame");

  ERR->SetFillStyle(3005);
  ERR->SetFillColor(12);
  ERR->SetLineColor(12);
  ERR->Draw("E2same");

  

  TLegend *pl2 = new TLegend(0.25,0.77,0.89,0.89);
  pl2->SetTextSize(0.030); 
  pl2->SetFillColor(0);
  TLegendEntry *ple2 = pl2->AddEntry(histo1, "QCD prediction from data-driven method",  "LP");
  ple2 = pl2->AddEntry(histo2, "QCD prediction from MC",  "F");
  pl2->Draw();

  if(save && !lepTauIso){
    c1->SaveAs("QCD_"+CHANNEL+"_prediction.png");
    //c1->SaveAs("QCD_"+CHANNEL+"_prediction.root");
    //c1->SaveAs("QCD_"+CHANNEL+"_prediction.C");
  }
  if(save && lepTauIso){
    c1->SaveAs("QCD_"+CHANNEL+"_prediction_tauIso.png");
    //c1->SaveAs("QCD_"+CHANNEL+"_prediction_tauIso.root");
    //c1->SaveAs("QCD_"+CHANNEL+"_prediction_tauIso.C");
  }
}
