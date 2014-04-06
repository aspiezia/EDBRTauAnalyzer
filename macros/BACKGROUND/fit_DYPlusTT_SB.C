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
  
  float XMassWidth=200;
  float XMassBin=15;
  float XMassMin=0;
  float XMassMax=XMassWidth;

  TH1F *histo_PRE = new TH1F("","",XMassBin,0,XMassWidth*XMassBin);
  TH1F *histo_DrY = new TH1F("","",XMassBin,0,XMassWidth*XMassBin);
  TH1F *histo_TTb = new TH1F("","",XMassBin,0,XMassWidth*XMassBin);
  TH1F *histo_OBS = new TH1F("","",XMassBin,0,XMassWidth*XMassBin);
  TH1F *histo_QCD = new TH1F("","",XMassBin,0,XMassWidth*XMassBin);
  TH1F *histo_WJe = new TH1F("","",XMassBin,0,XMassWidth*XMassBin);
  TH1F *histo_VVp = new TH1F("","",XMassBin,0,XMassWidth*XMassBin);
  TH1F *ERR       = new TH1F("","",XMassBin,0,XMassWidth*XMassBin);

  char *plot = "met";
  int bin=1; 
  float min=0; 
  float max=50000; 
  bool save=true;
  bool lepTauIso = true;
  int MET=100;
  
  char channel    [500]; sprintf(channel,    "MuoTau"); 
  char demoSB     [500]; sprintf(demoSB,     "demo/TreeSB1"); 
  char openTreeSB [500]; sprintf(openTreeSB, "%s%s",demoSB,channel);
  char demo       [500]; sprintf(demo,       "demo/TreeSB2"); 
  char openTree   [500]; sprintf(openTree,   "%s%s",demo,channel); 
  TString CHANNEL = channel;    

  TFile *file01 = TFile::Open("../RISULTATI/LepIso_010413/data_SL.root");       TTree *Tree01=(TTree*)file01->Get(openTreeSB);   TTree *Tree14=(TTree*)file01->Get(openTree);
  TFile *file02 = TFile::Open("../RISULTATI/LepIso_010413/DY100_SL.root");      TTree *Tree02=(TTree*)file02->Get(openTreeSB);   TTree *Tree15=(TTree*)file02->Get(openTree); 
  TFile *file03 = TFile::Open("../RISULTATI/LepIso_010413/DY70_SL.root");       TTree *Tree03=(TTree*)file03->Get(openTreeSB);   TTree *Tree16=(TTree*)file03->Get(openTree); 
  TFile *file04 = TFile::Open("../RISULTATI/LepIso_010413/DYM50_100_SL.root");  TTree *Tree04=(TTree*)file04->Get(openTreeSB);   TTree *Tree17=(TTree*)file04->Get(openTree); 
  TFile *file05 = TFile::Open("../RISULTATI/LepIso_010413/DYM50_70_SL.root");   TTree *Tree05=(TTree*)file05->Get(openTreeSB);   TTree *Tree18=(TTree*)file05->Get(openTree); 
  TFile *file06 = TFile::Open("../RISULTATI/LepIso_010413/QCD1000_SL.root");    TTree *Tree06=(TTree*)file06->Get(openTreeSB);   TTree *Tree19=(TTree*)file06->Get(openTree); 
  TFile *file07 = TFile::Open("../RISULTATI/LepIso_010413/QCD250_SL.root");     TTree *Tree07=(TTree*)file07->Get(openTreeSB);   TTree *Tree20=(TTree*)file07->Get(openTree); 
  TFile *file08 = TFile::Open("../RISULTATI/LepIso_010413/QCD500_SL.root");     TTree *Tree08=(TTree*)file08->Get(openTreeSB);   TTree *Tree21=(TTree*)file08->Get(openTree); 
  TFile *file09 = TFile::Open("../RISULTATI/LepIso_010413/TT_SL.root");         TTree *Tree09=(TTree*)file09->Get(openTreeSB);   TTree *Tree22=(TTree*)file09->Get(openTree);
  TFile *file10 = TFile::Open("../RISULTATI/LepIso_010413/WJets180_SL.root");   TTree *Tree10=(TTree*)file10->Get(openTreeSB);   TTree *Tree23=(TTree*)file10->Get(openTree); 
  TFile *file11 = TFile::Open("../RISULTATI/LepIso_010413/WW_SL.root");         TTree *Tree11=(TTree*)file11->Get(openTreeSB);   TTree *Tree24=(TTree*)file11->Get(openTree);  
  TFile *file12 = TFile::Open("../RISULTATI/LepIso_010413/WZ_SL.root");         TTree *Tree12=(TTree*)file12->Get(openTreeSB);   TTree *Tree25=(TTree*)file12->Get(openTree);  
  TFile *file13 = TFile::Open("../RISULTATI/LepIso_010413/ZZ_SL.root");         TTree *Tree13=(TTree*)file13->Get(openTreeSB);   TTree *Tree26=(TTree*)file13->Get(openTree); 

  if(save){
    ofstream myfile;
    char saveName [50]; sprintf(saveName, "BKG_%s_SB_met%i.txt",channel,MET);
    myfile.open(saveName); 
  }

  for(int i=0; i<XMassBin; i++){
    TH1F *data      = new TH1F("","",bin,min,max);           TH1F *data_SR      = new TH1F("","",bin,min,max);
    TH1F *DY100     = new TH1F("","",bin,min,max);	     TH1F *DY100_SR     = new TH1F("","",bin,min,max);
    TH1F *DY70      = new TH1F("","",bin,min,max);	     TH1F *DY70_SR      = new TH1F("","",bin,min,max);
    TH1F *DYM50_100 = new TH1F("","",bin,min,max);	     TH1F *DYM50_100_SR = new TH1F("","",bin,min,max);
    TH1F *DYM50_70  = new TH1F("","",bin,min,max);	     TH1F *DYM50_70_SR  = new TH1F("","",bin,min,max);
    TH1F *QCD1000   = new TH1F("","",bin,min,max);	     TH1F *TT_SR        = new TH1F("","",bin,min,max);
    TH1F *QCD250    = new TH1F("","",bin,min,max);	     TH1F *QCD1000_SR   = new TH1F("","",bin,min,max);
    TH1F *QCD500    = new TH1F("","",bin,min,max);	     TH1F *QCD250_SR    = new TH1F("","",bin,min,max);
    TH1F *TT        = new TH1F("","",bin,min,max);	     TH1F *QCD500_SR    = new TH1F("","",bin,min,max);
    TH1F *WJets180  = new TH1F("","",bin,min,max);	     TH1F *WJets180_SR  = new TH1F("","",bin,min,max);
    TH1F *WW        = new TH1F("","",bin,min,max);	     TH1F *WW_SR        = new TH1F("","",bin,min,max);
    TH1F *WZ        = new TH1F("","",bin,min,max);	     TH1F *WZ_SR        = new TH1F("","",bin,min,max);
    TH1F *ZZ        = new TH1F("","",bin,min,max);	     TH1F *ZZ_SR        = new TH1F("","",bin,min,max);

    char CUT [200];
    char CUT_SR [200];
    if(!lepTauIso){
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo")sprintf(CUT, "trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && MassVis>10)",
							MET,XMassMin,XMassMax);
      else if(CHANNEL=="EleMuo")                sprintf(CUT, "trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",MET,XMassMin,XMassMax);
      else                                      sprintf(CUT, "trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",
							MET,XMassMin,XMassMax);
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo")sprintf(CUT_SR,"trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && MassVis>10)",
							MET,XMassMin,XMassMax);
      else if(CHANNEL=="EleMuo")                sprintf(CUT_SR,"trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",MET,XMassMin,XMassMax);
      else                                      sprintf(CUT_SR,"trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",
							MET,XMassMin,XMassMax);
    } else {
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo")sprintf(CUT, "trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && MassVis>10)",
							MET,XMassMin,XMassMax);
      else if(CHANNEL=="EleMuo") 	        sprintf(CUT, "trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",MET,XMassMin,XMassMax);
      else if(CHANNEL=="EleTau")  	        sprintf(CUT, "trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && lepCorrPFIso<0.1)",
							MET,XMassMin,XMassMax);
      else                                      sprintf(CUT, "trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && lepCorrPFIso<0.2)",
							MET,XMassMin,XMassMax);
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo")sprintf(CUT_SR,"trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && MassVis>10)",
							MET,XMassMin,XMassMax);
      else if(CHANNEL=="EleMuo")                sprintf(CUT_SR,"trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",MET,XMassMin,XMassMax);
      else if(CHANNEL=="EleTau")  	        sprintf(CUT_SR,"trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && lepCorrPFIso<0.1)",
							MET,XMassMin,XMassMax);
      else   	                                sprintf(CUT_SR,"trigger*PUWeight*(met>%f && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && lepCorrPFIso<0.2)",
							MET,XMassMin,XMassMax);
    }
 
    //PLOTS IN SIDEBAND
    char input01[50]; sprintf(input01,  "%s>>h01(%i,%f,%f)", plot,bin,min,max);
    char input02[50]; sprintf(input02,  "%s>>h02(%i,%f,%f)", plot,bin,min,max);
    char input03[50]; sprintf(input03,  "%s>>h03(%i,%f,%f)", plot,bin,min,max);
    char input04[50]; sprintf(input04,  "%s>>h04(%i,%f,%f)", plot,bin,min,max);
    char input05[50]; sprintf(input05,  "%s>>h05(%i,%f,%f)", plot,bin,min,max);
    char input06[50]; sprintf(input06,  "%s>>h06(%i,%f,%f)", plot,bin,min,max);
    char input07[50]; sprintf(input07,  "%s>>h07(%i,%f,%f)", plot,bin,min,max);
    char input08[50]; sprintf(input08,  "%s>>h08(%i,%f,%f)", plot,bin,min,max);
    char input09[50]; sprintf(input09,  "%s>>h09(%i,%f,%f)", plot,bin,min,max);
    char input10[50]; sprintf(input10,  "%s>>h10(%i,%f,%f)", plot,bin,min,max);
    char input11[50]; sprintf(input11,  "%s>>h11(%i,%f,%f)", plot,bin,min,max);
    char input12[50]; sprintf(input12,  "%s>>h12(%i,%f,%f)", plot,bin,min,max);
    char input13[50]; sprintf(input13,  "%s>>h13(%i,%f,%f)", plot,bin,min,max);
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
  
    //PLOTS IN SIGNAL REGION
    char input14 [50]; sprintf(input14, "%s>>h14(%i,%f,%f)", plot,bin,min,max);
    char input15 [50]; sprintf(input15, "%s>>h15(%i,%f,%f)", plot,bin,min,max);
    char input16 [50]; sprintf(input16, "%s>>h16(%i,%f,%f)", plot,bin,min,max);
    char input17 [50]; sprintf(input17, "%s>>h17(%i,%f,%f)", plot,bin,min,max);
    char input18 [50]; sprintf(input18, "%s>>h18(%i,%f,%f)", plot,bin,min,max);
    char input19 [50]; sprintf(input19, "%s>>h19(%i,%f,%f)", plot,bin,min,max);
    char input20 [50]; sprintf(input20, "%s>>h20(%i,%f,%f)", plot,bin,min,max);
    char input21 [50]; sprintf(input21, "%s>>h21(%i,%f,%f)", plot,bin,min,max);
    char input22 [50]; sprintf(input22, "%s>>h22(%i,%f,%f)", plot,bin,min,max);
    char input23 [50]; sprintf(input23, "%s>>h23(%i,%f,%f)", plot,bin,min,max);
    char input24 [50]; sprintf(input24, "%s>>h24(%i,%f,%f)", plot,bin,min,max);
    char input25 [50]; sprintf(input25, "%s>>h25(%i,%f,%f)", plot,bin,min,max);
    char input26 [50]; sprintf(input26, "%s>>h26(%i,%f,%f)", plot,bin,min,max);
    Tree14->Draw(input14,CUT_SR,"E"); if(Tree14->Draw(input14,CUT_SR,"E")) {data_SR      = h14; }
    Tree15->Draw(input15,CUT_SR);     if(Tree15->Draw(input15,CUT_SR))     {DY100_SR     = h15; }
    Tree16->Draw(input16,CUT_SR);     if(Tree16->Draw(input16,CUT_SR))     {DY70_SR      = h16; }
    Tree17->Draw(input17,CUT_SR);     if(Tree17->Draw(input17,CUT_SR))     {DYM50_100_SR = h17; }
    Tree18->Draw(input18,CUT_SR);     if(Tree18->Draw(input18,CUT_SR))     {DYM50_70_SR  = h18; }
    Tree19->Draw(input19,CUT_SR);     if(Tree19->Draw(input19,CUT_SR))     {QCD1000_SR   = h19; }
    Tree20->Draw(input20,CUT_SR);     if(Tree20->Draw(input20,CUT_SR))     {QCD250_SR    = h20; }
    Tree21->Draw(input21,CUT_SR);     if(Tree21->Draw(input21,CUT_SR))     {QCD500_SR    = h21; }
    Tree22->Draw(input22,CUT_SR);     if(Tree22->Draw(input22,CUT_SR))     {TT_SR        = h22; }
    Tree23->Draw(input23,CUT_SR);     if(Tree23->Draw(input23,CUT_SR))     {WJets180_SR  = h23; }
    Tree24->Draw(input24,CUT_SR);     if(Tree24->Draw(input24,CUT_SR))     {WW_SR        = h24; }
    Tree25->Draw(input25,CUT_SR);     if(Tree25->Draw(input25,CUT_SR))     {WZ_SR        = h25; }
    Tree26->Draw(input26,CUT_SR);     if(Tree26->Draw(input26,CUT_SR))     {WZ_SR        = h26; }

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
    float w_WJets180  = ( 23.500*19702./9739464.0);

    double N_sb = data->Integral();
    double N_sb_err = sqrt(data->Integral());
    double f_sb_err = 0;  
    double f_sb = 0; 
    double alpha = 0; 
    double alpha_err = 0;

    double f_sb_num = (w_DY100*DY100->Integral()  +  w_DY70*DY70->Integral()  +  w_TT*TT->Integral()  +
		       w_DYM50_100*DYM50_100->Integral()  +  w_DYM50_70*DYM50_70->Integral());
    double f_sb_num_err = sqrt(w_DY100*w_DY100*DY100->Integral()  +  w_DY70*w_DY70*DY70->Integral()  +  w_TT*w_TT*TT->Integral()  +
			       w_DYM50_100*w_DYM50_100*DYM50_100->Integral()  +  w_DYM50_70*w_DYM50_70*DYM50_70->Integral()); 
    double f_sb_den = (1.9*(w_QCD1000*QCD1000->Integral()  +  w_QCD250*QCD250->Integral()  +  w_QCD500*QCD500->Integral())  +
    		       w_WZ*WZ->Integral()  +  w_WW*WW->Integral()  +  w_ZZ*ZZ->Integral()  +  w_WJets180*WJets180->Integral());
    double f_sb_den_err = sqrt(1.9*1.9*(w_QCD1000*w_QCD1000*QCD1000->Integral()  +  w_QCD250*w_QCD250*QCD250->Integral()  +  w_QCD500*w_QCD500*QCD500->Integral())  +
			       w_WZ*w_WZ*WZ->Integral()  +  w_WW*w_WW*WW->Integral()  +  w_ZZ*w_ZZ*ZZ->Integral()  +  w_WJets180*w_WJets180*WJets180->Integral()); 
    double alpha_num = (w_DY100*DY100_SR->Integral()  +  w_DY70*DY70_SR->Integral()  +  w_TT*TT_SR->Integral()  +
			w_DYM50_100*DYM50_100_SR->Integral()  +  w_DYM50_70*DYM50_70_SR->Integral());
    double alpha_num_err = sqrt(w_DY100*w_DY100*DY100_SR->Integral() + w_DY70*w_DY70*DY70_SR->Integral() + w_TT*w_TT*TT_SR->Integral()  +
				w_DYM50_100*w_DYM50_100*DYM50_100_SR->Integral() + w_DYM50_70*w_DYM50_70*DYM50_70_SR->Integral());
    double alpha_den = f_sb_num;
    double alpha_den_err = f_sb_num_err;
    
    if((f_sb_num+f_sb_den)!=0) {
      f_sb = f_sb_num/(f_sb_num+f_sb_den);
      f_sb_err = sqrt(f_sb_den*f_sb_den*f_sb_num_err*f_sb_num_err + f_sb_num*f_sb_num*f_sb_den_err*f_sb_den_err)/((f_sb_num+f_sb_den)*(f_sb_num+f_sb_den));
    }
    
    if(alpha_den!=0) {
      alpha = alpha_num/alpha_den;
      alpha_err = sqrt(alpha_den*alpha_den*alpha_num_err*alpha_num_err + alpha_num*alpha_num*alpha_den_err*alpha_den_err)/(alpha_den*alpha_den);
    }
    
    double N_DY = N_sb * f_sb * alpha;
    double N_DY_err = sqrt(N_sb_err*N_sb_err*f_sb*f_sb*alpha*alpha + f_sb_err*f_sb_err*N_sb*N_sb*alpha*alpha + alpha_err*alpha_err*N_sb*N_sb*f_sb*f_sb);
    
    cout<<endl;
    cout<<"Number of sideband events is "<<N_sb<<" +/- "<<N_sb_err<<endl;
    cout<<"Ratio between DY events in signal and sideband region is "<<alpha<<" +/- "<<alpha_err<<endl;
    cout<<"Fraction of bkg events that are from DY process is "<<f_sb<<" +/- "<<f_sb_err<<endl;
    cout<<"Number of predicted DY+TTbar events in the signal region (in "<<XMassMin<<" <M(ZH)< "<<XMassMax<<" GeV) is "<<N_DY<<" +/- "<<N_DY_err<<endl;
    cout<<"Number of DY+TTbar events from MC in signal region (in "<<XMassMin<<" <M(ZH)< "<<XMassMax<<" GeV)       is "<<alpha_num<<" +/- "<<alpha_num_err<<endl;
    cout<<"Number of observed events in data in signal region (in "<<XMassMin<<" <M(ZH)< "<<XMassMax<<" GeV)       is "<<data_SR->Integral()<<" +/- "<<sqrt(data_SR->Integral())<<endl;
    cout<<endl;
    if(save){
      myfile<<"Number of sideband events is "<<N_sb<<" +/- "<<N_sb_err<<endl;
      myfile<<"Ratio between DY events in signal and sideband region is "<<alpha<<" +/- "<<alpha_err<<endl;
      myfile<<"Fraction of bkg events that are from DY process is "<<f_sb<<" +/- "<<f_sb_err<<endl;
      myfile<<"Number of predicted DY+TTbar events in the signal region (in "<<XMassMin<<" <M(ZH)< "<<XMassMax<<" GeV) is "<<N_DY<<" +/- "<<N_DY_err<<endl;
      myfile<<"Number of DY+TTbar events from MC in signal region (in "<<XMassMin<<" <M(ZH)< "<<XMassMax<<" GeV)       is "<<alpha_num<<" +/- "<<alpha_num_err<<endl;
      myfile<<"Number of observed events in data in signal region (in "<<XMassMin<<" <M(ZH)< "<<XMassMax<<" GeV)       is "<<data_SR->Integral()<<" +/- "<<sqrt(data_SR->Integral())<<endl;
      myfile<<endl;
    }


    histo_PRE->SetBinContent(i+1,N_DY  +  
			  1.9*(w_QCD1000*QCD1000_SR->Integral() + w_QCD250*QCD250_SR->Integral() + w_QCD500*QCD500_SR->Integral()) +
			  w_WZ*WZ_SR->Integral() + w_WW*WW_SR->Integral() + w_ZZ*ZZ_SR->Integral() + w_WJets180*WJets180_SR->Integral()
			  );
    histo_PRE->SetBinError(i+1,sqrt(N_DY_err*N_DY_err + 
				 1.9*1.9*(w_QCD1000*w_QCD1000*QCD1000_SR->Integral() + w_QCD250*w_QCD250*QCD250_SR->Integral() + w_QCD500*w_QCD500*QCD500_SR->Integral()) + 
				 w_WZ*w_WZ*WZ_SR->Integral()  +  w_WW*w_WW*WW_SR->Integral()  +  w_ZZ*w_ZZ*ZZ_SR->Integral()  +  w_WJets180*w_WJets180*WJets180_SR->Integral()
				 ));
    histo_DrY->SetBinContent(i+1,w_DY100*DY100_SR->Integral()+w_DY70*DY70_SR->Integral()+w_DYM50_100*DYM50_100_SR->Integral()+w_DYM50_70*DYM50_70_SR->Integral());
    histo_DrY->SetBinError(i+1,sqrt(w_DY100*w_DY100*DY100_SR->Integral()+w_DY70*w_DY70*DY70_SR->Integral()+
				 w_DYM50_100*w_DYM50_100*DYM50_100_SR->Integral()+w_DYM50_70*w_DYM50_70*DYM50_70_SR->Integral()));
    histo_TTb->SetBinContent(i+1,w_TT*TT_SR->Integral());
    histo_TTb->SetBinError(i+1,sqrt(w_TT*w_TT*TT_SR->Integral()));
    histo_OBS->SetBinContent(i+1,data_SR->Integral());
    histo_OBS->SetBinError(i+1,sqrt(data_SR->Integral()));
    histo_QCD->SetBinContent(i+1,1.9*(w_QCD1000*QCD1000_SR->Integral()+w_QCD500*QCD500_SR->Integral()+w_QCD250*QCD250_SR->Integral()));
    histo_WJe->SetBinContent(i+1,w_WZ*WZ_SR->Integral() + w_WW*WW_SR->Integral() + w_ZZ*ZZ_SR->Integral());
    histo_VVp->SetBinContent(i+1,w_WJets180*WJets180_SR->Integral());

    ERR->SetBinContent(i+1,w_DY100*DY100_SR->Integral()+w_DY70*DY70_SR->Integral()+w_DYM50_100*DYM50_100_SR->Integral()+w_DYM50_70*DYM50_70_SR->Integral()+
		       w_TT*TT_SR->Integral()+
		       1.9*(w_QCD1000*QCD1000_SR->Integral()+w_QCD500*QCD500_SR->Integral()+w_QCD250*QCD250_SR->Integral())+
		       w_WZ*WZ_SR->Integral() + w_WW*WW_SR->Integral() + w_ZZ*ZZ_SR->Integral()+w_WJets180*WJets180_SR->Integral()
		       );
    ERR->SetBinError(i+1,sqrt(w_DY100*w_DY100*DY100_SR->Integral()+w_DY70*w_DY70*DY70_SR->Integral()+
			      w_DYM50_100*w_DYM50_100*DYM50_100_SR->Integral()+w_DYM50_70*w_DYM50_70*DYM50_70_SR->Integral()+
			      w_TT*w_TT*TT_SR->Integral()+ 
			      1.9*1.9*(w_QCD1000*w_QCD1000*QCD1000_SR->Integral() + w_QCD250*w_QCD250*QCD250_SR->Integral() + w_QCD500*w_QCD500*QCD500_SR->Integral()) + 
			      w_WZ*w_WZ*WZ_SR->Integral()  +  w_WW*w_WW*WW_SR->Integral()  +  w_ZZ*w_ZZ*ZZ_SR->Integral()  +  w_WJets180*w_WJets180*WJets180_SR->Integral()
			      ));

    XMassMin=XMassMin+XMassWidth;
    XMassMax=XMassMax+XMassWidth;
  }

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
  TH1D *RATIO = new TH1D("","",histo_OBS->GetNbinsX(),histo_OBS->GetXaxis()->GetXmin(),histo_OBS->GetXaxis()->GetXmax());
  for(int m=1; m<histo_OBS->GetNbinsX()+1; m++){ 
    float A    = histo_OBS->GetBinContent(m);
    float B    = histo_PRE->GetBinContent(m);
    float Aerr = histo_OBS->GetBinError(m);
    float Berr = histo_PRE->GetBinError(m);
    //RATIO->SetBinContent(m,A-B);
    //RATIO->SetBinError(m,sqrt(Aerr*Aerr + Berr*Berr));
    if(histo_OBS->GetBinContent(m)!=0 && histo_PRE->GetBinContent(m)!=0) {
      RATIO->SetBinContent(m,A/B);
      RATIO->SetBinError(m,sqrt(Aerr*Aerr/(B*B) + A*A*Berr*Berr/(B*B*B*B)));
    }
  }
  RATIO->SetMarkerStyle(21);
  RATIO->SetMinimum(-3);
  RATIO->SetMaximum(5);
  RATIO->SetLineColor(1);
  RATIO->SetLineWidth(2);
  RATIO->GetXaxis()->SetTitleOffset(0.9);
  RATIO->GetXaxis()->SetLabelSize(0.08);
  RATIO->GetXaxis()->SetTitleSize(0.09);
  RATIO->GetYaxis()->SetTitleOffset(0.35);
  RATIO->GetYaxis()->SetTitle("observation/prediction");
  RATIO->GetXaxis()->SetTitle("M(Z,H) [GeV] (SvFit)");
  RATIO->GetYaxis()->SetLabelSize(0.08);
  RATIO->GetYaxis()->SetTitleSize(0.09);
  TLine* line = new TLine(histo_OBS->GetXaxis()->GetXmin(),1,histo_OBS->GetXaxis()->GetXmax(),1);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  
  TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.98,0.32);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetTopMargin(0.05);
  c1_1->SetBottomMargin(0.2);
  c1_1->SetRightMargin(0.02);
  c1_1->SetLeftMargin(0.07);
  RATIO->Draw("E1");
  line->Draw("same");	
  c1->cd();
	
  TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.98,0.96);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetTopMargin(0.01);
  c1_2->SetBottomMargin(0.1);
  c1_2->SetRightMargin(0.02);
  c1_2->SetLeftMargin(0.07);
  
  histo_QCD->SetFillColor(kRed-2);
  histo_QCD->SetLineColor(kRed-2);
  histo_WJe->SetFillColor(kMagenta-5);
  histo_WJe->SetLineColor(kMagenta-5);
  histo_VVp->SetFillColor(kGreen+1);
  histo_VVp->SetLineColor(kGreen+1);

  histo_OBS->SetFillColor(1);
  histo_OBS->SetLineColor(1);
  histo_OBS->SetLineWidth(2); 
  histo_OBS->SetLineColor(1);
  histo_OBS->SetMarkerColor(1); 
  histo_OBS->SetMarkerStyle(20); 
  histo_OBS->SetMarkerSize(1.3);

  histo_TTb->SetFillColor(kBlue);
  histo_TTb->SetLineColor(kBlue);

  histo_DrY->SetFillColor(kBlue-9);
  histo_DrY->SetLineColor(kBlue-9);

  histo_PRE->SetLineWidth(2); 
  histo_PRE->SetLineColor(2);
  histo_PRE->SetMarkerColor(2); 
  histo_PRE->SetMarkerStyle(20); 
  histo_PRE->SetMarkerSize(1.3);

  TF1 *func1 = new TF1("fit1","gaus*landau",0,XMassWidth*XMassBin);
  func1->SetLineColor(kRed);
  func1->SetLineWidth(3);
  TF1 *func4 = new TF1("fit4","gaus*landau",0,XMassWidth*XMassBin);
  func4->SetLineColor(1);
  func4->SetLineWidth(3);
  //histo_PRE->Fit("fit1");
  //histo_OBS->Fit("fit4");

  THStack *hs = new THStack("hs","hs");
  hs->SetTitle();
  hs->Add(histo_DrY);
  hs->Add(histo_VVp);
  hs->Add(histo_TTb);
  hs->Add(histo_QCD);
  hs->Add(histo_WJe);
  hs->Draw("histo");
  hs->GetYaxis()->SetTitleSize(0.045);
  hs->GetXaxis()->SetTitleSize(0.045);
  hs->GetYaxis()->SetLabelSize(0.045);
  hs->GetXaxis()->SetLabelSize(0.045);
  hs->GetYaxis()->SetTitleOffset(0.8); 
  hs->GetYaxis()->SetTitle("Events");
  hs->GetXaxis()->SetTitle("M(Z,H) [GeV] (SvFit)"); 
  if(CHANNEL=="MuoMuo") hs->SetMaximum(60);
  if(CHANNEL=="EleEle") hs->SetMaximum(50);
  if(CHANNEL=="EleMuo") hs->SetMaximum(40);
  if(CHANNEL=="EleTau") hs->SetMaximum(100);
  if(CHANNEL=="MuoTau") hs->SetMaximum(100);
  hs->GetXaxis()->SetNdivisions(509);
  hs->SetMinimum(0);

  //histo_DrY->Draw("histo");
  histo_PRE->Draw("E1 same");
  histo_OBS->Draw("E1 same");
  
  ERR->SetFillStyle(3005);
  ERR->SetFillColor(12);
  ERR->SetLineColor(12);
  ERR->Draw("E2same");

  TLegend *pl2 = new TLegend(0.55,0.65,0.975,0.98);
  pl2->SetTextSize(0.030); 
  pl2->SetFillColor(0);
  TLegendEntry *ple2 = pl2->AddEntry(histo_PRE, "Prediction from data-driven method",  "LP");
  ple2 = pl2->AddEntry(histo_OBS, "data",  "LP");
  ple2 = pl2->AddEntry(histo_DrY, "Drell-Yan prediction from MC",  "F");
  ple2 = pl2->AddEntry(histo_TTb, "ttbar prediction from MC",  "F");
  ple2 = pl2->AddEntry(histo_QCD, "QCD prediction from MC",  "F");
  ple2 = pl2->AddEntry(histo_WJe, "WJets(pt180) prediction from MC",  "F");
  ple2 = pl2->AddEntry(histo_VVp, "Diboson prediction from MC",  "F");
  pl2->Draw();
 
  if(save){
    char saveName2 [50]; sprintf(saveName2, "BKG_%s_SB_met%i.png",channel,MET);
    c1->SaveAs(saveName2);
    //c1->SaveAs("BKG_"+CHANNEL+".root");
    //c1->SaveAs("BKG_"+CHANNEL+".C");
    myfile.close();
  }
}
