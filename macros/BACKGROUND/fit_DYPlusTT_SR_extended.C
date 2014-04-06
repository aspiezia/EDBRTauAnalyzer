void fit_DYPlusTT_SR_extended(){
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5); //title X location 
  gStyle->SetTitleY(0.96); //title Y location 
  gStyle->SetPaintTextFormat(".2f");
	
  int SignalMass = 1000;
  float XMassWidth=200; float XMassBin=15; float XMassMin=0; float XMassMax=XMassWidth; float MAX = XMassWidth*XMassBin;
  bool save=false; bool pull=true; bool resid=false;
	
  for(int i=0; i<1; i++){
    //TString FUNCTION = "Gaussian";         char function[500]; sprintf(function,"Gaussian"); 
    //TString FUNCTION = "Landau";           char function[500]; sprintf(function,"Landau"); 
    TString FUNCTION = "BreitWigner";      char function[500]; sprintf(function,"BreitWigner"); 
    //TString FUNCTION = "CrystalBall";      char function[500]; sprintf(function,"CrystalBall"); 
    //TString FUNCTION = "GaussXLandau";     char function[500]; sprintf(function,"GaussXLandau"); 
    //TString FUNCTION = "Voigtian";         char function[500]; sprintf(function,"Voigtian"); 
    //TString FUNCTION = "GaussXExpo";       char function[500]; sprintf(function,"GaussXExpo"); 
    //TString FUNCTION = "GaussPlusGauss";   char function[500]; sprintf(function,"GaussPlusGauss"); 
    //TString FUNCTION = "LogNormal";        char function[500]; sprintf(function,"LogNormal"); 
    //TString FUNCTION = "BifurGauss";       char function[500]; sprintf(function,"BifurGauss"); 
    ////TString FUNCTION = "BreitWignerXPol0"; char function[500]; sprintf(function,"BreitWignerXPol0"); 
    ////TString FUNCTION = "Argus";            char function[500]; sprintf(function,"Argus"); 
    ////TString FUNCTION = "keys";             char function[500]; sprintf(function,"keys"); 
    ////TString FUNCTION = "GaussXPoly";       char function[500]; sprintf(function,"GaussXPoly"); 
    ////TString FUNCTION = "GaussXGauss";      char function[500]; sprintf(function,"GaussXGauss"); 
		
    TH1F *histo_PRE1 = new TH1F("histo_PRE1","histo_PRE1",XMassBin,0,MAX);
    TH1F *histo_MC1  = new TH1F("histo_MC1", "histo_MC1", XMassBin,0,MAX); 
    TH1F *histo_PRE2 = new TH1F("histo_PRE2","histo_PRE2",XMassBin,0,MAX);
    TH1F *histo_MC2  = new TH1F("histo_MC2", "histo_MC2", XMassBin,0,MAX); 
    TH1F *histo_PRE3 = new TH1F("histo_PRE3","histo_PRE3",XMassBin,0,MAX);
    TH1F *histo_MC3  = new TH1F("histo_MC3", "histo_MC3", XMassBin,0,MAX); 
		
    float BkgYield = 0;
		
    TCanvas* c1 = new TCanvas("c1","c1",0,0,600,600); 
    BackgroundEstimation("MuoMuo", XMassWidth, XMassBin, XMassMin, XMassMax, histo_PRE1, histo_MC1, BkgYield);
    Fitting(c1, MAX, histo_PRE1, histo_MC1, FUNCTION, XMassWidth, XMassBin, XMassMin, XMassMax, "MuoMuo", function, SignalMass, BkgYield, save, pull, resid);
		
    TCanvas* c2 = new TCanvas("c2","c2",0,0,600,600); 
    BackgroundEstimation("EleMuo", XMassWidth, XMassBin, XMassMin, XMassMax, histo_PRE2, histo_MC2, BkgYield);
    Fitting(c2, MAX, histo_PRE2, histo_MC2, FUNCTION, XMassWidth, XMassBin, XMassMin, XMassMax, "EleMuo", function, SignalMass, BkgYield, save, pull, resid);
    		
    TCanvas* c3 = new TCanvas("c3","c3",0,0,600,600); 
    BackgroundEstimation("EleEle", XMassWidth, XMassBin, XMassMin, XMassMax, histo_PRE3, histo_MC3, BkgYield);
    Fitting(c3, MAX, histo_PRE3, histo_MC3, FUNCTION, XMassWidth, XMassBin, XMassMin, XMassMax, "EleEle", function, SignalMass, BkgYield, save, pull, resid);
		
    SignalMass+=500;
  }
}



void Fitting(TCanvas* c, float MAX, TH1F *histo_PRE, TH1F *histo_MC, TString FUNCTION, float XMassWidth, float XMassBin,
	     float XMassMin, float XMassMax, char *channel, char *function, int SignalMass, float BkgYield, bool save, bool pull, bool resid){  
  char demo       [500]; sprintf(demo,       "demo/Tree"); 
  char openTree   [500]; sprintf(openTree,   "%s%s",demo,channel);  
  RooRealVar XMassSVFit("XMassSVFit","XMassSVFit",0,MAX); 
  TTree *TreeSig;
  double WeightSig = 0;
  RooFormulaVar *wFunc;
  double MEAN1 = 0; double MEAN2 = 0;
  if(SignalMass==1000)      {
    TFile *fileSig = TFile::Open("../RISULTATI/SB2_050214/ZH1000_FL.root");  
    TreeSig=(TTree*)fileSig->Get(openTree); 
    WeightSig = 1.54*19.7/17958;
    wFunc = new RooFormulaVar("w","w","1.54*19.7/17958",XMassSVFit);
    MEAN1 = 1000; MEAN2 = 700;
  } else if(SignalMass==1500) {
    TFile *fileSig = TFile::Open("../RISULTATI/SB2_050214/ZH1500_FL.root");  
    TreeSig=(TTree*)fileSig->Get(openTree); 
    WeightSig = 0.22*19.7/14099;
    wFunc = new RooFormulaVar("w","w","0.22*19.7/14099",XMassSVFit);
    MEAN1 = 1500; MEAN2 = 1200;
  } else if(SignalMass==2000) {
    TFile *fileSig = TFile::Open("../RISULTATI/SB2_050214/ZH2000_FL.root");  
    TreeSig=(TTree*)fileSig->Get(openTree); 
    WeightSig = 0.03*19.7/13080;
    wFunc = new RooFormulaVar("w","w","0.03*19.7/13080",XMassSVFit);
    MEAN1 = 2000; MEAN2 = 1700;
  } else if(SignalMass==2500) {
    TFile *fileSig = TFile::Open("../RISULTATI/SB2_050214/ZH2500_FL.root");  
    TreeSig=(TTree*)fileSig->Get(openTree); 
    WeightSig = 0.01*19.7/13091;
    wFunc = new RooFormulaVar("w","w","0.01*19.7/13091",XMassSVFit);
    MEAN1 = 2500; MEAN2 = 2200;
  }  
	
  char Gauss1_name[50]; sprintf(Gauss1_name,  "SigGauss1%s",channel); 
  char Gauss2_name[50]; sprintf(Gauss2_name,  "SigGauss2%s",channel); 
  char mean1_name[50];  sprintf(mean1_name,   "SigMean1%s", channel); 
  char mean2_name[50];  sprintf(mean2_name,   "SigMean2%s", channel); 
  char sigma1_name[50]; sprintf(sigma1_name,  "SigSigma1%s",channel); 
  char sigma2_name[50]; sprintf(sigma2_name,  "SigSigma2%s",channel); 
  char nsig1_name[50];  sprintf(nsig1_name,   "SigNsig1%s", channel); 
  char nsig2_name[50];  sprintf(nsig2_name,   "SigNsig2%s", channel);  
  char signal_name[50]; sprintf(signal_name,  "signal");  
	
  TTree *SignalTree = TreeSig->CopyTree("met>50 && uncorrmet>30 && trigger==1");
  SignalTree->SetWeight(WeightSig);
  RooRealVar PUWeight("PUWeight","PUWeight",0,100);
  RooArgList VarSet(XMassSVFit,PUWeight);
  RooPlot* frame = XMassSVFit.frame();
  RooDataSet* signalSampleUnw = new RooDataSet("signalSampleUnw","signalSampleUnw",SignalTree,VarSet,0,"PUWeight");
  RooRealVar* w = (RooRealVar*) signalSampleUnw->addColumn(*wFunc);
  RooDataSet signalSample("signalSample","signalSample",signalSampleUnw,*signalSampleUnw->get(),0,"w") ;
  RooRealVar mean1(mean1_name,mean1_name, MEAN1, SignalMass-300,SignalMass+300); 
  RooRealVar mean2(mean2_name,mean2_name, MEAN2, SignalMass-750,SignalMass-100); 
  RooRealVar sigma1(sigma1_name,sigma1_name, 0, 500); 
  RooRealVar sigma2(sigma2_name,sigma2_name, 0, 500); 
  RooGaussModel Gauss1(Gauss1_name,Gauss1_name,XMassSVFit,mean1,sigma1);
  RooGaussModel Gauss2(Gauss2_name,Gauss2_name,XMassSVFit,mean2,sigma2);
  RooRealVar nsig1(nsig1_name,nsig1_name,0.,100.); 
  RooRealVar nsig2(nsig2_name,nsig2_name,0.,100.); 
  RooAddPdf signal(signal_name,signal_name,RooArgList(Gauss1, Gauss2), RooArgList(nsig1,nsig2));
  RooFitResult* r = signal.fitTo(signalSample, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  SignalTree->Draw("XMassSVFit>>data_sig(15,0,3000)","PUWeight*(met>50)");
  float SigYield = data_sig->Integral();   
	
  signalSample.plotOn(frame,RooFit::Name("signalSample"),RooFit::Binning(25),RooFit::LineColor(2),RooFit::MarkerColor(2),RooFit::LineWidth(2));  
  signal.plotOn(frame,RooFit::LineColor(2),RooFit::Name("signal"));
  Double_t chi2 = frame->chiSquare("signal","signalSample",6);
  RooArgSet* params = signal.getVariables();
	
  float mean1_value  = params->getRealValue(mean1_name); 
  float mean2_value  = params->getRealValue(mean2_name); 
  float sigma1_value = params->getRealValue(sigma1_name); 
  float sigma2_value = params->getRealValue(sigma2_name); 
  float nsig1_value  = params->getRealValue(nsig1_name); 
  float nsig2_value  = params->getRealValue(nsig2_name); 
	
  float mean1_error  = mean1->getError(); 
  float mean2_error  = mean2->getError(); 
  float sigma1_error = sigma1->getError(); 
  float sigma2_error = sigma2->getError(); 
  float nsig1_error  = nsig1->getError(); 
  float nsig2_error  = nsig2->getError(); 
	
  frame->Draw();
  frame->SetMinimum(0);
  cout<<endl;
  cout<<"SIGNAL MODEL --- m = "<<SignalMass<<" GeV"<<endl;
  params->Print("V");
  cout<<"  X) chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  TString MASS = "0";
  if(SignalMass==1000)      MASS="1000";
  else if(SignalMass==1500) MASS="1500";
  else if(SignalMass==2000) MASS="2000";
  else if(SignalMass==2500) MASS="2500";
	
  TString CHANNEL = channel; 
  if(save) c->SaveAs("fit_signal_"+CHANNEL+"_"+MASS+".png");
	
	
  TFile *fileDat = TFile::Open("../RISULTATI/SB2_050214/data_FL.root"); 
  TTree *TreeDat=(TTree*)fileDat->Get(openTree);  
  RooDataHist dataset1("s","s",XMassSVFit,histo_PRE);
  RooDataHist dataset2("s","s",XMassSVFit,histo_MC);
	
  RooPlot* frame1 = XMassSVFit.frame();
  if(FUNCTION=="Gaussian")  Gaussian(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
				     mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
				     mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="Landau") Landau(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
				mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
				mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="CrystalBall")     CrystalBall(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 5, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					      mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					      mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="Argus")  Argus(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 4, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
			       mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
			       mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="BreitWigner")  BreitWigner(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					   mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					   mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="Voigtian")  Voigtian(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 4, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
				     mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
				     mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  //if(FUNCTION=="keys")  KeysPdf(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
  //mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
  //mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussXLandau")  GaussXLandau(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 5, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					     mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					     mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussXGauss") GaussXGauss(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 5, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					  mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					  mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussXPoly") GaussXPoly(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 6, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussXExpo") GaussXExpo(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 4, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussPlusGauss") GaussPlusGauss(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 6, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
						mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
						mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="BifurGauss") BifurGauss(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="LogNormal") LogNormal(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
				      mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
				      mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="BreitWignerXPol0") BreitWignerXPol0(frame1, dataset1, XMassSVFit, XMassBin, MAX, 2, 4, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
						    mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
						    mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  frame1->SetTitle(0);
  frame1->SetXTitle("M(Z,H) [GeV]");
  frame1->SetMinimum(0);
  frame1->SetLineColor(2);
  frame1->Draw();
	
  RooPlot* frame2 = XMassSVFit.frame(); 
  if(FUNCTION=="Gaussian")  Gaussian(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
				     mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
				     mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="Landau") Landau(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
				mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
				mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="CrystalBall")     CrystalBall(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 5, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					      mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					      mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="Argus")  Argus(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 4, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
			       mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
			       mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="BreitWigner")  BreitWigner(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					   mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					   mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="Voigtian")  Voigtian(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
				     mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
				     mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  //if(FUNCTION=="keys")  KeysPdf(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
  //mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
  //mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussXLandau")  GaussXLandau(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 5, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					     mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					     mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussXGauss") GaussXGauss(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 5, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					  mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					  mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussXPoly") GaussXPoly(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 6, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussXExpo") GaussXExpo(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 4, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="GaussPlusGauss") GaussPlusGauss(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 6, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
						mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
						mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="BifurGauss") BifurGauss(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
					mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
					mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="LogNormal") LogNormal(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 3, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
				      mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
				      mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  if(FUNCTION=="BreitWignerXPol0") BreitWignerXPol0(frame2, dataset2, XMassSVFit, XMassBin, MAX, 1, 4, TreeDat, signal, channel, function, SignalMass, BkgYield, SigYield,
						    mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, 
						    mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error);
  frame2->Draw("same");
	
  SavePlot(c, histo_PRE, histo_MC, frame1, frame2, XMassSVFit, CHANNEL, FUNCTION, save, pull, resid, MAX);
}

void SavePlot(TCanvas* c, TH1F *histo_PRE, TH1F *histo_MC, RooPlot* frame1, RooPlot* frame2, RooRealVar XMassSVFit, 
	      TString CHANNEL, TString FUNCTION, bool save, bool pull, bool resid, float MAX){
	
  TLegend *pl2 = new TLegend(0.41,0.80,0.97,0.98);
  pl2->SetTextSize(0.043); 
  pl2->SetFillColor(0);
  TLegendEntry *ple2 = pl2->AddEntry(histo_PRE, "Prediction from data-driven method",  "LP");
  ple2 = pl2->AddEntry(histo_MC, "Prediction from MC",  "LP");
  pl2->Draw();
	
  c->cd();
  TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.98,0.32);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetTopMargin(0.05);
  c1_1->SetBottomMargin(0.2);
  c1_1->SetRightMargin(0.02);
  c1_1->SetLeftMargin(0.07);
  RooHist* hresid1 = frame1->residHist(); //Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hpull1  = frame1->pullHist();  //Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hresid2 = frame2->residHist(); //Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hpull2  = frame2->pullHist();  //Construct a histogram with the pulls of the data w.r.t the curve
  RooPlot* frame3 = XMassSVFit.frame(RooFit::Title("")); //Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot* frame4 = XMassSVFit.frame(RooFit::Title("")); //Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* frame5 = XMassSVFit.frame(RooFit::Title("")); //Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot* frame6 = XMassSVFit.frame(RooFit::Title("")); //Create a new frame to draw the pull distribution and add the distribution to the frame
  frame3->addPlotable(hresid1,"P");
  frame4->addPlotable(hpull1,"P");
  frame5->addPlotable(hresid2,"P");
  frame6->addPlotable(hpull2,"P");
  frame3->SetMinimum(-15);  frame3->SetMaximum(+15);
  frame4->SetMinimum(-5);   frame4->SetMaximum(+5);
  frame3->SetXTitle("M(Z,H) [GeV]");
  frame4->SetXTitle("M(Z,H) [GeV]");
  frame3->SetYTitle("residual");
  frame4->SetYTitle("pull");    
  frame3->SetTitle("");
  frame4->SetTitle("");    
  frame3->GetYaxis()->SetTitleSize(0.085);    frame4->GetYaxis()->SetTitleSize(0.085); 
  frame3->GetXaxis()->SetTitleSize(0.085);    frame4->GetXaxis()->SetTitleSize(0.085); 
  frame3->GetYaxis()->SetLabelSize(0.075);    frame4->GetYaxis()->SetLabelSize(0.075); 
  frame3->GetXaxis()->SetLabelSize(0.075);    frame4->GetXaxis()->SetLabelSize(0.075); 
  frame3->GetYaxis()->SetTitleOffset(0.4);    frame4->GetYaxis()->SetTitleOffset(0.3); 
  hresid1->SetLineColor(2);  hresid1->SetMarkerColor(2);  hresid1->SetLineWidth(2);
  hpull1->SetLineColor(2);   hpull1->SetMarkerColor(2);   hpull1->SetLineWidth(2);
  hresid2->SetLineColor(1);  hresid2->SetMarkerColor(1);  hresid2->SetLineWidth(2);
  hpull2->SetLineColor(1);   hpull2->SetMarkerColor(1);   hpull2->SetLineWidth(2);
  if(pull){
    frame4->Draw();
    frame6->Draw("same");
  }
  if(resid){
    frame3->Draw();
    frame5->Draw("same");
  }
	
  TLine* line = new TLine(0,0,MAX,0);
  line->SetLineColor(1);
  line->SetLineWidth(2);
  line->Draw("same");
	
  if(save){
    c->SaveAs("fit_DYPlusTT_"+CHANNEL+"_"+FUNCTION+"_SR.png");
    //c->SaveAs("DYPlusTT_"+CHANNEL+".root");
    //c->SaveAs("DYPlusTT_"+CHANNEL+".C");
  }
}

void SaveWorkspace(TTree *Tree1, RooAddPdf background, RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
		   float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		   float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error, 
		   float mean_value, float sigma_value, float nbkg_value, float mean_error, float sigma_error, float nbkg_error){
  char saveName [50]; sprintf(saveName, "simple-shapes_%s_%s_%i.root",channel,function,SignalMass);
  RooWorkspace *w = new RooWorkspace("w","w");
  TTree *SelectedTree = Tree1->CopyTree("met>50 && uncorrmet>30 && trigger==1");
  RooRealVar XMassSVFit("XMassSVFit","XMassSVFit",0,10000);
  RooDataSet data_obs("data_obs","data_obs",XMassSVFit,RooFit::Import(*SelectedTree)) ;
  w->var("XMassSVFit[0,3000]");
  w->importClassCode("RooAddPdf",kTRUE);
  w->importClassCode("RooDataSet",kTRUE);
  w->import(data_obs);
  w->import(signal);
  w->import(background);
  w->writeToFile(saveName);
  int data_obs_entries = SelectedTree->GetEntries();
  SaveDatacard(channel, function, SignalMass, data_obs_entries, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void SaveDatacard(char *channel, char *function, int SignalMass, int data_obs_entries, float BkgYield, float SigYield,
		  float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		  float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error, 
		  float mean_value, float sigma_value, float nbkg_value, float mean_error, float sigma_error, float nbkg_error){
  char saveName [50]; sprintf(saveName, "simple-shapes_%s_%s_%i.txt",channel,function,SignalMass);
  ofstream myfile;

  myfile.open(saveName); 
  myfile<<"imax 1"<<endl;
  myfile<<"jmax 1"<<endl;
  myfile<<"kmax *"<<endl;
  myfile<<"---------------"<<endl;
  myfile<<"shapes * "<<channel<<" simple-shapes_"<<channel<<"_"<<function<<"_"<<SignalMass<<".root w:$PROCESS"<<endl;
  myfile<<"---------------"<<endl;
  myfile<<"bin "<<channel<<endl;
  myfile<<"observation -1"<<endl;
  myfile<<"------------------------------"<<endl;
  myfile<<"bin        "<<channel<<"          "<<channel<<endl;
  myfile<<"process      signal     background"<<endl;
  myfile<<"process      0          1"<<endl;
  myfile<<"rate         "<<SigYield<<"   "<<BkgYield<<endl;
  myfile<<"--------------------------------"<<endl;
  myfile<<"BkgMean"  <<channel<<" param "<<mean_value  <<" "<<1+mean_error/mean_value  <<endl;
  myfile<<"BkgSigma" <<channel<<" param "<<sigma_value <<" "<<1+sigma_error/sigma_value <<endl;
  myfile<<"BkgNbkg"  <<channel<<" param "<<nbkg_value  <<" "<<1+nbkg_error/nbkg_value  <<endl;
  myfile<<"SigMean1" <<channel<<" param "<<mean1_value <<" "<<1+mean1_error/mean1_value <<endl;
  myfile<<"SigMean2" <<channel<<" param "<<mean2_value <<" "<<1+mean2_error/mean2_value <<endl;
  myfile<<"SigNsig1" <<channel<<" param "<<nsig1_value <<" "<<1+nsig1_error/nsig1_value <<endl;
  myfile<<"SigNsig2" <<channel<<" param "<<nsig2_value <<" "<<1+nsig2_error/nsig2_value <<endl;
  myfile<<"SigSigma1"<<channel<<" param "<<sigma1_value<<" "<<1+sigma1_error/sigma1_value<<endl;
  myfile<<"SigSigma2"<<channel<<" param "<<sigma2_value<<" "<<1+sigma2_error/sigma2_value<<endl;
  myfile.close();
	
}


void BackgroundEstimation(char *channel, float XMassWidth, float XMassBin, float XMassMin, float XMassMax, TH1F *histo_PRE, TH1F *histo_MC, float & BkgYield){
  char *plot = "met";
  bool lepTauIso = false;
  int bin=1; 
  float min=0; 
  float max=50000; 
  char demoSB     [500]; sprintf(demoSB,     "demo/TreeSB1"); 
  char openTreeSB [500]; sprintf(openTreeSB, "%s%s",demoSB,channel);
  char demo       [500]; sprintf(demo,       "demo/Tree"); 
  char openTree   [500]; sprintf(openTree,   "%s%s",demo,channel); 
  TString CHANNEL = channel; 
	
  TFile *file01 = TFile::Open("../RISULTATI/SB2_050214/data_FL.root");       TTree *Tree01=(TTree*)file01->Get(openTreeSB);   TTree *Tree14=(TTree*)file01->Get(openTree);
  TFile *file02 = TFile::Open("../RISULTATI/SB2_050214/DY100_FL.root");      TTree *Tree02=(TTree*)file02->Get(openTreeSB);   TTree *Tree15=(TTree*)file02->Get(openTree); 
  TFile *file03 = TFile::Open("../RISULTATI/SB2_050214/DY70_FL.root");       TTree *Tree03=(TTree*)file03->Get(openTreeSB);   TTree *Tree16=(TTree*)file03->Get(openTree); 
  TFile *file04 = TFile::Open("../RISULTATI/SB2_050214/DYM50_100_FL.root");  TTree *Tree04=(TTree*)file04->Get(openTreeSB);   TTree *Tree17=(TTree*)file04->Get(openTree); 
  TFile *file05 = TFile::Open("../RISULTATI/SB2_050214/DYM50_70_FL.root");   TTree *Tree05=(TTree*)file05->Get(openTreeSB);   TTree *Tree18=(TTree*)file05->Get(openTree); 
  TFile *file06 = TFile::Open("../RISULTATI/SB2_050214/QCD1000_FL.root");    TTree *Tree06=(TTree*)file06->Get(openTreeSB);   TTree *Tree19=(TTree*)file06->Get(openTree); 
  TFile *file07 = TFile::Open("../RISULTATI/SB2_050214/QCD250_FL.root");     TTree *Tree07=(TTree*)file07->Get(openTreeSB);   TTree *Tree20=(TTree*)file07->Get(openTree); 
  TFile *file08 = TFile::Open("../RISULTATI/SB2_050214/QCD500_FL.root");     TTree *Tree08=(TTree*)file08->Get(openTreeSB);   TTree *Tree21=(TTree*)file08->Get(openTree); 
  TFile *file09 = TFile::Open("../RISULTATI/SB2_050214/TT_FL.root");         TTree *Tree09=(TTree*)file09->Get(openTreeSB);   TTree *Tree22=(TTree*)file09->Get(openTree);
  TFile *file10 = TFile::Open("../RISULTATI/SB2_050214/WJets180_FL.root");   TTree *Tree10=(TTree*)file10->Get(openTreeSB);   TTree *Tree23=(TTree*)file10->Get(openTree); 
  TFile *file11 = TFile::Open("../RISULTATI/SB2_050214/WW_FL.root");         TTree *Tree11=(TTree*)file11->Get(openTreeSB);   TTree *Tree24=(TTree*)file11->Get(openTree);  
  TFile *file12 = TFile::Open("../RISULTATI/SB2_050214/WZ_FL.root");         TTree *Tree12=(TTree*)file12->Get(openTreeSB);   TTree *Tree25=(TTree*)file12->Get(openTree);  
  TFile *file13 = TFile::Open("../RISULTATI/SB2_050214/ZZ_FL.root");         TTree *Tree13=(TTree*)file13->Get(openTreeSB);   TTree *Tree26=(TTree*)file13->Get(openTree); 
	
  TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.98,0.96);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetTopMargin(0.01);
  c1_2->SetBottomMargin(0.1);
  c1_2->SetRightMargin(0.02);
  c1_2->SetLeftMargin(0.07);
	
  for(int i=0; i<XMassBin; i++){
    TH1F *data      = new TH1F("","",bin,min,max);       TH1F *data_SR      = new TH1F("","",bin,min,max);
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
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo")sprintf(CUT, "trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && MassVis>10)",XMassMin,XMassMax);
      else if(CHANNEL=="EleMuo")                sprintf(CUT, "trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",XMassMin,XMassMax);
      else                                      sprintf(CUT, "trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && charge==-1)",XMassMin,XMassMax);
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo")sprintf(CUT_SR,"trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && MassVis>10)",XMassMin,XMassMax);
      else if(CHANNEL=="EleMuo")                sprintf(CUT_SR,"trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",XMassMin,XMassMax);
      else                                      sprintf(CUT_SR,"trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && charge==-1)",XMassMin,XMassMax);
    } else {
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo")sprintf(CUT, "trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && MassVis>10)",XMassMin,XMassMax);
      else if(CHANNEL=="EleMuo") 	        sprintf(CUT, "trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",XMassMin,XMassMax);
      else if(CHANNEL=="EleTau")  	        sprintf(CUT, "trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && lepCorrPFIso<0.1 && charge==-1)",
							XMassMin,XMassMax);
      else                                      sprintf(CUT, "trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && lepCorrPFIso<0.2 && charge==-1)",
							XMassMin,XMassMax);
      if(CHANNEL=="EleEle" || CHANNEL=="MuoMuo")sprintf(CUT_SR,"trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && MassVis>10)",XMassMin,XMassMax);
      else if(CHANNEL=="EleMuo")                sprintf(CUT_SR,"trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f)",XMassMin,XMassMax);
      else if(CHANNEL=="EleTau")  	        sprintf(CUT_SR,"trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && lepCorrPFIso<0.1 && charge==-1)",
							XMassMin,XMassMax);
      else   	                                sprintf(CUT_SR,"trigger*PUWeight*(met>50 && uncorrmet>30 && XMassSVFit>%f && XMassSVFit<%f && lepCorrPFIso<0.2 && charge==-1)",
							XMassMin,XMassMax);
    }
		
    //PLOTS IN SIDEBAND
    char input01[50];   sprintf(input01,  "%s>>h01(%i,%f,%f)", plot,bin,min,max);
    char input02[50];   sprintf(input02,  "%s>>h02(%i,%f,%f)", plot,bin,min,max);
    char input03[50];   sprintf(input03,  "%s>>h03(%i,%f,%f)", plot,bin,min,max);
    char input04[50];   sprintf(input04,  "%s>>h04(%i,%f,%f)", plot,bin,min,max);
    char input05[50];   sprintf(input05,  "%s>>h05(%i,%f,%f)", plot,bin,min,max);
    char input06[50];   sprintf(input06,  "%s>>h06(%i,%f,%f)", plot,bin,min,max);
    char input07[50];   sprintf(input07,  "%s>>h07(%i,%f,%f)", plot,bin,min,max);
    char input08[50];   sprintf(input08,  "%s>>h08(%i,%f,%f)", plot,bin,min,max);
    char input09[50];   sprintf(input09,  "%s>>h09(%i,%f,%f)", plot,bin,min,max);
    char input10[50];   sprintf(input10,  "%s>>h10(%i,%f,%f)", plot,bin,min,max);
    char input11[50];   sprintf(input11,  "%s>>h11(%i,%f,%f)", plot,bin,min,max);
    char input12[50];   sprintf(input12,  "%s>>h12(%i,%f,%f)", plot,bin,min,max);
    char input13[50];   sprintf(input13,  "%s>>h13(%i,%f,%f)", plot,bin,min,max);
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
    cout<<endl;
		
    histo_PRE->SetBinContent(i+1,N_DY  +  
			     1.9*(w_QCD1000*QCD1000_SR->Integral() + w_QCD250*QCD250_SR->Integral() + w_QCD500*QCD500_SR->Integral()) +
			     w_WZ*WZ_SR->Integral() + w_WW*WW_SR->Integral() + w_ZZ*ZZ_SR->Integral() + w_WJets180*WJets180_SR->Integral()
			     );
    histo_PRE->SetBinError(i+1,sqrt(N_DY_err*N_DY_err + 
				    1.9*1.9*(w_QCD1000*w_QCD1000*QCD1000_SR->Integral() + w_QCD250*w_QCD250*QCD250_SR->Integral() + w_QCD500*w_QCD500*QCD500_SR->Integral()) + 
				    w_WZ*w_WZ*WZ_SR->Integral()  +  w_WW*w_WW*WW_SR->Integral()  +  w_ZZ*w_ZZ*ZZ_SR->Integral()  +  w_WJets180*w_WJets180*WJets180_SR->Integral()
				    ));
		
    histo_MC->SetBinContent(i+1,w_DY100*DY100_SR->Integral()+w_DY70*DY70_SR->Integral()+w_DYM50_100*DYM50_100_SR->Integral()+w_DYM50_70*DYM50_70_SR->Integral()+
			    w_TT*TT_SR->Integral()+
			    1.9*(w_QCD1000*QCD1000_SR->Integral()+w_QCD500*QCD500_SR->Integral()+w_QCD250*QCD250_SR->Integral())+
			    w_WZ*WZ_SR->Integral() + w_WW*WW_SR->Integral() + w_ZZ*ZZ_SR->Integral()+w_WJets180*WJets180_SR->Integral()
			    );
    histo_MC->SetBinError(i+1,sqrt(w_DY100*w_DY100*DY100_SR->Integral()+w_DY70*w_DY70*DY70_SR->Integral()+
				   w_DYM50_100*w_DYM50_100*DYM50_100_SR->Integral()+w_DYM50_70*w_DYM50_70*DYM50_70_SR->Integral()+
				   w_TT*w_TT*TT_SR->Integral()+ 
				   1.9*1.9*(w_QCD1000*w_QCD1000*QCD1000_SR->Integral() + w_QCD250*w_QCD250*QCD250_SR->Integral() + w_QCD500*w_QCD500*QCD500_SR->Integral()) + 
				   w_WZ*w_WZ*WZ_SR->Integral()  +  w_WW*w_WW*WW_SR->Integral()  +  w_ZZ*w_ZZ*ZZ_SR->Integral()  +  w_WJets180*w_WJets180*WJets180_SR->Integral()
				   ));
		
    XMassMin=XMassMin+XMassWidth;
    XMassMax=XMassMax+XMassWidth;
  }
	
  histo_PRE->SetLineWidth(2); 
  histo_PRE->SetLineColor(2);
  histo_PRE->SetMarkerColor(2); 
  histo_PRE->SetMarkerStyle(20); 
  histo_PRE->SetMarkerSize(1.3);
	
  histo_MC->SetLineColor(1);
  histo_MC->SetLineWidth(2);
  histo_MC->SetMarkerColor(1); 
  histo_MC->SetMarkerStyle(20); 
  histo_MC->SetMarkerSize(1.3);
  BkgYield = histo_PRE->Integral();
}


void Gaussian(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	      RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	      float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
	      float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  char bkg_name[50];        sprintf(bkg_name,         "Bkg%s",         channel); 
  char mean_name[50];       sprintf(mean_name,        "BkgMean%s",     channel); 
  char sigma_name[50];      sprintf(sigma_name,       "BkgSigma%s",    channel); 
  char nbkg_name[50];       sprintf(nbkg_name,        "BkgNbkg%s",     channel); 
  char background_name[50]; sprintf(background_name,  "background");  
  RooRealVar mean(mean_name,mean_name, 800,500,3000); 
  RooRealVar sigma(sigma_name,sigma_name, 700, 0,1000); 
  RooGaussModel bkg(bkg_name,bkg_name,XMassSVFit,mean,sigma);
  RooRealVar nbkg(nbkg_name,nbkg_name,100,0.,100000.); 
  RooAddPdf background(background_name,background_name,bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  float mean_value  = params->getRealValue(mean_name); 
  float sigma_value = params->getRealValue(sigma_name); 
  float nbkg_value  = params->getRealValue(nbkg_name);
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void Landau(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	    RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	    float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
	    float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 1000,700,2500); 
  RooRealVar sigma("sigma", "sigma", 100,0,2500); 
  RooLandau  bkg("bkg","bkg",XMassSVFit,mean,sigma);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void CrystalBall(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		 RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
		 float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		 float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 1000,700,2500); 
  RooRealVar sigma("sigma", "sigma", 100,0,2500); 
  RooRealVar alpha("alpha", "alpha", -1000,2500); 
  RooRealVar n("n", "n", 0,10); 
  RooCBShape bkg("bkg","bkg",XMassSVFit,mean,sigma,alpha,n);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void Argus(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	   RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	   float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
	   float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 7000,100,9000); 
  RooRealVar sigma("sigma", "sigma", -40,-100,-10); 
  RooRealVar alpha("alpha", "alpha", 10,0,250); 
  RooArgusBG bkg("bkg","bkg",XMassSVFit,mean,sigma,alpha);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void BreitWigner(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		 RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
		 float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		 float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  char bkg_name[50];        sprintf(bkg_name,         "Bkg%s",         channel); 
  char mean_name[50];       sprintf(mean_name,        "BkgMean%s",     channel); 
  char sigma_name[50];      sprintf(sigma_name,       "BkgSigma%s",    channel); 
  char nbkg_name[50];       sprintf(nbkg_name,        "BkgNbkg%s",     channel); 
  char background_name[50]; sprintf(background_name,  "background");  
  RooRealVar mean(mean_name,mean_name, 700,0,3000); 
  RooRealVar sigma(sigma_name,sigma_name, 700, 0,1000); 
  RooBreitWigner bkg(bkg_name,bkg_name,XMassSVFit,mean,sigma);
  RooRealVar nbkg(nbkg_name,nbkg_name,100,0.,100000.); 
  RooAddPdf background(background_name,background_name,bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  float mean_value  = params->getRealValue(mean_name); 
  float sigma_value = params->getRealValue(sigma_name); 
  float nbkg_value  = params->getRealValue(nbkg_name);
  float mean_error  = mean->getError();
  float sigma_error = sigma->getError();
  float nbkg_error  = nbkg->getError();
  //cout<<mean_value<<" "<<sigma_value<<" "<<nbkg_value<<endl;
  //cout<<mean->getError()<<" "<<sigma->getError()<<" "<<nbkg->getError()<<endl;
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void Voigtian(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	      RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	      float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
	      float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 700,1100); 
  RooRealVar width("width", "width", 200, 0,500); 
  RooRealVar sigma("sigma", "sigma", 100, 300); 
  RooVoigtian bkg("bkg","bkg",XMassSVFit,mean,width,sigma);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

//void KeysPdf(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1,
//	       RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
//             float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
//             float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
//  RooKeysPdf bkg("bkg","bkg",XMassSVFit,dataset,RooKeysPdf::MirrorBoth,1.6);
//  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
//  RooAddPdf background("background","background",bkg,nbkg);
//  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
//  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
//  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
//  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
//  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
//  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
//  Double_t chi2 = frame->chiSquare();
//  RooArgSet* params = background.getVariables();
//  cout<<endl;
//  cout<<"BACKGROUND MODEL"<<endl;
//  params->Print("V");
//  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
//  cout<<endl;
//  float mean_value=0; float sigma_value=0; float nbkg_value=0;
//  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
//}

void GaussXLandau(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		  RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	      float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		  float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  char bkgG_name[50];       sprintf(bkgG_name,        "BkgG%s",        channel); 
  char bkgL_name[50];       sprintf(bkgL_name,        "BkgL%s",        channel); 
  char bkg_name[50];        sprintf(bkg_name,         "Bkg%s",         channel); 
  char meanG_name[50];      sprintf(meanG_name,       "BkgMeanG%s",    channel); 
  char meanL_name[50];      sprintf(meanL_name,       "BkgMeanL%s",    channel); 
  char sigmaG_name[50];     sprintf(sigmaG_name,      "BkgSigmaG%s",   channel); 
  char sigmaL_name[50];     sprintf(sigmaL_name,      "BkgSigmaL%s",   channel); 
  char nbkg_name[50];       sprintf(nbkg_name,        "BkgNbkg%s",     channel); 
  char background_name[50]; sprintf(background_name,  "background");  
  RooRealVar meanG(meanG_name,meanG_name, 1000,700,2500); 
  RooRealVar meanL(meanL_name,meanL_name, -100,-80); 
  RooRealVar sigmaG(sigmaG_name,sigmaG_name, 100, 0,2500); 
  RooRealVar sigmaL(sigmaL_name,sigmaL_name, 0,2500); 
  RooLandau  Landau(bkgL_name,bkgL_name,XMassSVFit,meanL,sigmaL);
  RooGaussModel Gauss(bkgG_name,bkgG_name,XMassSVFit,meanG,sigmaG);
  RooFFTConvPdf bkg(bkg_name,bkg_name,XMassSVFit,Gauss,Landau);
  RooRealVar nbkg(nbkg_name,nbkg_name,100,0.,100000.); 
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void GaussXGauss(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		 RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	      float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		 float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean1("mean1", "mean1", 700,2000); 
  RooRealVar mean2("mean2", "mean2", 700,2000); 
  RooRealVar sigma1("sigma1", "sigma1", 0, 3000); 
  RooRealVar sigma2("sigma2", "sigma2", 0, 3000); 
  RooGaussModel Gauss1("Gauss1","Gauss1",XMassSVFit,mean1,sigma1);
  RooGaussModel Gauss2("Gauss2","Gauss2",XMassSVFit,mean2,sigma2);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooFFTConvPdf bkg("bkg","bkg",XMassSVFit,Gauss1,Gauss2);
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void GaussXExpo(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	      float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 0,3000); 
  RooRealVar sigma("sigma", "sigma", 150, 300); 
  RooRealVar alpha("alpha", "alpha", -1,1); 
  RooGaussModel Gauss("Gauss","Gauss",XMassSVFit,mean,sigma);
  RooExponential Expo("Expo","Expo",XMassSVFit,alpha);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooFFTConvPdf bkg("bkg","bkg",XMassSVFit,Gauss,Expo);
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void GaussXPoly(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	      float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 0,3000); 
  RooRealVar sigma("sigma", "sigma", 0, 3000); 
  RooRealVar c1("c1", "c1", -1000,1000); 
  RooRealVar c2("c2", "c2", -1000,1000); 
  RooRealVar c3("c3", "c3", -1000,1000); 
  RooGaussModel Gauss("Gauss","Gauss",XMassSVFit,mean,sigma);
  RooPolynomial Poly("Poly","Poly",XMassSVFit,RooArgList(c1,c2,c3));
  //RooChebychev Poly("bkg","background p.d.f.",XMassSVFit,RooArgList(c1,c2,c3));
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooFFTConvPdf bkg("bkg","bkg",XMassSVFit,Gauss,Poly);
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void GaussPlusGauss(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		    RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	      float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		    float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean1("mean1", "mean1", 700,2000); 
  RooRealVar mean2("mean2", "mean2", 700,2000); 
  RooRealVar sigma1("sigma1", "sigma1", 0, 3000); 
  RooRealVar sigma2("sigma2", "sigma2", 0, 3000); 
  RooGaussModel Gauss1("Gauss1","Gauss1",XMassSVFit,mean1,sigma1);
  RooGaussModel Gauss2("Gauss2","Gauss2",XMassSVFit,mean2,sigma2);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooRealVar nbkg2("nbkg2","signal fraction",100,0.,100000.); 
  RooAddPdf background("background", "gaussian plus exponential PDF",RooArgList(Gauss1, Gauss2), RooArgList(nbkg, nbkg2));
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void BifurGauss(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
		float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 0,3000); 
  RooRealVar sigma1("sigma1", "sigma1", 0, 400);
  RooRealVar sigma2("sigma2", "sigma2", 0, 700);
  RooBifurGauss bkg("bkg","bkg",XMassSVFit,mean,sigma1,sigma2);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void LogNormal(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	       RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
	       float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
	       float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 500,2000); 
  RooRealVar sigma("sigma", "sigma", 0, 3000);
  RooLognormal bkg("bkg","bkg",XMassSVFit,mean,sigma);
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}

void BreitWignerXPol0(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		      RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield,
		      float mean1_value, float mean2_value, float sigma1_value, float sigma2_value, float nsig1_value, float nsig2_value, 
		      float mean1_error, float mean2_error, float sigma1_error, float sigma2_error, float nsig1_error, float nsig2_error){
  RooRealVar mean("mean", "mean", 700,0,3000); 
  RooRealVar sigma("sigma", "sigma", 700, 0,1000); 
  RooRealVar c0("c0", "c0", 1, -10,10); 
  RooRealVar c1("c1", "c1", 1, -10,10); 
  RooRealVar c2("c2", "c2", 1, -10,10); 
  RooBreitWigner BreitWigner("BreitWigner","BreitWigner",XMassSVFit,mean,sigma);
  RooPolynomial Polynomial("Polynomial","Polynomial",XMassSVFit,RooArgList(c0));
  RooRealVar nbkg("nbkg","background fraction",100,0.,100000.); 
  RooFFTConvPdf bkg("bkg","bkg",XMassSVFit,BreitWigner,Polynomial);
  RooAddPdf background("background","background",bkg,nbkg);
  RooFitResult* r = background.fitTo(dataset, RooFit::Range(0,MAX), RooFit::Extended(kTRUE),RooFit::Save());
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2));  
  if(color==2)      background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kOrange),RooFit::Name("background")); 
  else if(color==1) background.plotOn(frame,RooFit::VisualizeError(*r,1,kFALSE),RooFit::FillColor(kGray),RooFit::Name("background")); 
  background.plotOn(frame,RooFit::LineColor(color),RooFit::Name("background"));
  dataset.plotOn(frame,RooFit::Name("dataset"),RooFit::Binning(XMassBin),RooFit::LineColor(color),RooFit::MarkerColor(color),RooFit::LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("background","dataset",nParam);
  RooArgSet* params = background.getVariables();
  cout<<endl;
  cout<<"BACKGROUND MODEL"<<endl;
  params->Print("V");
  cout<<"  "<<nParam+2<<") chi2/NDF = "<<chi2<<endl;
  cout<<endl;
  float mean_value=0; float sigma_value=0; float nbkg_value=0;
  float mean_error=0; float sigma_error=0; float nbkg_error=0;
  if(color==2) SaveWorkspace(Tree1, background, signal, channel, function, SignalMass, BkgYield, SigYield, mean1_value, mean2_value, sigma1_value, sigma2_value, nsig1_value, nsig2_value, mean1_error, mean2_error, sigma1_error, sigma2_error, nsig1_error, nsig2_error, mean_value, sigma_value, nbkg_value, mean_error, sigma_error, nbkg_error);
}















/*
//#include "TROOT.h"
//#include "TStyle.h"
//#include "TSystem.h"
//#include "TH1.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TCanvas.h"
//#include "TLegend.h"
//#include "TLatex.h"
//#include "TLine.h"
//#include "RooRealVar.h"
//#include "RooPlot.h"
//#include "RooAddPdf.h"
//#include "RooFormulaVar.h"
//#include "RooDataSet.h"
//#include "RooFitResult.h"
//#include "RooDataHist.h"
//#include "RooHist.h"
//#include "RooWorkspace.h"
//#include "RooGaussModel.h"
//#include "RooBreitWigner.h"
//#include "RooPolynomial.h"
//#include "RooFFTConvPdf.h"
//#include "RooLognormal.h"
//#include "RooBifurGauss.h"
//#include "RooExponential.h"
//#include "RooLandau.h"
//#include "RooKeysPdf.h"
//#include "RooVoigtian.h"
//#include "RooArgusBG.h"
//#include "RooCBShape.h"

void Fitting(TCanvas* c, float MAX, TH1F *histo_PRE, TH1F *histo_MC, TString FUNCTION, float XMassWidth, float XMassBin,
	     float XMassMin, float XMassMax, char *channel, char *function, int SignalMass, float BkgYield, bool save, bool pull, bool resid);
void SavePlot(TCanvas* c, TH1F *histo_PRE, TH1F *histo_MC, RooPlot* frame1, RooPlot* frame2, RooRealVar XMassSVFit, 
	      TString CHANNEL, TString FUNCTION, bool save, bool pull, bool resid, float MAX);
void SaveWorkspace(TTree *Tree1, RooAddPdf background, RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void SaveDatacard(char *channel, char *function, int SignalMass, int data_obs_entries, float BkgYield, float SigYield);
void BackgroundEstimation(char *channel, float XMassWidth, float XMassBin, float XMassMin, float XMassMax, TH1F *histo_PRE, TH1F *histo_MC, float & BkgYield);
void Gaussian(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	      RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void Landau(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	    RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void CrystalBall(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		 RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void Argus(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	   RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void BreitWigner(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		 RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void Voigtian(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	      RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
//void KeysPdf(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1,
//	     RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void GaussXLandau(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		  RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void GaussXGauss(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		 RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void GaussXExpo(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void GaussXPoly(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void GaussPlusGauss(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		    RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void BifurGauss(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void LogNormal(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
	       RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
void BreitWignerXPol0(RooPlot* frame, RooDataHist dataset, RooRealVar XMassSVFit, float XMassBin, float MAX, int color, int nParam, TTree *Tree1, 
		      RooAddPdf signal, char *channel, char *function, int SignalMass, float BkgYield, float SigYield);
*/
