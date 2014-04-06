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
  using namespace RooFit;
  
  bool save=false;
  bool pull=true;
  bool resid=false;
  int SignalMass = 1500;
  
  char channel    [500]; sprintf(channel,  "EleMuo"); 
  char demo       [500]; sprintf(demo,       "demo/Tree"); 
  char openTree   [500]; sprintf(openTree,   "%s%s",demo,channel); 
  TString CHANNEL = channel;


  RooRealVar XMassSVFit("XMassSVFit","XMassSVFit",SignalMass-750,SignalMass+750); 
  if(SignalMass==1000)      {
    TFile *fileSig = TFile::Open("../RISULTATI/SB2_050214/ZH1000_FL.root");  
    TTree *TreeSig=(TTree*)fileSig->Get(openTree); 
    double WeightSig = 1.54*19.7/17958;
    RooFormulaVar wFunc("w","w","1.54*19.7/17958",XMassSVFit);
    double MEAN1 = 1000; double MEAN2 = 700;
  } else if(SignalMass==1500) {
    TFile *fileSig = TFile::Open("../RISULTATI/SB2_050214/ZH1500_FL.root");  
    TTree *TreeSig=(TTree*)fileSig->Get(openTree); 
    double WeightSig = 0.22*19.7/14099;
    RooFormulaVar wFunc("w","w","0.22*19.7/14099",XMassSVFit);
    double MEAN1 = 1500; double MEAN2 = 1200;
  } else if(SignalMass==2000) {
    TFile *fileSig = TFile::Open("../RISULTATI/SB2_050214/ZH2000_FL.root");  
    TTree *TreeSig=(TTree*)fileSig->Get(openTree); 
    double WeightSig = 0.03*19.7/13080;
    RooFormulaVar wFunc("w","w","0.03*19.7/13080",XMassSVFit);
    double MEAN1 = 2000; double MEAN2 = 1700;
  } else if(SignalMass==2500) {
    TFile *fileSig = TFile::Open("../RISULTATI/SB2_050214/ZH2500_FL.root");  
    TTree *TreeSig=(TTree*)fileSig->Get(openTree); 
    double WeightSig = 0.01*19.7/13091;
    RooFormulaVar wFunc("w","w","0.01*19.7/13091",XMassSVFit);
    double MEAN1 = 2500; double MEAN2 = 2200;
  }

  TTree *SignalTree = TreeSig->CopyTree("met>50 && uncorrmet>30 && trigger==1");
  RooRealVar PUWeight("PUWeight","PUWeight",0,100);
  RooArgList VarSet(XMassSVFit,PUWeight);
  RooPlot* frame = XMassSVFit.frame();
  RooDataSet* signalSampleUnw = new RooDataSet("signalSampleUnw","signalSampleUnw",SignalTree,VarSet,0,"PUWeight");
  RooRealVar* w = (RooRealVar*) signalSampleUnw->addColumn(wFunc);
  RooDataSet signalSample("signalSample","signalSample",signalSampleUnw,*signalSampleUnw->get(),0,"w") ;

  RooRealVar mean1("mean1", "mean1", MEAN1, SignalMass-100,SignalMass+100); 
  RooRealVar mean2("mean2", "mean2", MEAN2, SignalMass-300,SignalMass); 
  RooRealVar sigma1("sigma1", "sigma1", 40, 80); 
  RooRealVar sigma2("sigma2", "sigma2", 40, 150); 
  RooGaussModel Gauss1("Gauss1","Gauss1",XMassSVFit,mean1,sigma1);
  RooGaussModel Gauss2("Gauss2","Gauss2",XMassSVFit,mean2,sigma2);
  RooRealVar fsig("fsig","fsig",0.,1.); 
  RooAddPdf signal("signal","signal",RooArgList(Gauss1, Gauss2), fsig);

  RooFitResult* r = signal.fitTo(signalSample, RooFit::Range(SignalMass-750,SignalMass+750), RooFit::Extended(kFALSE), RooFit::Save());

  signalSample->plotOn(frame,Name("signalSample"),Binning(50),LineColor(2),MarkerColor(2),LineWidth(2));  
  //signal.plotOn(frame,VisualizeError(*r,1,kFALSE),FillColor(kOrange),Name("signal")); 
  signal.plotOn(frame,LineColor(2),Name("signal"));
  signalSample->plotOn(frame,Name("signalSample"),Binning(50),LineColor(2),MarkerColor(2),LineWidth(2)); 
  Double_t chi2 = frame->chiSquare("signal","signalSample",5);
  RooArgSet* params = signal->getVariables();
  frame->Draw();
  frame->SetMinimum(0);
  cout<<endl;
  params->Print("V");
  cout<<"  X) chi2/NDF = "<<chi2<<endl;
  cout<<endl;

}
