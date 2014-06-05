{
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat("rme");
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5); //title X location 
  gStyle->SetTitleY(0.96); //title Y location 
  gStyle->SetPaintTextFormat(".2f");

  using namespace std;

  bool save = true;
  char channel[500]; 
  char efficiency1[500];  
  for(int i=0; i<8; i++){
    if(i==0){ sprintf(channel,   "_mm"); sprintf(efficiency1,   "Muo1%s",channel);}
    if(i==1){ sprintf(channel,   "_mm"); sprintf(efficiency1,   "Muo2%s",channel);}
    if(i==2){ sprintf(channel,   "_mt"); sprintf(efficiency1,   "Muo%s",channel);}
    if(i==3){ sprintf(channel,   "_em"); sprintf(efficiency1,   "Muo%s",channel);}
    if(i==4){ sprintf(channel,   "_ee"); sprintf(efficiency1,   "Ele1%s",channel);}
    if(i==5){ sprintf(channel,   "_ee"); sprintf(efficiency1,   "Ele2%s",channel);}
    if(i==6){ sprintf(channel,   "_et"); sprintf(efficiency1,   "Ele%s",channel);}
    if(i==7){ sprintf(channel,   "_em"); sprintf(efficiency1,   "Ele%s",channel);}

    char demo       [500];  sprintf(demo,      "demo/Tree"); 
    char openTree   [500];  sprintf(openTree,  "%s",demo); 
    TString EFFICIENCY = efficiency1; 

    TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
    TFile *file1 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/ZH1000.root");
    TFile *file2 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/ZH1500.root");
    TFile *file3 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/ZH2000.root");
    TFile *file4 = TFile::Open("../../RISULTATI/efficiency_220514_CorrPFIso/ZH2500.root");

    TTree *Tree1;  Tree1 = (TTree*)file1->Get(openTree); 
    TTree *Tree2;	 Tree2 = (TTree*)file2->Get(openTree); 
    TTree *Tree3;	 Tree3 = (TTree*)file3->Get(openTree); 
    TTree *Tree4;	 Tree4 = (TTree*)file4->Get(openTree); 
  
    char input1[500];
    sprintf(input1,  "gen%s>0.5",efficiency1);
  
    char variable1[500]; sprintf(variable1, "gen%s_Pt>>h1(25,0,500)",efficiency1); Tree1->Draw(variable1, input1, "E");
    char variable2[500]; sprintf(variable2, "gen%s_Pt>>h2(25,0,500)",efficiency1); Tree2->Draw(variable2, input1, "Esame");
    char variable3[500]; sprintf(variable3, "gen%s_Pt>>h3(25,0,500)",efficiency1); Tree3->Draw(variable3, input1, "Esame");
    char variable4[500]; sprintf(variable4, "gen%s_Pt>>h4(25,0,500)",efficiency1); Tree4->Draw(variable4, input1, "Esame");

    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h3->SetLineWidth(2);
    h4->SetLineWidth(2);
  
    h1->SetLineColor(1);
    h2->SetLineColor(2);
    h3->SetLineColor(kGreen+2);
    h4->SetLineColor(4);
  
    h1->SetMarkerColor(1);
    h2->SetMarkerColor(2);
    h3->SetMarkerColor(kGreen+2);
    h4->SetMarkerColor(4);

    h1->SetMarkerStyle(20); 
    h2->SetMarkerStyle(20); 
    h3->SetMarkerStyle(20); 
    h4->SetMarkerStyle(20); 
  
    h1->SetMarkerSize(1.3); 
    h2->SetMarkerSize(1.3); 
    h3->SetMarkerSize(1.3); 
    h4->SetMarkerSize(1.3); 

    float N1=1./h1->Integral();
    float N2=1./h2->Integral();
    float N3=1./h3->Integral();
    float N4=1./h4->Integral();

    h1->Scale(N1);
    h2->Scale(N2);
    h3->Scale(N3);
    h4->Scale(N4);

    h1->Draw("histo");
    c1->Update();
    TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    ps1->SetY1NDC(0.60); ps1->SetY2NDC(0.74); ps1->SetTextColor(1);
    ps1->SetX1NDC(0.70); ps1->SetX2NDC(0.89); 
    h4->Draw("histo sames");
    c1->Update();
    TPaveStats *ps4 = (TPaveStats*)h4->GetListOfFunctions()->FindObject("stats");
    ps4->SetY1NDC(0.45); ps4->SetY2NDC(0.59); ps4->SetTextColor(4);
    ps4->SetX1NDC(0.70); ps4->SetX2NDC(0.89); 

    //h1->Draw("histo");
    //h2->Draw("histo same");
    //h3->Draw("histo same");
    //h4->Draw("histo same");

    float min = h1->GetXaxis()->GetXmin();
    float max = h1->GetXaxis()->GetXmax();
    float bin = h1->GetXaxis()->GetNbins();
    h1->GetYaxis()->SetTitleSize(0.045);
    h1->GetXaxis()->SetTitleSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetTitleOffset(1.1); 
    h1->SetTitle("");
    h1->GetYaxis()->SetTitle(TString("Normalized Events / ")+TString::Format("%.2f",(max - min)/bin));
    h1->GetXaxis()->SetTitle("gen lepton pt [GeV]");
  
    TLegend *pl2 = new TLegend(0.60,0.77,0.89,0.89);
    pl2->SetTextSize(0.04); 
    pl2->SetFillColor(0);
    TLegendEntry *ple2 = pl2->AddEntry(h1, "M(X) = 1.0 TeV",  "LP");
    //ple2 = pl2->AddEntry(h2, "M(X) = 1.5 TeV",  "LP");
    //ple2 = pl2->AddEntry(h3, "M(X) = 2.0 TeV",  "LP");
    ple2 = pl2->AddEntry(h4, "M(X) = 2.5 TeV",  "LP");
    pl2->Draw();
  
    if(save){
      c1->SaveAs("pt_"+EFFICIENCY+".pdf");
      cout<<"Saving pt_"<<EFFICIENCY<<".pdf"<<endl;
    }
  }
}
