
void loopB02()
{
  //------------- Style --------------
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //----------------------------------
  
  TChain *HitTree = new TChain("B02Hits");
  TChain *EvtTree = new TChain("B02Evts");
  
  // adding new Tchain to read data without shielding
  TChain *HitTreewithoutshielding = new TChain("B02Hits");
  TChain *EvtTreewithoutshielding = new TChain("B02Evts");

  //----------------------------------
  // Input Monte Carlo simulation files
  //----------------------------------

  int nFiles = 200;
  string strID;
  for (int i=1;i<=nFiles;i++){
    TString s(Form("%d",i));
    if (i<10) strID = "00" + s;
    else if(i<100) strID = "0" + s;
    else if(i<1000) strID = "" + s;

   
   //HitTree->Add(Form("/home/cristianm/CCM/B02-files/output/muons_%s.root",strID.c_str())); // geometry's alexis
   //EvtTree->Add(Form("/home/cristianm/CCM/B02-files/output/muons_%s.root",strID.c_str())); // geometry's alexis
  
   HitTree->Add(Form("/home/cristianm/CCMOCT23/B02-files/outputwithshielding/muons_%s.root",strID.c_str())); // Complete shielding 
   EvtTree->Add(Form("/home/cristianm/CCMOCT23/B02-files/outputwithshielding/muons_%s.root",strID.c_str())); // Complete shielding 
  
   HitTreewithoutshielding->Add(Form("/home/cristianm/CCMOCT23/B02-files/outputwithoutshielding/muons_%s.root",strID.c_str())); // without shielding 
   EvtTreewithoutshielding->Add(Form("/home/cristianm/CCMOCT23/B02-files/outputwithoutshielding/muons_%s.root",strID.c_str())); // without shielding 
  
  }

  //-----------------------------------
  // Position histograms
  //-----------------------------------

  TH3F *hxyz = new TH3F("hxyz","",200,-150,150, 200,-150,150, 200,-150,150);
  hxyz->GetXaxis()->SetTitle("X (cm)");
  hxyz->GetYaxis()->SetTitle("Y (cm)");
  hxyz->GetZaxis()->SetTitle("Z (cm)");

  TH2F *hxy = new TH2F("hxy","",200,-150,150,200,-150,150);
  hxy->GetXaxis()->SetTitle("X (cm)");
  hxy->GetYaxis()->SetTitle("Y (cm)");

  TH2F *hxz = new TH2F("hxz","",200,-150,150,200,-150,150);
  hxz->GetXaxis()->SetTitle("X (cm)");
  hxz->GetYaxis()->SetTitle("Z (cm)");

  //-----------------------------------
  // Energy histograms
  //-----------------------------------

  double  loEd = 0.;
  double hiEd = 700;
  int     nbEd = 200;
  double loADC = 0.;
  double hiADC = 40000.;
  int    nbADC = 40000;

  TH1F *hEdep_mu = new TH1F("hEdep_mu","",nbEd,loEd,hiEd);
  hEdep_mu->SetLineColor(4);
  hEdep_mu->GetXaxis()->SetTitle("Edep (MeV)");
  hEdep_mu->GetYaxis()->SetTitle("events");
  TH1F *hEvis_mu = new TH1F("hEvis_mu","",nbEd,loEd,hiEd);
  hEvis_mu->SetLineColor(4);
  hEvis_mu->GetXaxis()->SetTitle("Evis (MeV)");
  hEvis_mu->GetYaxis()->SetTitle("events");
  TH1F *hEres_mu = new TH1F("hEres_mu","",nbEd,loEd,hiEd);
  hEres_mu->SetLineColor(4);
  hEres_mu->GetXaxis()->SetTitle("Eres (MeV)");
  hEres_mu->GetYaxis()->SetTitle("events");
  TH1F *hEadc_mu = new TH1F("hEadc_mu","",nbADC,loADC,hiADC);
  hEadc_mu->SetLineColor(4);
  hEadc_mu->GetXaxis()->SetTitle("Erec (PE)");
  hEadc_mu->GetYaxis()->SetTitle("events");

  TH1F *hEdep_bg = new TH1F("hEdep_bg","",nbEd,loEd,hiEd);
  hEdep_bg->SetLineColor(2);
  hEdep_bg->GetXaxis()->SetTitle("Edep (MeV)");
  hEdep_bg->GetYaxis()->SetTitle("events");
  TH1F *hEvis_bg = new TH1F("hEvis_bg","",nbEd,loEd,hiEd);
  hEvis_bg->SetLineColor(2);
  hEvis_bg->GetXaxis()->SetTitle("Evis (MeV)");
  hEvis_bg->GetYaxis()->SetTitle("events");
  TH1F *hEres_bg = new TH1F("hEres_bg","",nbEd,loEd,hiEd);
  hEres_bg->SetLineColor(2);
  hEres_bg->GetXaxis()->SetTitle("Eres (MeV)");
  hEres_bg->GetYaxis()->SetTitle("events");
  TH1F *hEadc_bg = new TH1F("hEadc_bg","",nbADC,loADC,hiADC);
  hEadc_bg->SetLineColor(2);
  hEadc_bg->GetXaxis()->SetTitle("Erec (PE)");
  hEadc_bg->GetYaxis()->SetTitle("events");


  TH1F *hEdep = new TH1F("hEdep","",nbEd,loEd,hiEd); // deposited energy histogram 
  hEdep->SetLineColor(3);
  hEdep->GetXaxis()->SetTitle("Edep (MeV)");
  hEdep->GetYaxis()->SetTitle("events");
  TH1F *hEvis = new TH1F("hEvis","",nbEd,loEd,hiEd); // EM saturation histogram
  hEvis->SetLineColor(3);
  hEvis->GetXaxis()->SetTitle("Evis (MeV)");
  hEvis->GetYaxis()->SetTitle("events");
  TH1F *hEres = new TH1F("hEres","",nbEd,loEd,hiEd);  // Gaussina resolution histogram 
  hEres->SetLineColor(3);
  hEres->GetXaxis()->SetTitle("Eres (MeV)");
  hEres->GetYaxis()->SetTitle("events");
  TH1F *hEadc = new TH1F("hEadc","",nbADC,loADC,hiADC); // ADC unit energy histogram 
  hEadc->SetLineColor(3);
  hEadc->GetXaxis()->SetTitle("Erec (PE)");
  hEadc->GetYaxis()->SetTitle("events");

  int NBmu = 200;
  double maxEd = 700.;
  double muonsbw = maxEd/NBmu;
  TH1F *muons = new TH1F("muons","",NBmu,0,maxEd);
  muons->GetXaxis()->SetTitle("Energy (MeV)");
  muons->GetYaxis()->SetTitle("events");
  
  //////////////////////////////////////////////////////////
 
  // added to data without shielding
  
  TH1F *muonswithoutshielding = new TH1F("muonswithoutshielding", "", NBmu, 0, maxEd);
  muonswithoutshielding->GetXaxis()->SetTitle("Energy (MeV)");
  muonswithoutshielding->GetYaxis()->SetTitle("events");

  //////////////////////////////////////////////////

  //-----------------------------------
  // Energy reponse parameters
  //-----------------------------------
	
  double E00 = 259.0;   // MeV
  double f =  0.08;   // frac
  double a0 = 0.0;    // PE
  double a1 = 70;    // PE/MeV
  double a2 = 0E-04;  // MeV^{-1}

  //-----------------------------------

  TRandom3 *r = new TRandom3();

  double  EevtBar;
  double  WevtBar;
  EvtTree->SetBranchAddress("EevtBar",&EevtBar);
  EvtTree->SetBranchAddress("WevtBar",&WevtBar);
  
  double  XhitBar, YhitBar, ZhitBar;
  HitTree->SetBranchAddress("XhitBar",&XhitBar);
  HitTree->SetBranchAddress("YhitBar",&YhitBar);
  HitTree->SetBranchAddress("ZhitBar",&ZhitBar);

   /////////////////////////////////////////////////////// 

  // added to data withoutshielding
  
  double  EevtBarwithoutshielding;
  double  WevtBarwithoutshielding;
  
  EvtTreewithoutshielding->SetBranchAddress("EevtBar",&EevtBarwithoutshielding);
  EvtTreewithoutshielding->SetBranchAddress("WevtBar",&WevtBarwithoutshielding);
  
  double  XhitBarwithoutshielding, YhitBarwithoutshielding, ZhitBarwithoutshielding;
  HitTreewithoutshielding->SetBranchAddress("XhitBar",&XhitBarwithoutshielding);
  HitTreewithoutshielding->SetBranchAddress("YhitBar",&YhitBarwithoutshielding);
  HitTreewithoutshielding->SetBranchAddress("ZhitBar",&ZhitBarwithoutshielding);

  ///////////////////////////////////////////////////

  double  GevtBar; //with Gaussian resolution function
  double  ADCevtBar; //in ADC with no-linearity

//----------------------------
//--Loop over EvtTree tree 
//----------------------------

  int nevent = EvtTree->GetEntries();
  cout<<"EvtTree entries: "<< nevent <<endl;
  for (int i=0;i<nevent;i++) {
     // read event
     EvtTree->GetEvent(i);
     double res = f;                           //
     double sig = res*sqrt(E00*WevtBar);     //
     double xg = r->Gaus(0, sig);              //
     GevtBar = xg+WevtBar;                     //
     double p1 = a0;                           //
     double p2 = a1;                           //
     double p3 = a2;                           //
     ADCevtBar = p1+p2*GevtBar/(1+p3*GevtBar); //
 
     if (EevtBar>0) hEdep_mu->Fill(EevtBar);
     if (WevtBar>0) hEvis_mu->Fill(WevtBar);
     if (WevtBar>0) muons->Fill(WevtBar); //** muons template **//

     if (GevtBar>0) hEres_mu->Fill(GevtBar);
     if (ADCevtBar>0) hEadc_mu->Fill(ADCevtBar);

  } //for i

//////////////////////////////////////////////////
// added to data without shielding 

int neventwithoutshielding = EvtTreewithoutshielding->GetEntries();
  cout<<"EvtTreewithoutshielding entries: "<< neventwithoutshielding <<endl;
  for (int i=0;i<neventwithoutshielding;i++) {
     // read event
     EvtTreewithoutshielding->GetEvent(i);
     //double res = f;                           //
     //double sig = res*sqrt(E00*WevtBar);     //
     //double xg = r->Gaus(0, sig);              //
     //GevtBar = xg+WevtBar;                     //
     //double p1 = a0;                           //
     //double p2 = a1;                           //
     //double p3 = a2;                           //
     //ADCevtBar = p1+p2*GevtBar/(1+p3*GevtBar); //
 
     //if (EevtBar>0) hEdep_mu->Fill(EevtBar);
     //if (WevtBar>0) hEvis_mu->Fill(WevtBar);
    if (WevtBarwithoutshielding>0) muonswithoutshielding->Fill(WevtBarwithoutshielding); //** muons template **//

     //if (GevtBar>0) hEres_mu->Fill(GevtBar);
     //if (ADCevtBar>0) hEadc_mu->Fill(ADCevtBar);

  }

///////////////////////////////////


//----------------------------
//-- Add background
//----------------------------
 
// Histograma de Ed simulado (200 bins de 0 a maxEx MeV) 
ifstream inpf;
inpf.open("background.dat");
double pv;
int nbins = 200;
double dE = maxEd/nbins;
double Emax = (nbins)*dE;

// Histograma de energías depositadas Ed
TH1F *hbackground = new TH1F("hbackground","",nbins,0,Emax);
hbackground->GetXaxis()->SetTitle("Energy (MeV)");
for (int i=1; i<=nbins; i++) {
    inpf >> pv;    // Reading input-file line
    hbackground->SetBinContent(i,pv);
} //for i
cout << "Integral hbackground: " << hbackground->Integral(1,nbins,"width") << "." << endl;


//int nbg = 17000;
int nbg = 17000;

for (int i=0;i<nbg;i++) {
     EevtBar = hbackground->GetRandom();
     WevtBar = hbackground->GetRandom();
     //WevtBar = EevtBar; // no Birks for now

     double res = f;                           //
     double sig = res*sqrt(E00*WevtBar);     //
     double xg = r->Gaus(0, sig);              //
     GevtBar = WevtBar+xg;                     //
     double p1 = a0;                           //
     double p2 = a1;                           //
     double p3 = a2;                           //
     ADCevtBar = p1+p2*GevtBar/(1+p3*GevtBar); //

     if (EevtBar>0) hEdep_bg->Fill(EevtBar);
     if (WevtBar>0) hEvis_bg->Fill(WevtBar);
     if (GevtBar>0) hEres_bg->Fill(GevtBar);
     if (ADCevtBar>0) hEadc_bg->Fill(ADCevtBar);

}

hEdep->Add(hEdep_mu); hEdep->Add(hEdep_bg); // deposited energy + background histogram 
hEvis->Add(hEvis_mu); hEvis->Add(hEvis_bg);  // Em saturation + background histogram 
hEres->Add(hEres_mu); hEres->Add(hEres_bg); // Gaussian resolution + background histogram 
hEadc->Add(hEadc_mu); hEadc->Add(hEadc_bg); // ADC units energy + background histogram

//----------------------------
//--Loop over HitTree tree 
//----------------------------

  nevent = HitTree->GetEntries();
  cout<<"HitTree entries: "<< nevent <<endl;
  //for (int i=0;i<nevent;i++) {
  for (int i=0;i<30000;i++) {

     // read event
     HitTree->GetEvent(i);
     hxyz->Fill(XhitBar,YhitBar,ZhitBar); // for fit

  } // for i

//////////////////////////////////////////
// added to data withoutshielding 

  neventwithoutshielding = HitTreewithoutshielding->GetEntries();
  cout<<"HitTreewithoutshielding entries: "<< neventwithoutshielding <<endl;
  //for (int i=0;i<nevent;i++) {
  for (int i=0;i<30000;i++) {

     // read event
     HitTreewithoutshielding->GetEvent(i);
     hxyz->Fill(XhitBarwithoutshielding,YhitBarwithoutshielding,ZhitBarwithoutshielding); // for fit

  }

//////////////////////////////////////

//----------------------------
//-- muons template data file 
//----------------------------

  ofstream myfile;
  myfile.open ("muons.dat");
  for (int i=0;i<NBmu; i++) 
  myfile << muons->GetBinContent(i+1)/(muons->Integral(1,NBmu)*muonsbw) << "\n";
  myfile.close();

  cout << "Integral muons: " << muons->Integral() << " entries." << endl;
  cout << "Integral muons (+ovflw): " << muons->Integral(0,NBmu+1) << " entries." << endl;

  //- Find muon peak
  int peakBin = 50;
  for (int i=30;i<NBmu;i++){
   if (muons->GetBinContent(i+1) > muons->GetBinContent(peakBin))
     peakBin = i+1;
  } //for i
  printf("Peak: %d  \n", peakBin);
  printf("Peak: %3.3f MeV \n", peakBin*maxEd/NBmu);
  double Epeak = peakBin*maxEd/NBmu;
  
  TLatex *latt = new TLatex;
  //latt->SetNDC();

  TCanvas *canvmu = new TCanvas("canvmu","",900,700);
  canvmu->cd(1);
  gPad->SetLogy(1);
  muons->Draw();
  latt->DrawLatex(Epeak, 1.05*muons->GetBinContent(peakBin), Form("%3.1f MeV",Epeak));
  canvmu->Print("muons.pdf");

  TFile *foutmuons = new TFile("muons.root", "recreate");
  muons->Write();
  foutmuons->Close();


///////////////////////////////////////////////////////  
  //  added to data without shielding 
  
  ofstream myfilewithoutshielding;
  myfilewithoutshielding.open ("muonswithoutshielding.dat");
  for (int i=0;i<NBmu; i++) 
  myfilewithoutshielding << muonswithoutshielding->GetBinContent(i+1)/(muonswithoutshielding->Integral(1,NBmu)*muonsbw) << "\n";
  myfilewithoutshielding.close();

  cout << "Integral muonswithoutshielding: " << muonswithoutshielding->Integral() << " entries." << endl;
  cout << "Integral muonswithoutshielding (+ovflw): " << muonswithoutshielding->Integral(0,NBmu+1) << " entries." << endl;
  
  TCanvas *canvmush = new TCanvas("canvmush","",900,700);
  canvmush->cd(1);
  gPad->SetLogy(1);
  muonswithoutshielding->Draw();
  //latt->DrawLatex(Epeak, 1.05*muons->GetBinContent(peakBin), Form("%3.1f MeV",Epeak));
  canvmush->Print("muonswithoutshielding.pdf");
  
  TFile *foutmuonswithoutshielding = new TFile("muonswithoutshielding.root", "recreate");
  muonswithoutshielding->Write();
  foutmuonswithoutshielding->Close();

////////////////////////////////////////////////////////


//----------------------------
//-- nonlinearity function def
//----------------------------

  TF1 *fMeVtoPE = new TF1("fMeVtoPE","[0]+[1]*x/(1+[2]*x)",0,hiEd);
  fMeVtoPE->SetParameters(a0,a1,a2);
  fMeVtoPE->GetXaxis()->SetTitle("E_{vis} (MeV)");
  fMeVtoPE->GetYaxis()->SetTitle("E_{rec} (PE)");
  fMeVtoPE->GetYaxis()->SetTitleOffset(1.4);

  //---------------------------
  //- Drawing section 
  //---------------------------

  TLatex *lat = new TLatex;
  lat->SetNDC();

  TH3F *frame3D = new TH3F("frame3D","",200,-150,150, 200,-150,150, 200,-150,150);
  frame3D->GetXaxis()->SetTitle("X (cm)");
  frame3D->GetYaxis()->SetTitle("Y (cm)");
  frame3D->GetZaxis()->SetTitle("Z (cm)");


  TH1F *frameEdep = new TH1F("frameEdep","",200,0,maxEd);
  
  TCanvas *canv0 = new TCanvas("canv0","",1*400,1*300);
  fMeVtoPE->Draw();
  double xx = 500;
  double yy = fMeVtoPE->Eval(xx);
  TLine *l0 = new TLine(0,0,hiEd,0); l0->Draw("same");
  TLine *lh = new TLine(0,yy,xx,yy); lh->Draw("same");
  TLine *lv = new TLine(xx,0,xx,yy); lv->Draw("same");
  lat->SetTextFont(42);
  lat->SetTextSize(0.035);
  lat->DrawLatex(0.6,0.50,"E_{PE} = a_{0} + #frac{a_{1} E_{vis}}{1 + a_{2} E_{vis}}");
  lat->DrawLatex(0.6,0.40,Form("a_{0} = %3.1f PE",a0));
  lat->DrawLatex(0.6,0.35,Form("a_{1} = %3.1f PE/MeV",a1));
  lat->DrawLatex(0.6,0.30,Form("a_{2} = %3.1e MeV^{-1}",a2));
  canv0->Print("MeVtoPE.pdf");

  TCanvas *Canv = new TCanvas("Canv","",2*300,2*300);
  Canv->Divide(2,2);
  
  Canv->cd(1);
  frame3D->Draw();
  hxyz->Draw("same");
  //TGeoTube *myCyl = new TGeoTube(0., 250, 250);
  //myCyl->Draw("same");
  Canv->cd(2);
  hEdep->SetLineStyle(3);
  hEdep->SetLineColor(4);
  gPad->SetLogy(1);
  hEdep->Draw();
  hEvis->Draw("same");
  Canv->cd(3);
  gPad->SetLogy(1);
  hEres->Draw();
  hEres_mu->Draw("same");
  hEres_bg->Draw("same");
  Canv->cd(4);
  gPad->SetLogy(0);
  gPad->SetGrid();
  int mbin = (int) 17500*nbADC/(hiADC-loADC);
  hEadc->SetMaximum(2.5*hEadc->GetBinContent(mbin));
  hEadc->Draw();
  hEadc_mu->Draw("same");
  hEadc_bg->Draw("same");
  hEadc->Draw("same");
  
  Canv->Print("loopPlot.pdf");
 
 
  
 
  TFile *outfile = new TFile("myHisto.root", "recreate");
  hEadc->Write();
  outfile->Close();
 
}
