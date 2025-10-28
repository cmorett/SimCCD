void B02()
{
  //------------- Style --------------
  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  //gStyle->SetOptStat(0);
  //----------------------------------
  
  //TFile *infile = 
  //new TFile("/home/alexis/GEANT4/work/bryanSim/B02/output/muonsconcrete_new_001.root");
  //TTree *HitTree = (TTree*)infile->Get("B02Hit");
  //TTree *EvtTree = (TTree*)infile->Get("B02Evt");

  TChain *HitTree = new TChain("B02Hits");
  TChain *EvtTree = new TChain("B02Evts");

  int nFiles = 10;
  string strID;
  for (int i=1;i<nFiles+1;i++){
    TString s(Form("%d",i));
    if (i<10) strID = "00" + s;
    else if(i<100) strID = "0" + s;
    else if(i<1000) strID = "" + s;

    HitTree->Add(Form("/home/cristianm/CCM/B02-files/output/muons_%s.root",strID.c_str()));
    EvtTree->Add(Form("/home/cristianm/CCM/B02-files/output/muons_%s.root",strID.c_str()));
  }

  double diam = 258.0;
  double heig = 225.0;

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

  TH1F *hEd = new TH1F("hEd","",100,0,1200);
  hEd->GetXaxis()->SetTitle("E dep (MeV)");
  hEd->GetYaxis()->SetTitle("events");

  TH1F *hEv = new TH1F("hEv","",100,0,1200);
  hEv->GetXaxis()->SetTitle("E vis  (MeV)");
  hEv->GetYaxis()->SetTitle("events");

  TH1F *hEg = new TH1F("hEg","",100,0,1200);
  hEg->GetXaxis()->SetTitle("E w/Res (MeV)");
  hEg->GetYaxis()->SetTitle("events");

  TH1F *hEa = new TH1F("hEa","",100,0,12000);
  hEa->GetXaxis()->SetTitle("E reco (PE)");
  hEa->GetYaxis()->SetTitle("events");

  HitTree->Draw("ZhitBar:YhitBar:XhitBar>>hxyz","evtId<10");
  //HitTree->Draw("YhitBar:XhitBar>>hxy","");
  //HitTree->Draw("ZhitBar:XhitBar>>hxz","");
  EvtTree->Draw("EevtBar >> hEd","EevtBar>2"); 
  EvtTree->Draw("WevtBar >> hEv","EevtBar>2"); 
  EvtTree->Draw("GevtBar >> hEg","EevtBar>2"); 
  //EvtTree->Draw("ADCevtBar >> hEa","EevtBar>2"); 
  EvtTree->Draw("0.0+GevtBar*10.0/(1+GevtBar/50000.) >> hEa","0.0+GevtBar*1.0/(1+GevtBar/5000.)>2");  

  //---------------------
  //----Drawing section 
  TH3F *frame3D = new TH3F("frame3D","",400,-150,150,400,-150,150,400,-150,150);
  frame3D->GetXaxis()->SetTitle("X (cm)");
  frame3D->GetYaxis()->SetTitle("Y (cm)");
  frame3D->GetZaxis()->SetTitle("Z (cm)");

  TCanvas *Canv = new TCanvas("Canv","",2*300,2*300);
  Canv->Divide(2,2);
  
  Canv->cd(1);
  frame3D->Draw();
  hxyz->Draw("same");
  Canv->cd(2);
  hEd->Draw();
  hEv->SetLineStyle(2);hEv->Draw("same");
  Canv->cd(3);
  hEg->Draw();
  Canv->cd(4);
  hEa->Draw();
  
  Canv->Print("plot.pdf");
}
