void Test_LvsThet(){
// TFile *file = new TFile("./root_files/muons_2M_vacuum_file.root");
// TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_2.root");
TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_nHG_0.root");
TTree *tree = (TTree*) file->Get("B02Evts");


TFile *file0 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_SIGMA_1.0_C_0.root");
TTree *tree0 = (TTree*) file0->Get("tree");


int NBL = 150;
double llo = 0.;
double lhi = 0.9;
int NBT = 150;
double tlo = 0.;
double thi = TMath::Pi()/2.0;
// TH1F *L = new TH1F("L", "Distance Distribution (CCD size: 0.9x0.6x0.0725 cm)", NB, tlow, thi);
// L->GetXaxis()->SetTitle("Distance (cm)");
TH2F *histLT = new TH2F("histLT", "Lvs#theta (Geant4)", NBL, llo, lhi, NBT, tlo, thi);
histLT->GetXaxis()->SetTitle(" Distance (cm)");
histLT->GetYaxis()->SetTitle("#theta (rad)");
histLT->SetStats(0);

TH2F *histLTpp = new TH2F("histLTpp", "Lvs#theta (SimPP)", NBL, llo, lhi, NBT, tlo, thi);
histLTpp->GetXaxis()->SetTitle(" Distance (cm)");
histLTpp->GetYaxis()->SetTitle("#theta (rad)");
histLTpp->SetStats(0);

// Fil histograms
double pl, pt, pedep, plpp, ptpp;
int nH;
tree->SetBranchAddress("LengthMuLAr", &pl);
tree->SetBranchAddress("thetaPri", &pt);
tree->SetBranchAddress("nHitBar", &nH);
tree->SetBranchAddress("EevtBar", &pedep);

tree0->SetBranchAddress("l", &plpp);
tree0->SetBranchAddress("thet", &ptpp);

int nentries = tree->GetEntries();
   for (int i = 0; i < nentries; i++) {
      tree->GetEntry(i);
    //   if (pl > 0 & nH >0){
    //   if (nH >0){
    //   if (pl > 0){
      if (pedep > 0){
        histLT->Fill(pl, pt);
      }
   }

nentries = tree0->GetEntries();
for (int i = 0; i < nentries; i++) {
    tree0->GetEntry(i);
    if (plpp > 0){
        histLTpp->Fill(plpp, ptpp);
    }
}

// tree->Draw("LengthMuLAr>>L", "LengthMuLAr>0 && thetaPri>25*TMath::Pi()/180");


TLine *line = new TLine(0.0725,0,0.0725,1.5);
line->SetLineStyle(2);
line->SetLineWidth(2);
line->SetLineColor(2);

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
histLT->Draw("colZ");
line->Draw("same");

TLegend *leg = new TLegend(0.5, 0.7, 0.8, 0.8);
leg->AddEntry(line, "Grosor de CCD: 0.0725 cm", "l");
leg->Draw();

canv->cd(2);
histLTpp->Draw("colZ");
line->Draw("same");

leg = new TLegend(0.5, 0.7, 0.8, 0.8);
leg->AddEntry(line, "Grosor de CCD: 0.0725 cm", "l");
leg->Draw();

//canv->cd(2);
//Lcut->Draw();
// func2->Draw("same");
canv->Print("Dis_L.pdf");

}
