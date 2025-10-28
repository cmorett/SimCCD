void Test_L(){
// TFile *file = new TFile("./root_files/muons_2M_vacuum_file.root");
// TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_2.root");
// TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_0.root");
TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_nHG_1.root");

TTree *tree = (TTree*) file->Get("B02Evts");


TFile *file0 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_SIGMA_1.0_C_0.root");
TTree *tree0 = (TTree*) file0->Get("tree");

TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_MeV.root"); // INFO ALL_CLUSTERS
// TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_400x700__MeV.root");	// INFO MUONS ONLY
TTree *tree_icn = (TTree*) f_icn->Get("tree");


int NB = 150;
double tlow = 0.;
double thi = 0.4;

// TH1F *L = new TH1F("L", "Distance Distribution (CCD size: 0.9x0.6x0.0725 cm)", NB, tlow, thi);
// L->GetXaxis()->SetTitle("Distance (cm)");
// TH1F *L = new TH1F("L", "Distribuci#acute{o}n de Longitudes con corte angular de 25^{o}", NB, tlow, thi);
TH1F *L = new TH1F("L", "Distribuci#acute{o}n de Longitudes (#theta > 20^{0})", NB, tlow, thi);
L->GetXaxis()->SetTitle("Distancia (cm)");
L->SetStats(0);
L->SetLineColor(2);

// TH1F *LPP = new TH1F("LPP", "Distribuci#acute{o}n de Longitudes con corte angular de 25^{o}", NB, tlow, thi);
TH1F *LPP = new TH1F("LPP", "Distribuci#acute{o}n de Longitudes", NB, tlow, thi);
LPP->GetXaxis()->SetTitle("Distancia (cm)");
LPP->SetStats(0);
LPP->SetLineStyle(1);
LPP->SetLineColor(4);

TH1F *L_ICN = new TH1F("L_ICN", "Distribuci#acute{o}n de Longitudes", NB, tlow, thi);
L_ICN->GetXaxis()->SetTitle("Distancia (cm)");
L_ICN->SetStats(0);
L_ICN->SetLineColor(1);


TH1F *Lcut = new TH1F("Lcut", "", NB, tlow, thi);

// Fill histograms //
// tree->Draw("LengthMuLAr>>L", "LengthMuLAr>0");  
// tree->Draw("LengthMuLAr>>L", "nHitBar>0");
// tree->Draw("LengthMuLAr>>L", "nHitBar>0 && LengthMuLAr>0");


tree0->Draw("l>>LPP", "l > 0");

tree->Draw("LengthMuLAr>>L", "LengthMuLAr>0 && thetaPri>20*TMath::Pi()/180");
// tree0->Draw("l>>LPP", "l > 0 && thet>25*TMath::Pi()/180");

tree_icn->Draw("l>>L_ICN", "l > 0");


//TLine *line = new TLine(0.0725,0,0.0725,30000);
TLine *line = new TLine(0.0725,0,0.0725,8000);
line->SetLineStyle(2);
line->SetLineWidth(2);

// LPP->Scale(0.94);
// LPP->Scale(0.1);
// L->Scale(10.);
L_ICN->Scale(6);

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
L->Draw("hist");
// LPP->Draw("hist same");
L_ICN->Draw("hist same");
line->Draw("same");
Lcut->Draw("same");

TLegend *leg = new TLegend(0.5, 0.7, 0.8, 0.8);
leg->AddEntry(L, "Simulaci#acute{o}n Geant4", "lp");
// leg->AddEntry(LPP, "Simulaci#acute{o}n ab initio", "lp");
leg->AddEntry(L_ICN, "Datos ICN", "l");
leg->AddEntry(line, "Grosor de CCD: 0.0725 cm", "l");
leg->Draw();

//canv->cd(2);
//Lcut->Draw();
// func2->Draw("same");
canv->Print("Dis_L.pdf");

}
