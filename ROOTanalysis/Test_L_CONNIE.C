void Test_L_CONNIE(){
// TFile *file = new TFile("./root_files/muons_2M_vacuum_file.root");
TFile *file = new TFile("./root_files/muons_1M_vacuum_CONNIE_420x1022_file.root");
TTree *tree = (TTree*) file->Get("B02Evts");

TFile *file0 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_CONNIE_NMUONS_1000000_PLANES_1.7_RADIO_7_CCDSIZE_420X1022_C_0.root");
TTree *tree0 = (TTree*) file0->Get("tree");


int NB = 200;
double tlow = 0;
double thi = 0.4;
TH1F *L = new TH1F("L", "Distance Distribution (CCD size: 1.533x0.63x0.068 cm)", NB, tlow, thi);
L->GetXaxis()->SetTitle("Distance (cm)");
//L->SetStats(0);

TH1F *LPP = new TH1F("LPP", "Distance Distribution (CCD size: 1.533x0.63x0.068 cm)", NB, tlow, thi);
LPP->GetXaxis()->SetTitle("Distance (cm)");
LPP->SetLineStyle(1);
LPP->SetLineColor(2);


TH1F *Lcut = new TH1F("Lcut", "", NB, tlow, thi);

// Fill histograms //
tree->Draw("LengthMuLAr>>L", "LengthMuLAr > 0");
tree0->Draw("l>>LPP", "l > 0");


//TLine *line = new TLine(0.0725,0,0.0725,30000);
TLine *line = new TLine(0.068,0,0.068,14000);
line->SetLineStyle(2);

LPP->Scale(1.0);

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
//canv->Divide(2,1);
canv->cd(1);
L->Draw();
LPP->Draw("hist same");
line->Draw("same");
// Lcut->Draw("same");

TLegend *leg = new TLegend(0.5, 0.7, 0.8, 0.8);
leg->AddEntry(L, "Distance distribution (GEANT4)", "lep");
leg->AddEntry(LPP, "Distance distribution (SimPP)", "lep");
leg->AddEntry(line, "CCD Thickness: 0.068 cm", "l");
leg->Draw();

//canv->cd(2);
//Lcut->Draw();
// func2->Draw("same");
canv->Print("Dis_L.pdf");

}
