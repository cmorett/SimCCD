void Test_edep(){

// TFile *f_geant = new TFile("./root_files/muons_2M_vacuum_file.root");
// TFile *f_geant = new TFile("./root_files/muons_1M_vacuum_250x529_file_m.root");
TFile *f_geant = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_nHG_1.root");
// TFile *file = new TFile("./root_files/muons_1M_vacuum_file.root");
TTree *tree_geant = (TTree*) f_geant->Get("B02Evts");
//TTree *tree = (TTree*) file->Get("B02Hits");

TFile *f_icn = new TFile("../../../../Simulacion_ab_initio//Edep_NSAMP324_MeV.root"); // INFO ALL_CLUSTERS
// TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_250x529_EXPOSURE_4504_MeV.root"); // INFO ALL_CLUSTERS
// TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_400x700__MeV.root");	// INFO MUONS ONLY
TTree *tree_icn = (TTree*) f_icn->Get("tree");

TFile *f_pp = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_SIGMA_1.0_C_0.root");
// TFile *f_pp = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5x1.5_RADIO_8_CCDSIZE_250x529_SIGMA_LV_1.0_.root");
TTree *tree_pp = (TTree*) f_pp->Get("tree");

TFile *f_conn = new TFile("Edep_CONNIE_NSAMP400_MeV.root");
TTree *tree_conn = (TTree*) f_conn->Get("tree");

// TFile *file3 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_400X525_C_0.root");
// TTree *tree3 = (TTree*) file3->Get("tree");

TLatex lat;

int NB = 150;
double tlow = 0;
double thi = 1000;

TH1F *edep_geant_birks = new TH1F("edep_geant_birks", "Energy Spectrum", NB, tlow, thi);
edep_geant_birks->GetXaxis()->SetTitle("Energy (MeV)");
edep_geant_birks->SetLineStyle(1);
edep_geant_birks->SetLineColor(3);
edep_geant_birks->SetStats(0);

//TH1F *edep0 = new TH1F("edep0", "Lab. Det. ", NB, tlow, thi);
TH1F *edep_icn = new TH1F("edep_icn", "Espectro de Energ#acute{i}as Depositada (ICN)", NB, tlow, thi);
// TH1F *edep_icn = new TH1F("edep_icn", "Espectro de Energ#acute{i}as Depositada (#theta > 20)", NB, tlow, thi);
edep_icn->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");
edep_icn->SetLineStyle(1);
edep_icn->SetLineColor(1);
edep_icn->SetStats(0);

TH1F *edep_g4 = new TH1F("edep_g4", "Distribuci#acute{o}n de Energ#acute{i}a Depositada", NB, tlow, thi);
edep_g4->SetStats(0);
edep_g4->SetLineStyle(1);
edep_g4->SetLineColor(2);
edep_g4->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");

TH1F *edep_conn = new TH1F("edep_conn", "Espectro de Energ#acute{i}as Depositada (CONNIE 2021-2022)", NB, tlow, thi);
edep_conn->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");
edep_conn->SetStats(0);
edep_conn->SetLineStyle(1);
edep_conn->SetLineColor(1);

TH1F *edep_pp = new TH1F("edep_pp", "Simulation PP", NB, tlow, thi);
edep_pp->SetStats(0);
edep_pp->SetLineStyle(1);
edep_pp->SetLineColor(4);

TH1F *edep_geant_scale = new TH1F("edep_geant_scale", "Distribuci#acute{o}n de Energ#acute{i}as Depositadas (#theta > 20^{o})", NB, tlow, thi);
edep_geant_scale->SetStats(0);
edep_geant_scale->SetLineStyle(1);
edep_geant_scale->SetLineColor(2);
edep_geant_scale->GetXaxis()->SetTitle("Energ#acute{i}a (KeV)");

// === Template for fit ==
int NBmu = 200;
double maxEd = 1000.;

double muonsbw = maxEd/NBmu;

TH1F *muons = new TH1F("muons","",NBmu,0,maxEd);
muons->GetXaxis()->SetTitle("Energy (KeV)");
muons->GetYaxis()->SetTitle("events");

// ============= GEANT4 histograms =========== //
tree_geant->Draw("WevtBar*1000>>edep_geant_birks", "nHitBar>0"); // GEANT4 INFO (BIRKS)
tree_geant->Draw("EevtBar*1000>>edep_g4", "EevtBar>0"); // GEANT4 INFO (NO BIRKS)
//tree->Draw("Ehitbar>>edep_cut", "thet>22*TMath::Pi()/180 & edep>0");
tree_geant->Draw("EevtBar*0.904*1000>>edep_geant_scale", "EevtBar>0 && thetaPri>20*TMath::Pi()/180"); // GEANT4 INFO (NO BIRKS)


// ============= EXPERIMENTAL histograms =========== //
tree_icn->Draw("edep*1000>>edep_icn", "edep>0"); // EXPERIMENTAL INFO in KeV
// tree0->Draw("edep>>edep0", "edep>0"); // EXPERIMENTAL INFO
tree_conn->Draw("edep*1000>>edep_conn", "edep>0"); // EXPERIMENTAL INFO CONNIE

edep_icn->SetMaximum(350);
edep_conn->SetMaximum(800);

// ============= SIM AB INITIO histograms =========== //
// tree_pp->Draw("edep*0.904*1000>>edep_pp", "edep>0 && thet>22*TMath::Pi()/180"); // SIM_AB_INITIO INFO
tree_pp->Draw("edep*0.904*1000>>edep_pp", "edep>0"); // SIM_AB_INITIO INFO


tree_geant->Draw("EevtBar*1000*0.92>>muons", "EevtBar>0");
// tree0->Draw("edep>>muons", "edep>0");

double cont_g4 = edep_g4->Integral();
double cont_g4scale = edep_geant_scale->Integral();
double cont_icn = edep_icn->Integral();
double cont_conn = edep_conn->Integral();
double cont_pp = edep_pp->Integral();

double cont_muons = muons->Integral();

cout << "Integral gean4: " << cont_g4 <<endl;
cout << "Integral muons: " << cont_muons <<endl;

// ========== Scale histograms =========== //
//edep->Scale(0.65, "");
// edep_pp->Scale(cont_g4/cont_pp);
//edep->SetLineColor(2);

// edep0->Scale(63.);
// edep0->Scale(50);
//edep0->SetLineColor(4);

// edep_g4->Scale(1);
//edep1->SetLineColor(1);

// edep_geant_scale->Scale(0.4);

// edep_icn->Scale(3.7);
// edep_icn->Scale(2.8); // Para 20°
edep_icn->Scale(6.4); // Para 25°
// edep_conn->Scale(3.0);

// ======================================= //



// ============ Create Canvas ============== //
TCanvas *canv = new TCanvas("canv","Edep", 2*800, 600);
canv->Divide(2,1);
canv->cd(1);
edep_geant_scale->Draw("hist same");
// edep_g4->Draw("hist"); 	// edepG4 no Birks
// edep_geant_birks->Draw("hist same");
// edep_pp->Draw("hist same"); 	// edepPP
// edep_conn->Draw("hist"); 	// CONNIE data
edep_icn->Draw("hist same"); 	// ICN data 

TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
// leg->AddEntry(edep, "SimG4-Birks: 0.09 cm/MeV", "lep");
// leg->AddEntry(edep_g4, "Simulaci#acute{o}n Geant4", "LP");
// leg->AddEntry(edep_icn, "Datos de todos los clusters (NSAMP324)", "LEP");
leg->AddEntry(edep_icn, "Datos de muones ICN (NSAMP324)", "LP");
// leg->AddEntry(edep_conn, "Datos de todos los clusters (NSAMP400)", "LEP");
// leg->AddEntry(edep_pp, "Simulaci#acute{o}n ab initio", "LP");
leg->AddEntry(edep_geant_scale, "Simulaci#acute{o}n de Geant4 (escalada: 0.904)", "LP");
leg->Draw();


//canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
canv->Print("Dis_edep_GEANT-PP.pdf");

}
