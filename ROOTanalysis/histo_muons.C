void histo_muons(){

// TFile *f_geant = new TFile("./root_files/muons_1M_vacuum_400x525_file.root");
TChain *f_geant_tree = new TChain("B02Evts");
// f_geant_tree->Add("./root_files/muons_1M_vacuum_250x529_file_m_old.root");

f_geant_tree->Add("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog.root");
// f_geant_tree->Add("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_0.root");
// f_geant_tree->Add("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_1.root");
// f_geant_tree->Add("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_2.root");
// f_geant_tree->Add("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_3.root");
// f_geant_tree->Add("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_4.root");


// TFile *f_geant = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog.root");
// TFile *f_geant = new TFile("./root_files/muons_1M_vacuum_300x529_file.root");
// TFile *file = new TFile("./root_files/muons_1M_vacuum_file.root");
// TTree *tree_geant = (TTree*) f_geant->Get("B02Evts");
//TTree *tree = (TTree*) file->Get("B02Hits");

// TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_400x700_MeV.root"); // INFO ALL_CLUSTERS
// TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_300x529_EXPOSURE_4504_MeV.root"); // INFO ALL_CLUSTERS
TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_250x529_EXPOSURE_4504_MeV.root"); // INFO ALL_CLUSTERS
TTree *tree_icn = (TTree*) f_icn->Get("tree");

TFile *f_conn = new TFile("Edep_CONNIE_NSAMP400_MeV.root");
TTree *tree_conn = (TTree*) f_conn->Get("tree");

TLatex lat;

int NB = 150;
double tlow = 0.;
double thi = 1000.;

TH1F *edep_icn = new TH1F("edep_icn", "Energy Spectrum (Sim. PP - All Clusters (ICN))", NB, tlow, thi);
edep_icn->GetXaxis()->SetTitle("Energy (MeV)");
edep_icn->SetLineStyle(1);
edep_icn->SetLineColor(1);
edep_icn->SetStats(0);

TH1F *edep_g4 = new TH1F("edep_g4", "Simulation no Birks", NB, tlow, thi);
// edep_g4->SetStats(0);
edep_g4->SetLineStyle(1);
edep_g4->SetLineColor(2);

TH1F *edep_conn = new TH1F("edep_conn", "CONNIE 2021-2022", NB, tlow, thi);
edep_conn->SetStats(0);
edep_conn->SetLineStyle(2);
edep_conn->SetLineColor(1);

TH1F *edep_geant_scale = new TH1F("edep_geant_scale", "Simulation GEANT4 scale", NB, tlow, thi);
edep_geant_scale->SetStats(0);
edep_geant_scale->SetLineStyle(2);
edep_geant_scale->SetLineColor(1);

// === Template for fit ==
int NBmu = 150;
double maxEd = 1500.;
double muonsbw = maxEd/NBmu;

TH1F *muons = new TH1F("muons","",NBmu,0,maxEd);
muons->GetXaxis()->SetTitle("Energy (KeV)");
muons->GetYaxis()->SetTitle("events");

f_geant_tree->Draw("EevtBar*1000>>muons", "EevtBar>0");


// ============= GEANT4 histograms =========== //
f_geant_tree->Draw("EevtBar*1000>>edep_g4", "EevtBar>0"); // GEANT4 INFO (NO BIRKS)
f_geant_tree->Draw("EevtBar*0.92*1000>>edep_geant_scale", "EevtBar>0"); // GEANT4 INFO (NO BIRKS)


// ============= EXPERIMENTAL histograms =========== //
tree_icn->Draw("edep*1000.>>edep_icn", "edep>0"); // EXPERIMENTAL INFO in KeV
edep_icn->SetMaximum(400);
// tree0->Draw("edep>>edep0", "edep>0"); // EXPERIMENTAL INFO
tree_conn->Draw("edep>>edep_conn", "edep>0"); // EXPERIMENTAL INFO CONNIE

double cont_g4 = edep_g4->Integral();
double cont_g4scale = edep_geant_scale->Integral();
double cont_icn = edep_icn->Integral();
double cont_conn = edep_conn->Integral();
double cont_muons = muons->Integral();

cout << "Integral muons: " << cont_muons <<endl;

// ============ Create Canvas ============== //
TCanvas *canv = new TCanvas("canv","Edep", 2*800, 600);
canv->Divide(1,1);
canv->cd(1);
// edep_g4->Draw("hist"); 	// edepG4 no Birks
muons->Draw("hist");
// edep_pp->Draw("hist same"); 	// edepPP
// // edep_conn->Draw("he0 same"); 	// CONNIE data
// edep_icn->Draw("hist"); 	// ICN data 
// edep_geant_scale->Draw("hist same");

// TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
// // leg->AddEntry(edep, "SimG4-Birks: 0.09 cm/MeV", "lep");
// leg->AddEntry(edep_g4, "Sim-GEANT4", "LEP");
// leg->AddEntry(edep_icn, "All Clusters (ICN-NSAMP324)", "LEP");
// //leg->AddEntry(edep2, "Datos CONNIE-NSAMP400", "LEP");
// leg->AddEntry(edep_pp, "Sim-PP", "LEP");
// leg->AddEntry(edep_geant_scale, "Sim-GEANT4 (0.92)", "LEP");
// leg->Draw();



TFile *fout = new TFile("ccdhisto.root", "recreate");
edep_icn->Write();
fout->Close();


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
  for (int i=50;i<NBmu;i++){
   if (muons->GetBinContent(i+1) > muons->GetBinContent(peakBin))
     peakBin = i+1;
  } //for i
  printf("Peak: %d  \n", peakBin);
  printf("Peak: %3.3f MeV \n", peakBin*maxEd/NBmu);
  double Epeak = peakBin*maxEd/NBmu;




//canv->cd(2);
//edep_cut->Draw();
// func2->Draw("same");
// canv->Print("Dis_edep_GEANT-PP.pdf");

}
