void Muonlengthdistribution(){

Double_t LMuLAr;
Double_t l;


TFile *file = TFile::Open("muons_nfile.root");
TTree *tree = (TTree*)file->Get("B02Evts");

TFile *file1 = TFile::Open("lengthDvalues.root");
TTree *tree1 = (TTree*)file1->Get("data");

gStyle->SetOptStat(0);  // 0 disables the stats box

TH1F *GeantHis = new TH1F("Geant4Histogram", "Length Distribution; Length (cm); Events", 300, 0, 250);
TH1F *firstpHis = new TH1F("FirstpHistogram", "Length Distribution; Length (cm); Events", 300, 0, 250);
tree->SetBranchAddress("LMuLAr", &LMuLAr); 
tree1->SetBranchAddress("l", &l); 


int nentries = tree->GetEntries();
for (int i = 0; i < nentries; i++){
  tree->GetEntry(i);
  if (LMuLAr > 0){
  GeantHis->Fill(LMuLAr);
  }
} 
int Nentries = tree1->GetEntries();
for (int i = 0; i < 15998; i++){
  tree1->GetEntry(i);
  firstpHis->Fill(l);
} 


TCanvas *c1 = new TCanvas("c1", "", 900,700);
   c1->SetGrid();
   GeantHis->Draw();
   firstpHis->SetLineColor(kRed);
   firstpHis->Draw("same");
   
   TLegend *legend = new TLegend(0.55, 0.6, 0.9, 0.8);
	legend->SetTextSize(0.02);
	legend->AddEntry(GeantHis, "Length Muon Distribution (Geant4)", "lp");
	legend->AddEntry(firstpHis, "Length Muon Distribution (first principles) ", "lp");
	legend->Draw();
	
c1->Print("LengthDistribution.pdf");



}

