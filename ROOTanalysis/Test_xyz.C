void Test_xyz(){
//TFile *file = new TFile("../build/muons_100K_vacuum_file.root");
TFile *file = new TFile("./root_files/muons_100K_vacuum_file.root");


// TTree *tree0 = (TTree*) file0->Get("tree");
TTree *tree = (TTree*) file->Get("B02Hits");
//tree->Print();

// int NB = 90;
int NBX = 90;
int NBY = 90;
int NBZ = 90;

double xlow = -0.5;
double xhi = 0.5;
double ylow = -0.5;
double yhi = 0.5;
double zlow = -0.1;
double zhi = 0.1;


TH2F *histXY = new TH2F("Hist_XY", "", NBX, xlow, xhi, NBY, ylow, xhi);
histXY->GetXaxis()->SetTitle("Side X (cm)");
histXY->GetYaxis()->SetTitle("Side Y (cm)");

TH2F *histXZ = new TH2F("Hist_XZ", "", NBX, xlow, xhi, NBZ, zlow, zhi);
histXZ->GetXaxis()->SetTitle("Side X (cm)");
histXZ->GetYaxis()->SetTitle("Side Z (cm)");

TH2F *histZY = new TH2F("Hist_ZY", "", NBY, ylow, yhi, NBZ, zlow, zhi);
histZY->GetXaxis()->SetTitle("Side Z (cm)");
histZY->GetYaxis()->SetTitle("Side Y (cm)");

TH3F *histXYZ = new TH3F("Hist_XYZ", "", NBX, xlow, xhi, NBY, ylow, xhi, NBZ, zlow, zhi);
histXYZ->GetXaxis()->SetTitle("Side X (cm)");
histXYZ->GetYaxis()->SetTitle("Side Y (cm)");
histXYZ->GetZaxis()->SetTitle("Side Z (cm)");

// Fill histograms //
double px, py, pz;
tree->SetBranchAddress("XhitBar", &px);
tree->SetBranchAddress("YhitBar", &py);
tree->SetBranchAddress("ZhitBar", &pz);

int nentries = tree->GetEntries();
   for (int i = 0; i < nentries; i++) {
      tree->GetEntry(i);
      histXY->Fill(px, py);
      histXZ->Fill(px, pz);
      histZY->Fill(py, pz);
      
      histXYZ->Fill(px, py, pz);
   }


// Create Canvas //
TCanvas *canv1 = new TCanvas("canv1","", 2*700, 600);
canv1->Divide(3,1);

canv1->cd(1);
histXY->Draw();

canv1->cd(2);
histXZ->Draw();

canv1->cd(3);
histZY->Draw();

TCanvas *canv = new TCanvas("canv","", 2*700, 600);
histXYZ->Draw();

canv1->Print("XYZ.pdf");

}
