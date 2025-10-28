void Test_thet(){
TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog_nHG_1.root");

// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_200000_PLANES_150x150_RADIO_100.root");
// TFile *file = new TFile("./treesROOT_Barra/Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_0.root");


// TTree *tree0 = (TTree*) file0->Get("tree");
TTree *tree = (TTree*) file->Get("B02Evts");
// tree->Print();


TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_MeV.root"); // INFO ALL_CLUSTERS
// TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_400x700__MeV.root");	// INFO MUONS ONLY
TTree *tree_icn = (TTree*) f_icn->Get("tree");

TFile *file0 = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_SIGMA_1.0_C_0.root");
TTree *tree0 = (TTree*) file0->Get("tree");

// int NB = 90;
int NB = 100;
double tlow = 0;
double thi = TMath::Pi()/2.0;
TH1F *theta_all = new TH1F("theta_all", "Distribuci#acute{o}n angular #theta de todos los muones simulados", NB, tlow, thi);
theta_all->GetXaxis()->SetTitle("#theta (rad)");
theta_all->SetLineColor(2);
theta_all->SetStats(0);

TH1F *theta_in = new TH1F("theta_in", "Distribuci#acute{o}n angular #theta de los muones que impactaron la CCD", NB, tlow, thi);
theta_in->GetXaxis()->SetTitle("#theta(rad)");
theta_in->SetLineColor(2);
// theta_in->SetStats(0);

TH1F *theta_icn = new TH1F("theta_icn", "Distribuci#acute{o}n angular #theta", NB, tlow, thi);
theta_icn->GetXaxis()->SetTitle("#theta(rad)");
theta_icn->SetLineColor(1);
theta_icn->SetStats(0);

TH1F *theta_cut = new TH1F("theta_cut", "Distribuci#acute{o}n angular #theta (#theta > 20^{o})", NB, tlow, thi);
theta_cut->GetXaxis()->SetTitle("#theta(rad)");
theta_cut->SetLineColor(2);
theta_cut->SetStats(0);
// TH1F *theta_incut = new TH1F("theta_incut", "", NB, tlow, thi);


TH1F *theta_pp_in = new TH1F("theta_pp_in", "Distribuci#acute{o}n angular #theta del ICN", NB, tlow, thi);
theta_pp_in->GetXaxis()->SetTitle("#theta(rad)");
theta_pp_in->SetLineColor(4);
// theta_cut->SetStats(0);

// Fill histograms //
// tree0->Draw("thet>>theta_all");
tree->Draw("thetaPri>>theta_all");

// tree0->Draw("thet>>theta_in", "l>0 ");
// tree->Draw("thetaPri>>theta_in", "nHitBar>0");
tree->Draw("thetaPri>>theta_in", "LengthMuLAr>0");
// tree->Draw("thetaPri>>theta_in", "nHitBar>0 && LengthMuLAr>0");


tree->Draw("thetaPri>>theta_cut", "nHitBar>0 && thetaPri > 20*TMath::Pi()/180");
// tree->Draw("thet>>theta_incut", "thet>22*TMath::Pi()/180");
tree_icn->Draw("thet>>theta_icn");

tree0->Draw("thet>>theta_pp_in", " l > 0 ");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0.01, 85*TMath::Pi()/180);
func1->SetLineColor(4);
func1->SetParameter(0, 5000);

// ====================== Función para la CCD =========================== //
// TF1 *func2 = new TF1("func2", "[0]*((12.169116/1)*sin(x)*(cos(x))^3+(4.06464/(2*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0.01, 85*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((3.042279/1)*sin(x)*(cos(x))^3+(0.25404/(1*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0.01, 85*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((3.042279/1)*sin(x)*(cos(x))^3+(2*(0.1389825 + 0.1150575)/(1*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0.01, 85*TMath::Pi()/180);
// double Ah = 0.54;
// double Al = 0.06525;
// double Ac = 0.0435;

double Ah = 0.2975625;
double Al = 0.05752875;
double Ac = 0.0271875;

TF1 *func2 = new TF1("func2", Form("[0]*((%f)*sin(x)*(cos(x))^3+(2*(%f + %f)/(1*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", Ah, Al, Ac), 0.01, 85*TMath::Pi()/180);
func2->SetLineColor(4);
// ====================== Funcion para la barra ========================  //
// TF1 *func2 = new TF1("func2", "[0]*((10)*sin(x)*(cos(x))^3+(22/(TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((10)*sin(x)*(cos(x))^3+(22/4)*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);
// TF1 *func2 = new TF1("func2", "[0]*((1.25)*sin(x)*(cos(x))^3+(11/(TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0, 90*TMath::Pi()/180);

func2->SetParameter(0, 2000);

// Fit functions //
theta_all->Fit("func1", "R");
theta_in->Fit("func2", "R");

double Prob1 = func1->GetProb();
double chi1 = func1->GetChisquare();
// double NFD1 = func1->GetN

double Prob2 = func2->GetProb();
double chi2 = func2->GetChisquare();

std::cout << "ChiSquare1 = "<< chi1 << std::endl;
std::cout << "Prob1 = "<< Prob1 << std::endl;

std::cout << "ChiSquare2 = "<< chi2 << std::endl;
std::cout << "Prob2 = "<< Prob2 << std::endl;

// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);
canv->cd(1);
theta_all->Draw();
func1->Draw("same");

TLegend *leg = new TLegend(0.5, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->AddEntry(func1, "Asin#theta cos^{2}#theta", "l");
leg->AddEntry(theta_all, "Datos Simulados", "f");
leg->Draw();

theta_icn->Scale(7.2);
canv->cd(2);
// theta_in->Draw("hist");
theta_cut->Draw("hist");
theta_icn->Draw("hist same");
// func2->Draw("same");
// theta_pp_in->Draw("hist same");

leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
// leg->SetHeader("", "C");
// leg->SetFillStyle(0);
// leg->AddEntry(func1, "A [A_{s} sin#theta cos^{3}#theta + #left(#frac{2(A_{c} + A_{l})}{#pi}#right)sin^{2}#theta cos^{2}#theta]", "l");
// leg->AddEntry((TObject*)0, "", "");
// leg->AddEntry((TObject*)0, "+ #left(#frac{A_c + A_l}{#pi}#right)sin^{2}#theta cos^{2}#theta]", " ");
// leg->AddEntry((TObject*)0, "", "");
leg->AddEntry(theta_all, "Simulaci#acute{o}n de Geant4", "f");
leg->AddEntry(theta_icn, "Datos ICN", "f");
// leg->AddEntry(theta_pp_in, "Sim-PP", "f");
leg->Draw();




canv->Print("Dis_thet.pdf");

}