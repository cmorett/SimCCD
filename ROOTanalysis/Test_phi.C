void Test_phi(){
//TFile *file = new TFile("Sim_ab_initio_NMUONS_300000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_400000.root");
// TFile *file = new TFile("Sim_ab_initio_NMUONS_10000_PLANES_3.0x3.0_RADIO_12_.root");
TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog.root");

// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_200000_PLANES_150x150_RADIO_100.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_500000_PLANES_150x150_RADIO_100(1).root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_100000_PLANES_150x150_RADIO_450_0.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_1.root");
// TFile *file = new TFile("Sim_ab_initio_Barra_NMUONS_300000_PLANES_150x150_RADIO_450_0.root");
TTree *tree = (TTree*) file->Get("B02Evts");

TFile *file_pp = new TFile("../../../../Simulacion_ab_initio/Sim_ab_initio_NMUONS_1000000_PLANES_1.5_RADIO_8_CCDSIZE_250X529_SIGMA_1.0_C_0.root");
TTree *tree_pp = (TTree*) file_pp->Get("tree");


TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_MeV.root"); // INFO ALL_CLUSTERS
// TFile *f_icn = new TFile("../../../../Simulacion_ab_initio/Edep_NSAMP324_400x700__MeV.root");	// INFO MUONS ONLY
TTree *tree_icn = (TTree*) f_icn->Get("tree");

int NB = 40;
// double tlow = - TMath::Pi() - 0.1;
double tlow = 0;
double thi = 2 * TMath::Pi() + 0.1;
TH1F *phi_all = new TH1F("phi_all", "Distribuci#acute{o}n angular #phi de todos los muones simulados", NB, tlow, thi);
phi_all->GetXaxis()->SetTitle("#phi (rad)");
phi_all->SetStats(0);
phi_all->SetLineColor(2);


TH1F *phi_in = new TH1F("phi_in", "Distribuci#acute{o}n angular #phi de los muones que impactaron la CCD", NB, tlow, thi);
phi_in->SetStats(0);
phi_in->GetXaxis()->SetTitle("#phi (rad)");
phi_in->SetLineColor(2);

TH1F *phi_cut = new TH1F("phi_cut", "Distribuci#acute{o}n angular #phi (para #theta > 20^{o})", NB, tlow, thi);
phi_cut->SetStats(0);
phi_cut->GetXaxis()->SetTitle("#phi (rad)");
phi_cut->SetLineColor(2);

TH1F *phi_icn = new TH1F("phi_icn", "Distribuci#acute{o}n angular #phi", NB, tlow, thi);
phi_icn->SetStats(0);
phi_icn->GetXaxis()->SetTitle("#phi (rad)");
phi_icn->SetLineColor(1);

// Fill histograms //
// tree->Draw("(phiPri + TMath::Pi())>>phi_all");
// tree->Draw("(phiPri + TMath::Pi())>>phi_in", " EevtBar > 0 ");

tree_pp->Draw("phi>>phi_all");
tree_pp->Draw("phi>>phi_in", " edep > 0 ");
tree_pp->Draw("phi>>phi_cut", " edep > 0 && thet>20*TMath::Pi()/180");
tree_icn->Draw("phi>>phi_icn");


// tree->Draw("phi>>theta_incut", "thet>22");

// Define fuctions //
TF1 *func1 = new TF1("func1", "[0]", 0.01, 2*TMath::Pi() - 0.1);
func1->SetLineColor(4);
// func1->SetParameter(0, 7500);

double Ah = 0.2975625;
double Al = 0.05752875;
double Ac = 0.0271875;

// ======================================= Funcionde ajuste para la CCD  =============================================== ###
// TF1 *func2 = new TF1("func2", Form("[0]*(((%f)/(1 * TMath::Pi())) + ((%f)/2)*abs(cos(x)) + ((%f)/2)*abs(sin(x)))", Ah, Al, Ac), 0.05,2*TMath::Pi() - 0.1); 
TF1 *func2 = new TF1("func2", Form("[0]*(((%f)/(1 * TMath::Pi())) + ((%f)/4)*abs(cos(x)) + ((%f)/4)*abs(sin(x)))", Ah, Al, Ac), 0.05,2*TMath::Pi() - 0.1); 
func2->SetLineColor(4);
// ======================================= Funcionde ajuste para la barra =============================================== ###
// TF1 *func2 = new TF1("func2", "[0]*((20/TMath::Pi()) + (5)*abs(cos(x)) + (1/2)*abs(sin(x)))",0.01,2*TMath::Pi() - 0.1);


// TF1 *func2 = new TF1("func2", "[0]*((8 * 0.8247739/(2 * TMath::Pi())) + (4 * 0.0499529/4)*(cos(x)) + (4 * 0.0376767715/4)*(sin(x)))", 0.01,2*TMath::Pi() - 0.1); 
// TF1 *func2 = new TF1("func2", "[0]*((8 * 0.8247739/(2 * TMath::Pi())) + (4 * 0.0499529/4)*abs(cos(x)))", 0.01,2*TMath::Pi() - 0.1); // EL mejor
// func2->SetParameter(0, 1600);
//func2->SetParameter(1, 1000);

// Fit functions //
phi_all->Fit("func1", "R");
phi_in->Fit("func2", "R");
// phi_cut->Fit("func2", "R");
double Prob1 = func1->GetProb();
double chi1 = func1->GetChisquare();

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
phi_all->Draw();
func1->Draw("same");

TLegend *leg = new TLegend(0.5, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
leg->AddEntry(func1, "A_{2}", "l");
leg->AddEntry(phi_all, "Datos Simulados", "f");
leg->Draw();


TLine *line1 = new TLine(TMath::Pi()/2,0,TMath::Pi()/2,2400);
line1->SetLineStyle(2);
line1->SetLineWidth(2);

TLine *line2 = new TLine(TMath::Pi(),0,TMath::Pi(),2500);
line2->SetLineStyle(2);
line2->SetLineWidth(2);

TLine *line3 = new TLine(3*TMath::Pi()/2,0,3*TMath::Pi()/2,2400);
line3->SetLineStyle(2);
line3->SetLineWidth(2);

canv->cd(2);
phi_cut->SetMinimum(0);
phi_in->SetMinimum(0);
// phi_in->Draw();
phi_cut->Draw("hist");
// func2->Draw("same"); 
phi_icn->Scale(6.5);
phi_icn->Draw("hist same");
line1->Draw("same");
line2->Draw("same");
line3->Draw("same");

leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
// leg->SetHeader("", "C");
// leg->SetFillStyle(0);
// leg->AddEntry(func1, "A [#frac{A_{s}}{#pi}+ #left(#frac{A_{l}}{4}#right)#left|cos#theta #right| + #left(#frac{A_{c}}{4}#right) #left|sin#theta #right| ]", "l");
// leg->AddEntry((TObject*)0, "", "");
// leg->AddEntry((TObject*)0, "+ #left(#frac{A_{c}}{4}#right) #left|sin#theta #right| ]", " ");
// leg->AddEntry((TObject*)0, "", "");
leg->AddEntry(phi_all, "Simulaci#acute{o}n de Geant4 ", "f");
leg->AddEntry(phi_icn, "Datos ICN ", "f");
leg->Draw();

canv->Print("Dis_phi.pdf");

}
