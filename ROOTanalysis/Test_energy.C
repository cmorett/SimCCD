
double Smith_Dull(double *lx , double *lpar){

    double E_mu = lx[0];
    double theta = TMath::Pi() * lpar[0]/180;

    //## ------------- Constantes físicas -------------- ##
    double A = 2 * pow(10,9);
    // double k = 8.0 / 3.0;
    double k = 2.645;
    double b = 0.771;
    double lambda_pi = 120.0;    // g/cm^2
    double y_0 = 1000.0;      // g/cm^2
    double r = 0.76;
    double b_mu = 0.8;   

    // # if units == 0:
    double a = 2.5;      //## MeV cm^2/g
    double m_mu = 105.7;    //## MeV/c^2 
    double m_pi = 139.6;    //## MeV/c^2

    double tau_mu_0 = 2.2 * pow(10,-6);   //## s
    double tau_0 = 2.6 * pow(10, -8);     //## s
    double rho_0 = 0.00129; //## g/cm^3
    double c = TMath::C()* 100; //## cm/s

    //### ---------------------- Parámetros ---------------------- ###
    double E_pi = (1 / r) * (E_mu + a * y_0 * ((1/TMath::Cos(theta)) - 0.1));
    double B_mu = (b_mu * m_mu * y_0)/(tau_mu_0 * rho_0 * c);
    // double B_mu = (b_mu * m_mu * y_0 * c)/(tau_mu_0 * rho_0);
    double P_mu = pow((0.1 * TMath::Cos(theta)) * (1 - (a * (y_0 *(1/TMath::Cos(theta)) - 100))/( r * E_pi)), ((B_mu)/((r * E_pi + 100 * a) * TMath::Cos(theta))));
    double j_pi = (m_pi * y_0)/(c * tau_0 * rho_0);
    // double j_pi = (m_pi * y_0 * c)/(tau_0 * rho_0);

    //## Intensidad diferencial
    //# C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    double C_1 = pow(E_pi, -k)* P_mu * lambda_pi * b * j_pi;
    double C_2 = E_pi * TMath::Cos(theta);
    double C_3 = b * j_pi;

    //# return (C_1 * np.sin(theta)) / (C_2 + C_3)
    double En_k = (C_1 ) / (C_2 + C_3);
    return En_k;
    }
/// ============================================================= ///

double Smith_Dull_Log(double *lx , double *lpar){

    double E_mu = pow(10, lx[0]);
    double theta = TMath::Pi() * lpar[0]/180;

    //## ------------- Constantes físicas -------------- ##
    double k = 8.0 / 3.0;
    double b = 0.771;
    double lambda_pi = 120;    // g/cm^2
    double y_0 = 1000;      // g/cm^2
    double r = 0.76;
    double b_mu = 0.8;   

    // # if units == 0:
    double a = 2.5;      //## MeV cm^2/g
    double m_mu = 105.7;    //## MeV/c^2 
    double m_pi = 139.6;    //## MeV/c^2

    double tau_mu_0 = 2.2 * pow(10,-6);   //## s
    double tau_0 = 2.6 * pow(10, -8);     //## s
    double rho_0 = 0.00129; //## g/cm^3
    double c = TMath::C()* 100; //## cm/s

    //### ---------------------- Parámetros ---------------------- ###
    double E_pi = (1.0 / r) * (E_mu + a * y_0 * ((1/TMath::Cos(theta)) - 0.1));
    double B_mu = (b_mu * m_mu * y_0)/(tau_mu_0 * rho_0 * c);
    // double B_mu = (b_mu * m_mu * y_0 * c)/(tau_mu_0 * rho_0);
    double P_mu = pow((0.1 * TMath::Cos(theta)) * (1.0 - (a * (y_0 *(1/TMath::Cos(theta)) - 100.0))/( r * E_pi)), ((B_mu)/((r * E_pi + 100.0 * a) * TMath::Cos(theta))));
    double j_pi = (m_pi * y_0)/(c * tau_0 * rho_0);
    // double j_pi = (m_pi * y_0 * c)/(tau_0 * rho_0);

    //## Intensidad diferencial
    //# C_1 = E_pi ** (-k) * P_mu * lambda_pi * b * j_pi
    double C_1 = pow(E_pi, -k)* P_mu * lambda_pi * b * j_pi;
    double C_2 = E_pi * TMath::Cos(theta);
    double C_3 = b * j_pi;

    //# return (C_1 * np.sin(theta)) / (C_2 + C_3)
    double En_k = (C_1 ) / (C_2 + C_3);
    double log_ek = En_k;

    return log_ek;
    // return log10(log_ek);
    // return pow(10, log_ek);
    // return pow(10, log10(log_ek));
    }
    


void Test_energy(){
TFile *file = new TFile("./root_files/muons_1M_vacuum_250x529_file_m_old_SDLog.root");
TTree *tree = (TTree*) file->Get("B02Evts");

// int NB = 100;
int NB = 70;
double tlow = 1;
double thi = 1 * pow(10,5);
double xbins[NB+1];
double low_bin;

double Delta_log = (log10(thi) - log10(tlow))/(NB -1);
for (int i = 0; i<NB;i++){
    low_bin = pow(10, log10(tlow) + i * Delta_log);
    // cout<<low_bin <<endl;
    xbins[i] = low_bin;
};

xbins[NB] = thi;


TH1F *thet_0_6 = new TH1F("thet_0_6", "", NB, xbins);
TH1F *thet_43_47= new TH1F("thet_43_47", "", NB, xbins);
TH1F *thet_73_77 = new TH1F("thet_73_77", "", NB, xbins);

thet_0_6->SetLineColor(1);
thet_43_47->SetLineColor(2);
thet_73_77->SetLineColor(3);


// Fill histograms //
tree->Draw("EevtPri*1000>>thet_0_6", "thetaPri>0 && thetaPri<(5*TMath::Pi()/180) && EevtBar > 0");
tree->Draw("EevtPri*1000>>thet_43_47", "thetaPri>(44*TMath::Pi()/180) && thetaPri<(46*TMath::Pi()/180) && EevtBar > 0");
tree->Draw("EevtPri*1000>>thet_73_77", "thetaPri>(73*TMath::Pi()/180) && thetaPri<(77*TMath::Pi()/180) && EevtBar > 0");

// tree->Draw("EevtPri*1000>>thet_0_6", "thetaPri>0 && thetaPri<(5*TMath::Pi()/180)");
// tree->Draw("EevtPri*1000>>thet_43_47", "thetaPri>(44*TMath::Pi()/180) && thetaPri<(46*TMath::Pi()/180)");
// tree->Draw("EevtPri*1000>>thet_73_77", "thetaPri>(74*TMath::Pi()/180) && thetaPri<(76*TMath::Pi()/180)");

double cont = 0;
double w = 0;
double val = 0;

// for (int i = 1; i < NB; i++){
//     cont = thet_0_6->GetBinContent(i);
//     w = thet_0_6->GetBinWidth(i);
//     val = cont / w;
//     thet_0_6->SetBinContent(i, val);
// }

// for (int i = 1; i < NB; i++){
//     cont = thet_43_47->GetBinContent(i);
//     w = thet_43_47->GetBinWidth(i);
//     val  = cont / w;
//     thet_43_47->SetBinContent(i, val);
// }

// for (int i = 1; i < NB; i++){
//     cont = thet_73_77->GetBinContent(i);
//     w = thet_73_77->GetBinWidth(i);
//     val  = cont / w;
//     thet_73_77->SetBinContent(i, val);
// }

double t0 = 0; // Degrees
double t45 = 45;
double t75 = 75;

// TF1 *Smith0 = new TF1("Distribuciones de Smith-Duller", Smith_Dull_Log, 1, 5, 1);
TF1 *Smith0 = new TF1("Distribuciones de Smith-Duller", Smith_Dull, tlow, thi, 1);
Smith0->SetParameter(0, t0);
Smith0->GetXaxis()->SetTitle("Energ#acute{i}a(MeV)");
double max0= Smith0->GetMaximum();
Smith0->SetLineColor(1);

TF1 *Smith45 = new TF1("Smith45", Smith_Dull, tlow, thi, 1);
Smith45->SetParameter(0, t45);
double max45= Smith45->GetMaximum();
Smith45->SetLineColor(2);

TF1 *Smith75 = new TF1("Smith75", Smith_Dull, tlow, thi, 1);
Smith75->SetParameter(0, t75);
double max75= Smith75->GetMaximum();
Smith75->SetLineColor(3);

double maxh0 = thet_0_6->GetMaximum();
double maxh45 = thet_43_47->GetMaximum();
double maxh75 = thet_73_77->GetMaximum();


/// ===== Se escalan a mano los histogramas al tamaño de la función ===== ///
thet_0_6->Scale(0.00000000024);
thet_43_47->Scale(0.000000000026);
thet_73_77->Scale(0.0000000000065);


// Create Canvas //
TCanvas *canv = new TCanvas("canv","", 2*700, 600);
canv->Divide(2,1);

canv->cd(1);
Smith0->Draw("L");
thet_0_6->Draw("hist same");

Smith45->Draw("L same");
thet_43_47->Draw("hist same");

Smith75->Draw("L same");
thet_73_77->Draw("hist same");

gPad->SetLogy(1);
gPad->SetLogx(1);


TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
// leg->SetTextAlign(11);
// leg->SetHeader("", "C");
leg->SetFillStyle(0);
leg->AddEntry(Smith0, "Curva te#acute{o}rica de #theta = 0", "l");
leg->AddEntry(Smith45, "Curva te#acute{o}rica para #theta = #frac{#pi}{4}", "l");
leg->AddEntry(Smith75, "Curva te#acute{o}rica para #theta = #frac{5#pi}{12}", "l");
leg->AddEntry(thet_0_6, "Datos Simulados", "f");
leg->Draw();

canv->Print("Dis_enpri.pdf");

}
