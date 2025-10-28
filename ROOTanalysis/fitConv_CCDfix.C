// G4Data from Smith-PlaneModel2 convolution anf fit 
// NOTE: Fitting values stored in LandauFittingValues file

bool doFit = true;
bool doPlot  = !doFit;
int rebinf = 1;
double KeVperbin = 6.666667; // 1000 KeV / 150 bins

//-------------------------------------
// The convolution function
//-------------------------------------

double f(double *x, double *par) {
    double pi = TMath::Pi();

    // Histograma de Ed simulado (200 bins de 0 a 700 MeV) 
    ifstream inp;
    inp.open("muons.dat");
    double p;
    int    Np    = 200;
    double Emin  = 0;              // Minimum energy
    double Emax  = 1000;           // Maximum energy
    double dE    = (Emax-Emin)/Np; // Energy interval
    double Ed;                     // Deposited energy
    //-------------------------------------------------------
    double Ea  = x[0];    // Function variable (E in PE)
    double f   = par[0];  // resolution at E0
    double Nmu = par[1];  // Muon normalization
    double Nb1 = par[2];  // Background 1 Normalizarion
	double Nb2 = par[3];  // Background 2 Normalizarion
    double a0  = par[4];  // MeV->PE, offset (PE)
    double a1  = par[5];  // MeV->PE linearity, (PE/MeV)
    double a2  = par[6];  // MeV->PE nonlinearity, (MeV^-1)
    double E0  = par[7];  // energy for which resol=f
    double ep1 = par[8];  // bkgd 1 decay const
    double ep2 = par[9];  // bkgd 2 decay const
    //-------------------------------------------------------
    // Visible energy in MeV. 
    // Inverse of Ea = a0 + a1*Ev/(1 + a2*Ev)
    double Ev = (Ea-a0)/(a1+a0*a2-a2*Ea); 

    // Convolution integral. Change units after
    double sig;
    double F = 0;        
    for (int i=1; i<=Np; i++) {
       inp >> p;         //- Read input-file line
       Ed  = (i-0.5)*dE; //- value at middle of bin
       sig = f*sqrt(Ev*E0);
       F += (Nmu*p + Nb1*(1./ep1)*exp(-Ed/ep1) + Nb2*(1./ep2)*exp(-Ed/ep2) )
           *(1./sqrt(2*pi*(pow(sig,2))))*exp(-(pow((Ev-Ed),2))/(2*(pow(sig,2))))*dE;
    } //for i
    double dEadEv = a1/pow(1 + a2*Ev,2);
    double z      = 1.0*F/dEadEv;  //change units to PE

    inp.close();
    return z * rebinf * KeVperbin;
}  // function f

//-------------------------------------
// Fit convolution function
//-------------------------------------
void fitConv_CCDfix() {

  //------------- Style --------------
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //----------------------------------

   // Histograma de Ed simulado (200 bins de 0 a maxEd MeV) 
   ifstream inpf;
   inpf.open("muons.dat");
   double pv;
   int    nbins = 200;
   double Emax  = 1000;
   double  dE   = Emax/nbins;


   // Histograma de energías depositadas Ed
   TH1F *hmuons = new TH1F("hmuons","",nbins,0,Emax);
   hmuons->GetXaxis()->SetTitle("Energy (KeV)");
   for (int i=1; i<=nbins; i++) {
       inpf >> pv;    // Reading input-file line
       hmuons->SetBinContent(i,pv);
   } //for i
   cout << "Integral hmuons: " << hmuons->Integral(1,nbins,"width") << "." << endl;

   TCanvas *c0 = new TCanvas("c0", "c0", 900,700);
   c0->SetGrid();
//    hmuons->Scale(0.014027296);
   hmuons->Draw();
   c0->Print("hmuons.pdf");

   //--CCD histogram
   //TFile *file = new TFile("ccm-data/bcm_t_nh_pe.root");
   TFile *file = new TFile("ccdhisto.root");
   //TFile *file = new TFile("ccm-data/from-Mayank/beamON_preBeam_bcm_nhits_previousEvent_selectingCosmicMuons.root");
   TH1F *hex = (TH1F*) file->FindObjectAny("edep_icn");
   cout<< hex->GetNbinsX() << endl;
   cout << "Integral Prompt_Energy: " << hex->Integral() << " entries (before rebin)." << endl;
//    hex->Rebin(rebinf); 
   cout << "Integral Prompt_Energy: " << hex->Integral() << " entries (after rebin)." << endl;
   TLatex *lat = new TLatex();
   lat->SetNDC();

   TCanvas *cx = new TCanvas("cx","testx", 900,700);
   cx->SetGrid();
   cx->cd();
   hex->GetXaxis()->SetRangeUser(0,1000);
   hex->SetMaximum(600);
   hex->Draw("hist");

   TF1 *f1 = new TF1();

//    hex->SetTitle("Energy histogram");
//    hex->SetXTitle("Energy (KeV)");
// //    hex->GetYaxis()->SetTitleOffset(1.1);
//    hex->SetYTitle("# Muons");
//    hex->SetLineColor(3);


   TCanvas *ct = new TCanvas("ct","test", 900,700);
   ct->SetGrid();
   ct->cd();
   hex->GetXaxis()->SetRangeUser(0,1000);
   hex->SetMaximum(300);
   hex->Draw("hist");


   // Data and Sim times (300x529 px)
//    double I0sim  = 101.2;
//    double nmusim = 113512; //1000000 simulados en total;
//    double Tsim   = 20969064.34; //sec 1M /  0.047689376 s^-1
//    double T      = 984670.02; //sec 766 * 1285.47 s
//    double eff = 1.0;

	// Data and Sim times 250x529 px)
 	double I0sim  = 101.2;
   	double nmusim = 113512; //1000000 simulados en total;
   	double Tsim   = 20969064.34; //sec 1M /  0.047689376 s^-1
   	double T      = 825914; //s
   	double eff = 1.0;

	TCanvas *c1 = new TCanvas("c1", "", 900,700);
	c1->SetGrid();

   if (doFit){
		cout<<"Hola"<<endl;
		int np = 10;           // Number of parameters
	double xm = 50;     // xmin to fit
	double xM = 700;    // xmax to fit

	// TF1 *f1 = new TF1("f1", f, xm, xM, n);
	f1 = new TF1("f1", f, xm, xM, np);

		// prebeam no muon selection
	double p0 = 0.05;    // Resolution
	double p1 = 4500.;    // Muon normalization
	double p2 = 7000.;    // Background 1
	double p3 = 4050.;    // Background 2
	double p4 = 0.00000e+00;    // PE offset
	double p5 = 1.0;    // No cambio de Unidades // por ahora
	double p6 = 0.00000e-04;    // PE scale (quadratic)
	double p7 = 220;    // E0
	double p8 = 70;    // exponente bkgd 1
	double p9 = 650;    // exponente bkgd 2


	f1->SetParameter(0, p0);    // Resolution
	f1->SetParameter(1, p1);    // Normalizing constant
	f1->SetParameter(2, p2);    // Background 1
	f1->SetParameter(3, p3);    // Background 2
	f1->FixParameter(4, p4);    // Energy scale offset
	f1->FixParameter(5, p5);    // Energy scale (linear)
	f1->FixParameter(6, p6);    // Energy scale (non-linear)
	f1->FixParameter(7, p7);    // E0
	f1->SetParameter(8, p8);    // exponente bkgd 1
	f1->SetParameter(9, p9);    // exponente bkgd 1

	// hex->Draw("hist");
	hex->Fit("f1", "R");   // Fitting in the specified range
	// hex->Draw("hist same");

	double r    = f1->GetParameter(0);
	double er   = f1->GetParError(0);
	double nmu  = f1->GetParameter(1);
	double enmu = f1->GetParError(1);
	double nb1  = f1->GetParameter(2);
	double enb1 = f1->GetParError(2);
	double nb2  = f1->GetParameter(3);
	double enb2 = f1->GetParError(3);
	double a0   = f1->GetParameter(4);
	double ea0  = f1->GetParError(4);
	double a1   = f1->GetParameter(5);
	double ea1  = f1->GetParError(5);
	double a2   = f1->GetParameter(6);
	double ea2  = f1->GetParError(6);
	double e0   = f1->GetParameter(7);
	double ee0  = f1->GetParError(7);
	double ep1  = f1->GetParameter(8);
	double eep1 = f1->GetParError(8);
	double ep2  = f1->GetParameter(9);
	double eep2 = f1->GetParError(9);

	double chi2 = f1->GetChisquare();
	int    ndf  = f1->GetNDF();
	double prob = TMath::Prob(chi2,ndf);

	// // Calculate I_0
	double I0  = I0sim*(nmu/nmusim)*(Tsim/T)*(1./eff);
	double eI0 = I0sim*(enmu/nmusim)*(Tsim/T)*(1./eff);

	TF1 *fm = new TF1("fm", f, xm, xM, np);
	fm->SetParameters(r,nmu,0*nb1, 0*nb2, a0,a1,a2,e0,ep1, ep2);
	fm->SetLineColor(4);

	TF1 *fb1 = new TF1("fb1", f, xm, xM, np);
	fb1->SetParameters(r,0*nmu,nb1, 0*nb2, a0,a1,a2,e0,ep1, ep2);
	fb1->SetLineColor(12);
	fb1->SetLineStyle(3);
	fb1->SetLineWidth(1);

	TF1 *fb2 = new TF1("fb2", f, xm, xM, np);
	fb2->SetParameters(r,0*nmu,0*nb1, nb2, a0,a1,a2,e0,ep1, ep2);
	fb2->SetLineColor(15);
	fb2->SetLineStyle(3);
	fb2->SetLineWidth(1);

	c1->cd();
	hex->SetMaximum(600);
	hex->SetLineColor(3);
	hex->Draw();
	// h->GetXaxis()->SetRangeUser(0,1000);
	f1->Draw("");
	fb1->Draw("same");
	fm->Draw("same");
	fb2->Draw("same");
	hex->Draw("hist same");

	lat->SetTextFont(42);
	lat->SetTextSize(0.034);
	lat->DrawLatex(0.15,0.85,Form("I^{CCD}_{0} = (%3.1f #pm %3.1f) m^{-2} s^{-1} sr^{-1}",I0,eI0));
	lat->DrawLatex(0.15,0.80,Form("N_{#mu} = (%6.0f #pm %4.0f)",nmu,enmu));
	lat->DrawLatex(0.15,0.76,Form("Resolution f = (%5.3f #pm %5.3f)",r,er));
	lat->DrawLatex(0.15,0.72,Form("N_{b1} = (%6.0f #pm %4.0f)",nb1,enb1));
	lat->DrawLatex(0.15,0.68,Form("#epsilon_1 = (%5.3f #pm %5.3f)",ep1,eep1));
	lat->DrawLatex(0.15,0.64,Form("N_{b2} = (%6.0f #pm %4.0f)",nb2,enb2));
	lat->DrawLatex(0.15,0.60,Form("#epsilon_2 = (%5.3f #pm %5.3f)",ep2,eep2));
	lat->DrawLatex(0.15,0.56,Form("a_{0} = (%5.3f #pm %5.3f) PE",a0,ea0));
	lat->DrawLatex(0.15,0.52,Form("a_{1} = (%5.3f #pm %5.3f) PE/MeV",a1,ea1));
	lat->DrawLatex(0.15,0.48,Form("a_{2} = (%6.5f #pm %6.5f) MeV^{-1}",a2,ea2));

	lat->DrawLatex(0.15,0.42,Form("#chi^{2}/ndf = %4.2f/%d",chi2,ndf));
	lat->DrawLatex(0.15,0.38,Form("Prob(#chi^{2}) = %6.4f",prob));

	TLegend *l = new TLegend(0.6, 0.6, 0.9, 0.75);
	l->SetTextSize(0.03);
	l->AddEntry(hex, "Data", "lp");
	l->AddEntry(f1, "Convolution Fit", "lp");
	l->Draw();
	// hex->GetXaxis()->SetRangeUser(0,1);

	c1->Print("ConvNonlinear.pdf");

        double funcInt = f1->Integral(xm,xM);
	cout << "Integral fitConv = " << funcInt << " muons" << endl;

 } //if doFit

} // fitConv
