// G4Data from Smith-PlaneModel2 convolution anf fit 
// NOTE: Fitting values stored in LandauFittingValues file

bool doFit = true;
bool doPlot  = !doFit;
int rebinf = 1;
double KeVperbin = 6.666667; // 1000 KeV / 150 bins

//-------------------------------------
// The convolution function
//-------------------------------------

// double f(double *x, double *par) {
//     double pi = TMath::Pi();

//     // Histograma de Ed simulado (200 bins de 0 a 700 MeV) 
//     ifstream inp;
//     inp.open("muons.dat");
//     double p;
//     int    Np    = 150;
//     double Emin  = 0;              // Minimum energy
//     double Emax  = 1500;           // Maximum energy
//     double dE    = (Emax-Emin)/Np; // Energy interval
//     double Ed;                     // Deposited energy
//     //-------------------------------------------------------
//     double Ea  = x[0];    // Function variable (E in PE)
//     double f   = par[0];  // resolution at E0
//     double Nmu = par[1];  // Muon normalization
//     double Nb1 = par[2];  // Background 1 Normalizarion
// 	double Nb2 = par[3];  // Background 2 Normalizarion
//     double a0  = par[4];  // MeV->PE, offset (PE)
//     double a1  = par[5];  // MeV->PE linearity, (PE/MeV)
//     double a2  = par[6];  // MeV->PE nonlinearity, (MeV^-1)
//     double E0  = par[7];  // energy for which resol=f
//     double ep1 = par[8];  // bkgd 1 decay const
//     double ep2 = par[9];  // bkgd 2 decay const
//     //-------------------------------------------------------
//     // Visible energy in MeV. 
//     // Inverse of Ea = a0 + a1*Ev/(1 + a2*Ev)
//     double Ev = (Ea-a0)/(a1+a0*a2-a2*Ea); 

//     // Convolution integral. Change units after
//     double sig;
//     double F = 0;        
//     for (int i=1; i<=Np; i++) {
//        inp >> p;         //- Read input-file line
//        Ed  = (i-0.5)*dE; //- value at middle of bin
//        sig = f*sqrt(Ev*E0);
//        F += (Nmu*p + Nb1*(1./ep1)*exp(-Ed/ep1) + Nb2*(1./ep2)*exp(-Ed/ep2) )
//            *(1./sqrt(2*pi*(pow(sig,2))))*exp(-(pow((Ev-Ed),2))/(2*(pow(sig,2))))*dE;
//     } //for i
//     double dEadEv = a1/pow(1 + a2*Ev,2);
//     double z      = 1.0*F/dEadEv;  //change units to PE

//     inp.close();
//     return z * rebinf * KeVperbin;
// }  // function f


double f(double *x, double *par) {
    double pi = TMath::Pi();

    // Histograma de Ed simulado (200 bins de 0 a 700 MeV) 
    ifstream inp;
    inp.open("muons.dat");
    double p;
    int    Np    = 150;
    double Emin  = 0;              // Minimum energy
    double Emax  = 1500;           // Maximum energy
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
	TGraph *fgr = new TGraph(Np);
    double sig;
    double F = 0;        
    for (int i=1; i<=Np; i++) {
       inp >> p;         //- Read input-file line
       Ed  = (i-0.5)*dE; //- value at middle of bin
	   F = (Nmu*p + Nb1*(1./ep1)*exp(-Ed/ep1) + Nb2*(1./ep2)*exp(-Ed/ep2) );
	   fgr->AddPoint(Ed, F);
    //    sig = f*sqrt(Ev*E0);
    } //for i
	F = fgr->Eval(Ev);
    double dEadEv = a1/pow(1 + a2*Ev,2);
    double z      = 1.0*F/dEadEv;  //change units to PE

    inp.close();
    return z * rebinf * KeVperbin;
}  // function f


//-------------------------------------
// Fit convolution function
//-------------------------------------
void fitConv_CCD_nores() {

  //------------- Style --------------
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //----------------------------------

   // Histograma de Ed simulado (200 bins de 0 a maxEd MeV) 
   ifstream inpf;
   inpf.open("muons.dat");
   double pv;
   int    nbins = 150;
   double Emax  = 1500;
   double  dE   = Emax/nbins;


   // Histograma de energías depositadas Ed
   TH1F *hmuons = new TH1F("hmuons","",nbins,0,Emax);
   hmuons->GetXaxis()->SetTitle("Energy (KeV)");
   for (int i=1; i<=nbins; i++) {
       inpf >> pv;    // Reading input-file line
       hmuons->SetBinContent(i,pv);
   } //for i
   cout << "Integral hmuons: " << hmuons->Integral(1,nbins,"width") << "." << endl;

   TCanvas *c0 = new TCanvas("c0", "", 900,700);
   c0->SetGrid();
//    hmuons->Scale(0.014027296);
   hmuons->Draw();
   c0->Print("hmuons.pdf");

   TCanvas *c1 = new TCanvas("c1", "", 900,700);
   c1->SetGrid();

   //--CCD histogram
   //TFile *file = new TFile("ccm-data/bcm_t_nh_pe.root");
   TFile *file = new TFile("ccdhisto.root");
   //TFile *file = new TFile("ccm-data/from-Mayank/beamON_preBeam_bcm_nhits_previousEvent_selectingCosmicMuons.root");
   TH1F *hex = (TH1F*) file->FindObjectAny("edep_icn");
//    TH1F *hex = (TH1F*) file->FindObjectAny("edep_conn");
   cout<< hex->GetNbinsX() << endl;
   cout << "Integral Prompt_Energy: " << hex->Integral() << " entries (before rebin)." << endl;
   hex->Rebin(rebinf); 
   cout << "Integral Prompt_Energy: " << hex->Integral() << " entries (after rebin)." << endl;
   TLatex *lat = new TLatex();
   lat->SetNDC();

   TF1 *f1 = new TF1();

//    hex->SetMaximum((rebinf/10)*8*hex->GetMaximum());
//    hex->SetMaximum(20000);
   hex->SetTitle("Energy histogram");
   hex->SetXTitle("Energy (KeV)");
   hex->GetYaxis()->SetTitleOffset(1.1);
   hex->SetYTitle("# Muons");
   hex->SetLineColor(3);
   hex->SetLineWidth(2);


   TCanvas *ct = new TCanvas("ct","test", 900,700);
   ct->cd();
   hex->GetXaxis()->SetLimits(0,1000);
   //hex->GetXaxis()->SetRangeUser(0,1000);
   hex->Draw("hist");


   // Data and Sim times (300x529 px)
//    double I0sim  = 101.2;
//    double nmusim = 113512; //1000000 simulados en total;
//    double Tsim   = 20969064.34; //sec 1M /  0.047689376 s^-1
//    double T      = 984670.02; //sec 766 * 1285.47 s
//    double eff = 1.0;

	// // Data and Sim times 250x529 px)
 	// double I0sim  = 101.2;
   	// double nmusim = 424665; //4000000 simulados en total;
   	// double Tsim   = 8.38761e+07; //sec 4M /  0.0476894 s^-1
   	// double T      = 825914; // Tiempo experimental s
   	// double eff = 1.0;

	double I0sim  = 101.2;
   	double nmusim = 105247; //1000000 simulados en total;
   	double Tsim   = 20969020; //sec 4M /  0.0476894 s^-1
   	double T      = 825914; // Tiempo experimental s
   	double eff = 1.0;

	// // Data and Sim CONNIE 420x1022 px)
	// double I0sim  = 101.2;
	// double nmusim = 268880; //1000000 simulados en total;
	// double Tsim   = 16325369.93; //sec 1M /  0.061254354 s^-1
	// double T      = 2434014 * 1; //s
	// double eff = 1.0;

	// Data and Sim CONNIE 420x700 px)
	// double I0sim  = 101.2;
	// double nmusim = 185579; //4000000 simulados en total;
	// double Tsim   = 16325357.9; //sec 1M /  0.061254354 s^-1
	// double T      = 1654349.759; //s
	// // double T      = 1067132; //s
	// double eff = 1.0;



   if (doFit){

	int np = 10;           // Number of parameters
	double xm = 50;     // xmin to fit
	double xM = 700;    // xmax to fit

	// TF1 *f1 = new TF1("f1", f, xm, xM, n);
	f1 = new TF1("f1", f, xm, xM, np);

		// prebeam no muon selection
		double p0 = 0.0;    // Resolution
		double p1 = 5000.;    // Muon normalization
		double p2 = 8075.;    // Background 1
		double p3 = 6178.;    // Background 2
		double p4 = 0.00000e+00;    // PE offset
		double p5 = 1.0;    // No cambio de Unidades // por ahora
		double p6 = 0.00000e-04;    // PE scale (quadratic)
		double p7 = 220;    // E0
		double p8 = 56.113;    // exponente bkgd 1
		double p9 = 332.344;    // exponente bkgd 2

	// ===== Parámetros para CONNIE ==== //
	// // double p0 = 0.031;    // Resolution
	// double p0 = 0.01;    // Resolution
	// double p1 = 20000.;    // Muon normalization
	// double p2 = 2075.;    // Background 1
	// double p3 = 3178.;    // Background 2
	// double p4 = 0.00000e+00;    // PE offset
	// double p5 = 1.0;    // No cambio de Unidades // por ahora
	// double p6 = 0.00000e-04;    // PE scale (quadratic)
	// double p7 = 220;    // E0
	// double p8 = 56.113;    // exponente bkgd 1
	// double p9 = 332.344;    // exponente bkgd 2


	f1->FixParameter(0, p0);    // Resolution
	f1->SetParameter(1, p1);    // Normalizing constant
	f1->SetParameter(2, p2);    // Background 1
	f1->SetParameter(3, p3);    // Background 2
	f1->FixParameter(4, p4);    // Energy scale offset
	f1->SetParameter(5, p5);    // Energy scale (linear)
	f1->FixParameter(6, p6);    // Energy scale (non-linear)
	f1->FixParameter(7, p7);    // E0
	f1->SetParameter(8, p8);    // exponente bkgd 1
	f1->SetParameter(9, p9);    // exponente bkgd 1

	hex->Fit("f1","R");   // Fitting in the specified range

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
	fb1->SetLineColor(6);
	fb1->SetLineStyle(3);
	// fb1->SetLineWidth(1);

	TF1 *fb2 = new TF1("fb2", f, xm, xM, np);
	fb2->SetParameters(r,0*nmu,0*nb1, nb2, a0,a1,a2,e0,ep1, ep2);
	fb2->SetLineColor(7);
	fb2->SetLineStyle(3);
	// fb2->SetLineWidth(1);

	c1->cd();
	hex->SetMaximum(800);
	// hex->SetLineColor(3);
	// f1->GetXaxis()->SetRangeUser(0,1000);
	// h->Draw("hist");
	hex->Draw("hist");
	f1->Draw("same");
	fb1->Draw("same");
	fm->Draw("same");
	fb2->Draw("same");

	lat->SetTextFont(42);
	lat->SetTextSize(0.034);
	lat->DrawLatex(0.15,0.85,Form("I^{CCM}_{0} = (%3.1f #pm %3.1f) m^{-2} s^{-1} sr^{-1}",I0,eI0));
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

	// c1->Print("ConvNonlinear_CONNIE_420x1022.pdf");
	c1->Print("ConvNonlinear_NORES_Log_ICN.pdf");
	// c1->Print("ConvNonlinear_NORES_Lin_ICN.pdf");

        double funcInt = f1->Integral(xm,xM);
	cout << "Integral fitConv = " << funcInt << " muons" << endl;

 } //if doFit


//-------------------------------------
// Plot convolution function
//-------------------------------------
  if (doPlot) {

	
    //     // prebeam no muon selection
	// double p0 = 0.029;    // Resolution
	// double p1 = 7998.;    // Muon normalization
	// double p2 = 15036.;    // Background 1
	// double p3 = 12287.;    // Background 2
	// double p4 = 0.00000e+00;    // PE offset
	// double p5 = 1.0;    // No cambio de Unidades // por ahora
	// double p6 = 0.00000e-04;    // PE scale (quadratic)
	// double p7 = 220;    // E0
	// double p8 = 58.533;    // exponente bkgd 1
	// double p9 = 452.959;    // exponente bkgd 2

	// ========== CCD size: 529x300 ============== //
	// double p0 = 0.031;    // Resolution
	// double p1 = 5710.;    // Muon normalization
	// double p2 = 10177.;    // Background 1
	// double p3 = 7823.;    // Background 2
	// double p4 = 0.00000e+00;    // PE offset
	// double p5 = 1.0;    // No cambio de Unidades // por ahora
	// double p6 = 0.00000e-04;    // PE scale (quadratic)
	// double p7 = 220;    // E0
	// double p8 = 61.118;    // exponente bkgd 1
	// double p9 = 390.377;    // exponente bkgd 2

	// // ========== CCD size: 529x250 ============== //
	double p0 = 0.0;    // Resolution
	double p1 = 5977.;    // Muon normalization
	double p2 = 9582.;    // Background 1
	double p3 = 2791;    // Background 2
	double p4 = 0.00000e+00;    // PE offset
	double p5 = 1.0;    // No cambio de Unidades // por ahora
	double p6 = 0.00000e-04;    // PE scale (quadratic)
	double p7 = 220;    // E0
	double p8 = 66.897;    // exponente bkgd 1
	double p9 = 397.646;    // exponente bkgd 2

	// ========== CONNIE size: 420x1022 ============== //
	// double p0 = 0.031;    // Resolution
	// double p1 = 15500.;    // Muon normalization
	// double p2 = 2075.;    // Background 1
	// double p3 = 3178.;    // Background 2
	// double p4 = 0.00000e+00;    // PE offset
	// double p5 = 1.0;    // No cambio de Unidades // por ahora
	// double p6 = 0.00000e-04;    // PE scale (quadratic)
	// double p7 = 220;    // E0
	// double p8 = 56.113;    // exponente bkgd 1
	// double p9 = 332.344;    // exponente bkgd 2



	int np = 10; //number of parameters
	double xm = 50;
	double xM = 700;

	f1 = new TF1("f1", f, xm, xM, np);
	f1->SetParameters(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9);
	f1->SetLineColor(2);

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

	// Calculate I_0
	double I0  = I0sim*(nmu/nmusim)*(Tsim/T)*(1./eff);
	double eI0 = I0sim*(enmu/nmusim)*(Tsim/T)*(1./eff);


	TF1 *fm = new TF1("fm", f, xm, xM, np);
	fm->SetParameters(p0, p1, 0*p2, 0*p3, p4, p5, p6, p7, p8, p9);
	fm->SetLineColor(4);

	TF1 *fb1 = new TF1("fb1", f, xm, xM, np);
	fb1->SetParameters(p0, 0*p1, p2, 0*p3, p4, p5, p6, p7, p8, p9);
	fb1->SetLineColor(6);
	fb1->SetLineStyle(2);

	TF1 *fb2 = new TF1("fb2", f, xm, xM, np);
	fb2->SetParameters(p0, 0*p1, 0*p2, p3, p4, p5, p6, p7, p8, p9);
	fb2->SetLineColor(7);
	fb2->SetLineStyle(2);


	c1->cd();
	//c1->SetLogy(1);
	hex->SetMaximum(800*rebinf);
	// h->SetMinimum(0);
	hex->Draw();
	// h->GetXaxis()->SetRangeUser(0,700);

	double xmax = f1->GetMaximumX();
	cout << xmax << endl;

	f1->Draw("l same");
	fm->Draw("same");
    fb1->Draw("same");
	fb2->Draw("same");

	lat->SetTextFont(42);
	lat->SetTextSize(0.034);
	lat->DrawLatex(0.15,0.85,Form("I^{CCD}_{0} = %3.1f m^{-2} s^{-1} sr^{-1}",I0));
	lat->DrawLatex(0.15,0.80,Form("N_{#mu} = %6.0f ",nmu));
	lat->DrawLatex(0.15,0.76,Form("Resolution f = %5.3f ",r));
	lat->DrawLatex(0.15,0.72,Form("N_{b1} = %6.0f ",nb1));
	lat->DrawLatex(0.15,0.68,Form("#epsilon_1 = %5.2f ",ep1));
	lat->DrawLatex(0.15,0.64,Form("N_{b2} = %6.0f ",nb2));
	lat->DrawLatex(0.15,0.60,Form("#epsilon_2 = %5.2f ",ep2));
	lat->DrawLatex(0.15,0.56,Form("a_{0} = %5.2f PE",a0));
	lat->DrawLatex(0.15,0.52,Form("a_{1} = %5.2f  PE/MeV",a1));
	lat->DrawLatex(0.15,0.48,Form("a_{2} = %6.5f  MeV^{-1}",a2));
	//lat->DrawLatex(0.15,0.51,Form("#chi^{2}/ndf = %4.2f/%d",chi2,ndf));
	//lat->DrawLatex(0.15,0.47,Form("Prob(#chi^{2}) = %6.4f",prob));


	TLegend *l = new TLegend(0.6, 0.6, 0.9, 0.75);
	l->SetTextSize(0.03);
	l->AddEntry(hex, "Data", "lp");
	l->AddEntry(f1, "Convolution model", "lp");
	l->Draw();

	c1->Print("ConvNonlinear_NORES.pdf");

 } // ifdoPlot

	// nonlinearity function
	TF1 *nonlin = new TF1("nonlin","[0]+[1]*x/(1+[2]*x)",0,Emax); 
        double a0 = f1->GetParameter(4);
        double a1 = f1->GetParameter(5);
        double a2 = f1->GetParameter(6);
	nonlin->SetParameters(a0,a1,a2);
	nonlin->GetXaxis()->SetTitle("Energy (MeV)");
	nonlin->GetYaxis()->SetTitleOffset(1.5);
	nonlin->GetYaxis()->SetTitle("Energy (PE)");
	// linear part
	TF1 *lin = new TF1("lin","[0]+[1]*x",0,Emax); 
	lin->SetParameters(a0,a1);
	lin->GetXaxis()->SetTitle("Energy (MeV)");
	lin->GetYaxis()->SetTitle("Energy (PE)");

	TCanvas *canv2 = new TCanvas("canv2","",400,300);
	canv2->cd();
	nonlin->Draw("");
	lin->SetLineColor(1);
	lin->Draw("same");
        double xx = 255.5;
        double yy = nonlin->Eval(xx);
        TLine *l0 = new TLine(0,0,Emax,0); l0->Draw("same");
        TLine *lh = new TLine(0,yy,xx,yy); lh->Draw("same");
        TLine *lv = new TLine(xx,0,xx,yy); lv->Draw("same");
        lat->SetTextFont(42);
        lat->SetTextSize(0.035);
        lat->DrawLatex(0.6,0.50,"E_{PE} = a_{0} + #frac{a_{1} E_{vis}}{1 + a_{2} E_{vis}}");
        lat->DrawLatex(0.6,0.40,Form("a_{0} = %3.1f PE",a0));
        lat->DrawLatex(0.6,0.35,Form("a_{1} = %3.1f PE/MeV",a1));
        lat->DrawLatex(0.6,0.30,Form("a_{2} = %3.1e MeV^{-1}",a2));
	canv2->Print("MeVtoPE_fitConv.pdf");


} // fitConv
