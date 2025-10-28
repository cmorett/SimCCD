// G4Data from Smith-PlaneModel2 convolution anf fit 
// NOTE: Fitting values stored in LandauFittingValues file

bool doFit = true;
bool doPlot  = !doFit;
int rebinf = 200;

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
    double Emax  = 700;           // Maximum energy
    double dE    = (Emax-Emin)/Np; // Energy interval
    double Ed;                     // Deposited energy
    //-------------------------------------------------------
    double Ea  = x[0];    // Function variable (E in PE)
    double f   = par[0];  // resolution at E0
    double Nmu = par[1];  // resolution at E0
    double Nbg = par[2];  // Background Normalizarion
    double a0  = par[3];  // MeV->PE, offset (PE)
    double a1  = par[4];  // MeV->PE linearity, (PE/MeV)
    double a2  = par[5];  // MeV->PE nonlinearity, (MeV^-1)
    double E0  = par[6];  // energy for which resol=f
    double eps = par[7];  // bkgd decay const
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
       sig = f*sqrt(Ed*E0);
       F += (Nmu*rebinf*p + Nbg*rebinf*(1./eps)*exp(-Ed/eps) )
           *(1./sqrt(2*pi*(sig**2)))*exp(-((Ev-Ed)**2)/(2*(sig**2)))*dE;
    } //for i
    double dEadEv = a1/pow(1 + a2*Ev,2);
    double z      = 1.0*F/dEadEv;  //change units to PE

    inp.close();
    return z;
}  // function f

//-------------------------------------
// Fit convolution function
//-------------------------------------
void fitConv_ccm() {

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
   double Emax  = 700;
   double  dE   = Emax/nbins;


   // Histograma de energías depositadas Ed
   TH1F *hmuons = new TH1F("hmuons","",nbins,0,Emax);
   hmuons->GetXaxis()->SetTitle("Energy (MeV)");
   for (int i=1; i<=nbins; i++) {
       inpf >> pv;    // Reading input-file line
       hmuons->SetBinContent(i,pv);
   } //for i
   cout << "Integral hmuons: " << hmuons->Integral(1,nbins,"width") << "." << endl;

   TCanvas *c0 = new TCanvas("c0", "", 900,700);
   c0->SetGrid();
   hmuons->Draw();
   c0->Print("hmuons.pdf");

   TCanvas *c1 = new TCanvas("c1", "", 900,700);
   c1->SetGrid();

   //--CCM histogram
   //TFile *file = new TFile("ccm-data/bcm_t_nh_pe.root");
   TFile *file = new TFile("ccm-data/from-Mayank/beamON_preBeam_bcm_nhits_previousEvent.root");
   //TFile *file = new TFile("ccm-data/from-Mayank/beamON_preBeam_bcm_nhits_previousEvent_selectingCosmicMuons.root");
   TH1F *h = (TH1F*) file->FindObjectAny("Prompt_Energy");
   cout << "Integral Prompt_Energy: " << h->Integral() << " entries (before rebin)." << endl;
   h->Rebin(rebinf); 
   cout << "Integral Prompt_Energy: " << h->Integral() << " entries (after rebin)." << endl;
   TLatex *lat = new TLatex();
   lat->SetNDC();

   TF1 *f1 = new TF1();

   h->SetMaximum((rebinf/10)*8*h->GetMaximum());
   h->SetTitle("Energy histogram");
   h->SetXTitle("Energy (PE)");
   h->GetYaxis()->SetTitleOffset(1.1);
   h->SetYTitle("# Muons");
   h->SetLineColor(3);

   if (doFit){

	int n = 8;           // Number of parameters
	double xm = 2500;     // xmin to fit
	double xM = 40000;    // xmax to fit

	//TF1 *f1 = new TF1("f1", f, xm, xM, n);
	f1 = new TF1("f1", f, xm, xM, n);

         
        //prebeam NO muon selection
	double p0 = 0.70000e-01;    // Resolution
	double p1 = 6.50000e+03;    // Muon normalization
	double p2 = 3.10000e+03;    // Background
	double p3 = 0.00000e+00;    // PE offset
	double p4 = 0.70000e+02;    // PE scale (linear)
	double p5 = 0.00000e-04;    // PE scale (quadratic)
	double p6 = 2.55000e+02;    // E0
	double p7 = 1.70000e+01;    // exponente bkgd

        /*
        //- prebeam WITH muon selection
	double p0 = 2.0000e-01;    // Resolution
	double p1 = 3.20000e+03;    // Muon normalization
	double p2 = 2.39202e+03;    // Background
	double p3 = 0.00000e+00;    // PE offset
	double p4 = 1.72000e+02;    // PE scale (linear)
	double p5 = 4.80000e-03;    // PE scale (quadratic)
	double p6 = 2.60000e+02;    // E0
	double p7 = 0.30000e+01;    // exponente bkgd
        */

	f1->SetParameter(0, p0);    // Resolution
	f1->SetParameter(1, p1);    // Normalizing constant
	f1->SetParameter(2, p2);    // Background
	f1->FixParameter(3, p3);    // PE offset
	f1->SetParameter(4, p4);    // PE scale (linear)
	f1->SetParameter(5, p5);    // PE scale (quadratic)
	f1->FixParameter(6, p6);    // E0
	f1->SetParameter(7, p7);    // exponente bkgd

	h->Fit("f1","R");   // Fitting in the specified range

	double r    = f1->GetParameter(0);
	double er   = f1->GetParError(0);
	double nmu  = f1->GetParameter(1);
	double enmu = f1->GetParError(1);
	double nbg  = f1->GetParameter(2);
	double enbg = f1->GetParError(2);
	double a0   = f1->GetParameter(3);
	double ea0  = f1->GetParError(3);
	double a1   = f1->GetParameter(4);
	double ea1  = f1->GetParError(4);
	double a2   = f1->GetParameter(5);
	double ea2  = f1->GetParError(5);
	double e0   = f1->GetParameter(6);
	double ee0  = f1->GetParError(6);
	double eps  = f1->GetParameter(7);
	double eeps = f1->GetParError(7);
	double chi2 = f1->GetChisquare();
	int    ndf  = f1->GetNDF();
	double prob = TMath::Prob(chi2,ndf);

	// Calculate I_0
	double I0sim  = 70;
	double nmusim = 308015;//393500;
	double Tsim   = 472.3; //sec
	double T      = 6.0; //sec
	double eff    = 1.0;
	double I0  = I0sim*(nmu/nmusim)*(Tsim/T)*(1./eff);
	double eI0 = I0sim*(enmu/nmusim)*(Tsim/T)*(1./eff);

	fm = new TF1("fm", f, xm, xM, n);
	fm->SetParameters(r,nmu,0,a0,a1,a2,e0,eps);
	fm->SetLineColor(4);

	TF1 *fb = new TF1("fb", f, xm, xM, n);
	fb->SetParameters(r,0*nmu,nbg,a0,a1,a2,e0,eps);
        fb->SetLineColor(13);
        fb->SetLineStyle(3);
        fb->SetLineWidth(1);

	//c1->SetLogy(1);
	h->SetMaximum(120);
	h->Draw();
	fm->Draw("same");
	fb->Draw("same");

	lat->SetTextFont(42);
	lat->SetTextSize(0.034);
	lat->DrawLatex(0.15,0.85,Form("I^{CCM}_{0} = (%3.1f #pm %3.1f) m^{-2} s^{-1} sr^{-1}",I0,eI0));
	lat->DrawLatex(0.15,0.80,Form("N_{#mu} = (%6.0f #pm %4.0f)",nmu,enmu));
	lat->DrawLatex(0.15,0.76,Form("Resolution f = (%5.3f #pm %5.3f)",r,er));
	lat->DrawLatex(0.15,0.72,Form("N_{bg} = (%6.0f #pm %4.0f)",nbg,enbg));
	lat->DrawLatex(0.15,0.68,Form("#epsilon = (%5.3f #pm %5.3f)",eps,eeps));
	lat->DrawLatex(0.15,0.64,Form("a_{0} = (%5.3f #pm %5.3f) PE",a0,ea0));
	lat->DrawLatex(0.15,0.60,Form("a_{1} = (%5.3f #pm %5.3f) PE/MeV",a1,ea1));
	lat->DrawLatex(0.15,0.56,Form("a_{2} = (%6.5f #pm %6.5f) MeV^{-1}",a2,ea2));
	lat->DrawLatex(0.15,0.51,Form("#chi^{2}/ndf = %4.2f/%d",chi2,ndf));
	lat->DrawLatex(0.15,0.47,Form("Prob(#chi^{2}) = %6.4f",prob));

	TLegend *l = new TLegend(0.6, 0.6, 0.9, 0.75);
	l->SetTextSize(0.03);
	l->AddEntry(h, "Data", "lp");
	l->AddEntry(f1, "Convolution Fit", "lp");
	l->Draw();
	h->GetXaxis()->SetRangeUser(0,40000);

	c1->Print("ConvNonlinear.pdf");

        double funcInt = f1->Integral(xm,xM);
	cout << "Integral fitConv = " << funcInt << " muons" << endl;

 } //if doFit


//-------------------------------------
// Plot convolution function
//-------------------------------------
  if (doPlot) {

	
        // prebeam no muon selection
	double p0 = 0.70000e-01;    // Resolution
	double p1 = 6.50000e+03;    // Muon normalization
	double p2 = 3.10000e+03;    // Background
	double p3 = 0.00000e+00;    // PE offset
	double p4 = 0.70000e+02;    // PE scale (linear)
	double p5 = 0.00000e-04;    // PE scale (quadratic)
	double p6 = 2.55000e+02;    // E0
	double p7 = 1.70000e+01;    // exponente bkgd

		

        /*
	// prebeam with muon selection
	double p0 = 2.0000e-01;    // Resolution
	double p1 = 3.20000e+03;    // Muon normalization
	double p2 = 2.00000e+03;    // Background
	double p3 = 0.00000e+00;    // PE offset
	double p4 = 1.35000e+02;    // PE scale (linear)
	double p5 = 2.87000e-03;    // PE scale (quadratic)
	double p6 = 2.60000e+02;    // E0
	double p7 = 0.30000e+01;    // exponente bkgd
        */	

	int np = 8;
	double xm = 2500;
	double xM = 40000;

	f1 = new TF1("f1", f, xm, xM, np);
	f1->SetParameters(p0, p1, p2, p3, p4, p5, p6, p7);
	f1->SetLineColor(2);

	double r    = f1->GetParameter(0);
	double er   = f1->GetParError(0);
	double nmu  = f1->GetParameter(1);
	double enmu = f1->GetParError(1);
	double nbg  = f1->GetParameter(2);
	double enbg = f1->GetParError(2);
	double a0   = f1->GetParameter(3);
	double ea0  = f1->GetParError(3);
	double a1   = f1->GetParameter(4);
	double ea1  = f1->GetParError(4);
	double a2   = f1->GetParameter(5);
	double ea2  = f1->GetParError(5);
	double e0   = f1->GetParameter(6);
	double ee0  = f1->GetParError(6);
	double eps  = f1->GetParameter(7);
	double eeps = f1->GetParError(7);
	double chi2 = f1->GetChisquare();
	int    ndf  = f1->GetNDF();
	double prob = TMath::Prob(chi2,ndf);

	// Calculate I_0
	double I0sim  = 70;
	double nmusim = 308015;//393500;
	double Tsim   = 472.3; //sec
	double T      = 6.0; //sec
	double eff    = 1.0;
	double I0  = I0sim*(nmu/nmusim)*(Tsim/T)*(1./eff);
	double eI0 = I0sim*(enmu/nmusim)*(Tsim/T)*(1./eff);

	fm = new TF1("fm", f, xm, xM, np);
	fm->SetParameters(p0, p1, 0*p2, p3, p4, p5, p6, p7);
	fm->SetLineColor(4);

	TF1 *fbg = new TF1("fbg", f, xm, xM, np);
	fbg->SetParameters(p0, 0*p1, p2, p3, p4, p5, p6, p7);
        fbg->SetLineColor(6);


	c1->cd();
	//c1->SetLogy(1);
	h->SetMaximum(120);
	h->Draw();
	h->GetXaxis()->SetRangeUser(0,40000);

	double xmax = f1->GetMaximumX();
	cout << xmax << endl;

	f1->Draw("l same");
        fm->Draw("same");
        fbg->Draw("same");

	lat->SetTextFont(42);
	lat->SetTextSize(0.034);
	lat->DrawLatex(0.15,0.85,Form("I^{CCM}_{0} = %3.1f m^{-2} s^{-1} sr^{-1}",I0));
	lat->DrawLatex(0.15,0.80,Form("N_{#mu} = %6.0f ",nmu));
	lat->DrawLatex(0.15,0.76,Form("Resolution f = %5.3f ",r));
	lat->DrawLatex(0.15,0.72,Form("N_{bg} = %6.0f ",nbg));
	lat->DrawLatex(0.15,0.68,Form("#epsilon = %5.2f ",eps));
	lat->DrawLatex(0.15,0.64,Form("a_{0} = %5.2f PE",a0));
	lat->DrawLatex(0.15,0.60,Form("a_{1} = %5.2f  PE/MeV",a1));
	lat->DrawLatex(0.15,0.56,Form("a_{2} = %6.5f  MeV^{-1}",a2));
	//lat->DrawLatex(0.15,0.51,Form("#chi^{2}/ndf = %4.2f/%d",chi2,ndf));
	//lat->DrawLatex(0.15,0.47,Form("Prob(#chi^{2}) = %6.4f",prob));


	TLegend *l = new TLegend(0.6, 0.6, 0.9, 0.75);
	l->SetTextSize(0.03);
	l->AddEntry(h, "Data", "lp");
	l->AddEntry(f1, "Convolution model", "lp");
	l->Draw();

	c1->Print("ConvNonlinear.pdf");

 } // ifdoPlot

	// nonlinearity function
	TF1 *nonlin = new TF1("nonlin","[0]+[1]*x/(1+[2]*x)",0,Emax); 
        double a0 = f1->GetParameter(3);
        double a1 = f1->GetParameter(4);
        double a2 = f1->GetParameter(5);
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
