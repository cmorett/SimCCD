void Test_edep_Chain(){

    TChain *chain = new TChain("tree");
    // chain->Add("/home/bruce/Documents/Programas/Simulacion_ab_initio/treesROOT_CCD/100k/Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_0.root");
    // chain->Add("Sim_ab_initio_NMUONS_50000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_.root");
    chain->Add("Sim_ab_initio_NMUONS_100000_PLANES_1x1_RADIO_5_CCDSIZE_400x600_SIGMA_LV_1_.root");
    chain->Add("Edep_allclusters_NSAMP324_KeV.root");


    // chain->Add("treesROOT_CCD/100k/Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_1.root");
    // chain->Add("treesROOT_CCD/100k/Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_2.root");
    // chain->Add("treesROOT_CCD/100k/Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_3.root");
    // chain->Add("treesROOT_CCD/100k/Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_4.root");
    // chain->Add("treesROOT_CCD/100k/Sim_ab_initio_NMUONS_100000_PLANES_3.0x3.0_RADIO_12_5.root");
    // chain->Draw("thet");


    // Sección de Energía depositada //
    int NB = 150;
    double tlow = 0;
    double thi = 1000;

    TH1F *edep = new TH1F("edep", "", NB, tlow, thi); // Histogram simulation all
    edep->GetXaxis()->SetTitle("Energy (KeV)");
    // edep->Scale(5);

    TH1F *edep_cut = new TH1F("edep_cut", "Sigma: 1", NB, tlow, thi); // Histogram simulation cut 
    edep_cut->GetXaxis()->SetTitle("Energy (KeV)");
    edep_cut->SetLineColor(1);

    TH1F *edep_exp = new TH1F("edep_exp", "", NB, tlow, thi); // Histogram simulation cut 
    edep_exp->GetXaxis()->SetTitle("Energy (KeV)");
    // edep_exp->Scale(5);

    chain->Draw("edep>>edep", "l>0");
    chain->Draw("edep>>edep_cut", "thet>22*TMath::Pi()/180 & edep>0");
    chain->Draw("edep>>edep_exp", "l<0");


    //scale hint1 to the pad coordinates
    double rightmax = 0.00168*edep->GetMaximum();
    double scale = gPad->GetUymax()/rightmax;
    // hint1->SetLineColor(kRed);
    edep->Scale(scale);
    edep->SetLineColor(2);
    // edep->SetLineStyle(0);
    // edep->SetFillStyle(0);
    // hint1->Draw("same");


    rightmax = 0.00168*edep_cut->GetMaximum();
    scale = gPad->GetUymax()/rightmax;
    // hint1->SetLineColor(kRed);
    edep_cut->Scale(scale);
    // hint1->Draw("same");

    // TGaxis*axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
    //                         gPad->GetUxmax(),gPad->GetUymax(),
    //                         0,rightmax,510,"+L");

    TCanvas *canv = new TCanvas("canv","", 2*700, 600);
    canv->Divide(2,1);
    canv->cd(1);
    edep->Draw();
    edep_exp->Draw("same");

    // axis->SetLineColor(2);
    // axis->SetLabelColor(2);
    // axis->Draw();
    // func1->Draw("same");

    canv->cd(2);
    edep_cut->Draw();
    edep_exp->Draw("same");

    // // Sección de theta //
    // // int NB = 90;
    // int NB = 300;
    // double tlow = 0;
    // double thi = TMath::Pi()/2.0;
    // TH1F *theta_all = new TH1F("theta_all", "", NB, tlow, thi);
    // TH1F *theta_in = new TH1F("theta_in", "", NB, tlow, thi);

    // chain->Draw("thet>>theta_all");
    // chain->Draw("thet>>theta_in", "l>0");

    // // Define fuctions //
    // TF1 *func1 = new TF1("func1", "[0]*sin(x)*(cos(x))^2", 0.01, 85*TMath::Pi()/180);
    // func1->SetParameter(0, 5000);

    // // ====================== Función para la CCD =========================== //
    // // TF1 *func2 = new TF1("func2", "[0]*((3.29909/1)*sin(x)*(cos(x))^3+(2.8041616/(1*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0.01, 85*TMath::Pi()/180);
    // TF1 *func2 = new TF1("func2", "[0]*((3.4275732/1)*sin(x)*(cos(x))^3+(2.289706464/(1*TMath::Pi()))*(sin(x))^2*(cos(x))^2)", 0.01, 85*TMath::Pi()/180);
    
    // theta_all->Fit("func1", "R");
    // theta_in->Fit("func2", "R");

    // TCanvas *canv = new TCanvas("canv","", 2*700, 600);
    // canv->Divide(2,1);
    // canv->cd(1);
    // theta_all->Draw();
    // func1->Draw("same");

    // canv->cd(2);
    // theta_in->Draw();
    // func2->Draw("same");

    // double Prob1 = func1->GetProb();
    // double chi1 = func1->GetChisquare();
    // // double NFD1 = func1->GetN

    // double Prob2 = func2->GetProb();
    // double chi2 = func2->GetChisquare();

    // std::cout << "ChiSquare1 = "<< chi1 << std::endl;
    // std::cout << "Prob1 = "<< Prob1 << std::endl;

    // std::cout << "ChiSquare2 = "<< chi2 << std::endl;
    // std::cout << "Prob2 = "<< Prob2 << std::endl;

    canv->Print("Dis_edep_chain_SIGMA_1.pdf");

}
