double Parametros_muons(){

    // ==== Datos de la CCD y simulación === //
    // double Nv = 250;    // px
    double Nv = 529;    // px
    double Nh = 250;    // px
    double px_size = 1.5 * pow(10, -5); // m
    double Nsamp = 324; 
    double tau = 50 * pow(10,-6); // seg
    double I0 = 101.2;


    double Nmuons_sim = 4000000; 
    double planes_size = 0.015; // m
    double Nmuons_t = (2 * TMath::Pi()/3) * I0 * pow(planes_size,2);
    // ==================================== //

    // ==== Número de extensiones === //
    // double Next1 = 388;
    // double Next2 = 378;

    double Next1 = 385;
    double Next2 = 386;
    double Next_total = Next1 + Next2;

    // double Next1 = 684;
    // double Next2 = 0;
    // double Next_total = Next1 + Next2;
    // ============================= //

    // ==== Tiempos === //
    double Tmax = Nh * Nv * Nsamp * tau; //seg
    double T_mean = 0.5 * Tmax;
    double T_exp = Next_total * T_mean;
    double T_sim = Nmuons_sim / Nmuons_t;
    // =============== //

    double dN_dtdA = (2 * TMath::Pi()/4) * I0;
    double Ah = Nh * Nv * pow(px_size, 2);
    double dN_dt = dN_dtdA * Ah; 
    double dNmuons = dN_dt * T_mean;

    double Totmuons = dNmuons * Next_total;

    cout << "Nh: " << Nh << " px" << endl;
    cout << "Nv: " << Nv << " px" << endl;
    cout << "Muons/time: " << Nmuons_t << " muons/seg " << endl;
    cout << "Ntotal_img: " << Next_total << endl;

    cout << "T_max: " << Tmax << " seg" << endl;
    cout << "T_mean: " << T_mean << " seg" << endl;
    cout << "T_exp: " << T_exp << " seg" << endl;
    cout << "T_sim: " << T_sim << " seg" <<endl;

    cout << "dN/dtdA: " << dN_dtdA << " m^-2 s^-1"<< endl;
    cout << "Ah: " << Ah << " m^2" << endl;
    cout << "dN/dt: " << dN_dt << " seg^-1" << endl;
    cout << "dNmuons: " << dNmuons <<  " muons" <<endl;
    cout << "Nmuons: " << Totmuons << " muons " << endl;

    return 0;

}