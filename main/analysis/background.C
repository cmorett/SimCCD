void background() 
{ //begin

//-----------------------------------
// background shape
//-----------------------------------

double Nbg1 = 1.0; // (MeV)
double eps1 = 20; // (MeV)

TF1 *fbg1 = new TF1("fbg1", "[0]*(1/[1])*exp(-x/[1])", 0,2000);
fbg1->SetLineColor(2);
fbg1->SetParameter(0, Nbg1);
fbg1->SetParameter(1, eps1);

//-----------------------------------
// Background histogram
//-----------------------------------

int NBb1 = 200;
double maxEd = 700.;
double backgroundbw = maxEd/NBb1;
TH1F *background = new TH1F("background","",NBb1,0,maxEd);
background->GetXaxis()->SetTitle("Energy (MeV)");
background->GetYaxis()->SetTitle("events");

//-----------------------------------
// Fill histogram
//-----------------------------------

int Nevbg1 = 1000000;

for (int i=0; i<Nevbg1; i++){
     double e = fbg1->GetRandom();  
     background->Fill(e);
}

//----------------------------
//-- background template data file 
//----------------------------

  ofstream myfile;
  myfile.open ("background.dat");
  for (int i=0;i<NBb1; i++) 
  myfile << background->GetBinContent(i+1)/(background->Integral(1,NBb1)*backgroundbw) << "\n";
  myfile.close();

  cout << "Integral background: " << background->Integral() << " entries." << endl;

  TCanvas *canvmu = new TCanvas("canvmu","",900,700);
  canvmu->cd(1);
  gPad->SetLogy(1);
  background->Draw();
  canvmu->Print("background.pdf");


//-----------------------------------
// draw
//-----------------------------------

TCanvas *canv1 = new TCanvas("canv1","",900,700);
canv1->cd();
background->Draw();
fbg1->SetParameter(0,Nevbg1*backgroundbw);
fbg1->Draw("same");
canv1->Print("background.pdf");


}// ebd
