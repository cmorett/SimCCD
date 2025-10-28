void twoHistograms()
{

TFile *f1 = new TFile("muons.root");
TFile *f2 = new TFile("muonswithoutshielding.root");


TCanvas *canvmuons = new TCanvas("canvmuons","",800,600);

canvmuons->cd(1);
gPad->SetLogy(1);
  
TH1F *h1 = new TH1F("h1", "muons", 200, 0, 700);
h1 = (TH1F*)f1->Get("muons");

TH1F *h2 = new TH1F("h2", "muonswithoutshielding", 200, 0, 700);
h2 = (TH1F*)f2->Get("muonswithoutshielding");
h2->SetLineColor(3);

h1->Draw();
h2->Draw("same");

 auto legend = new TLegend(0.2, 0.2, 0.7 , 0.3);
 legend->SetTextSize(0.03);
   //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(h1,"Deposited energy with shielding","l");
   legend->AddEntry(h2,"Deposited energy without shielding","l"); 
   legend->Draw();
canvmuons->SetGrid();
canvmuons->Print("twoHistograms.pdf"); 

};
