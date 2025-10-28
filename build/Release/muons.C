void muons()
{//begin


TFile *fin = new TFile("muons_output.root");

TTree *B02Evts = (TTree*) fin->Get("B02Evts");

TCanvas *canv = new TCanvas("canv","title",600,350);
canv->cd();
B02Evts->Draw("thetaPri");

canv->Print("muons.pdf");

}//end