#include "gethisto.C"

/// \author Javier Mariño Villadamigo

Double_t cal_NaI(Int_t ch) {return (-9.32 + ch*0.12081);}
Double_t cal_HPGe(Int_t ch) {return (0.27 + ch*0.17434);}



void cement() {

    //Style settings
	//gROOT->SetStyle("Default");
	//gStyle->SetCanvasColor(0);
	//gStyle->SetStatBorderSize(1);
	gStyle->SetOptStat("");
    gStyle->SetLabelFont(62, "XY");
    gStyle->SetLegendFont(62);
    gStyle->SetTitleFont(62, "XY");
    gStyle->SetOptFit(111);
	//gStyle->SetPalette(0);

    // Variables
    Int_t topCh = int(TMath::Power(2, 10));
	Int_t bins = topCh;
    Int_t xMin, xMax;
    Float_t xMin_en, xMax_en;
    Float_t new_end;

    Float_t normalization_factor = 15. / 29.; // 29 minutes for background and 15 minutes for cement sample
    xMin = 0, xMax = 30000;
    xMin_en = 20, xMax_en = 3000;
    
    
    


    TH1D *hist_NaI = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/samples/cement15min_521g.root", 0, bins, xMin, xMax);
    TH1D *hist_HPGe = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/samples/cement15min_521g.root", 1, bins, xMin, xMax);
    TH1D *hb_NaI = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/background29min.root", 0, bins, xMin, xMax);
    TH1D *hb_HPGe = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/background29min.root", 1, bins, xMin, xMax);
    
    Float_t new_beg_NaI = cal_NaI(xMin);
    Float_t new_beg_HPGe = cal_HPGe(xMax);
    Float_t new_end_NaI = cal_NaI(xMax);
    Float_t new_end_HPGe = cal_HPGe(xMax);

    TH1D *hcal_NaI = new TH1D("hcal_NaI", "Calibrated histo of cement (NaI)", bins, new_beg_NaI, new_end_NaI);
    TH1D *hcal_HPGe = new TH1D("hcal_HPGe", "Calibrated histo of cement (HPGe)", bins, new_beg_HPGe, new_end_HPGe);
    TH1D *hbcal_NaI = new TH1D("hbcal_NaI", "Calibrated histo of cement (NaI)", bins, new_beg_NaI, new_end_NaI);
    TH1D *hbcal_HPGe = new TH1D("hbcal_HPGe", "Calibrated histo of cement (HPGe)", bins, new_beg_HPGe, new_end_HPGe);

    for (Int_t i = 0; i < bins; i++)
    {
        hcal_NaI->SetBinContent(i, hist_NaI->GetBinContent(i));
        hcal_HPGe->SetBinContent(i+1, hist_HPGe->GetBinContent(i+1));
        hbcal_NaI->SetBinContent(i, hb_NaI->GetBinContent(i) * normalization_factor);
        hbcal_HPGe->SetBinContent(i, hb_HPGe->GetBinContent(i) * normalization_factor);
    }
    

    TCanvas *c1 = new TCanvas("c1", "Cement calibrated spectrum (NaI)");
    c1->cd();
    hcal_NaI->Rebin(4);
    hcal_NaI->SetTitle("");
    hcal_NaI->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_NaI->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_NaI->GetXaxis()->GetBinWidth(0)));
    hcal_NaI->GetXaxis()->SetLimits(20, xMax_en);
    hcal_NaI->GetXaxis()->SetRangeUser(20, xMax_en);
    hcal_NaI->GetYaxis()->SetTitleOffset(1.25);
    hcal_NaI->SetLineColor(kBlue);
    hcal_NaI->SetFillColorAlpha(kAzure+7, 0.4);
    hcal_NaI->Draw();

    hbcal_NaI->Rebin(4);
    hbcal_NaI->GetXaxis()->SetLimits(20, xMax_en);
    hbcal_NaI->GetXaxis()->SetRangeUser(20, xMax_en);
    hbcal_NaI->SetLineColor(kRed);
    hbcal_NaI->SetFillColorAlpha(kRed, 0.4);
    hbcal_NaI->Draw("SAME");


    TCanvas *c2 = new TCanvas("c2", "Cement calibrated spectrum (HPGe)");
    c2->cd();
    hcal_HPGe->Rebin(4);
    hcal_HPGe->SetTitle("");
    hcal_HPGe->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_NaI->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_NaI->GetXaxis()->GetBinWidth(0)));
    hcal_HPGe->GetXaxis()->SetLimits(0, xMax_en);
    hcal_HPGe->GetXaxis()->SetRangeUser(20, xMax_en);
    hcal_HPGe->GetYaxis()->SetTitleOffset(1.25);
    hcal_HPGe->SetLineColor(kBlue);
    hcal_HPGe->SetFillColorAlpha(kAzure+7, 0.4);
    hcal_HPGe->Draw();

    hbcal_HPGe->Rebin(4);
    hbcal_HPGe->GetXaxis()->SetLimits(0, xMax_en);
    hbcal_HPGe->GetXaxis()->SetRangeUser(20, xMax_en);
    hbcal_HPGe->SetLineColor(kRed);
    hbcal_HPGe->SetFillColorAlpha(kRed, 0.4);
    hbcal_HPGe->Draw("SAME");


    TCanvas *c3 = new TCanvas("c3", "Sub. background NaI");
    c3->cd();
    TH1D *hNaI = (TH1D*)hcal_NaI->Clone();
    hNaI->Add(hcal_NaI, hbcal_NaI, 1, -1); //substracting background from initial spectrum
    for (Int_t i = 0; i < hNaI->GetEntries(); i++) if (hNaI->GetBinContent(i) < 0) hNaI->SetBinContent(i,0);
    hNaI->Draw();

}
