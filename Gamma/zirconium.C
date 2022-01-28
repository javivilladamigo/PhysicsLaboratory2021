#include "gethisto.C"

/// \author Javier Mariño Villadamigo

Double_t cal_NaI(Int_t ch) {return (-9.32 + ch*0.12081);}
Double_t cal_HPGe(Int_t ch) {return (0.27 + ch*0.17434);}
Float_t eff_HPGe(Float_t E) {return (-128.4 + 8227 / TMath::Sqrt(E) + 10180 / E);}
Float_t eff_NaI(Float_t E) {return (1.069 + 56.62 / TMath::Sqrt(E));}


void zirconium() {

    //Style settings
	gROOT->SetStyle("Default");
	gStyle->SetCanvasColor(0);
    gStyle->SetStatColor(0);
	gStyle->SetStatBorderSize(1);
    gStyle->SetLabelFont(62, "XY");
    gStyle->SetLegendFont(62);
    gStyle->SetTitleFont(62, "XY");
    gStyle->SetOptFit(111);
	gStyle->SetPalette(0);
    gStyle->SetOptStat(0);

    // Variables
    Int_t topCh = int(TMath::Power(2, 14));
	Int_t bins = topCh;
    Int_t xMin, xMax;
    Float_t xMin_en, xMax_en;
    Float_t new_end;

    Float_t time = 15;
    Float_t normalization_factor = time / 29.; // 29 minutes for background and 15 minutes for zirconium sample
    xMin = 0, xMax = 30000;
    xMin_en = 40, xMax_en = 2000;
    
    
    


    TH1D *hist_NaI = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/samples/zirconium10min_1870g.root", 0, bins, xMin, xMax);
    TH1D *hist_HPGe = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/samples/zirconium10min_1870g.root", 1, bins, xMin, xMax);
    TH1D *hb_NaI = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/background29min.root", 0, bins, xMin, xMax);
    TH1D *hb_HPGe = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/background29min.root", 1, bins, xMin, xMax);
    
    Float_t new_beg_NaI = cal_NaI(xMin);
    Float_t new_beg_HPGe = cal_HPGe(xMin);
    Float_t new_end_NaI = cal_NaI(xMax);
    Float_t new_end_HPGe = cal_HPGe(xMax);

    TH1D *hcal_NaI = new TH1D("hcal_NaI", "Calibrated histo of zirconium (NaI)", bins, new_beg_NaI, new_end_NaI);
    TH1D *hcal_HPGe = new TH1D("hcal_HPGe", "Calibrated histo of zirconium (HPGe)", bins, new_beg_HPGe, new_end_HPGe);
    TH1D *hbcal_NaI = new TH1D("hbcal_NaI", "Calibrated histo of zirconium (NaI)", bins, new_beg_NaI, new_end_NaI);
    TH1D *hbcal_HPGe = new TH1D("hbcal_HPGe", "Calibrated histo of zirconium (HPGe)", bins, new_beg_HPGe, new_end_HPGe);

    for (Int_t i = 0; i < bins; i++)
    {
        hcal_NaI->SetBinContent(i, hist_NaI->GetBinContent(i));
        hcal_HPGe->SetBinContent(i+1, hist_HPGe->GetBinContent(i+1));
        hbcal_NaI->SetBinContent(i, hb_NaI->GetBinContent(i) * normalization_factor);
        hbcal_HPGe->SetBinContent(i, hb_HPGe->GetBinContent(i) * normalization_factor);
    }
    

    TCanvas *c1 = new TCanvas("c1", "zirconium calibrated spectrum (NaI)");
    c1->SetWindowSize(1102, 771);
    c1->SetLeftMargin(0.11749);
    c1->SetRightMargin(0.0344413);
    c1->SetBottomMargin(0.101053);
    c1->SetTopMargin(0.0336842);
    hcal_NaI->Rebin(8);
    hcal_NaI->SetTitle("");
    hcal_NaI->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_NaI->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_NaI->GetXaxis()->GetBinWidth(0)));
    hcal_NaI->GetXaxis()->SetRangeUser(xMin_en, xMax_en);
    hcal_NaI->GetYaxis()->SetTitleOffset(1.);
    hcal_NaI->SetLineColor(kBlue);
    hcal_NaI->SetFillColorAlpha(kAzure+7, 0.4);
    hcal_NaI->Draw();

    hbcal_NaI->Rebin(8);

    hbcal_NaI->GetXaxis()->SetRangeUser(xMin_en, xMax_en);
    hbcal_NaI->SetLineColor(kRed);
    hbcal_NaI->SetFillColorAlpha(kRed, 0.4);
    hbcal_NaI->Draw("SAME");

    TCanvas *c2 = new TCanvas("c2", "Sub. background NaI");
    c2->SetWindowSize(1102, 771);
    c2->SetLeftMargin(0.11749);
    c2->SetRightMargin(0.0344413);
    c2->SetBottomMargin(0.101053);
    c2->SetTopMargin(0.0336842);
    TH1D *hNaI = (TH1D*)hcal_NaI->Clone();
    hNaI->Add(hcal_NaI, hbcal_NaI, 1, -1); //substracting background from initial spectrum
    for (Int_t i = 0; i < hNaI->GetNbinsX(); i++) if (hNaI->GetBinContent(i) < 0) hNaI->SetBinContent(i,0);
    hNaI->Draw();


    TCanvas *c3 = new TCanvas("c3", "zirconium calibrated spectrum (HPGe)");
    c3->SetWindowSize(1102, 771);
    c3->SetLeftMargin(0.11749);
    c3->SetRightMargin(0.0344413);
    c3->SetBottomMargin(0.101053);
    c3->SetTopMargin(0.0336842);
    hcal_HPGe->Rebin(8);
    hcal_HPGe->SetTitle("");
    hcal_HPGe->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_HPGe->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_HPGe->GetXaxis()->GetBinWidth(0)));
    hcal_HPGe->GetXaxis()->SetRangeUser(xMin_en, xMax_en);
    hcal_HPGe->GetYaxis()->SetTitleOffset(1.);
    hcal_HPGe->SetLineColor(kBlue);
    hcal_HPGe->SetFillColorAlpha(kAzure+7, 0.4);
    hcal_HPGe->Draw();

    hbcal_HPGe->Rebin(8);
    hbcal_HPGe->GetXaxis()->SetRangeUser(xMin_en, xMax_en);
    hbcal_HPGe->SetLineColor(kRed);
    hbcal_HPGe->SetFillColorAlpha(kRed, 0.4);
    hbcal_HPGe->Draw("SAME");


    TCanvas *c4 = new TCanvas("c4", "Sub. background HPGe");
    c4->SetWindowSize(1102, 771);
    c4->SetLeftMargin(0.11749);
    c4->SetRightMargin(0.0344413);
    c4->SetBottomMargin(0.101053);
    c4->SetTopMargin(0.0336842);
    TH1D *hHPGe = (TH1D*)hcal_HPGe->Clone();
    hHPGe->Add(hcal_HPGe, hbcal_HPGe, 1, -1); //substracting background from initial spectrum
    for (Int_t i = 0; i < hHPGe->GetNbinsX(); i++)
    {
        if (hHPGe->GetBinContent(i) < 0)
        {
            
            hHPGe->SetBinContent(i, 0);
        }
    }
    hHPGe->Draw();


    // here begins the calculations

    Float_t N_NaI = hNaI->Integral(hNaI->FindFixBin(mean_NaI-3*sigma_NaI), hNaI->FindFixBin(mean_NaI+3*sigma_NaI));
    //use this line to integrate the gaussian of a peak you want in (-3sigma +3 sigma)






    /*
    // 1460 keV transition 
    Float_t P_gamma = 0.1067; // Branching ratio of 1460 keV transition

    Float_t sigma_NaI = 27.24;
    Float_t mean_NaI = 1451;

    Float_t N_NaI = hNaI->Integral(hNaI->FindFixBin(mean_NaI-3*sigma_NaI), hNaI->FindFixBin(mean_NaI+3*sigma_NaI));


    Float_t sigma_HPGe = 2.435;
    Float_t mean_HPGe = 1464;
    Float_t N_HPGe = hHPGe->Integral(hHPGe->FindFixBin(mean_NaI-3*sigma_NaI), hHPGe->FindFixBin(mean_NaI+3*sigma_NaI));
    cout << "integral NaI: " << N_NaI << endl;
    cout << "integral HPGe: " << N_HPGe << endl;

    // Geometry data
    Float_t distance_source_NaI = 2.5; // cm
    Float_t NaI_radius = 7.5/2.; // cm

    Float_t distance_source_HPGe = 7.5; // cm
    Float_t HPGe_radius = TMath::Sqrt(1200 / TMath::Pi()) / 10; // cm

    Float_t acceptancy_NaI = TMath::Pi() * TMath::Power(NaI_radius, 2) / TMath::Power(distance_source_NaI, 2);
    Float_t acceptancy_HPGe = TMath::Pi() * TMath::Power(HPGe_radius, 2) / TMath::Power(distance_source_HPGe, 2);


    cout << "eff: " << eff_NaI(mean_NaI) << endl;

    Float_t A_NaI = 4 * TMath::Pi() / acceptancy_NaI / P_gamma / (time * 60) * N_NaI / (eff_NaI(mean_NaI) / 100) / 8; // we have divided the efficiency by 100 to obtain a (0-1) number
    cout << "Mass activity for 40K in NaI detector: " << A_NaI << " Bq/g" << endl;
    Float_t A_HPGe = 4 * TMath::Pi() / acceptancy_HPGe / P_gamma / (time * 60) * N_HPGe / (eff_NaI(mean_HPGe) / 100) / 8; // we have divided the efficiency by 100 to obtain a (0-1) number
    cout << "Mass activity for 40K in HPGe detector: " << A_HPGe << " Bq/g" << endl;
    */
}