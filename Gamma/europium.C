#include "gethisto.C"

/// \author Javier Mariño Villadamigo

Double_t cal_NaI(Int_t ch) {return (-9.32 + ch*0.12081);}
Double_t cal_HPGe(Int_t ch) {return (0.27 + ch*0.17434);}



void europium() {

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
    Int_t topCh = int(TMath::Power(2, 11));
	Int_t bins = topCh;
    Int_t xMin, xMax;
    short chan;
    Float_t new_end;

    const int npeaks = 11;
    Float_t mean_vec[npeaks] = {};
    Float_t sigma_vec[npeaks] = {};
    Float_t area[npeaks] = {};



    xMin = (10 - 0.27) / 0.17434, xMax = 10000;
    chan = 1;
    

    TCanvas *c1 = new TCanvas("c1", "Eu calibrated spectrum");
    c1->SetLogy();
    TH1D *hist = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/HPGe_Eu20min.root", chan, bins, xMin, xMax);
    TH1D *bg = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/background29min.root", chan, bins, xMin, xMax);
    new_end = cal_HPGe(xMax);


    TH1D *hcal = new TH1D("hist_cal", "Europium spectrum", bins, cal_HPGe(xMin), cal_HPGe(xMax)); //create calibrated histo!!
    for (Int_t i = 0; i < bins; i++) hcal->SetBinContent(i+1, hist->GetBinContent(i+1)); //fill cal histo with same content!
    TH1D *hb = new TH1D("bg_cal", "Calibrated background", bins, cal_HPGe(xMin), cal_HPGe(xMax)); //create calibrated histo!!
    for (Int_t i = 0; i < bins; i++) hb->SetBinContent(i+1, bg->GetBinContent(i+1)); //fill cal histo with same content!

    hcal->Rebin(2);
    hcal->SetTitle("");
    hcal->GetXaxis()->SetTitle("Photon energy [keV] "); hcal->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal->GetXaxis()->GetBinWidth(0)));
    hcal->GetXaxis()->SetRangeUser(xMin, xMax);
    hcal->GetYaxis()->SetTitleOffset(1.25);
    hcal->SetLineColor(kBlue);
    hcal->SetFillColorAlpha(kAzure+7, 0.4);
    hcal->Draw();

    hb->Rebin(2);
    hb->SetLineColor(kRed);
    hb->SetFillColorAlpha(kRed, 0.4);
    hb->Draw("SAME");

    
    TCanvas *c2 = new TCanvas("c2", "Subtracted background spectrum");
    c2->SetLogy();
    TH1D *hclone = (TH1D*)hcal->Clone();
    hclone->Add(hclone, hb, 1, -1); //substracting background from initial spectrum
    hclone->Draw();
    for (Int_t i = 0; i<hclone->GetEntries(); i++) if (hclone->GetBinContent(i) < 0) hclone->SetBinContent(i, 0);

    Float_t xpeaks[npeaks] = {121.8, 344.3, 444.0, 778.9, 867.4, 964.0, 1089.7, 1112.1, 1212.9, 1299.1, 1408.0}; // 1089.7 has the contribution of also 1085.8
    Float_t intens[npeaks] = {141.0, 127.2, 15.00, 62.6, 20.54, 70.4, (48.7 + 8.26), 65.0, 6.67, 7.76, 100.};

    TF1 *fun1, *fun2, *clone_fun;
    for(Int_t i = 0; i<npeaks; i++)
    {
        Float_t xp = xpeaks[i];


        fun1 = new TF1("fun1name", "gaus", xp - 5, xp + 5);
        
        hclone->Fit(fun1,"RQN");


        fun2 = new TF1("fun2name","gaus(0) + exp(3)", fun1->GetParameter(1) - 2.5 * fun1->GetParameter(2), fun1->GetParameter(1) + 2.5 * fun1->GetParameter(2));
        
        fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
        fun2->SetParLimits(1, xpeaks[i] - 8, xpeaks[i] + 8);
        fun2->SetParLimits(2, 0, 20);
        fun2->SetLineColor(kRed);

        hclone->Fit(fun2,"RQN+");


        mean_vec[i] = fun2->GetParameter(1);
        sigma_vec[i] = fun2->GetParameter(2);

        clone_fun = (TF1*)fun2->Clone();
        clone_fun->SetRange(mean_vec[i] - 3 * sigma_vec[i], mean_vec[i] + 3 * sigma_vec[i]);

        clone_fun->SetLineColor(kRed);
        clone_fun->SetFillColorAlpha(kRed, 0.3);
        clone_fun->SetFillStyle(1001);

        clone_fun->Draw("SAME");

        area[i] = fun2->Integral(mean_vec[i] - 3 * sigma_vec[i], mean_vec[i] + 3 * sigma_vec[i]) / hclone->GetXaxis()->GetBinWidth(0);
        cout << "For" << xp << " " << area[i] * 100 << endl;
    }

    Float_t rel_efficiency[npeaks] = {}, urel_efficiency[npeaks] = {};

    for (Int_t i = 0; i<npeaks; i++)
    {
        rel_efficiency[i] = area[i] / area[npeaks-1] * 100 / intens[i];
        urel_efficiency[i] = 100 * 1 / rel_efficiency[i] * TMath::Sqrt(1 / area[i] + 1 / area[npeaks-1]);
    }


    TCanvas *c3 = new TCanvas("c3", "Relative efficiency");
    c3->SetLeftMargin(0.11749);
    c3->SetRightMargin(0.0544413);
    c3->SetBottomMargin(0.101053);
    c3->SetTopMargin(0.0336842);
    TH2F *h = new TH2F("h", ";Photon energy [keV];Relative efficiency [%]",  1000, 0, 1500, 1000, 0, 11);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->Draw();
    TGraphErrors *rel_ef_graph = new TGraphErrors(npeaks, mean_vec, rel_efficiency, sigma_vec, urel_efficiency); //calibration: channel vs energy of the found peaks
    rel_ef_graph->SetMarkerStyle(20);
    rel_ef_graph->SetMarkerColor(1);
    rel_ef_graph->Draw("P");
    
    TF1 *expo = new TF1("expo", "expo", 0, 1500);
    expo->SetLineColor(kRed);
    rel_ef_graph->Fit(expo, "RWW+");


}