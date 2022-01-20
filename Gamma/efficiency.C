#include "../../examples/data_analysis_examples/gethisto.C"

/// \author Javier Mariño Villadamigo

Double_t cal_NaI(Int_t ch) {return (-9.32 + ch*0.12081);}
Double_t cal_HPGe(Int_t ch) {return (0.27 + ch*0.17434);}


void efficiency(TString detector)
{   
    //Style settings
	gROOT->SetStyle("Default");
	gStyle->SetCanvasColor(0);
	gStyle->SetStatBorderSize(1);
	gStyle->SetOptFit(0000);
	gStyle->SetOptStat("");
	gStyle->SetPalette(0);

    // Activity data
    Float_t BR_Co1 = 0.9988;
    Float_t BR_Co2 = 0.9988 + 0.0012;
    Float_t BR_Am = 0.852 * 0.360;
    
    Float_t halflife_Co = 1925.28; // days
    Float_t halflife_Am = 158007.15; // days
    Float_t activity_Co_peak1 = 195000 * TMath::Exp( - TMath::Log(2) / halflife_Co * 2020) * BR_Co1; // a 1/07/2016: A = 195 kBq ---> 2020 days have passed
    Float_t activity_Co_peak2 = 195000 * TMath::Exp( - TMath::Log(2) / halflife_Co * 2020) * BR_Co2;
    Float_t activity_Am = 419000 * TMath::Exp( - TMath::Log(2) / halflife_Am * 163) * BR_Am;; // a 1/08/2021; A = 419 kBq ---> 163 days have passed

    Float_t emitted_Co_peak1 = activity_Co_peak1 * 10*60; // ten minutes
    Float_t emitted_Co_peak2 = activity_Co_peak2 * 10*60; // ten minutes
    Float_t emitted_Am = activity_Am * 10*60; // ten minutes

    cout << "Activity of Co at peak 1: " << activity_Co_peak1 << " and at peak 2: " << activity_Co_peak2 << endl;
    cout << "Activity of Am: " << activity_Am << endl;


    // Geometry data
    Float_t distance_source_NaI = 28.5; // cm
    Float_t NaI_radius = 7.5/2.; // cm

    Float_t distance_source_HPGe = 20; // cm
    Float_t HPGe_radius = TMath::Sqrt(1200 / TMath::Pi()) / 10; // cm

    Float_t acceptancy_NaI = TMath::Pi() * TMath::Power(NaI_radius, 2) / TMath::Power(distance_source_NaI, 2);
    Float_t acceptancy_HPGe = TMath::Pi() * TMath::Power(HPGe_radius, 2) / TMath::Power(distance_source_HPGe, 2);
    // In the limit of point-like source and long distance from the detector --> solid angle = A_detector / d^2


    Int_t topCh = int(TMath::Power(2, 12));
	Int_t bins = topCh;
    Int_t xMin, xMax;
    Float_t yMin = 0, yMax = 15;



    
    Float_t received_Co_peak1;
    Float_t received_Co_peak2;
    Float_t received_Am;

    Float_t efficiency_Co_peak1, uefficiency_Co_peak1;
    Float_t efficiency_Co_peak2, uefficiency_Co_peak2;
    Float_t efficiency_Am, uefficiency_Am;

    if (detector == "NaI")
    {   

        Float_t received_Co_peak1 = emitted_Co_peak1 * acceptancy_NaI;
        Float_t received_Co_peak2 = emitted_Co_peak2 * acceptancy_NaI;
        Float_t received_Am = emitted_Am * acceptancy_NaI;

        xMin = 0, xMax = 12000;
        short chan = 0;

        TCanvas *c1 = new TCanvas("c1", "Co calibrated spectrum");
        TH1D *hist_Co = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day1/10min_Co.root", chan, bins, xMin, xMax);

        Float_t new_end = cal_NaI(xMax);
        TH1D *hcal_Co = new TH1D("hist_cal", "NaI (Co source) spectrum", bins, cal_NaI(xMin), new_end); //create calibrated histo!!
        for(Int_t i = 0; i < bins; i++) hcal_Co->SetBinContent(i+1, hist_Co->GetBinContent(i+1)); //fill cal histo with same content!
        
        hcal_Co->Rebin(4);
        hcal_Co->SetTitle("");
        hcal_Co->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_Co->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_Co->GetXaxis()->GetBinWidth(0)));
        hcal_Co->GetXaxis()->SetRangeUser(xMin, xMax);
        hcal_Co->GetYaxis()->SetTitleOffset(1.25);
        hcal_Co->SetLineColor(kBlue);
        hcal_Co->SetFillColorAlpha(kAzure+7, 0.4);
        hcal_Co->Draw();


        TSpectrum *ss = new TSpectrum(2); //creating TSpectrum object
        TH1D *hclone_Co = (TH1D*)hcal_Co->Clone();
        TH1 *hb = ss->Background(hclone_Co,20,"same"); //calculating backgraund (TSpecrum class method)

        TCanvas *c2 = new TCanvas("c2", "Subtracted background spectrum");
        hclone_Co->Add(hclone_Co,hb,1,-1); //substracting background from initial spectrum
        TH1D *hfinal_Co = (TH1D*)hclone_Co->Clone();
        for (Int_t i = 0; i<hfinal_Co->GetEntries(); i++)
        {
            if (hfinal_Co->GetBinContent(i) < 0) hfinal_Co->SetBinContent(i, 0);
        }
        hfinal_Co->Draw();
        
        Float_t xpeaks[2] = {1173.228, 1332.514};

        Float_t mean_vec[2] = {};
        Float_t sigma_vec[2] = {};
        Float_t area[2] = {};


        TF1 *fun1, *fun2, *clone_fun;
        for(Int_t i=0; i<2; i++)
        {
            Float_t xp = xpeaks[i];
            //cout << " ************ Peak at: " << xp << endl;

            fun1 = new TF1("fun1name", "gaus", xp-120,xp+120);

            hfinal_Co->Fit(fun1,"RQN");

            fun2 = new TF1("fun2name","gaus", fun1->GetParameter(1)-(3.*fun1->GetParameter(2)), fun1->GetParameter(1)+(3.0*fun1->GetParameter(2)));
            fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
            fun2->SetLineColor(kRed);

            hfinal_Co->Fit(fun2,"RQN+");
            



            mean_vec[i] = fun2->GetParameter(1);
            sigma_vec[i] = fun2->GetParameter(2);

            clone_fun = (TF1*)fun2->Clone();
            clone_fun->SetRange(mean_vec[i]-3*sigma_vec[i], mean_vec[i]+3*sigma_vec[i]);

            clone_fun->SetLineColor(kRed);
            clone_fun->SetFillColorAlpha(kRed, 0.3);
            clone_fun->SetFillStyle(1001);

            clone_fun->Draw("SAME");

            area[i] = fun2->Integral(mean_vec[i] - 3 * sigma_vec[i], mean_vec[i] + 3 * sigma_vec[i]) / hfinal_Co->GetXaxis()->GetBinWidth(0);

        }

        efficiency_Co_peak1 = area[0] / received_Co_peak1;
        uefficiency_Co_peak1 = efficiency_Co_peak1 * TMath::Sqrt(1/area[0] + 1/received_Co_peak1);
        efficiency_Co_peak2 = area[1] / received_Co_peak2;
        uefficiency_Co_peak2 = efficiency_Co_peak2 * TMath::Sqrt(1/area[1] + 1/received_Co_peak2);




        // Am source

        xMin = 0;
        xMax = 1000;


        TCanvas *c3 = new TCanvas("c3", "Am calibrated spectrum");
        TH1D *hist_Am = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/day1/10min_Am.root", chan, bins, xMin, xMax);

        new_end = cal_NaI(xMax);
        TH1D *hcal_Am = new TH1D("hist_cal", "NaI (Am source) spectrum", bins, cal_NaI(xMin), cal_NaI(xMax)); //create calibrated histo!!
        for(Int_t i=0; i<bins; i++) hcal_Am->SetBinContent(i+1, hist_Am->GetBinContent(i+1)); //fill cal histo with same content!
        
        hcal_Am->Rebin(16);
        hcal_Am->SetTitle("");
        hcal_Am->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_Am->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_Am->GetXaxis()->GetBinWidth(0)));
        hcal_Am->GetXaxis()->SetRangeUser(xMin, xMax);
        hcal_Am->GetYaxis()->SetTitleOffset(1.25);
        hcal_Am->SetLineColor(kBlue);
        hcal_Am->SetFillColorAlpha(kAzure+7, 0.4);
        hcal_Am->Draw();


        ss = new TSpectrum(1); //creating TSpectrum object
        TH1D *hclone_Am = (TH1D*)hcal_Am->Clone();
        hb = ss->Background(hclone_Am,20,"same"); //calculating backgraund (TSpecrum class method)

        TCanvas *c4 = new TCanvas("c4", "Subtracted background spectrum");
        hclone_Am->Add(hclone_Am, hb, 1, -1); //substracting background from initial spectrum
        TH1D *hfinal_Am = (TH1D*)hclone_Am->Clone();
        for (Int_t i = 0; i<hfinal_Am->GetEntries(); i++)
        {
            if (hfinal_Am->GetBinContent(i) < 0) hfinal_Am->SetBinContent(i, 0);
        }
        hfinal_Am->Draw();
        

        
        Float_t mean_Am, sigma_Am, area_Am;
        Float_t xp = 59.54092;
        
        fun1 = new TF1("fun1name", "gaus", xp-120,xp+120);

        hfinal_Am->Fit(fun1,"RQN");

        fun2 = new TF1("fun2name","gaus", fun1->GetParameter(1)-(3.0*fun1->GetParameter(2)),fun1->GetParameter(1)+(3.0*fun1->GetParameter(2)));
        fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
        fun2->SetLineColor(kRed);

        hfinal_Am->Fit(fun2, "RQN+");
        



        mean_Am = fun2->GetParameter(1);
        sigma_Am = fun2->GetParameter(2);

        clone_fun = (TF1*)fun2->Clone();
        clone_fun->SetRange(mean_Am-3*sigma_Am, mean_Am+3*sigma_Am);

        clone_fun->SetLineColor(kRed);
        clone_fun->SetFillColorAlpha(kRed, 0.3);
        clone_fun->SetFillStyle(1001);

        clone_fun->Draw("SAME");

        area_Am = fun2->Integral(mean_Am - 3 * sigma_Am, mean_Am + 3 * sigma_Am) / hfinal_Am->GetXaxis()->GetBinWidth(0);
        efficiency_Am = area_Am / received_Am;
        uefficiency_Am = efficiency_Am*TMath::Sqrt(1/area_Am + 1/received_Am);


        Float_t efficiency[3] = {efficiency_Am, efficiency_Co_peak1, efficiency_Co_peak2};
        Float_t uefficiency[3] = {uefficiency_Am, uefficiency_Co_peak1, uefficiency_Co_peak2};
        Float_t energy[3]  = {59.54092, 1173.228, 1332.514};

        
        

    }
    
    if (detector == "HPGe")
    {   

        received_Co_peak1 = emitted_Co_peak1 * acceptancy_HPGe;
        received_Co_peak2 = emitted_Co_peak2 * acceptancy_HPGe;
        received_Am = emitted_Am * acceptancy_HPGe;

        xMin = 0, xMax = 10000;
        short chan = 1;









        TCanvas *c1 = new TCanvas("c1", "Co calibrated spectrum");
        TH1D *hist_Co = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/HPGe_Co10min.root", chan, bins, xMin, xMax);

        Float_t new_end = cal_HPGe(xMax);
        TH1D *hcal_Co = new TH1D("hist_cal", "HPGe (Co source) spectrum", bins, cal_HPGe(xMin), cal_HPGe(xMax)); //create calibrated histo!!
        for(Int_t i = 0; i < bins; i++) hcal_Co->SetBinContent(i+1, hist_Co->GetBinContent(i+1)); //fill cal histo with same content!
        
        hcal_Co->Rebin(4);
        hcal_Co->SetTitle("");
        hcal_Co->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_Co->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_Co->GetXaxis()->GetBinWidth(0)));
        hcal_Co->GetXaxis()->SetRangeUser(xMin, xMax);
        hcal_Co->GetYaxis()->SetTitleOffset(1.25);
        hcal_Co->SetLineColor(kBlue);
        hcal_Co->SetFillColorAlpha(kAzure+7, 0.4);
        hcal_Co->Draw();


        TSpectrum *ss = new TSpectrum(2); //creating TSpectrum object
        TH1D *hclone_Co = (TH1D*)hcal_Co->Clone();
        TH1 *hb = ss->Background(hclone_Co,20,"same"); //calculating backgraund (TSpecrum class method)

        TCanvas *c2 = new TCanvas("c2", "Subtracted background spectrum");
        hclone_Co->Add(hclone_Co,hb,1,-1); //substracting background from initial spectrum
        TH1D *hfinal_Co = (TH1D*)hclone_Co->Clone();
        for (Int_t i = 0; i<hfinal_Co->GetEntries(); i++)
        {
            if (hfinal_Co->GetBinContent(i) < 0) hfinal_Co->SetBinContent(i, 0);
        }
        hfinal_Co->Draw();
        
        Float_t xpeaks[2] = {1173.228, 1332.514};

        Float_t mean_vec[2] = {};
        Float_t sigma_vec[2] = {};
        Float_t area[2] = {};



        TF1 *fun1, *fun2, *clone_fun;
        for(Int_t i=0; i<2; i++)
        {
            Float_t xp = xpeaks[i];
            //cout << " ************ Peak at: " << xp << endl;

            fun1 = new TF1("fun1name", "gaus", xp-100,xp+100);

            hfinal_Co->Fit(fun1,"RQN");

            fun2 = new TF1("fun2name","gaus", fun1->GetParameter(1)-(3.*fun1->GetParameter(2)), fun1->GetParameter(1)+(3.0*fun1->GetParameter(2)));
            fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
            fun2->SetLineColor(kRed);

            hfinal_Co->Fit(fun2,"RQN+");
            



            mean_vec[i] = fun2->GetParameter(1);
            sigma_vec[i] = fun2->GetParameter(2);

            clone_fun = (TF1*)fun2->Clone();
            clone_fun->SetRange(mean_vec[i]-3*sigma_vec[i], mean_vec[i]+3*sigma_vec[i]);

            clone_fun->SetLineColor(kRed);
            clone_fun->SetFillColorAlpha(kRed, 0.3);
            clone_fun->SetFillStyle(1001);

            clone_fun->Draw("SAME");

            area[i] = fun2->Integral(mean_vec[i] - 3 * sigma_vec[i], mean_vec[i] + 3 * sigma_vec[i]) / hfinal_Co->GetXaxis()->GetBinWidth(0);

        }

        efficiency_Co_peak1 = area[0] / received_Co_peak1;
        uefficiency_Co_peak1 = efficiency_Co_peak1 * TMath::Sqrt(1/area[0] + 1/received_Co_peak1);
        efficiency_Co_peak2 = area[1] / received_Co_peak2;
        uefficiency_Co_peak2 = efficiency_Co_peak2 * TMath::Sqrt(1/area[1] + 1/received_Co_peak2);




        // Am source



        TCanvas *c3 = new TCanvas("c3", "Am calibrated spectrum");
        TH1D *hist_Am = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/HPGe_Am10min.root", chan, bins, xMin, xMax);

        new_end = cal_HPGe(xMax);
        TH1D *hcal_Am = new TH1D("hist_cal", "HPGe (Am source) spectrum", bins, cal_HPGe(xMin), cal_HPGe(xMax)); //create calibrated histo!!
        for(Int_t i=0; i<bins; i++) hcal_Am->SetBinContent(i+1, hist_Am->GetBinContent(i+1)); //fill cal histo with same content!
        
        hcal_Am->Rebin(16);
        hcal_Am->SetTitle("");
        hcal_Am->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_Am->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_Am->GetXaxis()->GetBinWidth(0)));
        hcal_Am->GetXaxis()->SetRangeUser(xMin, xMax);
        hcal_Am->GetYaxis()->SetTitleOffset(1.25);
        hcal_Am->SetLineColor(kBlue);
        hcal_Am->SetFillColorAlpha(kAzure+7, 0.4);
        hcal_Am->Draw();


        ss = new TSpectrum(1); //creating TSpectrum object
        TH1D *hclone_Am = (TH1D*)hcal_Am->Clone();
        hb = ss->Background(hclone_Am,20,"same"); //calculating backgraund (TSpecrum class method)

        TCanvas *c4 = new TCanvas("c4", "Subtracted background spectrum");
        hclone_Am->Add(hclone_Am, hb, 1, -1); //substracting background from initial spectrum
        TH1D *hfinal_Am = (TH1D*)hclone_Am->Clone();
        for (Int_t i = 0; i<hfinal_Am->GetEntries(); i++)
        {
            if (hfinal_Am->GetBinContent(i) < 0) hfinal_Am->SetBinContent(i, 0);
        }
        hfinal_Am->Draw();
        

        Float_t mean_Am, sigma_Am, area_Am;
        Float_t xp = 59.54092;
        
        fun1 = new TF1("fun1name", "gaus", xp-120,xp+120);

        hfinal_Am->Fit(fun1,"RQN");

        fun2 = new TF1("fun2name","gaus", fun1->GetParameter(1)-(3.0*fun1->GetParameter(2)),fun1->GetParameter(1)+(3.0*fun1->GetParameter(2)));
        fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
        fun2->SetLineColor(kRed);

        hfinal_Am->Fit(fun2, "RQN+");
        



        mean_Am = fun2->GetParameter(1);
        sigma_Am = fun2->GetParameter(2);

        clone_fun = (TF1*)fun2->Clone();
        clone_fun->SetRange(mean_Am-3*sigma_Am, mean_Am+3*sigma_Am);

        clone_fun->SetLineColor(kRed);
        clone_fun->SetFillColorAlpha(kRed, 0.3);
        clone_fun->SetFillStyle(1001);

        clone_fun->Draw("SAME");

        area_Am = fun2->Integral(mean_Am - 3 * sigma_Am, mean_Am + 3 * sigma_Am) / hfinal_Am->GetXaxis()->GetBinWidth(0);
        efficiency_Am = area_Am / received_Am;
        cout << "Efficiency Am: " << efficiency_Am << endl;
        uefficiency_Am = efficiency_Am*TMath::Sqrt(1/area_Am + 1/received_Am);
    }

    cout << "Efficiency Am: " << efficiency_Am << endl;
    Float_t efficiency[3] = {efficiency_Am*100, efficiency_Co_peak1*100, efficiency_Co_peak2*100};
    Float_t uefficiency[3] = {uefficiency_Am*100, uefficiency_Co_peak1*100, uefficiency_Co_peak2*100};
    Float_t energy[3] = {59.54092, 1173.228, 1332.514};

    TCanvas *c5 = new TCanvas("c5", "Intrinsic efficiency");
    c5->SetLeftMargin(0.11749);
    c5->SetRightMargin(0.0544413);
    c5->SetBottomMargin(0.101053);
    c5->SetTopMargin(0.0336842);
    TH2F *h = new TH2F("h", ";Photon energy [keV];Intrinsic efficiency [%]",  1000, 0, 1400, 1000, yMin, yMax);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->Draw();
    TGraphErrors *efficiency_graph = new TGraphErrors(3, energy, efficiency, 0, uefficiency); //calibration: channel vs energy of the found peaks
    efficiency_graph->SetMarkerStyle(20);
    efficiency_graph->SetMarkerColor(1);
    h->GetXaxis()->SetRangeUser(0, 1400);
    efficiency_graph->Draw("P");


    gSystem->Exec("mkdir preliminary_plots/");
	gSystem->Exec("mkdir preliminary_plots/efficiency");
	TString save_route = "preliminary_plots/efficiency/";
	c5->SaveAs(save_route + string(detector) + "_efficiencyCurve.pdf");
}