#include "gethisto.C"

/// \author Javier Mariño Villadamigo

Double_t cal_NaI(Int_t ch) {return (-9.32 + ch*0.12081);}
Double_t cal_HPGe(Int_t ch) {return (0.27 + ch*0.17434);}
Float_t eff_HPGe(Float_t E) {return (-128.4 + 8227 / TMath::Sqrt(E) + 10180 / E);}
Float_t eff_NaI(Float_t E) {return (1.069 + 56.62 / TMath::Sqrt(E));}


void pellet() {

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

    


    // Geometry data
    Float_t distance_source_NaI = 2.5; // cm
    Float_t NaI_radius = 7.5/2.; // cm

    Float_t distance_source_HPGe = 7.5; // cm
    Float_t HPGe_radius = TMath::Sqrt(1200 / TMath::Pi()) / 10; // cm

    Float_t acceptancy_NaI = TMath::Pi() * TMath::Power(NaI_radius, 2) / TMath::Power(distance_source_NaI, 2);
    Float_t acceptancy_HPGe = TMath::Pi() * TMath::Power(HPGe_radius, 2) / TMath::Power(distance_source_HPGe, 2);
    
    Float_t time = 15;
    Float_t mass = 250;
    Float_t normalization_factor = time / 29.; // 29 minutes for background and 15 minutes for Pellet sample
    xMin = 0, xMax = 30000;
    xMin_en = 40, xMax_en = 2000;
    


    TH1D *hist_NaI = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/samples/pellet15min_ash250g.root", 0, bins, xMin, xMax);
    TH1D *hist_HPGe = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/samples/pellet15min_ash250g.root", 1, bins, xMin, xMax);
    TH1D *hb_NaI = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/background29min.root", 0, bins, xMin, xMax);
    TH1D *hb_HPGe = getHistoForChannelFromTree("/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/background29min.root", 1, bins, xMin, xMax);
    
    Float_t new_beg_NaI = cal_NaI(xMin);
    Float_t new_beg_HPGe = cal_HPGe(xMin);
    Float_t new_end_NaI = cal_NaI(xMax);
    Float_t new_end_HPGe = cal_HPGe(xMax);

    TH1D *hcal_NaI = new TH1D("hcal_NaI", "Calibrated histo of Pellet (NaI)", bins, new_beg_NaI, new_end_NaI);
    TH1D *hcal_HPGe = new TH1D("hcal_HPGe", "Calibrated histo of Pellet (HPGe)", bins, new_beg_HPGe, new_end_HPGe);
    TH1D *hbcal_NaI = new TH1D("hbcal_NaI", "Calibrated histo of Pellet (NaI)", bins, new_beg_NaI, new_end_NaI);
    TH1D *hbcal_HPGe = new TH1D("hbcal_HPGe", "Calibrated histo of Pellet (HPGe)", bins, new_beg_HPGe, new_end_HPGe);

    for (Int_t i = 0; i < bins; i++)
    {
        hcal_NaI->SetBinContent(i, hist_NaI->GetBinContent(i));
        hcal_HPGe->SetBinContent(i+1, hist_HPGe->GetBinContent(i+1));
        hbcal_NaI->SetBinContent(i, hb_NaI->GetBinContent(i) * normalization_factor);
        hbcal_HPGe->SetBinContent(i, hb_HPGe->GetBinContent(i) * normalization_factor);
    }
    

    TCanvas *c1 = new TCanvas("c1", "Pellet calibrated spectrum (NaI)");
    c1->SetWindowSize(1102, 771);
    c1->SetLeftMargin(0.11749);
    c1->SetRightMargin(0.0344413);
    c1->SetBottomMargin(0.101053);
    c1->SetTopMargin(0.0336842);
    hcal_NaI->Rebin(16);
    hcal_NaI->SetTitle("");
    hcal_NaI->GetXaxis()->SetTitle("Photon energy [keV] "); hcal_NaI->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal_NaI->GetXaxis()->GetBinWidth(0)));
    hcal_NaI->GetXaxis()->SetRangeUser(xMin_en, xMax_en);
    hcal_NaI->GetYaxis()->SetTitleOffset(1.);
    hcal_NaI->SetLineColor(kBlue);
    hcal_NaI->SetFillColorAlpha(kAzure+7, 0.4);
    hcal_NaI->Draw();

    hbcal_NaI->Rebin(16);

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


    TCanvas *c3 = new TCanvas("c3", "Pellet calibrated spectrum (HPGe)");
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


    Float_t P_gamma[5] = {0.055, 0.44, 0.18, 0.36, 0.1067};
    Float_t real_xpeaks[5] = {186, 239, 295, 352, 1461};
    Float_t xpeaks[5] = {192, 246, 295, 352, 1461};

    /*
    Float_t mean_NaI[5] = {};
    Float_t sigma_NaI[5] = {};


    Float_t Nd_NaI[5] = {};
    Float_t A_NaI[5] = {};
    
    c2->cd();
    TF1 *fun1, *fun2, *clone_fun;
        for(Int_t i=0; i<5; i++)
        {
            Float_t xp = xpeaks[i];
            //cout << " ************ Peak at: " << xp << endl;

            fun1 = new TF1("fun1name", "gaus", xp-20,xp+20);

            hNaI->Fit(fun1,"RQN");

            fun2 = new TF1("fun2name","gaus+pol1(3)", fun1->GetParameter(1)-(1.5*fun1->GetParameter(2)), fun1->GetParameter(1)+(1.5*fun1->GetParameter(2)));
            fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
            fun2->SetLineColor(kRed);

            hNaI->Fit(fun2,"RQN+");
            



            mean_NaI[i] = fun2->GetParameter(1);
            sigma_NaI[i] = fun2->GetParameter(2);

            clone_fun = (TF1*)fun2->Clone();
            clone_fun->SetRange(mean_NaI[i]-3*sigma_NaI[i], mean_NaI[i]+3*sigma_NaI[i]);

            clone_fun->SetLineColor(kRed);
            clone_fun->SetFillColorAlpha(kRed, 0.3);
            clone_fun->SetFillStyle(1001);

            clone_fun->Draw("SAME");

            Nd_NaI[i] = fun2->Integral(mean_NaI[i] - 3 * sigma_NaI[i], mean_NaI[i] + 3 * sigma_NaI[i]) / hNaI->GetXaxis()->GetBinWidth(0);
            cout << "\n\n********** " << real_xpeaks[i] << " keV peak **********\n" << endl;
            cout << "Mean = " << mean_NaI[i] << " +/- " << fun2->GetParErrors()[1] << endl;
            cout << "Sigma = " << sigma_NaI[i] << " +/- " << fun2->GetParErrors()[2] << endl;
            cout << "Integral = " << Nd_NaI[i] << " +/- " << TMath::Sqrt(Nd_NaI[i]) << endl;
            cout << "BR = " << P_gamma[i] << endl;
            A_NaI[i] = 4 * TMath::Pi() / acceptancy_NaI / P_gamma[i] / (time * 60) * Nd_NaI[i] / (eff_NaI(mean_NaI[i]) / 100) / mass; // we have divided the efficiency by 100 to obtain a (0-1) number
            cout << "Mass activity in NaI detector: " << A_NaI[i] << " Bq/g" << endl;
        }


        Float_t BR_Cs = 0.946, BR_Po = 0.449;
        Float_t xp = 650;
        //cout << " ************ Peak at: " << xp << endl;

        fun1 = new TF1("fun1name", "gaus+gaus(3)", xp-80,xp+60);
        fun1->SetParameters(50, 609.317, 5, 50, 661.659, 5);
        fun1->SetParLimits(1, 600, 620);
        fun1->SetParLimits(4, 650, 670);
        hNaI->Fit(fun1,"RQN");

        fun2 = new TF1("fun2name","gaus+gaus(3)", 570, 710);
        fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2), fun1->GetParameter(3), fun1->GetParameter(4), fun1->GetParameter(5));
        fun2->SetLineColor(kRed);
        hNaI->Fit(fun2,"RQN+");



        clone_fun = (TF1*)fun2->Clone();
        clone_fun->SetRange(570, 710);

        clone_fun->SetLineColor(kRed);
        clone_fun->SetFillColorAlpha(kRed, 0.3);
        clone_fun->SetFillStyle(1001);

        clone_fun->Draw("SAME");

        

        TF1 *funPo = new TF1("funPo", "gaus", fun2->GetParameter(1) - 3 * fun1->GetParameter(2), fun2->GetParameter(1) + 3 * fun1->GetParameter(2)); // separate the gaussians to do the integration
        funPo->SetParameters(fun2->GetParameter(0), fun2->GetParameter(1), fun2->GetParameter(2));
        TF1 *funCs = new TF1("funPo", "gaus", fun2->GetParameter(4) - 3 * fun1->GetParameter(5), fun2->GetParameter(4) + 3 * fun1->GetParameter(5)); // separate the gaussians to do the integration
        funCs->SetParameters(fun2->GetParameter(3), fun2->GetParameter(4), fun2->GetParameter(5));


        Float_t Nd_Po_NaI = funPo->Integral(fun2->GetParameter(1) - 3 * fun1->GetParameter(2), fun2->GetParameter(1) + 3 * fun1->GetParameter(2)) / hNaI->GetXaxis()->GetBinWidth(0);
        Float_t Nd_Cs_NaI = funCs->Integral(fun2->GetParameter(4) - 3 * fun1->GetParameter(5), fun2->GetParameter(4) + 3 * fun1->GetParameter(5)) / hNaI->GetXaxis()->GetBinWidth(0);

        cout << "\n\n********** 609 keV peak **********\n" << endl;
        cout << "Mean = " << fun2->GetParameter(1) << " +/- " << fun2->GetParErrors()[1] << endl;
        cout << "Sigma = " << fun2->GetParameter(2) << " +/- " << fun2->GetParErrors()[2] << endl;
        cout << "Integral = " << Nd_Po_NaI << " +/- " << TMath::Sqrt(Nd_Po_NaI) << endl;
        cout << "BR = " << BR_Po << endl;
        Float_t A_Po_NaI = 4 * TMath::Pi() / acceptancy_NaI / BR_Po / (time * 60) * Nd_Po_NaI / (eff_NaI(fun2->GetParameter(1)) / 100) / mass; // we have divided the efficiency by 100 to obtain a (0-1) number
        cout << "Mass activity of Po in NaI detector: " << A_Po_NaI << " Bq/g" << endl;


        cout << "\n\n********** 662 keV peak **********\n" << endl;
        cout << "Mean = " << fun2->GetParameter(4) << " +/- " << fun2->GetParErrors()[4] << endl;
        cout << "Sigma = " << fun2->GetParameter(5) << " +/- " << fun2->GetParErrors()[5] << endl;
        cout << "Integral = " << Nd_Cs_NaI << " +/- " << TMath::Sqrt(Nd_Cs_NaI) << endl;
        cout << "BR = " << BR_Cs << endl;
        Float_t A_Cs_NaI = 4 * TMath::Pi() / acceptancy_NaI / BR_Cs / (time * 60) * Nd_Cs_NaI / (eff_NaI(fun2->GetParameter(4)) / 100) / mass; // we have divided the efficiency by 100 to obtain a (0-1) number
        cout << "Mass activity of Cs in NaI detector: " << A_Cs_NaI << " Bq/g" << endl;
        */




    Float_t mean_HPGe[5] = {};
    Float_t sigma_HPGe[5] = {};


    Float_t Nd_HPGe[5] = {};
    Float_t A_HPGe[5] = {};
    
    c4->cd();
    TF1 *fun1, *fun2, *clone_fun;
        for(Int_t i=0; i<5; i++)
        {
            Float_t xp = xpeaks[i];
            //cout << " ************ Peak at: " << xp << endl;

            if (i == 1)
            {
                fun1 = new TF1("fun1name", "gaus", 215, 257);
            }
            else if (i == 0)
            {
                fun1 = new TF1("fun1name", "gaus", 173, 201);
                fun1->SetParLimits(1, 180, 190);
                cout << "pinche puto" << endl;
            }
            else
            {
                fun1 = new TF1("fun1name", "gaus", xp-20,xp+20);
            }
            

            hHPGe->Fit(fun1,"RQNWW");

            fun2 = new TF1("fun2name","gaus+pol1(3)", fun1->GetParameter(1)-20, fun1->GetParameter(1)+20);

            fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
            
            fun2->SetParLimits(1, fun1->GetParameter(1) - 40, fun1->GetParameter(1) + 40);
            fun2->SetParLimits(2, 0, 20);
            
            fun2->SetLineColor(kRed);

            hHPGe->Fit(fun2,"RQNWW+");
            



            mean_HPGe[i] = fun2->GetParameter(1);
            sigma_HPGe[i] = fun2->GetParameter(2);

            clone_fun = (TF1*)fun2->Clone();
            clone_fun->SetRange(mean_HPGe[i]-3*sigma_HPGe[i], mean_HPGe[i]+3*sigma_HPGe[i]);

            clone_fun->SetLineColor(kRed);
            clone_fun->SetFillColorAlpha(kRed, 0.3);
            clone_fun->SetFillStyle(1001);

            clone_fun->Draw("SAME");

            Nd_HPGe[i] = fun2->Integral(mean_HPGe[i] - 3 * sigma_HPGe[i], mean_HPGe[i] + 3 * sigma_HPGe[i]) / hHPGe->GetXaxis()->GetBinWidth(0);

            cout << "\n\n********** " << real_xpeaks[i] << " keV peak **********\n" << endl;
            cout << "Mean = " << mean_HPGe[i] << " +/- " << fun2->GetParErrors()[1] << endl;
            cout << "Sigma = " << sigma_HPGe[i] << " +/- " << fun2->GetParErrors()[2] << endl;
            cout << "Integral = " << Nd_HPGe[i] << " +/- " << TMath::Sqrt(Nd_HPGe[i]) << endl;
            cout << "BR = " << P_gamma[i] << endl;
            A_HPGe[i] = 4 * TMath::Pi() / acceptancy_HPGe / P_gamma[i] / (time * 60) * Nd_HPGe[i] / (eff_HPGe(mean_HPGe[i]) / 100) / mass; // we have divided the efficiency by 100 to obtain a (0-1) number
            cout << "Mass activity in HPGe detector: " << A_HPGe[i] << " Bq/g" << endl;

        }

        Float_t BR_Cs = 0.946, BR_Po = 0.449;
        Float_t xp = 650;
        //cout << " ************ Peak at: " << xp << endl;

        fun1 = new TF1("fun1name", "gaus+gaus(3)", xp-80,xp+60);
        fun1->SetParameters(50, 609.317, 5, 50, 661.659, 5);
        fun1->SetParLimits(1, 600, 620);
        fun1->SetParLimits(4, 650, 670);
        hHPGe->Fit(fun1,"RQN");

        fun2 = new TF1("fun2name","gaus+gaus(3)", 570, 710);
        fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2), fun1->GetParameter(3), fun1->GetParameter(4), fun1->GetParameter(5));
        fun2->SetLineColor(kRed);
        hHPGe->Fit(fun2,"RQN+");



        clone_fun = (TF1*)fun2->Clone();
        clone_fun->SetRange(570, 710);

        clone_fun->SetLineColor(kRed);
        clone_fun->SetFillColorAlpha(kRed, 0.3);
        clone_fun->SetFillStyle(1001);

        clone_fun->Draw("SAME");

        

        TF1 *funPo = new TF1("funPo", "gaus", fun2->GetParameter(1) - 3 * fun1->GetParameter(2), fun2->GetParameter(1) + 3 * fun1->GetParameter(2)); // separate the gaussians to do the integration
        funPo->SetParameters(fun2->GetParameter(0), fun2->GetParameter(1), fun2->GetParameter(2));
        TF1 *funCs = new TF1("funPo", "gaus", fun2->GetParameter(4) - 3 * fun1->GetParameter(5), fun2->GetParameter(4) + 3 * fun1->GetParameter(5)); // separate the gaussians to do the integration
        funCs->SetParameters(fun2->GetParameter(3), fun2->GetParameter(4), fun2->GetParameter(5));


        Float_t Nd_Po_HPGe = funPo->Integral(fun2->GetParameter(1) - 3 * fun1->GetParameter(2), fun2->GetParameter(1) + 3 * fun1->GetParameter(2)) / hHPGe->GetXaxis()->GetBinWidth(0);
        Float_t Nd_Cs_HPGe = funCs->Integral(fun2->GetParameter(4) - 3 * fun1->GetParameter(5), fun2->GetParameter(4) + 3 * fun1->GetParameter(5)) / hHPGe->GetXaxis()->GetBinWidth(0);

        cout << "\n\n********** 609 keV peak **********\n" << endl;
        cout << "Mean = " << fun2->GetParameter(1) << " +/- " << fun2->GetParErrors()[1] << endl;
        cout << "Sigma = " << fun2->GetParameter(2) << " +/- " << fun2->GetParErrors()[2] << endl;
        cout << "Integral = " << Nd_Po_HPGe << " +/- " << TMath::Sqrt(Nd_Po_HPGe) << endl;
        cout << "BR = " << BR_Po << endl;
        Float_t A_Po_HPGe = 4 * TMath::Pi() / acceptancy_HPGe / BR_Po / (time * 60) * Nd_Po_HPGe / (eff_HPGe(fun2->GetParameter(1)) / 100) / mass; // we have divided the efficiency by 100 to obtain a (0-1) number
        cout << "Mass activity of Po in HPGe detector: " << A_Po_HPGe << " Bq/g" << endl;


        cout << "\n\n********** 662 keV peak **********\n" << endl;
        cout << "Mean = " << fun2->GetParameter(4) << " +/- " << fun2->GetParErrors()[4] << endl;
        cout << "Sigma = " << fun2->GetParameter(5) << " +/- " << fun2->GetParErrors()[5] << endl;
        cout << "Integral = " << Nd_Cs_HPGe << " +/- " << TMath::Sqrt(Nd_Cs_HPGe) << endl;
        cout << "BR = " << BR_Cs << endl;
        Float_t A_Cs_HPGe = 4 * TMath::Pi() / acceptancy_HPGe / BR_Cs / (time * 60) * Nd_Cs_HPGe / (eff_HPGe(fun2->GetParameter(4)) / 100) / mass; // we have divided the efficiency by 100 to obtain a (0-1) number
        cout << "Mass activity of Cs in HPGe detector: " << A_Cs_HPGe << " Bq/g" << endl;

}