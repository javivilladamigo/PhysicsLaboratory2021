#include "gethisto.C"

/// \author Javier Mariño Villadamigo

Double_t cal_NaI(Int_t ch) {return (-9.32 + ch*0.12081);}
Double_t cal_HPGe(Int_t ch) {return (0.27 + ch*0.17434);}



void europium() {

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
    Int_t topCh = int(TMath::Power(2, 15));
	Int_t bins = topCh;
    Int_t xMin, xMax;
    short chan;
    Float_t new_end;

    const int npeaks = 16;
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

    hcal->Rebin(16);
    hcal->SetTitle("");
    hcal->GetXaxis()->SetTitle("Photon energy [keV] "); hcal->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal->GetXaxis()->GetBinWidth(0)));
    hcal->GetXaxis()->SetRangeUser(xMin, xMax);
    hcal->GetYaxis()->SetTitleOffset(1.25);
    hcal->SetLineColor(kBlue);
    hcal->SetFillColorAlpha(kAzure+7, 0.4);
    hcal->Draw();

    hb->Rebin(16);
    hb->SetLineColor(kRed);
    hb->SetFillColorAlpha(kRed, 0.4);
    hb->Draw("SAME");

    
    TCanvas *c2 = new TCanvas("c2", "Subtracted background spectrum");
    c2->SetLogy();
    TH1D *hclone = (TH1D*)hcal->Clone();
    hclone->Add(hclone, hb, 1, -1); //substracting background from initial spectrum
    hclone->Draw();
    for (Int_t i = 0; i<hclone->GetNbinsX(); i++) if (hclone->GetBinContent(i) < 0) hclone->SetBinContent(i, 0);

    Float_t xpeaks[npeaks] = {121.8, 244.7, 344.3, 411.1, 444, 586.3, 688.7, 778.9, 867.4, 964, 1085.8, 1112.1, 1212.9, 1299.1, 1408, 1457.6}; // 1089.7 has the contribution of also 1085.8
    Float_t intens[npeaks] = {141, 36.6, 127.2, 10.71, 15, 2.24, 4.12, 62.6, 20.54, 70.4, 48.7, 65, 6.67, 7.76, 100, 2.52};
    Float_t vert_pos[npeaks] = {20000, 7000, 9000, 1500, 1500, 500, 500, 3000, 1000, 3000, 1500, 1500, 300, 300, 3000, 80};


    TF1 *fun1, *fun2, *clone_fun;
    TArrow *peak_ind;
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

        //clone_fun->Draw("SAME");
        

        if (i == 0)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 8000, 0.005, "|>");
        }
        else if (i == 1)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 4000, 0.005, "|>");
        }
        else if (i == 2)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 4000, 0.005, "|>");
        }
        else if (i == 3)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 700, 0.005, "|>");
        }
        else if (i == 4)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 700, 0.005, "|>");
        }
        else if (i == 5)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 200, 0.005, "|>");
        }
        else if (i == 6)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 200, 0.005, "|>");
        }
        else if (i == 7)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 1800, 0.005, "|>");
        }
        else if (i == 8)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 500, 0.005, "|>");
        }
        else if (i == 9)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 2000, 0.005, "|>");
        }
        else if (i == 10)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 800, 0.005, "|>");
        }
        else if (i == 11)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 800, 0.005, "|>");
        }
        else if (i == 12)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 200, 0.005, "|>");
        }
        else if (i == 13)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 200, 0.005, "|>");
        }
        else if (i == 14)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 2000, 0.005, "|>");
        }
        else if (i == 15)
        {
            peak_ind = new TArrow(xp, vert_pos[i], xp, vert_pos[i] - 50, 0.005, "|>");
        }

        peak_ind->SetAngle(40);
        peak_ind->SetLineWidth(2);
        peak_ind->Draw();

        area[i] = fun2->Integral(mean_vec[i] - 3 * sigma_vec[i], mean_vec[i] + 3 * sigma_vec[i]) / hclone->GetXaxis()->GetBinWidth(0);
        cout << "For" << xp << " " << area[i] * 100 << endl;
    }

    Float_t rel_efficiency[npeaks] = {}, urel_efficiency[npeaks] = {};
    Float_t areas_Attilio[npeaks] = {73728.2, 13109.9, 33317.4, 2607.5, 3306.6, 398.5, 662, 7865, 2417.8, 7609.5, 5571.2, 6204.1, 619.7, 619.6, 7628.3, 156.9};
    
    for (Int_t i = 0; i<npeaks; i++)
    {
        rel_efficiency[i] = areas_Attilio[i] / areas_Attilio[npeaks-2] * 10000 / intens[i];
        
        urel_efficiency[i] = rel_efficiency[i] * TMath::Sqrt(1 / areas_Attilio[i] + 1 / areas_Attilio[npeaks-2]);
        cout << urel_efficiency[i] << endl;
    }





    //Float_t rel_ef_Attilio[npeaks] = {6.8546, 4.6956, 3.4336, 3.1916, 2.8897, 2.3321, 2.1063, 1.6470, 1.5431, 1.4169, 1.4996, 1.2512, 1.2179, 1.0467, 1.0000, 0.8162};



    TCanvas *c3 = new TCanvas("c3", "Relative efficiency");
    c3->SetLogy();
    c3->SetWindowSize(1102, 771);
    c3->SetLeftMargin(0.11749);
    c3->SetRightMargin(0.0344413);
    c3->SetBottomMargin(0.101053);
    c3->SetTopMargin(0.0336842);
    TH2F *h = new TH2F("h", ";Photon energy [keV];Efficiency [%]",  1000, 0, 1500, 1000, 2, 1000);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->Draw();

    TGraphErrors *rel_ef_graph = new TGraphErrors(npeaks, xpeaks, rel_efficiency, 0, urel_efficiency); //calibration: channel vs energy of the found peaks
    rel_ef_graph->SetMarkerStyle(20);
    rel_ef_graph->SetMarkerColor(1);
    rel_ef_graph->Draw("P");
    
    TF1 *square_law = new TF1("square_law", "[0] + [1] / TMath::Sqrt(x) + [2] / x", 0, 1500);
    square_law->SetLineColor(kRed);
    square_law->SetLineWidth(2);
    rel_ef_graph->Fit(square_law, "R+");
    gPad->Update();

    TPaveStats *st = (TPaveStats*)rel_ef_graph->FindObject("stats");
    st->SetX1NDC(0.216364); //new x start position
    st->SetY1NDC(0.439678); //new y start position
    st->SetX2NDC(0.577273); //new x end position
    st->SetY2NDC(0.647453); //new y end position

   
    TF1 *abs_eff = new TF1("abs_eff", "-0.5963 + 112.2 / TMath::Sqrt(x)", 0, 1500);
    abs_eff->SetLineColor(kBlue);
    abs_eff->SetLineWidth(2);
    abs_eff->Draw("SAME");
   
    TF1 *normalized = new TF1("normalized", "square_law * (2.39257 / 100)", 0, 1500);
    normalized->SetLineColor(kGreen + 2);
    normalized->SetLineWidth(2);
    normalized->Draw("SAME");
    
    cout << "The normalized parameters are" << endl;

    cout << "p0 = " << square_law->GetParameter(0) * (2.39257 / 100) << " +/- " << square_law->GetParErrors()[0] * (2.39257 / 100) << endl;
    cout << "p1 = " << square_law->GetParameter(1) * (2.39257 / 100) << " +/- " << square_law->GetParErrors()[1] * (2.39257 / 100) << endl;
    cout << "p2 = " << square_law->GetParameter(2) * (2.39257 / 100) << " +/- " << square_law->GetParErrors()[2] * (2.39257 / 100) << endl;

    TLegend *legend = new TLegend(0.586364,0.430976,0.84, 0.638751);
    legend->SetTextSize(0.025);
    legend->AddEntry(rel_ef_graph, "HPGe (Eu source) rel. efficiency", "p");
    legend->AddEntry(square_law, "Fit (Eu source) to rel. eff.","l");
    legend->AddEntry(abs_eff, "Fit (241Am + 60Co) to abs. eff.","l");
    legend->AddEntry(normalized, "Normalized absolute efficiency","l");
    legend->Draw("SAME");

   

    gSystem->Exec("mkdir preliminary_plots/");
	gSystem->Exec("mkdir preliminary_plots/efficiency");
	TString save_route = "preliminary_plots/efficiency/";
	c3->SaveAs(save_route + "Eu_efficiency.pdf");


}