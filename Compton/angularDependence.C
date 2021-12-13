/// \author Javier Mariño Villadamigo


void angularDependence()
{   

    //Style settings
    gStyle->SetPalette(1);

    gStyle->SetTitleFont(62, "XYZ"); // title of axes
    gStyle->SetTitleSize(0.035, "XYZ"); // size of axes titles

    gStyle->SetTitleFont(62, "t"); // size of axes titles
    gStyle->SetTitleFontSize(0.08); // size of axes titles

    gStyle->SetLabelFont(62, "XYZ"); // labels of axes
    gStyle->SetLabelSize(0.045, "XYZ"); // size of axes labels

    gStyle->SetTitleOffset(1.4, "X"); // offset of x-axis
    gStyle->SetTitleOffset(1.2, "Y"); // offset of y-axis
    gStyle->SetTitleOffset(1.8, "Z"); // offset of z-axis
    gStyle->SetTextFont(62); // general texts

    gStyle->SetOptStat(0); // clearing TPaveStats




    Double_t deltaTheta = (TMath::ATan(4. / 3.757) + TMath::ATan(8. / 26.)) / TMath::Pi() * 180; //TMath::ATan(4. / 26.) / TMath::Pi() * 180; // geometry definition for 260 mm of distance between SCATTERER and DETECTOR


    // peaks from the calibrated and filtered data (this was done by FitPanel manually, since it is more appropiate to fit visually in order to not take spurious contributions)
    Double_t angles[6] = {0, 20, 40, 60, 90, 100};
    Double_t sigma_angles[6];
    for (Int_t i = 0; i<6; i++)
    {
        sigma_angles[i] = deltaTheta;
    }
    Double_t Ee_vec[6] = {35.6702, 55.1499, 93.5697, 155.887, 226.388, 246.219};
    Double_t sigma_Ee_vec[6] = {10.3801, 30.3350, 30.5068, 52.1386, 52.9440, 39.3669};
    Double_t Ef_vec[6] = {458.064, 441.250, 396.069, 349.047, 279.245, 267.993};
    Double_t sigma_Ef_vec[6] = {21.2928, 27.2681, 41.5107, 35.5887, 34.2632, 30.7702};

    TCanvas *c1 = new TCanvas("c1", "Angular dependence");
    c1->SetWindowSize(1920, 643);
    c1->Divide(2, 1);
    c1->cd(1)->SetMargin(0.088873, 0.04458308, 0.100197, 0.0991254);
    c1->cd(1);
    TH2F *h = new TH2F("h", ";Scattering angle [deg];Energy [keV]",  1000, -5, 115, 1000, -10, 550); // initialising the histogram with a 10-day margin
    h->Draw();


    TGraphErrors *graph_Ef = new TGraphErrors(6, angles, Ef_vec, 0, sigma_Ef_vec);
    graph_Ef->SetMarkerStyle(21);
    graph_Ef->SetMarkerColor(2);
    graph_Ef->Draw("PSAME");

    TF1 *Ef_th = new TF1("Ef_th", "511 / (2 - TMath::Cos(x / 180 * TMath::Pi()))", 0, 110);
    Ef_th->SetLineColor(2);
    Ef_th->Draw("SAME");

    TF1 *Ef_th_fit = new TF1("Ee_th", "[0] / (2 - TMath::Cos((x) / 180 * TMath::Pi()))", 0, 110);
    Ef_th_fit->SetLineColor(2);
    Ef_th_fit->SetLineStyle(2);
    graph_Ef->Fit(Ef_th_fit, "EX0");
    Ef_th_fit->Draw("SAME");


    
    TGraphErrors *graph_Ee = new TGraphErrors(6, angles, Ee_vec, 0, sigma_Ee_vec);
    graph_Ee->SetMarkerStyle(21);
    graph_Ee->SetMarkerColor(4);
    graph_Ee->Draw("P");

    TF1 *Ee_th = new TF1("Ee_th", "511 * (1 - 1 / (2 - TMath::Cos(x / 180 * TMath::Pi())))", 0, 110);
    Ee_th->SetLineColor(4);
    Ee_th->Draw("SAME");

    TF1 *Ee_th_fit = new TF1("Ee_th", "[0] * (1 - 1 / (2 - TMath::Cos(x / 180 * TMath::Pi())))", 0, 110);
    Ee_th_fit->SetLineColor(4);
    Ee_th_fit->SetLineStyle(2);
    graph_Ee->Fit(Ee_th_fit, "EX0");
    Ee_th_fit->Draw("SAME");

    TLatex text3;
	text3.SetTextSize(0.070);
    text3.DrawLatex(70, 470, Form("#color[2]{E_{#gamma} = %1.0f(%2.0f)}", Ef_th_fit->GetParameter(0), Ef_th_fit->GetParErrors()[0]));

    TLatex text4;
	text4.SetTextSize(0.070);
    text4.DrawLatex(70, 50, Form("#color[4]{E_{#gamma} = %1.0f(%2.0f)}", Ee_th_fit->GetParameter(0), Ee_th_fit->GetParErrors()[0]));

    TLegend *legend = new TLegend(0.129062, 0.315072, 0.477733, 0.602371);
    legend->SetTextSize(0.034);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    //legend->AddEntry(h1,"Histogram filled with random numbers","f");
    legend->AddEntry(graph_Ef, "Photon energy data", "p"); //l draws the line;
    legend->AddEntry(Ef_th, "Theoretical photon energy", "l");
    legend->AddEntry(Ef_th_fit, "Fit photon energy", "l");
    legend->AddEntry(graph_Ee, "Electron energy data", "p"); // p draws the polymarker (see https://root.cern.ch/doc/master/classTLegend.html)
    legend->AddEntry(Ee_th, "Theoretical electron energy", "l");
    legend->AddEntry(Ee_th_fit, "Fit electron energy", "l");
    legend->Draw();







    



    Double_t res_Ee[6] = {};
    Double_t res_Ef[6] = {};
    Double_t chi2_Ee = 0, chi2_Ef = 0;
    for (Int_t i = 0; i<6; i ++)
    {
        res_Ee[i] = Ee_vec[i] - Ee_th_fit->Eval(angles[i]);
        chi2_Ee += TMath::Power((Ee_vec[i] - Ee_th_fit->Eval(angles[i])) / sigma_Ee_vec[i], 2);
        res_Ef[i] = (Ef_vec[i] - Ef_th_fit->Eval(angles[i]));
        chi2_Ef += TMath::Power((Ef_vec[i] - Ef_th_fit->Eval(angles[i])) / sigma_Ef_vec[i], 2);
    }


    cout << "chi2/NDf photon: " << chi2_Ef/5 << " chi2/NDf electron: " << chi2_Ee/5 << endl;
    c1->cd(2)->SetMargin(0.088873, 0.04458308, 0.100197, 0.0991254);
    c1->cd(2);
    TH2F *h2 = new TH2F("h2", ";Scattering angle [deg];Residuals [keV]",  1000, -5, 115, 1000, -30, 60); // initialising the histogram with a 10-day margin
    h2->Draw();

    TGraph *graph_res_Ee = new TGraph(6, angles, res_Ee);
    graph_res_Ee->SetMarkerStyle(21);
    graph_res_Ee->SetMarkerColor(2);
    graph_res_Ee->Draw("P");

    TGraph *graph_res_Ef = new TGraph(6, angles, res_Ef);
    graph_res_Ef->SetMarkerStyle(21);
    graph_res_Ef->SetMarkerColor(4);
    graph_res_Ef->Draw("P");

    TLine *zero = new TLine(0, 0, 110, 0);
    zero->SetLineStyle(9);
    zero->SetLineColor(1);
    zero->SetLineWidth(2);
    zero->Draw("same");

    TLatex text;
	text.SetTextSize(0.075);
    text.DrawLatex(37, -10, Form("#color[2]{#chi^{2}/NDf = %1.1f}", chi2_Ef/5));

    TLatex text2;
	text2.SetTextSize(0.075);
    text2.DrawLatex(37, -20, Form("#color[4]{#chi^{2}/NDf = %1.1f}", chi2_Ee/5));

    TString save_route = "preliminary_plots/angularAnalysis/";
	mkdir(save_route, 0777);
    c1->SaveAs(save_route + "angularDependence.pdf");

}