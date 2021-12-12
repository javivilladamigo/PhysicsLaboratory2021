#include "../../examples/data_analysis_examples/gethisto.C"

/// \author Javier Mariño Villadamigo

void fractionEvents(TString detector_to_calibrate = "TAGGER") {
	
	//Style settings
	gROOT->SetStyle("Default");
	gStyle->SetCanvasColor(0);
	gStyle->SetStatBorderSize(1);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	gStyle->SetStatColor(0);
	gStyle->SetPalette(0);
	

	Int_t topCh = int(TMath::Power(2, 15));
	Int_t bins = topCh;
	Int_t xMin = 0.;
	Int_t xMax = topCh;
	Int_t inferior_limit = 7640;
	Int_t left_limit = 0;
	Int_t right_limit = 9600;
	
	short chan;
	const char *route;
	if (detector_to_calibrate == "TAGGER")
	{
		route = "data/day1/511keV_fraction/TAGGER_Na.root";
		chan = 0;
	}
	

	TCanvas *c1 = new TCanvas("c1", "c1");
	gPad->SetMargin(0.12, 0.04, 0.09, 0.02); // left right bottom top

	TH1D *hist = getHistoForChannelFromTree(route, chan, bins, xMin, xMax);
	hist->Rebin(32); //rebinin for a better visualization and peak finding
	hist->SetTitle("");
	hist->SetName("22Na spectrum");
	hist->SetLineColor(kBlue);
	hist->SetFillColorAlpha(kAzure + 7, 0.4);
	hist->GetXaxis()->SetTitle("channel "); hist->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist->GetXaxis()->GetBinWidth(0)));

	hist->GetYaxis()->SetTitleOffset(1.5);
	


	
	hist->GetXaxis()->SetRangeUser(xMin, 25000);

	Double_t rangeFit[2] = {6500, 10000};
	TF1 *gaus_plus_bg = new TF1("Gaussian + BG", "gaus(0)+pol1(3)", rangeFit[0], rangeFit[1]);

	gaus_plus_bg->SetParameters(1, 8250);
	gaus_plus_bg->SetParLimits(1, 7500, 9000);
	gaus_plus_bg->SetParLimits(2, 0, 3000);

	hist->Fit(gaus_plus_bg, "R+");
	Double_t offset = gaus_plus_bg->GetParameter(3);
	Double_t slope = gaus_plus_bg->GetParameter(4);

	//TF1 *linear = new TF1("linear", "[0]+x*[1]", rangeFit[0], rangeFit[1]);
	TF1 *gaus_fit = new TF1("Gaussian Fit", "gaus", rangeFit[0], rangeFit[1]);

	gaus_fit->SetParameter(0, gaus_plus_bg->GetParameter(0));
	gaus_fit->SetParameter(1, gaus_plus_bg->GetParameter(1));
	gaus_fit->SetParameter(2, gaus_plus_bg->GetParameter(2));


	//linear->SetParameter(0, offset);
	//linear->SetParameter(1, slope);

	//linear->SetLineColor(kOrange);
	gaus_fit->SetLineColor(kRed);
	gaus_fit->SetFillColorAlpha(kRed, 0.3);
	gaus_fit->SetFillStyle(1001);
	//linear->SetLineWidth(4);
	//linear->Draw("SAME");
	gaus_fit->Draw("SAME");

	



	Double_t A_511 = gaus_fit->Integral(rangeFit[0], rangeFit[1]) / hist->GetXaxis()->GetBinWidth(0);
	Double_t A_tot = hist->GetEntries();
	Double_t fraction = A_511/A_tot;
	Double_t ufraction = fraction*TMath::Sqrt(1/A_511 + 1/A_tot);

	//grshade->Draw("SAME");

	cout << A_511 << " " << A_tot << endl;
	cout << "A(511)/Atot = " << fraction << " +/- " << ufraction << endl;
	TLatex text;
	text.SetTextSize(0.075);
    text.DrawLatex(11320.38, 30000, Form("#frac{#color[2]{A_{511}}}{#color[4]{A_{tot}}} #color[1]{= %1.4f(%1.0f)}", fraction, 1e4*ufraction));
    




	
	TLegend *legend = new TLegend(0.532951, 0.806316, 0.951289, 0.966316);
    legend->SetShadowColor(0);
    legend->AddEntry(gaus_plus_bg, "Gaussian + BG fit","l"); //l draws the line;
    legend->AddEntry(gaus_fit, "Gaussian fit","l"); // p draws the polymarker (see https://root.cern.ch/doc/master/classTLegend.html)
	legend->SetTextSize(0.0526316);
	legend->SetTextFont(62);
    legend->Draw();
    
	TString save_route = "preliminary_plots/fraction_511keV_events/";
	mkdir(save_route, 0777);
	c1->SaveAs(save_route + string(detector_to_calibrate) + "_fraction_511keV_events.pdf");

	/*
	std::ofstream data;
	if (detector_to_calibrate == "TAGGER") {data.open("calibration_data.csv", std::ofstream::out); data << "," << Form("%2.2f", xval[0]) << " keV," << Form("%2.2f", xval[1]) << " keV," << Form("%3.0f", xval[2]) << " keV," << Form("%4.2f", xval[3]) << " keV," << "a+b*ch" << endl; data << "#DETECTOR," << "Resolution(%)," << "Resolution(%)," << "Resolution(%)," << "Resolution(%)," << "a +/- ua," << "b +/- ub" << endl;}
	else if (detector_to_calibrate != "TAGGER") data.open("calibration_data.csv", std::ofstream::app);
	data << detector_to_calibrate << "," << Form("%2.2f,", resolution[0]) << Form("%2.2f,", resolution[1]) << Form("%2.2f,", resolution[2]) << Form("%2.2f,", resolution[3]) << Form("%1.2f +/- %1.2f,", a, ua) << Form("%1.5f +/- %1.5f,", b, ub) << endl;
	data.close();
	*/
}

