#include "../../examples/data_analysis_examples/gethisto.C"

/// \author Javier Mariño Villadamigo

void fractionEvents(TString detector_to_calibrate = "TAGGER") {
	
	//Style settings
	gROOT->SetStyle("Default");
	gStyle->SetCanvasColor(0);
	gStyle->SetStatBorderSize(1);
	gStyle->SetOptFit(0);
	//gStyle->SetOptStat(0);
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
	hist->GetXaxis()->SetTitle("channel"); hist->GetYaxis()->SetTitle("counts");

	hist->GetYaxis()->SetTitleOffset(1.5);
	
	

	TH1D *h_substracted = (TH1D*)hist->Clone();
	for (Int_t i=0; i<hist->GetNbinsX(); i++)
	{	
		if (hist->GetBinContent(i) >= inferior_limit) h_substracted->SetBinContent(i, hist->GetBinContent(i) - inferior_limit);
		else if (hist->GetBinContent(i) < inferior_limit) h_substracted->SetBinContent(i, 0);
	}

	TH1D *integral = (TH1D*)hist->Clone();
	TH1D *white = (TH1D*)integral->Clone();

	for (Int_t i=0; i<integral->GetNbinsX(); i++)
	{	
		if (integral->GetBinContent(i) >= inferior_limit) white->SetBinContent(i, inferior_limit);
		else if (integral->GetBinContent(i) < inferior_limit) white->SetBinContent(i, integral->GetBinContent(i));
	}


	integral->GetXaxis()->SetRangeUser(left_limit, right_limit);
	white->GetXaxis()->SetRangeUser(left_limit, right_limit);
	//integral->GetYaxis()->SetRangeUser(inferior_limit, hist->GetNbinsY());
	integral->SetLineWidth(0);
	white->SetLineWidth(0);
	white->SetFillColorAlpha(0, 1);
	integral->SetFillColorAlpha(kRed, 0.5);

	
	hist->GetXaxis()->SetRangeUser(xMin, 25000);
	hist->Draw();
	integral->Draw("SAME");
	white->Draw("SAME");
	
	Double_t A_511 = h_substracted->Integral(left_limit, right_limit);
	Double_t A_tot = hist->GetEntries();

	Double_t fraction = A_511/A_tot;
	Double_t ufraction = fraction*TMath::Sqrt(1/A_511+1/A_tot);

	cout << "A(511)/Atot = " << fraction << " +/- " << ufraction << endl;
	TLatex text;
	text.SetTextSize(0.075);
    text.DrawLatex(11320.38, 30000, Form("#frac{#color[2]{A_{511}}}{#color[4]{A_{tot}}} #color[1]{= %1.5f(%2.0f)}", fraction, ufraction*1e5));
    

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

