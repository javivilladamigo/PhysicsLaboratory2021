#include "../../examples/data_analysis_examples/gethisto.C"

/// \author Javier Mariño Villadamigo


   


void calibration(TString detector_to_calibrate, Int_t npeaks = 3) {
	//Style settings
	gROOT->SetStyle("Default");
	gStyle->SetCanvasColor(0);
	gStyle->SetStatBorderSize(1);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(0);

	Int_t topCh = int(TMath::Power(2, 15));
	Int_t bins = topCh;
	Int_t xMin = 700.;
	Int_t xMax = topCh;
	
	short chan;
	const char *route;
	if (detector_to_calibrate == "TAGGER")
	{
		route = "data/day1/calibration/TAGGER_calibration_Na&Am.root";
		chan = 0;
	}
	else if (detector_to_calibrate == "SCATTERER")
	{
		route = "data/day1/calibration/SCATTERER_calibration_Na&Am.root";
		chan = 1;
	}
	else if (detector_to_calibrate == "DETECTOR")
	{
		route = "data/day1/calibration/DETECTOR_calibration_Na&Am.root";
		chan = 2;
	}
	

	TCanvas *c1 = new TCanvas("c1", "calibration_canvas");
	c1->SetLogy();

	TH1D *hist = getHistoForChannelFromTree(route, chan, bins, 0, xMax);
	
	hist->Rebin(16); //rebinin for a better visualization and peak finding
	hist->SetTitle("");
	hist->GetXaxis()->SetTitle("channel "); hist->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist->GetXaxis()->GetBinWidth(0)));
	hist->GetXaxis()->SetRangeUser(xMin, xMax);
	hist->GetYaxis()->SetTitleOffset(1.25);
	hist->SetLineColor(kBlue);
	hist->SetFillColorAlpha(kAzure+7, 0.4);

	

	TSpectrum *ss = new TSpectrum(npeaks); //creating TSpectrum object
	Int_t nfound = ss->Search(hist, 10, "goff", 0.005); //searching peaks (TSpecrum class method)
	//cout << "\n --------- Found " << nfound << " peaks" << endl << endl << endl;
	Double_t *xpeaks;
	xpeaks = (Double_t*) ss->GetPositionX(); //saving the position of the found peaks in an array

	Double_t xp = 0.;
	Float_t mean = 0., sigma = 0.;
	Float_t resolution[3]={0., 0., 0.};
	Float_t chi2 = 0.;
	Float_t Ndf = 0.;

	//Sort xpeaks array
	int n = sizeof(*xpeaks)/sizeof(xpeaks[0]);
	std::sort(xpeaks,xpeaks+nfound);

	hist->Draw();
	TF1 *fun1, *fun2, *fit;
	for(Int_t i=0;i<nfound;i++)
	{
		xp = xpeaks[i];
		//cout << " ************ Peak at: " << xp << endl;

		fun1 = new TF1("fun1name", "gaus", xp-120,xp+120);

		hist->Fit(fun1,"RQNP");

		fun2 = new TF1("fun2name","gaus(0)+pol1(3)", fun1->GetParameter(1)-(1.5*fun1->GetParameter(2)),fun1->GetParameter(1)+(1.5*fun1->GetParameter(2)));
		fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
		fun2->SetLineColor(kRed);

		hist->Fit(fun2,"RP+");

		mean  = fun2->GetParameter(1);
		sigma = fun2->GetParameter(2);

		fit = hist->GetFunction("fun2name");

		chi2 = fit->GetChisquare();
		Ndf = fit->GetNDF();

		resolution[i]=((2.35*sigma)/mean)*100;
		//cout << "\nRESULTS:  Mean at: " << mean << " with sigma " << sigma << " , FWHM of " << 2.35*sigma << " , and resolution of " << resolution[i] << " %" << endl << endl;
		//cout << "\nThe fit shows a Chi2/Ndf = " << chi2/Ndf << endl;
	}
		

	//Calibration
	Float_t xval[3];
	if(npeaks == 1) {xval[0] = 59.54092;} else if(npeaks == 2) {xval[0] = 59.54092; xval[1] = 511;} else if(npeaks == 3) {xval[0] = 59.54092; xval[1] = 511; xval[2] = 1274.537;};
	Float_t peaks[npeaks];
	for(Int_t i = 1; i<npeaks; i++) peaks[i] = xpeaks[i];
	TGraph *cal_graph = new TGraph(npeaks, peaks, xval); //calibration: channel vs energy of the found peaks
	cal_graph->SetNameTitle("cal_graph","");
	cal_graph->SetMarkerStyle(34);
	cal_graph->SetMarkerSize(1.5);
	cal_graph->SetMarkerColor(4);
	cal_graph->SetLineWidth(0);
	cal_graph->GetXaxis()->SetTitle("channel");
	
	cal_graph->GetYaxis()->SetTitle("Photon energy [keV]");
	cal_graph->GetYaxis()->SetTitleOffset(1.25);
	TF1 *cal_fn = new TF1("cal_fn","[0]*x+[1]"); //defining a first order polynomial function and fit the data
	cal_fn->SetLineColor(kRed);
	cal_graph->Fit(cal_fn, "F"); // specify "F" option to use Minuit for linear fits (by default, polN, are fitted by the linear fitter)
	Double_t a = cal_fn->GetParameter(1); Double_t ua = cal_fn->GetParErrors()[1];
	Double_t b = cal_fn->GetParameter(0); Double_t ub = cal_fn->GetParErrors()[0];
	TCanvas *cal_can = new TCanvas("cal_can","cal_can");
	
	cal_graph->Draw();
	cal_graph->GetYaxis()->SetRangeUser(-100, 1400);
	cal_graph->GetXaxis()->SetRangeUser(-1000, xMax);


	//Re-scale the histogram X-axis with the calculated calibration
	Float_t new_end = cal_fn->Eval(xMax);
	TH1F *hcal = new TH1F("hcal", "", hist->GetNbinsX(), 0, new_end); //create calibrated histo!!
	for(Int_t i=0;i<bins;i++) hcal->SetBinContent(i, hist->GetBinContent(i)); //fill cal histo with same content!
	TCanvas *cal_histo_can = new TCanvas("cal_histo_can","cal_histo_can");
	cal_histo_can->SetLogy();
	hcal->SetLineColor(kBlue);
	hcal->SetFillColorAlpha(kAzure+7, 0.4);
	hcal->GetXaxis()->SetTitle("Photon energy [keV]");
	hcal->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hcal->GetXaxis()->GetBinWidth(0)));
	hcal->GetYaxis()->SetTitleOffset(1.25);

	

	TString save_route = "preliminary_plots/calibration/";
	mkdir(save_route, 0777);

	hist->GetXaxis()->SetRangeUser(0, 25000);
	hist->GetXaxis()->SetNdivisions(508, kTRUE);
	c1->SaveAs(save_route + string(detector_to_calibrate) + "_peakFit.pdf");
	c1->Close();

	cal_can->SaveAs(save_route + string(detector_to_calibrate) + "_calcan.pdf");
	cal_can->Close();

	hcal->GetXaxis()->SetRangeUser(0, cal_fn->Eval(25000));
	hcal->Draw();
	cal_histo_can->SaveAs(save_route + string(detector_to_calibrate) + "_calSpectrum.pdf");


	std::ofstream data;
	if (detector_to_calibrate == "TAGGER") {data.open("calibration_data.csv", std::ofstream::out); data << "," << Form("%2.2f", xval[0]) << " keV," << Form("%2.2f", xval[1]) << " keV," << Form("%3.0f", xval[2]) << " keV," << "a+b*ch" << endl; data << "#DETECTOR," << "Resolution(%)," << "Resolution(%)," << "Resolution(%)," << "a," << "ua," << "b," << "ub" << endl;}
	else if (detector_to_calibrate != "TAGGER") data.open("calibration_data.csv", std::ofstream::app);
	data << detector_to_calibrate << "," << Form("%2.2f,", resolution[0]) << Form("%2.2f,", resolution[1]) << Form("%2.2f,", resolution[2]) << Form("%1.2f,", a) << Form("%1.2f,", ua) << Form("%1.5f,", b) << Form("%1.5f,", ub) << endl;
	data.close();

}

void detectorCalibration()
{
	TString routes[3] = {"TAGGER", "SCATTERER", "DETECTOR"};
	for (Int_t i = 0; i<3; i++) calibration(routes[i]);
	
}