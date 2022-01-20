#include "../../examples/data_analysis_examples/gethisto.C"

/// \author Javier Mariño Villadamigo


void calibrationCo(TString detector_to_calibrate, Int_t npeaks = 2, Float_t *meanCo1 = 0, Float_t *meanCo2 = 0, Float_t *sigmaCo1 = 0, Float_t *sigmaCo2 = 0) {

	Int_t topCh = int(TMath::Power(2, 15));
	Int_t bins = topCh;
	Int_t xMin;
	Int_t xMax;
	if (detector_to_calibrate == "NaI") {xMin = 9000; xMax = 12000;}
	if (detector_to_calibrate == "HPGe") {xMin = 6500; xMax = 8000;}

	
	short chan;
	const char *route;
	if (detector_to_calibrate == "NaI")
	{
		route = "/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day1/10min_Co.root";
		chan = 0;
	}
	else if (detector_to_calibrate == "HPGe")
	{
		route = "/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/HPGe_Co10min.root";
		chan = 1;
	}

	TH1D *hist = getHistoForChannelFromTree(route, chan, bins, 0, xMax);
	
	hist->Rebin(16); //rebinin for a better visualization and peak finding
	hist->SetTitle("");
	hist->GetXaxis()->SetTitle("channel "); hist->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist->GetXaxis()->GetBinWidth(0)));
	hist->GetXaxis()->SetRangeUser(xMin, xMax);
	hist->GetYaxis()->SetTitleOffset(1.25);
	hist->SetLineColor(kBlue);
	hist->SetFillColorAlpha(kAzure+7, 0.4);

	hist->Draw();

	TSpectrum *ss = new TSpectrum(npeaks); //creating TSpectrum object
	Int_t nfound = ss->Search(hist, 10, "goff", 0.1); //searching peaks (TSpecrum class method)
	//cout << "\n --------- Found " << nfound << " peaks" << endl << endl << endl;
	Double_t *xpeaks;
	xpeaks = (Double_t*) ss->GetPositionX(); //saving the position of the found peaks in an array
	Double_t xp = 0.;
	Float_t mean[2] = {};
	Float_t sigma[2] = {};
	Float_t resolution[2]={0., 0.};
	Float_t chi2 = 0.;
	Float_t Ndf = 0.;

	//Sort xpeaks array
	int n = sizeof(*xpeaks)/sizeof(xpeaks[0]);
	std::sort(xpeaks,xpeaks+nfound);


	TF1 *fun1, *fun2, *fit;
	for(Int_t i=0; i<nfound; i++)
	{
		xp = xpeaks[i];
		//cout << " ************ Peak at: " << xp << endl;

		fun1 = new TF1("fun1name", "gaus", xp-120,xp+120);

		hist->Fit(fun1,"RQN");

		fun2 = new TF1("fun2name","gaus(0)+pol1(3)", fun1->GetParameter(1)-(2.0*fun1->GetParameter(2)),fun1->GetParameter(1)+(2.0*fun1->GetParameter(2)));
		fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
		fun2->SetLineColor(kRed);

		hist->Fit(fun2,"RQ+");
		
		mean[i] = fun2->GetParameter(1);
		sigma[i] = fun2->GetParameter(2);
	}
	*meanCo1 = mean[0];
	*meanCo2 = mean[1];
	*sigmaCo1 = sigma[0];
	*sigmaCo2 = sigma[1];
}
void calibrationAm(TString detector_to_calibrate, Int_t npeaks = 1, Float_t *meanAm = 0, Float_t *sigmaAm = 0) {

	Int_t topCh = int(TMath::Power(2, 15));
	Int_t bins = topCh;
	Int_t xMin = 100.;
	Int_t xMax;
	if (detector_to_calibrate == "NaI") xMax = 15000;
	if (detector_to_calibrate == "HPGe") xMax = 3000;
	
	short chan;
	const char *route;
	if (detector_to_calibrate == "NaI")
	{
		route = "/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day1/10min_Am.root";
		chan = 0;
	}
	else if (detector_to_calibrate == "HPGe")
	{
		route = "/Users/javi/Documents/Padova/1/PhysicsLaboratory/PhysicsLaboratory2021/Gamma/data/full_data/day2/HPGe_Am10min.root";
		chan = 1;
	}

	TH1D *hist = getHistoForChannelFromTree(route, chan, bins, 0, xMax);
	
	hist->Rebin(32); //rebinin for a better visualization and peak finding
	hist->SetTitle("");
	hist->GetXaxis()->SetTitle("channel "); hist->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist->GetXaxis()->GetBinWidth(0)));
	hist->GetXaxis()->SetRangeUser(xMin, xMax);
	hist->GetYaxis()->SetTitleOffset(1.25);
	hist->SetLineColor(kBlue);
	hist->SetFillColorAlpha(kAzure+7, 0.4);



	TSpectrum *ss = new TSpectrum(npeaks); //creating TSpectrum object
	Int_t nfound = ss->Search(hist, 10, "goff", 0.1); //searching peaks (TSpecrum class method)
	//cout << "\n --------- Found " << nfound << " peaks" << endl << endl << endl;
	Double_t *xpeaks;
	xpeaks = (Double_t*) ss->GetPositionX(); //saving the position of the found peaks in an array

	Double_t xp = 0.;

	//Sort xpeaks array
	int n = sizeof(*xpeaks)/sizeof(xpeaks[0]);
	std::sort(xpeaks,xpeaks+nfound);

	TF1 *fun1, *fun2, *fit;
	for(Int_t i=0; i<nfound; i++)
	{
		xp = xpeaks[i];
		//cout << " ************ Peak at: " << xp << endl;

		fun1 = new TF1("fun1name", "gaus", xp-120,xp+120);

		hist->Fit(fun1,"RQN");

		fun2 = new TF1("fun2name","gaus(0)+pol1(3)", fun1->GetParameter(1)-(2.0*fun1->GetParameter(2)),fun1->GetParameter(1)+(2.0*fun1->GetParameter(2)));
		fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
		fun2->SetLineColor(kRed);

		hist->Fit(fun2,"RQN+");

		*meanAm = fun2->GetParameter(1);
		*sigmaAm = fun2->GetParameter(2);
	}
}



void detectorCalibration(TString detector_to_calibrate)
{
	Int_t xMin, xMax, yMin, yMax;
	
	if (detector_to_calibrate == "NaI") {xMin = -200; xMax = 12000; yMin = -50; yMax = 1600;}
	if (detector_to_calibrate == "HPGe") {xMin = -200; xMax = 8000; yMin = -50; yMax = 1400;}

	//Style settings
	gROOT->SetStyle("Default");
	gStyle->SetCanvasColor(0);
	gStyle->SetStatBorderSize(1);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(0);

	Int_t npeaks = 3;
	Float_t meanCo1, meanCo2, sigmaCo1, sigmaCo2;
	Float_t meanAm, sigmaAm;

	calibrationCo(detector_to_calibrate, 2, &meanCo1, &meanCo2, &sigmaCo1, &sigmaCo2);
	calibrationAm(detector_to_calibrate, 1, &meanAm, &sigmaAm);
	cout << "Peaks for Co found in: " << meanCo1 << ", " << meanCo2 << " ch." << endl;
	cout << "Peaks for Am found in: " << meanAm << " ch." << endl;

	Float_t xpeaks[3] = {meanAm, meanCo1, meanCo2};
	Float_t sigma_peaks[3] = {sigmaAm, sigmaCo1, sigmaCo2};
	Float_t xval[3] = {59.54092, 1173.228, 1332.514};
	Float_t chi2 = 0;

	//Calibration

	TCanvas *cal_can = new TCanvas("cal_can","cal_can");
    TH2F *h = new TH2F("h", ";channel;Photon energy [keV]",  1000, xMin, xMax, 1000, yMin, yMax);
	cal_can->SetLeftMargin(0.11749);
	cal_can->SetRightMargin(0.0544413);
	cal_can->SetBottomMargin(0.101053);
	cal_can->SetTopMargin(0.0336842);
	h->GetYaxis()->SetTitleOffset(1.4);
    h->Draw();

	
	TGraphErrors *cal_graph = new TGraphErrors(npeaks, xpeaks, xval, 0, sigma_peaks); //calibration: channel vs energy of the found peaks
	cal_graph->SetNameTitle("cal_graph","");
	cal_graph->SetMarkerStyle(20);
	cal_graph->SetMarkerColor(1);
	

	TF1 *cal_fn = new TF1("cal_fn","[0]*x+[1]"); //defining a first order polynomial function and fit the data
	cal_fn->SetLineColor(kRed);


	cal_graph->Draw("P");
	cal_graph->Fit(cal_fn);
	chi2 = cal_fn->GetChisquare();

	Double_t a = cal_fn->GetParameter(1); Double_t ua = cal_fn->GetParErrors()[1];
	Double_t b = cal_fn->GetParameter(0); Double_t ub = cal_fn->GetParErrors()[0];

	TLatex text;
	text.SetTextSize(0.05);
	if (detector_to_calibrate == "NaI") text.DrawLatex(1000., 1200., Form("E [keV] =  %1.0f(%1.0f) + %1.2f(%1.0f) #upoint E_{ch}", a, ua, b, ub*100));
    if (detector_to_calibrate == "HPGe") text.DrawLatex(700., 1050., Form("E [keV] =  %1.0f(%1.0f) + %1.3f(%1.0f) #upoint E_{ch}", a, ua, b, ub*1000));

	gSystem->Exec("mkdir preliminary_plots/");
	gSystem->Exec("mkdir preliminary_plots/calibration");
	TString save_route = "preliminary_plots/calibration/";
	cal_can->SaveAs(save_route + string(detector_to_calibrate) + "_calcan.pdf");	

	std::ofstream data;
	if (detector_to_calibrate == "NaI") {data.open("calibration_data.csv", std::ofstream::out);  data << "a+b*ch" << endl; data << "#DETECTOR," << "a," << "ua," << "b," << "ub," << "chi2" << endl;}
	else if (detector_to_calibrate == "HPGe") data.open("calibration_data.csv", std::ofstream::app);
	data << detector_to_calibrate << "," 	<< Form("%1.2f,", a) << Form("%1.2f,", ua) << Form("%1.5f,", b) << Form("%1.5f,", ub) << Form("%1.10f", chi2) << endl;
	data.close();

}