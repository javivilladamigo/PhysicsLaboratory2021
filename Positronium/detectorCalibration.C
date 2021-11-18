#include "../../examples/data_analysis_examples/gethisto.C"


// in order to do the calibration of several detectors
// modify the corresponding channel (ch1) in the route variable (i.e. const char* route = "data/day1/ch1Calibration.root")
// and always keep channel 1 in the call "TH1D *hist = getHistoForChannelFromTree(route, 1, bins, xMin, xMax)"

void detectorCalibration(const char* route = "data/day1/ch4Calibration.root", Int_t npeaks = 2) {
	
	//Style settings
	gROOT->SetStyle("Default");
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	gStyle->SetPalette(0);

	Int_t bins = 4*4096;
	Int_t xMin = 0;
	Int_t xMax = 4*4096;
	
	TCanvas *c1 = new TCanvas("c1", "calibration_canvas");
	TH1D *hist = getHistoForChannelFromTree(route, 1, bins, xMin, xMax);
	hist->Rebin(4); //rebinin for a better visualization and peak finding

	hist->SetTitle("");
	hist->GetXaxis()->SetTitle("channel"); hist->GetYaxis()->SetTitle("counts");
	hist->GetYaxis()->SetTitleOffset(1.5);

	hist->Draw();
	

	TSpectrum *ss = new TSpectrum(npeaks); //creating TSpectrum object
	Int_t nfound = ss->Search(hist,10,"goff", 0.1); //searching peaks (TSpecrum class method)
	cout << "\n --------- Found " << nfound << " peaks" << endl << endl << endl;
	Double_t *xpeaks;
	xpeaks = (Double_t*) ss->GetPositionX(); //saving the position of the found peaks in an array

	TH1D *hfinal = (TH1D*)hist->Clone();

	Double_t xp = 0.;
	Float_t mean = 0., sigma = 0.;
	Float_t offset = 0, slope = 0;

	Float_t resolution[3]={0.,0.,0.};

	//Sort xpeaks array
	int n = sizeof(*xpeaks)/sizeof(xpeaks[0]);
	std::sort(xpeaks,xpeaks+nfound);



	TF1 *fun2, *fun1, *fit;
	for(Int_t i=0;i<nfound;i++)
	{
		xp = xpeaks[i];
		cout << " ************ Peak at: " << xp << endl;

		fun1 = new TF1("fun1name","gaus",xp-120,xp+120);

		hfinal->Fit(fun1,"RQN");

		fun2 = new TF1("fun2name","gaus(0)+pol1(3)", fun1->GetParameter(1)-(2.5*fun1->GetParameter(2)),fun1->GetParameter(1)+(2.5*fun1->GetParameter(2)));
		fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
		fun2->SetLineColor(kRed);
		hfinal->Fit(fun2,"R+");

		mean  = fun2->GetParameter(1);
		sigma = fun2->GetParameter(2);


		fit = hfinal->GetFunction("fun2name");
		Double_t chi2 = fit->GetChisquare();
		Double_t Ndf = fit->GetNDF();

		resolution[i]=((2.35*sigma)/mean)*100;
		cout << "\nRESULTS:  Mean at: " << mean << " with sigma " << sigma << " , FWHM of " << 2.35*sigma << " , and resolution of " << resolution[i] << " %" << endl << endl;
		cout << "\nThe fit shows a Chi2/Ndf = " << chi2/Ndf << endl;

		}
	

	//Calibration
	Float_t xval[2];
	if(npeaks==1) xval[0] = 511; else if(npeaks==2) {xval[0] = 511; xval[1] = 1275;};
	Float_t peaks[npeaks];
	for(Int_t i=0;i<npeaks;i++) peaks[i] = xpeaks[i];
	TGraph *cal_graph = new TGraph(npeaks,peaks,xval); //calibration: channel vs energy of the found peaks
	cal_graph->SetNameTitle("cal_graph","");
	cal_graph->SetMarkerStyle(21);
	cal_graph->SetMarkerColor(4);
	TF1 *cal_fn = new TF1("cal_fn","[0]*x+[1]"); //defining a first order polynomial function and fit the data
	cal_graph->Fit(cal_fn);
	TCanvas *cal_can = new TCanvas("cal_can","cal_can");
	cal_graph->Draw();

	cout << endl << "The calibration has offset: " << cal_fn->GetParameter(1) << " and slope: " << cal_fn->GetParameter(0) << endl;
	//Re-escale the histogram X-axis with the calculated calibration
	Float_t new_end = cal_fn->Eval(bins);
	TH1F *hcal = new TH1F("hcal","",hist->GetNbinsX(),0,new_end); //create calibrated histo!!
	for(Int_t i=0;i<bins;i++) hcal->SetBinContent(i+1,hist->GetBinContent(i+1)); //fill cal histo with same content!
	TCanvas *cal_histo_can = new TCanvas("cal_histo_can","cal_histo_can");
	hcal->Draw();
	hcal->GetXaxis()->SetTitle("Photon energy [keV]");
	hcal->GetXaxis()->CenterTitle();

}


