#include "../../examples/data_analysis_examples/gethisto.C"



void detectorCalibration(const char* route = "/Users/javi/OneDrive/Documentos/Padova/1/PhysicsLaboratory/Gr28/TAC_calibration_7ns.root", Int_t npeaks = 2) {
	
	//Style settings
	gROOT->SetStyle("Plain");
	gStyle->SetOptFit(111);
	gStyle->SetPalette(1);

	Int_t bins = 1024;

	TH1D *hist = getHistoForChannelFromTree(route, 3, bins, 0, 25);
	TCanvas *c1 = new TCanvas("c1", "calibration_canvas");

	hist->Draw();
	hist->Rebin(4); //rebinin for a better visualization and peak finding

	TSpectrum *ss = new TSpectrum(npeaks); //creating TSpectrum object
	Int_t nfound = ss->Search(hist,10,"new",0.1); //searching peaks (TSpecrum class method)
	cout << "\n --------- Found " << nfound << " peaks" << endl << endl << endl;
	Double_t *xpeaks;
	xpeaks = (Double_t*) ss->GetPositionX(); //saving the position of the found peaks in an array 


	
}


