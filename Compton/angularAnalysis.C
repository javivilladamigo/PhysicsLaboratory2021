/// \author Javier Mariño Villadamigo

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};

Double_t cal_TAGGER(Int_t ch) {return (42.42 + ch*0.05906);}
Double_t cal_SCATTERER(Int_t ch) {return (42.66 + ch*0.05179);}
Double_t cal_DETECTOR(Int_t ch) {return (42.68 + ch*0.05240);}



TH1D *selectionOfEvents(TH1D *hist, Double_t xMin, Double_t xMax, Int_t CI, Double_t *Ei_min, Double_t *Ei_max) {
    hist->GetXaxis()->SetRangeUser(xMin, xMax);
    TSpectrum *ss = new TSpectrum(1); //creating TSpectrum object
    Int_t nfound = ss->Search(hist, 10, "goff", 0.005); //searching peaks (TSpectrum class method)
    Double_t *xpeaks;
	xpeaks = (Double_t*) ss->GetPositionX(); //saving the position of the found peaks in an array
    Double_t xp = 0.;
	Float_t mean = 0., sigma = 0.;
    TF1 *fun1, *fun2;
	for(Int_t i=0;i<nfound;i++)
	{
		xp = xpeaks[i];

		fun1 = new TF1("fun1name", "gaus", xp-120,xp+120);

		hist->Fit(fun1,"RQNP");

		fun2 = new TF1("fun2name","gaus(0)+pol1(3)", fun1->GetParameter(1)-(1.5*fun1->GetParameter(2)),fun1->GetParameter(1)+(1.5*fun1->GetParameter(2)));
		fun2->SetParameters(fun1->GetParameter(0), fun1->GetParameter(1), fun1->GetParameter(2));
		fun2->SetLineColor(kRed);

		hist->Fit(fun2,"RN+");

		mean  = fun2->GetParameter(1);
		sigma = fun2->GetParameter(2);
	}
    TH1D *h_filtered = (TH1D*)hist->Clone();
    *Ei_min = mean-CI*sigma;
    *Ei_max = mean+CI*sigma;

    cout << "Selection of TAGGER; from: " << *Ei_min << " keV to " << *Ei_max << " keV. " << endl;
    for (Int_t i=0; i<h_filtered->GetNbinsX(); i++)
    {
        if (i > hist->FindFixBin(*Ei_min) && i < hist->FindFixBin(*Ei_max))
        {
            h_filtered->SetBinContent(i, hist->GetBinContent(i));
        }
        else
        {
            h_filtered->SetBinContent(i, 0);
        }
    }
    return h_filtered;
}

void angularAnalysis(Int_t theta)
{   

    //Style settings
    Int_t MyPalette[100];
    Double_t Red[]    = {0., 0.0, 1.0, 1.0, 1.0};
    Double_t Green[]  = {0., 0.0, 0.0, 1.0, 1.0};
    Double_t Blue[]   = {0., 1.0, 0.0, 0.0, 1.0};
    Double_t Length[] = {0., .25, .50, .75, 1.0};
    Int_t FI = TColor::CreateGradientColorTable(5, Length, Red, Green, Blue, 100);
    for (int i=0;i<100;i++) MyPalette[i] = FI+i;
    gStyle->SetPalette(1);

    gStyle->SetTitleFont(62, "XYZ"); // title of axes
    gStyle->SetTitleSize(0.035, "XYZ"); // size of axes titles

    gStyle->SetTitleFont(62, "t"); // size of axes titles
    gStyle->SetTitleFontSize(0.08); // size of axes titles

    gStyle->SetLabelFont(62, "XYZ"); // labels of axes
    gStyle->SetLabelSize(0.045, "XYZ"); // size of axes labels

    gStyle->SetTitleOffset(1.4, "X"); // offset of x-axis
    gStyle->SetTitleOffset(2.5, "Y"); // offset of y-axis
    gStyle->SetTitleOffset(1.8, "Z"); // offset of z-axis
    gStyle->SetTextFont(62); // general texts

    gStyle->SetOptStat(0); // clearing TPaveStats


    
    // declaring variables
    Int_t topCh = int(TMath::Power(2, 12));
	Int_t bins = topCh;
	Int_t xMin = 0.;
	Int_t xMax = int(TMath::Power(2, 15));
    Double_t deltaTheta = (TMath::ATan(4. / 3.757) + TMath::ATan(8. / 26.)) / TMath::Pi() * 180; //TMath::ATan(4. / 26.) / TMath::Pi() * 180; // geometry definition for 260 mm of distance between SCATTERER and DETECTOR
    Double_t Ei, Ee, Ef; // energies of initial photon, electron and final phton
    Double_t Ei_min, Ei_max; // selection of energies in TAGGER


    // declaring calibrated histograms
    TH1D *hist_TAGGER_cal = new TH1D("tagger_cal", "Tagger spectrum", bins, 0, cal_TAGGER(xMax));
    TH1D *hist_SCATTERER_cal = new TH1D("scatterer_cal", "Scatterer spectrum", bins, 0, cal_SCATTERER(xMax));
    TH1D *hist_DETECTOR_cal = new TH1D("detector_cal", "Detector spectrum", bins, 0, cal_DETECTOR(xMax));
    
    
    // reading the file
    TString fileName = "data/day2/angular/D2D3=260mm/" + TString::Itoa(theta, 10) + "deg_fixed.root";
    const char *route = fileName;
    slimport_data_t indata_TAGGER, indata_SCATTERER, indata_DETECTOR;
    TFile *infile = new TFile(route);
    TTree *intree= (TTree*)infile->Get("acq_tree_0");
    TBranch *TAGGER_inbranch = intree->GetBranch("acq_ch0");
    TAGGER_inbranch->SetAddress(&indata_TAGGER.timetag);
    TBranch *SCATTERER_inbranch = intree->GetBranch("acq_ch1");
    SCATTERER_inbranch->SetAddress(&indata_SCATTERER.timetag);
    TBranch *DETECTOR_inbranch = intree->GetBranch("acq_ch2");
    DETECTOR_inbranch->SetAddress(&indata_DETECTOR.timetag);

    for (int i=0; i<DETECTOR_inbranch->GetEntries(); i++)
    {
        TAGGER_inbranch->GetEntry(i);
        hist_TAGGER_cal->Fill(cal_TAGGER(indata_TAGGER.qlong)-32.42); // here I am just shifting a little bit the spectrum to recalibrate (since the ordinate of the calibration (~40+-30) does not seem very accurate)
        SCATTERER_inbranch->GetEntry(i);
        hist_SCATTERER_cal->Fill(cal_SCATTERER(indata_SCATTERER.qlong)-32.66);
        DETECTOR_inbranch->GetEntry(i);
        hist_DETECTOR_cal->Fill(cal_DETECTOR(indata_DETECTOR.qlong)-32.68);
    }
    
    // declaring filtered histograms
    TH1D *hist_TAGGER_filtered = selectionOfEvents(hist_TAGGER_cal, 500, 570, 3, &Ei_min, &Ei_max); // filtering the events in TAGGER within 3 sigma from the 511 peak
    TH1D *hist_SCATTERER_filtered = (TH1D*)hist_SCATTERER_cal->Clone();
    TH1D *hist_DETECTOR_filtered = new TH1D("hist_DETECTOR_filtered", "Detector spectrum", bins, 0, cal_SCATTERER(xMax));
    
    // kinematically filtering the events in the SCATTERER from those in TAGGER
    Double_t Ee_min = Ei_min*(1 - 1 / (2 - TMath::Cos((theta - deltaTheta) / 180 * TMath::Pi()))); // electron energy corresponding to the minimum Ei (and theta-deltaTheta)
    Double_t Ee_max = Ei_max*(1 - 1 / (2 - TMath::Cos((theta + deltaTheta) / 180 * TMath::Pi()))); // electron energy corresponding to the maximum Ei (and theta+deltaTheta)
    if (theta -  deltaTheta <= 0) Ee_min = 0; // fixing the distribution of the cosine
    // (because when the angular aperture is centered in 0, theta+-deltaTheta no longer give the min and max value of cos respectively;
    // both are minima and cos(0) is the maximum)

    cout << "Selection of SCATTERER; from: " << Ee_min << " keV to " << Ee_max << " keV." << endl;

    for (Int_t i = 0; i<hist_SCATTERER_filtered->GetNbinsX(); i++)
    {
        if (i > hist_SCATTERER_cal->FindFixBin(Ee_min) && i < hist_SCATTERER_cal->FindFixBin(Ee_max))
        {
            hist_SCATTERER_filtered->SetBinContent(i, hist_SCATTERER_cal->GetBinContent(i));
        }
        else
        {
            hist_SCATTERER_filtered->SetBinContent(i, 0);
        }
    }
    
    TH2D *hist_correlation = new TH2D("hist_correlation", "Detector vs. Scatterer correlation", hist_SCATTERER_filtered->GetNbinsX(), 0, cal_DETECTOR(xMax), hist_DETECTOR_filtered->GetNbinsX(), 0, cal_DETECTOR(xMax));
    // coincidence between tagger and scatterer selected spectra
	for (int i=0; i<SCATTERER_inbranch->GetEntries(); i++)
    {
        TAGGER_inbranch->GetEntry(i);
        Ei = cal_TAGGER(indata_TAGGER.qlong);
        SCATTERER_inbranch->GetEntry(i);
        Ee = cal_SCATTERER(indata_SCATTERER.qlong);
        DETECTOR_inbranch->GetEntry(i);
        Ef = cal_DETECTOR(indata_DETECTOR.qlong);

        if ((Ei<Ei_max && Ei>Ei_min) && (Ee<Ee_max && Ee>Ee_min) && (indata_TAGGER.timetag == indata_SCATTERER.timetag) && (indata_SCATTERER.timetag == indata_DETECTOR.timetag))
        {
            hist_DETECTOR_filtered->Fill(Ef);
            hist_correlation->Fill(Ee, Ef);
        }
        
        
    }

    // format and plotting
    TCanvas *c1 = new TCanvas ("c1", "Events unfiltered");
    c1->SetWindowSize(1920, 582);
    c1->Divide(3, 1);

    c1->cd(1);
    gPad->SetMargin(0.15881, 0.00358308, 0.100197, 0.0991254); // left right bottom top
    hist_TAGGER_cal->Rebin(32);
    hist_TAGGER_cal->GetXaxis()->SetNdivisions(505, kTRUE);
    hist_TAGGER_cal->SetLineColor(kBlue);
	hist_TAGGER_cal->SetFillColorAlpha(kAzure+7, 0.4);
    hist_TAGGER_cal->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist_TAGGER_cal->GetXaxis()->GetBinWidth(0)));
    hist_TAGGER_cal->GetXaxis()->SetTitle("Energy [keV]");
    hist_TAGGER_cal->Draw();

    c1->cd(2);
    gPad->SetMargin(0.15881, 0.00358308, 0.100197, 0.0991254); // left right bottom top
    hist_SCATTERER_cal->Rebin(32);
    hist_SCATTERER_cal->GetXaxis()->SetNdivisions(505, kTRUE);
    hist_SCATTERER_cal->SetLineColor(kBlue);
	hist_SCATTERER_cal->SetFillColorAlpha(kAzure+7, 0.4);
    hist_SCATTERER_cal->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist_SCATTERER_cal->GetXaxis()->GetBinWidth(0)));
    hist_SCATTERER_cal->GetXaxis()->SetTitle("Energy [keV]");
    hist_SCATTERER_cal->Draw();

    c1->cd(3);
    gPad->SetMargin(0.15881, 0.00358308, 0.100197, 0.0991254); // left right bottom top
    hist_DETECTOR_cal->Rebin(32);
    hist_DETECTOR_cal->GetXaxis()->SetNdivisions(505, kTRUE);
    hist_DETECTOR_cal->SetLineColor(kBlue);
	hist_DETECTOR_cal->SetFillColorAlpha(kAzure+7, 0.4);
    hist_DETECTOR_cal->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist_DETECTOR_cal->GetXaxis()->GetBinWidth(0)));
    hist_DETECTOR_cal->GetXaxis()->SetTitle("Energy [keV]");
    hist_DETECTOR_cal->Draw();

    TCanvas *c2 = new TCanvas("c2", "Events filtered");
    c2->SetWindowSize(1920, 582);
    c2->Divide(3, 1);

    c2->cd(1);
    gPad->SetMargin(0.15881, 0.00358308, 0.100197, 0.0991254); // left right bottom top
    hist_TAGGER_filtered->Rebin(32);
    hist_TAGGER_filtered->GetXaxis()->SetNdivisions(505, kTRUE);
    hist_TAGGER_filtered->SetLineColor(kBlue);
	hist_TAGGER_filtered->SetFillColorAlpha(kAzure+7, 0.4);
    hist_TAGGER_filtered->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist_TAGGER_filtered->GetXaxis()->GetBinWidth(0)));
    hist_TAGGER_filtered->GetXaxis()->SetTitle("Energy [keV]");
    hist_TAGGER_filtered->Draw();

    c2->cd(2);
    gPad->SetMargin(0.15881, 0.00358308, 0.100197, 0.0991254); // left right bottom top
    hist_SCATTERER_filtered->Rebin(32);
    hist_SCATTERER_filtered->GetXaxis()->SetNdivisions(505, kTRUE);
    hist_SCATTERER_filtered->SetLineColor(kBlue);
	hist_SCATTERER_filtered->SetFillColorAlpha(kAzure+7, 0.4);
    hist_SCATTERER_filtered->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist_SCATTERER_filtered->GetXaxis()->GetBinWidth(0)));
    hist_SCATTERER_filtered->GetXaxis()->SetTitle("Energy [keV]");
    hist_SCATTERER_filtered->Draw();

    c2->cd(3);
    gPad->SetMargin(0.15881, 0.00358308, 0.100197, 0.0991254); // left right bottom top
    hist_DETECTOR_filtered->Rebin(32);
    hist_DETECTOR_filtered->GetXaxis()->SetNdivisions(505, kTRUE);
    hist_DETECTOR_filtered->SetLineColor(kBlue);
	hist_DETECTOR_filtered->SetFillColorAlpha(kAzure+7, 0.4);
    hist_DETECTOR_filtered->GetYaxis()->SetTitle(Form("counts / %1.2f keV", hist_DETECTOR_filtered->GetXaxis()->GetBinWidth(0)));
    hist_DETECTOR_filtered->GetXaxis()->SetTitle("Energy [keV]");
    hist_DETECTOR_filtered->Draw();


    TCanvas *c3 = new TCanvas("c3", "Correlated Events");
    c3->SetWindowSize(1080, 1080);
    c3->cd();
    gPad->SetMargin(0.15881, 0.05, 0.100197, 0.0991254); // left right bottom top
    hist_correlation->GetYaxis()->SetTitle("Detector energy [keV]");
    hist_correlation->GetXaxis()->SetTitle("Scatterer energy [keV]");
    hist_correlation->GetZaxis()->SetTitle(Form("counts / %1.2f keV", hist_DETECTOR_filtered->GetXaxis()->GetBinWidth(0)));
    hist_correlation->SetMinimum(14);
    hist_correlation->Rebin2D(16, 16);
    hist_correlation->Draw("lego2");
    hist_correlation->GetXaxis()->SetRangeUser(0, 520);
    hist_correlation->GetYaxis()->SetRangeUser(0, 520);
    cout << hist_correlation->GetEntries() << endl;


    TString save_route = "preliminary_plots/angularAnalysis/";
	mkdir(save_route, 0777);
	c1->SaveAs(save_route + theta + ".pdf");
}