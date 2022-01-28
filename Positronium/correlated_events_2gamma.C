#include <stdexcept>

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};

Double_t cal1(Int_t ch) {return (-18.894 + ch*0.119077);}
Double_t cal2(Int_t ch) {return (-28.3349 + ch*0.115687);}
Double_t TAC_cal(Int_t ch) {return 0.01179*ch-140.2 + 30;} // 0.01179*ch-140.2 real one

void correlated_events_2gamma(const char *name_file = "data/day3/2gamma_filteredEvents.root", int bins = 16384, double xMin = 128, double xMax = 16384) {
    
    gROOT->SetStyle("Default");
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
	gStyle->SetStatBorderSize(1);
    gStyle->SetTitleFontSize(.09);

    // variables
    Double_t angular_aperture = TMath::ATan(5./15.);

    Double_t Tmin = 0; Double_t Tmax = 200;
    Double_t Emin = 406; Double_t Emax = 615;

    slimport_data_t indataA, indataB, indataD;
    TFile *infile = new TFile(name_file);
    TTree *intree = (TTree*)infile->Get("acq_tree_0");

    TBranch *inbranchA = intree->GetBranch("acq_ch0");
    TBranch *inbranchB = intree->GetBranch("acq_ch1");
	TBranch *inbranchD = intree->GetBranch("acq_ch3");

    inbranchA->SetAddress(&indataA.timetag);
    inbranchB->SetAddress(&indataB.timetag);
    inbranchD->SetAddress(&indataD.timetag);

    TCanvas *c1 = new TCanvas("c1", "Detector #1");
    TH1D *histA = new TH1D("histA", "Detector #1", bins, cal1(xMin), cal1(xMax));


    TCanvas *c2 = new TCanvas("c2", "Detector #2");
    TH1D *histB = new TH1D("histB", "Detector #2", bins, cal2(xMin), cal2(xMax));

    TCanvas *c4 = new TCanvas("c4", "TAC");
    c4->SetLogy();
    TH1D *histD = new TH1D("histD", "TAC", bins, Tmin, Tmax);
    TH1D *histD_clone = (TH1D*)histD->Clone();

    TCanvas *c5 = new TCanvas("c5", "sum_spectra");
    TH1D *hist_sum = new TH1D("hist_sum", "", bins, 0, 4000);

    bool goodToGo = 0;
    // by now we should check we have the same number of entries in every branch
    if (inbranchA->GetEntries() == inbranchB->GetEntries() && inbranchB->GetEntries() == inbranchD->GetEntries())
    {
        cout << "All channels have the same number of entries." << endl;
        goodToGo = 1;

    }
    else {
        throw std::runtime_error( "Channels do not have the same number of entries." );
    }

    if (goodToGo == 1) {
        for (Int_t i = 0; i<inbranchA->GetEntries(); i++) // parsing the file
        {
            inbranchA->GetEntry(i);
            double E0 = cal1(indataA.qlong); // pass to energy
            

            inbranchB->GetEntry(i);
            double E1 = cal2(indataB.qlong); // pass to energy
            

            inbranchD->GetEntry(i);
            double T = TAC_cal(indataD.qlong);
            
            if (indataA.timetag == indataB.timetag && indataB.timetag == indataD.timetag)
            {
                if (E0 > Emin && E0 < Emax && Tmin<T && T<Tmax) histA->Fill(E0);
                if (E1 > Emin && E1 < Emax && Tmin<T && T<Tmax) histB->Fill(E1);
                
                if (((E0 > Emin && E0 < Emax) || (E1 > Emin && E1 < Emax)) && Tmin<T && T<Tmax)
                {
                    histD->Fill(T); // get the TAC entry for that possible coincidences
                    hist_sum->Fill(E0+E1); // fill the sum histogram if any of the events is in the desired energy range
                    if ((E0+E1)<1104.73 && (E0+E1) > 956.563 && T<260)
                    {
                        histD_clone->Fill(T); // filling the timetags of the events that are later going to be counted as coincidences
                    }
                }
            }      
        }
        c1->cd();
        histA->Rebin(32);
        histA->GetYaxis()->SetTitle("counts"); histA->GetXaxis()->SetTitle("Photon energy [keV]");
        histA->Draw();

        c2->cd();
        histB->Rebin(32);
        histB->GetYaxis()->SetTitle("counts"); histB->GetXaxis()->SetTitle("Photon energy [keV]");
        histB->Draw();

        c4->cd();
        histD->Rebin(128);
        histD->GetYaxis()->SetTitle("counts"); histD->GetXaxis()->SetTitle("Delay [ns]");
        histD_clone->Rebin(128);
        histD_clone->SetFillColorAlpha(kBlue, 0.5);

        // filtering the TAC
        for (Int_t i = 0; i<histD->GetNbinsX(); i++)
        {
            if (i < histD->FindFixBin(0))
            {
                histD_clone->SetBinContent(i, 0);
            }
        }

        histD->Draw();
        histD_clone->Draw("SAME");
        TF1 *expo_free = new TF1("expo_freename", "[0]+[1]*TMath::Exp([2]*x+[3])", 22, 100);
        expo_free->SetParameters(2, 80, -5, 0);
        expo_free->SetLineColor(kRed);
        expo_free->SetLineWidth(8.);
        histD_clone->Fit("expo_freename", "R+");
        Double_t lifetime1 = -1/expo_free->GetParameter(2); // in ns
        Double_t error_lifetime1 = expo_free->GetParErrors()[2]/(lifetime1*lifetime1);
        //Double_t lifetime_2 = -expo_free->GetParameter(5); // in ns
        cout << "First lifetime comp: " << lifetime1 << endl;
        //cout << "Second lifetime comp: " << lifetime_2 << endl;
        Double_t chi2_expo = expo_free->GetChisquare();
		Double_t Ndf_expo = expo_free->GetNDF();
        cout << "chi2/NDf : " << chi2_expo/Ndf_expo << endl;


        c5->cd();

        hist_sum->GetYaxis()->SetTitleOffset(1.3);
        hist_sum->GetXaxis()->SetTitleSize(0.045);
        hist_sum->GetYaxis()->SetTitleSize(0.045);

        hist_sum->Rebin(16);
        hist_sum->GetYaxis()->SetTitle("counts"); hist_sum->GetXaxis()->SetTitle("Photon energy sum [keV]");
        hist_sum->SetLineColor(kBlue);
        hist_sum->SetLineWidth(2);
        //hist_sum->GetYaxis()->SetRangeUser(0, 40000); 
        hist_sum->GetXaxis()->SetRangeUser(450, 1350);
        hist_sum->Draw();



        Double_t limits[2] = {923.58, 1132.81};
        TF1 *gauss_plus_bg = new TF1("gauss_plus_bg","gaus(0)+pol1(3)", limits[0], limits[1]);
		gauss_plus_bg->SetParameters(3e5, 1030, 24);
        gauss_plus_bg->SetParLimits(1, 950, 1100);
		gauss_plus_bg->SetLineColor(kRed);
        gauss_plus_bg->SetLineWidth(8.);
        hist_sum->Fit("gauss_plus_bg", "R+");

        Double_t mean = gauss_plus_bg->GetParameter(1);
        Double_t sigma = gauss_plus_bg->GetParameter(2);
        TF1 *fit = hist_sum->GetFunction("gauss_plus_bg");
        Double_t chi2 = fit->GetChisquare();
		Double_t Ndf = fit->GetNDF();

        Double_t nb_of_events = hist_sum->Integral(hist_sum->FindFixBin(mean-3*sigma), hist_sum->FindFixBin(mean+3*sigma));
        TH1D *colored_hist = (TH1D*)hist_sum->Clone();
        colored_hist->GetXaxis()->SetRangeUser(mean-3*sigma, mean+3*sigma);
        colored_hist->SetLineColor(kRed);
        colored_hist->SetLineWidth(4);
        colored_hist->SetFillColorAlpha(kRed, 0.4);
        colored_hist->Draw("SAME");
        
        TLatex *text = new TLatex;
        text->DrawLatex(1100, 3e3, Form("#color[2]{#chi^{2}/Ndf = %1.2f}", chi2/Ndf));
        text->SetTextSize(0.06);


        TLatex *text2 = new TLatex;
        text2->DrawLatex(1100, 3e3, Form("#color[2]{Integral (#pm 3#sigma) = %.0f events}", nb_of_events));
        text2->SetTextSize(0.06);


    }
}
