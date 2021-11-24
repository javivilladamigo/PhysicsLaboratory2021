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
Double_t cal3(Int_t ch) {return (-25.3076 + ch*0.110853);}
Double_t TAC_cal(Int_t ch) {return ((ch-12040)/84.74);} // real one 

void correlated_events(const char *name_file = "data/day3/3gamma_filteredEvents.root", int bins = 16384, double xMin = 128, double xMax = 16384) {
    
    gROOT->SetStyle("Default");

    // variables
    Double_t angular_aperture = TMath::ATan(5./15.);
    
    Double_t Tmin = -30; Double_t Tmax = 240;
    Double_t Emin = 236; Double_t Emax = 445;

    slimport_data_t indataA, indataB, indataC, indataD;
    TFile *infile = new TFile(name_file);
    TTree *intree = (TTree*)infile->Get("acq_tree_0");

    TBranch *inbranchA = intree->GetBranch("acq_ch0");
    TBranch *inbranchB = intree->GetBranch("acq_ch1");
	TBranch *inbranchC = intree->GetBranch("acq_ch2");
	TBranch *inbranchD = intree->GetBranch("acq_ch3");
    

    inbranchA->SetAddress(&indataA.timetag);
    inbranchB->SetAddress(&indataB.timetag);
    inbranchC->SetAddress(&indataC.timetag);
    inbranchD->SetAddress(&indataD.timetag);


    TCanvas *c1 = new TCanvas("c1", "Detector #1");
    TH1D *histA = new TH1D("histA", "Detector #1", bins, cal1(xMin), cal1(xMax));


    TCanvas *c2 = new TCanvas("c2", "Detector #2");
    TH1D *histB = new TH1D("histB", "Detector #2", bins, cal2(xMin), cal2(xMax));

    TCanvas *c3 = new TCanvas("c3", "Detector #3");
    TH1D *histC = new TH1D("histC", "Detector #3", bins, cal3(xMin), cal3(xMax));

    TCanvas *c4 = new TCanvas("c4", "TAC");
    TH1D *histD = new TH1D("histD", "TAC", bins, Tmin, Tmax);

    TCanvas *c5 = new TCanvas("c5", "sum_spectra");
    TH1D *hist_sum = new TH1D("hist_sum", "", bins, 0, 4000);

    bool goodToGo = 0;
    // by now we should check we have the same number of entries in every branch
    if (inbranchA->GetEntries() == inbranchB->GetEntries() && inbranchC->GetEntries() == inbranchD->GetEntries() && inbranchC->GetEntries() == inbranchA->GetEntries())
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
            

            inbranchC->GetEntry(i);
            double E2 = cal3(indataC.qlong); // pass to energy
            

            inbranchD->GetEntry(i);
            double T = TAC_cal(indataD.qlong);
            
            if ((indataA.timetag == indataB.timetag && indataB.timetag == indataC.timetag) && (indataC.timetag == indataD.timetag))
            {
                if (E0 > Emin && E0 < Emax && Tmin<T && T<Tmax) histA->Fill(E0);
                if (E1 > Emin && E1 < Emax && Tmin<T && T<Tmax) histB->Fill(E1);
                if (E2 > Emin && E2 < Emax && Tmin<T && T<Tmax) histC->Fill(E2);

                if ((((E0 > 250 && E0 < 430) || (E1 > 250 && E1 < 430) || (E2 > 250 && E2 < 430)) && Tmin<T) && T<Tmax)
                {
                    histD->Fill(T); // get the TAC entry for that possible coincidences
                    hist_sum->Fill(E0+E1+E2); // fill the sum histogram if any of the events is in the desired energy range
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

        c3->cd();
        histC->Rebin(32);
        histC->GetYaxis()->SetTitle("counts"); histC->GetXaxis()->SetTitle("Photon energy [keV]");
        histC->Draw();

        c4->cd();
        histD->Rebin(32);
        histD->GetYaxis()->SetTitle("counts"); histD->GetXaxis()->SetTitle("Delay [ns]");
        histD->Draw();

        c5->cd();
        hist_sum->Rebin(8);
        hist_sum->GetYaxis()->SetTitle("counts"); hist_sum->GetXaxis()->SetTitle("Photon energy sum [keV]");
        hist_sum->SetFillColorAlpha(kRed, 0.3);
        hist_sum->GetXaxis()->SetRangeUser(722, 1322);
        hist_sum->Draw();

        TF1 *gauss_plus_bg = new TF1("gauss_plus_bg","gaus(0)+pol1(3)", 900, 1122);
		gauss_plus_bg->SetParameters(4, 1022, 24);
		gauss_plus_bg->SetLineColor(kBlue);
        hist_sum->Fit("gauss_plus_bg", "R+");

        TF1 *fit = hist_sum->GetFunction("gauss_plus_bg");
        Double_t chi2 = fit->GetChisquare();
		Double_t Ndf = fit->GetNDF();

        TLatex text;
        text.DrawLatex(1100, 6.96053, Form("#color[4]{#chi^{2}/Ndf = %1.1f}", chi2/Ndf));
        text.SetTextSize(0.0842105);
    }
}