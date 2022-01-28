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
Double_t TAC_cal(Int_t ch) {return 0.01179*ch-140.2 + 30;} // 0.01179*ch-140.2 real one

void correlated_events(const char *name_file = "data/day3/3gamma_filteredEvents.root", int bins = 16384, double xMin = 128, double xMax = 16384) {
    
    gROOT->SetStyle("Default");
    //gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
	gStyle->SetStatBorderSize(1);

    // variables
    Double_t angular_aperture = TMath::ATan(5./15.);
    
    Double_t Tmin = 0; Double_t Tmax = 600;
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
    TH1D *histA = new TH1D("histA", "", bins, cal1(xMin), cal1(xMax));
    TH1D *unfiltered_D1 = (TH1D*)histA->Clone();


    TCanvas *c2 = new TCanvas("c2", "Detector #2");
    TH1D *histB = new TH1D("histB", "", bins, cal2(xMin), cal2(xMax));
    TH1D *unfiltered_D2 = (TH1D*)histB->Clone();


    TCanvas *c3 = new TCanvas("c3", "Detector #3");
    TH1D *histC = new TH1D("histC", "", bins, cal3(xMin), cal3(xMax));
    TH1D *unfiltered_D3 = (TH1D*)histC->Clone();
    

    TCanvas *c4 = new TCanvas("c4", "TAC");
    TH1D *histD = new TH1D("histD", "", bins, Tmin, Tmax);
    TH1D *histD_clone = (TH1D*)histD->Clone();

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
            unfiltered_D1->Fill(E0);

            inbranchB->GetEntry(i);
            double E1 = cal2(indataB.qlong); // pass to energy
            unfiltered_D2->Fill(E1);
            

            inbranchC->GetEntry(i);
            double E2 = cal3(indataC.qlong); // pass to energy
            unfiltered_D3->Fill(E2);
            

            inbranchD->GetEntry(i);
            double T = TAC_cal(indataD.qlong);
            
            if ((indataA.timetag == indataB.timetag && indataB.timetag == indataC.timetag) && (indataC.timetag == indataD.timetag))
            {
                if (E0 > Emin && E0 < Emax && Tmin<T && T<Tmax) histA->Fill(E0);
                if (E1 > Emin && E1 < Emax && Tmin<T && T<Tmax) histB->Fill(E1);
                if (E2 > Emin && E2 < Emax && Tmin<T && T<Tmax) histC->Fill(E2);

                if ((((E0 > Emin && E0 < Emax) || (E1 > Emin && E1 < Emax) || (E2 > Emin && E2 < Emax))) && Tmin<T && T<Tmax)
                {
                    histD->Fill(T); // get the TAC entry for that possible coincidences
                    hist_sum->Fill(E0+E1+E2); // fill the sum histogram if any of the events is in the desired energy range
                    if ((E0+E1+E2)<(1017.63+3*26.1516) && (E0+E1+E2) > (1017.63-3*26.1516))
                    {
                        histD_clone->Fill(T); // filling the timetags of the events that are later going to be counted as coincidences
                    }
                    
                }
            }         
        }

        c1->cd();
        c1->SetLogy();
        histA->Rebin(64);
        unfiltered_D1->Rebin(64);
        unfiltered_D1->GetYaxis()->SetTitle("counts"); unfiltered_D1->GetXaxis()->SetTitle("Photon energy [keV]");
        unfiltered_D1->GetYaxis()->SetRangeUser(1., 1e5);
        histA->SetLineColor(kBlue);
        histA->SetFillColorAlpha(kBlue, 0.4);
        unfiltered_D1->Draw();
        histA->Draw("SAME");
        
        
        

        c2->cd();
        c2->SetLogy();
        histB->Rebin(64);
        unfiltered_D2->Rebin(64);
        unfiltered_D2->GetYaxis()->SetTitle("counts"); unfiltered_D2->GetXaxis()->SetTitle("Photon energy [keV]");
        unfiltered_D2->GetYaxis()->SetRangeUser(1., 1e5);
        histB->SetLineColor(kBlue);
        histB->SetFillColorAlpha(kBlue, 0.4);
        unfiltered_D2->Draw();
        histB->Draw("SAME");

        c3->cd();
        c3->SetLogy();
        histC->Rebin(64);
        unfiltered_D3->Rebin(64);
        unfiltered_D3->GetYaxis()->SetTitle("counts"); unfiltered_D3->GetXaxis()->SetTitle("Photon energy [keV]");
        unfiltered_D3->GetYaxis()->SetRangeUser(1., 1e5);
        histC->SetLineColor(kBlue);
        histC->SetFillColorAlpha(kBlue, 0.4);
        unfiltered_D3->Draw();
        histC->Draw("SAME");

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
                histD_clone->Fill(i);
            }
        }
        histD->GetXaxis()->SetRangeUser(0, 200);
        histD->SetLineColor(kBlue);
        histD->SetLineWidth(2.);
        histD->Draw();
        //histD_clone->Draw("SAME");

        TF1 *expo_free = new TF1("expo_freename", "expo+expo(2)", 25, 200);
        expo_free->SetParLimits(1, -0.025, -0.005);
        expo_free->SetParLimits(3, -1, -0.025);

        expo_free->SetLineColor(kRed);
        expo_free->SetLineWidth(4.);

        histD->Fit("expo_freename", "R+");
        Double_t lifetime1 = -1/expo_free->GetParameter(1); // in ns
        Double_t error_lifetime1 = expo_free->GetParErrors()[1]/(expo_free->GetParameter(1)*expo_free->GetParameter(1));
        Double_t lifetime2 = -1/expo_free->GetParameter(3); // in ns
        Double_t error_lifetime2 = expo_free->GetParErrors()[3]/(expo_free->GetParameter(3)*expo_free->GetParameter(3));
        cout << endl << "First lifetime comp: " << lifetime1 << " +- " << error_lifetime1 << endl;
        cout << "Second lifetime comp: " << lifetime2 << " +- " << error_lifetime2 << endl;
        Double_t chi2_expo = expo_free->GetChisquare();
		Double_t Ndf_expo = expo_free->GetNDF();


        cout << "Lifetime fit has a chi2/Ndf: " << chi2_expo/Ndf_expo << endl;

        TLatex text3;
        text3.DrawLatex(100, 100, Form("#color[2]{#chi^{2}/Ndf = %1.2f}", chi2_expo/Ndf_expo));
        text3.SetTextSize(0.0671141);

        TLatex text4;
        text4.DrawLatex(100, 90, Form("#color[2]{#tau_{1} = %1.1f(%1.1f) ns}", lifetime1, error_lifetime1));
        text4.SetTextSize(0.0671141);

        TLatex text5;
        text5.DrawLatex(100, 80, Form("#color[2]{#tau_{2} = %1.2f(%2.0f) ns}", lifetime2, error_lifetime2*100));
        text5.SetTextSize(0.0671141);


        c5->cd();
        hist_sum->GetYaxis()->SetTitleOffset(0.9);
        hist_sum->GetXaxis()->SetTitleSize(0.045);
        hist_sum->GetYaxis()->SetTitleSize(0.045);

        hist_sum->Rebin(64);
        hist_sum->GetYaxis()->SetTitle("counts"); hist_sum->GetXaxis()->SetTitle("Photon energy sum [keV]");
        hist_sum->SetLineColor(kBlue);
        hist_sum->SetLineWidth(2);
        //hist_sum->GetYaxis()->SetRangeUser(0, 80); 
        hist_sum->GetXaxis()->SetRangeUser(380, 2400);
        hist_sum->Draw();

        TF1 *gauss_plus_bg = new TF1("gauss_plus_bg","gaus(0)+pol1(3)", Emin*3, Emax*3);
		gauss_plus_bg->SetParameters(4, 1022, 23);
        gauss_plus_bg->SetParLimits(1, 950, 1100);
        gauss_plus_bg->SetParLimits(2, 0, 10000);
		gauss_plus_bg->SetLineColor(kRed);
        gauss_plus_bg->SetLineWidth(8.);
        hist_sum->Fit("gauss_plus_bg", "RQ");
        

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
        
        TLatex text;
        text.DrawLatex(1100, 6.96053, Form("#color[2]{#chi^{2}/Ndf = %1.2f}", chi2/Ndf));
        text.SetTextSize(0.06);

        TLatex text2;
        text2.DrawLatex(1100, 6.0, Form("#color[2]{Integral (#pm 3#sigma) = %.0f events}", nb_of_events));
        text2.SetTextSize(0.06);

    }
}