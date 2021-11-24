


struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};



void TAC_calibration(){
    gROOT->SetStyle("Default");
    gStyle->SetOptStat(0); // clearing the f***ing TPaveStats
    gStyle->SetOptFit(111);
    gStyle->SetPalette(1);

    // desde canal 12000 hasta 18000
    const int length = 6;
    TString delays[length] = {"0ns", "3ns", "7ns", "15ns", "24ns", "56ns"};

    Double_t means[length] = {};
    Double_t stdDev[length] = {};
    
    for(Int_t i = 0; i<length; i++) {

        TString fileName = "data/day2/TAC_calibration_" + delays[i] + ".root";

        slimport_data_t indata;
        TFile *infile = new TFile(fileName);

        TTree *intree = (TTree*)infile->Get("acq_tree_0");
        

        TBranch *inbranch = intree->GetBranch("acq_ch3");
        inbranch->SetAddress(&indata.timetag);

        TCanvas *c1 = new TCanvas("c1", "c1");
        TH1D *hist = new TH1D("hist", "TAC", 40000, 0, 40000);

        for (Int_t ent = 0; ent<inbranch->GetEntries(); ent++)
        {
            inbranch->GetEntry(ent);
            hist->Fill(indata.qlong);
        }
        

        

        hist->Rebin(64);
        hist->Draw();

        // fitting
        TF1 *fun1, *fun2;
        Double_t mean, sigma;
        fun1 = new TF1("fun1name","gaus", 10000, 20000);
        hist->Fit(fun1,"RQN");
        mean  = fun1->GetParameter(1);
        sigma = fun1->GetParameter(2);
        fun2 = new TF1("fun2name","gaus",mean-(2.5*sigma),mean+(2.5*sigma));
        fun2->SetLineColor(kRed);
        hist->Fit(fun2,"R+");

        c1->SetLogy();
        c1->SaveAs(Form("preliminary_plots/TAC_calibration/%i.pdf", i));

        means[i] = fun2->GetParameter(1);
        stdDev[i] = fun2->GetParameter(2);
    }

    Double_t delay[length] = {0, 3, 7, 15, 24, 56};

    TCanvas *cal_can = new TCanvas("cal_can","cal_can");
    TGraph *cal_graph = new TGraph(length, means, delay);
    cal_graph->SetMarkerStyle(21);
	cal_graph->SetMarkerColor(4);
    TF1 *cal_fn = new TF1("cal_fn","[0]*x+[1]"); //defining a first order polynomial function and fit the data
	cal_graph->Fit(cal_fn);
    cal_graph->GetXaxis()->SetTitle("delay [ch]");
    cal_graph->GetYaxis()->SetTitle("delay [ns]");
    cal_graph->Draw();
	
}