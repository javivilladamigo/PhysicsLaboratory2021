#include <string>
#include <iostream>

//
// Simple root macro for gamma (or others) spectrum analysis 
//
// Run it from the command line:
// user@host:~/path$ root -l 'Mytool.C("spectrum_137Cs_60Co.Spe",3,300,4000)'
//   Don't forget the comas!!
// Run it inside a root sesion:
// root [0] .x Mytool.C("spectrum_137Cs_60Co.Spe",3,300,4000)
// 
// Arguments: 
// fileName = nome from the file to be analyzed (including file extension)
// npeaks = number of peaks to be searched for
// xMin,xMax = searching range (and drawing) in the histogram (in channels)
//

void Mytool(char* fileName,Int_t npeaks=3,Float_t xMin=400,Float_t xMax=1600) {

	//Style settings
	gROOT->SetStyle("Plain");
	gStyle->SetOptFit(111);
	gStyle->SetPalette(1);

	ifstream *file=new ifstream(fileName); //lee un un ficheiro, QUERO LER TODOS!!!!!

        cout << "\n \t ++++++++++ Searching for " << npeaks << " peaks..." << endl << endl;

        Int_t bins = 8192; //number of bins/channels of your histogram 

	Int_t data; //variable where to store data

	TH1I *hist = new TH1I("hist","hist",bins,0,bins); //data histogram

        for(Int_t i=0;i<12;i++) file->ignore(256,'\n'); //skipping the header of the ascii file (12 lines in this case)

	//Reading nbins/lines with data
        for(Int_t i=0;i<bins;i++) { 

          *file >> data;
          hist->SetBinContent(i+1,data);

        }

	TCanvas *can = new TCanvas("name","title");
        can->Divide(1,2);
        can->cd(1);
	hist->Draw();
        hist->Rebin(4); //rebinin for a better visualization and peak finding
        hist->GetXaxis()->SetRangeUser(xMin,xMax);

        TSpectrum *ss = new TSpectrum(npeaks); //creating TSpectrum object
        TH1I *hclone = (TH1I*)hist->Clone();
        can->cd(2);
        hclone->Draw();
        Int_t nfound = ss->Search(hist,10,"new",0.1); //searching peaks (TSpecrum class method)
        cout << "\n --------- Found " << nfound << " peaks" << endl << endl << endl;
        Double_t *xpeaks;
        xpeaks = (Double_t*) ss->GetPositionX(); //saving the position of the found peaks in an array 
        TH1 *hb = ss->Background(hclone,20,"same"); //calculating backgraund (TSpecrum class method)
        can->Update();

        hclone->Add(hclone,hb,1,-1); //substracting background from initial spectrum
        TCanvas *fin = new TCanvas("fin","fin");
        TH1I *hfinal = (TH1I*)hclone->Clone();
        hfinal->Draw();

        Double_t xp = 0.;
        Float_t mean = 0., sigma = 0.;

	Float_t resolution[3]={0.,0.,0.};

        //Sort xpeaks array
        int n = sizeof(xpeaks)/sizeof(xpeaks[0]);
        std::sort(xpeaks,xpeaks+nfound);

	//2-setp gauss fit to each of the found peaks, calculating the relative resolution afterwards
        TF1 *fun1, *fun2;
        for(Int_t i=0;i<nfound;i++)
        {   

          xp = xpeaks[i];
          cout << " ************ Peak at: " << xp << endl;
          fun1 = new TF1("fun1name","gaus",xp-120,xp+120);
          hfinal->Fit(fun1,"RQN");
          mean  = fun1->GetParameter(1);
          sigma = fun1->GetParameter(2);
          fun2 = new TF1("fun2name","gaus",mean-(1.5*sigma),mean+(1.5*sigma));
          hfinal->Fit(fun2,"R+");
          mean  = fun2->GetParameter(1);
          sigma = fun2->GetParameter(2);

	  resolution[i]=((2.35*sigma)/mean)*100;

          cout << "\nRESULTS:  Mean at: " << mean << " with sigma " << sigma << " , FWHM of " << 2.35*sigma << " , and resolution of " << resolution[i] << " %" << endl << endl;
        
        }

	//Fit of the resolution values for each peak to a function Res~1/sqrt(E), and evaluating it at 1 MeV
	Float_t xval[3];
	if(npeaks==1) xval[0] = 0.662; else if(npeaks==2) { xval[0] = 1.173; xval[1] = 1.332;} else if(npeaks==3) { xval[0] = 0.662; xval[1] = 1.173; xval[2] = 1.332;} 
	TGraph *graph = new TGraph(npeaks,xval,resolution);
	TCanvas *c_graph = new TCanvas("c_graph","c_graph");
	graph->SetMarkerStyle(20);
	graph->Draw("AP");
	TF1 *reso = new TF1("reso","[0]/x+[1]/sqrt(x)+[2]");
	graph->Fit("reso");
	Float_t res_at_one = reso->Eval(1.);
	cout << "\n\n ************** RESOLUTION @ 1MeV:  " << res_at_one << " % ************** " << endl << endl;

	//Calibration
	Float_t peaks[npeaks];
	for(Int_t i=0;i<npeaks;i++) peaks[i] = xpeaks[i];
	TGraph *cal_graph = new TGraph(npeaks,peaks,xval); //calibration: channel vs energy of the found peaks
	cal_graph->SetNameTitle("cal_graph","Calibration");
	cal_graph->SetMarkerStyle(21);
	cal_graph->SetMarkerColor(4);
	TF1 *cal_fn = new TF1("cal_fn","[0]*x+[1]"); //defining a first order polynomial function and fit the data
	cal_graph->Fit(cal_fn);
	TCanvas *cal_can = new TCanvas("cal_can","cal_can");
	cal_graph->Draw();

	//Re-escale the histogram X-axis with the calculated calibration
	Float_t new_end = cal_fn->Eval(bins);
	TH1F *hcal = new TH1F("hcal","hcal",hist->GetNbinsX(),0,new_end); //create calibrated histo!!
	for(Int_t i=0;i<bins;i++) hcal->SetBinContent(i+1,hist->GetBinContent(i+1)); //fill cal histo with same content!
	TCanvas *cal_histo_can = new TCanvas("cal_histo_can","cal_histo_can");
	hcal->Draw();
	hcal->GetXaxis()->SetTitle("Photon energy [MeV]");
	hcal->GetXaxis()->CenterTitle();


}
