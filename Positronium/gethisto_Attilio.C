//Digitizer data from the LAB

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};

TH1D* getHistoFromTree(const char *name_file, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch("acq_ch3");
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h_spectrum->Fill(indata.qlong);
	}
	// return
	return h_spectrum;
}

TH1D* getHistoForChannelFromTree(const char *name_file, short chan, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h_spectrum->Fill(indata.qlong);
	}
	// return
	return h_spectrum;
}


TH1D* getHistoWithFilter(const char *name_file, int numBins, double minX, double maxX, double lowThr = 0, double highThr = 999999) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		if (indata.qlong>lowThr && indata.qlong<highThr) {
			h_spectrum->Fill(indata.qlong);
		}
	}
	// return
	return h_spectrum;
}


TGraph *getSignal(const char *name_file, int nSamples=250, short chan=0, int nrEv=1){
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TGraph *graph = new TGraph();
	
	//Setting the desired event
	inbranch->GetEntry(nrEv);

	//Looping over the samples
	for (int i=0; i<nSamples; ++i){
		graph->SetPoint(i, i, indata.samples[i]);
	}
	return graph;
}

TH1D* getRefinedHistoFromTree(const char *name_file, int numBins, double minX, double maxX, double lowThr=0, double higThr=999999) {
	// variables
	slimport_data_t indata0, indata1, indata2, indata3;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch0 = intree->GetBranch("acq_ch0");
	TBranch *inbranch1=intree->GetBranch("acq_ch1");
	TBranch *inbranch2=intree->GetBranch("acq_ch2");
	TBranch *inbranch3=intree->GetBranch("acq_ch3");
	inbranch0->SetAddress(&indata0.timetag);
	inbranch1->SetAddress(&indata1.timetag);
    	inbranch2->SetAddress(&indata2.timetag);
    	inbranch3->SetAddress(&indata3.timetag);    
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch3->GetEntries(); i++) {
		inbranch0->GetEntry(i);
		double k0=indata0.qlong;
		double E0=-18.9+k0*0.119077;
		
		inbranch1->GetEntry(i);
		double k1=indata1.qlong;
		double E1=-28.3+k1*0.115687;
		
		inbranch2->GetEntry(i);
		double k2=indata2.qlong;
		double E2=-25.3+k2*0.110853;
		
		inbranch3->GetEntry(i);
		double T;
		T=((indata3.qlong-11250)/85.)*1000;
		
		if(E0+E1>950 && E0+E1<1090)//this limits are just for fun, they should be adjusted
		{
		inbranch3->GetEntry(i);
			if(indata3.qlong>lowThr && indata3.qlong<higThr){
			h_spectrum->Fill(T);
			}
		}
	}
	// return
	return h_spectrum;
}

