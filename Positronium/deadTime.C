
Double_t nonParalyzableDeadTime(Double_t true_rate, Double_t observed_rate) {
    Double_t deadTime = (true_rate-observed_rate)/(true_rate*observed_rate);
    return deadTime;
}

Double_t standardDeadTime(Double_t true_rate, Double_t observed_rate) {
    Double_t deadTime = (1-observed_rate/true_rate)*100;
    return deadTime;
}

void deadTime() {


    // 2 gammas with trigger on D1&D2
    Double_t counts_D1andD2_2gammas_scaler = 56852;
    Double_t time_D1andD2_2gammas_scaler = 10.;

    Double_t counts_D1andD2_2gammas_VERDI = 664007;
    Double_t time_D1andD2_2gammas_VERDI = 121.;


    // 2 gammas with trigger on (D1&D2)&D4
    Double_t counts_D1andD2andD4_2gammas_scaler = 608;
    Double_t time_D1andD2andD4_2gammas_scaler = 10.;

    Double_t counts_D1andD2andD4_2gammas_VERDI = 10048;
    Double_t time_D1andD2andD4_2gammas_VERDI = 183.;


    // 3 gammas (trigger on D1&D2D3)
    Double_t counts_D1andD2andD3_3gammas_scaler = 3278;
    Double_t time_D1andD2andD3_3gammas_scaler = 600.;

    Double_t counts_D1andD2andD3_3gammas_VERDI = 398987;
    Double_t time_D1andD2andD3_3gammas_VERDI = 20*3600+22*60+17;

    gROOT->SetStyle("Default");
    gStyle->SetOptStat(0); // clearing the f***ing TPaveStats
    gStyle->SetOptFit(111);
    gStyle->SetPalette(1);

    Double_t x_D1andD2_2gammas[2] = {1. , 1.};
    Double_t D1andD2_2gammas[2] = {counts_D1andD2_2gammas_scaler/time_D1andD2_2gammas_scaler, counts_D1andD2_2gammas_VERDI/time_D1andD2_2gammas_VERDI};
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->SetLogy();
    TGraph *D1andD2_2gammas_graph = new TGraph(2, x_D1andD2_2gammas, D1andD2_2gammas);
    D1andD2_2gammas_graph->SetMarkerStyle(21);
	D1andD2_2gammas_graph->SetMarkerColor(4);
    D1andD2_2gammas_graph->SetLineWidth(0);
    D1andD2_2gammas_graph->GetYaxis()->SetRangeUser(1, 10000);
    D1andD2_2gammas_graph->GetXaxis()->SetLimits(0, 4);
    D1andD2_2gammas_graph->GetXaxis()->SetTitle("Arrangements");
    D1andD2_2gammas_graph->GetYaxis()->SetTitle("Rate");
    D1andD2_2gammas_graph->Draw();

    D1andD2_2gammas_graph->Draw();


    Double_t x_D1andD2andD4_2gammas[2] = {2. , 2.};
    Double_t D1andD2andD4_2gammas[2] = {counts_D1andD2andD4_2gammas_scaler/time_D1andD2andD4_2gammas_scaler, counts_D1andD2andD4_2gammas_VERDI/time_D1andD2andD4_2gammas_VERDI};
    TGraph *D1andD2andD4_2gammas_graph = new TGraph(2, x_D1andD2andD4_2gammas, D1andD2andD4_2gammas);
    D1andD2andD4_2gammas_graph->SetMarkerStyle(21);
	D1andD2andD4_2gammas_graph->SetMarkerColor(2);
    D1andD2andD4_2gammas_graph->SetLineWidth(0);

    D1andD2andD4_2gammas_graph->Draw("PSAME");


    Double_t x_D1andD2andD3_3gammas[2] = {3. , 3.};
    Double_t D1andD2andD3_3gammas[2] = {counts_D1andD2andD3_3gammas_scaler/time_D1andD2andD3_3gammas_scaler, counts_D1andD2andD3_3gammas_VERDI/time_D1andD2andD3_3gammas_VERDI};
    TGraph *D1andD2andD3_3gammas_graph = new TGraph(2, x_D1andD2andD3_3gammas, D1andD2andD3_3gammas);
    D1andD2andD3_3gammas_graph->SetMarkerStyle(21);
	D1andD2andD3_3gammas_graph->SetMarkerColor(3);
    D1andD2andD3_3gammas_graph->SetLineWidth(0);
    D1andD2andD3_3gammas_graph->Draw("PSAME");
    

    TLegend *legend = new TLegend(0.72,0.656842,0.982808, 0.995);
    legend->AddEntry(D1andD2_2gammas_graph, "2 gammas (D1&D2)","p");
    legend->AddEntry(D1andD2andD4_2gammas_graph, "2 gammas (D1&D2)&D4","p");
    legend->AddEntry(D1andD2andD3_3gammas_graph, "3gammas","p");
    legend->Draw();


    TCanvas *c2 = new TCanvas("c2", "c2");
    TH2D *hist = new TH2D("hist", "hist", 1000, 0, 4, 1000, 0, 10);
    hist->GetXaxis()->SetTitle("arrangements");
    hist->GetYaxis()->SetTitle("dead time [%]");
    hist->Draw();

    Double_t deadTime_2gammas_D1andD2 = standardDeadTime(counts_D1andD2_2gammas_scaler/time_D1andD2_2gammas_scaler, counts_D1andD2_2gammas_VERDI/time_D1andD2_2gammas_VERDI);
    cout << "Dead time for 2 gammas configuration and D1&D2: " << deadTime_2gammas_D1andD2 << " per cent" << endl;

    Double_t deadTime_2gammas_D1andD2andD4 = standardDeadTime(counts_D1andD2andD4_2gammas_scaler/time_D1andD2andD4_2gammas_scaler, counts_D1andD2andD4_2gammas_VERDI/time_D1andD2andD4_2gammas_VERDI);
    cout << "Dead time for 2 gammas configuration and D1&D2&D4: " << deadTime_2gammas_D1andD2andD4 << " per cent" << endl;

    Double_t deadTime_3gammas_D1andD2andD3 = standardDeadTime(counts_D1andD2andD3_3gammas_scaler/time_D1andD2andD3_3gammas_scaler, counts_D1andD2andD3_3gammas_VERDI/time_D1andD2andD3_3gammas_VERDI);
    cout << "Dead time for 3 gammas configuration and D1&D2&D3: " << deadTime_3gammas_D1andD2andD3 << " per cent" << endl;

    Double_t x3[3] = {1,2,3};
    Double_t deadTimes[3] = {deadTime_2gammas_D1andD2, deadTime_2gammas_D1andD2andD4, deadTime_3gammas_D1andD2andD3};
    

    TMarker *dT_2g_D1D2_marker = new TMarker(x3[0], deadTimes[0],21);
    TMarker *dT_2g_D1D2D4_marker = new TMarker(x3[1], deadTimes[1],20);
    TMarker *dT_3g_D1D2D3_marker = new TMarker(x3[2], deadTimes[2],22);

    dT_2g_D1D2_marker->Draw();
    dT_2g_D1D2D4_marker->Draw();
    dT_3g_D1D2D3_marker->Draw();

    TLegend *legend2 = new TLegend(0.72,0.656842,0.982808, 0.995);       
    legend2->AddEntry(dT_2g_D1D2_marker, "2 gammas (D1&D2)","p");
    legend2->AddEntry(dT_2g_D1D2D4_marker, "2 gammas (D1&D2)&D4","p");
    legend2->AddEntry(dT_3g_D1D2D3_marker, "3gammas","p");
    legend2->Draw();


}

