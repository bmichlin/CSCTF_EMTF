#include <math.h>

TH1F* getHist(TFile* file, TString hname){
	TIter next(file->GetListOfKeys());
	TKey *key;
	while((key = (TKey*)next())){
		TDirectory* dir = file->GetDirectory(key->GetName());
		if(dir->FindKey(hname) != 0){
			TH1F *h = (TH1F*)dir->Get(hname);
			return h;
		}
	}
	//cout << hname << endl;
	return NULL;
}

double roundDouble(double value){

	double val = value*1000;
	double num = (int)val;
	double output = num/1000;
	//cout << "val, div, output: " << val << ", " << div << ", " << output << endl;
	return output;
}

void printTable(double *value, double *error){

	cout << "+-----------------------------------------------+" << endl;
	cout << "|\tChamber\t|\t\tCorrection\t|" << endl;
	for(int i =0; i < 36; i++){
		cout << "|\t" << i+1 << "\t|\t" << roundDouble(value[i]) << "\t+/- "  << " " << roundDouble(error[i]) << "\t|" << endl;
	}
	cout << "+-----------------------------------------------+" << endl;

}

void ChiSquarePlots(TString filename, int runNumber, int station, int ring, int numChambersInRing, int endcap){

double inputChiValues[36] = {-0.0809471, 
0.0220545, 
0.00219286, 
0.0556113, 
0.0475148, 
0.158269, 
0.0465512, 
0.155559, 
0.124524, 
0.15363, 
0.0248635, 
0.131979, 
0.12787, 
0.246296, 
0.160614, 
0.116149, 
0.0593604, 
0.129642, 
-0.0931235, 
-0.152684, 
-0.199359, 
-0.195445, 
-0.16198, 
-0.187031, 
-0.237379, 
-0.198947, 
-0.171845, 
-0.166389, 
-0.156949, 
-0.185, 
-0.1477, 
-0.10938, 
-0.065182, 
-0.0415139, 
0.0245867, 
0.0239319};

double inputChiErrors[36] = {0.0448514, 
0.043628, 
0.0416011, 
0.0404225, 
0.0390291, 
0.038018, 
0.0368951, 
0.0363442, 
0.0361403, 
0.0362033, 
0.0365041, 
0.0372255, 
0.0380233, 
0.0401305, 
0.0407451, 
0.0427709, 
0.0440956, 
0.0465246, 
0.0475629, 
0.0482525, 
0.0487611, 
0.0489634, 
0.0488342, 
0.0486882, 
0.0483945, 
0.0482302, 
0.048064, 
0.0479439, 
0.0477545, 
0.0474891, 
0.0473868, 
0.0472509, 
0.0472151, 
0.0470214, 
0.0465479, 
0.0457038};


//	double inputChiValues[36] = {-0.428327, 
		/*{{{*/
/*		-0.440006, 
		-0.351737, 
		-0.0330415, 
		-0.0844775, 
		0.340695, 
		0.0968491, 
		0.229901, 
		0.396996, 
		0.369743, 
		0.177973, 
		1.21075, 
		0.455162, 
		0.530616, 
		0.643555, 
		0.56129, 
		0.506935, 
		0.528321, 
		0.206211, 
		0.105172, 
		0.151464, 
		0.0979389, 
		-0.233295, 
		-0.206126, 
		-0.203086, 
		-0.314195, 
		-0.475122, 
		-0.280711, 
		-0.551098, 
		-0.611501, 
		-0.544541, 
		-0.61929, 
		-0.573261, 
		-0.572481, 
		-0.531746, 
		-0.384724};
		*/
		/*}}}*/

	//double inputChiErrors[36] = {0.105823, 
		/*{{{*/
	/*
	0.102682, 
		0.101035, 
		0.0961825, 
		0.0950807, 
		0.0902251, 
		0.0885011, 
		0.08595, 
		0.084719, 
		0.0837659, 
		0.0830174, 
		0.0826727, 
		0.0822851, 
		0.0821744, 
		0.0824333, 
		0.0830177, 
		0.0834895, 
		0.0840564, 
		0.0856365, 
		0.0866472, 
		0.0865308, 
		0.0863351, 
		0.0859117, 
		0.085956, 
		0.0866027, 
		0.0876516, 
		0.0880984, 
		0.0891224, 
		0.0902623, 
		0.0911986, 
		0.0950467, 
		0.0970292, 
		0.100203, 
		0.101305, 
		0.105165, 
		0.106045};
	*/	
		/*}}}*/
	//Load File
	TFile *tfile = TFile::Open(filename);	

	//define a few variables
	if(endcap == 1) TString cap = "pos";
	if(endcap == -1) TString cap = "neg";

	if(endcap == 1) TString sign = "+";
	if(endcap == -1) TString sign = "-";

	//Define vectors that need to be used for plots	
	vector<double> chi_low_x; //These are the chi corrected values
	vector<double> chi_low_y;
	vector<double> chi_low_ex;
	vector<double> chi_low_ey;

	vector<double> chi_high_x; //These are the chi corrected values
	vector<double> chi_high_y;
	vector<double> chi_high_ex;
	vector<double> chi_high_ey;

	vector<double> chi_RelBX_x; 
	vector<double> chi_RelBX_y;
	vector<double> chi_RelBX_ex;
	vector<double> chi_RelBX_ey;

	vector<double> chi_corr_x;
	vector<double> chi_corr_y;
	vector<double> chi_corr_ex;
	vector<double> chi_corr_ey;

	vector<double>  clbg_x; //clbg=CSCTF LCT BX Graph with Chi Correction
	vector<double>  clbg_y;
	vector<double>  clbg_ex;
	vector<double>  clbg_ey;

	vector<double> tof_bx_x; //These follow Pieter's TOF correction & Chi Correction
	vector<double> tof_bx_y;
	vector<double> tof_bx_ex;
	vector<double> tof_bx_ey;

	vector<double> IncAngle_x; // This is for the angle of the incident muon
	vector<double> IncAngle_y;
	vector<double> IncAngle_ex;
	vector<double> IncAngle_ey;

	vector<double> IncAngleNeighbors_x; // This is for the angle of the incident muon
	vector<double> IncAngleNeighbors_y;
	vector<double> IncAngleNeighbors_ex;
	vector<double> IncAngleNeighbors_ey;

	for(int x = 1; x < numChambersInRing+1; x++){

		int highChamber = 0; //These chambers are generic and can be used for any of the plots
		int lowChamber = 0;
		if(x == numChambersInRing){
			int lowChamber =  x;
			int highChamber = 1;
		}else{
			int lowChamber = x;
			int highChamber = x+1; 
		}


		TString chiIndBX_low = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,station,ring,lowChamber,endcap);
		TString chiIndBX_high = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,station,ring,highChamber,endcap);

		//Adding errors in quadrature
		chi_low_x.push_back(lowChamber);
		chi_low_y.push_back((getHist(tfile,chiIndBX_low)->GetMean()) + inputChiValues[lowChamber-1] );
		chi_low_ex.push_back((0.5));
		chi_low_ey.push_back(  sqrt( (getHist(tfile,chiIndBX_low)->GetMeanError()*getHist(tfile,chiIndBX_low)->GetMeanError()) + (inputChiErrors[lowChamber-1]* inputChiErrors[lowChamber-1]) )  );

		chi_high_x.push_back(lowChamber);
		chi_high_y.push_back((getHist(tfile,chiIndBX_high)->GetMean()) + inputChiValues[highChamber-1]);
		chi_high_ex.push_back((0.5));
		chi_high_ey.push_back( sqrt( (getHist(tfile,chiIndBX_high)->GetMeanError()*getHist(tfile,chiIndBX_high)->GetMeanError()) + (inputChiErrors[highChamber-1]*inputChiErrors[highChamber-1])   ));

		TString chiRelBX = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,endcap);
		chi_RelBX_x.push_back(lowChamber);
		chi_RelBX_y.push_back( ((getHist(tfile,chiIndBX_low)->GetMean()) + inputChiValues[lowChamber-1] )-(getHist(tfile,chiIndBX_high)->GetMean() + inputChiValues[highChamber-1]) );
		chi_RelBX_ex.push_back((0.5));
		double er_ey1 = sqrt( (getHist(tfile,chiIndBX_low)->GetMeanError()*getHist(tfile,chiIndBX_low)->GetMeanError()) + (inputChiErrors[lowChamber-1]* inputChiErrors[lowChamber-1]) );
		double er_ey2 = sqrt( (getHist(tfile,chiIndBX_high)->GetMeanError()*getHist(tfile,chiIndBX_high)->GetMeanError()) + (inputChiErrors[highChamber-1]*inputChiErrors[highChamber-1]) );
		double er_tot = sqrt( (er_ey1*er_ey1) + (er_ey2*er_ey2) );
		chi_RelBX_ey.push_back(er_tot );



		chi_corr_x.push_back(lowChamber);
		chi_corr_y.push_back(inputChiValues[lowChamber-1]);
		chi_corr_ex.push_back(0.5);
		chi_corr_ey.push_back(inputChiErrors[lowChamber-1]);


		TString bxHist = TString("BX_")+cap+TString::Format("ME_%d%d_%d_run%d",station,ring,x,runNumber);

		clbg_x.push_back(x);
		clbg_y.push_back((getHist(tfile,bxHist)->GetMean()) + inputChiValues[lowChamber-1]);
		clbg_ex.push_back( (0.5));
		clbg_ey.push_back(  sqrt( (getHist(tfile, bxHist)->GetMeanError()*getHist(tfile, bxHist)->GetMeanError()) + (inputChiErrors[lowChamber-1]*inputChiErrors[lowChamber-1])) );


		TString bxTOF = TString("TOF_BX_")+cap+TString::Format("_ME%d%d_%d",station,ring,x);

		tof_bx_x.push_back(x);
		tof_bx_y.push_back((getHist(tfile,bxTOF)->GetMean()) + inputChiValues[lowChamber-1]);
		tof_bx_ex.push_back( (0.5) );
		tof_bx_ey.push_back( sqrt( (getHist(tfile,bxTOF)->GetMeanError()*getHist(tfile,bxTOF)->GetMeanError()) + (inputChiErrors[lowChamber-1]*inputChiErrors[lowChamber-1])) );


		TString incAngle = TString("IncidentAngle_")+cap+TString::Format("_ME%d%d_%d",station,ring,x);

		IncAngle_x.push_back(x);
		IncAngle_y.push_back(getHist(tfile, incAngle)->GetMean());
		IncAngle_ex.push_back(0.5);
		IncAngle_ey.push_back(getHist(tfile, incAngle)->GetMeanError());

		TString incAngle = TString("IncidentAngleNeighbors_")+cap+TString::Format("_ME%d%d_%d",station,ring,x);

		IncAngleNeighbors_x.push_back(x);
		IncAngleNeighbors_y.push_back(getHist(tfile, incAngle)->GetMean());
		IncAngleNeighbors_ex.push_back(0.5);
		IncAngleNeighbors_ey.push_back(getHist(tfile, incAngle)->GetMeanError());
	}

	//Turn vectors in to arrays of data
	double* array_chi_low_x = &chi_low_x[0];
	double* array_chi_low_y = &chi_low_y[0];
	double* array_chi_low_ex = &chi_low_ex[0];
	double* array_chi_low_ey = &chi_low_ey[0];

	double* array_chi_high_x = &chi_high_x[0];
	double* array_chi_high_y = &chi_high_y[0];
	double* array_chi_high_ex = &chi_high_ex[0];
	double* array_chi_high_ey = &chi_high_ey[0];

	double* array_chi_Rel_BX_x = &chi_RelBX_x[0];
	double* array_chi_Rel_BX_y = &chi_RelBX_y[0];
	double* array_chi_Rel_BX_ex = &chi_RelBX_ex[0];
	double* array_chi_Rel_BX_ey = &chi_RelBX_ey[0];

	double* array_chi_corr_x = &chi_corr_x[0];
	double* array_chi_corr_y = &chi_corr_y[0];
	double* array_chi_corr_ex = &chi_corr_ex[0];
	double* array_chi_corr_ey = &chi_corr_ey[0];

	double* array_x = &clbg_x[0];
	double* array_y = &clbg_y[0];
	double* array_ex = &clbg_ex[0];
	double* array_ey = &clbg_ey[0];

	double* array_tof_bx_x = &tof_bx_x[0];
	double* array_tof_bx_y = &tof_bx_y[0];
	double* array_tof_bx_ex = &tof_bx_ex[0];
	double* array_tof_bx_ey = &tof_bx_ey[0];

double* array_incAngle_x = &IncAngle_x[0];
double* array_incAngle_y = &IncAngle_y[0];
double* array_incAngle_ex = &IncAngle_ex[0];
double* array_incAngle_ey = &IncAngle_ey[0];

double* array_IncAngleNeighbors_x = &IncAngleNeighbors_x[0];
double* array_IncAngleNeighbors_y = &IncAngleNeighbors_y[0];
double* array_IncAngleNeighbors_ex = &IncAngleNeighbors_ex[0];
double* array_IncAngleNeighbors_ey = &IncAngleNeighbors_ey[0];



	//Fill plots with arrays
	TGraphErrors* chi_indBX_low_graph = new TGraphErrors(numChambersInRing, array_chi_low_x, array_chi_low_y, array_chi_low_ex, array_chi_low_ey);
	TGraphErrors* chi_indBX_high_graph = new TGraphErrors(numChambersInRing, array_chi_high_x, array_chi_high_y, array_chi_high_ex, array_chi_high_ey);
	TGraphErrors* chi_RelBX_graph = new TGraphErrors(numChambersInRing,array_chi_Rel_BX_x,array_chi_Rel_BX_y,array_chi_Rel_BX_ex,array_chi_Rel_BX_ey);
	TGraphErrors* chi_corr_graph = new TGraphErrors(numChambersInRing,array_chi_corr_x,array_chi_corr_y,array_chi_corr_ex,array_chi_corr_ey);
	TGraphErrors* graph = new TGraphErrors(numChambersInRing, array_x,array_y, array_ex, array_ey);
	TGraphErrors* tof_bx_graph = new TGraphErrors(numChambersInRing,array_tof_bx_x,array_tof_bx_y,array_tof_bx_ex,array_tof_bx_ey);

TGraph2DErrors* angleBX_graph = new TGraph2DErrors(numChambersInRing, array_tof_bx_x,array_tof_bx_y, array_IncAngleNeighbors_x,array_IncAngleNeighbors_y,array_IncAngleNeighbors_ex,array_IncAngleNeighbors_ey);
	//Stylistic Choices for Plots (and Plotting)	

	TString plotTitle = TString("BX (with ChiSquare Correction) for ME")+sign+Form("%d%d",station,ring);

/*{{{*/
char* xlabel[36]={"Chambers 1 & 2",
"Chambers 2 & 3",
"Chambers 3 & 4",
"Chambers 4 & 5",
"Chambers 5 & 6",
"Chambers 6 & 7",
"Chambers 7 & 8",
"Chambers 8 & 9",
"Chambers 9 & 10",
"Chambers 10 & 11",
"Chambers 11 & 12",
"Chambers 12 & 13",
"Chambers 13 & 14",
"Chambers 14 & 15",
"Chambers 15 & 16",
"Chambers 16 & 17",
"Chambers 17 & 18",
"Chambers 18 & 19",
"Chambers 19 & 20",
"Chambers 20 & 21",
"Chambers 21 & 22",
"Chambers 22 & 23",
"Chambers 23 & 24",
"Chambers 24 & 25",
"Chambers 25 & 26",
"Chambers 26 & 27",
"Chambers 27 & 28",
"Chambers 28 & 29",
"Chambers 29 & 30",
"Chambers 30 & 31",
"Chambers 31 & 32",
"Chambers 32 & 33",
"Chambers 33 & 34",
"Chambers 34 & 35",
"Chambers 35 & 36",
"Chambers 36 & 1"};
/*}}}*/

TH1F *h = new TH1F("h","test",36,4.5,39.5);
for(k=1;k<numChambersInRing;k++) h->GetXaxis()->SetBinLabel(k,xlabel[k-1]);

	TCanvas* mycanvas5 = new TCanvas();
	mycanvas5->SetGrid();
	chi_RelBX_graph->SetMarkerColor(8);
	chi_RelBX_graph->SetMarkerStyle(21);
	chi_RelBX_graph->SetMarkerSize(1);
	TMultiGraph *mg2 = new TMultiGraph();
	chi_indBX_low_graph->SetMarkerColor(2);
	chi_indBX_low_graph->SetMarkerStyle(21);
	chi_indBX_low_graph->SetMarkerSize(1);
	chi_indBX_high_graph->SetMarkerColor(4);
	chi_indBX_high_graph->SetMarkerStyle(21);
	chi_indBX_high_graph->SetMarkerSize(1);
	mg2->SetTitle(plotTitle+"; Chamber Number (n); BX (25ns)");
	mg2->Add(chi_indBX_high_graph);
	mg2->Add(chi_indBX_low_graph);
	mg2->Add(chi_RelBX_graph);
	//h->Draw(); //This will put in the horrible x-axis
	mg2->Draw("AP");

	//Fit the relative difference
	//gStyle->SetOptFit(0001);
	//chi_RelBX_graph->Fit("pol0","F");

	TLegend* bleg = new TLegend(0.5545977,0.3813559,0.8994253,0.5233051,NULL,"brNDC");
	bleg->SetFillStyle(1);
	bleg->SetBorderSize(1);
	bleg->AddEntry(chi_indBX_high_graph, "n+1", "lp");
	bleg->AddEntry(chi_indBX_low_graph, "n", "lp");
	bleg->AddEntry(chi_RelBX_graph, "Relative BX (BX_{n}-BX_{n+1})", "lp");
	bleg->Draw("same");


	chi_corr_graph->SetTitle(TString("Suggested Timing Correction to ME")+sign+Form("%d%d",station,ring));
	chi_corr_graph->SetMarkerColor(2);
	chi_corr_graph->SetMarkerStyle(21);
	chi_corr_graph->SetMarkerSize(1);
	chi_corr_graph->GetXaxis()->SetTitle("Chamber Number");
	chi_corr_graph->GetYaxis()->SetTitle("BX (25ns)");
	TCanvas* mycanvas6 = new TCanvas();
	mycanvas6->SetGrid();
	//TAxis *axis = graph->GetYaxis();
	//axis->SetRangeUser(4,8);
	chi_corr_graph->Draw("AP");


	graph->SetTitle(TString("CSCTF LCT BX With Chi Correction ME")+sign+Form("%d%d",station,ring));
	graph->SetMarkerColor(2);
	graph->SetMarkerStyle(21);
	graph->SetMarkerSize(1);
	graph->GetXaxis()->SetTitle("Chamber Number");
	graph->GetYaxis()->SetTitle("BX (25ns)");
	TCanvas* mycanvas7 = new TCanvas();
	mycanvas7->SetGrid();
	TAxis *axis = graph->GetYaxis();
	axis->SetRangeUser(4,8);
	graph->Draw("AP");

	tof_bx_graph->SetTitle(TString("CSCTF LCT BX With TOF and Chi Corrections ME")+sign+Form("%d%d",station,ring));
	tof_bx_graph->SetMarkerColor(2);
	tof_bx_graph->SetMarkerStyle(21);
	tof_bx_graph->SetMarkerSize(1);
	tof_bx_graph->GetXaxis()->SetTitle("Chamber Number");
	tof_bx_graph->GetYaxis()->SetTitle("BX (25ns)");
	TCanvas* mycanvas8 = new TCanvas();
	mycanvas8->SetGrid();
	TAxis *axis = tof_bx_graph->GetYaxis();
	axis->SetRangeUser(4,8);
	tof_bx_graph->Draw("AP");
	gStyle->SetOptFit(0001);
	//cout << "" << endl;
	//cout << "Fit information from the CSCTF LCT BX With Pieter's Corrections plot: " << endl;
	tof_bx_graph->Fit("pol0","F"); //The F option uses the minuit fitter.  Adding Q ("FQ") makes the fit not print.
	//cout << "" << endl;

	angleBX_graph->SetTitle(TString("Angle By BX with TOF and Chi Corrections for ME")+sign+Form("%d%d",station,ring));
	angleBX_graph->SetMarkerColor(2);
	angleBX_graph->SetMarkerStyle(21);
	angleBX_graph->SetMarkerSize(1);
	//angleBX_graph->GetXaxis()->SetTitle("Chamber Number");
	//angleBX_graph->GetYaxis()->SetTitle("BX (25ns)");
	TCanvas* mycanvas11 = new TCanvas();
	mycanvas11->SetGrid();
	angleBX_graph->Draw("err p0");

	mycanvas5->SaveAs(TString::Format("WalkAround_Chi_ME%d%d.pdf",station,ring));
	mycanvas6->SaveAs(TString::Format("Chi_Timing_Correction_ME%d%d.pdf",station,ring));
	mycanvas7->SaveAs(TString::Format("BX_By_Chamber_Chi_ME%d%d.pdf",station,ring));
	mycanvas8->SaveAs(TString::Format("BX_By_Chamber_TOF_Chi_ME%d%d.pdf",station,ring));

	printTable(inputChiValues, inputChiErrors);
}
