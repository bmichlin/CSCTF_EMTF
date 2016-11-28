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

void MakePlots(TString filename, int runNumber, int station, int ring, int numChambersInRing, int endcap){
	//Load File
	TFile *tfile = TFile::Open(filename);	

	//define a few variables
	if(endcap == 1) TString cap = "pos";
	if(endcap == -1) TString cap = "neg";

	if(endcap == 1) TString sign = "+";
	if(endcap == -1) TString sign = "-";

	//Define vectors that need to be used for plots	
	vector<double>  clbg_x; //clbg=CSCTF LCT BX Graph
	vector<double>  clbg_y;
	vector<double>  clbg_ex;
	vector<double>  clbg_ey;

	vector<double>  xpos_x; 
	vector<double>  xpos_y;
	vector<double>  xpos_ex;
	vector<double>  xpos_ey;

	vector<double>  ypos_x; 
	vector<double>  ypos_y;
	vector<double>  ypos_ex;
	vector<double>  ypos_ey;

	vector<double>  zpos_x; 
	vector<double>  zpos_y;
	vector<double>  zpos_ex;
	vector<double>  zpos_ey;

	vector<double> tof_bx_x; //These follow Pieter's TOF correction
	vector<double> tof_bx_y;
	vector<double> tof_bx_ex;
	vector<double> tof_bx_ey;

	vector<double> WA_RelBX_x; //This is for RelBX
	vector<double> WA_RelBX_y;
	vector<double> WA_RelBX_ex;
	vector<double> WA_RelBX_ey;

	vector<double> WA_indBX_low_x; //This is for the individual WalkAround BX (ie: the individual, NOT relative quantities).
	vector<double> WA_indBX_low_y;
	vector<double> WA_indBX_low_ex;
	vector<double> WA_indBX_low_ey;

	vector<double> WA_indBX_high_x; //This is for the individual WalkAround BX (ie: the individual, NOT relative quantities).
	vector<double> WA_indBX_high_y;
	vector<double> WA_indBX_high_ex;
	vector<double> WA_indBX_high_ey;

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


		TString bxHist = TString("BX_")+cap+TString::Format("ME_%d%d_%d_run%d",station,ring,x,runNumber);

		clbg_x.push_back(x);    
		clbg_y.push_back((getHist(tfile,bxHist)->GetMean()));
		clbg_ex.push_back( (0.5));
		clbg_ey.push_back(  (getHist(tfile, bxHist)->GetMeanError()));

		TString xpos = TString("xCoordinate_")+cap+TString::Format("_ME%d%d_%d",station,ring,x);
		xpos_x.push_back(x);
		xpos_y.push_back((getHist(tfile,xpos)->GetMean()));
		xpos_ex.push_back(0.5);
		xpos_ey.push_back(getHist(tfile,xpos)->GetMeanError());

		TString ypos = TString("yCoordinate_")+cap+TString::Format("_ME%d%d_%d",station,ring,x);
		ypos_x.push_back(x);
		ypos_y.push_back((getHist(tfile,ypos)->GetMean()));
		ypos_ex.push_back(0.5);
		ypos_ey.push_back(getHist(tfile,ypos)->GetMeanError());

		TString zpos = TString("zCoordinate_")+cap+TString::Format("_ME%d%d_%d",station,ring,x);
		zpos_x.push_back(x);
		zpos_y.push_back((getHist(tfile,zpos)->GetMean()));
		zpos_ex.push_back(0.5);
		zpos_ey.push_back(getHist(tfile,zpos)->GetMeanError());


		TString bxTOF = TString("TOF_BX_")+cap+TString::Format("_ME%d%d_%d",station,ring,x);

		tof_bx_x.push_back(x);
		tof_bx_y.push_back((getHist(tfile,bxTOF)->GetMean()));
		tof_bx_ex.push_back( (0.5) );
		tof_bx_ey.push_back((getHist(tfile,bxTOF)->GetMeanError()));

		TString waRelBX = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,endcap);

		WA_RelBX_x.push_back(lowChamber);
		WA_RelBX_y.push_back((getHist(tfile,waRelBX)->GetMean()));
		WA_RelBX_ex.push_back((0.5));
		WA_RelBX_ey.push_back((getHist(tfile,waRelBX)->GetMeanError()));



		TString waIndBX_low = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,station,ring,lowChamber,endcap);
		TString waIndBX_high = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,station,ring,highChamber,endcap);
		WA_indBX_low_x.push_back(lowChamber);
		WA_indBX_low_y.push_back((getHist(tfile,waIndBX_low)->GetMean()));
		WA_indBX_low_ex.push_back((0.5));
		WA_indBX_low_ey.push_back((getHist(tfile,waIndBX_low)->GetMeanError()));

		WA_indBX_high_x.push_back(lowChamber);
		WA_indBX_high_y.push_back((getHist(tfile,waIndBX_high)->GetMean()));
		WA_indBX_high_ex.push_back((0.5));
		WA_indBX_high_ey.push_back((getHist(tfile,waIndBX_high)->GetMeanError()));


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
	double* array_x = &clbg_x[0];
	double* array_y = &clbg_y[0];
	double* array_ex = &clbg_ex[0];
	double* array_ey = &clbg_ey[0];

	double* array_tof_bx_x = &tof_bx_x[0];
	double* array_tof_bx_y = &tof_bx_y[0];
	double* array_tof_bx_ex = &tof_bx_ex[0];
	double* array_tof_bx_ey = &tof_bx_ey[0];

	double* array_WA_Rel_BX_x = &WA_RelBX_x[0];
	double* array_WA_Rel_BX_y = &WA_RelBX_y[0];
	double* array_WA_Rel_BX_ex = &WA_RelBX_ex[0];
	double* array_WA_Rel_BX_ey = &WA_RelBX_ey[0];

	double* array_wa_indbx_low_x = &WA_indBX_low_x[0];
	double* array_wa_indbx_low_y = &WA_indBX_low_y[0];
	double* array_wa_indbx_low_ex = &WA_indBX_low_ex[0];
	double* array_wa_indbx_low_ey = &WA_indBX_low_ey[0];

	double* array_wa_indbx_high_x = &WA_indBX_high_x[0];
	double* array_wa_indbx_high_y = &WA_indBX_high_y[0];
	double* array_wa_indbx_high_ex = &WA_indBX_high_ex[0];
	double* array_wa_indbx_high_ey = &WA_indBX_high_ey[0];

	double* array_incAngle_x = &IncAngle_x[0];
	double* array_incAngle_y = &IncAngle_y[0];
	double* array_incAngle_ex = &IncAngle_ex[0];
	double* array_incAngle_ey = &IncAngle_ey[0];

	double* array_IncAngleNeighbors_x = &IncAngleNeighbors_x[0];
	double* array_IncAngleNeighbors_y = &IncAngleNeighbors_y[0];
	double* array_IncAngleNeighbors_ex = &IncAngleNeighbors_ex[0];
	double* array_IncAngleNeighbors_ey = &IncAngleNeighbors_ey[0];

	double* array_xpos_x = &xpos_x[0]; 
	double* array_xpos_y = &xpos_y[0];
	double* array_xpos_ex = &xpos_ex[0];
	double* array_xpos_ey = &xpos_ey[0];

	double* array_ypos_x = &ypos_x[0]; 
	double* array_ypos_y = &ypos_y[0];
	double* array_ypos_ex = &ypos_ex[0];
	double* array_ypos_ey = &ypos_ey[0];

	double* array_zpos_x = &zpos_x[0]; 
	double* array_zpos_y = &zpos_y[0];
	double* array_zpos_ex = &zpos_ex[0];
	double* array_zpos_ey = &zpos_ey[0];


	//Fill plots with arrays
	TGraphErrors* graph = new TGraphErrors(numChambersInRing, array_x,array_y, array_ex, array_ey);
	TGraphErrors* tof_bx_graph = new TGraphErrors(numChambersInRing,array_tof_bx_x,array_tof_bx_y,array_tof_bx_ex,array_tof_bx_ey);
	TGraphErrors* wa_RelBX_graph = new TGraphErrors(numChambersInRing,array_WA_Rel_BX_x,array_WA_Rel_BX_y,array_WA_Rel_BX_ex,array_WA_Rel_BX_ey);
	TGraphErrors* wa_indBX_low_graph = new TGraphErrors(numChambersInRing,array_wa_indbx_low_x,array_wa_indbx_low_y,array_wa_indbx_low_ex,array_wa_indbx_low_ey);
	TGraphErrors* wa_indBX_high_graph = new TGraphErrors(numChambersInRing,array_wa_indbx_high_x,array_wa_indbx_high_y,array_wa_indbx_high_ex,array_wa_indbx_high_ey);
	TGraphErrors* IncAngle_graph = new TGraphErrors(numChambersInRing, array_incAngle_x,array_incAngle_y,array_incAngle_ex,array_incAngle_ey);
	TGraphErrors* IncAngleNeighbors_graph = new TGraphErrors(numChambersInRing, array_IncAngleNeighbors_x,array_IncAngleNeighbors_y,array_IncAngleNeighbors_ex,array_IncAngleNeighbors_ey);
	TGraphErrors* xpos_graph = new TGraphErrors(numChambersInRing, array_xpos_x, array_xpos_y, array_xpos_ex, array_xpos_ey);
	TGraphErrors* ypos_graph = new TGraphErrors(numChambersInRing, array_ypos_x, array_ypos_y, array_ypos_ex, array_ypos_ey);
	TGraphErrors* zpos_graph = new TGraphErrors(numChambersInRing, array_zpos_x, array_zpos_y, array_zpos_ex, array_zpos_ey);


	//Stylistic Choices for Plots (and Plotting)	
	graph->SetTitle(TString("CSCTF LCT BX ME")+sign+Form("%d%d",station,ring));
	graph->SetMarkerColor(2);
	graph->SetMarkerStyle(21);
	graph->SetMarkerSize(1);	
	graph->GetXaxis()->SetTitle("Chamber Number");
	graph->GetYaxis()->SetTitle("BX (25ns)");	
	TCanvas* mycanvas = new TCanvas();
	mycanvas->SetGrid();
	TAxis *axis = graph->GetYaxis();
	axis->SetRangeUser(4,8);
	graph->Draw("AP");
	gStyle->SetOptFit(0001);
	cout << "" << endl;
	cout << "Fit information from the CSCTF LCT BX plot: " << endl;
	graph->Fit("pol0","F"); //The F option uses the minuit fitter.  Adding Q ("FQ") makes the fit not print.
	cout << "" << endl;


	tof_bx_graph->SetTitle(TString("CSCTF LCT BX With TOF Corrections ME")+sign+Form("%d%d",station,ring));
	tof_bx_graph->SetMarkerColor(2);
	tof_bx_graph->SetMarkerStyle(21);
	tof_bx_graph->SetMarkerSize(1);
	tof_bx_graph->GetXaxis()->SetTitle("Chamber Number");
	tof_bx_graph->GetYaxis()->SetTitle("BX (25ns)");
	TCanvas* mycanvas2 = new TCanvas();
	mycanvas2->SetGrid();
	TAxis *axis = tof_bx_graph->GetYaxis();
	axis->SetRangeUser(4,8);
	tof_bx_graph->Draw("AP");
	gStyle->SetOptFit(0001);
	cout << "" << endl;
	cout << "Fit information from the CSCTF LCT BX With TOF Corrections plot: " << endl;
	tof_bx_graph->Fit("pol0","F"); //The F option uses the minuit fitter.  Adding Q ("FQ") makes the fit not print.
	cout << "" << endl;


	wa_RelBX_graph->SetTitle(TString("BX for ME")+sign+Form("%d%d",station,ring));
	wa_RelBX_graph->SetMarkerColor(2);
	wa_RelBX_graph->SetMarkerStyle(21);
	wa_RelBX_graph->SetMarkerSize(1);
	wa_RelBX_graph->GetXaxis()->SetTitle("Chamber Number");
	wa_RelBX_graph->GetYaxis()->SetTitle("BX (25ns)");
	TCanvas* mycanvas3 = new TCanvas();
	mycanvas3->SetGrid();
	TAxis *axis = wa_RelBX_graph->GetYaxis();
	axis->SetRangeUser(-4,4);
	wa_RelBX_graph->Draw("AP");


	TCanvas* mycanvas4 = new TCanvas();
	mycanvas4->SetGrid();
	TMultiGraph *mg = new TMultiGraph();
	wa_indBX_low_graph->SetMarkerColor(2);
	wa_indBX_low_graph->SetMarkerStyle(21);
	wa_indBX_low_graph->SetMarkerSize(1);
	wa_indBX_high_graph->SetMarkerColor(4);
	wa_indBX_high_graph->SetMarkerStyle(21);
	wa_indBX_high_graph->SetMarkerSize(1);
	TString plotTitle = TString("Individual BX for ME")+sign+Form("%d%d",station,ring);
	mg->SetTitle(plotTitle+"; Chamber Number (n); BX (25ns)");
	mg->Add(wa_indBX_high_graph);
	mg->Add(wa_indBX_low_graph);
	mg->Draw("AP");
	TLegend* aleg = new TLegend(0.8232759,0.3771186,0.8994253,0.4978814,NULL,"brNDC");
	aleg->SetFillStyle(1);
	aleg->SetBorderSize(1);
	aleg->AddEntry(wa_indBX_high_graph, "n+1", "lp");
	aleg->AddEntry(wa_indBX_low_graph, "n", "lp");
	aleg->Draw("same");


	TCanvas* mycanvas1 = new TCanvas();
	mycanvas1->SetGrid();
	wa_RelBX_graph->SetMarkerColor(8);
	TMultiGraph *mg2 = new TMultiGraph();
	wa_indBX_low_graph->SetMarkerColor(2);
	wa_indBX_low_graph->SetMarkerStyle(21);
	wa_indBX_low_graph->SetMarkerSize(1);
	wa_indBX_high_graph->SetMarkerColor(4);
	wa_indBX_high_graph->SetMarkerStyle(21);
	wa_indBX_high_graph->SetMarkerSize(1);
	mg2->SetTitle(plotTitle+"; Chamber Number (n); BX (25ns)");
	mg2->Add(wa_indBX_high_graph);
	mg2->Add(wa_indBX_low_graph);
	mg2->Add(wa_RelBX_graph);
	mg2->Draw("AP");
	TLegend* bleg = new TLegend(0.5545977,0.3813559,0.8994253,0.5233051,NULL,"brNDC");
	bleg->SetFillStyle(1);
	bleg->SetBorderSize(1);
	bleg->AddEntry(wa_indBX_high_graph, "n+1", "lp");
	bleg->AddEntry(wa_indBX_low_graph, "n", "lp");
	bleg->AddEntry(wa_RelBX_graph, "Relative BX (BX_{n}-BX_{n+1})", "lp");
	bleg->Draw("same");

	IncAngle_graph->SetTitle(TString("Angle of Incident Muon for ME")+sign+Form("%d%d",station,ring));
	IncAngle_graph->SetMarkerColor(2);
	IncAngle_graph->SetMarkerStyle(21);
	IncAngle_graph->SetMarkerSize(1);
	IncAngle_graph->GetXaxis()->SetTitle("Chamber Number");
	IncAngle_graph->GetYaxis()->SetTitle("Angle (radian)");
	TCanvas* mycanvas10 = new TCanvas();
	mycanvas10->SetGrid();
	TAxis *axis10 = IncAngle_graph->GetYaxis();
	axis10->SetRangeUser(0,1.6);
	IncAngle_graph->Draw("AP");

	IncAngleNeighbors_graph->SetTitle(TString("Angle of Incident Muon for Neighboring Chambers in ME")+sign+Form("%d%d",station,ring));
	IncAngleNeighbors_graph->SetMarkerColor(2);
	IncAngleNeighbors_graph->SetMarkerStyle(21);
	IncAngleNeighbors_graph->SetMarkerSize(1);
	IncAngleNeighbors_graph->GetXaxis()->SetTitle("Chamber Number");
	IncAngleNeighbors_graph->GetYaxis()->SetTitle("Angle (radian)");
	TCanvas* mycanvas11 = new TCanvas();
	mycanvas11->SetGrid();
	TAxis *axis11 = IncAngleNeighbors_graph->GetYaxis();
	axis11->SetRangeUser(0,1.6);
	IncAngleNeighbors_graph->Draw("AP");

	xpos_graph->SetTitle(TString("x Coordinate for ME")+sign+Form("%d%d",station,ring));
	xpos_graph->SetMarkerColor(2);
	xpos_graph->SetMarkerStyle(21);
	xpos_graph->SetMarkerSize(1);
	xpos_graph->GetXaxis()->SetTitle("Chamber Number");
	xpos_graph->GetYaxis()->SetTitle("x-Position");
	TCanvas* mycanvas12 = new TCanvas();
	mycanvas12->SetGrid();
	xpos_graph->Draw("AP");

	ypos_graph->SetTitle(TString("y Coordinate for ME")+sign+Form("%d%d",station,ring));
	ypos_graph->SetMarkerColor(2);
	ypos_graph->SetMarkerStyle(21);
	ypos_graph->SetMarkerSize(1);
	ypos_graph->GetXaxis()->SetTitle("Chamber Number");
	ypos_graph->GetYaxis()->SetTitle("y-Position");
	TCanvas* mycanvas13 = new TCanvas();
	mycanvas13->SetGrid();
	ypos_graph->Draw("AP");

	zpos_graph->SetTitle(TString("z Coordinate for ME")+sign+Form("%d%d",station,ring));
	zpos_graph->SetMarkerColor(2);
	zpos_graph->SetMarkerStyle(21);
	zpos_graph->SetMarkerSize(1);
	zpos_graph->GetXaxis()->SetTitle("Chamber Number");
	zpos_graph->GetYaxis()->SetTitle("z-Position");
	TCanvas* mycanvas14 = new TCanvas();
	mycanvas14->SetGrid();
	zpos_graph->Draw("AP");

	TCanvas* mycanvas15 = new TCanvas();
	mycanvas15->SetGrid();
	getHist(tfile,"WalkAround_Rings_ME42_21_and_ME42_22_Overlap_BX_ME42_21_Endcap_1")->Draw();


	TCanvas* mycanvas16 = new TCanvas();
	mycanvas16->SetGrid();
	getHist(tfile,"WalkAroundwithTOF_Rings_ME42_21_and_ME42_22_Overlap_BX_ME42_21_Endcap_1")->Draw();

	TCanvas* mycanvas17 = new TCanvas();
	mycanvas17->SetGrid();
	getHist(tfile,"WalkAroundwithTOF_Rings_ME42_20_and_ME42_21_Overlap_BX_ME42_21_Endcap_1")->Draw();

	TCanvas* mycanvas18 = new TCanvas();
	mycanvas18->SetGrid();
	getHist(tfile,"WalkAround_Rings_ME42_20_and_ME42_21_Overlap_BX_ME42_21_Endcap_1")->Draw();

	mycanvas->SaveAs(TString::Format("BX_By_Chamber_ME%d%d.pdf",station,ring));
	mycanvas1->SaveAs(TString::Format("WalkAround_ME%d%d.pdf",station,ring));
	mycanvas2->SaveAs(TString::Format("BX_By_Chamber_TOF_ME%d%d.pdf",station,ring));
	mycanvas3->SaveAs(TString::Format("WalkAround_relBX_ME%d%d.pdf",station,ring));
	mycanvas4->SaveAs(TString::Format("WalkAround_High_Low_ME%d%d.pdf",station,ring));
	mycanvas10->SaveAs(TString::Format("IncidentAngle_ME%d%d.pdf",station,ring));
	mycanvas11->SaveAs(TString::Format("IncidentAngleNieghbors_ME%d%d.pdf",station,ring));
	mycanvas12->SaveAs(TString::Format("XCoordinate_By_Chamber_ME%d%d.pdf",station,ring));
	mycanvas13->SaveAs(TString::Format("YCoordinate_By_Chamber_ME%d%d.pdf",station,ring));
	mycanvas14->SaveAs(TString::Format("ZCoordinate_By_Chamber_ME%d%d.pdf",station,ring));
	mycanvas15->SaveAs("WalkAround_Rings_ME42_21_and_ME42_22_Overlap_BX_ME42_21_Endcap_1.pdf");
	mycanvas16->SaveAs("WalkAroundwithTOF_Rings_ME42_21_and_ME42_22_Overlap_BX_ME42_21_Endcap_1.pdf");
	mycanvas17->SaveAs("WalkAroundwithTOF_Rings_ME42_20_and_ME42_21_Overlap_BX_ME42_21_Endcap_1.pdf");
	mycanvas18->SaveAs("WalkAround_Rings_ME42_20_and_ME42_21_Overlap_BX_ME42_21_Endcap_1.pdf");
}
