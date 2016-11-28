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

void Matrix_Sync(TString filename, int runNumber, int station, int ring, const int numChambersInRing, int endcap, double* inputChiValues, double* inputChiErrors){
    
    //Load File
    TFile *tfile = TFile::Open(filename);
    //Creates the array of mean time differences marray and the array of corresponding errors sarray
    double marray [numChambersInRing];
    double sarray [numChambersInRing];
    for(int i = 1; i<numChambersInRing+1; i++){
        int highChamber = 0; //These chambers are generic and can be used for any of the plots
        int lowChamber = 0;
        if(i == numChambersInRing){
            int lowChamber =  i;
            int highChamber = 1;
        }else{
            int lowChamber = i;
            int highChamber = i+1;
        }
        TString waRelBX = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,endcap);
        marray[i-1] = (getHist(tfile,waRelBX)->GetMean());
        sarray[i-1] = (getHist(tfile,waRelBX)->GetMeanError());
    }
    
    //Instantiate vectors s2, square of each of the errors, and m, the mean time differences
    vector<double> s2;
    vector<double> m;
    for(int i = 0; i < 36; i++){
        s2.push_back(sarray[i]);
        m.push_back(marray[i]);
    }
    //Square each entry
    for (int i=0;i<=35;i++){
        s2[i] = s2[i]*s2[i];
    }
    
    //Create A matrix and b vector
    TMatrixD A(36,36);
    TVectorD b(36);
    //Make a Matrix of all zeros
    for(int i =0; i < 35; i++){
        for(int j =0; j < 36; j++){
            A(i,j) = 0.000;
        }
    }
    //Fill in the matrix A
    A(0,0) = s2[35] + s2[0];
    A(0,1) = -s2[35];
    A(0,35) = -s2[0];
    A(35,0) = -s2[34];
    A(35,34) = -s2[35];
    A(35,35) = s2[34] + s2[35];
    for (int i=1;i<35;i++){
        A(i,i-1) = -s2[i];
        A(i,i) = s2[i-1] + s2[i];
        A(i,i+1) = -s2[i-1];
    }
    //Fill in the vector b
    b[0] = s2[0]*m[35] - s2[35]*m[0];
    for (int i=1;i<36;i++){
        b[i] = s2[i]*m[i-1] - s2[i-1]*m[i];
    }
    
    //In the timing note x_pseudo is referred to as the Greek letter delta
    TMatrixD Aw = A;
    TVectorD yw = b;
    TDecompSVD svd(Aw);
    
    //Solve Ax=b via TDecompSVD pseudoinverse method
    const TMatrixD pseudo1  = svd.Invert();
    TVectorD x_pseudo = yw;
    x_pseudo *= pseudo1;
    
    double temp = (m[35]+x_pseudo[35]-x_pseudo[0]);
    double chi2 = (temp*temp)/(s2[35]);
    for(int i=0;i<35;i++){
        if(s2[i] != 0){
            temp = (m[i]+x_pseudo[i]-x_pseudo[i+1]);
            chi2 = chi2 + (temp*temp)/(s2[i]);
        }
    }
    std::cout << "Minimum chi-squared value: " << std::endl;
    printf("%.20f\n\n",chi2);
    
    
    TVectorD x_errors(36);
    for(int i = 0; i < 36; i++){
        for(int j = 0; j < 36; j++){
            int plusOne = 0; //These chambers are generic and can be used for any of the plots
            int minusOne = 0;
            if(j == 0){
                int minusOne =  35;
                int plusOne = 1;
            }
            else if(j==35){
                int minusOne = 34;
                int plusOne = 0;
            }
            else{
                int minusOne = j-1;
                int plusOne = j+1;
            }
            x_errors(i) = x_errors(i) + s2[j]*(pseudo1(i,plusOne)*s2[plusOne]-pseudo1(i,j)*s2[minusOne])*(pseudo1(i,plusOne)*s2[plusOne]-pseudo1(i,j)*s2[minusOne]);
        }
    }
    x_errors.Sqrt();
    
    for(int i =0; i<numChambersInRing; i++){
        inputChiValues[i] = x_pseudo(i);
        inputChiErrors[i] = x_errors(i);
    }
}

double roundDouble(double value){

	double val = value*1000;
	double num = (int)val;
	double output = num/1000;
	//cout << "val, div, output: " << val << ", " << div << ", " << output << endl;
	return output;
}

void printTable(double *value, double *error, int runNumber,  int endcap, int station, int ring){
ofstream writeFile(TString::Format("Correction_Table_%d_ME%d_%d%d.txt",runNumber, endcap, station, ring),ios::app);
	cout << "+-----------------------------------------------+" << endl;
	cout << "|\tChamber\t|\t\tCorrection\t|" << endl;
writeFile << "+-----------------------------------------------+\n";
writeFile << "|\tChamber\t|\t\tCorrection\t|\n";
	for(int i =0; i < 36; i++){
		cout << "|\t" << i+1 << "\t|\t" << roundDouble(value[i]) << "\t+/- "  << " " << roundDouble(error[i]) << "\t|" << endl;
writeFile << "|\t" << i+1 << "\t|\t" << roundDouble(value[i]) << "\t+/- "  << " " << roundDouble(error[i]) << "\t|\n";	
}
	cout << "+-----------------------------------------------+" << endl;
writeFile << "+-----------------------------------------------+\n";

}

void ChiSquarePlots(TString filename, int runNumber, int station, int ring, int numChambersInRing, int endcap){

    double inputChiValues[36];

    double inputChiErrors[36];
    
    Matrix_Sync(filename, runNumber, station, ring, numChambersInRing, endcap, inputChiValues, inputChiErrors);
    
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
	TGraphErrors* chi_RelBX_graph_fit = new TGraphErrors(numChambersInRing,array_chi_Rel_BX_x,array_chi_Rel_BX_y,array_chi_Rel_BX_ex,array_chi_Rel_BX_ey);
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
for(int k=1;k<numChambersInRing;k++) h->GetXaxis()->SetBinLabel(k,xlabel[k-1]);

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
	bleg->SetTextSize(0.037);
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
	gStyle->SetOptFit(0001);
	graph->Fit("pol0","F");
	graph->GetFunction("pol0")->SetLineColor(2);

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

	chi_RelBX_graph_fit->SetTitle(TString("Relative BX with Chi Corrections for ME")+sign+Form("%d%d",station,ring));
	chi_RelBX_graph_fit->SetMarkerColor(8);
	chi_RelBX_graph_fit->SetMarkerStyle(21);
	chi_RelBX_graph_fit->SetMarkerSize(1);
	chi_RelBX_graph_fit->GetXaxis()->SetTitle("Chamber Number");
	chi_RelBX_graph_fit->GetYaxis()->SetTitle("BX (25ns)");
	TCanvas* mycanvas9 = new TCanvas();
	mycanvas9->SetGrid();
	TAxis *axis = chi_RelBX_graph_fit->GetYaxis();
	axis->SetRangeUser(-4,4);

	chi_RelBX_graph_fit->Draw("AP");
	gStyle->SetOptFit(0001);
	chi_RelBX_graph_fit->Fit("pol0","F");
	chi_RelBX_graph_fit->GetFunction("pol0")->SetLineColor(8);

	mycanvas5->SaveAs(cap+TString::Format("Matrix_WalkAround_Chi_ME%d%d_%d.pdf",station,ring,runNumber));
	mycanvas6->SaveAs(cap+TString::Format("Matrix_Chi_Timing_Correction_ME%d%d_%d.pdf",station,ring,runNumber));
	mycanvas7->SaveAs(cap+TString::Format("Matrix_BX_By_Chamber_Chi_ME%d%d_%d.pdf",station,ring,runNumber));
	mycanvas8->SaveAs(cap+TString::Format("Matrix_BX_By_Chamber_TOF_Chi_ME%d%d_%d.pdf",station,ring,runNumber));
	mycanvas9->SaveAs(cap+TString::Format("Matrix_WalkAround_RelBX_withChi_ME%d%d_%d.pdf",station,ring,runNumber));

	printTable(inputChiValues, inputChiErrors, runNumber, endcap, station, ring);
}
