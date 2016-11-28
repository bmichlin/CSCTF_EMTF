#include "Math/GSLMinimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TH2.h"
#include <iomanip>
/*
#include "TH1.h"
#include "TFile.h"
#include "TKey.h"


TString filename;
int endcap, station, ring, chamber, InnumChambersInRing;


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
*/

double chisquare(const double *xx){

	const Double_t x0 = xx[0];
	/*{{{*/
	const Double_t x1 = xx[1];
	const Double_t x2 = xx[2];
	const Double_t x3 = xx[3];
	const Double_t x4 = xx[4];
	const Double_t x5 = xx[5];
	const Double_t x6 = xx[6];
	const Double_t x7 = xx[7];
	const Double_t x8 = xx[8];
	const Double_t x9 = xx[9];
	const Double_t x10 = xx[10];
	const Double_t x11 = xx[11];
	const Double_t x12 = xx[12];
	const Double_t x13 = xx[13];
	const Double_t x14 = xx[14];
	const Double_t x15 = xx[15];
	const Double_t x16 = xx[16];
	const Double_t x17 = xx[17];
	const Double_t x18 = xx[18];
	const Double_t x19 = xx[19];
	const Double_t x20 = xx[20];
	const Double_t x21 = xx[21];
	const Double_t x22 = xx[22];
	const Double_t x23 = xx[23];
	const Double_t x24 = xx[24];
	const Double_t x25 = xx[25];
	const Double_t x26 = xx[26];
	const Double_t x27 = xx[27];
	const Double_t x28 = xx[28];
	const Double_t x29 = xx[29];
	const Double_t x30 = xx[30];
	const Double_t x31 = xx[31];
	const Double_t x32 = xx[32];
	const Double_t x33 = xx[33];
	/*}}}1*/
	const Double_t x34 = xx[34];
	const Double_t x35 = xx[35];

	//my variable names are poor. xx[0] corresponds to chamber 1. xx[int+1]==chamber number.
	//I did this so that I could follow root tutorial

	//double chisquare is defined in the folded lines below

	double chisquare =pow(  ( x0 - x1 + 0.137255 )  / 0.0621269 ,2)  +  
pow(  ( x1 - x2 + -0.0588235 )  / 0.0570672 ,2)  +  
pow(  ( x2 - x3 + 0.129032 )  / 0.0882473 ,2)  +  
pow(  ( x3 - x4 + -0.2 )  / 0.178885 ,2)  +  
pow(  ( x4 - x5 + -0.1 )  / 0.0974679 ,2)  +  
//pow(  ( x5 - x6 + 0 )  / 0 ,2)  +  
pow(  ( x6 - x7 + 0.0769231 )  / 0.0739053 ,2)  +  
//pow(  ( x7 - x8 + 0 )  / 0 ,2)  +  
pow(  ( x8 - x9 + -0.111111 )  / 0.0740741 ,2)  +  
//pow(  ( x9 - x10 + 0 )  / 0 ,2)  +  
//pow(  ( x10 - x11 + 0 )  / 0 ,2)  +  
//pow(  ( x11 - x12 + 0 )  / 0 ,2)  +  
//pow(  ( x12 - x13 + 0 )  / 0 ,2)  +  
//pow(  ( x13 - x14 + 0 )  / 0 ,2)  +  
pow(  ( x14 - x15 + 0.2 )  / 0.189737 ,2)  +  
pow(  ( x15 - x16 + 0 )  / 0.176777 ,2)  +  
pow(  ( x16 - x17 + 0.103448 )  / 0.0565523 ,2)  +  
pow(  ( x17 - x18 + -0.133333 )  / 0.0877707 ,2)  +  
pow(  ( x18 - x19 + 0.0322581 )  / 0.0966006 ,2)  +  
//pow(  ( x19 - x20 + 0 )  / 0 ,2)  +  
pow(  ( x20 - x21 + -0.25 )  / 0.216506 ,2)  +  
pow(  ( x21 - x22 + 0.333333 )  / 0.136083 ,2)  +  
//pow(  ( x22 - x23 + 0 )  / 0 ,2)  +  
//pow(  ( x23 - x24 + 0 )  / 0 ,2)  +  
pow(  ( x24 - x25 + 0 )  / 0.235702 ,2)  +  
//pow(  ( x25 - x26 + 0 )  / 0 ,2)  +  
pow(  ( x26 - x27 + 0.333333 )  / 0.272166 ,2)  +  
pow(  ( x27 - x28 + 1 )  / 0 ,2)  +  
//pow(  ( x28 - x29 + 0 )  / 0.282843 ,2)  +  
pow(  ( x29 - x30 + 0.333333 )  / 0.272166 ,2)  +  
pow(  ( x30 - x31 + 0.2 )  / 0.126491 ,2)  +  
//pow(  ( x31 - x32 + 0 )  / 0 ,2)  +  
pow(  ( x32 - x33 + 0.107143 )  / 0.0922962 ,2)  +  
pow(  ( x33 - x34 + 0.2 )  / 0.126491 ,2)  +  
pow(  ( x34 - x35 + -0.0833333 )  / 0.0815788 ,2)  +  
pow(  ( x35 - x0 + 0.047619 )  / 0.0818214 ,2); 
	/*
	   TFile *tfile = TFile::Open(filename);
	   double chisquare;
	   for(int chamber = 1; chamber < 37; chamber++){
	   int highChamber = 0; //These chambers are generic and can be used for any of the plots
	   int lowChamber = 0;
	   if(chamber == InnumChambersInRing){
	   int lowChamber =  chamber;
	   int highChamber = 1;
	   }else{
	   int lowChamber = chamber;
	   int highChamber = chamber+1;
	   }

	   TString waRelBX = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,endcap);
	   if(chamber != InnumChambersInRing){
	   if(getHist(tfile, waRelBX)->GetMeanError() != 0){

	   string First = Form("x%d",chamber-1);
	   string Second = Form("x%d",chamber);

	   chisquare += pow( ( First.c_str()- Second.c_str() + getHist(tfile, waRelBX)->GetMean()  )  / getHist(tfile, waRelBX)->GetMeanError() ,2);
	   }
	//myString += TString("pow(  ( ") + TString::Format("x%d",chamber-1) + TString(" - ") + TString::Format("x%d",chamber) + TString(" + ") + TString("getHist(tfile, waRelBX)->GetMean())") + TString("  / ") + TString("getHist(tfile, waRelBX)->GetMeanError()") + TString(" ,2)  +  ");

	}else{
	TString waRelBX = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,endcap);
	//	myString += TString("pow(  ( ") +  TString::Format("x%d",chamber-1) + TString(" - ") + TString("x0") + TString(" + ") + TString(",getHist(tfile, waRelBX)->GetMean())") + TString("  / ") + TString("getHist(tfile, waRelBX)->GetMeanError())") +TString(" ,2);  ");


	}
	}

	cout << chisquare << endl;

*/


	return chisquare;
}
///
void MinimizeBFGS(int numChambersInRing, double stepSize, double variableVal/*, TString filename, int runNumber, int station, int ring, int numChambersInRing, int endcap*/){
	//Minimizer
	/*	ROOT::Math::GSLMinimizer min(ROOT::Math::kVectorBFGS);

		min.SetMaxFunctionCalls(10000);
		min.SetMaxIterations(10000);
		min.SetTolerance(0.01);

		ROOT::Math::Functor f(&chisquare, numChambersInRing);

		vector<double> stepVector;
		vector<double> variableVector;

		for(int i =0; i < numChambersInRing; i++){
		stepVector.push_back(stepSize);
		variableVector.push_back(variableVal);
		}

		double* step = &stepVector[0];
		double* variable = &variableVector[0];

		min.SetFunction(f);
		min.SetVariable(0, "x0", variable[0], step[0]);
		min.SetVariable(1, "x1", variable[1], step[1]);
		min.SetVariable(2, "x2", variable[2], step[2]);
		min.SetVariable(3, "x3", variable[3], step[3]);
		min.SetVariable(4, "x4", variable[4], step[4]);
		min.SetVariable(5, "x5", variable[5], step[5]);
		min.SetVariable(6, "x6", variable[6], step[6]);
		min.SetVariable(7, "x7", variable[7], step[7]);
		min.SetVariable(8, "x8", variable[8], step[8]);
		min.SetVariable(9, "x9", variable[9], step[9]);
		min.SetVariable(10, "x10", variable[10], step[10]);
		min.SetVariable(11, "x11", variable[11], step[11]);
		min.SetVariable(12, "x12", variable[12], step[12]);
		min.SetVariable(13, "x13", variable[13], step[13]);
		min.SetVariable(14, "x14", variable[14], step[14]);
		min.SetVariable(15, "x15", variable[15], step[15]);
		min.SetVariable(16, "x16", variable[16], step[16]);
		min.SetVariable(17, "x17", variable[17], step[17]);
		if(numChambersInRing != 18){	
		min.SetVariable(18, "x18", variable[18], step[18]);
		min.SetVariable(19, "x19", variable[19], step[19]);
		min.SetVariable(20, "x20", variable[20], step[20]);
		min.SetVariable(21, "x21", variable[21], step[21]);
		min.SetVariable(22, "x22", variable[22], step[22]);
		min.SetVariable(23, "x23", variable[23], step[23]);
		min.SetVariable(24, "x24", variable[24], step[24]);
		min.SetVariable(25, "x25", variable[25], step[25]);
		min.SetVariable(26, "x26", variable[26], step[26]);
		min.SetVariable(27, "x27", variable[27], step[27]);
		min.SetVariable(28, "x28", variable[28], step[28]);
		min.SetVariable(29, "x29", variable[29], step[29]);
		min.SetVariable(30, "x30", variable[30], step[30]);
		min.SetVariable(31, "x31", variable[31], step[31]);
		min.SetVariable(32, "x32", variable[32], step[32]);
		min.SetVariable(33, "x33", variable[33], step[33]);
		min.SetVariable(34, "x34", variable[34], step[34]);
		min.SetVariable(35, "x35", variable[35], step[35]);
		}	

		min.Minimize();
		const double *xs = min.X();
		cout << "Optimal Values (that minimize chi-square): " << endl;
		for(int j = 0; j < numChambersInRing; j++){
		cout << "Chamber " << j+1 << ": " << xs[j] << endl;
		}
		cout << "Value of minimum: " << chisquare(xs) << endl; 
		*/
}
///
void MinimizeMinuit(/*TString Infilename, int InStation, int InRing, int InEndcap,*/ int numChambersInRing, int searchLow, int searchHigh, double stepSize, int printFlag){

	//int searchLow = -1000;
	//int searchHigh = 1000;	

	TH2D *minPlot = new TH2D("minPlot","minPlot", 1000, searchLow, searchHigh, 200000, 1.8806528818530438179e+20, 1.8806528818530458179e+20); 
	minPlot->SetMarkerColor(kRed);	
	double minValue = 9999999999999999999999999999999999999999999;
	double pOfMin;
	for(double p = searchLow; p < searchHigh; p+= stepSize){
		//Minimizer
		ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );

		min.SetMaxFunctionCalls(1000000);
		min.SetMaxIterations(1000000);
		min.SetTolerance(0.001);
		min.SetPrintLevel(1);
		min.SetPrintLevel(0);
		/*
		   filename = Infilename;
		   station = InStation;
		   ring = InRing;
		   endcap = InEndcap;
		   InnumChambersInRing = numChambersInRing;  
		   */

		ROOT::Math::Functor f(&chisquare, numChambersInRing);

		vector<double> stepVector;
		vector<double> variableVector;

		for(int i =0; i < numChambersInRing; i++){
			stepVector.push_back(stepSize);
			variableVector.push_back(p);
			//variableVector.push_back(variableVal);
		}

		double* step = &stepVector[0];
		double* variable = &variableVector[0];

		min.SetFunction(f);
		min.SetVariable(0, "x0", variable[0], step[0]);
		/*{{{1*/	
		min.SetVariable(1, "x1", variable[1], step[1]);
		min.SetVariable(2, "x2", variable[2], step[2]);
		min.SetVariable(3, "x3", variable[3], step[3]);
		min.SetVariable(4, "x4", variable[4], step[4]);
		min.SetVariable(5, "x5", variable[5], step[5]);
		min.SetVariable(6, "x6", variable[6], step[6]);
		min.SetVariable(7, "x7", variable[7], step[7]);
		min.SetVariable(8, "x8", variable[8], step[8]);
		min.SetVariable(9, "x9", variable[9], step[9]);
		min.SetVariable(10, "x10", variable[10], step[10]);
		min.SetVariable(11, "x11", variable[11], step[11]);
		min.SetVariable(12, "x12", variable[12], step[12]);
		min.SetVariable(13, "x13", variable[13], step[13]);
		min.SetVariable(14, "x14", variable[14], step[14]);
		min.SetVariable(15, "x15", variable[15], step[15]);
		min.SetVariable(16, "x16", variable[16], step[16]);
		min.SetVariable(17, "x17", variable[17], step[17]);
		if(numChambersInRing != 18){	
			min.SetVariable(18, "x18", variable[18], step[18]);
			min.SetVariable(19, "x19", variable[19], step[19]);
			min.SetVariable(20, "x20", variable[20], step[20]);
			min.SetVariable(21, "x21", variable[21], step[21]);
			min.SetVariable(22, "x22", variable[22], step[22]);
			min.SetVariable(23, "x23", variable[23], step[23]);
			min.SetVariable(24, "x24", variable[24], step[24]);
			min.SetVariable(25, "x25", variable[25], step[25]);
			min.SetVariable(26, "x26", variable[26], step[26]);
			min.SetVariable(27, "x27", variable[27], step[27]);
			min.SetVariable(28, "x28", variable[28], step[28]);
			min.SetVariable(29, "x29", variable[29], step[29]);
			min.SetVariable(30, "x30", variable[30], step[30]);
			min.SetVariable(31, "x31", variable[31], step[31]);
			min.SetVariable(32, "x32", variable[32], step[32]);
			min.SetVariable(33, "x33", variable[33], step[33]);
			min.SetVariable(34, "x34", variable[34], step[34]);
			min.SetVariable(35, "x35", variable[35], step[35]);
		}	
		/*}}}*/	

		min.Minimize();
		//if(printFlag !=0) cout << min.Minimize() << endl;
		const double *xs = min.X();
		const double *error = min.Errors();
		double runningSum = 0;
//cout << setprecision(25) << chisquare(xs)*pow(10,18) << endl;
cout << "Minimum: " << setprecision(20) << chisquare(xs)*pow(10,18) << endl;	
	minPlot->Fill(p, chisquare(xs)*pow(10,18));
		if(minValue > chisquare(xs)){
			minValue = chisquare(xs);	
			pOfMin = p;
		}
		//cout << "Optimal Chamber Adjustment Values: " << endl;
		//for(int j = 0; j < numChambersInRing; j++){
		//	if(printFlag !=0)cout << "Chamber " << j+1 << ": " << xs[j] << endl;
		//	if(printFlag !=0)cout << "Chamber error" << j+1 << ": " << error[j] << endl;
		//	runningSum += xs[j];	
		//	if(printFlag == 0) cout << xs[j] << ", " << endl;	
		//}
		//if(printFlag == 0) cout << "Now Printing Errors" << endl;
		//   for(int j = 0; j < numChambersInRing; j++){
		//            if(printFlag == 0) cout << error[j] << ", " << endl;
		//  }
		//if(printFlag !=0)cout << "Value of minimum: " << setprecision(15) << 1000000000000000*chisquare(xs) << endl; 
		//if(printFlag !=0)cout << "1-2: " << xs[0]-xs[1] << endl;
		//if(printFlag !=0)cout << "19-20: " << xs[18]-xs[19] << endl;
		//if(printFlag !=0)cout << "36-1: " << xs[35]-xs[0] << endl;
		//if(printFlag !=0)cout << "" << endl;
		//cout << "Average adjustment: "<< runningSum/36 << endl;
	}
	cout << "The true minimum at seedvalue " << setprecision(20) << pOfMin << " is: " << minValue << endl;


	//Minimizer
	ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );

	min.SetMaxFunctionCalls(1000000);
	min.SetMaxIterations(1000000);
	min.SetTolerance(0.001);
	min.SetPrintLevel(1);
	min.SetPrintLevel(1);
	/*
	   filename = Infilename;
	   station = InStation;
	   ring = InRing;
	   endcap = InEndcap;
	   InnumChambersInRing = numChambersInRing;  
	   */

	ROOT::Math::Functor f(&chisquare, numChambersInRing);

	vector<double> stepVector;
	vector<double> variableVector;

	for(int i =0; i < numChambersInRing; i++){
		stepVector.push_back(stepSize);
		variableVector.push_back(pOfMin);
		//variableVector.push_back(variableVal);
	}

	double* step = &stepVector[0];
	double* variable = &variableVector[0];

	min.SetFunction(f);
	min.SetVariable(0, "x0", variable[0], step[0]);
	/*{{{1*/	
	min.SetVariable(1, "x1", variable[1], step[1]);
	min.SetVariable(2, "x2", variable[2], step[2]);
	min.SetVariable(3, "x3", variable[3], step[3]);
	min.SetVariable(4, "x4", variable[4], step[4]);
	min.SetVariable(5, "x5", variable[5], step[5]);
	min.SetVariable(6, "x6", variable[6], step[6]);
	min.SetVariable(7, "x7", variable[7], step[7]);
	min.SetVariable(8, "x8", variable[8], step[8]);
	min.SetVariable(9, "x9", variable[9], step[9]);
	min.SetVariable(10, "x10", variable[10], step[10]);
	min.SetVariable(11, "x11", variable[11], step[11]);
	min.SetVariable(12, "x12", variable[12], step[12]);
	min.SetVariable(13, "x13", variable[13], step[13]);
	min.SetVariable(14, "x14", variable[14], step[14]);
	min.SetVariable(15, "x15", variable[15], step[15]);
	min.SetVariable(16, "x16", variable[16], step[16]);
	min.SetVariable(17, "x17", variable[17], step[17]);
	if(numChambersInRing != 18){	
		min.SetVariable(18, "x18", variable[18], step[18]);
		min.SetVariable(19, "x19", variable[19], step[19]);
		min.SetVariable(20, "x20", variable[20], step[20]);
		min.SetVariable(21, "x21", variable[21], step[21]);
		min.SetVariable(22, "x22", variable[22], step[22]);
		min.SetVariable(23, "x23", variable[23], step[23]);
		min.SetVariable(24, "x24", variable[24], step[24]);
		min.SetVariable(25, "x25", variable[25], step[25]);
		min.SetVariable(26, "x26", variable[26], step[26]);
		min.SetVariable(27, "x27", variable[27], step[27]);
		min.SetVariable(28, "x28", variable[28], step[28]);
		min.SetVariable(29, "x29", variable[29], step[29]);
		min.SetVariable(30, "x30", variable[30], step[30]);
		min.SetVariable(31, "x31", variable[31], step[31]);
		min.SetVariable(32, "x32", variable[32], step[32]);
		min.SetVariable(33, "x33", variable[33], step[33]);
		min.SetVariable(34, "x34", variable[34], step[34]);
		min.SetVariable(35, "x35", variable[35], step[35]);
	}	
	/*}}}*/	

	min.Minimize();
	//if(printFlag !=0) cout << min.Minimize() << endl;
	setprecision(6);
	const double *xs = min.X();
	const double *error = min.Errors();
	double runningSum = 0;
	cout << "Optimal Chamber Adjustment Values: " << endl;
	for(int j = 0; j < numChambersInRing; j++){
	      if(printFlag !=0)cout << "Chamber " << j+1 << ": " << xs[j] << endl;
	      if(printFlag !=0)cout << "Chamber error" << j+1 << ": " << error[j] << endl;
	      runningSum += xs[j];    
	      if(printFlag == 0) cout << setprecision(6) <<  xs[j] << ", " << endl;       
	}
	if(printFlag == 0) cout << "Now Printing Errors" << endl;
	   for(int j = 0; j < numChambersInRing; j++){
	            if(printFlag == 0) cout << setprecision(6) << error[j] << ", " << endl;
	  }
	if(printFlag !=0)cout << "Value of minimum: " << setprecision(6) << chisquare(xs) << endl; 
	if(printFlag !=0)cout << "1-2: " << xs[0]-xs[1] << endl;
	if(printFlag !=0)cout << "19-20: " << xs[18]-xs[19] << endl;
	if(printFlag !=0)cout << "36-1: " << xs[35]-xs[0] << endl;
	if(printFlag !=0)cout << "" << endl;
	cout << "Average adjustment: "<< runningSum/36 << endl;

	        cout << "Optimal Chamber Adjustment Values About 0: " << endl;
        for(int j = 0; j < numChambersInRing; j++){
              if(printFlag !=0)cout << "Chamber " << j+1 << ": " << xs[j] << endl;
              if(printFlag !=0)cout << "Chamber error" << j+1 << ": " << error[j] << endl;
              if(printFlag == 0) cout << setprecision(6) <<  xs[j]-(runningSum/36) << ", " << endl;
        }
        if(printFlag == 0) cout << "Now Printing Errors" << endl;
           for(int j = 0; j < numChambersInRing; j++){
                    if(printFlag == 0) cout << setprecision(6) << error[j] << ", " << endl;
          }

	minPlot->Draw();
	minPlot->SaveAs("minPlot.pdf");
}
