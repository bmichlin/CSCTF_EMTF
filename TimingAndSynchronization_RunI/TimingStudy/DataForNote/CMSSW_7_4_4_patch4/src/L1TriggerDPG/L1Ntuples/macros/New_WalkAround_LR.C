
#include "L1Ntuple.h"
#include<exception>
#include<TFile.h>
#include<TChain.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<TGraphErrors.h>
#include<TCanvas.h>
#include<TLegend.h>
#include<vector>
#include<iostream>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TLatex.h>
#include <iostream>
#include <fstream>
#include <TGraph2D.h>

// ---- CP error ----
#include <TGraphAsymmErrors.h>
#include <cmath>

#include "TEfficiency.h"//for Clopper Pearson
//--------

class WalkAround : public L1Ntuple
{
	public :
		//constructor
		WalkAround(std::string filename) : L1Ntuple(filename) {}
		WalkAround() {}
		~WalkAround() {}
		//main function macro : arguments can be adapted to needs
		void run(int runNumber, Long64_t numEvents, int EventSelection, int stationInput, int ringInput, int numChambersInRing, int ME2Trigger, int ME3Trigger, int forwardOnly, int PileUp50ns);
		//EventSelection = 0, all have precise topology, but precise to each step (TOF and Chi have different topologies).
		//EventSelection = 1, all have the same precise topology
		//forwardOnly selects only forward muons.  1= yes, 0/other = no
	private :

		//your private methods can be declared here
};

using namespace std;

// ------------------------------------------------------------------------------
// global variables
// ------------------------------------------------------------------------------
/*{{{*/
double Pi  = TMath::Pi();
double speedLight = 299.792; //in mm/ns. 
//Input Chamber parameters
//Members should be accessed using (Vector[Station][ring])
double chamberWidthBottom[5][5] = {{},{0,311,740,859},{0, 751, 895},{0, 835, 895},{0,903,895}};
double chamberWidthTop[5][5] = {{},{0,613,1078,1192},{0, 1534, 1530},{0, 1534, 1530},{0,1534,1530}};
double chamberPhi[5][5] = {{},{0,10,10,10},{0, 20, 10},{0, 20, 10},{0,20,10}};
double numWireGroups[5][5] = {{},{0,48,48,48},{0, 112, 64},{0, 96, 64},{0,96,64}};
double numHalfStrips[5][5] = {{},{0,224,160,128},{0, 160, 160},{0, 160, 160},{0,160,160}};
double stripWidthTop[5][5] = {{},{0,7.6,10.4,14.9},{0, 15.6, 16.0},{0, 15.6, 16.0},{0,15.6,16.0}};
double stripWidthBottom[5][5] = {{},{0,3.15,6.6,11.1},{0, 6.8, 8.5},{0, 7.8, 8.5},{0,8.6,8.5}};
double lengthOfStrip[5][5] = {{},{0,1680,1800,1900},{0, 2065, 3380},{0, 1845, 3380},{0,1665,3380}};
double stripPhi[5][5] = {{},{0,2.96,2.33,2.16},{0, 4.65, 2.33},{0, 4.65, 2.33},{0,4.65,2.33}};
/*}}}*/
// ------------------------------------------------------------------------------
// methods
// ------------------------------------------------------------------------------
//All supplementary methods 
/*{{{*/

double dist(double x1, double y1, double z1, double x2, double y2, double z2){
	double dist = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
	return dist;
}

bool contains(vector<int> list, int object){
	bool contain = false;
	for(int i =0; i < list.size(); i++){
		if(list[i] == object) contain = true;
	}
	return contain;
}

int findIndex(vector<int> list, int object){
	int index = -999;
	int indexCounter = 0;	
	if(contains(list,object)){
		for(int i = 0; i < list.size(); i++){
			if(list[i] == object){
				index = i;
				indexCounter++;
			}
		}
	}
	if(indexCounter > 1) ////cout << "The object appears more than once in the vector" << endl;
	return index;
}

bool vector_contains(vector < vector<int> > dupeList, int object){
	int dupeCounter = 0;
	for(int i = 0; i < dupeList.size(); i++){
		for(int j = 0; j < dupeList[i].size(); j++){
			if(object == dupeList[i][j] && object != dupeList[i][0]){
				dupeCounter++;			
			}
		}
	}
	if(dupeCounter > 0){
		return true;
	}else{
		return false;
	}
}

bool vecvec_contains(vector < vector<int> > vecvec, int object){
	bool out = false;
	for(int i = 0; i < vecvec.size(); i++){
		for(int j = 0; j < vecvec[i].size(); j++){
			if(object == vecvec[i][j])
				out = true;
		}
	}
	return out;
}

vector<int> neighborToChamber(vector < vector<int> > neighborList, int chamber){
	//If order = 1, then the input chamber was the low chamber and the ouput chamber is the high chamber.
	//If order = 2, then vice versa.	
	vector<int> info;
	int neighborChamber = 0;
	int order = 0;
	for(int i =0; i<neighborList.size(); i++){
		if(neighborList[i][0] == chamber){
			neighborChamber = neighborList[i][1];
			order = 1;
			info.push_back(neighborChamber);
			info.push_back(order);		
		}
		if(neighborList[i][1] == chamber){
			neighborChamber = neighborList[i][0];
			order = 2;	
			info.push_back(neighborChamber);
			info.push_back(order);
		}
	}
	return info;
}

void print_contents(vector<int> vec){
	int size = vec.size();
	//if(size == 0) //cout << "" << endl; 
	//if(size == 1) //cout << vec[0] << endl;
	//if(size == 2) //cout << vec[0] << ", " << vec[1] << endl;
	//if(size == 3) //cout << vec[0] << ", " << vec[1] << ", " << vec[2]  << endl;
	//if(size == 4) //cout << vec[0] << ", " << vec[1] << ", " << vec[2]  << ", " << vec[3] << endl;
	//if(size == 5) //cout << vec[0] << ", " << vec[1] << ", " << vec[2]  << ", " << vec[3] << ", " << vec[4] << endl;
	//if(size == 6) //cout << vec[0] << ", " << vec[1] << ", " << vec[2]  << ", " << vec[3] << ", " << vec[4] << ", " << vec[5]  << endl;

}


bool vectorvector_contains(vector < vector<int> > PairList, int object1, int object2){
	vector< int > forwards;
	forwards.push_back(object1);
	forwards.push_back(object2);
	vector< int > backwards;
	backwards.push_back(object2);
	backwards.push_back(object1);

	int counter = 0;

	for(int i = 0; i < PairList.size(); i++){
		if(PairList[i] == forwards || PairList[i] == backwards){
			counter++;
		}
	}
	if(counter > 0){
		return true;
	}else{
		return false;
	}
}

//The current version of this method uses the actual global phi and global eta.  To use the bits, uncomment the lines below.
double get_TOF(int globalEta1, int globalPhi1, int station1, int ring1, int sector1, int globalEta2, int globalPhi2, int station2, int ring2, int sector2){
	// for z position, I am using average station position from p204 of muTDR.  I am not worrying about even/odd stations.  ME11 is shifted compared to ME1/2 and ME1/3 and position is taken from p104 of muTDR.  No other station is treated differently.  Units of z are mm.

	////cout << "" << endl;
	////cout << "For ME1/1:  eta, phi, station, ring, sector: " << globalEta1 << ", " << globalPhi1 << ", " << station1 << ", " << ring1 << ", " << sector1 << endl;
	////cout << "For Other:  eta, phi, station, ring, sector: " << globalEta2 << ", " << globalPhi2 << ", " << station2 << ", " << ring2 << ", " << sector2 << endl;
	////cout << "" << endl;

	double z1, z2;
	if(station1 == 1 && ring1 == 1) z1 = 6000;
	if(station1 == 1 && ring1 != 1) z1 = 7070;
	if(station1 == 2) z1 = 8222.5;
	if(station1 == 3) z1 = 9407.5; //ALSO taken from p104
	if(station1 == 4) z1 = 10322.5;

	if(station2 == 1 && ring2 == 1) z2 = 6000;
	if(station2 == 1 && ring2 != 1) z2 = 7070;
	if(station2 == 2) z2 = 8222.5;
	if(station2 == 3) z2 = 9407.5; //ALSO taken from p104
	if(station2 == 4) z2 = 10322.5;

	//Now I need to go from eta and phi bits to a usable form
	//Each 62deg sector is split up in to 4096 equally spaced regions.
	//Note that for a track to be made, the LCTs are automatically in the same sector.
	//Eta goes form 0.9 to 2.4 in equally spaced bins of 126
	/*
	//Find phi (in radians) from bit
	int phiPoint1;
	if(sector1 == 1 || sector1 == 7) phiPoint1 = 15;
	if(sector1 == 2 || sector1 == 8) phiPoint1 = 75;
	if(sector1 == 3 || sector1 == 9) phiPoint1 = 135;
	if(sector1 == 4 || sector1 == 10) phiPoint1 = 195;
	if(sector1 == 5 || sector1 == 11) phiPoint1 = 255;
	if(sector1 == 6 || sector1 == 12) phiPoint1 = 315;
	int phiPoint2;
	if(sector2 == 1 || sector2 == 7) phiPoint2 = 15;
	if(sector2 == 2 || sector2 == 8) phiPoint2 = 75;
	if(sector2 == 3 || sector2 == 9) phiPoint2 = 135;
	if(sector2 == 4 || sector2 == 10) phiPoint2 = 195;
	if(sector2 == 5 || sector2 == 11) phiPoint2 = 255;
	if(sector2 == 6 || sector2 == 12) phiPoint2 = 315;

	double value1 = (double)62/4096;
	double value2 = (double)1/180;

	double phiAddition1 = value1*globalPhi1;
	double phiAddition2 = value1*globalPhi2;

	double phi1deg = phiAddition1+phiPoint1;
	double phi2deg = phiAddition2+phiPoint2;	
	//now to radians
	double phi1 = phi1deg*Pi*(value2);
	double phi2 = phi2deg*Pi*(value2);


	if(phi1 > (2*Pi)) phi1 = phi1-(2*Pi);
	if(phi2 > (2*Pi)) phi2 = phi2-(2*Pi);
	////cout << "Phi1 (phiPoint, phiAddition, phideg, phi, z): " << phiPoint1 << ", " << phiAddition1 << ", " << phi1deg << ", " << phi1 << ", " << z1 << endl;


	//Find eta from bit

	double eta1 = (1.5/126)*globalEta1;
	double eta2 = (1.5/126)*globalEta2;
	*/
	//Find theta from eta
	double theta1 = 2*atan(exp(-1*globalEta1));
	double theta2 = 2*atan(exp(-1*globalEta2));

	//Find r from z and theta
	double r1 = (z1)*(1/cos(theta1));
	double r2 = (z2)*(1/cos(theta2));

	//Find x and y from r, theta, phi
	double x1 = r1*sin(theta1)*cos(globalPhi1);
	double x2 = r2*sin(theta2)*cos(globalPhi2);
	double y1 = r1*sin(theta1)*sin(globalPhi1);
	double y2 = r2*sin(theta2)*sin(globalPhi2);
	////cout << "For ME1/1 (eta, theta, phi, r): " << globalEta1 << ", " << theta1 << ", " << globalPhi1 << ", " << r1 << endl;
	////cout << "For Other (eta, theta, phi, r): " << globalEta2 << ", " << theta2 << ", " << globalPhi2 << ", " << r2 << endl;
	////cout << "" << endl;
	//Find the differences
	double dx = x2-x1;
	double dy = y2-y1;
	double dz = z2-z1;
	////cout << "For ME1/1 (x,y,z): " << x1 << ", " << y1 << ", " << z1 << endl;
	////cout << "For Other (x,y,z): " << x2 << ", " << y2 << ", " << z2 << endl;
	//Square the differences
	double dxSqr = dx*dx;
	double dySqr = dy*dy;
	double dzSqr = dz*dz;

	//Take sqrt of the argument
	double arg = dxSqr+dySqr+dzSqr;	
	double distance = sqrt(arg);

	//TOF=d/v.  NOTE: distance is in mm. c is in mm/ns.  BX is in ns
	double tof = distance/speedLight;
	////cout << "distance, tof: " << distance << ", " << tof << endl;
	return tof;
}

double get_TOF_Bx_Correction(double tof, int bx1, int bx){
	//The TOF correction is rather coarse, broken up by 12.5ns.  
	//ie: if tof <12.5 ns, then there is no correction. if 12.5 < tof < 25 then it moves by one bin. etc.

	int correction = (int) tof/12.5; //This works.  See test() below.
	double corBX1, corBX;
	if(bx1 > bx){
		corBX1 = bx1-correction;
		corBX = bx;
	}else{
		corBX = bx-correction;
		corBX1 = bx1; 
	}	
	double correctedBX = corBX1-corBX;
	return correctedBX;

}

int numLCTs_inTrack(L1Analysis::L1AnalysisCSCTFDataFormat *csctf_, int trk){
	int lctCounter = 0;
	bool go = true;
	for(int i = 0; go == true; i++){
		if(csctf_->trLctChamber[trk][i] > 0){
			lctCounter++;
		}else{
			go = false;
		}  
	}
	return lctCounter;
}


/*double getz(int globalEta1, int globalPhi1, int station1, int ring1, int endcap){
// for z position, I am using average station position from p204 of muTDR.  I am not worrying about even/odd stations.  ME11 is shifted compared to ME1/2 and ME1/3 and position is taken from p104 of muTDR.  No other station is treated differently.  Units of z are mm.

double z1, z2;
if(station1 == 1 && ring1 == 1) z1 = 6000;
if(station1 == 1 && ring1 != 1) z1 = 7070;
if(station1 == 2) z1 = 8222.5;
if(station1 == 3) z1 = 9407.5; //ALSO taken from p104
if(station1 == 4) z1 = 10322.5;

if(endcap>1) z1=-z1;
return z1;
}
*/

double getz(int globalEta1, int globalPhi1, int station1, int ring1, int endcap, int chamber){
	// for z position, I am using average station position from p204 of muTDR.  I am not worrying about even/odd stations.  
	//ME11 is shifted compared to ME1/2 and ME1/3 and position is taken from p104 of muTDR.  
	//No other station is treated differently.  Units of z are mm. 

	double z1, z2;
	if(endcap==1){
		if(station1==1 && ring1==1 && chamber%2==1) return 6152;
		else if(station1==1 && ring1==1 && chamber%2==0) return 5861;
		else if(station1==1 && ring1==2 && chamber%2==1) return 7115;
		else if(station1==1 && ring1==2 && chamber%2==0) return 6844;
		else if(station1==1 && ring1==3) return 6940;
		else if(station1==2 && chamber%2==1) return 8398;
		else if(station1==2 && chamber%2==0) return 8150;
		else if(station1==3 && chamber%2==1) return 9219;
		else if(station1==3 && chamber%2==0) return 9467;
		else if(station1==4 && chamber%2==1) return 10127;
		else if(station1==4 && chamber%2==0) return 10374;
	}
	else{
		if(station1==1 && ring1==1 && chamber%2==1) return -6153;
		else if(station1==1 && ring1==1 && chamber%2==0) return -5861;
		else if(station1==1 && ring1==2 && chamber%2==1) return -7116;
		else if(station1==1 && ring1==2 && chamber%2==0) return -6845;
		else if(station1==1 && ring1==3) return -6942;
		else if(station1==2 && chamber%2==1) return -8405;
		else if(station1==2 && chamber%2==0) return -8157;
		else if(station1==3 && chamber%2==1) return -9231;
		else if(station1==3 && chamber%2==0) return -9478;
		else if(station1==4 && chamber%2==1) return -10126;
		else if(station1==4 && chamber%2==0) return -10373;
	}
	return 0;
}


double gety(int globalEta1, int globalPhi1, int station1, int ring1, int sector1){
	// for z position, I am using average station position from p204 of muTDR.  I am not worrying about even/odd stations.  ME11 is shifted compared to ME1/2 and ME1/3 and position is taken from p104 of muTDR.  No other station is treated differently.  Units of z are mm.


	int phiPoint1;
	if(sector1 == 1 || sector1 == 7) phiPoint1 = 15;
	if(sector1 == 2 || sector1 == 8) phiPoint1 = 75;
	if(sector1 == 3 || sector1 == 9) phiPoint1 = 135;
	if(sector1 == 4 || sector1 == 10) phiPoint1 = 195;
	if(sector1 == 5 || sector1 == 11) phiPoint1 = 255;
	if(sector1 == 6 || sector1 == 12) phiPoint1 = 315;


	double value1 = (double)62/4096;
	double value2 = (double)1/180;

	double phiAddition1 = value1*globalPhi1;

	double phi1deg = phiAddition1+phiPoint1;
	//now to radians
	double phi1 = phi1deg*TMath::Pi()*(value2);


	if(phi1 > (2*TMath::Pi())) phi1 = phi1-(2*TMath::Pi());


	//Find eta from bit

	double eta1 = (1.5/126)*globalEta1+0.9;


	double z1, z2;
	if(station1 == 1 && ring1 == 1) z1 = 6000;
	if(station1 == 1 && ring1 != 1) z1 = 7070;
	if(station1 == 2) z1 = 8222.5;
	if(station1 == 3) z1 = 9407.5; //ALSO taken from p104
	if(station1 == 4) z1 = 10322.5;

	double theta1 = 2*atan(exp(-1*eta1));

	//Find r from z and theta
	double r1 = (z1)*(1/cos(theta1));

	//Find x and y from r, theta, phi
	double x1 = r1*sin(theta1)*cos(phi1);
	double y1 = r1*sin(theta1)*sin(phi1);
	return y1;
}

double getx(int globalEta1, int globalPhi1, int station1, int ring1, int sector1){
	// for z position, I am using average station position from p204 of muTDR.  I am not worrying about even/odd stations.  ME11 is shifted compared to ME1/2 and ME1/3 and position is taken from p104 of muTDR.  No other station is treated differently.  Units of z are mm.

	//Now I need to go from eta and phi bits to a usable form
	//Each 62deg sector is split up in to 4096 equally spaced regions.
	//Note that for a track to be made, the LCTs are automatically in the same sector.
	//Eta goes form 0.9 to 2.4 in equally spaced bins of 126
	//Find phi (in radians) from bit
	int phiPoint1;
	if(sector1 == 1 || sector1 == 7) phiPoint1 = 15;
	if(sector1 == 2 || sector1 == 8) phiPoint1 = 75;
	if(sector1 == 3 || sector1 == 9) phiPoint1 = 135;
	if(sector1 == 4 || sector1 == 10) phiPoint1 = 195;
	if(sector1 == 5 || sector1 == 11) phiPoint1 = 255;
	if(sector1 == 6 || sector1 == 12) phiPoint1 = 315;


	double value1 = (double)62/4096;
	double value2 = (double)1/180;

	double phiAddition1 = value1*globalPhi1;

	double phi1deg = phiAddition1+phiPoint1;
	//now to radians
	double phi1 = phi1deg*TMath::Pi()*(value2);


	if(phi1 > (2*TMath::Pi())) phi1 = phi1-(2*TMath::Pi());


	//Find eta from bit

	double eta1 = (1.5/126)*globalEta1+0.9;



	double z1, z2;
	if(station1 == 1 && ring1 == 1) z1 = 6000;
	if(station1 == 1 && ring1 != 1) z1 = 7070;
	if(station1 == 2) z1 = 8222.5;
	if(station1 == 3) z1 = 9407.5; //ALSO taken from p104
	if(station1 == 4) z1 = 10322.5;

	double theta1 = 2*atan(exp(-1*eta1));

	//Find r from z and theta
	double r1 = (z1)*(1/cos(theta1));

	//Find x and y from r, theta, phi
	//      std::////cout<<"z1 "<<z1<<" eta1 "<<eta1<<" theta1 "<<theta1<<" r1 "<<r1<<" radius "<<r1*sin(theta1)<<std::endl;
	double x1 = r1*sin(theta1)*cos(phi1);
	double y1 = r1*sin(theta1)*sin(phi1);
	//      std::////cout<<x1<<std::endl;
	return x1;
}
double degToRad(double degree){
	double rad = (degree*Pi)/180;
	return rad;
}

double signalPropTime_BAD(int ring, int station, int strip, int wire){

	double yposA = numWireGroups[station][ring]*(wire-1);
	double yposB = (chamberWidthTop[station][ring]-chamberWidthBottom[station][ring])/(2*tan(degToRad(chamberPhi[station][ring])/2)); 
	double ypos = yposA/yposB;

	double w1 = chamberWidthBottom[station][ring]/(2*tan(degToRad(stripPhi[station][ring])/2));

	double stripThickness = (w1+ypos)*tan(degToRad(stripPhi[station][ring]));

	double lctStrip = strip;
	if(station == 1 && ring == 1 && strip >= 128){
		//Then ME11A
		lctStrip = 96 - strip;
	} //ME11B is normal

	double signalDistance = stripThickness*(lctStrip-1); //in mm

	double signalTime = signalDistance/speedLight; // in ns

	double BXtime = signalTime/25; //This gives the propogation time in BX units.
	cout << "" << endl;
	cout << "numWireGroups, wire, yposA                            : " << numWireGroups[station][ring] << ", " << wire << ", " << yposA << ", " << endl;
	cout << "chamberWidthTop, chamberWidthBottom, chamberPhi, yposB: " << chamberWidthTop[station][ring] << ", " << chamberWidthBottom[station][ring] << ", " << chamberPhi[station][ring] << ", " << yposB << endl;
	cout << "yposA/yposB, ypos                                     : " << yposA/yposB << ", " << ypos << endl;
	cout << "chamberWidthBottom, stripPhi, W1                      : " << chamberWidthBottom[station][ring] << ", " << stripPhi[station][ring] << ", " << w1 << endl;
	cout << "W1, ypos, stripPhi, stripThickness                    : " << w1 << ", " << ypos << ", " << stripPhi[station][ring] << ", " << stripThickness << endl;
	cout << "strip, lctStrip, station, ring                        : " << strip << ", " << lctStrip << ", " << station << ", " << ring << endl;
	cout << "" << endl;
	cout << "Ypos, W1, S_T, signal Distance, signalTime, BXtime    : " << ypos << ", " << w1 << ", " << stripThickness << ", " << signalDistance << ", " << signalTime << ", " << BXtime << endl;

	return BXtime;
}
double signalPropTime(int ring, int station, int strip, int wire){
	double height = (chamberWidthTop[station][ring]-chamberWidthBottom[station][ring])/(2*tan(degToRad(chamberPhi[station][ring])/2));
	double h1 = (numWireGroups[station][ring]/height)*(wire+1);
	double w1 = 2*((h1+chamberWidthBottom[station][ring]/(2*tan(degToRad(chamberPhi[station][ring])/2))))*tan(degToRad(chamberPhi[station][ring])/2);
	double st = w1/numHalfStrips[station][ring];
	double lctStrip = strip;
	double sp = 999;
	if(station == 1 && ring == 1 && strip >= 128){
		//Then ME11A
		//sp = st*(128-lctStrip); //This is for AFEB on the right	
		sp=st*(lctStrip-128); //This is for AFEB on the left
	}else{ //ME11B is normal
		//sp = st*(numHalfStrips[station][ring]-lctStrip); //AFEB on right
		sp=st*(lctStrip + 1);	//AFEB on left
	}	
	double tbx = (sp)/(25*speedLight);

	//cout << "Signal Distance        : " << sp << endl;
	//cout << "Signal Propagation Time: " << sp/speedLight << endl;
	//cout << "BX Time                : " << tbx << endl;

	return tbx;
}

double signalPropTimePrint(int ring, int station, int strip, int wire){
	double height = (chamberWidthTop[station][ring]-chamberWidthBottom[station][ring])/(2*tan(degToRad(chamberPhi[station][ring])/2));
	double h1 = (numWireGroups[station][ring]/height)*(wire+1);
	double w1 = 2*((h1+chamberWidthBottom[station][ring]/(2*tan(degToRad(chamberPhi[station][ring])/2))))*tan(degToRad(chamberPhi[station][ring])/2);
	double st = w1/numHalfStrips[station][ring];
	double lctStrip = strip;
	double sp = 999;

	cout << "" << endl;
	if(station == 1 && ring == 1 && strip >= 128){
		//Then ME11A
		sp = st*(128-lctStrip); //This is for AFEB on the right 
		//sp=st*(lctStrip-128); //This is for AFEB on the left
	}else{ //ME11B is normal
		sp = st*(numHalfStrips[station][ring]-lctStrip); //AFEB on right
		//sp=st*(lctStrip + 1); //AFEB on left
	}
	double tbx = (sp)/(25*speedLight);

	cout << "Signal Distance        : " << sp << endl;
	cout << "Signal Propagation Time: " << sp/speedLight << endl;
	cout << "BX Time                : " << tbx << endl;

	return tbx;
}

/*}}}*/
void WalkAround::run(int runNumber, Long64_t numEvents, int EventSelection, int stationInput, int ringInput, int numChambersInRing, int ME2Trigger, int ME3Trigger, int forwardOnly, int PileUp50ns) 
{

	////cout << "runNumber: " << runNumber << endl;	
	TFile *BAM = new TFile(Form("output_WalkAround_ME%d%d_%d_%d_Trig%d%d.root",stationInput,ringInput,runNumber,EventSelection,ME2Trigger,ME3Trigger),"RECREATE");

	//Plot style settings
	/*{{{*/
	gROOT->Clear();
	gStyle->SetOptStat(210);  
	gStyle->SetHistFillColor(4);
	gStyle->SetHistFillStyle(3004);
	/*}}}*/
	//--------------------------------------------------------------------------
	// Declare Plots 
	//--------------------------------------------------------------------------

	////cout << "Begin Creating Plots" << endl;
	/*{{{*/
	//Checkpoint	
	TDirectory *Misc = BAM->mkdir("Misc");	
	Misc->cd();

	TH1* NumLCTs					= new TH1D("NumLCTs", "NumLCTs", 10, 0, 10);	

	TObjArray HListWire(0);
	int wireStationCounter = 0;
	while(wireStationCounter < 4){
		wireStationCounter++;
		int wireRingCounter = 0;
		while(wireRingCounter < 3){
			wireRingCounter++;		
			TString histoWire = TString::Format("MEplus%d%d_WireGroup",wireStationCounter,wireRingCounter);
			TH2F *myhistWire = ((TH2F *)(HListWire.FindObject(histoWire)));
			if(!myhistWire){
				myhistWire = new TH2F(histoWire,histoWire, 36,0,36,150,0,150);
				HListWire.Add(myhistWire);
			}
			TString histoWireNeg = TString::Format("MEneg%d%d_WireGroup",wireStationCounter,wireRingCounter);
			TH2F *myhistWireNeg = ((TH2F *)(HListWire.FindObject(histoWireNeg)));
			if(!myhistWireNeg){
				myhistWireNeg = new TH2F(histoWireNeg,histoWireNeg, 36,0,36,150,0,150);
				HListWire.Add(myhistWireNeg);
			}
		}
	}

	TDirectory *PositiveEndcap = BAM->mkdir("PositiveEndcap");
	PositiveEndcap->cd();

	TH1* ME11_pos_BX                                  = new TH1D(Form("ME11_pos_BX_%d",runNumber),"ME11_pos_BX",7,3,10);
	TH1* ME12_pos_BX                                  = new TH1D(Form("ME12_pos_BX_%d",runNumber),"ME12_pos_BX",7,3,10);
	TH1* ME13_pos_BX                                  = new TH1D(Form("ME13_pos_BX_%d",runNumber),"ME13_pos_BX",7,3,10);
	TH1* ME21_pos_BX                                  = new TH1D(Form("ME21_pos_BX_%d",runNumber),"ME21_pos_BX",7,3,10);
	TH1* ME22_pos_BX                                  = new TH1D(Form("ME22_pos_BX_%d",runNumber),"ME22_pos_BX",7,3,10);
	TH1* ME31_pos_BX                                  = new TH1D(Form("ME31_pos_BX_%d",runNumber),"ME31_pos_BX",7,3,10);
	TH1* ME32_pos_BX                                  = new TH1D(Form("ME32_pos_BX_%d",runNumber),"ME32_pos_BX",7,3,10);
	TH1* ME41_pos_BX                                  = new TH1D(Form("ME41_pos_BX_%d",runNumber),"ME41_pos_BX",7,3,10);
	TH1* ME42_pos_BX                                  = new TH1D(Form("ME42_pos_BX_%d",runNumber),"ME42_pos_BX",7,3,10);

	TDirectory *NegativeEndcap = BAM->mkdir("NegativeEndcap");
	NegativeEndcap->cd();

	TH1* ME11_neg_BX                                  = new TH1D(Form("ME11_neg_BX_%d",runNumber),"ME11_neg_BX",7,3,10);
	TH1* ME12_neg_BX                                  = new TH1D(Form("ME12_neg_BX_%d",runNumber),"ME12_neg_BX",7,3,10);
	TH1* ME13_neg_BX                                  = new TH1D(Form("ME13_neg_BX_%d",runNumber),"ME13_neg_BX",7,3,10);
	TH1* ME21_neg_BX                                  = new TH1D(Form("ME21_neg_BX_%d",runNumber),"ME21_neg_BX",7,3,10);
	TH1* ME22_neg_BX                                  = new TH1D(Form("ME22_neg_BX_%d",runNumber),"ME22_neg_BX",7,3,10);
	TH1* ME31_neg_BX                                  = new TH1D(Form("ME31_neg_BX_%d",runNumber),"ME31_neg_BX",7,3,10);
	TH1* ME32_neg_BX                                  = new TH1D(Form("ME32_neg_BX_%d",runNumber),"ME32_neg_BX",7,3,10);
	TH1* ME41_neg_BX                                  = new TH1D(Form("ME41_neg_BX_%d",runNumber),"ME41_neg_BX",7,3,10);
	TH1* ME42_neg_BX                                  = new TH1D(Form("ME42_neg_BX_%d",runNumber),"ME42_neg_BX",7,3,10);

	TDirectory *Strip = BAM->mkdir("Strip");
	Strip->cd();
	TObjArray HListStrip(0);
	int chamberCounterStrip = 0;
	while(chamberCounterStrip < 36){
		chamberCounterStrip++;
		TString histonameStrip = TString::Format("Strip_BX_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterStrip);
		TH2F *myhistStrip = ((TH2F *)(HListStrip.FindObject(histonameStrip)));
		if(!myhistStrip){
			myhistStrip = new TH2F(histonameStrip, histonameStrip,225,0,225,7,3,10);
			HListStrip.Add(myhistStrip);
		}
		TString histonameStrip2 = TString::Format("Strip_BX_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterStrip);
		TH2F *myhistStrip2 = ((TH2F *)(HListStrip.FindObject(histonameStrip2)));
		if(!myhistStrip2){
			myhistStrip2 = new TH2F(histonameStrip2, histonameStrip2,225,0,225,7,3,10);
			HListStrip.Add(myhistStrip2);
		}

	}

	TDirectory *Strip_TOF = BAM->mkdir("Strip_TOF");
	Strip_TOF->cd();
	TObjArray HListStrip_TOF(0);
	int chamberCounterStrip_TOF = 0;
	while(chamberCounterStrip_TOF < 36){
		chamberCounterStrip_TOF++;
		TString histonameStrip_TOF = TString::Format("Strip_TOF_BX_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterStrip_TOF);
		TH2F *myhistStrip_TOF = ((TH2F *)(HListStrip_TOF.FindObject(histonameStrip_TOF)));
		if(!myhistStrip_TOF){
			myhistStrip_TOF = new TH2F(histonameStrip_TOF, histonameStrip_TOF,225,0,225,700,3,10);
			HListStrip_TOF.Add(myhistStrip_TOF);
		}
		TString histonameStrip_TOF2 = TString::Format("Strip_TOF_BX_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterStrip_TOF);
		TH2F *myhistStrip_TOF2 = ((TH2F *)(HListStrip_TOF.FindObject(histonameStrip_TOF2)));
		if(!myhistStrip_TOF2){
			myhistStrip_TOF2 = new TH2F(histonameStrip_TOF2, histonameStrip_TOF2,225,0,225,700,3,10);
			HListStrip_TOF.Add(myhistStrip_TOF2);
		}

	}


	TDirectory *NeighborStripTOF_Rings_ME11 = BAM->mkdir("NeighborStripTOF_Rings_ME11");
	TDirectory *NeighborStripTOF_Rings_ME12 = BAM->mkdir("NeighborStripTOF_Rings_ME12");
	TDirectory *NeighborStripTOF_Rings_ME13 = BAM->mkdir("NeighborStripTOF_Rings_ME13");
	TDirectory *NeighborStripTOF_Rings_ME21 = BAM->mkdir("NeighborStripTOF_Rings_ME21");
	TDirectory *NeighborStripTOF_Rings_ME22 = BAM->mkdir("NeighborStripTOF_Rings_ME22");
	TDirectory *NeighborStripTOF_Rings_ME31 = BAM->mkdir("NeighborStripTOF_Rings_ME31");
	TDirectory *NeighborStripTOF_Rings_ME32 = BAM->mkdir("NeighborStripTOF_Rings_ME32");
	TDirectory *NeighborStripTOF_Rings_ME41 = BAM->mkdir("NeighborStripTOF_Rings_ME41");
	TDirectory *NeighborStripTOF_Rings_ME42 = BAM->mkdir("NeighborStripTOF_Rings_ME42");

	TObjArray NeighborStripTOF(0);
	for(int endcap = 0; endcap < 2; endcap++){
		int endcapActual = 0;
		if(endcap == 0) endcapActual = 1;
		if(endcap == 1) endcapActual = -1;
		if(endcapActual == 0) continue;
		int stationPlotCounterWA = 0;
		while(stationPlotCounterWA < 4){
			stationPlotCounterWA++;
			int ringPlotCounterWA = 0;
			while(ringPlotCounterWA < 3){
				ringPlotCounterWA++;
				int chamberPlotCounterWA = 0;
				int chamberPlotCounterAdditionWA = 0;
				if(stationPlotCounterWA == 1 || ringPlotCounterWA == 2){
					while(chamberPlotCounterWA < 36){
						chamberPlotCounterWA++;

						if(chamberPlotCounterWA != 36)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 36)  chamberPlotCounterAdditionWA = 1;

						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 1) NeighborStripTOF_Rings_ME11->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 2) NeighborStripTOF_Rings_ME12->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 3) NeighborStripTOF_Rings_ME13->cd();
						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 2) NeighborStripTOF_Rings_ME22->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 2) NeighborStripTOF_Rings_ME32->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 2) NeighborStripTOF_Rings_ME42->cd();

						TString histoname = TString::Format("NeighborStripTOF_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH2F *myhist = ((TH2F *)(NeighborStripTOF.FindObject(histoname)));
						if(!myhist){
							myhist = new TH2F(histoname, histoname, 225,0,225,700,3,10);
							NeighborStripTOF.Add(myhist);
						}
						TString histoname_ChamberA = TString::Format("NeighborStripTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH2F *myhist_ChamberA = ((TH2F *)(NeighborStripTOF.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH2F(histoname_ChamberA, histoname_ChamberA, 225,0,225,700,3,10);
							NeighborStripTOF.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("NeighborStripTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH2F *myhist_ChamberB = ((TH2F *)(NeighborStripTOF.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH2F(histoname_ChamberB, histoname_ChamberB, 225,0,225,700,3,10);
							NeighborStripTOF.Add(myhist_ChamberB);
						}
					}
				}else{
					while(chamberPlotCounterWA < 18){
						chamberPlotCounterWA++;
						if(chamberPlotCounterWA != 18)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 18)  chamberPlotCounterAdditionWA = 1;
						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 1) NeighborStripTOF_Rings_ME21->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 1) NeighborStripTOF_Rings_ME31->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 1) NeighborStripTOF_Rings_ME41->cd();

						TString histoname = TString::Format("NeighborStripTOF_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA,endcapActual);
						TH2F *myhist = ((TH2F *)(NeighborStripTOF.FindObject(histoname)));
						TString histoname_ChamberA = TString::Format("NeighborStripTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH2F *myhist_ChamberA = ((TH2F *)(NeighborStripTOF.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH2F(histoname_ChamberA, histoname_ChamberA, 225,0,225,700,3,10);
							NeighborStripTOF.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("NeighborStripTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH2F *myhist_ChamberB = ((TH2F *)(NeighborStripTOF.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH2F(histoname_ChamberB, histoname_ChamberB, 225,0,225,700,3,10);
							NeighborStripTOF.Add(myhist_ChamberB);
						}
						if(!myhist){
							myhist = new TH2F(histoname, histoname, 225,0,225,700,3,10);
							NeighborStripTOF.Add(myhist);
						}
					}
				}
			}
		}
	}

	TDirectory *NeighborStripTOF_And_Prop_Rings_ME11 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME11");
	TDirectory *NeighborStripTOF_And_Prop_Rings_ME12 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME12");
	TDirectory *NeighborStripTOF_And_Prop_Rings_ME13 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME13");
	TDirectory *NeighborStripTOF_And_Prop_Rings_ME21 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME21");
	TDirectory *NeighborStripTOF_And_Prop_Rings_ME22 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME22");
	TDirectory *NeighborStripTOF_And_Prop_Rings_ME31 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME31");
	TDirectory *NeighborStripTOF_And_Prop_Rings_ME32 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME32");
	TDirectory *NeighborStripTOF_And_Prop_Rings_ME41 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME41");
	TDirectory *NeighborStripTOF_And_Prop_Rings_ME42 = BAM->mkdir("NeighborStripTOF_And_Prop_Rings_ME42");

	TObjArray NeighborStripTOF_And_Prop(0);
	for(int endcap = 0; endcap < 2; endcap++){
		int endcapActual = 0;
		if(endcap == 0) endcapActual = 1;
		if(endcap == 1) endcapActual = -1;
		if(endcapActual == 0) continue;
		int stationPlotCounterWA = 0;
		while(stationPlotCounterWA < 4){
			stationPlotCounterWA++;
			int ringPlotCounterWA = 0;
			while(ringPlotCounterWA < 3){
				ringPlotCounterWA++;
				int chamberPlotCounterWA = 0;
				int chamberPlotCounterAdditionWA = 0;
				if(stationPlotCounterWA == 1 || ringPlotCounterWA == 2){
					while(chamberPlotCounterWA < 36){
						chamberPlotCounterWA++;

						if(chamberPlotCounterWA != 36)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 36)  chamberPlotCounterAdditionWA = 1;

						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 1) NeighborStripTOF_And_Prop_Rings_ME11->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 2) NeighborStripTOF_And_Prop_Rings_ME12->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 3) NeighborStripTOF_And_Prop_Rings_ME13->cd();
						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 2) NeighborStripTOF_And_Prop_Rings_ME22->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 2) NeighborStripTOF_And_Prop_Rings_ME32->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 2) NeighborStripTOF_And_Prop_Rings_ME42->cd();

						TString histoname = TString::Format("NeighborStripTOF_And_Prop_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH2F *myhist = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(histoname)));
						if(!myhist){
							myhist = new TH2F(histoname, histoname, 225,0,225,700,3,10);
							NeighborStripTOF_And_Prop.Add(myhist);
						}
						TString histoname_ChamberA = TString::Format("NeighborStripTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH2F *myhist_ChamberA = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH2F(histoname_ChamberA, histoname_ChamberA, 225,0,225,700,3,10);
							NeighborStripTOF_And_Prop.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("NeighborStripTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH2F *myhist_ChamberB = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH2F(histoname_ChamberB, histoname_ChamberB, 225,0,225,700,3,10);
							NeighborStripTOF_And_Prop.Add(myhist_ChamberB);
						}
					}
				}else{
					while(chamberPlotCounterWA < 18){
						chamberPlotCounterWA++;
						if(chamberPlotCounterWA != 18)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 18)  chamberPlotCounterAdditionWA = 1;
						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 1) NeighborStripTOF_And_Prop_Rings_ME21->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 1) NeighborStripTOF_And_Prop_Rings_ME31->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 1) NeighborStripTOF_And_Prop_Rings_ME41->cd();

						TString histoname = TString::Format("NeighborStripTOF_And_Prop_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA,endcapActual);
						TH2F *myhist = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(histoname)));
						TString histoname_ChamberA = TString::Format("NeighborStripTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH2F *myhist_ChamberA = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH2F(histoname_ChamberA, histoname_ChamberA, 225,0,225,700,3,10);
							NeighborStripTOF_And_Prop.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("NeighborStripTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH2F *myhist_ChamberB = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH2F(histoname_ChamberB, histoname_ChamberB, 225,0,225,700,3,10);
							NeighborStripTOF_And_Prop.Add(myhist_ChamberB);
						}
						if(!myhist){
							myhist = new TH2F(histoname, histoname, 225,0,225,700,3,10);
							NeighborStripTOF_And_Prop.Add(myhist);
						}
					}
				}
			}
		}
	}




	TDirectory *TOF = BAM->mkdir("TOF");
	TOF->cd();
	TObjArray HListTOF(0);
	int chamberCounterTOF = 0;
	while(chamberCounterTOF < 36){
		chamberCounterTOF++;
		TString histonameTOF = TString::Format("TOF_BX_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterTOF);
		TH1F *myhistTOF = ((TH1F *)(HListTOF.FindObject(histonameTOF)));
		if(!myhistTOF){
			myhistTOF = new TH1F(histonameTOF, histonameTOF, 700,3,10);
			HListTOF.Add(myhistTOF);
		}
		TString histonameTOF2 = TString::Format("TOF_BX_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterTOF);
		TH1F *myhistTOF2 = ((TH1F *)(HListTOF.FindObject(histonameTOF2)));
		if(!myhistTOF2){
			myhistTOF2 = new TH1F(histonameTOF2, histonameTOF2, 700,3,10);
			HListTOF.Add(myhistTOF2);
		}

	}



	TDirectory *TOFforNeighbors = BAM->mkdir("TOFforNeighbors");
	TOFforNeighbors->cd();
	TObjArray HListTOFforNeighbors(0);
	int chamberCounterTOFforNeighbors = 0;
	while(chamberCounterTOFforNeighbors < 36){
		chamberCounterTOFforNeighbors++;
		TString histonameTOFforNeighbors = TString::Format("TOFforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterTOFforNeighbors);
		TH1F *myhistTOFforNeighbors = ((TH1F *)(HListTOFforNeighbors.FindObject(histonameTOFforNeighbors)));
		if(!myhistTOFforNeighbors){
			myhistTOFforNeighbors = new TH1F(histonameTOFforNeighbors, histonameTOFforNeighbors, 700,3,10);
			HListTOFforNeighbors.Add(myhistTOFforNeighbors);
		}
		TString histonameTOFforNeighbors2 = TString::Format("TOFforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterTOFforNeighbors);
		TH1F *myhistTOFforNeighbors2 = ((TH1F *)(HListTOFforNeighbors.FindObject(histonameTOFforNeighbors2)));
		if(!myhistTOFforNeighbors2){
			myhistTOFforNeighbors2 = new TH1F(histonameTOFforNeighbors2, histonameTOFforNeighbors2, 700,3,10);
			HListTOFforNeighbors.Add(myhistTOFforNeighbors2);
		}

	}

	TDirectory *TOF_And_PropforNeighbors = BAM->mkdir("TOF_And_PropforNeighbors");
	TOF_And_PropforNeighbors->cd();
	TObjArray HListTOF_And_PropforNeighbors(0);
	int chamberCounterTOF_And_PropforNeighbors = 0;
	while(chamberCounterTOF_And_PropforNeighbors < 36){
		chamberCounterTOF_And_PropforNeighbors++;
		TString histonameTOF_And_PropforNeighbors = TString::Format("TOF_And_PropforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterTOF_And_PropforNeighbors);
		TH1F *myhistTOF_And_PropforNeighbors = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(histonameTOF_And_PropforNeighbors)));
		if(!myhistTOF_And_PropforNeighbors){
			myhistTOF_And_PropforNeighbors = new TH1F(histonameTOF_And_PropforNeighbors, histonameTOF_And_PropforNeighbors, 700,3,10);
			HListTOF_And_PropforNeighbors.Add(myhistTOF_And_PropforNeighbors);
		}
		TString histonameTOF_And_PropforNeighbors2 = TString::Format("TOF_And_PropforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterTOF_And_PropforNeighbors);
		TH1F *myhistTOF_And_PropforNeighbors2 = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(histonameTOF_And_PropforNeighbors2)));
		if(!myhistTOF_And_PropforNeighbors2){
			myhistTOF_And_PropforNeighbors2 = new TH1F(histonameTOF_And_PropforNeighbors2, histonameTOF_And_PropforNeighbors2, 700,3,10);
			HListTOF_And_PropforNeighbors.Add(myhistTOF_And_PropforNeighbors2);
		}

	}


	TDirectory *TOFCorrection = BAM->mkdir("TOFCorrection");
	TOFCorrection->cd();
	TObjArray HListTOFCorrection(0);
	int chamberCounterTOFCorrection = 0;
	while(chamberCounterTOFCorrection < 36){
		chamberCounterTOFCorrection++;
		TString histonameTOFCorrection = TString::Format("TOFCorrection_BX_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterTOFCorrection);
		TH1F *myhistTOFCorrection = ((TH1F *)(HListTOFCorrection.FindObject(histonameTOFCorrection)));
		if(!myhistTOFCorrection){
			myhistTOFCorrection = new TH1F(histonameTOFCorrection, histonameTOFCorrection, 2000,-1,1);
			HListTOFCorrection.Add(myhistTOFCorrection);
		}
		TString histonameTOFCorrection2 = TString::Format("TOFCorrection_BX_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterTOFCorrection);
		TH1F *myhistTOFCorrection2 = ((TH1F *)(HListTOFCorrection.FindObject(histonameTOFCorrection2)));
		if(!myhistTOFCorrection2){
			myhistTOFCorrection2 = new TH1F(histonameTOFCorrection2, histonameTOFCorrection2, 2000,-1,1);
			HListTOFCorrection.Add(myhistTOFCorrection2);
		}

	}

	TDirectory *xCoordinate = BAM->mkdir("xCoordinate");
	xCoordinate->cd();
	TObjArray HListxCoordinate(0);
	int chamberCounterxCoordinate = 0;
	while(chamberCounterxCoordinate < 36){
		chamberCounterxCoordinate++;
		TString histonamexCoordinate = TString::Format("xCoordinate_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterxCoordinate);
		TH1F *myhistxCoordinate = ((TH1F *)(HListxCoordinate.FindObject(histonamexCoordinate)));
		if(!myhistxCoordinate){
			myhistxCoordinate = new TH1F(histonamexCoordinate, histonamexCoordinate, 50000,-25000,25000);
			HListxCoordinate.Add(myhistxCoordinate);
		}
		TString histonamexCoordinate2 = TString::Format("xCoordinate_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterxCoordinate);
		TH1F *myhistxCoordinate2 = ((TH1F *)(HListxCoordinate.FindObject(histonamexCoordinate2)));
		if(!myhistxCoordinate2){
			myhistxCoordinate2 = new TH1F(histonamexCoordinate2, histonamexCoordinate2, 50000,-25000,25000);
			HListxCoordinate.Add(myhistxCoordinate2);
		}

	}

	TDirectory *yCoordinate = BAM->mkdir("yCoordinate");
	yCoordinate->cd();
	TObjArray HListyCoordinate(0);
	int chamberCounteryCoordinate = 0;
	while(chamberCounteryCoordinate < 36){
		chamberCounteryCoordinate++;
		TString histonameyCoordinate = TString::Format("yCoordinate_pos_ME%d%d_%d",stationInput,ringInput,chamberCounteryCoordinate);
		TH1F *myhistyCoordinate = ((TH1F *)(HListyCoordinate.FindObject(histonameyCoordinate)));
		if(!myhistyCoordinate){
			myhistyCoordinate = new TH1F(histonameyCoordinate, histonameyCoordinate, 50000,-25000,25000);
			HListyCoordinate.Add(myhistyCoordinate);
		}
		TString histonameyCoordinate2 = TString::Format("yCoordinate_neg_ME%d%d_%d",stationInput,ringInput,chamberCounteryCoordinate);
		TH1F *myhistyCoordinate2 = ((TH1F *)(HListyCoordinate.FindObject(histonameyCoordinate2)));
		if(!myhistyCoordinate2){
			myhistyCoordinate2 = new TH1F(histonameyCoordinate2, histonameyCoordinate2, 50000,-25000,25000);
			HListyCoordinate.Add(myhistyCoordinate2);
		}

	}
	TDirectory *zCoordinate = BAM->mkdir("zCoordinate");
	zCoordinate->cd();
	TObjArray HListzCoordinate(0);
	int chamberCounterzCoordinate = 0;
	while(chamberCounterzCoordinate < 36){
		chamberCounterzCoordinate++;
		TString histonamezCoordinate = TString::Format("zCoordinate_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterzCoordinate);
		TH1F *myhistzCoordinate = ((TH1F *)(HListzCoordinate.FindObject(histonamezCoordinate)));
		if(!myhistzCoordinate){
			myhistzCoordinate = new TH1F(histonamezCoordinate, histonamezCoordinate, 50000,-25000,25000);
			HListzCoordinate.Add(myhistzCoordinate);
		}
		TString histonamezCoordinate2 = TString::Format("zCoordinate_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterzCoordinate);
		TH1F *myhistzCoordinate2 = ((TH1F *)(HListzCoordinate.FindObject(histonamezCoordinate2)));
		if(!myhistzCoordinate2){
			myhistzCoordinate2 = new TH1F(histonamezCoordinate2, histonamezCoordinate2, 50000,-25000,25000);
			HListzCoordinate.Add(myhistzCoordinate2);
		}

	}


	TDirectory *PropCorrection = BAM->mkdir("PropCorrection");
	PropCorrection->cd();
	TObjArray HListPropCorrection(0);
	int chamberCounterPropCorrection = 0;
	while(chamberCounterPropCorrection < 36){
		chamberCounterPropCorrection++;
		TString histonamePropCorrection = TString::Format("PropCorrection_BX_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterPropCorrection);
		TH1F *myhistPropCorrection = ((TH1F *)(HListPropCorrection.FindObject(histonamePropCorrection)));
		if(!myhistPropCorrection){
			myhistPropCorrection = new TH1F(histonamePropCorrection, histonamePropCorrection, 2000,-1,1);
			HListPropCorrection.Add(myhistPropCorrection);
		}
		TString histonamePropCorrection2 = TString::Format("PropCorrection_BX_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterPropCorrection);
		TH1F *myhistPropCorrection2 = ((TH1F *)(HListPropCorrection.FindObject(histonamePropCorrection2)));
		if(!myhistPropCorrection2){
			myhistPropCorrection2 = new TH1F(histonamePropCorrection2, histonamePropCorrection2, 2000,-1,1);
			HListPropCorrection.Add(myhistPropCorrection2);
		}

	}



	TDirectory *IncidentAngle = BAM->mkdir("IncidentAngle");
	IncidentAngle->cd();
	TObjArray HListIncidentAngle(0);
	int chamberCounterIncidentAngle = 0;
	while(chamberCounterIncidentAngle < 36){
		chamberCounterIncidentAngle++;
		TString histonameIncidentAngle = TString::Format("IncidentAngle_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterIncidentAngle);
		TH1F *myhistIncidentAngle = ((TH1F *)(HListIncidentAngle.FindObject(histonameIncidentAngle)));
		if(!myhistIncidentAngle){
			myhistIncidentAngle = new TH1F(histonameIncidentAngle, histonameIncidentAngle, 400,-2,2);
			HListIncidentAngle.Add(myhistIncidentAngle);
		}
		TString histonameIncidentAngle2 = TString::Format("IncidentAngle_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterIncidentAngle);
		TH1F *myhistIncidentAngle2 = ((TH1F *)(HListIncidentAngle.FindObject(histonameIncidentAngle2)));
		if(!myhistIncidentAngle2){
			myhistIncidentAngle2 = new TH1F(histonameIncidentAngle2, histonameIncidentAngle2, 400,-2,2);
			HListIncidentAngle.Add(myhistIncidentAngle2);
		}

	}

	TDirectory *IncidentAngleNeighbors = BAM->mkdir("IncidentAngleNeighbors");
	IncidentAngleNeighbors->cd();
	TObjArray HListIncidentAngleNeighbors(0);
	int chamberCounterIncidentAngleNeighbors = 0;
	while(chamberCounterIncidentAngleNeighbors < 36){
		chamberCounterIncidentAngleNeighbors++;
		TString histonameIncidentAngleNeighbors = TString::Format("IncidentAngleNeighbors_pos_ME%d%d_%d",stationInput,ringInput,chamberCounterIncidentAngleNeighbors);
		TH1F *myhistIncidentAngleNeighbors = ((TH1F *)(HListIncidentAngleNeighbors.FindObject(histonameIncidentAngleNeighbors)));
		if(!myhistIncidentAngleNeighbors){
			myhistIncidentAngleNeighbors = new TH1F(histonameIncidentAngleNeighbors, histonameIncidentAngleNeighbors, 400,-2,2);
			HListIncidentAngleNeighbors.Add(myhistIncidentAngleNeighbors);
		}
		TString histonameIncidentAngleNeighbors2 = TString::Format("IncidentAngleNeighbors_neg_ME%d%d_%d",stationInput,ringInput,chamberCounterIncidentAngleNeighbors);
		TH1F *myhistIncidentAngleNeighbors2 = ((TH1F *)(HListIncidentAngleNeighbors.FindObject(histonameIncidentAngleNeighbors2)));
		if(!myhistIncidentAngleNeighbors2){
			myhistIncidentAngleNeighbors2 = new TH1F(histonameIncidentAngleNeighbors2, histonameIncidentAngleNeighbors2, 400,-2,2);
			HListIncidentAngleNeighbors.Add(myhistIncidentAngleNeighbors2);
		}

	}


	TDirectory *BX_Prop_ByChamber_posME1 = BAM->mkdir("BX_Prop_ByChamber_posME1");
	TDirectory *BX_Prop_ByChamber_posME2 = BAM->mkdir("BX_Prop_ByChamber_posME2");
	TDirectory *BX_Prop_ByChamber_posME3 = BAM->mkdir("BX_Prop_ByChamber_posME3");
	TDirectory *BX_Prop_ByChamber_posME4 = BAM->mkdir("BX_Prop_ByChamber_posME4");


	TObjArray HList_prop(0);
	int stationPlotCounterA = 0;
	while(stationPlotCounterA < 4){
		stationPlotCounterA++;
		int ringPlotCounterA = 0;
		if(stationPlotCounterA == 1){
			BX_Prop_ByChamber_posME1->cd();
			while(ringPlotCounterA < 3){
				ringPlotCounterA++;
				int chamberPlotCounterA = 0;
				while(chamberPlotCounterA < 36){
					chamberPlotCounterA++;
					TString histoname = TString::Format("BX_Prop_posME_%d%d_%d_run%d",(stationPlotCounterA),(ringPlotCounterA),(chamberPlotCounterA),(runNumber));
					TH1F *myhist = ((TH1F *)(HList_prop.FindObject(histoname)));
					if(!myhist){
						myhist = new TH1F(histoname, histoname, 700,3,10);
						HList_prop.Add(myhist);
					}

					TString histoname2 = TString::Format("BX_Prop_posME_%d%d_%d_run%d_Conditions",(stationPlotCounterA),(ringPlotCounterA),(chamberPlotCounterA),(runNumber));
					TH1F *myhist2 = ((TH1F *)(HList_prop.FindObject(histoname2)));
					if(!myhist2){
						myhist2 = new TH1F(histoname2, histoname2, 700,3,10);
						HList_prop.Add(myhist2);
					}

				}
			}
		}else{                        if(stationPlotCounterA == 2) BX_Prop_ByChamber_posME2->cd();
			if(stationPlotCounterA == 3) BX_Prop_ByChamber_posME3->cd();
			if(stationPlotCounterA == 4) BX_Prop_ByChamber_posME4->cd();
			while(ringPlotCounterA < 2){                                
				ringPlotCounterA++;
				int chamberPlotCounterA = 0;
				while(chamberPlotCounterA < 36){
					chamberPlotCounterA++;
					TString histoname = TString::Format("BX_Prop_posME_%d%d_%d_run%d",(stationPlotCounterA),(ringPlotCounterA),(chamberPlotCounterA),(runNumber));
					TH1F *myhist = ((TH1F *)(HList_prop.FindObject(histoname)));
					if(!myhist){
						myhist = new TH1F(histoname, histoname,700,3,10);
						HList_prop.Add(myhist);
					}
					TString histoname3 = TString::Format("BX_Prop_posME_%d%d_%d_run%d_Conditions",(stationPlotCounterA),(ringPlotCounterA),(chamberPlotCounterA),(runNumber));
					TH1F *myhist3 = ((TH1F *)(HList_prop.FindObject(histoname3)));
					if(!myhist3){
						myhist3 = new TH1F(histoname3, histoname3,700,3,10);
						HList_prop.Add(myhist3);
					}
				}
			}
		}
	}


	TDirectory *BX_Prop_ByChamber_negME1 = BAM->mkdir("BX_Prop_ByChamber_negME1");
	TDirectory *BX_Prop_ByChamber_negME2 = BAM->mkdir("BX_Prop_ByChamber_negME2");
	TDirectory *BX_Prop_ByChamber_negME3 = BAM->mkdir("BX_Prop_ByChamber_negME3");
	TDirectory *BX_Prop_ByChamber_negME4 = BAM->mkdir("BX_Prop_ByChamber_negME4");     

	TObjArray HList_propneg(0);
	int stationPlotCounter_Propneg = 0;
	while(stationPlotCounter_Propneg < 4){
		stationPlotCounter_Propneg++;
		int ringPlotCounter_Propneg = 0;
		if(stationPlotCounter_Propneg == 1){
			BX_Prop_ByChamber_negME1->cd();
			while(ringPlotCounter_Propneg < 3){
				ringPlotCounter_Propneg++;
				int chamberPlotCounter_Propneg = 0;
				while(chamberPlotCounter_Propneg < 36){
					chamberPlotCounter_Propneg++;
					TString histoname = TString::Format("BX_Prop_negME_%d%d_%d_run%d",(stationPlotCounter_Propneg),(ringPlotCounter_Propneg),(chamberPlotCounter_Propneg),(runNumber));
					TH1F *myhist = ((TH1F *)(HList_propneg.FindObject(histoname)));
					if(!myhist){
						myhist = new TH1F(histoname, histoname, 7,3,10);
						HList_propneg.Add(myhist);
					}
					TString histoname2 = TString::Format("BX_Prop_negME_%d%d_%d_run%d_Conditions",(stationPlotCounter_Propneg),(ringPlotCounter_Propneg),(chamberPlotCounter_Propneg),(runNumber));
					TH1F *myhist2 = ((TH1F *)(HList_propneg.FindObject(histoname2)));
					if(!myhist2){
						myhist2 = new TH1F(histoname2, histoname2, 7,3,10);
						HList_propneg.Add(myhist2);
					}

				}
			}
		}else{
			if(stationPlotCounter_Propneg == 2) BX_Prop_ByChamber_negME2->cd();
			if(stationPlotCounter_Propneg == 3) BX_Prop_ByChamber_negME3->cd();
			if(stationPlotCounter_Propneg == 4) BX_Prop_ByChamber_negME4->cd();
			while(ringPlotCounter_Propneg < 2){
				ringPlotCounter_Propneg++;
				int chamberPlotCounter_Propneg = 0;
				while(chamberPlotCounter_Propneg < 36){
					chamberPlotCounter_Propneg++;
					TString histoname = TString::Format("BX_Prop_negME_%d%d_%d_run%d",(stationPlotCounter_Propneg),(ringPlotCounter_Propneg),(chamberPlotCounter_Propneg),(runNumber));
					TH1F *myhist = ((TH1F *)(HList_propneg.FindObject(histoname)));
					if(!myhist){
						myhist = new TH1F(histoname, histoname,7,3,10);
						HList_propneg.Add(myhist);
					}
					TString histoname3 = TString::Format("BX_Prop_negME_%d%d_%d_run%d_Conditions",(stationPlotCounter_Propneg),(ringPlotCounter_Propneg),(chamberPlotCounter_Propneg),(runNumber));
					TH1F *myhist3 = ((TH1F *)(HList_propneg.FindObject(histoname3)));
					if(!myhist3){
						myhist3 = new TH1F(histoname3, histoname3,7,3,10);
						HList_propneg.Add(myhist3);
					}
				}
			}
		}
	}


	TDirectory *BX_ByChamber_posME1 = BAM->mkdir("BX_ByChamber_posME1");
	TDirectory *BX_ByChamber_posME2 = BAM->mkdir("BX_ByChamber_posME2");
	TDirectory *BX_ByChamber_posME3 = BAM->mkdir("BX_ByChamber_posME3");
	TDirectory *BX_ByChamber_posME4 = BAM->mkdir("BX_ByChamber_posME4");

	TDirectory *BX_ByChamber_negME1 = BAM->mkdir("BX_ByChamber_negME1");
	TDirectory *BX_ByChamber_negME2 = BAM->mkdir("BX_ByChamber_negME2");
	TDirectory *BX_ByChamber_negME3 = BAM->mkdir("BX_ByChamber_negME3");
	TDirectory *BX_ByChamber_negME4 = BAM->mkdir("BX_ByChamber_negME4");

	TObjArray HList(0);
	int stationPlotCounter = 0;
	while(stationPlotCounter < 4){
		stationPlotCounter++;
		int ringPlotCounter = 0;
		if(stationPlotCounter == 1){
			BX_ByChamber_posME1->cd();
			while(ringPlotCounter < 3){
				ringPlotCounter++;
				int chamberPlotCounter = 0;
				while(chamberPlotCounter < 36){
					chamberPlotCounter++;
					TString histoname = TString::Format("BX_posME_%d%d_%d_run%d",(stationPlotCounter),(ringPlotCounter),(chamberPlotCounter),(runNumber));
					TH1F *myhist = ((TH1F *)(HList.FindObject(histoname)));
					if(!myhist){
						myhist = new TH1F(histoname, histoname, 7,3,10);
						HList.Add(myhist);
					}
					TString histoname2 = TString::Format("BX_posME_%d%d_%d_run%d_Conditions",(stationPlotCounter),(ringPlotCounter),(chamberPlotCounter),(runNumber));
					TH1F *myhist2 = ((TH1F *)(HList.FindObject(histoname2)));
					if(!myhist2){
						myhist2 = new TH1F(histoname2, histoname2, 7,3,10);
						HList.Add(myhist2);
					}

				}		
			}
		}else{
			if(stationPlotCounter == 2) BX_ByChamber_posME2->cd();
			if(stationPlotCounter == 3) BX_ByChamber_posME3->cd();
			if(stationPlotCounter == 4) BX_ByChamber_posME4->cd();
			while(ringPlotCounter < 2){
				ringPlotCounter++;
				int chamberPlotCounter = 0;
				while(chamberPlotCounter < 36){
					chamberPlotCounter++;
					TString histoname = TString::Format("BX_posME_%d%d_%d_run%d",(stationPlotCounter),(ringPlotCounter),(chamberPlotCounter),(runNumber));
					TH1F *myhist = ((TH1F *)(HList.FindObject(histoname)));
					if(!myhist){
						myhist = new TH1F(histoname, histoname,7,3,10);
						HList.Add(myhist);
					}	
					TString histoname3 = TString::Format("BX_posME_%d%d_%d_run%d_Conditions",(stationPlotCounter),(ringPlotCounter),(chamberPlotCounter),(runNumber));
					TH1F *myhist3 = ((TH1F *)(HList.FindObject(histoname3)));
					if(!myhist3){
						myhist3 = new TH1F(histoname3, histoname3,7,3,10);
						HList.Add(myhist3);
					}
				}		
			}
		}
	}

	TObjArray HListneg(0);
	int stationPlotCounterneg = 0;
	while(stationPlotCounterneg < 4){
		stationPlotCounterneg++;
		int ringPlotCounterneg = 0;
		if(stationPlotCounterneg == 1){
			BX_ByChamber_negME1->cd();
			while(ringPlotCounterneg < 3){
				ringPlotCounterneg++;
				int chamberPlotCounterneg = 0;
				while(chamberPlotCounterneg < 36){
					chamberPlotCounterneg++;
					TString histoname = TString::Format("BX_negME_%d%d_%d_run%d",(stationPlotCounterneg),(ringPlotCounterneg),(chamberPlotCounterneg),(runNumber));
					TH1F *myhist = ((TH1F *)(HListneg.FindObject(histoname)));
					if(!myhist){
						myhist = new TH1F(histoname, histoname, 7,3,10);
						HListneg.Add(myhist);
					}
					TString histoname2 = TString::Format("BX_negME_%d%d_%d_run%d_Conditions",(stationPlotCounterneg),(ringPlotCounterneg),(chamberPlotCounterneg),(runNumber));
					TH1F *myhist2 = ((TH1F *)(HListneg.FindObject(histoname2)));
					if(!myhist2){
						myhist2 = new TH1F(histoname2, histoname2, 7,3,10);
						HListneg.Add(myhist2);
					}

				}		
			}
		}else{
			if(stationPlotCounterneg == 2) BX_ByChamber_negME2->cd();
			if(stationPlotCounterneg == 3) BX_ByChamber_negME3->cd();
			if(stationPlotCounterneg == 4) BX_ByChamber_negME4->cd();
			while(ringPlotCounterneg < 2){
				ringPlotCounterneg++;
				int chamberPlotCounterneg = 0;
				while(chamberPlotCounterneg < 36){
					chamberPlotCounterneg++;
					TString histoname = TString::Format("BX_negME_%d%d_%d_run%d",(stationPlotCounterneg),(ringPlotCounterneg),(chamberPlotCounterneg),(runNumber));
					TH1F *myhist = ((TH1F *)(HListneg.FindObject(histoname)));
					if(!myhist){
						myhist = new TH1F(histoname, histoname,7,3,10);
						HListneg.Add(myhist);
					}	
					TString histoname3 = TString::Format("BX_negME_%d%d_%d_run%d_Conditions",(stationPlotCounterneg),(ringPlotCounterneg),(chamberPlotCounterneg),(runNumber));
					TH1F *myhist3 = ((TH1F *)(HListneg.FindObject(histoname3)));
					if(!myhist3){
						myhist3 = new TH1F(histoname3, histoname3,7,3,10);
						HListneg.Add(myhist3);
					}
				}		
			}
		}
	}

	//Make the walk-around plots for each ring

	TDirectory *WalkAround_Rings_ME11 = BAM->mkdir("WalkAround_Rings_ME11");
	TDirectory *WalkAround_Rings_ME12 = BAM->mkdir("WalkAround_Rings_ME12");
	TDirectory *WalkAround_Rings_ME13 = BAM->mkdir("WalkAround_Rings_ME13");
	TDirectory *WalkAround_Rings_ME21 = BAM->mkdir("WalkAround_Rings_ME21");
	TDirectory *WalkAround_Rings_ME22 = BAM->mkdir("WalkAround_Rings_ME22");
	TDirectory *WalkAround_Rings_ME31 = BAM->mkdir("WalkAround_Rings_ME31");
	TDirectory *WalkAround_Rings_ME32 = BAM->mkdir("WalkAround_Rings_ME32");
	TDirectory *WalkAround_Rings_ME41 = BAM->mkdir("WalkAround_Rings_ME41");
	TDirectory *WalkAround_Rings_ME42 = BAM->mkdir("WalkAround_Rings_ME42");

	TObjArray WalkAround(0);
	for(int endcap = 0; endcap < 2; endcap++){
		int endcapActual = 0;
		if(endcap == 0) endcapActual = 1;
		if(endcap == 1) endcapActual = -1;
		if(endcapActual == 0) continue;
		int stationPlotCounterWA = 0;
		while(stationPlotCounterWA < 4){
			stationPlotCounterWA++;
			int ringPlotCounterWA = 0;
			while(ringPlotCounterWA < 3){
				ringPlotCounterWA++;
				int chamberPlotCounterWA = 0;
				int chamberPlotCounterAdditionWA = 0;
				if(stationPlotCounterWA == 1 || ringPlotCounterWA == 2){
					while(chamberPlotCounterWA < 36){
						chamberPlotCounterWA++;

						if(chamberPlotCounterWA != 36)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 36)  chamberPlotCounterAdditionWA = 1;

						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 1) WalkAround_Rings_ME11->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 2) WalkAround_Rings_ME12->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 3) WalkAround_Rings_ME13->cd();
						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 2) WalkAround_Rings_ME22->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 2) WalkAround_Rings_ME32->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 2) WalkAround_Rings_ME42->cd();

						TString histoname = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist = ((TH1F *)(WalkAround.FindObject(histoname)));
						if(!myhist){
							myhist = new TH1F(histoname, histoname, 20,-10,10);
							WalkAround.Add(myhist);
						}
						TString histoname_ChamberA = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH1F *myhist_ChamberA = ((TH1F *)(WalkAround.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH1F(histoname_ChamberA, histoname_ChamberA, 20,-10,10);
							WalkAround.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist_ChamberB = ((TH1F *)(WalkAround.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH1F(histoname_ChamberB, histoname_ChamberB, 20,-10,10);
							WalkAround.Add(myhist_ChamberB);
						}
					}
				}else{
					while(chamberPlotCounterWA < 18){
						chamberPlotCounterWA++;
						if(chamberPlotCounterWA != 18)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 18)  chamberPlotCounterAdditionWA = 1;

						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 1) WalkAround_Rings_ME21->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 1) WalkAround_Rings_ME31->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 1) WalkAround_Rings_ME41->cd();

						TString histoname = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA,endcapActual);
						TH1F *myhist = ((TH1F *)(WalkAround.FindObject(histoname)));
						TString histoname_ChamberA = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH1F *myhist_ChamberA = ((TH1F *)(WalkAround.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH1F(histoname_ChamberA, histoname_ChamberA, 20,-10,10);
							WalkAround.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist_ChamberB = ((TH1F *)(WalkAround.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH1F(histoname_ChamberB, histoname_ChamberB, 20,-10,10);
							WalkAround.Add(myhist_ChamberB);
						}
						if(!myhist){
							myhist = new TH1F(histoname, histoname, 20,-10,10);
							WalkAround.Add(myhist);
						}
					}
				}
			}
		}
	}

	//Make the walk-around plots for each ring including signal propogation time

	TDirectory *WalkAroundProp_Rings_ME11 = BAM->mkdir("WalkAroundProp_Rings_ME11");
	TDirectory *WalkAroundProp_Rings_ME12 = BAM->mkdir("WalkAroundProp_Rings_ME12");
	TDirectory *WalkAroundProp_Rings_ME13 = BAM->mkdir("WalkAroundProp_Rings_ME13");
	TDirectory *WalkAroundProp_Rings_ME21 = BAM->mkdir("WalkAroundProp_Rings_ME21");
	TDirectory *WalkAroundProp_Rings_ME22 = BAM->mkdir("WalkAroundProp_Rings_ME22");
	TDirectory *WalkAroundProp_Rings_ME31 = BAM->mkdir("WalkAroundProp_Rings_ME31");
	TDirectory *WalkAroundProp_Rings_ME32 = BAM->mkdir("WalkAroundProp_Rings_ME32");
	TDirectory *WalkAroundProp_Rings_ME41 = BAM->mkdir("WalkAroundProp_Rings_ME41");
	TDirectory *WalkAroundProp_Rings_ME42 = BAM->mkdir("WalkAroundProp_Rings_ME42");

	TObjArray WalkAroundProp(0);
	for(int endcap = 0; endcap < 2; endcap++){
		int endcapActual = 0;
		if(endcap == 0) endcapActual = 1;
		if(endcap == 1) endcapActual = -1;
		if(endcapActual == 0) continue;
		int stationPlotCounterWA = 0;
		while(stationPlotCounterWA < 4){
			stationPlotCounterWA++;
			int ringPlotCounterWA = 0;
			while(ringPlotCounterWA < 3){
				ringPlotCounterWA++;
				int chamberPlotCounterWA = 0;
				int chamberPlotCounterAdditionWA = 0;
				if(stationPlotCounterWA == 1 || ringPlotCounterWA == 2){
					while(chamberPlotCounterWA < 36){
						chamberPlotCounterWA++;

						if(chamberPlotCounterWA != 36)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 36)  chamberPlotCounterAdditionWA = 1;

						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 1) WalkAroundProp_Rings_ME11->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 2) WalkAroundProp_Rings_ME12->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 3) WalkAroundProp_Rings_ME13->cd();
						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 2) WalkAroundProp_Rings_ME22->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 2) WalkAroundProp_Rings_ME32->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 2) WalkAroundProp_Rings_ME42->cd();

						TString histoname = TString::Format("WalkAroundProp_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist = ((TH1F *)(WalkAroundProp.FindObject(histoname)));
						if(!myhist){
							myhist = new TH1F(histoname, histoname, 1400,-10,10);
							WalkAroundProp.Add(myhist);
						}
						TString histoname_ChamberA = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH1F *myhist_ChamberA = ((TH1F *)(WalkAroundProp.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH1F(histoname_ChamberA, histoname_ChamberA, 1400,-10,10);
							WalkAroundProp.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist_ChamberB = ((TH1F *)(WalkAroundProp.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH1F(histoname_ChamberB, histoname_ChamberB, 1400,-10,10);
							WalkAroundProp.Add(myhist_ChamberB);

						}
					}
				}else{
					while(chamberPlotCounterWA < 18){
						chamberPlotCounterWA++;
						if(chamberPlotCounterWA != 18)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 18)  chamberPlotCounterAdditionWA = 1;

						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 1) WalkAroundProp_Rings_ME21->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 1) WalkAroundProp_Rings_ME31->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 1) WalkAroundProp_Rings_ME41->cd();

						TString histoname = TString::Format("WalkAroundProp_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA,endcapActual);
						TH1F *myhist = ((TH1F *)(WalkAroundProp.FindObject(histoname)));
						TString histoname_ChamberA = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH1F *myhist_ChamberA = ((TH1F *)(WalkAroundProp.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH1F(histoname_ChamberA, histoname_ChamberA, 1400,-10,10);
							WalkAroundProp.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist_ChamberB = ((TH1F *)(WalkAroundProp.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH1F(histoname_ChamberB, histoname_ChamberB, 1400,-10,10);
							WalkAroundProp.Add(myhist_ChamberB);
						}
						if(!myhist){
							myhist = new TH1F(histoname, histoname, 1400,-10,10);
							WalkAroundProp.Add(myhist);
						}
					}
				}
			}
		}
	}




	//Make the walk-around plots for each ring using TOF

	TDirectory *WalkAroundwithTOF_Rings_ME11 = BAM->mkdir("WalkAroundwithTOF_Rings_ME11");
	TDirectory *WalkAroundwithTOF_Rings_ME12 = BAM->mkdir("WalkAroundwithTOF_Rings_ME12");
	TDirectory *WalkAroundwithTOF_Rings_ME13 = BAM->mkdir("WalkAroundwithTOF_Rings_ME13");
	TDirectory *WalkAroundwithTOF_Rings_ME21 = BAM->mkdir("WalkAroundwithTOF_Rings_ME21");
	TDirectory *WalkAroundwithTOF_Rings_ME22 = BAM->mkdir("WalkAroundwithTOF_Rings_ME22");
	TDirectory *WalkAroundwithTOF_Rings_ME31 = BAM->mkdir("WalkAroundwithTOF_Rings_ME31");
	TDirectory *WalkAroundwithTOF_Rings_ME32 = BAM->mkdir("WalkAroundwithTOF_Rings_ME32");
	TDirectory *WalkAroundwithTOF_Rings_ME41 = BAM->mkdir("WalkAroundwithTOF_Rings_ME41");
	TDirectory *WalkAroundwithTOF_Rings_ME42 = BAM->mkdir("WalkAroundwithTOF_Rings_ME42");

	TObjArray WalkAroundwithTOF(0);
	for(int endcap = 0; endcap < 2; endcap++){
		int endcapActual = 0;
		if(endcap == 0) endcapActual = 1;
		if(endcap == 1) endcapActual = -1;
		if(endcapActual == 0) continue;
		int stationPlotCounterWA = 0;
		while(stationPlotCounterWA < 4){
			stationPlotCounterWA++;
			int ringPlotCounterWA = 0;
			while(ringPlotCounterWA < 3){
				ringPlotCounterWA++;
				int chamberPlotCounterWA = 0;
				int chamberPlotCounterAdditionWA = 0;
				if(stationPlotCounterWA == 1 || ringPlotCounterWA == 2){
					while(chamberPlotCounterWA < 36){
						chamberPlotCounterWA++;

						if(chamberPlotCounterWA != 36)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 36)  chamberPlotCounterAdditionWA = 1;
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 1) WalkAroundwithTOF_Rings_ME11->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 2) WalkAroundwithTOF_Rings_ME12->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 3) WalkAroundwithTOF_Rings_ME13->cd();
						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 2) WalkAroundwithTOF_Rings_ME22->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 2) WalkAroundwithTOF_Rings_ME32->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 2) WalkAroundwithTOF_Rings_ME42->cd();

						TString histoname = TString::Format("WalkAroundwithTOF_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist = ((TH1F *)(WalkAroundwithTOF.FindObject(histoname)));
						if(!myhist){
							myhist = new TH1F(histoname, histoname, 8000,-10,10);
							WalkAroundwithTOF.Add(myhist);
						}
						TString histoname_ChamberA = TString::Format("WalkAroundwithTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH1F *myhist_ChamberA = ((TH1F *)(WalkAroundwithTOF.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH1F(histoname_ChamberA, histoname_ChamberA, 2000,5,10);
							WalkAroundwithTOF.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("WalkAroundwithTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist_ChamberB = ((TH1F *)(WalkAroundwithTOF.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH1F(histoname_ChamberB, histoname_ChamberB, 2000,5,10);
							WalkAroundwithTOF.Add(myhist_ChamberB);
						}
					}
				}else{
					while(chamberPlotCounterWA < 18){
						chamberPlotCounterWA++;
						if(chamberPlotCounterWA != 18)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 18)  chamberPlotCounterAdditionWA = 1;

						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 1) WalkAroundwithTOF_Rings_ME21->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 1) WalkAroundwithTOF_Rings_ME31->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 1) WalkAroundwithTOF_Rings_ME41->cd();

						TString histoname = TString::Format("WalkAroundwithTOF_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA,endcapActual);
						TH1F *myhist = ((TH1F *)(WalkAroundwithTOF.FindObject(histoname)));
						TString histoname_ChamberA = TString::Format("WalkAroundwithTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH1F *myhist_ChamberA = ((TH1F *)(WalkAroundwithTOF.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH1F(histoname_ChamberA, histoname_ChamberA, 2000,5,10);
							WalkAroundwithTOF.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("WalkAroundwithTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist_ChamberB = ((TH1F *)(WalkAroundwithTOF.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH1F(histoname_ChamberB, histoname_ChamberB, 2000,5,10);
							WalkAroundwithTOF.Add(myhist_ChamberB);
						}
						if(!myhist){
							myhist = new TH1F(histoname, histoname, 2000,5,10);
							WalkAroundwithTOF.Add(myhist);
						}
					}
				}
			}
		}
	}


	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME11 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME11");
	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME12 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME12");
	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME13 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME13");
	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME21 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME21");
	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME22 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME22");
	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME31 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME31");
	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME32 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME32");
	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME41 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME41");
	TDirectory *WalkAroundwithTOF_And_Prop_Rings_ME42 = BAM->mkdir("WalkAroundwithTOF_And_Prop_Rings_ME42");

	TObjArray WalkAroundwithTOF_And_Prop(0);
	for(int endcap = 0; endcap < 2; endcap++){
		int endcapActual = 0;
		if(endcap == 0) endcapActual = 1;
		if(endcap == 1) endcapActual = -1;
		if(endcapActual == 0) continue;
		int stationPlotCounterWA = 0;
		while(stationPlotCounterWA < 4){
			stationPlotCounterWA++;
			int ringPlotCounterWA = 0;
			while(ringPlotCounterWA < 3){
				ringPlotCounterWA++;
				int chamberPlotCounterWA = 0;
				int chamberPlotCounterAdditionWA = 0;
				if(stationPlotCounterWA == 1 || ringPlotCounterWA == 2){
					while(chamberPlotCounterWA < 36){
						chamberPlotCounterWA++;

						if(chamberPlotCounterWA != 36)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 36)  chamberPlotCounterAdditionWA = 1;
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 1) WalkAroundwithTOF_And_Prop_Rings_ME11->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 2) WalkAroundwithTOF_And_Prop_Rings_ME12->cd();
						if(stationPlotCounterWA == 1 && ringPlotCounterWA == 3) WalkAroundwithTOF_And_Prop_Rings_ME13->cd();
						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 2) WalkAroundwithTOF_And_Prop_Rings_ME22->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 2) WalkAroundwithTOF_And_Prop_Rings_ME32->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 2) WalkAroundwithTOF_And_Prop_Rings_ME42->cd();

						TString histoname = TString::Format("WalkAroundwithTOF_And_Prop_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(histoname)));
						if(!myhist){
							myhist = new TH1F(histoname, histoname, 8000,-10,10);
							WalkAroundwithTOF_And_Prop.Add(myhist);
						}
						TString histoname_ChamberA = TString::Format("WalkAroundwithTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH1F *myhist_ChamberA = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH1F(histoname_ChamberA, histoname_ChamberA, 2000,5,10);
							WalkAroundwithTOF_And_Prop.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("WalkAroundwithTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist_ChamberB = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH1F(histoname_ChamberB, histoname_ChamberB, 2000,5,10);
							WalkAroundwithTOF_And_Prop.Add(myhist_ChamberB);
						}
					}
				}else{
					while(chamberPlotCounterWA < 18){
						chamberPlotCounterWA++;
						if(chamberPlotCounterWA != 18)  chamberPlotCounterAdditionWA = chamberPlotCounterWA + 1;
						if(chamberPlotCounterWA == 18)  chamberPlotCounterAdditionWA = 1;

						if(stationPlotCounterWA == 2 && ringPlotCounterWA == 1) WalkAroundwithTOF_And_Prop_Rings_ME21->cd();
						if(stationPlotCounterWA == 3 && ringPlotCounterWA == 1) WalkAroundwithTOF_And_Prop_Rings_ME31->cd();
						if(stationPlotCounterWA == 4 && ringPlotCounterWA == 1) WalkAroundwithTOF_And_Prop_Rings_ME41->cd();

						TString histoname = TString::Format("WalkAroundwithTOF_And_Prop_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA,endcapActual);
						TH1F *myhist = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(histoname)));
						TString histoname_ChamberA = TString::Format("WalkAroundwithTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, endcapActual);
						TH1F *myhist_ChamberA = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(histoname_ChamberA)));
						if(!myhist_ChamberA){
							myhist_ChamberA = new TH1F(histoname_ChamberA, histoname_ChamberA, 2000,5,10);
							WalkAroundwithTOF_And_Prop.Add(myhist_ChamberA);
						}
						TString histoname_ChamberB = TString::Format("WalkAroundwithTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, stationPlotCounterWA, ringPlotCounterWA, chamberPlotCounterAdditionWA, endcapActual);
						TH1F *myhist_ChamberB = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(histoname_ChamberB)));
						if(!myhist_ChamberB){
							myhist_ChamberB = new TH1F(histoname_ChamberB, histoname_ChamberB, 2000,5,10);
							WalkAroundwithTOF_And_Prop.Add(myhist_ChamberB);
						}
						if(!myhist){
							myhist = new TH1F(histoname, histoname, 2000,5,10);
							WalkAroundwithTOF_And_Prop.Add(myhist);
						}
					}
				}
			}
		}
	}




	TDirectory *IntermediateProcesses = BAM->mkdir("IntermediateProcesses");
	IntermediateProcesses->cd();


	/*}}}*/

	//cout << "Checkpoint 1" << endl;
	// ------------------------------------------------------------------
	// Loop over the events
	// ------------------------------------------------------------------
	//Global Counters
	/*{{{*/
	int eventCounter = 0;
	int me11trkCounterTotal = 0;
	int dupeCounterTotal = 0;
	int tofandneighborCounter = 0;
	int tofCounter = 0;
	int chiCounter = 0;
	int chiAndtofCounter = 0;
	int tofForNeighborsPlotCounter = 0;
	int chiPlotCounter = 0;
	int chiAndtofPlotCounter = 0;
	int backwardsCounter = 0;	
	int totalCounter = 0;

	double maxPropTime = 0;
	double minPropTime = 999;

	/*}}}*/

	if((numEvents == -1 && eventCounter < GetEntries()) || (numEvents != -1 && eventCounter < numEvents)){
		if(numEvents == -1) numEvents = GetEntries();						
		cout << "This is printing for Run " << runNumber << endl;		
		cout << numEvents << " to process ..." << endl;
		cout << "" << endl;		
		for(Long64_t e = 0; e < numEvents; e++){
			eventCounter++;
			int dupeCounter = 0;	
			//Load the e-th event
			Long64_t ientry = LoadTree(e);
			if(ientry < 0) break;
			//cout << "Checkpoint 7" << endl;
			//cout << "Event: " << e << endl;
			//try{	
			GetEntry(e); //The spash run is breaking here
			//}catch(const length_error& le){
			//	cout << "Length error caught: " << le.what() << endl;		
			//	continue;
			//}
			//cout << "Checkpoint 8" << endl;
			if(e % 10000 == 0){			
				cout << "Event Number: " << e << endl; 
				//cout << "Event Number: " << e << "\r" << std::flush;			
			}
			//if(e <= 557100) continue;
	//		cout << "Event Number: " << e << endl;

			//Event Counters			
			int waTOFcounter=0;
			int NeighborTOFcounter = 0;
			int me11posCounter = 0;
			int me11negCounter = 0;
			int neighborCounter = 0;
			//Find duplicates.  duplicates is the index of the duplicate LCT.
			/*{{{*/
			//There are no duplicate LCTs in a given track, but quite often the duplicates are part
			//of different tracks.  This is fine for most of the code.
			//However, any loop which just looks at raw variables (not just trk variables)
			//needs to be fixed.

			//print everything about these events
			for(int i = 0; i < csctf_->lctBx.size(); i++){
				//cout << "LCT Index " << i << " (Station, Ring, Chamber, GlobalEta, GlobalPhi, BX, Endcap): " << endl;
				//cout <<  csctf_->lctStation[i] << ", " << csctf_->lctRing[i] << ", " << csctf_->lctChamber[i] << ", " << csctf_->lctglobalEta[i] << ", " << csctf_->lctglobalPhi[i] << ", " << csctf_->lctBx[i] << ", " << csctf_->lctEndcap[i] << endl;
				//cout << "" << endl;
			}


			//cout << "Checkpoint 1" << endl;
			//Begin with finding the duplicate LCTs
			vector< vector<int> > duplicateList;
			vector<int> duplicates;
			for(int i = 0; i < csctf_->lctBx.size(); i++){
				vector<int> duplicateVector;
				//////cout << "LCT " << i << " information: " << endl;
				//////cout <<  csctf_->lctStation[i] << ", " << csctf_->lctRing[i] << ", " << csctf_->lctChamber[i] << ", " << csctf_->lctglobalEta[i] << ", " << csctf_->lctglobalPhi[i] << ", " << csctf_->lctBx[i] << endl;
				for(int j = 0; j < csctf_->lctBx.size(); j ++){
					//Say that two LCTs are the same if they have the same station, ring, chamber, BX, eta, and phi.
					if(i != j){
						if(csctf_->lctStation[i]   != csctf_->lctStation[j]) continue;
						if(csctf_->lctRing[i]      != csctf_->lctRing[j]) continue;
						if(csctf_->lctChamber[i]   != csctf_->lctChamber[j]) continue;
						if(csctf_->lctglobalEta[i] != csctf_->lctglobalEta[j]) continue;
						if(csctf_->lctglobalPhi[i] != csctf_->lctglobalPhi[j]) continue;
						if(csctf_->lctBx[i]        != csctf_->lctBx[j]) continue;
						if(csctf_->lctEndcap[i]        != csctf_->lctEndcap[j]) continue;
						//////cout << "" << endl;
						//////cout << "There was a duplicate between " << i << " and " << j << endl;						
						dupeCounter++;
						dupeCounterTotal++;
						duplicateVector.push_back(i);
						duplicateVector.push_back(j); //shouldn't only one be added to the duplicate List? otherwise you throw out the original and the duplicate...?
						duplicates.push_back(i);	
						//LCTCheck
						////cout << "Index " << i << "has been added to the duplicates list and will not be used." << endl;
					}
				}
				duplicateList.push_back(duplicateVector);
			}

			//cout << "Checkpoint 2" << endl;
			//cout << "" << endl;
			//cout << "The Unique LCTs in the event" << endl;
			for(int i = 0; i < csctf_->lctBx.size(); i++){
				if(contains(duplicates,i)) continue;
				//cout << "LCT Index " << i << " information: " << endl;
				//cout <<  csctf_->lctStation[i] << ", " << csctf_->lctRing[i] << ", " << csctf_->lctChamber[i] << ", " << csctf_->lctglobalEta[i] << ", " << csctf_->lctglobalPhi[i] << ", " << csctf_->lctBx[i] << ", " << csctf_->lctEndcap[i] << endl;
				if(csctf_->lctStation[i] == 1 && csctf_->lctRing[i] == 1 && csctf_->lctEndcap[i] == 1) me11posCounter++;
				if(csctf_->lctStation[i] == 1 && csctf_->lctRing[i] == 1 && csctf_->lctEndcap[i] == -1) me11negCounter++;			
			}

			//cout << "Checkpoint 3" << endl;
			NumLCTs->Fill(csctf_->lctBx.size());

			if(me11posCounter > 2 || me11negCounter > 2) ////cout << "There are too many LCTs in ME11 for TOF to be accurate...skipping this event" << endl;
			if(me11posCounter > 2 || me11negCounter > 2) continue;
			/*}}}*/
			//LCTCheck and catagorize LCTs (ie: make station index vectors and histograms of basic information))
			/*{{{*/			
			//Let's look at the events and print some info about events that pass my various criteria.
			//Chi/WA 1 unique LCT in ME1/1 per chamber
			//TOF: 3 unique LCTs. 1 in each station ME1,2,3
			//Chi+TOF: 4 total unique LCTs. 1 in ME2,3. 2 in ME1 (neighbors).
			//For now, I only care about ME11 as shown in if statement below
			bool statOne = false;
			bool statTwo = false;
			bool statThree = false;
			bool statFour = false;
			vector<int> statOneIndex;
			vector<int> statOneChambers;		
			vector<int> statTwoIndex;
			vector<int> statThreeIndex;
			vector<int> statFourIndex;

			vector< vector<int> > neighborList;
			for(int j = 0; j < csctf_->lctBx.size(); j ++){
				if(contains(duplicates, j) == true) continue; //Make sure not to count the duplicates
				//Adding Quality Cut Here:  Renjie says to use Quality > 10
				//if(csctf_->lctQuality[j] <= 10) continue;	

				if(csctf_->lctStation[j] == stationInput && csctf_->lctRing[j] == ringInput){
					statOne = true;
					statOneIndex.push_back(j);
					statOneChambers.push_back(csctf_->lctChamber[j]);			
				}
				if(csctf_->lctStation[j] == 2){
					statTwo = true;
					statTwoIndex.push_back(j);
				}
				if(csctf_->lctStation[j] == 3){
					statThree = true;
					statThreeIndex.push_back(j);
				}
				if(csctf_->lctStation[j] == 4){
					statFour = true;
					statFourIndex.push_back(j);
				}
			}
			//cout << "Checkpoint 4" << endl;
			//A quick look at the Cartesion points
			for(int l = 0; l < csctf_->lctBx.size(); l++){
				if(contains(duplicates, l) == true) continue;
				if(csctf_->lctStation[l] != stationInput) continue;
				if(csctf_->lctRing[l] != ringInput) continue;
				double xpos = getx(csctf_->lctglobalEta[l],csctf_->lctglobalPhi[l],stationInput,ringInput,csctf_->lctSector[l]);
				double ypos = gety(csctf_->lctglobalEta[l],csctf_->lctglobalPhi[l],stationInput,ringInput,csctf_->lctSector[l]);
				double zpos = getz(csctf_->lctglobalEta[l],csctf_->lctglobalPhi[l],stationInput,ringInput,csctf_->lctEndcap[l], csctf_->lctChamber[l]);
				if(csctf_->lctEndcap[l] == 1){
					TString getcoordx = TString::Format("xCoordinate_pos_ME%d%d_%d",stationInput,ringInput,csctf_->lctChamber[l]);
					TH1F *thexhist = ((TH1F *)(HListxCoordinate.FindObject(getcoordx)));
					thexhist->Fill(xpos);
					TString getcoordy = TString::Format("yCoordinate_pos_ME%d%d_%d",stationInput,ringInput,csctf_->lctChamber[l]);
					TH1F *theyhist = ((TH1F *)(HListyCoordinate.FindObject(getcoordy)));
					theyhist->Fill(ypos);
					TString getcoordz = TString::Format("zCoordinate_pos_ME%d%d_%d",stationInput,ringInput,csctf_->lctChamber[l]);
					TH1F *thezhist = ((TH1F *)(HListzCoordinate.FindObject(getcoordz)));
					thezhist->Fill(zpos);
				}else{
					TString getcoordx = TString::Format("xCoordinate_neg_ME%d%d_%d",stationInput,ringInput,csctf_->lctChamber[l]);
					TH1F *thexhist = ((TH1F *)(HListxCoordinate.FindObject(getcoordx)));
					thexhist->Fill(xpos);
					TString getcoordy = TString::Format("yCoordinate_neg_ME%d%d_%d",stationInput,ringInput,csctf_->lctChamber[l]);
					TH1F *theyhist = ((TH1F *)(HListyCoordinate.FindObject(getcoordy)));
					theyhist->Fill(ypos);
					TString getcoordz = TString::Format("zCoordinate_neg_ME%d%d_%d",stationInput,ringInput,csctf_->lctChamber[l]);
					TH1F *thezhist = ((TH1F *)(HListzCoordinate.FindObject(getcoordz)));
					thezhist->Fill(zpos);
				}			
			}
			//cout << "Checkpoint 5" << endl;
			//At this point, I have vectors of indices for each chamber.  All LCTs are unique
			//Need to see if Chamber 1 has at least 2 LCTs and if they are neighbors
			//cout << "The LCTs have now been catagorized." << endl;
			//cout << "Station 1 LCTs:" << endl;			
			//cout << "Index: " << endl;
			print_contents(statOneIndex);
			//cout << "Chambers: " << endl;
			print_contents(statOneChambers);
			////cout << "" << endl;
			////cout << "Station 2 LCTs: " << endl;
			////cout << "Index: " << endl;
			//print_contents(statTwoIndex);
			////cout << "" << endl;
			////cout << "Station 3 LCTs: " << endl;
			////cout << "Index: " << endl;
			//print_contents(statThreeIndex);
			vector<int> statOneIndexPos;			
			vector<int> statOneIndexNeg;			
			vector<int> statTwoIndexPos;			
			vector<int> statTwoIndexNeg;			
			vector<int> statThreeIndexPos;			
			vector<int> statThreeIndexNeg;			
			//A quick catagorization
			for(int i =0; i < statOneIndex.size(); i++){
				if(csctf_->lctEndcap[statOneIndex[i]] == 1) statOneIndexPos.push_back(statOneIndex[i]); 
				if(csctf_->lctEndcap[statOneIndex[i]] == -1) statOneIndexNeg.push_back(statOneIndex[i]); 
			}
			for(int i =0; i < statTwoIndex.size(); i++){
				if(csctf_->lctEndcap[statTwoIndex[i]] == 1) statTwoIndexPos.push_back(statTwoIndex[i]);
				if(csctf_->lctEndcap[statTwoIndex[i]] == -1) statTwoIndexNeg.push_back(statTwoIndex[i]);
			}
			for(int i =0; i < statThreeIndex.size(); i++){
				if(csctf_->lctEndcap[statThreeIndex[i]] == 1) statThreeIndexPos.push_back(statThreeIndex[i]);
				if(csctf_->lctEndcap[statThreeIndex[i]] == -1) statThreeIndexNeg.push_back(statThreeIndex[i]);
			}



			for(int i = 0; i < statOneIndex.size(); i++){
				vector<int> neighbors;
				//cout << "Searching for neighbors..." << endl;
				if(contains(statOneChambers,statOneChambers[i]+1)){
					if(csctf_->lctEndcap[statOneIndex[findIndex(statOneChambers,statOneChambers[i])]] != csctf_->lctEndcap[statOneIndex[findIndex(statOneChambers,statOneChambers[i]+1)]] )continue;	
					neighbors.push_back(statOneIndex[findIndex(statOneChambers,statOneChambers[i])]);
					neighbors.push_back(statOneIndex[findIndex(statOneChambers,statOneChambers[i]+1)]);
					neighborList.push_back(neighbors);
					////cout << "Added Neighbors: " << endl;
					print_contents(neighbors);
					neighborCounter++;				
				}
				if(statOneChambers[i] == numChambersInRing){
					if(contains(statOneChambers, 1)){
						if(csctf_->lctEndcap[statOneIndex[findIndex(statOneChambers,numChambersInRing)]] != csctf_->lctEndcap[statOneIndex[findIndex(statOneChambers,1)]] )continue;	
						neighbors.push_back(statOneIndex[findIndex(statOneChambers,numChambersInRing)]);
						neighbors.push_back(statOneIndex[findIndex(statOneChambers,1)]);
						neighborList.push_back(neighbors);
						////cout << "Added Neighbors: " << endl;
						print_contents(neighbors);					
						neighborCounter++;					
					}
				}
			}
			//cout << "Checkpoint 6" << endl;

			for(int n =0; n < neighborList.size(); n++){
				//cout << "Neighbors List: " << endl;
				for(int m = 0; m < neighborList[n].size(); m++){
					//cout << "neighborList Entry: " << neighborList[n][m] << endl; 				
					//print_contents(neighborList[n]); 
				}
				//cout << "Neighbor Chambers Check: " << csctf_->lctChamber[neighborList[n][0]] << ", " << csctf_->lctChamber[neighborList[n][1]] << endl;
			}
			//cout << "Checkpoint 6.1" << endl;
			//Now I know whether or not the LCTs  are neighbors
			//Checking my LCT number requirements (this will tell me the number of inputs I expect to see in the histograms)		
			if(statOne && statTwo && statThree){
				tofCounter++;
				////cout << "This event has an LCT in Station 1,2,3.  TOF may be calculated." << endl;			
			}
			if(neighborList.size() > 0){
				chiCounter++;
				////cout << "This event has neighboring LCTs in ME11.  A WalkAround may be performed. " << endl;
			}
			if(statOne && statTwo && statThree){
				if(neighborList.size() >0){
					chiAndtofCounter++;
					////cout << "This event has an LCT in Station 1,2,3 and neighbors in ME11.  TOF in conjunction with the WA may be performed." << endl;
				}
			}
			//cout << "Checkpoint 6.2" << endl;
			///
			//Checkpoint
			for(int i=0; i<csctf_->lctBx.size(); i++){
			//  cout << "Filling the ME11_pos_BX plot" << endl;
				if(contains(duplicates, i) == true) continue; //Make sure not to count the duplicates
//cout << "Duplicate cut passed" << endl;			
	//if(csctf_->lctChamber[i] > numChambersInRing) ////cout << "Chamber Number out of range" << endl;
				//positive endcap only

				//if(csctf_->lctQuality[i] <= 10) continue;
//cout << "Quality cut passed" << endl;			
//cout << "endcap, station, ring: " << csctf_->lctEndcap[i] <<  ", " << csctf_->lctStation[i] << ", " << csctf_->lctRing[i] << endl;
				if(csctf_->lctEndcap[i] == 1){
//cout << "Positive endcap cut passed" << endl;				
	//cout << "Checkpoint 0" << endl;
					TString getWirePosHist = TString::Format("MEplus%d%d_WireGroup",csctf_->lctStation[i],csctf_->lctRing[i]);
					TH2F *theWirePosHist = ((TH2F *)(HListWire.FindObject(getWirePosHist)));
					theWirePosHist->Fill(csctf_->lctChamber[i],csctf_->lctwireGroup[i]);

					if(csctf_->lctStation[i] == stationInput && csctf_->lctRing[i] == ringInput){
						TString gethist = TString::Format("Strip_BX_pos_ME%d%d_%d",(csctf_->lctStation[i]),(csctf_->lctRing[i]),(csctf_->lctChamber[i]));
						TH2F *theHist = ((TH2F *)(HListStrip.FindObject(gethist)));
						theHist->Fill(csctf_->lctstripNum[i],csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 1 && csctf_->lctRing[i] == 1){
						ME11_pos_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 1 && csctf_->lctRing[i] == 2){
						ME12_pos_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 1 && csctf_->lctRing[i] == 3){
						ME13_pos_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 2 && csctf_->lctRing[i] == 1){
						ME21_pos_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 2 && csctf_->lctRing[i] == 2){
						ME22_pos_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 3 && csctf_->lctRing[i] == 1){
						ME31_pos_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 3 && csctf_->lctRing[i] == 2){
						ME32_pos_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 4 && csctf_->lctRing[i] == 1){
						ME41_pos_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 4 && csctf_->lctRing[i] == 2){
						ME42_pos_BX->Fill(csctf_->lctBx[i]);
					}
				}
				//negative endcap only
				if(csctf_->lctEndcap[i] == -1){
					TString getWireNegHist = TString::Format("MEneg%d%d_WireGroup",csctf_->lctStation[i],csctf_->lctRing[i]);
					TH2F *theWireNegHist = ((TH2F *)(HListWire.FindObject(getWireNegHist)));
					theWireNegHist->Fill(csctf_->lctChamber[i],csctf_->lctwireGroup[i]);					

					if(csctf_->lctStation[i] == stationInput && csctf_->lctRing[i] == ringInput){
						TString gethist = TString::Format("Strip_BX_neg_ME%d%d_%d",(csctf_->lctStation[i]),(csctf_->lctRing[i]),(csctf_->lctChamber[i]));
						TH2F *theHist = ((TH2F *)(HListStrip.FindObject(gethist)));
						theHist->Fill(csctf_->lctstripNum[i],csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 1 && csctf_->lctRing[i] == 1){
						ME11_neg_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 1 && csctf_->lctRing[i] == 2){
						ME12_neg_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 1 && csctf_->lctRing[i] == 3){
						ME13_neg_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 2 && csctf_->lctRing[i] == 1){
						ME21_neg_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 2 && csctf_->lctRing[i] == 2){
						ME22_neg_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 3 && csctf_->lctRing[i] == 1){
						ME31_neg_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 3 && csctf_->lctRing[i] == 2){
						ME32_neg_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 4 && csctf_->lctRing[i] == 1){
						ME41_neg_BX->Fill(csctf_->lctBx[i]);
					}
					if(csctf_->lctStation[i] == 4 && csctf_->lctRing[i] == 2){
						ME42_neg_BX->Fill(csctf_->lctBx[i]);
					}
				}
				//BX by chamber 
				if(csctf_->lctEndcap[i] == -1){
					TString gethist = TString::Format("BX_negME_%d%d_%d_run%d",(csctf_->lctStation[i]),(csctf_->lctRing[i]),(csctf_->lctChamber[i]),(runNumber));
					TH1F *theHist = ((TH1F *)(HListneg.FindObject(gethist)));
					theHist->Fill(csctf_->lctBx[i]);
				}else{ //pos endcap
					TString gethist = TString::Format("BX_posME_%d%d_%d_run%d",(csctf_->lctStation[i]),(csctf_->lctRing[i]),(csctf_->lctChamber[i]),(runNumber));
					TH1F *theHist = ((TH1F *)(HList.FindObject(gethist)));
					theHist->Fill(csctf_->lctBx[i]);
					//BX by chamber with conditions
				}
				//BX by chamber 
				//TODO: signalProp must change from statOne to statTwo, etc. As necessary!			
				//cout << "Checkpoint 6.3" << endl;

			}//end LCT loop

			//loop for signal propagation time
			vector<int> tempVector;
			if(stationInput == 1){
				for(int m = 0; m < statOneIndex.size(); m++){
					tempVector.push_back( statOneIndex[m]);
				}  
			}
			if(stationInput == 2){
				for(int m = 0; m < statOneIndex.size(); m++){
					tempVector.push_back( statOneIndex[m]);
				}
			}
			if(stationInput == 3){
				for(int m = 0; m < statTwoIndex.size(); m++){
					tempVector.push_back( statTwoIndex[m]);
				}
			}
			if(stationInput == 4){
				for(int m = 0; m < statThreeIndex.size(); m++){
					tempVector.push_back( statThreeIndex[m]);
				}
			}

			//cout << "Check A" << endl;
			for(int q = 0; q < tempVector.size(); q++){
				//cout << "tempVector[q]  : " << tempVector[q] << endl;
				//cout << "statOneIndex[q]: " << statOneIndex[q] << endl;
			}
			for(int n = 0; n < tempVector.size(); n++){
				if(csctf_->lctEndcap[tempVector[n]] == 1){ 
					TString gethistprop = TString::Format("BX_Prop_posME_%d%d_%d_run%d",(csctf_->lctStation[n]),(csctf_->lctRing[n]),(csctf_->lctChamber[n]),(runNumber));
					//	cout << "the next line is gethistprop" << endl; 
					//cout << gethistprop << endl;

					TH1F *theHistprop = ((TH1F *)(HList_prop.FindObject(gethistprop)));
					/*
					 * cout << "Ring (temp)                              : " << csctf_->lctRing[tempVector[n]] << endl;
					 cout << "Ring (statOneIndex)                      : " << csctf_->lctRing[statOneIndex[n]] << endl;
					 cout << "Station (temp)                           : " << csctf_->lctStation[tempVector[n]] << endl;
					 cout << "Station (statOneIndex)                   : " << csctf_->lctStation[statOneIndex[n]] << endl;
					 cout << "Strip (temp)                             : " << csctf_->lctstripNum[tempVector[n]] << endl;
					 cout << "Strip (statOneIndex)                     : " << csctf_->lctstripNum[statOneIndex[n]] << endl;
					 cout << "Wire (temp)                              : " << csctf_->lctwireGroup[tempVector[n]] << endl;
					 cout << "Wire (statOneIndex)                      : " << csctf_->lctwireGroup[statOneIndex[n]] << endl;
					 cout << "Ring, Station, Strip, Wire (temp)        : " << csctf_->lctRing[tempVector[n]] << ", " << csctf_->lctStation[tempVector[n]] << ", " << csctf_->lctstripNum[tempVector[n]] << ", " << csctf_->lctwireGroup[tempVector[n]] << endl;
					 cout << "Ring, Station, Strip, Wire (statOneIndex): " << csctf_->lctRing[statOneIndex[n]] << ", " << csctf_->lctStation[statOneIndex[n]] << ", " << csctf_->lctstripNum[statOneIndex[n]] << ", " << csctf_->lctwireGroup[statOneIndex[n]] << endl;
					 cout << "signalPropTime (temp)                    : " << signalPropTime(csctf_->lctRing[tempVector[n]],csctf_->lctStation[tempVector[n]],csctf_->lctstripNum[tempVector[n]],csctf_->lctwireGroup[tempVector[n]]) << endl;
					 cout << "signalPropTime (statOneIndex)            : " << signalPropTime(csctf_->lctRing[statOneIndex[n]],csctf_->lctStation[statOneIndex[n]],csctf_->lctstripNum[statOneIndex[n]],csctf_->lctwireGroup[statOneIndex[n]]) << endl;
					//csctf_->lctBx[tempVector[i]]
					cout << "Now printing BXs" << endl;				

					for(int b = 0; b < csctf_->lctBx.size(); b++){
					cout << "BX " << b << " : " << csctf_->lctBx[b] << endl;
					}
					cout << " csctf_->lctBx.size(): " <<  csctf_->lctBx.size() << endl;

					cout << "BX (StatOneIndex)                        : " << csctf_->lctBx[statOneIndex[n]] << endl;
					cout << "BX (temp)                                : " << csctf_->lctBx[tempVector[n]] << endl;
					*/
					theHistprop->Fill(csctf_->lctBx[tempVector[n]] - signalPropTime(csctf_->lctRing[tempVector[n]],csctf_->lctStation[tempVector[n]],csctf_->lctstripNum[tempVector[n]],csctf_->lctwireGroup[tempVector[n]]) );

				}
				else{
					TString gethistprop = TString::Format("BX_Prop_negME_%d%d_%d_run%d",(csctf_->lctStation[n]),(csctf_->lctRing[n]),(csctf_->lctChamber[n]),(runNumber));
					//cout << "the next line is gethistprop" << endl; 
					//cout << gethistprop << endl;

					TH1F *theHistprop = ((TH1F *)(HList_propneg.FindObject(gethistprop)));
					/*				
									cout << "Ring (temp)                              : " << csctf_->lctRing[tempVector[n]] << endl;
									cout << "Ring (statOneIndex)                      : " << csctf_->lctRing[statOneIndex[n]] << endl;
									cout << "Station (temp)                           : " << csctf_->lctStation[tempVector[n]] << endl;
									cout << "Station (statOneIndex)                   : " << csctf_->lctStation[statOneIndex[n]] << endl;
									cout << "Strip (temp)                             : " << csctf_->lctstripNum[tempVector[n]] << endl;
									cout << "Strip (statOneIndex)                     : " << csctf_->lctstripNum[statOneIndex[n]] << endl;
									cout << "Wire (temp)                              : " << csctf_->lctwireGroup[tempVector[n]] << endl;
									cout << "Wire (statOneIndex)                      : " << csctf_->lctwireGroup[statOneIndex[n]] << endl;
									cout << "Ring, Station, Strip, Wire (temp)        : " << csctf_->lctRing[tempVector[n]] << ", " << csctf_->lctStation[tempVector[n]] << ", " << csctf_->lctstripNum[tempVector[n]] << ", " << csctf_->lctwireGroup[tempVector[n]] << endl;
									cout << "Ring, Station, Strip, Wire (statOneIndex): " << csctf_->lctRing[statOneIndex[n]] << ", " << csctf_->lctStation[statOneIndex[n]] << ", " << csctf_->lctstripNum[statOneIndex[n]] << ", " << csctf_->lctwireGroup[statOneIndex[n]] << endl;
									cout << "signalPropTime (temp)                    : " << signalPropTime(csctf_->lctRing[tempVector[n]],csctf_->lctStation[tempVector[n]],csctf_->lctstripNum[tempVector[n]],csctf_->lctwireGroup[tempVector[n]]) << endl;
									cout << "signalPropTime (statOneIndex)            : " << signalPropTime(csctf_->lctRing[statOneIndex[n]],csctf_->lctStation[statOneIndex[n]],csctf_->lctstripNum[statOneIndex[n]],csctf_->lctwireGroup[statOneIndex[n]]) << endl;
									cout << "BX (StatOneIndex)                        : " << csctf_->lctBx[statOneIndex[n]] << endl;
									cout << "BX (temp)                                : " << csctf_->lctBx[tempVector[n]] << endl;
									*/
					theHistprop->Fill(csctf_->lctBx[tempVector[n]] - signalPropTime(csctf_->lctRing[tempVector[n]],csctf_->lctStation[tempVector[n]],csctf_->lctstripNum[tempVector[n]],csctf_->lctwireGroup[tempVector[n]]) );

				}			

				/*
				   cout << "the next line is gethistprop" << endl;	
				   cout << gethistprop << endl;

				   TH1F *theHistprop = ((TH1F *)(HList_prop.FindObject(gethistprop)));
				   cout << "Ring (temp)                              : " << csctf_->lctRing[tempVector[n]] << endl;
				   cout << "Ring (statOneIndex)                      : " << csctf_->lctRing[statOneIndex[n]] << endl;
				   cout << "Station (temp)                           : " << csctf_->lctStation[tempVector[n]] << endl;
				   cout << "Station (statOneIndex)                   : " << csctf_->lctStation[statOneIndex[n]] << endl;
				   cout << "Strip (temp)                             : " << csctf_->lctstripNum[tempVector[n]] << endl;
				   cout << "Strip (statOneIndex)                     : " << csctf_->lctstripNum[statOneIndex[n]] << endl;
				   cout << "Wire (temp)                              : " << csctf_->lctwireGroup[tempVector[n]] << endl;
				   cout << "Wire (statOneIndex)                      : " << csctf_->lctwireGroup[statOneIndex[n]] << endl;
				   cout << "Ring, Station, Strip, Wire (temp)        : " << csctf_->lctRing[tempVector[n]] << ", " << csctf_->lctStation[tempVector[n]] << ", " << csctf_->lctstripNum[tempVector[n]] << ", " << csctf_->lctwireGroup[tempVector[n]] << endl;
				   cout << "Ring, Station, Strip, Wire (statOneIndex): " << csctf_->lctRing[statOneIndex[n]] << ", " << csctf_->lctStation[statOneIndex[n]] << ", " << csctf_->lctstripNum[statOneIndex[n]] << ", " << csctf_->lctwireGroup[statOneIndex[n]] << endl;
				   cout << "signalPropTime (temp)                    : " << signalPropTime(csctf_->lctRing[tempVector[n]],csctf_->lctStation[tempVector[n]],csctf_->lctstripNum[tempVector[n]],csctf_->lctwireGroup[tempVector[n]]) << endl;
				   cout << "signalPropTime (statOneIndex)            : " << signalPropTime(csctf_->lctRing[statOneIndex[n]],csctf_->lctStation[statOneIndex[n]],csctf_->lctstripNum[statOneIndex[n]],csctf_->lctwireGroup[statOneIndex[n]]) << endl;

				   theHistprop->Fill(csctf_->lctBx[tempVector[n]] - signalPropTime(csctf_->lctRing[tempVector[n]],csctf_->lctStation[tempVector[n]],csctf_->lctstripNum[tempVector[n]],csctf_->lctwireGroup[tempVector[n]]) );
				   */
			}
			//cout << "Check B" << endl;
			//BX by chamber with conditions

			//}
			//}

			//}//end LCT loop
			//cout << "Checkpoint 7" << endl;
			/*}}}*/
//Run Pieter's TOF code with my event selection
/*{{{*/

double deltaY=0;
double deltaT1a=0;
double deltaT2a=0;
double deltaT1b=0;
double deltaT2b=0;
double deltaT3a=0;
double deltaT3b=0;
double Segx1=0;
double Segy1=0;
double Segz1=0;
double Segx2=0;
double Segy2=0;
double Segz2=0;
double Segx3=0;
double Segy3=0;
double Segz3=0;
int ALCT_Bx2=10;
bool saveInfo=false;
bool hastwo=false;
bool hasthree=false;
int TOFchamber = 0;
int TOFendcap = 0;
int TOFendcap2 = 0;
int TOFendcap3 = 0;		
double correctedBX = 0;
int bottom = 0;

vector< vector<double> > NeighborTOFcorrectedBX;

//BAM: Added requirement to ensure same endcap 11/7/2014

//Corrected TOF algorithm 11/11/2014
for(int i = 0; i < statOneIndex.size(); i++){
	//Event topology.  This is the same topology as later, the only difference is that here, I do not require the 2 hits in ME11 to be neighbors.  In which case, TOF will be calculated for each point.	
	if(EventSelection==0){
		if(statOneIndex.size() == 0     || statTwoIndex.size() == 0 || statThreeIndex.size() == 0)continue;
		if(statOneIndexPos.size() > 2   || statOneIndexNeg.size() > 2) continue;
		if(statTwoIndexPos.size() > 1   || statTwoIndexNeg.size() > 1) continue;
		if(statThreeIndexPos.size() > 1 || statThreeIndexNeg.size() > 1) continue;	
	}
	if(EventSelection == 1){
		if(statOneIndex.size() == 0     || statTwoIndex.size() == 0 || statThreeIndex.size() == 0)continue;
		if(statOneIndexPos.size() > 2   || statOneIndexNeg.size() > 2) continue;
		if(statTwoIndexPos.size() > 1   || statTwoIndexNeg.size() > 1) continue;
		if(statThreeIndexPos.size() > 1 || statThreeIndexNeg.size() > 1) continue;
		if(neighborList.size() == 0) continue;
	}

	if((TOFendcap2 == 0 && TOFendcap3 == 0) || (TOFendcap2 != 0 && TOFendcap2 == csctf_->lctEndcap[statOneIndex[i]]) || (TOFendcap3 !=0 && TOFendcap3 == csctf_->lctEndcap[statOneIndex[i]])){
		//////cout << "LCT In Station 1 (for TOF)" << endl;					
		//////cout << "Checkpoint 1" << endl;
		saveInfo=true;
		ALCT_Bx2=csctf_->lctBx[statOneIndex[i]];
		Segx1=getx(csctf_->lctglobalEta[statOneIndex[i]],csctf_->lctglobalPhi[statOneIndex[i]],csctf_->lctStation[statOneIndex[i]],csctf_->lctRing[statOneIndex[i]],csctf_->lctSector[statOneIndex[i]]);
		Segy1=gety(csctf_->lctglobalEta[statOneIndex[i]],csctf_->lctglobalPhi[statOneIndex[i]],csctf_->lctStation[statOneIndex[i]],csctf_->lctRing[statOneIndex[i]],csctf_->lctSector[statOneIndex[i]]);
		Segz1=getz(csctf_->lctglobalEta[statOneIndex[i]],csctf_->lctglobalPhi[statOneIndex[i]],csctf_->lctStation[statOneIndex[i]],csctf_->lctRing[statOneIndex[i]],csctf_->lctEndcap[statOneIndex[i]], csctf_->lctChamber[statOneIndex[i]]);
		TOFchamber = csctf_->lctChamber[statOneIndex[i]];                        
		TOFendcap = csctf_->lctEndcap[statOneIndex[i]];
		//////cout << "Checkpoint 2" << endl;

	}

	for(int j = 0; j < statTwoIndex.size(); j++){
		//////cout << "Checkpoint 3" << endl;
		if((TOFendcap == 0 && TOFendcap3 == 0) || (TOFendcap != 0 && TOFendcap == csctf_->lctEndcap[statTwoIndex[j]]) || (TOFendcap3 !=0 && TOFendcap3 == csctf_->lctEndcap[statTwoIndex[j]])){
			//////cout << "Checkpoint 4" << endl;
			//////cout << "LCT In Station 2 (for TOF)" << endl;				

			//Make sure that the triggering LCT has 5 <= BX <= 7
			//For LR 228354 the triggering LCT is ME2
			//if(csctf_->lctBx[statTwoIndex[j]] < 5 || csctf_->lctBx[statTwoIndex[j]] > 7) cout << "The event should be skipped. ME2 BX: " << csctf_->lctBx[statTwoIndex[j]] << endl;
			//if((csctf_->lctBx[statTwoIndex[j]] < 5 || csctf_->lctBx[statTwoIndex[j]] > 7) && ME2Trigger == 1) continue;
			if(TOFendcap == 1 && statTwoIndexPos.size() == 1 && statThreeIndexPos.size() == 1){
				if((csctf_->lctBx[statTwoIndexPos[0]] < 5 || csctf_->lctBx[statTwoIndexPos[0]] > 7 ) && ME2Trigger == 1) continue;
			}else{
				if(TOFendcap == -1 && statTwoIndexNeg.size() == 1 && statThreeIndexNeg.size() == 1){
					if((csctf_->lctBx[statTwoIndexNeg[0]] < 5 || csctf_->lctBx[statTwoIndexNeg[0]] > 7 ) && ME2Trigger == 1) continue;
				}else{
					continue;
				}
			}

			hastwo=true;
			Segx2=getx(csctf_->lctglobalEta[statTwoIndex[j]],csctf_->lctglobalPhi[statTwoIndex[j]],csctf_->lctStation[statTwoIndex[j]],csctf_->lctRing[statTwoIndex[j]],csctf_->lctSector[statTwoIndex[j]]);
			Segy2=gety(csctf_->lctglobalEta[statTwoIndex[j]],csctf_->lctglobalPhi[statTwoIndex[j]],csctf_->lctStation[statTwoIndex[j]],csctf_->lctRing[statTwoIndex[j]],csctf_->lctSector[statTwoIndex[j]]);
			Segz2=getz(csctf_->lctglobalEta[statTwoIndex[j]],csctf_->lctglobalPhi[statTwoIndex[j]],csctf_->lctStation[statTwoIndex[j]],csctf_->lctRing[statTwoIndex[j]],csctf_->lctEndcap[statTwoIndex[j]], csctf_->lctChamber[statTwoIndex[j]]);

			TOFendcap2 = csctf_->lctEndcap[statTwoIndex[j]];	
			//////cout << "Checkpoint 5" << endl;
		}

		for(int k =0; k < statThreeIndex.size(); k++){
			//////cout << "Checkpoint 6" << endl;
			if((TOFendcap2 == 0 && TOFendcap == 0) || (TOFendcap2 != 0 && TOFendcap2 == csctf_->lctEndcap[statThreeIndex[k]]) || (TOFendcap !=0 && TOFendcap == csctf_->lctEndcap[statThreeIndex[k]])){				
				//////cout << "Checkpoint 7" << endl;
				//////cout << "LCT In Station 3 (for TOF)" << endl;				
				//if((csctf_->lctBx[statThreeIndex[k]] < 5 || csctf_->lctBx[statThreeIndex[k]] > 7) && ME3Trigger == 1) continue;

				if(TOFendcap == 1 && statTwoIndexPos.size() == 1 && statThreeIndexPos.size() == 1){
					if((csctf_->lctBx[statThreeIndexPos[0]] < 5 || csctf_->lctBx[statThreeIndexPos[0]] > 7 ) && ME3Trigger == 1) continue;
				}else{
					if(TOFendcap == -1 && statTwoIndexNeg.size() == 1 && statThreeIndexNeg.size() == 1){
						if((csctf_->lctBx[statThreeIndexNeg[0]] < 5 || csctf_->lctBx[statThreeIndexNeg[0]] > 7 ) && ME3Trigger == 1) continue;
					}else{
						continue;
					}
				}

				hasthree=true;
				Segx3=getx(csctf_->lctglobalEta[statThreeIndex[k]],csctf_->lctglobalPhi[statThreeIndex[k]],csctf_->lctStation[statThreeIndex[k]],csctf_->lctRing[statThreeIndex[k]],csctf_->lctSector[statThreeIndex[k]]);
				Segy3=gety(csctf_->lctglobalEta[statThreeIndex[k]],csctf_->lctglobalPhi[statThreeIndex[k]],csctf_->lctStation[statThreeIndex[k]],csctf_->lctRing[statThreeIndex[k]],csctf_->lctSector[statThreeIndex[k]]);
				Segz3=getz(csctf_->lctglobalEta[statThreeIndex[k]],csctf_->lctglobalPhi[statThreeIndex[k]],csctf_->lctStation[statThreeIndex[k]],csctf_->lctRing[statThreeIndex[k]],csctf_->lctEndcap[statThreeIndex[k]], csctf_->lctChamber[statThreeIndex[k]]);
				TOFendcap3 = csctf_->lctEndcap[statThreeIndex[k]];			
				//////cout << "Checkpoint 8" << endl;
			}			

			//Now move on to part two of TOF.  Use the information gathered above:
			//////cout << "Checkpoint 9" << endl;

			if(saveInfo && hastwo && hasthree) {
				//////cout << "Checkpoint 10" << endl;

				deltaY=Segy3-Segy2;
				deltaT1a=(sqrt(Segx2*Segx2+Segy2*Segy2+Segz2*Segz2)-sqrt(Segx1*Segx1+Segy1*Segy1+Segz1*Segz1))/29.98/25;
				deltaT1b=(sqrt((Segx2-Segx1)*(Segx2-Segx1)+(Segy2-Segy1)*(Segy2-Segy1)+(Segz2-Segz1)*(Segz2-Segz1)))/29.98/25;

				deltaT2a=(sqrt(Segx3*Segx3+Segy3*Segy3+Segz3*Segz3)-sqrt(Segx1*Segx1+Segy1*Segy1+Segz1*Segz1))/29.98/25;
				deltaT2b=(sqrt((Segx3-Segx1)*(Segx3-Segx1)+(Segy3-Segy1)*(Segy3-Segy1)+(Segz3-Segz1)*(Segz3-Segz1)))/29.98/25;

				deltaT3a=(sqrt((Segx3+Segx2)/2*(Segx3+Segx2)/2+(Segy3+Segy2)/2*(Segy3+Segy2)/2+(Segz3+Segz2)/2*(Segz3+Segz2)/2)-sqrt(Segx1*Segx1+Segy1*Segy1+Segz1*Segz1))/29.98/25;
				deltaT3b=(sqrt(((Segx3+Segx2)/2-Segx1)*((Segx3+Segx2)/2-Segx1)+((Segy3+Segy2)/2-Segy1)*((Segy3+Segy2)/2-Segy1)+((Segz3+Segz2)/2-Segz1)*((Segz3+Segz2)/2-Segz1)))/29.98/25;


				//Find incident angle of muon
				double sideA = 0;
				double sideB = 0;
				double incAngle = 0;
				if(deltaY >= 0){
					sideA = dist(0,Segy3, Segz3, Segx3, Segy3, Segz3);
					sideB = dist(0,Segy1, Segz1, 0, Segy3, Segz3);
					incAngle = atan(sideA/sideB);
				}else{
					sideA = dist(0,Segy2, Segz2, Segx2, Segy2, Segz2);
					sideB = dist(0,Segy1, Segz1, 0, Segy2, Segz2);
					incAngle = atan(sideA/sideB);
				}
				//////cout << "Checkpoint 11" << endl;

				//I should be able to put the Propogation time here:
				double propCorr = signalPropTime(csctf_->lctRing[statOneIndex[i]],csctf_->lctStation[statOneIndex[i]],csctf_->lctstripNum[statOneIndex[i]],csctf_->lctwireGroup[statOneIndex[i]]);
				if(TOFendcap == -1){
					TString PropCorrectionstringB = TString::Format("PropCorrection_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
					TH1F *thePropHistB = ((TH1F *)(HListPropCorrection.FindObject(PropCorrectionstringB)));
					thePropHistB->Fill(propCorr);
				}else{
					TString PropCorrectionstringA = TString::Format("PropCorrection_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
					TH1F *thePropHistA = ((TH1F *)(HListPropCorrection.FindObject(PropCorrectionstringA)));
					thePropHistA->Fill(propCorr);
				}
				///
				if(stationInput == 1){
					//Make incident angle plots and TOF correction plots
					/*{{{*/
					if(TOFendcap == -1){
						if((csctf_->lctBx[statTwoIndexNeg[0]] < 5 || csctf_->lctBx[statTwoIndexNeg[0]] > 7) && ME2Trigger == 1) cout << "There is a problem!" << endl;
						if((csctf_->lctBx[statThreeIndexNeg[0]] < 5 || csctf_->lctBx[statThreeIndexNeg[0]] > 7) && ME3Trigger == 1) cout << "There is a problem!" << endl;
						TString incAngleStringA = TString::Format("IncidentAngle_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
						TH1F *theAngleHistA = ((TH1F *)(HListIncidentAngle.FindObject(incAngleStringA)));
						theAngleHistA->Fill(incAngle);			
						//////cout << "Checkpoint 12" << endl;
						if(deltaY > 100 && deltaY < 1500 && ALCT_Bx2 > 0){
							//////cout << "Checkpoint 13" << endl;
							if(forwardOnly == 1) continue;										
							correctedBX = ALCT_Bx2-(deltaT2a+deltaT2b)/10; //The factor of 10 is there so that we are in mm and ns
							TString TOFstringA = TString::Format("TOF_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theHistA = ((TH1F *)(HListTOF.FindObject(TOFstringA)));
							theHistA->Fill(correctedBX);
							//cout << "Checkpoint 3" << endl;

							TString gethist = TString::Format("Strip_TOF_BX_neg_ME%d%d_%d",(stationInput),(ringInput),(TOFchamber));
							TH2F *theHist = ((TH2F *)(HListStrip_TOF.FindObject(gethist)));
							theHist->Fill(csctf_->lctstripNum[statOneIndex[i]],correctedBX);
							//These are the backwards muons
							backwardsCounter++;
							totalCounter++;									
						}        
						if(deltaY<-100 && deltaY>-1500 && ALCT_Bx2>0){
							//////cout << "Checkpoint 14" << endl;
							correctedBX = ALCT_Bx2-(deltaT1a-deltaT1b)/10;
							TString TOFstringA = TString::Format("TOF_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theHistA = ((TH1F *)(HListTOF.FindObject(TOFstringA)));
							theHistA->Fill(correctedBX);
							//cout << "Checkpoint 4" << endl;

							TString gethist = TString::Format("Strip_TOF_BX_neg_ME%d%d_%d",(stationInput),(ringInput),(TOFchamber));
							TH2F *theHist = ((TH2F *)(HListStrip_TOF.FindObject(gethist)));
							theHist->Fill(csctf_->lctstripNum[statOneIndex[i]],correctedBX);			
							totalCounter++;									
						}
					}else{
						if((csctf_->lctBx[statTwoIndexPos[0]] < 5 || csctf_->lctBx[statTwoIndexPos[0]] > 7) && ME2Trigger == 1) cout << "There is a problem!" << endl;
						if((csctf_->lctBx[statThreeIndexPos[0]] < 5 || csctf_->lctBx[statThreeIndexPos[0]] > 7) && ME3Trigger == 1) cout << "There is a problem!" << endl;
						TString incAngleStringA = TString::Format("IncidentAngle_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
						TH1F *theAngleHistA = ((TH1F *)(HListIncidentAngle.FindObject(incAngleStringA)));
						theAngleHistA->Fill(incAngle);
						if(deltaY > 100 && deltaY < 1500 && ALCT_Bx2 > 0){
							//////cout << "Checkpoint 15" << endl;
							if(forwardOnly == 1) continue;										
							correctedBX = ALCT_Bx2-(deltaT2a+deltaT2b)/10; //The factor of 10 is there so that we are in mm and ns
							TString TOFstringA = TString::Format("TOF_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theHistA = ((TH1F *)(HListTOF.FindObject(TOFstringA)));
							theHistA->Fill(correctedBX);
							//cout << "Checkpoint 5" << endl;

							TString gethist = TString::Format("Strip_TOF_BX_pos_ME%d%d_%d",(stationInput),(ringInput),(TOFchamber));
							TH2F *theHist = ((TH2F *)(HListStrip_TOF.FindObject(gethist)));
							theHist->Fill(csctf_->lctstripNum[statOneIndex[i]],correctedBX);
							//These are backwards muons
							backwardsCounter++;
							totalCounter++;									
						} 
						if(deltaY<-100 && deltaY>-1500 && ALCT_Bx2>0){
							//////cout << "Checkpoint 16" << endl;
							correctedBX = ALCT_Bx2-(deltaT1a-deltaT1b)/10;
							TString TOFstringA = TString::Format("TOF_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theHistA = ((TH1F *)(HListTOF.FindObject(TOFstringA)));
							theHistA->Fill(correctedBX);            
							//cout << "Checkpoint 6" << endl;

							TString gethist = TString::Format("Strip_TOF_BX_pos_ME%d%d_%d",(stationInput),(ringInput),(TOFchamber));
							TH2F *theHist = ((TH2F *)(HListStrip_TOF.FindObject(gethist)));
							theHist->Fill(csctf_->lctstripNum[statOneIndex[i]],correctedBX); 
							totalCounter++;									
						}
					}
					/*}}}*/
					//////cout << "Checkpoint 17" << endl;
					//Now only for the events that also have neighboring LCTs in ME11.  I am removing the tolerance on deltaY
					if(neighborList.size() > 0 && vecvec_contains(neighborList, statOneIndex[i])){
						//Since this is only for neighbors, I want to impliment my event topology. 
						//Require that there are only two (2) LCTs in ME11. 1 LCT in ME2 and 1 LCT in ME3
						if(statOneIndexPos.size() > 2   || statOneIndexNeg.size() > 2) continue;
						if(statTwoIndexPos.size() > 1   || statTwoIndexNeg.size() > 1) continue;
						if(statThreeIndexPos.size() > 1 || statThreeIndexNeg.size() > 1) continue;									

						//////cout << "Checkpoint 18" << endl;
						tofForNeighborsPlotCounter++;			
						NeighborTOFcounter++;	
						////cout << "TOFforNeighbors_BX is being filled with corrected BX for ME11_" << TOFchamber << endl;
						////cout << "" << endl;					
						//cout << "Filling TOFforNeighbors_BX plots" << endl;

						if(TOFendcap == -1){
							TString incAngleStringA = TString::Format("IncidentAngleNeighbors_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theAngleHistA = ((TH1F *)(HListIncidentAngleNeighbors.FindObject(incAngleStringA)));
							theAngleHistA->Fill(incAngle);
							//////cout << "Checkpoint 19" << endl;
							if(deltaY > 100 && deltaY < 1500 && ALCT_Bx2 > 0){
								//////cout << "Checkpoint 20" << endl;
								if(forwardOnly == 1) continue;										
								correctedBX = ALCT_Bx2-(deltaT2a+deltaT2b)/10; //The factor of 10 is there so that we are in mm and ns
								double correction = (deltaT2a+deltaT2b)/10;
								//cout << "Loop Alpha" << endl;
								//cout << "Chamber Number: " << TOFchamber << endl;						
								TString TOFforNeighborsstringA = TString::Format("TOFforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistA = ((TH1F *)(HListTOFforNeighbors.FindObject(TOFforNeighborsstringA)));
								theHistA->Fill(correctedBX);
								TString TOFCorrectionstringA = TString::Format("TOFCorrection_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistB = ((TH1F *)(HListTOFCorrection.FindObject(TOFCorrectionstringA)));
								theHistB->Fill(correction);			
								vector<double> ntcb;
								ntcb.push_back(statOneIndex[i]);
								ntcb.push_back(correctedBX);
								NeighborTOFcorrectedBX.push_back(ntcb);			

								TString TOF_And_PropforNeighborsstringA = TString::Format("TOF_And_PropforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistPropA = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(TOF_And_PropforNeighborsstringA)));
								theHistPropA->Fill(correctedBX - propCorr);							


							}
							if(deltaY<-100 && deltaY>-1500 && ALCT_Bx2>0){
								//////cout << "Checkpoint 21" << endl;
								correctedBX = ALCT_Bx2-(deltaT1a-deltaT1b)/10;
								double correction = (deltaT1a-deltaT1b)/10;					
								//cout << "Loop Beta" << endl;
								//cout << "Chamber Number: " << TOFchamber << endl;						
								TString TOFforNeighborsstringA = TString::Format("TOFforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistA = ((TH1F *)(HListTOFforNeighbors.FindObject(TOFforNeighborsstringA)));
								theHistA->Fill(correctedBX);
								TString TOFCorrectionstringA = TString::Format("TOFCorrection_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistB = ((TH1F *)(HListTOFCorrection.FindObject(TOFCorrectionstringA)));
								theHistB->Fill(correction);					
								vector<double> ntcb;
								ntcb.push_back(statOneIndex[i]);
								ntcb.push_back(correctedBX);
								NeighborTOFcorrectedBX.push_back(ntcb);

								TString TOF_And_PropforNeighborsstringA = TString::Format("TOF_And_PropforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistPropA = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(TOF_And_PropforNeighborsstringA)));
								theHistPropA->Fill(correctedBX - propCorr);
							}
						}else{
							TString incAngleStringA = TString::Format("IncidentAngleNeighbors_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theAngleHistA = ((TH1F *)(HListIncidentAngleNeighbors.FindObject(incAngleStringA)));
							theAngleHistA->Fill(incAngle);
							if(deltaY > 100 && deltaY < 1500 && ALCT_Bx2 > 0){
								if(forwardOnly == 1) continue;										
								//////cout << "Checkpoint 22" << endl;
								correctedBX = ALCT_Bx2-(deltaT2a+deltaT2b)/10; 
								double correction = (deltaT2a+deltaT2b)/10;						
								//cout << "Loop Gamma" << endl;
								//cout << "Chamber Number: " << TOFchamber << endl;						
								TString TOFforNeighborsstringA = TString::Format("TOFforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistA = ((TH1F *)(HListTOFforNeighbors.FindObject(TOFforNeighborsstringA)));
								theHistA->Fill(correctedBX);
								TString TOFCorrectionstringA = TString::Format("TOFCorrection_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistB = ((TH1F *)(HListTOFCorrection.FindObject(TOFCorrectionstringA)));
								theHistB->Fill(correction);
								vector<double> ntcb;
								ntcb.push_back(statOneIndex[i]);
								ntcb.push_back(correctedBX);
								NeighborTOFcorrectedBX.push_back(ntcb);					

								TString TOF_And_PropforNeighborsstringA = TString::Format("TOF_And_PropforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistPropA = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(TOF_And_PropforNeighborsstringA)));
								theHistPropA->Fill(correctedBX - propCorr);
							}
							if(deltaY<-100 && deltaY>-1500 && ALCT_Bx2>0){
								correctedBX = ALCT_Bx2-(deltaT1a-deltaT1b)/10;
								//////cout << "Checkpoint 23" << endl;
								double correction = (deltaT1a-deltaT1b)/10;						
								//cout << "Loop Delta" << endl;
								//cout << "Chamber Number: " << TOFchamber << endl;						
								TString TOFforNeighborsstringA = TString::Format("TOFforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistA = ((TH1F *)(HListTOFforNeighbors.FindObject(TOFforNeighborsstringA)));
								theHistA->Fill(correctedBX);
								TString TOFCorrectionstringA = TString::Format("TOFCorrection_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistB = ((TH1F *)(HListTOFCorrection.FindObject(TOFCorrectionstringA)));
								theHistB->Fill(correction);
								vector<double> ntcb;
								ntcb.push_back(statOneIndex[i]);
								ntcb.push_back(correctedBX);
								NeighborTOFcorrectedBX.push_back(ntcb);	

								TString TOF_And_PropforNeighborsstringA = TString::Format("TOF_And_PropforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistPropA = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(TOF_And_PropforNeighborsstringA)));
								theHistPropA->Fill(correctedBX);				
							}
						}
						//cout << "" << endl;
					}
				} //close ME1 loop

				if(stationInput == 4){
					if(TOFendcap == -1){
						if((csctf_->lctBx[statTwoIndexNeg[0]] < 5 || csctf_->lctBx[statTwoIndexNeg[0]] > 7) && ME2Trigger == 1) cout << "There is a problem!" << endl;
						if((csctf_->lctBx[statThreeIndexNeg[0]] < 5 || csctf_->lctBx[statThreeIndexNeg[0]] > 7) && ME3Trigger == 1) cout << "There is a problem!" << endl;								

						TString incAngleStringA = TString::Format("IncidentAngle_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
						TH1F *theAngleHistA = ((TH1F *)(HListIncidentAngle.FindObject(incAngleStringA)));
						theAngleHistA->Fill(incAngle);			
						//////cout << "Checkpoint 12" << endl;
						if(deltaY > 100 && deltaY < 1500 && ALCT_Bx2 > 0){
							if(forwardOnly == 1) continue;										
							//////cout << "Checkpoint 13" << endl;
							correctedBX = ALCT_Bx2-(deltaT2a-deltaT2b)/10; //The factor of 10 is there so that we are in mm and ns
							TString TOFstringA = TString::Format("TOF_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theHistA = ((TH1F *)(HListTOF.FindObject(TOFstringA)));
							theHistA->Fill(correctedBX);

							TString gethist = TString::Format("Strip_TOF_BX_neg_ME%d%d_%d",(stationInput),(ringInput),(TOFchamber));
							TH2F *theHist = ((TH2F *)(HListStrip_TOF.FindObject(gethist)));
							theHist->Fill(csctf_->lctstripNum[statOneIndex[i]],correctedBX);
						}        
						if(deltaY<-100 && deltaY>-1500 && ALCT_Bx2>0){
							//////cout << "Checkpoint 14" << endl;
							correctedBX = ALCT_Bx2-(deltaT1a+deltaT1b)/10;
							TString TOFstringA = TString::Format("TOF_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theHistA = ((TH1F *)(HListTOF.FindObject(TOFstringA)));
							theHistA->Fill(correctedBX);  

							TString gethist = TString::Format("Strip_TOF_BX_neg_ME%d%d_%d",(stationInput),(ringInput),(TOFchamber));
							TH2F *theHist = ((TH2F *)(HListStrip_TOF.FindObject(gethist)));
							theHist->Fill(csctf_->lctstripNum[statOneIndex[i]],correctedBX);			
						}
					}else{
						if((csctf_->lctBx[statTwoIndexPos[0]] < 5 || csctf_->lctBx[statTwoIndexPos[0]] > 7) && ME2Trigger == 1) cout << "There is a problem!" << endl;
						if((csctf_->lctBx[statThreeIndexPos[0]] < 5 || csctf_->lctBx[statThreeIndexPos[0]] > 7) && ME3Trigger == 1) cout << "There is a problem!" << endl;								

						TString incAngleStringA = TString::Format("IncidentAngle_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
						TH1F *theAngleHistA = ((TH1F *)(HListIncidentAngle.FindObject(incAngleStringA)));
						theAngleHistA->Fill(incAngle);
						if(deltaY > 100 && deltaY < 1500 && ALCT_Bx2 > 0){
							if(forwardOnly == 1) continue;										
							//////cout << "Checkpoint 15" << endl;
							correctedBX = ALCT_Bx2-(deltaT2a-deltaT2b)/10; //The factor of 10 is there so that we are in mm and ns
							TString TOFstringA = TString::Format("TOF_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theHistA = ((TH1F *)(HListTOF.FindObject(TOFstringA)));
							theHistA->Fill(correctedBX);

							TString gethist = TString::Format("Strip_TOF_BX_pos_ME%d%d_%d",(stationInput),(ringInput),(TOFchamber));
							TH2F *theHist = ((TH2F *)(HListStrip_TOF.FindObject(gethist)));
							theHist->Fill(csctf_->lctstripNum[statOneIndex[i]],correctedBX);
						} 
						if(deltaY<-100 && deltaY>-1500 && ALCT_Bx2>0){
							//////cout << "Checkpoint 16" << endl;
							correctedBX = ALCT_Bx2+(deltaT1a-deltaT1b)/10;
							TString TOFstringA = TString::Format("TOF_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theHistA = ((TH1F *)(HListTOF.FindObject(TOFstringA)));
							theHistA->Fill(correctedBX);             

							TString gethist = TString::Format("Strip_TOF_BX_pos_ME%d%d_%d",(stationInput),(ringInput),(TOFchamber));
							TH2F *theHist = ((TH2F *)(HListStrip_TOF.FindObject(gethist)));
							theHist->Fill(csctf_->lctstripNum[statOneIndex[i]],correctedBX);
						}
					}
					//////cout << "Checkpoint 17" << endl;
					//Now only for the events that also have neighboring LCTs in ME11.  I am removing the tolerance on deltaY
					if(neighborList.size() > 0 && vecvec_contains(neighborList, statOneIndex[i])){
						//Since this is only for neighbors, I want to impliment my event topology. 
						//Require that there are only two (2) LCTs in ME11. 1 LCT in ME2 and 1 LCT in ME3
						if(statOneIndexPos.size() > 2   || statOneIndexNeg.size() > 2) continue;
						if(statTwoIndexPos.size() > 1   || statTwoIndexNeg.size() > 1) continue;
						if(statThreeIndexPos.size() > 1 || statThreeIndexNeg.size() > 1) continue;									

						//////cout << "Checkpoint 18" << endl;
						tofForNeighborsPlotCounter++;			
						NeighborTOFcounter++;	
						////cout << "TOFforNeighbors_BX is being filled with corrected BX for ME11_" << TOFchamber << endl;
						////cout << "" << endl;					
						//cout << "Filling TOFforNeighbors_BX plots" << endl;

						if(TOFendcap == -1){
							TString incAngleStringA = TString::Format("IncidentAngleNeighbors_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theAngleHistA = ((TH1F *)(HListIncidentAngleNeighbors.FindObject(incAngleStringA)));
							theAngleHistA->Fill(incAngle);
							//////cout << "Checkpoint 19" << endl;
							if(deltaY > 100 && deltaY < 1500 && ALCT_Bx2 > 0){
								if(forwardOnly == 1) continue;										
								//////cout << "Checkpoint 20" << endl;
								correctedBX = ALCT_Bx2-(deltaT2a-deltaT2b)/10; //The factor of 10 is there so that we are in mm and ns
								double correction = (deltaT2a-deltaT2b)/10;
								//cout << "Loop Alpha" << endl;
								//cout << "Chamber Number: " << TOFchamber << endl;						
								TString TOFforNeighborsstringA = TString::Format("TOFforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistA = ((TH1F *)(HListTOFforNeighbors.FindObject(TOFforNeighborsstringA)));
								theHistA->Fill(correctedBX);
								TString TOFCorrectionstringA = TString::Format("TOFCorrection_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistB = ((TH1F *)(HListTOFCorrection.FindObject(TOFCorrectionstringA)));
								theHistB->Fill(correction);			
								vector<double> ntcb;
								ntcb.push_back(statOneIndex[i]);
								ntcb.push_back(correctedBX);
								NeighborTOFcorrectedBX.push_back(ntcb);										

								TString TOF_And_PropforNeighborsstringA = TString::Format("TOF_And_PropforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistPropA = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(TOF_And_PropforNeighborsstringA)));
								theHistPropA->Fill(correctedBX - propCorr);
							}
							if(deltaY<-100 && deltaY>-1500 && ALCT_Bx2>0){
								//////cout << "Checkpoint 21" << endl;
								correctedBX = ALCT_Bx2-(deltaT1a+deltaT1b)/10;
								double correction = (deltaT1a+deltaT1b)/10;					
								//cout << "Loop Beta" << endl;
								//cout << "Chamber Number: " << TOFchamber << endl;						
								TString TOFforNeighborsstringA = TString::Format("TOFforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistA = ((TH1F *)(HListTOFforNeighbors.FindObject(TOFforNeighborsstringA)));
								theHistA->Fill(correctedBX);
								TString TOFCorrectionstringA = TString::Format("TOFCorrection_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistB = ((TH1F *)(HListTOFCorrection.FindObject(TOFCorrectionstringA)));
								theHistB->Fill(correction);					
								vector<double> ntcb;
								ntcb.push_back(statOneIndex[i]);
								ntcb.push_back(correctedBX);
								NeighborTOFcorrectedBX.push_back(ntcb);

								TString TOF_And_PropforNeighborsstringA = TString::Format("TOF_And_PropforNeighbors_BX_neg_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistPropA = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(TOF_And_PropforNeighborsstringA)));
								theHistPropA->Fill(correctedBX - propCorr);
							}
						}else{
							TString incAngleStringA = TString::Format("IncidentAngleNeighbors_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
							TH1F *theAngleHistA = ((TH1F *)(HListIncidentAngleNeighbors.FindObject(incAngleStringA)));
							theAngleHistA->Fill(incAngle);
							if(deltaY > 100 && deltaY < 1500 && ALCT_Bx2 > 0){
								if(forwardOnly == 1) continue;										
								//////cout << "Checkpoint 22" << endl;
								correctedBX = ALCT_Bx2-(deltaT2a-deltaT2b)/10; 
								double correction = (deltaT2a-deltaT2b)/10;						
								//cout << "Loop Gamma" << endl;
								//cout << "Chamber Number: " << TOFchamber << endl;						
								TString TOFforNeighborsstringA = TString::Format("TOFforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistA = ((TH1F *)(HListTOFforNeighbors.FindObject(TOFforNeighborsstringA)));
								theHistA->Fill(correctedBX);
								TString TOFCorrectionstringA = TString::Format("TOFCorrection_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistB = ((TH1F *)(HListTOFCorrection.FindObject(TOFCorrectionstringA)));
								theHistB->Fill(correction);
								vector<double> ntcb;
								ntcb.push_back(statOneIndex[i]);
								ntcb.push_back(correctedBX);
								NeighborTOFcorrectedBX.push_back(ntcb);	

								TString TOF_And_PropforNeighborsstringA = TString::Format("TOF_And_PropforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistPropA = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(TOF_And_PropforNeighborsstringA)));
								theHistPropA->Fill(correctedBX - propCorr);				
							}
							if(deltaY<-100 && deltaY>-1500 && ALCT_Bx2>0){
								correctedBX = ALCT_Bx2-(deltaT1a+deltaT1b)/10;
								//////cout << "Checkpoint 23" << endl;
								double correction = (deltaT1a+deltaT1b)/10;						
								//cout << "Loop Delta" << endl;
								//cout << "Chamber Number: " << TOFchamber << endl;						
								TString TOFforNeighborsstringA = TString::Format("TOFforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistA = ((TH1F *)(HListTOFforNeighbors.FindObject(TOFforNeighborsstringA)));
								theHistA->Fill(correctedBX);
								TString TOFCorrectionstringA = TString::Format("TOFCorrection_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistB = ((TH1F *)(HListTOFCorrection.FindObject(TOFCorrectionstringA)));
								theHistB->Fill(correction);
								vector<double> ntcb;
								ntcb.push_back(statOneIndex[i]);
								ntcb.push_back(correctedBX);
								NeighborTOFcorrectedBX.push_back(ntcb);	

								TString TOF_And_PropforNeighborsstringA = TString::Format("TOF_And_PropforNeighbors_BX_pos_ME%d%d_%d",stationInput,ringInput,TOFchamber);
								TH1F *theHistPropA = ((TH1F *)(HListTOF_And_PropforNeighbors.FindObject(TOF_And_PropforNeighborsstringA)));
								theHistPropA->Fill(correctedBX - propCorr);				
							}
						}
						//cout << "" << endl;
					}
				} //close ME4 loop

				//////cout << "Checkpoint 24" << endl;
				deltaY=0;
				deltaT1a=0;
				deltaT2a=0;
				deltaT1b=0;
				deltaT2b=0;
				deltaT3a=0;
				deltaT3b=0;
				correctedBX = 0;


			}
			//////cout << "Checkpoint 25" << endl;
			Segx3=0;
			Segy3=0;
			Segz3=0;
			hasthree=false;
			TOFendcap3 = 0; 
		}
		//////cout << "Checkpoint 26" << endl;
		Segx2=0;
		Segy2=0;
		Segz2=0;
		hastwo=false;
		TOFendcap2 = 0;
	}
	//////cout << "Checkpoint 27" << endl;
	Segx1=0;
	Segy1=0;
	Segz1=0;
	ALCT_Bx2=0;
	TOFchamber = 0;
	TOFendcap = 0;
	saveInfo=false;
}
/*}}}*/
//Now input the TOF corrected Neighbor values for the walk around
/*{{{*/
int neighborFillerCounter = 0;				
for(int c = 0; c < NeighborTOFcorrectedBX.size(); c++){
	//This should only fill once per event with my current topology
	if(neighborFillerCounter > 0) continue;
	chiAndtofPlotCounter++;
	waTOFcounter++;

	vector<int> info = neighborToChamber(neighborList, NeighborTOFcorrectedBX[c][0]);
	int lowerIndex = 0;
	int upperIndex = 0;

	if(info[1] == 1){
		lowerIndex = NeighborTOFcorrectedBX[c][0];
		upperIndex = info[0];
	}
	if(info[1] == 2){
		lowerIndex = info[0];
		upperIndex = NeighborTOFcorrectedBX[c][0];
	}

	double lowerBX;
	double upperBX;
	double lowerBXProp, upperBXProp;

	if(  PileUp50ns == 1 ){
		if(upperBX > 7 || upperBX < 5) continue;
		if(lowerBX > 7 || lowerBX < 5) continue;
	}
	for(int k =0; k < NeighborTOFcorrectedBX.size(); k++){
		if(NeighborTOFcorrectedBX[k][0] == lowerIndex) lowerBX = NeighborTOFcorrectedBX[k][1];
		if(NeighborTOFcorrectedBX[k][0] == upperIndex) upperBX = NeighborTOFcorrectedBX[k][1];	
	}


	double relBX = lowerBX - upperBX;

	double lowerBXPropCorr = signalPropTime(csctf_->lctRing[lowerIndex],csctf_->lctStation[lowerIndex],csctf_->lctstripNum[lowerIndex],csctf_->lctwireGroup[lowerIndex]);
	double upperBXPropCorr = signalPropTime(csctf_->lctRing[upperIndex],csctf_->lctStation[upperIndex],csctf_->lctstripNum[upperIndex],csctf_->lctwireGroup[upperIndex]);

	double lowprop  = lowerBX - lowerBXPropCorr;
	double highprop = upperBX - upperBXPropCorr;
	double relBXProp = lowprop - highprop; 

	if(lowerBXPropCorr < 0){
		cout << "THERE IS A BIG PROBLEM HERE: THE MIN PROP TIME < 0" << endl;

		cout << "Chamber, Ring, Station, Strip, Wire: " << csctf_->lctChamber[lowerIndex] << ", " << csctf_->lctRing[lowerIndex] << ", " << csctf_->lctStation[lowerIndex] << ", " << csctf_->lctstripNum[lowerIndex] << ", " << csctf_->lctwireGroup[lowerIndex] << endl;
		cout << "Prop Time: " << lowerBXPropCorr << endl;
	}

	if(lowerBXPropCorr < minPropTime) minPropTime = lowerBXPropCorr;
	if(lowerBXPropCorr > maxPropTime) maxPropTime = lowerBXPropCorr;
	if(upperBXPropCorr > maxPropTime) maxPropTime = upperBXPropCorr;
	if(upperBXPropCorr < minPropTime) minPropTime = upperBXPropCorr;

	//print some info
	if(csctf_->lctChamber[upperIndex] == 19 || csctf_->lctChamber[lowerIndex] == 19){				
		if(csctf_->lctEndcap[lowerIndex] == 1){
			cout << "" << endl;
			cout << "ME" << csctf_->lctStation[lowerIndex] << csctf_->lctRing[lowerIndex] << "/" << csctf_->lctChamber[lowerIndex] << " to ME" << csctf_->lctStation[upperIndex] << csctf_->lctRing[upperIndex] << "/" << csctf_->lctChamber[upperIndex] << " in endcap: " << csctf_->lctEndcap[lowerIndex] << endl;
			cout << "BX for upper: " << csctf_->lctBx[upperIndex] << endl;
			cout << "BX for lower: " << csctf_->lctBx[lowerIndex] << endl; 
			cout << "" << endl;
			cout << "ME2 and ME3 BX: " << csctf_->lctBx[statTwoIndexPos[0]] << ", " << csctf_->lctBx[statThreeIndexPos[0]] << endl;
		}

		cout << "" << endl;
	}
	TString gethist = TString::Format("WalkAroundwithTOF_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex], csctf_->lctEndcap[lowerIndex]);
	////cout << gethist << endl;					
	TH1F *theHist = ((TH1F *)(WalkAroundwithTOF.FindObject(gethist)));
	theHist->Fill(relBX);
	////cout << "Filled Plot: " << gethist << endl;

	TString gethistA = TString::Format("WalkAroundwithTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex], csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctEndcap[lowerIndex]);
	////cout << gethistA << endl;			
	TH1F *theHistA = ((TH1F *)(WalkAroundwithTOF.FindObject(gethistA)));
	theHistA->Fill(lowerBX);
	////cout << "Filled Plot: " << gethistA << endl;
	TString gethistB = TString::Format("WalkAroundwithTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex],csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[upperIndex], csctf_->lctEndcap[lowerIndex]);
	////cout << gethistB << endl;					
	TH1F *theHistB = ((TH1F *)(WalkAroundwithTOF.FindObject(gethistB)));
	theHistB->Fill(upperBX);

	TString getPropHist1 = TString::Format("WalkAroundwithTOF_And_Prop_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex], csctf_->lctEndcap[lowerIndex]);
	TH1F *thePropHist1 = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(getPropHist1)));
	thePropHist1->Fill(relBXProp);

	TString getPropHistA1 = TString::Format("WalkAroundwithTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex], csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctEndcap[lowerIndex]);
	TH1F *thePropHistA1 = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(getPropHistA1)));
	thePropHistA1->Fill(lowprop);
	TString getPropHistB1 = TString::Format("WalkAroundwithTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex],csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[upperIndex], csctf_->lctEndcap[lowerIndex]);
	TH1F *thePropHistB1 = ((TH1F *)(WalkAroundwithTOF_And_Prop.FindObject(getPropHistB1)));
	thePropHistB1->Fill(highprop);


	//********//
	TString getStripHist = TString::Format("NeighborStripTOF_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex], csctf_->lctEndcap[lowerIndex]);
	TH2F *theStripHist = ((TH2F *)(NeighborStripTOF.FindObject(getStripHist)));
	theStripHist->Fill((csctf_->lctstripNum[lowerIndex]-csctf_->lctstripNum[upperIndex]),relBX);

	TString getStripHistA = TString::Format("NeighborStripTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex], csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctEndcap[lowerIndex]);
	TH2F *theStripHistA = ((TH2F *)(NeighborStripTOF.FindObject(getStripHistA)));
	theStripHistA->Fill(csctf_->lctstripNum[lowerIndex],lowerBX);

	TString getStripHistB = TString::Format("NeighborStripTOF_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex],csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[upperIndex], csctf_->lctEndcap[lowerIndex]);
	TH2F *theStripHistB = ((TH2F *)(NeighborStripTOF.FindObject(getStripHistB)));
	theStripHistB->Fill(csctf_->lctstripNum[upperIndex],upperBX);

	neighborFillerCounter++;

	TString getPropHist = TString::Format("NeighborStripTOF_And_Prop_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex], csctf_->lctEndcap[lowerIndex]);
	TH2F *thePropHist = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(getPropHist)));
	thePropHist->Fill((csctf_->lctstripNum[lowerIndex]-csctf_->lctstripNum[upperIndex]),relBXProp);

	//cout << "Ring, Station, Strip, Wire: " << csctf_->lctRing[upperIndex]<< ", " << csctf_->lctStation[upperIndex]<< ", " << csctf_->lctstripNum[upperIndex]<< ", " << csctf_->lctwireGroup[upperIndex] << endl;
	//cout << "signalProp: " << signalPropTime(csctf_->lctRing[upperIndex],csctf_->lctStation[upperIndex],csctf_->lctstripNum[upperIndex],csctf_->lctwireGroup[upperIndex]) << endl;
	//cout << " upperBXProp, uppBXProp-upp

	TString getPropHistA = TString::Format("NeighborStripTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex], csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctEndcap[lowerIndex]);
	TH2F *thePropHistA = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(getPropHistA)));
	thePropHistA->Fill(csctf_->lctstripNum[lowerIndex],lowerBX-lowerBXPropCorr);

	TString getPropHistB = TString::Format("NeighborStripTOF_And_Prop_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[lowerIndex], csctf_->lctStation[upperIndex], csctf_->lctRing[upperIndex], csctf_->lctChamber[upperIndex],csctf_->lctStation[lowerIndex], csctf_->lctRing[lowerIndex], csctf_->lctChamber[upperIndex], csctf_->lctEndcap[lowerIndex]);
	TH2F *thePropHistB = ((TH2F *)(NeighborStripTOF_And_Prop.FindObject(getPropHistB)));
	thePropHistB->Fill(csctf_->lctstripNum[upperIndex],upperBX-upperBXPropCorr);

}
/*}}}*/
//Run the WalkAround without any TOF requirements
/*{{{*/
//CheckpointA
//At this point I have collected all sorts of data and have a good understanding of these events.
//The next step is to look at synchronizing the chambers.
//This will be done via a "walk-around" for each ring and each station. 
//Then a similar method will be applied to neighboring chambers from different stations.
//It is important to note that all stations/rings/chambers are already in time except for ME1/1 and ME4/2, so
//All other chambers may be used as a reference.

//I will begin by just looking for LCTs in neighboring chambers.
//Duplicates must also be taken in to account so as not to skew the statistics
//TODO It may be better to limit the number of LCTs in the event to 2 or 3 to ensure that the LCTs came from the same Muon.
//It also may be worth it to look at LCTs within the same track (assuming that they are in neighboring chambers)			

//LCTs are being double counted. Need to create a unique identifier and make sure that each LCT is only used once.


//quick loop to see how many unique LCTs are in ME11		
int numLCTsME11 = 0;
for(int i = 0; i<csctf_->lctBx.size(); i++){
	//if(vector_contains(duplicateList, i) == true) continue; 
	if(contains(duplicates, i) == true) continue; 
	if(csctf_->lctStation[i] != 1) continue; 
	if(csctf_->lctRing[i] != 1) continue;
	numLCTsME11++;
}

vector< vector<int> > usedLCTpairs_Ring;
for(int i=0; i<csctf_->lctBx.size(); i++){
	for(int j = 0; j < csctf_->lctBx.size(); j++){
		if(i == j) continue;
		if(contains(duplicates, i) == true) continue; //Make sure not to count the duplicates
		//if(vector_contains(duplicateList, j) == true) continue; //Make sure not to count the duplicates
		if(csctf_->lctEndcap[i] != csctf_->lctEndcap[j]) continue; //Make sure that the LCTs are in the same endcap					
		if(csctf_->lctStation[i] != csctf_->lctStation[j]) continue; //Make sure that the LCTs are in the same station
		if(csctf_->lctRing[i] != csctf_->lctRing[j]) continue; // Make sure that the LCTs are in the same ring
		if(csctf_->lctChamber[i] == csctf_->lctChamber[j]) continue; //Make sure that the LCTs are NOT in the same chamber	

		//if(csctf_->lctQuality[j] <= 10 || csctf_->lctQuality[i] <= 10) continue;

		//topology implimentation
		if(EventSelection == 1){						
			if(statOneIndex.size() == 0     || statTwoIndex.size() == 0 || statThreeIndex.size() == 0)continue;
			if(statOneIndexPos.size() > 2   || statOneIndexNeg.size() > 2) continue;
			if(statTwoIndexPos.size() > 1   || statTwoIndexNeg.size() > 1) continue;
			if(statThreeIndexPos.size() > 1 || statThreeIndexNeg.size() > 1) continue;
			if(neighborList.size() == 0) continue;						
		}
		//LCTCheck
		//if(numLCTsME11 != 2) continue;				
		int upper = 0;
		int lower = 0;
		int upperInt = 0;
		int lowerInt = 0;

		vector< int > LCTpair;
		LCTpair.push_back(i);		
		LCTpair.push_back(j);
		if(vectorvector_contains(usedLCTpairs_Ring, i, j) == true) continue;

		if(csctf_->lctChamber[i] > csctf_->lctChamber[j]){
			upper = csctf_->lctChamber[i]; 
			lower = csctf_->lctChamber[j];
			upperInt = i;
			lowerInt = j;					
		}else{
			upper = csctf_->lctChamber[j];
			lower = csctf_->lctChamber[i];
			upperInt = j;
			lowerInt = i;

		}

		if(csctf_->lctStation[i] == 1 || csctf_->lctRing[i] == 2){ //These chambers go from 1-36
			if((upper == lower+1) || (upper == 36 && lower == 1)){
				///
				//Adding the ME2,3 Trigger requirement to check endcap
				if(csctf_->lctEndcap[i] == 1){
					if(statTwoIndexPos.size() == 1){
						if((csctf_->lctBx[statTwoIndexPos[0]] < 5 || csctf_->lctBx[statTwoIndexPos[0]] > 7) && ME2Trigger == 1) continue;
					}else{
						continue;
					}
					if(statThreeIndexPos.size() == 1){
						if((csctf_->lctBx[statThreeIndexPos[0]] < 5 || csctf_->lctBx[statThreeIndexPos[0]] > 7) && ME3Trigger == 1) continue;
					}else{
						continue;
					}
				}else{
					if(statTwoIndexNeg.size() == 1){
						if((csctf_->lctBx[statTwoIndexNeg[0]] < 5 || csctf_->lctBx[statTwoIndexNeg[0]] > 7) && ME2Trigger == 1) continue;
					}else{
						continue;
					}
					if(statThreeIndexNeg.size() == 1){
						if((csctf_->lctBx[statThreeIndexNeg[0]] < 5 || csctf_->lctBx[statThreeIndexNeg[0]] > 7) && ME3Trigger == 1) continue;
					}else{
						continue;
					}
				}							
				///

				usedLCTpairs_Ring.push_back(LCTpair);							
				int relBX = csctf_->lctBx[lowerInt] - csctf_->lctBx[upperInt];
				//Checkpoint
				chiPlotCounter++;

				//////
				double lowerBX = csctf_->lctBx[lowerInt];
				double upperBX = csctf_->lctBx[upperInt];
				if(  PileUp50ns == 1 ){
					if(upperBX > 7 || upperBX < 5) continue;					
					if(lowerBX > 7 || lowerBX < 5) continue;
				}
				/*
				   if(upperBX == 8 || upperBX == 4){
				   cout << "There is a problem, here! B" << endl;
				   cout << "upperBX, lowerBX: " << upperBX << ", " << lowerBX << endl;
				   if(upperBX <= 7) cout << "upperBX <= 7" << endl;				
				   if(upperBX >= 5) cout << "upperBX >= 5" << endl;				
				   if(lowerBX <= 7) cout << "lowerBX <= 7" << endl;				
				   if(lowerBX <= 5) cout << "lowerBX >= 5" << endl;				
				   }
				   if(lowerBX == 8 || lowerBX == 4){
				   cout << "There is a problem, here! C" << endl;
				   cout << "upperBX, lowerBX: " << upperBX << ", " << lowerBX << endl;
				   if(upperBX <= 7) cout << "upperBX <= 7" << endl;                
				   if(upperBX >= 5) cout << "upperBX >= 5" << endl;
				   if(lowerBX <= 7) cout << "lowerBX <= 7" << endl;
				   if(lowerBX <= 5) cout << "lowerBX >= 5" << endl;
				   }
				//This is where it should not get to...
				*/
				double lowerBXPropCorr = signalPropTime(csctf_->lctRing[lowerInt],csctf_->lctStation[lowerInt],csctf_->lctstripNum[lowerInt],csctf_->lctwireGroup[lowerInt]);
				double upperBXPropCorr = signalPropTime(csctf_->lctRing[upperInt],csctf_->lctStation[upperInt],csctf_->lctstripNum[upperInt],csctf_->lctwireGroup[upperInt]);

				double lowprop  = lowerBX - lowerBXPropCorr;
				double highprop = upperBX - upperBXPropCorr;
				double relBXProp = lowprop - highprop;

				if(lowerBXPropCorr < minPropTime) minPropTime = lowerBXPropCorr;
				if(lowerBXPropCorr > maxPropTime) maxPropTime = lowerBXPropCorr;
				if(upperBXPropCorr > maxPropTime) maxPropTime = upperBXPropCorr;
				if(upperBXPropCorr < minPropTime) minPropTime = upperBXPropCorr;
				/////

				if(upper == lower+1){
					TString gethist = TString::Format(	           "WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctEndcap[i]); 
					TH1F *theHist = ((TH1F *)(WalkAround.FindObject(gethist)));
					//////cout << "gethist, theHist, relBX: " << gethist << ", " << theHist << ", " << relBX << endl;
					theHist->Fill(relBX);
					TString gethistA = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctEndcap[i]);
					TH1F *theHistA = ((TH1F *)(WalkAround.FindObject(gethistA)));
					theHistA->Fill(csctf_->lctBx[lowerInt]);
					TString gethistB = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctEndcap[i]);
					TH1F *theHistB = ((TH1F *)(WalkAround.FindObject(gethistB)));
					theHistB->Fill(csctf_->lctBx[upperInt]);


					///prop

					TString getPropHist1 = TString::Format(                 "WalkAroundProp_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctEndcap[i]);
					TH1F *thePropHist1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHist1)));
					thePropHist1->Fill(relBXProp);

					TString getPropHistA1 = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctEndcap[i]);
					TH1F *thePropHistA1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHistA1)));
					thePropHistA1->Fill(lowprop);
					TString getPropHistB1 = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctEndcap[i]);
					TH1F *thePropHistB1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHistB1)));
					thePropHistB1->Fill(highprop);

				}else{
					TString gethist = TString::Format(                 "WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctEndcap[i]); 
					TH1F *theHist = ((TH1F *)(WalkAround.FindObject(gethist)));
					theHist->Fill(relBX);

					TString gethistA = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctEndcap[i]);
					TH1F *theHistA = ((TH1F *)(WalkAround.FindObject(gethistA)));
					theHistA->Fill(csctf_->lctBx[lowerInt]);
					TString gethistB = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctEndcap[i]);
					TH1F *theHistB = ((TH1F *)(WalkAround.FindObject(gethistB)));
					theHistB->Fill(csctf_->lctBx[upperInt]);

					///prop

					TString getPropHist1 = TString::Format(                 "WalkAroundProp_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctEndcap[i]);
					TH1F *thePropHist1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHist1)));
					thePropHist1->Fill(relBXProp);
					TString getPropHistA1 = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctEndcap[i]);
					TH1F *thePropHistA1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHistA1)));
					thePropHistA1->Fill(lowprop);
					TString getPropHistB1 = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctEndcap[i]);
					TH1F *thePropHistB1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHistB1)));
					thePropHistB1->Fill(highprop);


				}	
			}
		}else{ //These chambers go from 1-18
			if((upper == lower+1) || (upper == 18 && lower == 1)){ 
				///
				//Adding the ME2,3 Trigger requirement to check endcap
				if(csctf_->lctEndcap[i] == 1){
					if(statTwoIndexPos.size() == 1){
						if((csctf_->lctBx[statTwoIndexPos[0]] < 5 || csctf_->lctBx[statTwoIndexPos[0]] > 7) && ME2Trigger == 1) continue;
					}else{
						continue;
					}
					if(statThreeIndexPos.size() == 1){
						if((csctf_->lctBx[statThreeIndexPos[0]] < 5 || csctf_->lctBx[statThreeIndexPos[0]] > 7) && ME3Trigger == 1) continue;
					}else{
						continue;
					}
				}else{
					if(statTwoIndexNeg.size() == 1){
						if((csctf_->lctBx[statTwoIndexNeg[0]] < 5 || csctf_->lctBx[statTwoIndexNeg[0]] > 7) && ME2Trigger == 1) continue;
					}else{
						continue;
					}
					if(statThreeIndexNeg.size() == 1){
						if((csctf_->lctBx[statThreeIndexNeg[0]] < 5 || csctf_->lctBx[statThreeIndexNeg[0]] > 7) && ME3Trigger == 1) continue;
					}else{
						continue;
					}
				}
				///					
				usedLCTpairs_Ring.push_back(LCTpair);							
				int relBX = csctf_->lctBx[i] - csctf_->lctBx[j];

				//
				double lowerBX = csctf_->lctBx[lowerInt];
				double upperBX = csctf_->lctBx[upperInt];
				if(  PileUp50ns == 1 ){
					if(upperBX > 7 || upperBX < 5) continue;
					if(lowerBX > 7 || lowerBX < 5) continue;
				}
				double lowerBXPropCorr = signalPropTime(csctf_->lctRing[lowerInt],csctf_->lctStation[lowerInt],csctf_->lctstripNum[lowerInt],csctf_->lctwireGroup[lowerInt]);
				double upperBXPropCorr = signalPropTime(csctf_->lctRing[upperInt],csctf_->lctStation[upperInt],csctf_->lctstripNum[upperInt],csctf_->lctwireGroup[upperInt]);

				double lowprop  = lowerBX - lowerBXPropCorr;
				double highprop = upperBX - upperBXPropCorr;
				double relBXProp = lowprop - highprop;
				//
				if(lowerBXPropCorr < minPropTime) minPropTime = lowerBXPropCorr;
				if(lowerBXPropCorr > maxPropTime) maxPropTime = lowerBXPropCorr;
				if(upperBXPropCorr > maxPropTime) maxPropTime = upperBXPropCorr;
				if(upperBXPropCorr < minPropTime) minPropTime = upperBXPropCorr;
				if(upper == lower+1){
					TString gethist = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctEndcap[i]); 
					TH1F *theHist = ((TH1F *)(WalkAround.FindObject(gethist)));
					//////cout << "gethist, theHist, relBX: " << gethist << ", " << theHist << ", " << relBX << endl;
					theHist->Fill(relBX);
					TString gethistA = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctEndcap[i]);
					TH1F *theHistA = ((TH1F *)(WalkAround.FindObject(gethistA)));
					theHistA->Fill(csctf_->lctBx[lowerInt]);
					TString gethistB = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctEndcap[i]);
					TH1F *theHistB = ((TH1F *)(WalkAround.FindObject(gethistB)));
					theHistB->Fill(csctf_->lctBx[upperInt]);							

					///prop

					TString getPropHist1 = TString::Format(                 "WalkAroundProp_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctEndcap[i]);
					TH1F *thePropHist1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHist1)));
					thePropHist1->Fill(relBXProp);

					TString getPropHistA1 = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctEndcap[i]);
					TH1F *thePropHistA1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHistA1)));
					thePropHistA1->Fill(lowprop);
					TString getPropHistB1 = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctStation[j], csctf_->lctRing[j], upper, csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctEndcap[i]);
					TH1F *thePropHistB1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHistB1)));
					thePropHistB1->Fill(highprop);

				}else{
					TString gethist = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctEndcap[i]);                            
					TH1F *theHist = ((TH1F *)(WalkAround.FindObject(gethist)));
					theHist->Fill(relBX);
					TString gethistA = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctEndcap[i]);
					TH1F *theHistA = ((TH1F *)(WalkAround.FindObject(gethistA)));
					theHistA->Fill(csctf_->lctBx[lowerInt]);
					TString gethistB = TString::Format("WalkAround_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctEndcap[i]);
					TH1F *theHistB = ((TH1F *)(WalkAround.FindObject(gethistB)));
					theHistB->Fill(csctf_->lctBx[upperInt]);


					///prop

					TString getPropHist1 = TString::Format(                 "WalkAroundProp_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctEndcap[i]);
					TH1F *thePropHist1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHist1)));
					thePropHist1->Fill(relBXProp);
					TString getPropHistA1 = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctStation[i], csctf_->lctRing[i], lower, csctf_->lctEndcap[i]);
					TH1F *thePropHistA1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHistA1)));
					thePropHistA1->Fill(lowprop);
					TString getPropHistB1 = TString::Format("WalkAroundProp_Rings_ME%d%d_%d_and_ME%d%d_%d_Overlap_BX_ME%d%d_%d_Endcap_%d",csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctStation[j], csctf_->lctRing[j], lower, csctf_->lctStation[i], csctf_->lctRing[i], upper, csctf_->lctEndcap[i]);
					TH1F *thePropHistB1 = ((TH1F *)(WalkAroundProp.FindObject(getPropHistB1)));
					thePropHistB1->Fill(highprop);

				}
			}
		}
	}
}
/*}}}*/
//commented out WalkAround Disk and Endcap Code
/*{{{*/
/*
//Now I will look for LCTs in neighboring rings.  For now, the requirement will be strict.  
//Requiring that the chambers in concentric rings match exactly. 
//This restriction may be loosened later.
vector< vector<int> > usedLCTpairs_Station;
for(int i=0; i<csctf_->lctBx.size(); i++){
for(int j = 0; j < csctf_->lctBx.size(); j++){
if(i == j) continue;
if(vector_contains(duplicateList, i) == true) continue; //Make sure not to count the duplicates
if(vector_contains(duplicateList, j) == true) continue; //Make sure not to count the duplicates
if(csctf_->lctEndcap[i] != csctf_->lctEndcap[j]) continue; //Make sure that the LCTs are in the same endcap
if(csctf_->lctStation[i] != csctf_->lctStation[j]) continue; //Make sure that the LCTs are in the same station
if(csctf_->lctRing[i] == csctf_->lctRing[j]) continue; // Make sure that the LCTs are NOT the same ring
int outerRing = 0;
int innerRing = 0;
int innerChamber = 0;
int outerChamber = 0;

vector< int > LCTpair;
LCTpair.push_back(i);
LCTpair.push_back(j);
if(vectorvector_contains(usedLCTpairs_Station, i, j) == true) continue;

if(csctf_->lctRing[i] > csctf_->lctRing[j]){
outerRing = csctf_->lctRing[i];
innerRing = csctf_->lctRing[j];

innerChamber = csctf_->lctChamber[j];
outerChamber = csctf_->lctChamber[i];
}else{
outerRing = csctf_->lctRing[j];
innerRing = csctf_->lctRing[i];

outerChamber = csctf_->lctChamber[j];
innerChamber = csctf_->lctChamber[i];					
}
if(csctf_->lctStation[i] == 1){
if(innerRing == 1 && outerRing == 2){
if(csctf_->lctChamber[i] == csctf_->lctChamber[j]){
usedLCTpairs_Station.push_back(LCTpair);								
int relBX = csctf_->lctBx[i] - csctf_->lctBx[j];
TString gethist = TString::Format("WalkAround_Stations_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",1, innerRing, csctf_->lctChamber[i], 1, outerRing, csctf_->lctChamber[i], csctf_->lctEndcap[i]);	
TH1F *theHist = ((TH1F *)(WalkAround2.FindObject(gethist)));
theHist->Fill(relBX);						
}
}
if(innerRing == 2 && outerRing == 3){
if(csctf_->lctChamber[i] == csctf_->lctChamber[j]){
usedLCTpairs_Station.push_back(LCTpair);								
int relBX = csctf_->lctBx[i] - csctf_->lctBx[j];
TString gethist = TString::Format("WalkAround_Stations_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",1, innerRing, csctf_->lctChamber[i], 1, outerRing, csctf_->lctChamber[i], csctf_->lctEndcap[i]); 
TH1F *theHist = ((TH1F *)(WalkAround2.FindObject(gethist)));
theHist->Fill(relBX);                                           
}
}
}else{
int TwoChamber = innerChamber*2;
if((outerChamber == TwoChamber) || (outerChamber == TwoChamber - 1)){
usedLCTpairs_Station.push_back(LCTpair);							
int relBX = csctf_->lctBx[i] - csctf_->lctBx[j];
////cout << "InnerChamber, TwoChamber,  OuterChamber, OuterChamber-1: " << innerChamber << ", " << TwoChamber << ", " << outerChamber << ", " << outerChamber-1 << endl;
TString gethist = TString::Format("WalkAround_Stations_RelBX_ME%d%d_%d_to_ME%d%d_%d_or_%d_Endcap_%d",csctf_->lctStation[i], innerRing, innerChamber, csctf_->lctStation[j], outerRing, TwoChamber-1, TwoChamber, csctf_->lctEndcap[i]);
TH1F *theHist = ((TH1F *)(WalkAround2.FindObject(gethist)));
////cout << gethist << endl;
theHist->Fill(relBX);
}
}

}
}
*/

//Now walk around Disks

/*
   vector< vector<int> > usedLCTpairs_Disks;
   for(int i = 0; i < csctf_->lctBx.size(); i++){
   for(int j = 0; j < csctf_->lctBx.size(); j++){
   if(i == j) continue;
   if(vector_contains(duplicateList, i) == true) continue; //Make sure not to count the duplicates
   if(vector_contains(duplicateList, j) == true) continue; //Make sure not to count the duplicates
   if(csctf_->lctStation[i] == csctf_->lctStation[j]) continue; //different stations
   if(csctf_->lctSector[i] != csctf_->lctSector[j]) continue; //same sector

   vector< int > LCTpair;
   LCTpair.push_back(i);
   LCTpair.push_back(j);
   if(vectorvector_contains(usedLCTpairs_Disks, i, j) == true) continue;                        

   int lowInt;
   int highInt;
   if(csctf_->lctStation[i] < csctf_->lctStation[j]){
   lowInt = i;
   highInt = j;
   }else{
   lowInt = j;
   highInt = i;
   }
   if(csctf_->lctStation[lowInt]+1 != csctf_->lctStation[highInt]) continue; //make sure that the stations are neighboring
   usedLCTpairs_Disks.push_back(LCTpair);                                                       
   int relBX = csctf_->lctBx[i] - csctf_->lctBx[j];

   TString gethist = TString::Format("WalkAround_Disks_RelBX_Sector_%d_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",csctf_->lctSector[i],csctf_->lctStation[lowInt],csctf_->lctRing[lowInt],csctf_->lctChamber[lowInt],csctf_->lctStation[highInt],csctf_->lctRing[highInt],csctf_->lctChamber[highInt],csctf_->lctEndcap[i]);
   TH1F *theHist = ((TH1F *)(WalkAroundDisks.FindObject(gethist)));
   theHist->Fill(relBX);
   }
   }

*/
/*}}}*/
}//Event loop
}//Event Number if Statement
//Print final counters and write to root output file
/*{{{*/
//////cout << "There were " << me11trkCounterTotal << " tracks with an LCT from ME1/1." << endl;	
//////cout << "There were " << dupeCounterTotal << " events with at least one duplicate LCT." << endl;	
////cout << "(Event Counters) tofCounter, chiCounter, chiAndtofCounter: " << tofCounter << ", " << chiCounter << ", " << chiAndtofCounter << endl;
//TOF for neighbor and chi&TOF counters should be the same.
////cout << "(Plot Counters) NieghborPlotCounter, tofForNeighborsCounter, chiAndtofPlotCounter: " << chiPlotCounter << ", " << tofForNeighborsPlotCounter << ", " << chiAndtofPlotCounter << endl;
// Keep track of trigger count for failed muons
cout << "" << endl;
cout << "" << endl;
cout << "" << endl;
cout << "" << endl;
cout << "Min and Max Prop Time: " << minPropTime << ", " << maxPropTime << endl;
cout << "" << endl;;	

cout << "Number of Backwards Muons: " << backwardsCounter << endl;
cout << "Total number of muons    : " << totalCounter << endl;	
HList.Write();
HListTOF.Write();
HListneg.Write();
HList_prop.Write();
HList_propneg.Write();
WalkAround.Write();
WalkAroundProp.Write();
HListTOFforNeighbors.Write();
WalkAroundwithTOF.Write();
HListTOFCorrection.Write();
HListIncidentAngle.Write();
HListIncidentAngleNeighbors.Write();
HListStrip.Write();	
HListStrip_TOF.Write();	
NeighborStripTOF.Write();	
HListWire.Write();
NeighborStripTOF_And_Prop.Write();
HListPropCorrection.Write();
WalkAroundwithTOF_And_Prop.Write();
HListTOF_And_PropforNeighbors.Write();
HListxCoordinate.Write();
HListyCoordinate.Write();
HListzCoordinate.Write();	
BAM->Write();
/*}}}*/
} 

