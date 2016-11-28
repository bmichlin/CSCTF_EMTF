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
void MakeChiSquare(TString filename, int runNumber, int station, int ring, int numChambersInRing, int endcap){
TFile *tfile = TFile::Open(filename);
for(int chamber = 1; chamber < 37; chamber++){
		int highChamber = 0; //These chambers are generic and can be used for any of the plots
		int lowChamber = 0;
		if(chamber == numChambersInRing){
			int lowChamber =  chamber;
			int highChamber = 1;
		}else{
			int lowChamber = chamber;
			int highChamber = chamber+1;
		}
		TString waRelBX = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,endcap);
		if(chamber != numChambersInRing){
			cout << "pow(  ( " << Form("x%d",chamber-1) << " - " << Form("x%d",chamber) << " + " << getHist(tfile, waRelBX)->GetMean() << " )  / " << getHist(tfile, waRelBX)->GetMeanError() << " ,2)  +  " << endl;   
//cout << "Chamber Number, RelBX, Error on RelBX: " << chamber << ", " << getHist(tfile, waRelBX)->GetMean() << ", " << getHist(tfile, waRelBX)->GetMeanError() << endl;
		}else{
		 TString waRelBX = TString::Format("WalkAround_Rings_RelBX_ME%d%d_%d_to_ME%d%d_%d_Endcap_%d",station,ring,lowChamber,station,ring,highChamber,endcap);
		 cout << "pow(  ( " << Form("x%d",chamber-1) << " - " << "x0" << " + " << getHist(tfile, waRelBX)->GetMean() << " )  / " << getHist(tfile, waRelBX)->GetMeanError() << " ,2);  " << endl; 
//cout << "Chamber Number, RelBX, Error on RelBX: " << chamber << ", " << getHist(tfile, waRelBX)->GetMean() << ", " << getHist(tfile, waRelBX)->GetMeanError() << endl;
		}

	}
}
