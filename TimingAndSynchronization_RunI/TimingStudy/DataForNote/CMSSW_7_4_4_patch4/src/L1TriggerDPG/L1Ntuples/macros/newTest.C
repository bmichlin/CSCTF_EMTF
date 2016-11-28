void newTest(){

	//TString arr[] = {"root://eoscms///store/user/bmichlin/Cosmics/Complete_DataSets_For_Synchronization_Note/150715_113113/0000/CRAFT_CERN_Complete_163.root","root://eoscms///store/user/bmichlin/Cosmics/Complete_DataSets_For_Synchronization_Note/150715_113113/0000/CRAFT_CERN_Complete_164.root"};

	//for(int i = 0; i < arr.size(); i++){
	//TString command = TString::Format("WalkAround *macro = new WalkAround(\"%s\")");


	gROOT->ProcessLine(".x initL1Analysis.C+");
	gROOT->ProcessLine(".L New_output_WalkAround.C+");
	gROOT->ProcessLine("WalkAround *macro = new WalkAround(\"root://eoscms///store/user/bmichlin/Cosmics/Complete_DataSets_For_Synchronization_Note/150715_113113/0000/CRAFT_CERN_Complete_163.root\")");	
	gROOT->ProcessLine("macro->run(1,666666,-1,0,1,1,36,1,1,0)");
	gROOT->ProcessLine("WalkAround *macro2 = new WalkAround(\"root://eoscms///store/user/bmichlin/Cosmics/Complete_DataSets_For_Synchronization_Note/150715_113113/0000/CRAFT_CERN_Complete_164.root\")");
	gROOT->ProcessLine("macro2->run(2,666666,-1,0,1,1,36,1,1,0)");
	gROOT->ProcessLine(".q");


	//}
}
