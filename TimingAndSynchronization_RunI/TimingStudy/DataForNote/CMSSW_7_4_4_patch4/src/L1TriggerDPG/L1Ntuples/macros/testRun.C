void testRun(){
	gROOT->ProcessLine(".x initL1Analysis.C+");
	gROOT->ProcessLine(".L New_WalkAround.C+");
	gROOT->ProcessLine("WalkAround *macro = new WalkAround(\"root://eoscms///store/user/bmichlin/Cosmics/Individual_CRUZET_237205/150902_125618/0000/Combined_237205.root\")");
	gROOT->ProcessLine("macro->run(237205, -1, 0, 1, 1, 36, 1, 1, 0, 0)");
	gROOT->ProcessLine(".q");
}
