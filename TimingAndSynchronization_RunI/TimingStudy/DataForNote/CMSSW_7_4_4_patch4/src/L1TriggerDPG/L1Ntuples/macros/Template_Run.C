void asdfasdf(){
	gROOT->ProcessLine(".x initL1Analysis.C+");
	gROOT->ProcessLine(".L New_WalkAround.C+");
	gROOT->ProcessLine("WalkAround *macro_RUNNUMBER_11 = new WalkAround(\"FILE\")");
	gROOT->ProcessLine("macro_RUNNUMBER_11->run(RUNNUMBER, -1, 0, 1, 1, 36, 1, 1, 0, 0)");
        gROOT->ProcessLine("WalkAround *macro_RUNNUMBER_12 = new WalkAround(\"FILE\")");
        gROOT->ProcessLine("macro_RUNNUMBER_12->run(RUNNUMBER, -1, 0, 1, 2, 36, 1, 1, 0, 0)");
        gROOT->ProcessLine("WalkAround *macro_RUNNUMBER_42 = new WalkAround(\"FILE\")");
        gROOT->ProcessLine("macro_RUNNUMBER_42->run(RUNNUMBER, -1, 0, 4, 2, 36, 1, 1, 0, 0)");
	gROOT->ProcessLine(".q");
}
