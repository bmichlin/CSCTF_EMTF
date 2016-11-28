void asdf(){
	gROOT->ProcessLine(".x initL1Analysis.C+");
	gROOT->ProcessLine(".L New_output_WalkAround.C+");
	gROOT->ProcessLine("WalkAround *macro = new WalkAround(\"InputFile\")");
	gROOT->ProcessLine("macro->run(filenum,runnum,-1,0,station,ring,chambers,trig,trig,0)");
	gROOT->ProcessLine(".q");
}
