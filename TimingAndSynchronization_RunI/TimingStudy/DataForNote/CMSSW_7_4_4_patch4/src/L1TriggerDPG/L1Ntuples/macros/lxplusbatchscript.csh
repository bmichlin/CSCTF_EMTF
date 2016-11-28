# Lxplus Batch Job Script

cd /afs/cern.ch/work/b/bmichlin/private/DataForNote/CMSSW_7_4_4_patch4/src/
cmsenv
cd L1TriggerDPG/L1Ntuples/macros
root -l
.x initL1Analysis.C+
.L New_output_WalkAround.C+
WalkAround *macro = new WalkAround("CRUZET_Complete_Cosmics_1486.root")
macro->run(0,111111,-1,0,1,1,36,0,0,0)
.q
cd .
