# Lxplus Batch Job Script

cd /afs/cern.ch/work/b/bmichlin/private/DataForNote/CMSSW_7_4_4_patch4/src/
cmsenv
cd L1TriggerDPG/L1Ntuples/macros/
root -b -q temp_1486.C >& run_bash.log
