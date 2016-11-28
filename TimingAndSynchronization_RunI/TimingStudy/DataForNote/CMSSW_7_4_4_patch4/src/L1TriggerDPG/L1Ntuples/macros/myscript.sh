# Actually, let's do something more useful than a simple hello world... this will print the input arguments passed to the script
echo "Here there are all the input arguments"
echo $@

# If you are curious, you can have a look at the tweaked PSet. This however won't give you any information...
echo "================= PSet.py file =================="
#cat PSet.py

# This is what you need if you want to look at the tweaked parameter set!!
echo "================= Dumping PSet ===================="
#python -c "import PSet; print PSet.process.dumpPython()"

# Ok, let's stop fooling around and execute the job:
#cmsRun -j FrameworkJobReport.xml -p PSet.py

echo "======RUN ROOT======"
root -b -q temp_163.C
