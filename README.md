# vlq-1lepDnn-RDF
RDataFrame VLQ pair -> single lepton analysis

Replication of step1 for LJMet-Slimmer-1lepDnn analysis with RDataFrame structure for increased readability, and decreased run-time + 'sit-and-wait' time.
callRDF.C is used as the command line script that will then call runRDF.C, which houses the main step1RDF_forLJMet.cpp and step1RDF_forLJMet.h files, as well as various .cc files containing functions used within the .cpp. They were moved into .cc scripts to maintain some of the integrity of .cpp and not suffocate the analysis code in functions. The input root files then go through a series of definitions and cuts, with an option to save a snapshot at both preselection and post step1 analysis. This results in a .root file containing the defined branches for plotting and/or furthur analysis.

TO RUN:
If you intend to save a snapshot, make sure that the final .Snapshot at the very end of the script is uncommented. There's an optional place to start the analysis from a preselection .root file, which cuts down on time significatly. Thus, if you already have preselection root files, I reccomend utilizing that area. If not, the analysis takes anywhere from 30 minutes to 3 hours depending on the size of your input. To run on the command line, use the following:

root -l callRDF.C\(\"Muon(OR)Electron\",\"testNumber\"\)

Additionally, we're in the process of adding additional files for condor job submission.
Happy analyzing!
