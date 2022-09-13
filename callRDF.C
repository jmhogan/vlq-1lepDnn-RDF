// callRDF.C with lwtnn includes
void callRDF(TString channel, TString testNum)
{
	gSystem->Load("~/nobackup/YOURWORKINGAREA/CMSSW_11_0_0/lib/slc7_amd64_gcc820/liblwtnnlwtnn.so");
	gROOT->ProcessLine(".x runRDF.C(\""+channel+"\",\""+testNum+"\")");
};
