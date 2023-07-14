#include "step1RDF_forLJMet.cpp"
#include "cleanJet.cc"
#include "dnnPrep.cc"
#include "W_t_reco.cc"
#include "TBPrime.cc"
#include "utilities.cc"

void runRDF(TString testNum)
{
	// TPrime root files
//	std::vector<std::string> inputFiles = {"root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv7/TprimeTprime_M-1200_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/110000/84EFA4DB-56CC-5C46-9B7C-2286B546B23B.root","root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv7/TprimeTprime_M-1200_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/110000/5B3E2711-DA8D-FF4E-8E2E-D1E1F57E4BEE.root"};

	// TTToSemiLeptonic root files
//	std::vector<std::string> inputFiles = {"root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/70000/DDE55425-6106-4B4E-958C-87545D734E2B.root","root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv7/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/70000/DBB9294E-FBDF-8A48-B666-CD4F7187C317.root"}; // Incomplete file list

	// Data root files
  std::vector<std::string> inputFiles = {"root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL18NanoAODv9/TprimeTprime_M-1500_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/447AD74F-034B-FA42-AD05-CD476A98C43D.root","root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL18NanoAODv9/TprimeTprime_M-1500_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/8FC43306-01B6-A343-848C-1309279D568E.root","root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL18NanoAODv9/TprimeTprime_M-1500_TuneCP5_13TeV-madgraph-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/BD29BF55-6295-0044-8621-191CD4D5CF72.root"};

  rdf t(inputFiles[0],"preselTree_"+testNum,"finalselTree_"+testNum); // names get set to class members, should be known w/o passing
  int year = 2018;
  std::cout << "Number of Root Files: " << inputFiles.size() << std::endl;
  
  t.step1RDF_forLJMet(inputFiles, testNum ,year);
};
