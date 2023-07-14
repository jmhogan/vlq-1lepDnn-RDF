// --------------------------------------------------------------------------------------- //
// Implimentation of RDataFrame in C++.					                   //
// Replication of step1.cc for VLQ analysis				                   //
// To Run on Command Line:   root -l callRDF.C\(\"Muon(OR)Electron\",\"testNumber\"\)      //
// --------------------------------------------------------------------------------------- //

#define rdf_cxx
#include "step1RDF_forLJMet.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include <TFile.h>
#include <TChain.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH3.h>
#include <algorithm> // std::sort
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <sstream>
#include <chrono> // for high_resolution_clock

using namespace ROOT::VecOps;
void rdf::step1RDF_forLJMet(std::vector<std::string> sample, TString testNum, int year)
{
  ROOT::EnableImplicitMT();
  //	const std::string samplesBasePath = "root://cmsxrootd.fnal.gov/";
  TStopwatch time;
  time.Start();
  bool isNominal = isNominal;
  TString region = "Signal";
  if(isSig == true){region = "Signal";} // TPrimeTPrime or BPrimeBPrime
  else if(isTT == true){region = "TTToSemiLeptonic";} // TTToSemiLeptonic
  std::cout<< "Region: " << region << std::endl;
  
  // --------------------------------------------------------------------------------------------------------------------
  // 							LAMBDA FXNS
  // --------------------------------------------------------------------------------------------------------------------
  
  // ----------------------------------------------------
  //   		DEACY/GENTTBAR CALCULATOR:
  // ----------------------------------------------------
  auto decayModeSelection_genTTbarMassCalc = [region](unsigned int nGenPart, ROOT::VecOps::RVec<int>& GenPart_pdgId, ROOT::VecOps::RVec<float>& GenPart_mass, ROOT::VecOps::RVec<float>& GenPart_pt, ROOT::VecOps::RVec<float>& GenPart_phi, ROOT::VecOps::RVec<float>& GenPart_eta, ROOT::VecOps::RVec<int>& GenPart_genPartIdxMother, ROOT::VecOps::RVec<int>& GenPart_status)
    {
      int returnVar = 0;
      if(region == "Signal")
	{
	  std::vector<int> tPrimeID;
	  std::vector<int> bPrimeID;
	  std::vector<int> listofQuarkIDs;
	  std::vector<int> listofBosonIDs;
	  std::vector<unsigned int> quarks;
	  std::vector<unsigned int> bosons;
	  
	  bool isBWBW = false;
	  bool isTZTZ = false;
	  bool isTHTH = false;
	  bool isTZTH = false;
	  bool isTZBW = false;
	  bool isTHBW = false;
	  
	  bool isTWTW = false;
	  bool isBZBZ = false;
	  bool isBHBH = false;
	  bool isBZBH = false;
	  bool isBZTW = false;
	  bool isBHTW = false;
	  
	  int decayMode = 0;
	  
	  tPrimeID.clear();
	  bPrimeID.clear();
	  listofQuarkIDs.clear();
	  listofBosonIDs.clear();
	  quarks.clear();
	  bosons.clear();
	  
	  for(unsigned int p = 0; p < nGenPart; p++)
	    {
	      int id=GenPart_pdgId[p];
	      // find T' and B' particles
	      if(abs(id) != 8000001 && abs(id) != 8000002){continue;}
	      bool hasTdaughter = false;
	      vector<unsigned int> daughters;
	      daughters.clear();
	      for(unsigned int  dau = 0; dau < nGenPart; dau++)
		{
		  if(GenPart_genPartIdxMother[dau]!=p){continue;}
		  daughters.push_back(dau);
		  if(abs(id) == 8000001 && abs(GenPart_pdgId[dau]) == 8000001){hasTdaughter = true;}
		  if(abs(id) == 8000002 && abs(GenPart_pdgId[dau]) == 8000002){hasTdaughter = true;}
		}
	      if(hasTdaughter){continue;}
	      int mother = GenPart_genPartIdxMother[p];
	      int mother_id = GenPart_pdgId[mother];
	      if(abs(id) == 8000001)
		{
		  if(abs(mother_id) == 8000001){tPrimeID.push_back(GenPart_pdgId[mother]);}
		  else{tPrimeID.push_back(GenPart_pdgId[p]);}
		}
	      if(abs(id) == 8000002)
		{
		  if(abs(mother_id) == 8000002){bPrimeID.push_back(GenPart_pdgId[mother]);}
		  else{bPrimeID.push_back(GenPart_pdgId[p]);}
		}
	      for(unsigned int j = 0; j < daughters.size(); j++)
		{
		  unsigned int d = daughters.at(j);
		  int dauId = GenPart_pdgId[d];
		  if(abs(dauId) == 5 || abs(dauId) == 6)
		    {
		      quarks.push_back(d);
		      listofQuarkIDs.push_back(dauId);
		    }
		  else if(abs(dauId) > 22 && abs(dauId) < 26)
		    {
		      bosons.push_back(d);
		      listofBosonIDs.push_back(dauId);
		    }
		  else{continue;}
		}
	    }
	  
	  if(tPrimeID.size() > 0 && bPrimeID.size() > 0) {std::cout << "Found both T' and B' " << std::endl;}
	  if(listofQuarkIDs.size() != 0 && listofQuarkIDs.size() != 2)
	    {
	      std::cout << "More/less than 2 quarks stored: " << listofQuarkIDs.size() << std::endl;
	      for(unsigned int i = 0; i < listofQuarkIDs.size(); i++){std::cout << "quark " << i << " = " << listofQuarkIDs.at(i) << std::endl;}
	      int test = listofQuarkIDs.at(0)*listofQuarkIDs.at(1);
	      int sign = -1;
	      if(test > 0){sign = 1;}
	      if(sign > 0)
		{
		  if(listofQuarkIDs.size() == 4)
		    {
		      std::swap(listofQuarkIDs.at(2),listofQuarkIDs.at(3));
		      std::swap(quarks.at(2),quarks.at(3));
		    }
		  std::swap(listofQuarkIDs.at(1),listofQuarkIDs.at(2));
		  std::swap(quarks.at(1),quarks.at(2));
		  test = listofQuarkIDs.at(0)*listofQuarkIDs.at(1);
		  sign = -1;
		  if(test > 0){sign = 1;}
		  if(sign < 0){std::cout << "Signs are fixed!" << std::endl;}
		}
	      if(listofQuarkIDs.size() > 3 && abs(listofQuarkIDs.at(3)) == 6)
		{
		  std::swap(listofQuarkIDs.at(2),listofQuarkIDs.at(3));
		  std::swap(quarks.at(2),quarks.at(3));
		}
	      if(listofQuarkIDs.size() > 2 && abs(listofQuarkIDs.at(2)) == 6)
		{
		  std::swap(listofQuarkIDs.at(1),listofQuarkIDs.at(2));
		  std::swap(quarks.at(1),quarks.at(2));
		}
	    }
	  if(listofBosonIDs.size() != 0 && listofBosonIDs.size() != 2)
	    {
	      std::cout << "More/less than 2 bosons stored: " << listofBosonIDs.size() << std::endl;
	    }
	  // tag the decay chains according to ID'd quarks and bosons.
	  
	  // TPrime Decay Mode Selector
	  if(tPrimeID.size() > 1 && bPrimeID.size() == 0)
	    {
	      if(abs(listofQuarkIDs.at(0)) == 5 && abs(listofQuarkIDs.at(1)) == 5)
		{
		  if(abs(listofBosonIDs.at(0)) == 24 && abs(listofBosonIDs.at(1)) == 24)
		    {
		      isBWBW = true;
		      decayMode = 1; // BWBW ID!
		    }
		}
	      // 2 t quarks, check for Z's and H's
	      else if(abs(listofQuarkIDs.at(0)) == 6 && abs(listofQuarkIDs.at(1)) == 6)
		{
		  if(listofBosonIDs.at(0) == 23 && listofBosonIDs.at(1) == 23)
		    {
		      isTZTZ = true;
		      decayMode = 2; // TZTZ ID!
		    }
		  else if(listofBosonIDs.at(0) == 25 && listofBosonIDs.at(1) == 25)
		    {
		      isTHTH = true;
		      decayMode = 3; // THTH ID!
		    }
		  else if(listofBosonIDs.at(0) == 25 && listofBosonIDs.at(1) == 23)
		    {
		      isTZTH = true;
		      decayMode = 4; //TZTH ID!
		    }
		  else if(listofBosonIDs.at(0) == 23 && listofBosonIDs.at(1) == 25)
		    {
		      isTZTH = true;
		      decayMode = 4; // TZTH ID!
		    }
		  else
		    {
		      std::cout << "2 t daughters didn't match tZtZ, tHtH, or tZtH" << listofBosonIDs.at(0) << ", " << listofBosonIDs.at(1) << std::endl;
		    }
		}
	      // t-b pairs, check for correlating bosons in the right spots
	      else if(abs(listofQuarkIDs.at(0)) == 6 && abs(listofQuarkIDs.at(1)) == 5)
		{
		  if(listofBosonIDs.at(0) == 23 && abs(listofBosonIDs.at(1)) == 24)
		    {
		      isTZBW = true;
		      decayMode = 5; // TZBW ID!
		    }
		  else if(listofBosonIDs.at(0) == 25 && abs(listofBosonIDs.at(1)) == 24)
		    {
		      isTHBW = true;
		      decayMode = 6; // THBW ID!
		    }
		  else{std::cout<< "t - b pair didn't match Z/H - W pair" << listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
		}
	      // b-t pairs, check for correlating bosons in the right spots
	      else if(abs(listofQuarkIDs.at(1)) == 6 && abs(listofQuarkIDs.at(0)) == 5)
		{
		  if(listofBosonIDs.at(1) == 23 && abs(listofBosonIDs.at(0)) == 24)
		    {
		      isTZBW = true;
		      decayMode = 5; // TZBW ID!
		    }
		  else if(listofBosonIDs.at(1) == 25 && abs(listofBosonIDs.at(0)) == 24)
		    {
		      isTHBW = true;
		      decayMode = 6; //THBW ID!
		    }
		  else{std::cout<< "b - t pair didn't match W - Z/H pair" << listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
		}
	      // error messages if we found something else entirely
	      else
		{
		  std::cout << "T' daughters didn't match a recognized pattern" << std::endl;
		  for(size_t i = 0; i < listofQuarkIDs.size(); i++)
		    {
		      std::cout << "quark " << i << " = " << listofQuarkIDs.at(i) << std::endl;
		    }
		  for(size_t i = 0; i < listofBosonIDs.size(); i++)
		    {
		      std::cout << "boson " << i << " = " << listofBosonIDs.at(i) << std::endl;
		    }
		  decayMode = -1;
		}
	    }
	  // BPrime Decay Mode Selector
	  if(bPrimeID.size() > 1 && tPrimeID.size() == 0)
	    {
	      // 2 t quarks, check for matching W's
	      if(abs(listofQuarkIDs.at(0)) == 6 && abs(listofQuarkIDs.at(1)) == 6)
		{
		  if(abs(listofBosonIDs.at(0)) == 24 && abs(listofBosonIDs.at(1)) == 24)
		    {
		      isTWTW = true;
		      decayMode = 1; // TWTW ID!
		    }
		  else{std::cout<< "2 t daughters didn't match tWtW: " <<listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
		}
	      // 2 b quarks, check for Z's and H's
	      else if(abs(listofQuarkIDs.at(0)) == 5 && abs(listofQuarkIDs.at(1)) == 5)
		{
		  if(listofBosonIDs.at(0) == 23 && listofBosonIDs.at(1) == 23)
		    {
		      isBZBZ = true;
		      decayMode = 2; // BZBZ ID!
		    }
		  else if(listofBosonIDs.at(0) == 25 && listofBosonIDs.at(1) == 25)
		    {
		      isBHBH = true;
		      decayMode = 3; // BHBH ID!
		    }
		  else if(listofBosonIDs.at(0) == 25 && listofBosonIDs.at(1) == 23)
		    {
		      isBZBH = true;
		      decayMode = 4; // BZBH ID!
		    }
		  else if(listofBosonIDs.at(0) == 23 && listofBosonIDs.at(1) == 25)
		    {
		      isBZBH = true;
		      decayMode = 4; //BZBH ID!
		    }
		  else
		    {
		      std::cout << "2 b daughters didn't match bZbZ, bHbH, or bZbH" << listofBosonIDs.at(0) << ", " << listofBosonIDs.at(1) << std::endl;
		    }
		}
	      // b-t pairs, check for correlating bosons in the right spots
	      else if(abs(listofQuarkIDs.at(0)) == 5 && abs(listofQuarkIDs.at(1)) == 6)
		{
		  if(listofBosonIDs.at(0) == 23 && abs(listofBosonIDs.at(1)) == 24)
		    {
		      isBZTW = true;
		      decayMode = 5; // BZTW ID!
		    }
		  else if(listofBosonIDs.at(0) == 25 && abs(listofBosonIDs.at(1)) == 24)
		    {
		      isBHTW = true;
		      decayMode = 6; // BHTW ID!
		    }
		  else{std::cout<< "b - t pair didn't match Z/H - W pair" << listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
		}
	      // t-b pairs, check for correlating bosons in the right spots
	      else if(abs(listofQuarkIDs.at(1)) == 5 && abs(listofQuarkIDs.at(0)) == 6)
		{
		  if(listofBosonIDs.at(1) == 23 && abs(listofBosonIDs.at(0)) == 24)
		    {
		      isBZTW = true;
		      decayMode = 5; // BZTW ID!
		    }
		  else if(listofBosonIDs.at(1) == 25 && abs(listofBosonIDs.at(0)) == 24)
		    {
		      isBHTW = true;
		      decayMode = 6; // BHTW ID!
		    }
		  else{std::cout<< "t - b pair didn't match W - Z/H pair" << listofBosonIDs.at(0)<<", "<<listofBosonIDs.at(1) << std::endl;}
		}
	      // error messages if we found something else entirely
	      else
		{
		  std::cout << "B' daughters didn't match a recognized pattern" << std::endl;
		  for(size_t i = 0; i < listofQuarkIDs.size(); i++)
		    {
		      std::cout << "quark " << i << " = " << listofQuarkIDs.at(i) << std::endl;
		    }
		  for(size_t i = 0; i < listofBosonIDs.size(); i++)
		    {
		      std::cout << "boson " << i << " = " << listofBosonIDs.at(i) << std::endl;
		    }
		  decayMode = -1;
		}
	    }
	  returnVar = decayMode;
	}
      
      else if(region == "TTToSemiLeptonic")
	{
	  int genTTbarMass = -999;
	  double topPtWeight = 1.0;
	  TLorentzVector top, antitop;
	  bool gottop = false;
	  bool gotantitop = false;
	  bool gottoppt = false;
	  bool gotantitoppt = false;
	  float toppt, antitoppt;
	  for(unsigned int p = 0; p < nGenPart; p++)
	    {
	      int id = GenPart_pdgId[p];
	      if (abs(id) != 6){continue;}
	      if (GenPart_mass[p] < 10){continue;}
	      int motherid = GenPart_pdgId[GenPart_genPartIdxMother[p]];
	      if(abs(motherid) != 6)
		{
		  if (!gottop && id == 6)
		    {
		      top.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
		      gottop = true;
		    }
		  if (!gotantitop && id == -6)
		    {
		      antitop.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
		      gotantitop = true;
		    }
		}
	      if(GenPart_status[p] == 62)
		{
		  if (!gottoppt && id == 6)
		    {
		      toppt = GenPart_pt[p];
		      gottoppt = true;
		    }
		  if (!gotantitoppt && id == -6)
		    {
		      antitoppt = GenPart_pt[p];
		      gotantitoppt = true;
		    }
		}
	    }
	  if(gottop && gotantitop){genTTbarMass = (top+antitop).M();}
	  if(gottoppt && gotantitoppt)
	    {
	      float SFtop = TMath::Exp(0.0615-0.0005*toppt);
	      float SFantitop = TMath::Exp(0.0615-0.0005*antitoppt);
	      topPtWeight = TMath::Sqrt(SFtop*SFantitop);
	    }
	  returnVar = genTTbarMass;
	}
      return returnVar;
    };
  
  // ----------------------------------------------------
  //     minMleppJet VECTOR RETURN + NJETSDEEPFLAV
  // ----------------------------------------------------
  auto minMleppJet_calc = [isNominal](ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, TLorentzVector lepton_lv, ROOT::VecOps::RVec<float> Jet_btag)
    {
      float ind_MinMlj = -1; // This gets changed into int in .Define()
      float minMleppJet = 1e8;
      ROOT::VecOps::RVec<int> theJetBTag (jet_pt.size(),0);
      float NJetsDeepFlavwithSF = 0;
      TLorentzVector jet_lv;
      
      for(unsigned int ijet=0; ijet < jet_pt.size(); ijet++)
	{
	  jet_lv.SetPtEtaPhiM(jet_pt.at(ijet),jet_eta.at(ijet),jet_phi.at(ijet),jet_mass.at(ijet));
	  if(Jet_btag[ijet] > 0.2783){theJetBTag.at(ijet) = 1;} // BTagged or not
	  else if(Jet_btag[ijet] < 0.2783){theJetBTag.at(ijet) = 0;}
	  
	  if((lepton_lv + jet_lv).M() < minMleppJet)
	    {
	      minMleppJet = fabs((lepton_lv + jet_lv).M());
	      ind_MinMlj = ijet;
	    }
	  if(isNominal && theJetBTag.at(ijet) == 1){NJetsDeepFlavwithSF += 1;}
	}
      ROOT::VecOps::RVec<float> minMlj = {minMleppJet,ind_MinMlj,NJetsDeepFlavwithSF};
      return minMlj;
    };
  
  
  // -------------------------------------------------------
  //               Flags and First Filter 
  // -------------------------------------------------------
  auto rdf = ROOT::RDataFrame("Events" , sample); // Initial data
  auto METfilters = rdf.Filter("Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_goodVertices == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_eeBadScFilter == 1 && Flag_globalSuperTightHalo2016Filter == 1 && Flag_BadPFMuonFilter == 1 && Flag_ecalBadCalibFilter == 1", "MET Filters")
    .Filter("MET_pt > 50", "Pass MET > 50")
    .Filter("nJet > 0 && nFatJet > 0", "Event has jets");
    
  auto LepDefs = METfilters.Define("TightMu", "abs(Muon_eta) < 2.4 && Muon_tightId == true && Muon_miniIsoId >= 3 && Muon_pt > 30")
    .Define("TightEl", "abs(Electron_eta) < 2.5 && Electron_mvaFall17V2noIso_WP90 == true && Electron_miniPFRelIso_all < 0.1 && Electron_pt > 30")
    .Define("VetoMu", "abs(Muon_eta) < 2.4 && Muon_looseId == true && Muon_miniIsoId >= 1 && Muon_pt > 10 && TightMu == false")
    .Define("VetoEl", "abs(Electron_eta) < 2.5 && Electron_mvaFall17V2noIso_WPL == true && Electron_miniPFRelIso_all < 0.4 && Electron_pt > 10 && TightEl == false")
    .Define("nVetoLep", "(int) (Sum(VetoMu)+Sum(VetoEl))")
    .Define("nTightMu", "(int) (Sum(TightMu))")
    .Define("nTightEl", "(int) (Sum(TightEl))")
    .Define("TMuon_pt", "Muon_pt[TightMu == true]")
    .Define("TMuon_eta", "Muon_eta[TightMu == true]")
    .Define("TMuon_phi", "Muon_phi[TightMu == true]")
    .Define("TMuon_mass", "Muon_mass[TightMu == true]")
    .Define("TElectron_pt", "Electron_pt[TightEl == true]")
    .Define("TElectron_eta", "Electron_eta[TightEl == true]")
    .Define("TElectron_phi", "Electron_phi[TightEl == true]")
    .Define("TElectron_mass", "Electron_mass[TightEl == true]")
    .Define("TMuon_P4", "fVectorConstructor(TMuon_pt,TMuon_eta,TMuon_phi,TMuon_mass)")
    .Define("TElectron_P4", "fVectorConstructor(TElectron_pt,TElectron_eta,TElectron_phi,TElectron_mass)")
    .Define("TMuon_jetIdx", "Muon_jetIdx[TightMu == true]")
    .Define("TElectron_jetIdx", "Electron_jetIdx[TightEl == true]");

  auto LepSelect = LepDefs.Define("isMu","nMuon > 0 && nTightMu == 1 && (nElectron == 0 || nTightEl == 0) && nVetoLep == 0 && (HLT_Mu50 == 1 || HLT_Mu15_IsoVVVL_PFHT450 == 1)") \
    .Define("isEl","nElectron > 0 && nTightEl == 1 && (nMuon == 0 || nTightMu == 0) && nVetoLep==0 && (HLT_Ele35_WPTight_Gsf == 1 || HLT_Ele15_IsoVVVL_PFHT450 == 1)") \
    .Filter("isMu || isEl","Event is either muon or electron");

  auto LepAssign = LepSelect.Define("assignleps","assign_leps(isMu,isEl,TightMu,TightEl,Muon_pt,Muon_eta,Muon_phi,Muon_mass,Muon_miniPFRelIso_all,Electron_pt,Electron_eta,Electron_phi,Electron_mass,Electron_miniPFRelIso_all)") \
    .Define("lepton_pt","assignleps[0]")				\
    .Define("lepton_eta","assignleps[1]")				\
    .Define("lepton_phi","assignleps[2]")				\
    .Define("lepton_mass","assignleps[3]")				\
    .Define("lepton_miniIso","assignleps[4]");
  
  auto JetCleaner = LepAssign.Define("Jet_P4", "fVectorConstructor(Jet_pt,Jet_eta,Jet_phi,Jet_mass)")
    .Define("cleanedJets", "cleanJets(Jet_P4,Jet_rawFactor,TMuon_P4,TMuon_jetIdx,TElectron_P4,TElectron_jetIdx)")
    .Define("cleanedJet_pt", "cleanedJets[0]")
    .Define("cleanedJet_eta", "cleanedJets[1]")
    .Define("cleanedJet_phi", "cleanedJets[2]")
    .Define("cleanedJet_mass", "cleanedJets[3]")
    .Define("cleanedJet_rawFactor", "cleanedJets[4]")
    .Define("DR_lepJets","DeltaR_VecAndFloat(cleanedJet_eta,cleanedJet_phi,lepton_eta,lepton_phi)")
    .Define("ptrel_lepJets","ptRel(cleanedJet_pt,cleanedJet_eta,cleanedJet_phi,cleanedJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)") \
    .Define("goodcleanJets", "cleanedJet_pt > 30 && abs(cleanedJet_eta) < 2.4 && Jet_jetId > 1 && (DR_lepJets > 0.4 || ptrel_lepJets > 20)")
    .Define("NJets_central", "(int) Sum(goodcleanJets)")
    .Define("gcJet_pt", "cleanedJet_pt[goodcleanJets == true]")
    .Define("gcJet_eta", "cleanedJet_eta[goodcleanJets == true]")
    .Define("gcJet_phi", "cleanedJet_phi[goodcleanJets == true]")
    .Define("gcJet_mass", "cleanedJet_mass[goodcleanJets == true]")
    .Define("gcJet_DeepFlav", "Jet_btagDeepFlavB[goodcleanJets == true]")
    .Define("gcJet_DeepFlavM", "gcJet_DeepFlav > 0.2783")
    .Define("NJets_DeepFlavM", "(int) Sum(gcJet_DeepFlavM)")
    .Define("DR_lepFatJets","DeltaR_VecAndFloat(FatJet_eta,FatJet_phi,lepton_eta,lepton_phi)")
    .Define("ptrel_lepFatJets","ptRel(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_mass,lepton_pt,lepton_eta,lepton_phi,lepton_mass)") \
    .Define("goodcleanFatJets", "FatJet_pt > 200 && abs(FatJet_eta) < 2.4 && FatJet_jetId > 1 && (DR_lepFatJets > 0.8 || ptrel_lepFatJets > 20)")
    .Define("NFatJets_central", "(int) Sum(goodcleanFatJets)")
    .Define("gcFatJet_pt", "FatJet_pt[goodcleanFatJets == true]")
    .Define("gcFatJet_eta", "FatJet_eta[goodcleanFatJets == true]")
    .Define("gcFatJet_phi", "FatJet_phi[goodcleanFatJets == true]")
    .Define("gcFatJet_mass", "FatJet_mass[goodcleanFatJets == true]")
    .Define("gcFatJet_sdmass", "FatJet_msoftdrop[goodcleanFatJets == true]");
  

  // ---------------------------------------------------------
  // 	  HT Calculation and Final Preselection Cut
  // ---------------------------------------------------------
  auto HT_calc = JetCleaner.Define("AK4HT","Sum(gcJet_pt)")
    .Filter("AK4HT > 510","AK4 HT Pass")
    .Filter("NFatJets_central > 2", "3 AK8s Pass");

  // ---------------------------------------------------------
  // 		Post Preselection Analysis
  // ---------------------------------------------------------
  auto postPresel = HT_calc.Define("decayMode_or_genTTbarMass",decayModeSelection_genTTbarMassCalc,{"nGenPart","GenPart_pdgId","GenPart_mass", \
  	"GenPart_pt","GenPart_phi","GenPart_eta",			\
  	"GenPart_genPartIdxMother","GenPart_status"})			\
    .Define("lepton_lv","lvConstructor(lepton_pt,lepton_eta,lepton_phi,lepton_mass)")\
    .Define("AK4HTpMETpLepPt","AK4HT + lepton_pt + MET_pt") \
    .Define("dnnJ","FatJet_deepTag_QCDothers[goodcleanFatJets == true]")	\
    .Define("int_dnnT","(FatJet_deepTag_TvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_TvsQCD)") \
    .Define("dnnT","int_dnnT[goodcleanFatJets == true]")			\
    .Define("dnnH","FatJet_deepTag_H[goodcleanFatJets == true]")		\
    .Define("int_dnnZ","(FatJet_deepTag_ZvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_ZvsQCD)") \
    .Define("dnnZ","int_dnnZ[goodcleanFatJets == true]")			\
    .Define("int_dnnW","(FatJet_deepTag_WvsQCD * FatJet_deepTag_QCD) / (1 - FatJet_deepTag_WvsQCD)") \
    .Define("dnnW","int_dnnW[goodcleanFatJets == true]")			\
    .Define("int_dnnB","(FatJet_deepTag_QCD - FatJet_deepTag_QCDothers)") \
    .Define("dnnB","int_dnnB[goodcleanFatJets == true]")			\
    .Define("dnnLargest","maxFxn(dnnJ,dnnT,dnnH,dnnZ,dnnW,dnnB)")	\
    .Define("nJ_DeepAK8","Sum(dnnLargest == 0)")			\
    .Define("nT_DeepAK8","Sum(dnnLargest == 1)")			\
    .Define("nH_DeepAK8","Sum(dnnLargest == 2)")			\
    .Define("nZ_DeepAK8","Sum(dnnLargest == 3)")			\
    .Define("nW_DeepAK8","Sum(dnnLargest == 4)")			\
    .Define("nB_DeepAK8","Sum(dnnLargest == 5)")			\
    .Define("int_tau21","(FatJet_tau2 / FatJet_tau1)")			\
    .Define("tau21","int_tau21[goodcleanFatJets == true]")			\
    .Define("tau21_1","tau21[0]")					\
    .Define("tau21_2","tau21[1]")					\
    .Define("tau21_3","tau21[2]")					\
    .Define("minDR_ptRel_lead_lepAK8","minDR_ptRel_lead_calc(gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,lepton_lv)")\
    .Define("minDR_lepAK8","minDR_ptRel_lead_lepAK8[0]")		\
    .Define("ptRel_lepAK8","minDR_ptRel_lead_lepAK8[1]")		\
    .Define("minDR_leadAK8otherAK8","minDR_ptRel_lead_lepAK8[2]")	\
    .Define("minDR_ptRel_lead_lepJets","minDR_ptRel_lead_calc(gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass,lepton_lv)")\
    .Define("minDR_lepJet","minDR_ptRel_lead_lepJets[0]")		\
    .Define("ptRel_lepJet","minDR_ptRel_lead_lepJets[1]")		\
    .Define("DR_lepAK8s","DeltaR_VecAndFloat(gcFatJet_eta,gcFatJet_phi,lepton_eta,lepton_phi)")\
    .Define("W_lv","lpNu_WCalc(MET_pt,MET_phi,lepton_lv)") \
    .Define("minMlj",minMleppJet_calc,{"gcJet_pt","gcJet_eta","gcJet_phi","gcJet_mass", "lepton_lv","gcJet_DeepFlav"})				\
    .Define("W_dRLep","dR_Wt_Calc(W_lv,lepton_lv)")			\
    .Define("minMleppJet","minMlj[0]")					\
    .Define("ind_MinMlj","(int) minMlj[1]")				\
    .Define("NJetsDeepFlavwithSF","(int) minMlj[2]")		\
    .Define("isLeptonic","isLeptonic_X(minMleppJet)")			\
    .Define("t_lv","lpNu_t_Calc(isLeptonic,gcJet_pt,gcJet_eta,gcJet_phi,gcJet_mass,W_lv,minMleppJet,ind_MinMlj)")\
    .Define("t_pt","t_lv[0]")						\
    .Define("t_eta","t_lv[1]")						\
    .Define("t_phi","t_lv[2]")						\
    .Define("t_mass","t_lv[3]")						\
    .Define("t_dRWb","t_lv[4]")						\
    .Define("top_lv","top_lvConstructor(t_pt,t_eta,t_phi,t_mass)")	\
    .Define("tj_vec","three_jet(top_lv,W_lv,isLeptonic,gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,dnnT,dnnH,dnnZ,dnnW,dnnB,dnnLargest,gcFatJet_sdmass)")\
    .Define("Tprime1_DeepAK8_Mass","tj_vec[0]")				\
    .Define("Tprime2_DeepAK8_Mass","tj_vec[1]")				\
    .Define("Tprime1_DeepAK8_pt","tj_vec[2]")				\
    .Define("Tprime2_DeepAK8_pt","tj_vec[3]")				\
    .Define("Tprime1_DeepAK8_eta","tj_vec[4]")				\
    .Define("Tprime2_DeepAK8_eta","tj_vec[5]")				\
    .Define("Tprime1_DeepAK8_Phi","tj_vec[6]")				\
    .Define("Tprime2_DeepAK8_Phi","tj_vec[7]")				\
    .Define("Tprime1_DeepAK8_deltaR","tj_vec[8]")			\
    .Define("Tprime2_DeepAK8_deltaR","tj_vec[9]")			\
    .Define("Bprime1_DeepAK8_Mass","tj_vec[10]")			\
    .Define("Bprime2_DeepAK8_Mass","tj_vec[11]")			\
    .Define("Bprime1_DeepAK8_pt","tj_vec[12]")				\
    .Define("Bprime2_DeepAK8_pt","tj_vec[13]")				\
    .Define("Bprime1_DeepAK8_eta","tj_vec[14]")				\
    .Define("Bprime2_DeepAK8_eta","tj_vec[15]")				\
    .Define("Bprime1_DeepAK8_Phi","tj_vec[16]")				\
    .Define("Bprime2_DeepAK8_Phi","tj_vec[17]")				\
    .Define("Bprime1_DeepAK8_deltaR","tj_vec[18]")			\
    .Define("Bprime2_DeepAK8_deltaR","tj_vec[19]")			\
    .Define("leptonicTprimeJetIDs_DeepAK8","(int) tj_vec[23]")		\
    .Define("leptonicBprimeJetIDs_DeepAK8","(int) tj_vec[24]")		\
    .Define("hadronicTprimeJetIDs1_DeepAK8","(int) tj_vec[25]")		\
    .Define("hadronicTprimeJetIDs2_DeepAK8","(int) tj_vec[26]")		\
    .Define("hadronicBprimeJetIDs1_DeepAK8","(int) tj_vec[27]")		\
    .Define("hadronicBprimeJetIDs2_DeepAK8","(int) tj_vec[28]")		\
    .Define("TPrime1_lv","top_lvConstructor(Tprime1_DeepAK8_pt,Tprime1_DeepAK8_eta,Tprime1_DeepAK8_Phi,Tprime1_DeepAK8_Mass)") \
    .Define("TPrime2_lv","top_lvConstructor(Tprime2_DeepAK8_pt,Tprime2_DeepAK8_eta,Tprime2_DeepAK8_Phi,Tprime2_DeepAK8_Mass)") \
    .Define("BPrime1_lv","top_lvConstructor(Bprime1_DeepAK8_pt,Bprime1_DeepAK8_eta,Bprime1_DeepAK8_Phi,Bprime1_DeepAK8_Mass)") \
    .Define("BPrime2_lv","top_lvConstructor(Bprime2_DeepAK8_pt,Bprime2_DeepAK8_eta,Bprime2_DeepAK8_Phi,Bprime2_DeepAK8_Mass)") \
    .Define("TTagged","three_jet_TTag(top_lv,W_lv,isLeptonic,gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,dnnLargest, \
  				    				   hadronicTprimeJetIDs1_DeepAK8,hadronicTprimeJetIDs2_DeepAK8)")\
    .Define("validT","TTagged[0]")					\
    .Define("TtaggedDecay","TTagged[1]")				\
    .Define("BTagged","three_jet_BTag(top_lv,W_lv,isLeptonic,gcFatJet_pt,gcFatJet_eta,gcFatJet_phi,gcFatJet_mass,dnnLargest, \
  				    				   hadronicBprimeJetIDs1_DeepAK8,hadronicBprimeJetIDs2_DeepAK8)")\
    .Define("validB","BTagged[0]")					\
    .Define("BtaggedDecay","BTagged[1]");

  // -------------------------------------------------
  // 		Save Snapshot to file
  // -------------------------------------------------
  auto dfFinal = postPresel;
  std::cout << "-------------------------------------------------" << std::endl << ">>> Saving " << region << " Snapshot..." << std::endl;
  TString outputFile = "RDF_"+region+"_test"+testNum+".root";
  const char* stdOutputFile = outputFile;
  std::cout << "Output File: " << outputFile << std::endl << "-------------------------------------------------" << std::endl;

  auto colNames = dfFinal.GetColumnNames();
  vector<std::string> snapCol;
  int i = 0;
  for (auto &&ColName : colNames)
  {
    TString colName = ColName;
    if ((!colName.Contains("P4") && colName != "cleanedJets" && !colName.BeginsWith("L1") && !colName.BeginsWith("Gen") && !colName.BeginsWith("Soft") && !colName.BeginsWith("fixed") && !colName.BeginsWith("Sub") && !colName.BeginsWith("LHE") && !colName.BeginsWith("Raw") && !colName.BeginsWith("Calo") && !colName.BeginsWith("Chs") && !colName.BeginsWith("Corr") && !colName.BeginsWith("Fsr") && !colName.BeginsWith("Iso") && !colName.BeginsWith("Tau") && !colName.BeginsWith("SV") && !colName.BeginsWith("Puppi") && !colName.BeginsWith("Jet_") && !colName.BeginsWith("FatJet_") && !colName.BeginsWith("Photon") && !colName.BeginsWith("Low") && !colName.BeginsWith("HLT") && !colName.BeginsWith("Muon") && !colName.BeginsWith("Electron") && !colName.BeginsWith("boosted") && !colName.BeginsWith("Flag")))
    {
      std::string name = colName.Data();
      snapCol.push_back(name);
      i++;
    }
  }
  cout << "Number of Columns in Snapshot: " << i << endl;

  dfFinal.Snapshot("Events", stdOutputFile, snapCol);

  time.Stop();
  time.Print();

  cout << "Cut statistics:" << endl;
  dfFinal.Report()->Print();

  cout << "Adding Counter tree to the file:" << endl;
  auto rdf_runs = ROOT::RDataFrame("Runs", sample); 
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "UPDATE";
  rdf_runs.Snapshot("Runs", stdOutputFile, rdf_runs.GetColumnNames(), opts);

  cout << "Done!" << endl;

}
