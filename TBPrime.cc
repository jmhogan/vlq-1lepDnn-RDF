// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// Fxn to return any and all float TPrime and BPrime variables needed for plotting
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
RVec<float> three_jet(TLorentzVector top_lv, TLorentzVector Wlv, int isLeptonic, RVec<float>& ak8_pt, RVec<float>& ak8_eta, RVec<float>& ak8_phi, RVec<float>& ak8_mass, RVec<float>& dnn_T_DeepAK8Calc_PtOrdered, RVec<float>& dnn_H_DeepAK8Calc_PtOrdered, RVec<float>& dnn_Z_DeepAK8Calc_PtOrdered, RVec<float>& dnn_W_DeepAK8Calc_PtOrdered, RVec<float>& dnn_B_DeepAK8Calc_PtOrdered, RVec<int>& dnn_largest_DeepAK8Calc_PtOrdered, RVec<float>& theJetAK8SoftDropCorr_PtOrdered)
{
	TLorentzVector jet_lv;
	std::vector<pair<TLorentzVector,int>> jets_lv;
	float deltaR_leptonicparticle_AK8_PtOrdered = 0;
	bool isLeptonic_W = false;
	bool isLeptonic_t = false;
	for(unsigned int ijet=0; ijet < ak8_pt.size(); ijet++)
	{
		jet_lv.SetPtEtaPhiM(ak8_pt.at(ijet),ak8_eta.at(ijet),ak8_phi.at(ijet),ak8_mass.at(ijet));
		if(isLeptonic == 0)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(Wlv);
			isLeptonic_W = true;
		}
		if(isLeptonic == 1)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(top_lv);
			isLeptonic_t = true;
		}
		// Get 3 highest-pT jets that are not close to t/W (deltaR > .8) and store AK8 index and 4-vector
		if(jets_lv.size() >= 3){continue;}
		if(jet_lv.DeltaR(top_lv) > 0.8 and isLeptonic_t) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
		if(jet_lv.DeltaR(Wlv) > 0.8 and isLeptonic_W) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
	}
	
	float highPtAK8Jet1_SoftDropCorrectedMass = -999;
	float highPtAK8Jet2_SoftDropCorrectedMass = -999;
	float highPtAK8Jet3_SoftDropCorrectedMass = -999;
	
	float Tprime1_DeepAK8_Mass = -999;
	float Tprime2_DeepAK8_Mass = -999;
	float Tprime1_DeepAK8_Pt = -9999;
	float Tprime2_DeepAK8_Pt = -9999;
	float Tprime1_DeepAK8_Eta = 9;
	float Tprime2_DeepAK8_Eta = 9;
	float Tprime1_DeepAK8_Phi = 9;
	float Tprime2_DeepAK8_Phi = 9;
	float Tprime1_DeepAK8_deltaR = -9;
	float Tprime2_DeepAK8_deltaR = -9;
	float Bprime1_DeepAK8_Mass = -999;
	float Bprime2_DeepAK8_Mass = -999;
	float Bprime1_DeepAK8_Pt = -9999;
	float Bprime2_DeepAK8_Pt = -9999;
	float Bprime1_DeepAK8_Eta = 9;
	float Bprime2_DeepAK8_Eta = 9;
	float Bprime1_DeepAK8_Phi = 9;
	float Bprime2_DeepAK8_Phi = 9;
	float Bprime1_DeepAK8_deltaR = -9;
	float Bprime2_DeepAK8_deltaR = -9;
	
	float probSum_DeepAK8_decay = -999;
	float probSum_DeepAK8_four = -999;
	
	float leptonicTprimeJetIDs_DeepAK8 = -1; //Changed into int within .Define()
	float leptonicBprimeJetIDs_DeepAK8 = -1;
	
	bool validTDecay_DeepAK8 = false;
	bool validBDecay_DeepAK8 = false;
	
	RVec<float> hadronicTprimeJetIDs_DeepAK8 (2,0);
	RVec<float> hadronicBprimeJetIDs_DeepAK8 (2,0);
	
	// ----------------------------------------------------------------------------
	// VLQ Decay -- 3 AK8 jets away from leptonic particle
	// ----------------------------------------------------------------------------
	if(jets_lv.size() > 3) {std::cout << "Problem: > 3 AK8s for Tprime reco" << std::endl;}
	int npass_ThreeJets = 0;
	if(jets_lv.size() == 3)
	{
		npass_ThreeJets = npass_ThreeJets + 1;
		// Sums of tag probabilities
		probSum_DeepAK8_decay = 0;
		probSum_DeepAK8_four = 0;
		for (unsigned int i = 0; i < jets_lv.size(); i++)
		{
			// "decay" weighted sum of probabilities including B
			probSum_DeepAK8_decay += dnn_B_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + 2*dnn_W_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + 3*dnn_T_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + 4*dnn_Z_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + 5*dnn_H_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second);
			
			// "four" raw sum of heavy particle probabilities
			probSum_DeepAK8_four += dnn_W_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + dnn_T_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + dnn_Z_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second) + dnn_H_DeepAK8Calc_PtOrdered.at(jets_lv.at(i).second);
		}
		
		// First guess at leptonic T particles
		if(isLeptonic_W) {leptonicTprimeJetIDs_DeepAK8 = 4;}
		if(isLeptonic_t) {leptonicTprimeJetIDs_DeepAK8 = 1;}
		
		// Save masses of the 3 jets for plotting
		highPtAK8Jet1_SoftDropCorrectedMass = theJetAK8SoftDropCorr_PtOrdered.at(jets_lv.at(0).second);
		highPtAK8Jet2_SoftDropCorrectedMass = theJetAK8SoftDropCorr_PtOrdered.at(jets_lv.at(1).second);
		highPtAK8Jet3_SoftDropCorrectedMass = theJetAK8SoftDropCorr_PtOrdered.at(jets_lv.at(2).second);
		
		// ----------------------------------------------------------------------------
		// DeepAK8 SECTION -- TT
		// ----------------------------------------------------------------------------
		
		// get the tags
		float jet1_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(0).second);
		float jet2_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(1).second);
		float jet3_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(2).second);
		// pair up the jet tag with the pT index 0,1,2 and sort by tag (orders J, T, H, Z, W, B)
		std::vector <pair<float,float>> decayJets_DeepAK8;
		decayJets_DeepAK8.push_back(std::make_pair(jet1_DeepAK8,0));
		decayJets_DeepAK8.push_back(std::make_pair(jet2_DeepAK8,1));
		decayJets_DeepAK8.push_back(std::make_pair(jet3_DeepAK8,2));
		std::sort(decayJets_DeepAK8.begin(),decayJets_DeepAK8.end());
		
		// Start forming 4 particle groups
		TLorentzVector Tprime1_DeepAK8_lv;
		TLorentzVector Tprime2_DeepAK8_lv;
		TLorentzVector Bprime1_DeepAK8_lv;
		TLorentzVector Bprime2_DeepAK8_lv;
		validBDecay_DeepAK8 = false;
		validTDecay_DeepAK8 = false;
		if(isLeptonic_t)
		{
			float massDiff1=(top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first).M()-(jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first).M();
			float massDiff2=(top_lv+jets_lv.at(decayJets_DeepAK8.at(2).second).first).M()-(jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first).M();
			float massDiff3=(top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first).M()-(jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first).M();
			
			if(decayJets_DeepAK8.at(0).first==2 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tH bW, BB -> tW bH
				validTDecay_DeepAK8 = true;
				leptonicTprimeJetIDs_DeepAK8 = 2; // assign the H with the leptonic top
				hadronicTprimeJetIDs_DeepAK8 = {4,5};      // assign the b & W as hadronic
				Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first; // decayJets.second gives the jets_lv index to get 4-vec
				Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
				Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				
				validBDecay_DeepAK8 = true;
				leptonicBprimeJetIDs_DeepAK8 = 4; // assign the W with leptonic top
				hadronicBprimeJetIDs_DeepAK8 = {2,5};      // assign bH hadronic
				Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first; // decayJets.second gives the jets_lv index to get 4-vec
				Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
				Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
				Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
			
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==2)
			{ // TTbar --> tH and tH
				validTDecay_DeepAK8 = true;
				leptonicTprimeJetIDs_DeepAK8 = 2; // assign an H with leptonic top
				hadronicTprimeJetIDs_DeepAK8 = {1,2};      // assign tH hadronic
				// options (lepTop + H1) - (T0 + H2) OR (lepTop + H2) - (T0 + H1) checking smallest
				if(massDiff1 < massDiff2)
				{ // (lepTop + H1) wins
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // (lepTop + H2) wins
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				}
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==3)
			{ // TTbar --> tH and tZ
				validTDecay_DeepAK8 = true;
				// options (lepTop + H1) - (T0 + Z2) OR (lepTop + Z2) - (T0 + H1)
				if(massDiff1 < massDiff2)
				{ // (lepTop + H1) wins
					leptonicTprimeJetIDs_DeepAK8 = 2; // tH
					hadronicTprimeJetIDs_DeepAK8 = {1,3}; // tZ
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // (lepTop + Z2) wins
					leptonicTprimeJetIDs_DeepAK8 = 3;
					hadronicTprimeJetIDs_DeepAK8 = {1,2};
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				}
			}
			else if(decayJets_DeepAK8.at(0).first==3 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tZ bW, BB -> tW bZ
				validTDecay_DeepAK8 = true;
				leptonicTprimeJetIDs_DeepAK8 = 3; // tZ
				hadronicTprimeJetIDs_DeepAK8 = {4,5}; // bW
				Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
				Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
				Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				
				validBDecay_DeepAK8 = true;
				leptonicBprimeJetIDs_DeepAK8 = 4; // assign the W with the leptonic top
				hadronicBprimeJetIDs_DeepAK8 = {3,5};      // assign the b & Z as hadronic
				Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first; // decayJets.second gives the jets_lv index to get 4-vec
				Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
				Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
				Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==3 && decayJets_DeepAK8.at(2).first==3)
			{ // TTbar --> tZ tZ
				validTDecay_DeepAK8 = true;
				leptonicTprimeJetIDs_DeepAK8 = 3; // tZ
				hadronicTprimeJetIDs_DeepAK8 = {1,3}; // tZ
				// options (lepTop + Z1) - (T0 + Z2) OR (lepTop + Z2) - (T0 + Z1)
				if(massDiff1 < massDiff2)
				{ // (lepTop + Z1) wins
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // (lepTop + Z2) wins
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				}
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==4)
			{ // BB/XX -> tW tW, jets t W W
				validBDecay_DeepAK8 = true;
				leptonicBprimeJetIDs_DeepAK8 = 4; // tW
				hadronicBprimeJetIDs_DeepAK8 = {1,4}; // tW
				// options (lepTop + W1) - (T0 + W2) OR (lepTop + W2) - (T0 + W1)
				if(massDiff1 < massDiff2)
				{ // (lepTop + W1) wins
					Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // (lepTop + W2) wins
					Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
				}
			}
			
			if(!validTDecay_DeepAK8)
			{ // Not a valid T decay combination
				if(massDiff3 < massDiff1 and massDiff3 < massDiff2)
				{  // lepTop + 0 is best
					leptonicTprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(0).first;
					hadronicTprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(1).first, decayJets_DeepAK8.at(2).first};
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else if(massDiff1 < massDiff2 and massDiff1 < massDiff3)
				{ // lepTop + 1 is best
					leptonicTprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(1).first;
					hadronicTprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(0).first, decayJets_DeepAK8.at(2).first};
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // lepTop + 2 is best
					leptonicTprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(2).first;
					hadronicTprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(0).first, decayJets_DeepAK8.at(1).first};
					Tprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				}
			}
			if(!validBDecay_DeepAK8)
			{ // Not a valid B decay combination
				if(massDiff3 < massDiff1 and massDiff3 < massDiff2)
				{  // lepTop + 0 is best
					leptonicBprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(0).first;
					hadronicBprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(1).first, decayJets_DeepAK8.at(2).first};
					Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else if(massDiff1 < massDiff2 and massDiff1 < massDiff3)
				{ // lepTop + 1 is best
					leptonicBprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(1).first;
					hadronicBprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(0).first, decayJets_DeepAK8.at(2).first};
					Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // lepTop + 2 is best
					leptonicBprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(2).first;
					hadronicBprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(0).first, decayJets_DeepAK8.at(1).first};
					Bprime1_DeepAK8_lv = top_lv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Bprime1_DeepAK8_deltaR = top_lv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				}
			}
		}
		else
		{ // isLeptonic_W
			float massDiff1=(Wlv+jets_lv.at(decayJets_DeepAK8.at(1).second).first).M()-(jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first).M();
			float massDiff2=(Wlv+jets_lv.at(decayJets_DeepAK8.at(2).second).first).M()-(jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first).M();
			float massDiff3=(Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first).M()-(jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first).M();
			
			if(decayJets_DeepAK8.at(0).first==4 && decayJets_DeepAK8.at(1).first==5 && decayJets_DeepAK8.at(2).first==5)
			{ // bW bW
				validTDecay_DeepAK8 = true;
				leptonicTprimeJetIDs_DeepAK8 = 5; // bW
				hadronicTprimeJetIDs_DeepAK8 = {4,5}; // bW
				// options (lepW + b1) - (W0 + b2) OR (lepW + b2) - (W0 + b1)
				if(massDiff1 < massDiff2)
				{ // (lepW + b1) wins
					Tprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // (lepW + b2) wins
					Tprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
				}
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==3 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tZ, BB -> tW bZ
				validTDecay_DeepAK8 = true;
				leptonicTprimeJetIDs_DeepAK8 = 5; // bW
				hadronicTprimeJetIDs_DeepAK8 = {1,3}; // tZ
				Tprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
				Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
				Tprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
				
				validBDecay_DeepAK8 = true;
				leptonicBprimeJetIDs_DeepAK8 = 4; // tW
				hadronicBprimeJetIDs_DeepAK8 = {3,5}; // bZ
				Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
				Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
				Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tH, BB -> tW bH
				validTDecay_DeepAK8 = true;
				leptonicTprimeJetIDs_DeepAK8 = 5; // bW
				hadronicTprimeJetIDs_DeepAK8 = {1,2}; // tH
				Tprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
				Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
				Tprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
				
				validBDecay_DeepAK8 = true;
				leptonicBprimeJetIDs_DeepAK8 = 4; // tW
				hadronicBprimeJetIDs_DeepAK8 = {2,5}; // bH
				Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
				Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
				Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
				Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(2).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==1 && decayJets_DeepAK8.at(2).first==4)
			{ // BB -> tW tW, jets t t W
				validBDecay_DeepAK8 = true;
				leptonicBprimeJetIDs_DeepAK8 = 1; // tW
				hadronicBprimeJetIDs_DeepAK8 = {1,4}; // tW
				// options (lepW + t0) - (W2 + t1) OR (lepW + t1) - (W2 + t0)
				if(massDiff3 < massDiff1)
				{ // (lepW + t0) wins
					Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // (lepW + t1) wins
					Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
			}
			if(!validTDecay_DeepAK8)
			{ // not a valid grouping
				if(massDiff3 < massDiff1 and massDiff3 < massDiff2)
				{ // lepW + 0 wins
					leptonicTprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(0).first;
					hadronicTprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(1).first, decayJets_DeepAK8.at(2).first};
					Tprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else if(massDiff1 < massDiff2 and massDiff1 < massDiff3)
				{ // lepW + 1 wins
					leptonicTprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(1).first;
					hadronicTprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(0).first, decayJets_DeepAK8.at(2).first};
					Tprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // lepW + 2 wins
					leptonicTprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(2).first;
					hadronicTprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(0).first, decayJets_DeepAK8.at(1).first};
					Tprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Tprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Tprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Tprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
			}
			if(!validBDecay_DeepAK8)
			{ // not a valid grouping
				if(massDiff3 < massDiff1 and massDiff3 < massDiff2)
				{ // lepW + 0 wins
					leptonicBprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(0).first;
					hadronicBprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(1).first, decayJets_DeepAK8.at(2).first};
					Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(0).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(1).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(0).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(1).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else if(massDiff1 < massDiff2 and massDiff1 < massDiff3)
				{ // lepW + 1 wins
					leptonicBprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(1).first;
					hadronicBprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(0).first, decayJets_DeepAK8.at(2).first};
					Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(1).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
				else
				{ // lepW + 2 wins
					leptonicBprimeJetIDs_DeepAK8 = decayJets_DeepAK8.at(2).first;
					hadronicBprimeJetIDs_DeepAK8 = {decayJets_DeepAK8.at(0).first, decayJets_DeepAK8.at(1).first};
					Bprime1_DeepAK8_lv = Wlv+jets_lv.at(decayJets_DeepAK8.at(2).second).first;
					Bprime2_DeepAK8_lv = jets_lv.at(decayJets_DeepAK8.at(0).second).first+jets_lv.at(decayJets_DeepAK8.at(1).second).first;
					Bprime1_DeepAK8_deltaR = Wlv.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
					Bprime2_DeepAK8_deltaR = jets_lv.at(decayJets_DeepAK8.at(0).second).first.DeltaR(jets_lv.at(decayJets_DeepAK8.at(2).second).first);
				}
			}
		}
		if(Tprime1_DeepAK8_lv.M() != -999)
		{
			Tprime1_DeepAK8_Mass = Tprime1_DeepAK8_lv.M();
			Tprime2_DeepAK8_Mass = Tprime2_DeepAK8_lv.M();
			Tprime1_DeepAK8_Pt = Tprime1_DeepAK8_lv.Pt();
			Tprime2_DeepAK8_Pt = Tprime2_DeepAK8_lv.Pt();
			Tprime1_DeepAK8_Eta = Tprime1_DeepAK8_lv.Eta();
			Tprime2_DeepAK8_Eta = Tprime2_DeepAK8_lv.Eta();
			Tprime1_DeepAK8_Phi = Tprime1_DeepAK8_lv.Phi();
			Tprime2_DeepAK8_Phi = Tprime2_DeepAK8_lv.Phi();
		}
		if(Bprime1_DeepAK8_lv.M() != -999)
		{
			Bprime1_DeepAK8_Mass = Bprime1_DeepAK8_lv.M();
			Bprime2_DeepAK8_Mass = Bprime2_DeepAK8_lv.M();
			Bprime1_DeepAK8_Pt = Bprime1_DeepAK8_lv.Pt();
			Bprime2_DeepAK8_Pt = Bprime2_DeepAK8_lv.Pt();
			Bprime1_DeepAK8_Eta = Bprime1_DeepAK8_lv.Eta();
			Bprime2_DeepAK8_Eta = Bprime2_DeepAK8_lv.Eta();
			Bprime1_DeepAK8_Phi = Bprime1_DeepAK8_lv.Phi();
			Bprime2_DeepAK8_Phi = Bprime2_DeepAK8_lv.Phi();
		}
	}
	RVec<float> TandBPrimeVec = {Tprime1_DeepAK8_Mass,Tprime2_DeepAK8_Mass,Tprime1_DeepAK8_Pt,Tprime2_DeepAK8_Pt,Tprime1_DeepAK8_Eta,Tprime2_DeepAK8_Eta,Tprime1_DeepAK8_Phi,Tprime2_DeepAK8_Phi,Tprime1_DeepAK8_deltaR,Tprime2_DeepAK8_deltaR,Bprime1_DeepAK8_Mass,Bprime2_DeepAK8_Mass,Bprime1_DeepAK8_Pt,Bprime2_DeepAK8_Pt,Bprime1_DeepAK8_Eta,Bprime2_DeepAK8_Eta,Bprime1_DeepAK8_Phi,Bprime2_DeepAK8_Phi,Bprime1_DeepAK8_deltaR,Bprime2_DeepAK8_deltaR,highPtAK8Jet1_SoftDropCorrectedMass,highPtAK8Jet2_SoftDropCorrectedMass,highPtAK8Jet3_SoftDropCorrectedMass,leptonicTprimeJetIDs_DeepAK8,leptonicBprimeJetIDs_DeepAK8,hadronicTprimeJetIDs_DeepAK8[0],hadronicTprimeJetIDs_DeepAK8[1],hadronicBprimeJetIDs_DeepAK8[0],hadronicBprimeJetIDs_DeepAK8[1]};
	return TandBPrimeVec;
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
// 									TTAGGING FXN
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
RVec<int> three_jet_TTag(TLorentzVector top_lv, TLorentzVector Wlv, int isLeptonic, RVec<float>& ak8_pt, RVec<float>& ak8_eta, RVec<float>& ak8_phi, RVec<float>& ak8_mass, RVec<int>& dnn_largest_DeepAK8Calc_PtOrdered, int hadronicTprimeJetIDs1_DeepAK8, int hadronicTprimeJetIDs2_DeepAK8)
{
	RVec<int> validTTagged (2,0);
	int val = -99;
	int tag = -99;
	TString tagCut = "";
	TLorentzVector jet_lv;
	std::vector<pair<TLorentzVector,int>> jets_lv;
	float deltaR_leptonicparticle_AK8_PtOrdered = 0;
	bool isLeptonic_W = false;
	bool isLeptonic_t = false;
	for(unsigned int ijet=0; ijet < ak8_pt.size(); ijet++)
	{
		jet_lv.SetPtEtaPhiM(ak8_pt.at(ijet),ak8_eta.at(ijet),ak8_phi.at(ijet),ak8_mass.at(ijet));
		if(isLeptonic == 0)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(Wlv);
			isLeptonic_W = true;
		}
		if(isLeptonic == 1)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(top_lv);
			isLeptonic_t = true;
		}
		// Get 3 highest-pT jets that are not close to t/W (deltaR > .8) and store AK8 index and 4-vector
		if(jets_lv.size() >= 3){continue;}
		if(jet_lv.DeltaR(top_lv) > 0.8 and isLeptonic_t) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
		if(jet_lv.DeltaR(Wlv) > 0.8 and isLeptonic_W) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
	}
	
	bool validTDecay_DeepAK8 = false; // [0,y]
	bool taggedBWBW_DeepAK8 = false; // [x,0]
	bool taggedTHBW_DeepAK8 = false; // [x,1]
	bool taggedTHTH_DeepAK8 = false; // [x,2]
	bool taggedTZBW_DeepAK8 = false; // [x,3]
	bool taggedTZTH_DeepAK8 = false; // [x,4]
	bool taggedTZTZ_DeepAK8 = false; // [x,5]
	
	// ----------------------------------------------------------------------------
	// VLQ Decay -- 3 AK8 jets away from leptonic particle
	// ----------------------------------------------------------------------------
	if(jets_lv.size() > 3) {std::cout << "Problem: > 3 AK8s for Tprime reco" << std::endl;}
	if(jets_lv.size() == 3)
	{
		// ----------------------------------------------------------------------------
		// DeepAK8 SECTION -- TT
		// ----------------------------------------------------------------------------
		
		// get the tags
		int jet1_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(0).second);
		int jet2_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(1).second);
		int jet3_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(2).second);
		// pair up the jet tag with the pT index 0,1,2 and sort by tag (orders J, T, H, Z, W, B)
		std::vector <pair<int,int>> decayJets_DeepAK8;
		decayJets_DeepAK8.push_back(std::make_pair(jet1_DeepAK8,0));
		decayJets_DeepAK8.push_back(std::make_pair(jet2_DeepAK8,1));
		decayJets_DeepAK8.push_back(std::make_pair(jet3_DeepAK8,2));
		std::sort(decayJets_DeepAK8.begin(),decayJets_DeepAK8.end());
		
		// Start forming 4 particle groups
		validTDecay_DeepAK8 = false;
		if(isLeptonic_t)
		{
			if(decayJets_DeepAK8.at(0).first==2 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tH bW, BB -> tW bH
				validTDecay_DeepAK8 = true; val = 0;
				taggedTHBW_DeepAK8 = true; tag = 1;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==2)
			{ // TTbar --> tH and tH
				validTDecay_DeepAK8 = true; val = 0;
				taggedTHTH_DeepAK8 = true; tag = 2;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==3)
			{ // TTbar --> tH and tZ
				validTDecay_DeepAK8 = true; val = 0;
				taggedTZTH_DeepAK8 = true; tag = 4;
			}
			else if(decayJets_DeepAK8.at(0).first==3 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tZ bW, BB -> tW bZ
				validTDecay_DeepAK8 = true; val = 0;
				taggedTZBW_DeepAK8 = true; tag = 3;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==3 && decayJets_DeepAK8.at(2).first==3)
			{ // TTbar --> tZ tZ
				validTDecay_DeepAK8 = true; val = 0;
				taggedTZTZ_DeepAK8 = true; tag = 5;
			}
		}
		else
		{ // isLeptonic_W
			if(decayJets_DeepAK8.at(0).first==4 && decayJets_DeepAK8.at(1).first==5 && decayJets_DeepAK8.at(2).first==5)
			{ // bW bW
				validTDecay_DeepAK8 = true; val = 0;
				taggedBWBW_DeepAK8 = true; tag = 0;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==3 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tZ, BB -> tW bZ
				validTDecay_DeepAK8 = true; val = 0;
				taggedTZBW_DeepAK8 = true; tag = 3;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tH, BB -> tW bH
				validTDecay_DeepAK8 = true; val = 0;
				taggedTHBW_DeepAK8 = true; tag = 1;
			}
		}
	}
	if(!validTDecay_DeepAK8)
	{
		// signal categories for only hadronic VLQ valid
		if(hadronicTprimeJetIDs1_DeepAK8 == 1 && hadronicTprimeJetIDs2_DeepAK8 == 2) {tag = 9;} // 'tH'
		else if(hadronicTprimeJetIDs1_DeepAK8 == 1 && hadronicTprimeJetIDs2_DeepAK8 == 3) {tag = 10;} // 'tZ'
		else if(hadronicTprimeJetIDs1_DeepAK8 == 4 && hadronicTprimeJetIDs2_DeepAK8 == 5) {tag = 11;} // 'bW'
		// if whichSig == 'TT' and 'tH' not in tag and 'tZ' not in tag and 'bW' not in tag:
		else{tag = -5;}
	}
	validTTagged = {val,tag};
	return validTTagged;
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
//  									 BTAGGING FXN 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------
RVec<int> three_jet_BTag(TLorentzVector top_lv, TLorentzVector Wlv, int isLeptonic, RVec<float>& ak8_pt, RVec<float>& ak8_eta, RVec<float>& ak8_phi, RVec<float>& ak8_mass, RVec<int>& dnn_largest_DeepAK8Calc_PtOrdered, int hadronicBprimeJetIDs1_DeepAK8, int hadronicBprimeJetIDs2_DeepAK8)
{
	RVec<int> validBTagged (2,0);
	int val = -99;
	int tag = -99;
	TLorentzVector jet_lv;
	std::vector<pair<TLorentzVector,int>> jets_lv;
	float deltaR_leptonicparticle_AK8_PtOrdered = 0;
	bool isLeptonic_W = false;
	bool isLeptonic_t = false;
	for(unsigned int ijet=0; ijet < ak8_pt.size(); ijet++)
	{
		jet_lv.SetPtEtaPhiM(ak8_pt.at(ijet),ak8_eta.at(ijet),ak8_phi.at(ijet),ak8_mass.at(ijet));
		if(isLeptonic == 0)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(Wlv);
			isLeptonic_W = true;
		}
		if(isLeptonic == 1)
		{
			deltaR_leptonicparticle_AK8_PtOrdered = jet_lv.DeltaR(top_lv);
			isLeptonic_t = true;
		}
		// Get 3 highest-pT jets that are not close to t/W (deltaR > .8) and store AK8 index and 4-vector
		if(jets_lv.size() >= 3){continue;}
		if(jet_lv.DeltaR(top_lv) > 0.8 and isLeptonic_t) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
		if(jet_lv.DeltaR(Wlv) > 0.8 and isLeptonic_W) {jets_lv.push_back(std::make_pair(jet_lv,ijet));}
	}
	
	bool validBDecay_DeepAK8 = false; // [1,y]
	bool taggedTZBW_DeepAK8 = false; // [x,3] T and B
	bool taggedTWTW_DeepAK8 = false; // [x,6] B
	bool taggedBZTW_DeepAK8 = false; // [x,7] B
	bool taggedBHTW_DeepAK8 = false; // [x,8] B
	
	// ----------------------------------------------------------------------------
	// VLQ Decay -- 3 AK8 jets away from leptonic particle
	// ----------------------------------------------------------------------------
	if(jets_lv.size() > 3) {std::cout << "Problem: > 3 AK8s for Tprime reco" << std::endl;}
	if(jets_lv.size() == 3)
	{
		// ----------------------------------------------------------------------------
		// DeepAK8 SECTION -- TT
		// ----------------------------------------------------------------------------
		
		// get the tags
		int jet1_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(0).second);
		int jet2_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(1).second);
		int jet3_DeepAK8 = dnn_largest_DeepAK8Calc_PtOrdered.at(jets_lv.at(2).second);
		// pair up the jet tag with the pT index 0,1,2 and sort by tag (orders J, T, H, Z, W, B)
		std::vector <pair<int,int>> decayJets_DeepAK8;
		decayJets_DeepAK8.push_back(std::make_pair(jet1_DeepAK8,0));
		decayJets_DeepAK8.push_back(std::make_pair(jet2_DeepAK8,1));
		decayJets_DeepAK8.push_back(std::make_pair(jet3_DeepAK8,2));
		std::sort(decayJets_DeepAK8.begin(),decayJets_DeepAK8.end());
		
		// Start forming 4 particle groups
		validBDecay_DeepAK8 = false;
		if(isLeptonic_t)
		{
			if(decayJets_DeepAK8.at(0).first==2 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tH bW, BB -> tW bH
				validBDecay_DeepAK8 = true; val = 1;
				taggedBHTW_DeepAK8 = true; tag = 8;
			}
			else if(decayJets_DeepAK8.at(0).first==3 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> tZ bW, BB -> tW bZ
				validBDecay_DeepAK8 = true; val = 1;
				taggedBZTW_DeepAK8 = true; tag = 7;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==4 && decayJets_DeepAK8.at(2).first==4)
			{ // BB/XX -> tW tW, jets t W W
				validBDecay_DeepAK8 = true; val = 1;
				taggedTWTW_DeepAK8 = true; tag = 6;
			}
		}
		else
		{ // isLeptonic_W
			if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==3 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tZ, BB -> tW bZ
				validBDecay_DeepAK8 = true; val = 1;
				taggedBZTW_DeepAK8 = true; tag = 7;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==2 && decayJets_DeepAK8.at(2).first==5)
			{ // TT -> bW tH, BB -> tW bH
				validBDecay_DeepAK8 = true; val = 1;
				taggedBHTW_DeepAK8 = true; tag = 8;
			}
			else if(decayJets_DeepAK8.at(0).first==1 && decayJets_DeepAK8.at(1).first==1 && decayJets_DeepAK8.at(2).first==4)
			{ // BB -> tW tW, jets t t W
				validBDecay_DeepAK8 = true; val = 1;
				taggedTWTW_DeepAK8 = true; tag = 6;
			}
		}
	}
	if(!validBDecay_DeepAK8)
	{
		// signal categories for only hadronic VLQ valid
		if(hadronicBprimeJetIDs1_DeepAK8 == 1 && hadronicBprimeJetIDs2_DeepAK8 == 4) {tag = 12;} // 'tW'
		else if(hadronicBprimeJetIDs1_DeepAK8 == 3 && hadronicBprimeJetIDs2_DeepAK8 == 5) {tag = 13;} // 'bZ'
		else if(hadronicBprimeJetIDs1_DeepAK8 == 2 && hadronicBprimeJetIDs2_DeepAK8 == 5){tag = 14;} // 'bH'
		// if(whichSig == "BB" and 'bH' not in tag and 'bZ' not in tag and 'tW' not in tag)
		else{tag = -5;}
	}
	validBTagged = {val,tag};
	return validBTagged;
};
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------



