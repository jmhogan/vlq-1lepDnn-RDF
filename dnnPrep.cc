
// --------------------------------------------
// 	   MINDR & PTREL CALCULATOR
// --------------------------------------------
ROOT::VecOps::RVec<float> minDR_ptRel_lead_calc(ROOT::VecOps::RVec<float>& jet_pt, ROOT::VecOps::RVec<float>& jet_eta, ROOT::VecOps::RVec<float>& jet_phi, ROOT::VecOps::RVec<float>& jet_mass, TLorentzVector lepton_lv)
{
	TLorentzVector jet_TLV, leadJet;
	float deltaR_lepJets = 0;
	float ptRel_lepJets = 0;
	float minDR_lepJets = 1000;
	float minDR_leadJetotherJet = 1000;
	leadJet.SetPtEtaPhiM(jet_pt[0],jet_eta[0],jet_phi[0],jet_mass[0]);
	if(jet_pt.size() < 1) {minDR_lepJets = -99.0;}
	if(jet_pt.size() < 2) {minDR_leadJetotherJet = -99.0;}
	for(int i = 0; i < jet_pt.size(); i++)
	{
		jet_TLV.SetPtEtaPhiM(jet_pt[i],jet_eta[i],jet_phi[i],jet_mass[i]);
		deltaR_lepJets = lepton_lv.DeltaR(jet_TLV);
		if(deltaR_lepJets < minDR_lepJets)
		{
			minDR_lepJets = lepton_lv.DeltaR(jet_TLV);
			ptRel_lepJets = lepton_lv.P()*(jet_TLV.Vect().Cross(lepton_lv.Vect()).Mag()/jet_TLV.P()/lepton_lv.P());
		}
		if(i > 0)
		{
			float tempDR = leadJet.DeltaR(jet_TLV);
			if(tempDR < minDR_leadJetotherJet){minDR_leadJetotherJet = tempDR;}
		}
	}
	ROOT::VecOps::RVec<float> returnVec = {minDR_lepJets,ptRel_lepJets,minDR_leadJetotherJet};
	return returnVec;
};
