Events->Draw("FatJet_particleNetMD_QCD","(FatJet_particleNetMD_Xbb/(FatJet_particleNetMD_Xbb+FatJet_particleNetMD_QCD))==Max$(FatJet_particleNetMD_Xbb/(FatJet_particleNetMD_Xbb+FatJet_particleNetMD_QCD))&&Sum$(FatJet_pt>250 && abs(FatJet_eta) < 2.5 && FatJet_msoftdrop > 40)>=3")
python stack_plots.py -b
