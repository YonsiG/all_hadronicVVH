#ifndef _MAKEHIST_H_
#define _MAKEHIST_H_

#include <iostream>
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"
#include "TLegend.h"

class makeHists
{
public:
  TFile *hf;

  TH1D *weight_Scale;
  TH1D *cutflow;

  TH1D *fatjet_btagger;

  TH1D *HiggsM_QJ_DeltaR0_01LL;
  TH1D *HiggsM_QJ_DeltaR0_01MM;
  TH1D *HiggsM_QJ_DeltaR0_01TT;

  TH1D *HiggsM_QJ_DeltaR01_025LL;
  TH1D *HiggsM_QJ_DeltaR01_025MM;
  TH1D *HiggsM_QJ_DeltaR01_025TT;

  TH1D *HiggsM_QJ_DeltaR025_05LL;
  TH1D *HiggsM_QJ_DeltaR025_05MM;
  TH1D *HiggsM_QJ_DeltaR025_05TT;

  TH2D *QJ_DeltaR_vs_HiggsMLL;
  TH2D *QJ_DeltaR_vs_HiggsMMM;
  TH2D *QJ_DeltaR_vs_HiggsMTT;

  TH2D *QJ_DeltaR1_vs_HiggsMLL;
  TH2D *QJ_DeltaR1_vs_HiggsMMM;
  TH2D *QJ_DeltaR1_vs_HiggsMTT;

  TH2D *QJ_DeltaR2_vs_HiggsMLL;
  TH2D *QJ_DeltaR2_vs_HiggsMMM;
  TH2D *QJ_DeltaR2_vs_HiggsMTT;

  TH2D *QJ_DeltaR1_vs_QJ_DeltaR2LL;
  TH2D *QJ_DeltaR1_vs_QJ_DeltaR2MM;
  TH2D *QJ_DeltaR1_vs_QJ_DeltaR2TT;

  TH1D *QJ_DeltaRLL;
  TH1D *QJ_DeltaRMM;
  TH1D *QJ_DeltaRTT;

  TH1D *QJ_DeltaR1LL;
  TH1D *QJ_DeltaR1MM;
  TH1D *QJ_DeltaR1TT;

  TH1D *QJ_DeltaR2LL;
  TH1D *QJ_DeltaR2MM;
  TH1D *QJ_DeltaR2TT;

  TH1D *total_Higgspt;
  TH1D *cut_HiggsptLL;
  TH1D *cut_HiggsptMM;
  TH1D *cut_HiggsptTT;
  TH1D *total_quarkDeltaR;
  TH1D *cut_quarkDeltaRLL;
  TH1D *cut_quarkDeltaRMM;
  TH1D *cut_quarkDeltaRTT;
  TH1D *cut_Higgspt_2passLL;
  TH1D *cut_Higgspt_2passMM;
  TH1D *cut_Higgspt_2passTT;
  TH1D *cut_Higgspt_1passLL;
  TH1D *cut_Higgspt_1passMM;
  TH1D *cut_Higgspt_1passTT;
  TH1D *cut_Higgspt_0passLL;
  TH1D *cut_Higgspt_0passMM;
  TH1D *cut_Higgspt_0passTT;

  TH1D *recoHiggs_MassLL;
  TH1D *recoHiggs_msoftdropLL;
  TH1D *recoHiggs_PtLL;
  TH1D *recoHiggs_EtaLL;
  TH1D *recoHiggs_Eta_2passLL;
  TH1D *recoHiggs_Eta_1passLL;
  TH1D *recoHiggs_Eta_0passLL;
  TH1D *recoHiggs_Pt_2passLL;
  TH1D *recoHiggs_Pt_1passLL;
  TH1D *recoHiggs_Pt_0passLL;
  TH1D *recoHiggs_msoftdrop_2passLL;
  TH1D *recoHiggs_msoftdrop_1passLL;
  TH1D *recoHiggs_msoftdrop_0passLL;

  TH1D *recoHiggs_MassMM;
  TH1D *recoHiggs_msoftdropMM;
  TH1D *recoHiggs_PtMM;
  TH1D *recoHiggs_EtaMM;
  TH1D *recoHiggs_Eta_2passMM;
  TH1D *recoHiggs_Eta_1passMM;
  TH1D *recoHiggs_Eta_0passMM;
  TH1D *recoHiggs_Pt_2passMM;
  TH1D *recoHiggs_Pt_1passMM;
  TH1D *recoHiggs_Pt_0passMM;
  TH1D *recoHiggs_msoftdrop_2passMM;
  TH1D *recoHiggs_msoftdrop_1passMM;
  TH1D *recoHiggs_msoftdrop_0passMM;

  TH1D *recoHiggs_MassTT;
  TH1D *recoHiggs_msoftdropTT;
  TH1D *recoHiggs_PtTT;
  TH1D *recoHiggs_EtaTT;
  TH1D *recoHiggs_Eta_2passTT;
  TH1D *recoHiggs_Eta_1passTT;
  TH1D *recoHiggs_Eta_0passTT;
  TH1D *recoHiggs_Pt_2passTT;
  TH1D *recoHiggs_Pt_1passTT;
  TH1D *recoHiggs_Pt_0passTT;
  TH1D *recoHiggs_msoftdrop_2passTT;
  TH1D *recoHiggs_msoftdrop_1passTT;
  TH1D *recoHiggs_msoftdrop_0passTT;

  TH1D *nBquark;
  TH1D *double_Bquark_Pt;
  TH1D *Bquark_Pt;
  TH1D *anti_Bquark_Pt;
  TH1D *leading_Bquark_Pt;
  TH1D *sub_Bquark_Pt;
  TH1D *leading_Bquark_Eta;
  TH1D *sub_Bquark_Eta;
  TH1D *Bquark_DeltaR;
  TH2D *double_Bquark_Pt_vs_Higgs_Pt;
  TH2D *leading_Bquark_Pt_vs_Higgs_Pt;
  TH2D *sub_Bquark_Pt_vs_Higgs_Pt;
  TH2D *Bquark_DeltaR_vs_Higgs_Pt;

  //  Th1D *leading_Bquark_Phi;
  //  Th1D *sub_Bquark_Phi;
  void createHists(const char *fileName);
  void saveHists();
};
#endif
