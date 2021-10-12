#include "makeHists.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TLegend.h"

using namespace std;
void makeHists::createHists(const char *fileName)
{
    hf = new TFile(fileName, "RECREATE");

    weight_Scale = new TH1D("weight_Scale", "weight_Scale", 1, 0, 1);

    cutflow = new TH1D("cutflow", "cutflow", 10, 0, 10);

    fatjet_btagger = new TH1D("fatjet_btagger", "fatjet_btagger", 100, 0, 1);
    fatjet_btagger->Sumw2();

    HiggsM_QJ_DeltaR0_01LL = new TH1D("HiggsM_QJ_DeltaR0_01LL", "HiggsM_QJ_DeltaR0_01LL", 100, 75, 175);
    HiggsM_QJ_DeltaR0_01LL->Sumw2();

    HiggsM_QJ_DeltaR0_01MM = new TH1D("HiggsM_QJ_DeltaR0_01MM", "HiggsM_QJ_DeltaR0_01MM", 100, 75, 175);
    HiggsM_QJ_DeltaR0_01MM->Sumw2();

    HiggsM_QJ_DeltaR0_01TT = new TH1D("HiggsM_QJ_DeltaR0_01TT", "HiggsM_QJ_DeltaR0_01TT", 100, 75, 175);
    HiggsM_QJ_DeltaR0_01TT->Sumw2();

    HiggsM_QJ_DeltaR01_025LL = new TH1D("HiggsM_QJ_DeltaR01_025LL", "HiggsM_QJ_DeltaR01_025LL", 100, 75, 175);
    HiggsM_QJ_DeltaR01_025LL->Sumw2();

    HiggsM_QJ_DeltaR01_025MM = new TH1D("HiggsM_QJ_DeltaR01_025MM", "HiggsM_QJ_DeltaR01_025MM", 100, 75, 175);
    HiggsM_QJ_DeltaR01_025MM->Sumw2();

    HiggsM_QJ_DeltaR01_025TT = new TH1D("HiggsM_QJ_DeltaR01_025TT", "HiggsM_QJ_DeltaR01_025TT", 100, 75, 175);
    HiggsM_QJ_DeltaR01_025TT->Sumw2();

    HiggsM_QJ_DeltaR025_05LL = new TH1D("HiggsM_QJ_DeltaR025_05LL", "HiggsM_QJ_DeltaR025_05LL", 100, 75, 175);
    HiggsM_QJ_DeltaR025_05LL->Sumw2();

    HiggsM_QJ_DeltaR025_05MM = new TH1D("HiggsM_QJ_DeltaR025_05MM", "HiggsM_QJ_DeltaR025_05MM", 100, 75, 175);
    HiggsM_QJ_DeltaR025_05MM->Sumw2();

    HiggsM_QJ_DeltaR025_05TT = new TH1D("HiggsM_QJ_DeltaR025_05TT", "HiggsM_QJ_DeltaR025_05TT", 100, 75, 175);
    HiggsM_QJ_DeltaR025_05TT->Sumw2();

    QJ_DeltaR_vs_HiggsMLL = new TH2D("QJ_DeltaR_vs_HiggsMLL", "QJ_DeltaR_vs_HiggsMLL", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR_vs_HiggsMLL->Sumw2();

    QJ_DeltaR_vs_HiggsMMM = new TH2D("QJ_DeltaR_vs_HiggsMMM", "QJ_DeltaR_vs_HiggsMMM", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR_vs_HiggsMMM->Sumw2();

    QJ_DeltaR_vs_HiggsMTT = new TH2D("QJ_DeltaR_vs_HiggsMTT", "QJ_DeltaR_vs_HiggsMTT", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR_vs_HiggsMTT->Sumw2();

    QJ_DeltaR1_vs_HiggsMLL = new TH2D("QJ_DeltaR1_vs_HiggsMLL", "QJ_DeltaR1_vs_HiggsMLL", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR1_vs_HiggsMLL->Sumw2();

    QJ_DeltaR1_vs_HiggsMMM = new TH2D("QJ_DeltaR1_vs_HiggsMMM", "QJ_DeltaR1_vs_HiggsMMM", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR1_vs_HiggsMMM->Sumw2();

    QJ_DeltaR1_vs_HiggsMTT = new TH2D("QJ_DeltaR1_vs_HiggsMTT", "QJ_DeltaR1_vs_HiggsMTT", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR1_vs_HiggsMTT->Sumw2();

    QJ_DeltaR2_vs_HiggsMLL = new TH2D("QJ_DeltaR2_vs_HiggsMLL", "QJ_DeltaR2_vs_HiggsMLL", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR2_vs_HiggsMLL->Sumw2();

    QJ_DeltaR2_vs_HiggsMMM = new TH2D("QJ_DeltaR2_vs_HiggsMMM", "QJ_DeltaR2_vs_HiggsMMM", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR2_vs_HiggsMMM->Sumw2();

    QJ_DeltaR2_vs_HiggsMTT = new TH2D("QJ_DeltaR2_vs_HiggsMTT", "QJ_DeltaR2_vs_HiggsMTT", 50, 0, 0.5, 100, 75, 175);
    QJ_DeltaR2_vs_HiggsMTT->Sumw2();

    QJ_DeltaR1_vs_QJ_DeltaR2LL = new TH2D("QJ_DeltaR1_vs_QJ_DeltaR2LL", "QJ_DeltaR1_vs_QJ_DeltaR2LL", 50, 0, 0.5, 50, 0, 0.5);
    QJ_DeltaR1_vs_QJ_DeltaR2LL->Sumw2();

    QJ_DeltaR1_vs_QJ_DeltaR2MM = new TH2D("QJ_DeltaR1_vs_QJ_DeltaR2MM", "QJ_DeltaR1_vs_QJ_DeltaR2MM", 50, 0, 0.5, 50, 0, 0.5);
    QJ_DeltaR1_vs_QJ_DeltaR2MM->Sumw2();

    QJ_DeltaR1_vs_QJ_DeltaR2TT = new TH2D("QJ_DeltaR1_vs_QJ_DeltaR2TT", "QJ_DeltaR1_vs_QJ_DeltaR2TT", 50, 0, 0.5, 50, 0, 0.5);
    QJ_DeltaR1_vs_QJ_DeltaR2TT->Sumw2();

    QJ_DeltaRLL = new TH1D("QJ_DeltaRLL", "QJ_DeltaRLL", 50, 0, 2);
    QJ_DeltaRLL->Sumw2();

    QJ_DeltaRMM = new TH1D("QJ_DeltaRMM", "QJ_DeltaRMM", 50, 0, 2);
    QJ_DeltaRMM->Sumw2();

    QJ_DeltaRTT = new TH1D("QJ_DeltaRTT", "QJ_DeltaRTT", 50, 0, 2);
    QJ_DeltaRTT->Sumw2();

    QJ_DeltaR1LL = new TH1D("QJ_DeltaR1LL", "QJ_DeltaR1LL", 50, 0, 2);
    QJ_DeltaR1LL->Sumw2();

    QJ_DeltaR1MM = new TH1D("QJ_DeltaR1MM", "QJ_DeltaR1MM", 50, 0, 2);
    QJ_DeltaR1MM->Sumw2();

    QJ_DeltaR1TT = new TH1D("QJ_DeltaR1TT", "QJ_DeltaR1TT", 50, 0, 2);
    QJ_DeltaR1TT->Sumw2();

    QJ_DeltaR2LL = new TH1D("QJ_DeltaR2LL", "QJ_DeltaR2LL", 50, 0, 2);
    QJ_DeltaR2LL->Sumw2();

    QJ_DeltaR2MM = new TH1D("QJ_DeltaR2MM", "QJ_DeltaR2MM", 50, 0, 2);
    QJ_DeltaR2MM->Sumw2();

    QJ_DeltaR2TT = new TH1D("QJ_DeltaR2TT", "QJ_DeltaR2TT", 50, 0, 2);
    QJ_DeltaR2TT->Sumw2();

    recoHiggs_MassLL = new TH1D("recoHiggs_MassLL", "recoHiggs_MassLL", 100, 75, 175);
    recoHiggs_MassLL->Sumw2();

    recoHiggs_msoftdropLL = new TH1D("recoHiggs_msoftdropLL", "recoHiggs_msoftdropLL", 100, 75, 175);
    recoHiggs_msoftdropLL->Sumw2();

    recoHiggs_PtLL = new TH1D("recoHiggs_PtLL", "recoHiggs_PtLL", 100, 0, 2000);
    recoHiggs_PtLL->Sumw2();

    recoHiggs_EtaLL = new TH1D("recoHiggs_EtaLL", "recoHiggs_EtaLL", 100, -2.5, 2.5);
    recoHiggs_EtaLL->Sumw2();

    double ybin[31];
    for (int ibin = 0; ibin < 20; ibin++)
    {
        ybin[ibin] = ibin * 50;
    }
    for (int ibin = 20; ibin < 31; ibin++)
    {
        ybin[ibin] = ibin * 100 - 1000;
    }

    total_Higgspt = new TH1D("total_Higgspt", "total_Higgspt", 30, ybin);
    total_Higgspt->Sumw2();

    cut_HiggsptLL = new TH1D("cut_HiggsptLL", "cut_HiggsptLL", 30, ybin);
    cut_HiggsptLL->Sumw2();

    cut_Higgspt_2passLL = new TH1D("cut_Higgspt_2passLL", "cut_Higgspt_2passLL", 30, ybin);
    cut_Higgspt_2passLL->Sumw2();

    cut_Higgspt_1passLL = new TH1D("cut_Higgspt_1passLL", "cut_Higgspt_1passLL", 30, ybin);
    cut_Higgspt_1passLL->Sumw2();

    cut_Higgspt_0passLL = new TH1D("cut_Higgspt_0passLL", "cut_Higgspt_0passLL", 30, ybin);
    cut_Higgspt_0passLL->Sumw2();

    recoHiggs_Eta_2passLL = new TH1D("recoHiggs_Eta_2passLL", "recoHiggs_Eta_2passLL", 100, -2.5, 2.5);
    recoHiggs_Eta_2passLL->Sumw2();

    recoHiggs_Eta_1passLL = new TH1D("recoHiggs_Eta_1passLL", "recoHiggs_Eta_1passLL", 100, -2.5, 2.5);
    recoHiggs_Eta_1passLL->Sumw2();

    recoHiggs_Eta_0passLL = new TH1D("recoHiggs_Eta_0passLL", "recoHiggs_Eta_0passLL", 100, -2.5, 2.5);
    recoHiggs_Eta_0passLL->Sumw2();

    recoHiggs_Pt_2passLL = new TH1D("recoHiggs_Pt_2passLL", "recoHiggs_Pt_2passLL", 100, 0, 2000);
    recoHiggs_Pt_2passLL->Sumw2();

    recoHiggs_Pt_1passLL = new TH1D("recoHiggs_Pt_1passLL", "recoHiggs_Pt_1passLL", 100, 0, 2000);
    recoHiggs_Pt_1passLL->Sumw2();

    recoHiggs_Pt_0passLL = new TH1D("recoHiggs_Pt_0passLL", "recoHiggs_Pt_0passLL", 100, 0, 2000);
    recoHiggs_Pt_0passLL->Sumw2();

    recoHiggs_msoftdrop_2passLL = new TH1D("recoHiggs_msoftdrop_2passLL", "recoHiggs_msoftdrop_2passLL", 100, 75, 175);
    recoHiggs_msoftdrop_2passLL->Sumw2();

    recoHiggs_msoftdrop_1passLL = new TH1D("recoHiggs_msoftdrop_1passLL", "recoHiggs_msoftdrop_1passLL", 100, 75, 175);
    recoHiggs_msoftdrop_1passLL->Sumw2();

    recoHiggs_msoftdrop_0passLL = new TH1D("recoHiggs_msoftdrop_0passLL", "recoHiggs_msoftdrop_0passLL", 100, 75, 175);
    recoHiggs_msoftdrop_0passLL->Sumw2();

    double xbin[14] = {0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};

    total_quarkDeltaR = new TH1D("total_quarkDeltaR", "total_quarkDeltaR", 13, xbin);
    total_quarkDeltaR->Sumw2();

    cut_quarkDeltaRLL = new TH1D("cut_quarkDeltaRLL", "cut_quarkDeltaRLL", 13, xbin);
    cut_quarkDeltaRLL->Sumw2();

    cut_quarkDeltaRMM = new TH1D("cut_quarkDeltaRMM", "cut_quarkDeltaRMM", 13, xbin);
    cut_quarkDeltaRMM->Sumw2();

    cut_quarkDeltaRTT = new TH1D("cut_quarkDeltaRTT", "cut_quarkDeltaRTT", 13, xbin);
    cut_quarkDeltaRTT->Sumw2();

    recoHiggs_MassMM = new TH1D("recoHiggs_MassMM", "recoHiggs_MassMM", 100, 75, 175);
    recoHiggs_MassMM->Sumw2();

    recoHiggs_msoftdropMM = new TH1D("recoHiggs_msoftdropMM", "recoHiggs_msoftdropMM", 100, 75, 175);
    recoHiggs_msoftdropMM->Sumw2();

    recoHiggs_PtMM = new TH1D("recoHiggs_PtMM", "recoHiggs_PtMM", 100, 0, 2000);
    recoHiggs_PtMM->Sumw2();

    recoHiggs_EtaMM = new TH1D("recoHiggs_EtaMM", "recoHiggs_EtaMM", 100, -2.5, 2.5);
    recoHiggs_EtaMM->Sumw2();

    cut_HiggsptMM = new TH1D("cut_HiggsptMM", "cut_HiggsptMM", 30, ybin);
    cut_HiggsptMM->Sumw2();

    cut_Higgspt_2passMM = new TH1D("cut_Higgspt_2passMM", "cut_Higgspt_2passMM", 30, ybin);
    cut_Higgspt_2passMM->Sumw2();

    cut_Higgspt_1passMM = new TH1D("cut_Higgspt_1passMM", "cut_Higgspt_1passMM", 30, ybin);
    cut_Higgspt_1passMM->Sumw2();

    cut_Higgspt_0passMM = new TH1D("cut_Higgspt_0passMM", "cut_Higgspt_0passMM", 30, ybin);
    cut_Higgspt_0passMM->Sumw2();

    recoHiggs_Eta_2passMM = new TH1D("recoHiggs_Eta_2passMM", "recoHiggs_Eta_2passMM", 100, -2.5, 2.5);
    recoHiggs_Eta_2passMM->Sumw2();

    recoHiggs_Eta_1passMM = new TH1D("recoHiggs_Eta_1passMM", "recoHiggs_Eta_1passMM", 100, -2.5, 2.5);
    recoHiggs_Eta_1passMM->Sumw2();

    recoHiggs_Eta_0passMM = new TH1D("recoHiggs_Eta_0passMM", "recoHiggs_Eta_0passMM", 100, -2.5, 2.5);
    recoHiggs_Eta_0passMM->Sumw2();

    recoHiggs_Pt_2passMM = new TH1D("recoHiggs_Pt_2passMM", "recoHiggs_Pt_2passMM", 100, 0, 2000);
    recoHiggs_Pt_2passMM->Sumw2();

    recoHiggs_Pt_1passMM = new TH1D("recoHiggs_Pt_1passMM", "recoHiggs_Pt_1passMM", 100, 0, 2000);
    recoHiggs_Pt_1passMM->Sumw2();

    recoHiggs_Pt_0passMM = new TH1D("recoHiggs_Pt_0passMM", "recoHiggs_Pt_0passMM", 100, 0, 2000);
    recoHiggs_Pt_0passMM->Sumw2();

    recoHiggs_msoftdrop_2passMM = new TH1D("recoHiggs_msoftdrop_2passMM", "recoHiggs_msoftdrop_2passMM", 100, 75, 175);
    recoHiggs_msoftdrop_2passMM->Sumw2();

    recoHiggs_msoftdrop_1passMM = new TH1D("recoHiggs_msoftdrop_1passMM", "recoHiggs_msoftdrop_1passMM", 100, 75, 175);
    recoHiggs_msoftdrop_1passMM->Sumw2();

    recoHiggs_msoftdrop_0passMM = new TH1D("recoHiggs_msoftdrop_0passMM", "recoHiggs_msoftdrop_0passMM", 100, 75, 175);
    recoHiggs_msoftdrop_0passMM->Sumw2();

    recoHiggs_MassTT = new TH1D("recoHiggs_MassTT", "recoHiggs_MassTT", 100, 75, 175);
    recoHiggs_MassTT->Sumw2();

    recoHiggs_msoftdropTT = new TH1D("recoHiggs_msoftdropTT", "recoHiggs_msoftdropTT", 100, 75, 175);
    recoHiggs_msoftdropTT->Sumw2();

    recoHiggs_PtTT = new TH1D("recoHiggs_PtTT", "recoHiggs_PtTT", 100, 0, 2000);
    recoHiggs_PtTT->Sumw2();

    recoHiggs_EtaTT = new TH1D("recoHiggs_EtaTT", "recoHiggs_EtaTT", 100, -2.5, 2.5);
    recoHiggs_EtaTT->Sumw2();

    cut_HiggsptTT = new TH1D("cut_HiggsptTT", "cut_HiggsptTT", 30, ybin);
    cut_HiggsptTT->Sumw2();

    cut_Higgspt_2passTT = new TH1D("cut_Higgspt_2passTT", "cut_Higgspt_2passTT", 30, ybin);
    cut_Higgspt_2passTT->Sumw2();

    cut_Higgspt_1passTT = new TH1D("cut_Higgspt_1passTT", "cut_Higgspt_1passTT", 30, ybin);
    cut_Higgspt_1passTT->Sumw2();

    cut_Higgspt_0passTT = new TH1D("cut_Higgspt_0passTT", "cut_Higgspt_0passTT", 30, ybin);
    cut_Higgspt_0passTT->Sumw2();

    recoHiggs_Eta_2passTT = new TH1D("recoHiggs_Eta_2passTT", "recoHiggs_Eta_2passTT", 100, -2.5, 2.5);
    recoHiggs_Eta_2passTT->Sumw2();

    recoHiggs_Eta_1passTT = new TH1D("recoHiggs_Eta_1passTT", "recoHiggs_Eta_1passTT", 100, -2.5, 2.5);
    recoHiggs_Eta_1passTT->Sumw2();

    recoHiggs_Eta_0passTT = new TH1D("recoHiggs_Eta_0passTT", "recoHiggs_Eta_0passTT", 100, -2.5, 2.5);
    recoHiggs_Eta_0passTT->Sumw2();

    recoHiggs_Pt_2passTT = new TH1D("recoHiggs_Pt_2passTT", "recoHiggs_Pt_2passTT", 100, 0, 2000);
    recoHiggs_Pt_2passTT->Sumw2();

    recoHiggs_Pt_1passTT = new TH1D("recoHiggs_Pt_1passTT", "recoHiggs_Pt_1passTT", 100, 0, 2000);
    recoHiggs_Pt_1passTT->Sumw2();

    recoHiggs_Pt_0passTT = new TH1D("recoHiggs_Pt_0passTT", "recoHiggs_Pt_0passTT", 100, 0, 2000);
    recoHiggs_Pt_0passTT->Sumw2();

    recoHiggs_msoftdrop_2passTT = new TH1D("recoHiggs_msoftdrop_2passTT", "recoHiggs_msoftdrop_2passTT", 100, 75, 175);
    recoHiggs_msoftdrop_2passTT->Sumw2();

    recoHiggs_msoftdrop_1passTT = new TH1D("recoHiggs_msoftdrop_1passTT", "recoHiggs_msoftdrop_1passTT", 100, 75, 175);
    recoHiggs_msoftdrop_1passTT->Sumw2();

    recoHiggs_msoftdrop_0passTT = new TH1D("recoHiggs_msoftdrop_0passTT", "recoHiggs_msoftdrop_0passTT", 100, 75, 175);
    recoHiggs_msoftdrop_0passTT->Sumw2();

    nBquark = new TH1D("nBquark", "nBquark", 5, 0, 5);

    double_Bquark_Pt = new TH1D("double_Bquark_Pt", "double_Bquark_Pt", 100, 0, 2000);
    double_Bquark_Pt->Sumw2();

    Bquark_Pt = new TH1D("Bquark_Pt", "Bquark_Pt", 100, 0, 2000);
    Bquark_Pt->Sumw2();

    anti_Bquark_Pt = new TH1D("anti_Bquark_Pt", "anti_Bquark_Pt", 100, 0, 2000);
    anti_Bquark_Pt->Sumw2();

    leading_Bquark_Pt = new TH1D("leading_Bquark_Pt", "leading_Bquark_Pt", 100, 0, 2000);
    leading_Bquark_Pt->Sumw2();

    sub_Bquark_Pt = new TH1D("sub_Bquark_Pt", "sub_Bquark_Pt", 50, 0, 1000);
    sub_Bquark_Pt->Sumw2();

    leading_Bquark_Eta = new TH1D("leading_Bquark_Eta", "leading_Bquark_Eta", 100, -2.5, 2.5);
    leading_Bquark_Eta->Sumw2();

    sub_Bquark_Eta = new TH1D("sub_Bquark_Eta", "sub_Bquark_Eta", 100, -2.5, 2.5);
    sub_Bquark_Eta->Sumw2();

    Bquark_DeltaR = new TH1D("Bquark_DeltaR", "Bquark_DeltaR", 50, 0, 5);
    Bquark_DeltaR->Sumw2();

    leading_Bquark_Pt_vs_Higgs_Pt = new TH2D("leading_Bquark_Pt_vs_Higgs_Pt", "leading_Bquark_Pt_vs_Higgs_Pt", 100, 0, 2000, 100, 0, 2000);
    leading_Bquark_Pt_vs_Higgs_Pt->Sumw2();

    sub_Bquark_Pt_vs_Higgs_Pt = new TH2D("sub_Bquark_Pt_vs_Higgs_Pt", "sub_Bquark_Pt_vs_Higgs_Pt", 50, 0, 1000, 100, 0, 2000);
    sub_Bquark_Pt_vs_Higgs_Pt->Sumw2();

    double_Bquark_Pt_vs_Higgs_Pt = new TH2D("double_Bquark_Pt_vs_Higgs_Pt", "double_Bquark_Pt_vs_Higgs_Pt", 100, 0, 2000, 100, 0, 2000);
    double_Bquark_Pt_vs_Higgs_Pt->Sumw2();

    Bquark_DeltaR_vs_Higgs_Pt = new TH2D("Bquark_DeltaR_vs_Higgs_Pt", "Bquark_DeltaR_vs_Higgs_Pt", 50, 0, 5, 100, 0, 2000);
    Bquark_DeltaR_vs_Higgs_Pt->Sumw2();
}

void makeHists::saveHists()
{
    hf->cd();
    hf->Write();
    hf->Close();
}
