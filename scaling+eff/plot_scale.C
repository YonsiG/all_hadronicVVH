#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMath.h"
#include "TH1D.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH3F.h"
#include <TApplication.h>
#include <TEnv.h>
#include <TComplex.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>

using namespace std;

int main()
{
    char infileName[100] = "../outfiles/TTJets_SingleLeptFromTbar/TTJets_SingleLeptFromTbar_selected.root";
    char outfileName[100] = "../outfiles/TTJets_SingleLeptFromTbar/TTJets_SingleLeptFromTbar_scaled.root";

    TFile *inputFile = new TFile(infileName);
    TH1D *cutflow = (TH1D *)inputFile->Get("cutflow");
    TH1D *weight_Scale = (TH1D *)inputFile->Get("weight_Scale");

    double scaleNum = 1 / weight_Scale->GetBinContent(1);

    TFile *outputFile = new TFile(outfileName, "RECREATE");
    TH1D *fatjet_btagger;

    fatjet_btagger = (TH1D *)inputFile->Get("fatjet_btagger")->Clone();
    fatjet_btagger->Scale(scaleNum);

    TH1D *HiggsM_QJ_DeltaR0_01LL = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR0_01LL")->Clone();
    HiggsM_QJ_DeltaR0_01LL->Scale(scaleNum);
    TH1D *HiggsM_QJ_DeltaR0_01MM = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR0_01MM")->Clone();
    HiggsM_QJ_DeltaR0_01MM->Scale(scaleNum);
    TH1D *HiggsM_QJ_DeltaR0_01TT = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR0_01TT")->Clone();
    HiggsM_QJ_DeltaR0_01TT->Scale(scaleNum);

    TH1D *HiggsM_QJ_DeltaR01_025LL = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR01_025LL")->Clone();
    HiggsM_QJ_DeltaR01_025LL->Scale(scaleNum);
    TH1D *HiggsM_QJ_DeltaR01_025MM = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR01_025MM")->Clone();
    HiggsM_QJ_DeltaR01_025MM->Scale(scaleNum);
    TH1D *HiggsM_QJ_DeltaR01_025TT = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR01_025TT")->Clone();
    HiggsM_QJ_DeltaR01_025TT->Scale(scaleNum);

    TH1D *HiggsM_QJ_DeltaR025_05LL = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR025_05LL")->Clone();
    HiggsM_QJ_DeltaR025_05LL->Scale(scaleNum);
    TH1D *HiggsM_QJ_DeltaR025_05MM = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR025_05MM")->Clone();
    HiggsM_QJ_DeltaR025_05MM->Scale(scaleNum);
    TH1D *HiggsM_QJ_DeltaR025_05TT = (TH1D *)inputFile->Get("HiggsM_QJ_DeltaR025_05TT")->Clone();
    HiggsM_QJ_DeltaR025_05TT->Scale(scaleNum);

    TH2D *QJ_DeltaR_vs_HiggsMLL = (TH2D *)inputFile->Get("QJ_DeltaR_vs_HiggsMLL")->Clone();
    QJ_DeltaR_vs_HiggsMLL->Scale(scaleNum);
    TH2D *QJ_DeltaR_vs_HiggsMMM = (TH2D *)inputFile->Get("QJ_DeltaR_vs_HiggsMMM")->Clone();
    QJ_DeltaR_vs_HiggsMMM->Scale(scaleNum);
    TH2D *QJ_DeltaR_vs_HiggsMTT = (TH2D *)inputFile->Get("QJ_DeltaR_vs_HiggsMTT")->Clone();
    QJ_DeltaR_vs_HiggsMTT->Scale(scaleNum);

    TH2D *QJ_DeltaR1_vs_HiggsMLL = (TH2D *)inputFile->Get("QJ_DeltaR1_vs_HiggsMLL")->Clone();
    QJ_DeltaR1_vs_HiggsMLL->Scale(scaleNum);
    TH2D *QJ_DeltaR1_vs_HiggsMMM = (TH2D *)inputFile->Get("QJ_DeltaR1_vs_HiggsMMM")->Clone();
    QJ_DeltaR1_vs_HiggsMMM->Scale(scaleNum);
    TH2D *QJ_DeltaR1_vs_HiggsMTT = (TH2D *)inputFile->Get("QJ_DeltaR1_vs_HiggsMTT")->Clone();
    QJ_DeltaR1_vs_HiggsMTT->Scale(scaleNum);

    TH2D *QJ_DeltaR2_vs_HiggsMLL = (TH2D *)inputFile->Get("QJ_DeltaR2_vs_HiggsMLL")->Clone();
    QJ_DeltaR2_vs_HiggsMLL->Scale(scaleNum);
    TH2D *QJ_DeltaR2_vs_HiggsMMM = (TH2D *)inputFile->Get("QJ_DeltaR2_vs_HiggsMMM")->Clone();
    QJ_DeltaR2_vs_HiggsMMM->Scale(scaleNum);
    TH2D *QJ_DeltaR2_vs_HiggsMTT = (TH2D *)inputFile->Get("QJ_DeltaR2_vs_HiggsMTT")->Clone();
    QJ_DeltaR2_vs_HiggsMTT->Scale(scaleNum);

    TH2D *QJ_DeltaR1_vs_QJ_DeltaR2LL = (TH2D *)inputFile->Get("QJ_DeltaR1_vs_QJ_DeltaR2LL")->Clone();
    QJ_DeltaR1_vs_QJ_DeltaR2LL->Scale(scaleNum);
    TH2D *QJ_DeltaR1_vs_QJ_DeltaR2MM = (TH2D *)inputFile->Get("QJ_DeltaR1_vs_QJ_DeltaR2MM")->Clone();
    QJ_DeltaR1_vs_QJ_DeltaR2MM->Scale(scaleNum);
    TH2D *QJ_DeltaR1_vs_QJ_DeltaR2TT = (TH2D *)inputFile->Get("QJ_DeltaR1_vs_QJ_DeltaR2TT")->Clone();
    QJ_DeltaR1_vs_QJ_DeltaR2TT->Scale(scaleNum);

    TH1D *QJ_DeltaRLL = (TH1D *)inputFile->Get("QJ_DeltaRLL")->Clone();
    QJ_DeltaRLL->Scale(scaleNum);
    TH1D *QJ_DeltaRMM = (TH1D *)inputFile->Get("QJ_DeltaRMM")->Clone();
    QJ_DeltaRMM->Scale(scaleNum);
    TH1D *QJ_DeltaRTT = (TH1D *)inputFile->Get("QJ_DeltaRTT")->Clone();
    QJ_DeltaRTT->Scale(scaleNum);

    TH1D *QJ_DeltaR1LL = (TH1D *)inputFile->Get("QJ_DeltaR1LL")->Clone();
    QJ_DeltaR1LL->Scale(scaleNum);
    TH1D *QJ_DeltaR1MM = (TH1D *)inputFile->Get("QJ_DeltaR1MM")->Clone();
    QJ_DeltaR1MM->Scale(scaleNum);
    TH1D *QJ_DeltaR1TT = (TH1D *)inputFile->Get("QJ_DeltaR1TT")->Clone();
    QJ_DeltaR1TT->Scale(scaleNum);

    TH1D *QJ_DeltaR2LL = (TH1D *)inputFile->Get("QJ_DeltaR2LL")->Clone();
    QJ_DeltaR2LL->Scale(scaleNum);
    TH1D *QJ_DeltaR2MM = (TH1D *)inputFile->Get("QJ_DeltaR2MM")->Clone();
    QJ_DeltaR2MM->Scale(scaleNum);
    TH1D *QJ_DeltaR2TT = (TH1D *)inputFile->Get("QJ_DeltaR2TT")->Clone();
    QJ_DeltaR2TT->Scale(scaleNum);

    TH1D *total_Higgspt = (TH1D *)inputFile->Get("total_Higgspt")->Clone();
    total_Higgspt->Scale(scaleNum);
    TH1D *cut_HiggsptLL = (TH1D *)inputFile->Get("cut_HiggsptLL")->Clone();
    cut_HiggsptLL->Scale(scaleNum);
    TH1D *cut_HiggsptMM = (TH1D *)inputFile->Get("cut_HiggsptMM")->Clone();
    cut_HiggsptMM->Scale(scaleNum);
    TH1D *cut_HiggsptTT = (TH1D *)inputFile->Get("cut_HiggsptTT")->Clone();
    cut_HiggsptTT->Scale(scaleNum);
    TH1D *total_quarkDeltaR = (TH1D *)inputFile->Get("total_quarkDeltaR")->Clone();
    total_quarkDeltaR->Scale(scaleNum);
    TH1D *cut_quarkDeltaRLL = (TH1D *)inputFile->Get("cut_quarkDeltaRLL")->Clone();
    cut_quarkDeltaRLL->Scale(scaleNum);
    TH1D *cut_quarkDeltaRMM = (TH1D *)inputFile->Get("cut_quarkDeltaRMM")->Clone();
    cut_quarkDeltaRMM->Scale(scaleNum);
    TH1D *cut_quarkDeltaRTT = (TH1D *)inputFile->Get("cut_quarkDeltaRTT")->Clone();
    cut_quarkDeltaRTT->Scale(scaleNum);
    TH1D *cut_Higgspt_2passLL = (TH1D *)inputFile->Get("cut_Higgspt_2passLL")->Clone();
    cut_Higgspt_2passLL->Scale(scaleNum);
    TH1D *cut_Higgspt_2passMM = (TH1D *)inputFile->Get("cut_Higgspt_2passMM")->Clone();
    cut_Higgspt_2passMM->Scale(scaleNum);
    TH1D *cut_Higgspt_2passTT = (TH1D *)inputFile->Get("cut_Higgspt_2passTT")->Clone();
    cut_Higgspt_2passTT->Scale(scaleNum);
    TH1D *cut_Higgspt_1passLL = (TH1D *)inputFile->Get("cut_Higgspt_1passLL")->Clone();
    cut_Higgspt_1passLL->Scale(scaleNum);
    TH1D *cut_Higgspt_1passMM = (TH1D *)inputFile->Get("cut_Higgspt_1passMM")->Clone();
    cut_Higgspt_1passMM->Scale(scaleNum);
    TH1D *cut_Higgspt_1passTT = (TH1D *)inputFile->Get("cut_Higgspt_1passTT")->Clone();
    cut_Higgspt_1passTT->Scale(scaleNum);
    TH1D *cut_Higgspt_0passLL = (TH1D *)inputFile->Get("cut_Higgspt_0passLL")->Clone();
    cut_Higgspt_0passLL->Scale(scaleNum);
    TH1D *cut_Higgspt_0passMM = (TH1D *)inputFile->Get("cut_Higgspt_0passMM")->Clone();
    cut_Higgspt_0passMM->Scale(scaleNum);
    TH1D *cut_Higgspt_0passTT = (TH1D *)inputFile->Get("cut_Higgspt_0passTT")->Clone();
    cut_Higgspt_0passTT->Scale(scaleNum);

    TH1D *Higgs_msoftdrop_afterpurityLL = (TH1D *)inputFile->Get("Higgs_msoftdrop_afterpurityLL")->Clone();
    Higgs_msoftdrop_afterpurityLL->Scale(scaleNum);
    TH1D *recoHiggs_MassLL = (TH1D *)inputFile->Get("recoHiggs_MassLL")->Clone();
    recoHiggs_MassLL->Scale(scaleNum);
    TH1D *recoHiggs_msoftdropLL = (TH1D *)inputFile->Get("recoHiggs_msoftdropLL")->Clone();
    recoHiggs_msoftdropLL->Scale(scaleNum);
    TH1D *recoHiggs_PtLL = (TH1D *)inputFile->Get("recoHiggs_PtLL")->Clone();
    recoHiggs_PtLL->Scale(scaleNum);
    TH1D *recoHiggs_EtaLL = (TH1D *)inputFile->Get("recoHiggs_EtaLL")->Clone();
    recoHiggs_EtaLL->Scale(scaleNum);
    TH1D *recoHiggs_Eta_2passLL = (TH1D *)inputFile->Get("recoHiggs_Eta_2passLL")->Clone();
    recoHiggs_Eta_2passLL->Scale(scaleNum);
    TH1D *recoHiggs_Eta_1passLL = (TH1D *)inputFile->Get("recoHiggs_Eta_1passLL")->Clone();
    recoHiggs_Eta_1passLL->Scale(scaleNum);
    TH1D *recoHiggs_Eta_0passLL = (TH1D *)inputFile->Get("recoHiggs_Eta_0passLL")->Clone();
    recoHiggs_Eta_0passLL->Scale(scaleNum);
    TH1D *recoHiggs_Pt_2passLL = (TH1D *)inputFile->Get("recoHiggs_Pt_2passLL")->Clone();
    recoHiggs_Pt_2passLL->Scale(scaleNum);
    TH1D *recoHiggs_Pt_1passLL = (TH1D *)inputFile->Get("recoHiggs_Pt_1passLL")->Clone();
    recoHiggs_Pt_1passLL->Scale(scaleNum);
    TH1D *recoHiggs_Pt_0passLL = (TH1D *)inputFile->Get("recoHiggs_Pt_0passLL")->Clone();
    recoHiggs_Pt_0passLL->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_2passLL = (TH1D *)inputFile->Get("Higgs_msoftdrop_2passLL")->Clone();
    Higgs_msoftdrop_2passLL->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_1passLL = (TH1D *)inputFile->Get("Higgs_msoftdrop_1passLL")->Clone();
    Higgs_msoftdrop_1passLL->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_0passLL = (TH1D *)inputFile->Get("Higgs_msoftdrop_0passLL")->Clone();
    Higgs_msoftdrop_0passLL->Scale(scaleNum);

    TH1D *Higgs_msoftdrop_afterpurityMM = (TH1D *)inputFile->Get("Higgs_msoftdrop_afterpurityMM")->Clone();
    Higgs_msoftdrop_afterpurityMM->Scale(scaleNum);
    TH1D *recoHiggs_MassMM = (TH1D *)inputFile->Get("recoHiggs_MassMM")->Clone();
    recoHiggs_MassMM->Scale(scaleNum);
    TH1D *recoHiggs_msoftdropMM = (TH1D *)inputFile->Get("recoHiggs_msoftdropMM")->Clone();
    recoHiggs_msoftdropMM->Scale(scaleNum);
    TH1D *recoHiggs_PtMM = (TH1D *)inputFile->Get("recoHiggs_PtMM")->Clone();
    recoHiggs_PtMM->Scale(scaleNum);
    TH1D *recoHiggs_EtaMM = (TH1D *)inputFile->Get("recoHiggs_EtaMM")->Clone();
    recoHiggs_EtaMM->Scale(scaleNum);
    TH1D *recoHiggs_Eta_2passMM = (TH1D *)inputFile->Get("recoHiggs_Eta_2passMM")->Clone();
    recoHiggs_Eta_2passMM->Scale(scaleNum);
    TH1D *recoHiggs_Eta_1passMM = (TH1D *)inputFile->Get("recoHiggs_Eta_1passMM")->Clone();
    recoHiggs_Eta_1passMM->Scale(scaleNum);
    TH1D *recoHiggs_Eta_0passMM = (TH1D *)inputFile->Get("recoHiggs_Eta_0passMM")->Clone();
    recoHiggs_Eta_0passMM->Scale(scaleNum);
    TH1D *recoHiggs_Pt_2passMM = (TH1D *)inputFile->Get("recoHiggs_Pt_2passMM")->Clone();
    recoHiggs_Pt_2passMM->Scale(scaleNum);
    TH1D *recoHiggs_Pt_1passMM = (TH1D *)inputFile->Get("recoHiggs_Pt_1passMM")->Clone();
    recoHiggs_Pt_1passMM->Scale(scaleNum);
    TH1D *recoHiggs_Pt_0passMM = (TH1D *)inputFile->Get("recoHiggs_Pt_0passMM")->Clone();
    recoHiggs_Pt_0passMM->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_2passMM = (TH1D *)inputFile->Get("Higgs_msoftdrop_2passMM")->Clone();
    Higgs_msoftdrop_2passMM->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_1passMM = (TH1D *)inputFile->Get("Higgs_msoftdrop_1passMM")->Clone();
    Higgs_msoftdrop_1passMM->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_0passMM = (TH1D *)inputFile->Get("Higgs_msoftdrop_0passMM")->Clone();
    Higgs_msoftdrop_0passMM->Scale(scaleNum);

    TH1D *Higgs_msoftdrop_afterpurityTT = (TH1D *)inputFile->Get("Higgs_msoftdrop_afterpurityTT")->Clone();;
    Higgs_msoftdrop_afterpurityTT->Scale(scaleNum);
    TH1D *recoHiggs_MassTT = (TH1D *)inputFile->Get("recoHiggs_MassTT")->Clone();;
    recoHiggs_MassTT->Scale(scaleNum);
    TH1D *recoHiggs_msoftdropTT = (TH1D *)inputFile->Get("recoHiggs_msoftdropTT")->Clone();;
    recoHiggs_msoftdropTT->Scale(scaleNum);
    TH1D *recoHiggs_PtTT = (TH1D *)inputFile->Get("recoHiggs_PtTT")->Clone();;
    recoHiggs_PtTT->Scale(scaleNum);
    TH1D *recoHiggs_EtaTT = (TH1D *)inputFile->Get("recoHiggs_EtaTT")->Clone();;
    recoHiggs_EtaTT->Scale(scaleNum);
    TH1D *recoHiggs_Eta_2passTT = (TH1D *)inputFile->Get("recoHiggs_Eta_2passTT")->Clone();;
    recoHiggs_Eta_2passTT->Scale(scaleNum);
    TH1D *recoHiggs_Eta_1passTT = (TH1D *)inputFile->Get("recoHiggs_Eta_1passTT")->Clone();;
    recoHiggs_Eta_1passTT->Scale(scaleNum);
    TH1D *recoHiggs_Eta_0passTT = (TH1D *)inputFile->Get("recoHiggs_Eta_0passTT")->Clone();;
    recoHiggs_Eta_0passTT->Scale(scaleNum);
    TH1D *recoHiggs_Pt_2passTT = (TH1D *)inputFile->Get("recoHiggs_Pt_2passTT")->Clone();;
    recoHiggs_Pt_2passTT->Scale(scaleNum);
    TH1D *recoHiggs_Pt_1passTT = (TH1D *)inputFile->Get("recoHiggs_Pt_1passTT")->Clone();;
    recoHiggs_Pt_1passTT->Scale(scaleNum);
    TH1D *recoHiggs_Pt_0passTT = (TH1D *)inputFile->Get("recoHiggs_Pt_0passTT")->Clone();;
    recoHiggs_Pt_0passTT->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_2passTT = (TH1D *)inputFile->Get("Higgs_msoftdrop_2passTT")->Clone();;
    Higgs_msoftdrop_2passTT->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_1passTT = (TH1D *)inputFile->Get("Higgs_msoftdrop_1passTT")->Clone();;
    Higgs_msoftdrop_1passTT->Scale(scaleNum);
    TH1D *Higgs_msoftdrop_0passTT = (TH1D *)inputFile->Get("Higgs_msoftdrop_0passTT")->Clone();;
    Higgs_msoftdrop_0passTT->Scale(scaleNum);

    TH1D *nBquark = (TH1D *)inputFile->Get("nBquark")->Clone();;
    nBquark->Scale(scaleNum);
    TH1D *double_Bquark_Pt = (TH1D *)inputFile->Get("double_Bquark_Pt")->Clone();;
    double_Bquark_Pt->Scale(scaleNum);
    TH1D *Bquark_Pt = (TH1D *)inputFile->Get("Bquark_Pt")->Clone();;
    Bquark_Pt->Scale(scaleNum);
    TH1D *anti_Bquark_Pt = (TH1D *)inputFile->Get("anti_Bquark_Pt")->Clone();;
    anti_Bquark_Pt->Scale(scaleNum);
    TH1D *leading_Bquark_Pt = (TH1D *)inputFile->Get("leading_Bquark_Pt")->Clone();;
    leading_Bquark_Pt->Scale(scaleNum);
    TH1D *sub_Bquark_Pt = (TH1D *)inputFile->Get("sub_Bquark_Pt")->Clone();;
    sub_Bquark_Pt->Scale(scaleNum);
    TH1D *leading_Bquark_Eta = (TH1D *)inputFile->Get("leading_Bquark_Eta")->Clone();;
    leading_Bquark_Eta->Scale(scaleNum);
    TH1D *sub_Bquark_Eta = (TH1D *)inputFile->Get("sub_Bquark_Eta")->Clone();;
    sub_Bquark_Eta->Scale(scaleNum);
    TH1D *Bquark_DeltaR = (TH1D *)inputFile->Get("Bquark_DeltaR")->Clone();;
    Bquark_DeltaR->Scale(scaleNum);
    TH2D *double_Bquark_Pt_vs_Higgs_Pt = (TH2D *)inputFile->Get("double_Bquark_Pt_vs_Higgs_Pt")->Clone();;
    double_Bquark_Pt_vs_Higgs_Pt->Scale(scaleNum);
    TH2D *leading_Bquark_Pt_vs_Higgs_Pt = (TH2D *)inputFile->Get("leading_Bquark_Pt_vs_Higgs_Pt")->Clone();;
    leading_Bquark_Pt_vs_Higgs_Pt->Scale(scaleNum);
    TH2D *sub_Bquark_Pt_vs_Higgs_Pt = (TH2D *)inputFile->Get("sub_Bquark_Pt_vs_Higgs_Pt")->Clone();;
    sub_Bquark_Pt_vs_Higgs_Pt->Scale(scaleNum);
    TH2D *Bquark_DeltaR_vs_Higgs_Pt = (TH2D *)inputFile->Get("Bquark_DeltaR_vs_Higgs_Pt")->Clone();;
    Bquark_DeltaR_vs_Higgs_Pt->Scale(scaleNum);

    inputFile->Close();
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();
}
