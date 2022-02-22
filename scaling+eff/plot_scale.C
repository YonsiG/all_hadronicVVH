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

void scale(TString fileName)
{
    TString infileName_s = "../../outfiles/" + fileName + "_selected.root";
    TString outfileName_s = "../../outfiles/" + fileName + "_scaled.root";
    TFile *inputFile = new TFile(infileName_s);
    TH1D *weight_Scale = (TH1D *)inputFile->Get("weight_Scale");
    double scaleNum = 1 / weight_Scale->GetBinContent(1);

    TFile *outputFile = new TFile(outfileName_s, "RECREATE");
    char plotname[50];
    TH1D *cutflow[14];
    TH1D *number_of_fatjets[14];
    TH1D *number_of_jets[14];
    TH1D *number_of_central_jets[14];
    TH1D *fatjet_btag_score;
    TH1D *first_fatjet_btag_score;
    TH1D *second_fatjet_btag_score;
    TH1D *third_fatjet_btag_score;
    TH1D *first_fatjet_btag_score_2match;
    TH1D *second_fatjet_btag_score_2match;
    TH1D *third_fatjet_btag_score_2match;
    TH1D *first_fatjet_btag_score_1match;
    TH1D *second_fatjet_btag_score_1match;
    TH1D *third_fatjet_btag_score_1match;
    TH1D *first_fatjet_btag_score_0match;
    TH1D *second_fatjet_btag_score_0match;
    TH1D *third_fatjet_btag_score_0match;
    TH1D *VBFJet_leadingPt_2match;
    TH1D *VBFJet_leadingPt_1match;
    TH1D *VBFJet_leadingPt_0match;
    TH1D *VBFJet_subleadingPt_2match;
    TH1D *VBFJet_subleadingPt_1match;
    TH1D *VBFJet_subleadingPt_0match;
    TH1D *VBFJet_Mjj_2match;
    TH1D *VBFJet_Mjj_1match;
    TH1D *VBFJet_Mjj_0match;
    TH1D *VBFJet_DeltaEta_2match_total;
    TH1D *VBFJet_DeltaEta_1match_total;
    TH1D *VBFJet_DeltaEta_0match_total;
    TH1D *VBFJet_DeltaEta_2match[14];
    TH1D *VBFJet_DeltaEta_1match[14];
    TH1D *VBFJet_DeltaEta_0match[14];
    TH1D *matched_VBFJet_Pt_total;
    TH1D *matched_VBFJet_Eta_total;
    TH1D *matched_VBFJet_qgl_total;
    TH1D *matched_VBFJet_Pt[14];
    TH1D *matched_VBFJet_Eta[14];
    TH1D *matched_VBFJet_qgl[14];
    TH1D *fatjet_msoftdrop[14][3];
    TH1D *fatjet_pt[14][3];
    TH1D *fatjet_eta[14][3];
    TH1D *fatjet_WvsQCD[14][3];
    TH1D *fatjet_mass[14][3];
    TH1D *fatjet_Xbb_modified[14][3];
    TH1D *fatjet_Xcc[14][3];
    TH1D *fatjet_Xqq[14][3];
    TH1D *fatjet_QCD[14][3];
    TH1D *VBF_max_mass[14];

    for (int icategory = 0; icategory < 14; icategory++)
    {
        sprintf(plotname, "cutflow%i", icategory);
        cutflow[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "number_of_fatjets%i", icategory);
        number_of_fatjets[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "number_of_jets%i", icategory);
        number_of_jets[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "number_of_central_jets%i", icategory);
        number_of_central_jets[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "VBFJet_DeltaEta_2match%i", icategory);
        VBFJet_DeltaEta_2match[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "VBFJet_DeltaEta_1match%i", icategory);
        VBFJet_DeltaEta_1match[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "VBFJet_DeltaEta_0match%i", icategory);
        VBFJet_DeltaEta_0match[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "matched_VBFJet_Pt%i", icategory);
        matched_VBFJet_Pt[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "matched_VBFJet_Eta%i", icategory);
        matched_VBFJet_Eta[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "matched_VBFJet_qgl%i", icategory);
        matched_VBFJet_qgl[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "VBF_max_mass%i", icategory);
        VBF_max_mass[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();

        cutflow[icategory]->Scale(scaleNum);
        number_of_fatjets[icategory]->Scale(scaleNum);
        number_of_jets[icategory]->Scale(scaleNum);
        number_of_central_jets[icategory]->Scale(scaleNum);
        VBFJet_DeltaEta_2match[icategory]->Scale(scaleNum);
        VBFJet_DeltaEta_1match[icategory]->Scale(scaleNum);
        VBFJet_DeltaEta_0match[icategory]->Scale(scaleNum);
        matched_VBFJet_Pt[icategory]->Scale(scaleNum);
        matched_VBFJet_Eta[icategory]->Scale(scaleNum);
        matched_VBFJet_qgl[icategory]->Scale(scaleNum);
        VBF_max_mass[icategory]->Scale(scaleNum);

        for (int ifatjet=0; ifatjet<3; ifatjet++){
            sprintf(plotname, "fatjet_msoftdrop%i_%i", icategory, ifatjet);
            fatjet_msoftdrop[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
            sprintf(plotname, "fatjet_pt%i_%i", icategory, ifatjet);
            fatjet_pt[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
            sprintf(plotname, "fatjet_eta%i_%i", icategory, ifatjet);
            fatjet_eta[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
            sprintf(plotname, "fatjet_WvsQCD%i_%i", icategory, ifatjet);
            fatjet_WvsQCD[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
//            sprintf(plotname, "fatjet_mass%i_%i", icategory, ifatjet);
 //           fatjet_mass[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
            sprintf(plotname, "fatjet_Xbb_modified%i_%i", icategory, ifatjet);
            fatjet_Xbb_modified[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
            sprintf(plotname, "fatjet_Xcc%i_%i", icategory, ifatjet);
            fatjet_Xcc[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
            sprintf(plotname, "fatjet_Xqq%i_%i", icategory, ifatjet);
            fatjet_Xqq[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
            sprintf(plotname, "fatjet_QCD%i_%i", icategory, ifatjet);
            fatjet_QCD[icategory][ifatjet] = (TH1D *)inputFile->Get(plotname)->Clone();
            fatjet_msoftdrop[icategory][ifatjet]->Scale(scaleNum);
            fatjet_pt[icategory][ifatjet]->Scale(scaleNum);
            fatjet_eta[icategory][ifatjet]->Scale(scaleNum);
            fatjet_WvsQCD[icategory][ifatjet]->Scale(scaleNum);
//            fatjet_mass[icategory][ifatjet]->Scale(scaleNum);
            fatjet_Xbb_modified[icategory][ifatjet]->Scale(scaleNum);
            fatjet_Xcc[icategory][ifatjet]->Scale(scaleNum);
            fatjet_Xqq[icategory][ifatjet]->Scale(scaleNum);
            fatjet_QCD[icategory][ifatjet]->Scale(scaleNum);
        }
    }

    fatjet_btag_score = (TH1D *)inputFile->Get("fatjet_btag_score")->Clone();
    fatjet_btag_score->Scale(scaleNum);
    first_fatjet_btag_score = (TH1D *)inputFile->Get("first_fatjet_btag_score")->Clone();
    first_fatjet_btag_score->Scale(scaleNum);
    second_fatjet_btag_score = (TH1D *)inputFile->Get("second_fatjet_btag_score")->Clone();
    second_fatjet_btag_score->Scale(scaleNum);
    third_fatjet_btag_score = (TH1D *)inputFile->Get("third_fatjet_btag_score")->Clone();
    third_fatjet_btag_score->Scale(scaleNum);

    first_fatjet_btag_score_2match = (TH1D *)inputFile->Get("first_fatjet_btag_score_2match")->Clone();
    first_fatjet_btag_score_2match->Scale(scaleNum);
    second_fatjet_btag_score_2match = (TH1D *)inputFile->Get("second_fatjet_btag_score_2match")->Clone();
    second_fatjet_btag_score_2match->Scale(scaleNum);
    third_fatjet_btag_score_2match = (TH1D *)inputFile->Get("third_fatjet_btag_score_2match")->Clone();
    third_fatjet_btag_score_2match->Scale(scaleNum);

    first_fatjet_btag_score_1match = (TH1D *)inputFile->Get("first_fatjet_btag_score_1match")->Clone();
    first_fatjet_btag_score_1match->Scale(scaleNum);
    second_fatjet_btag_score_1match = (TH1D *)inputFile->Get("second_fatjet_btag_score_1match")->Clone();
    second_fatjet_btag_score_1match->Scale(scaleNum);
    third_fatjet_btag_score_1match = (TH1D *)inputFile->Get("third_fatjet_btag_score_1match")->Clone();
    third_fatjet_btag_score_1match->Scale(scaleNum);

    first_fatjet_btag_score_0match = (TH1D *)inputFile->Get("first_fatjet_btag_score_0match")->Clone();
    first_fatjet_btag_score_0match->Scale(scaleNum);
    second_fatjet_btag_score_0match = (TH1D *)inputFile->Get("second_fatjet_btag_score_0match")->Clone();
    second_fatjet_btag_score_0match->Scale(scaleNum);
    third_fatjet_btag_score_0match = (TH1D *)inputFile->Get("third_fatjet_btag_score_0match")->Clone();
    third_fatjet_btag_score_0match->Scale(scaleNum);

    VBFJet_leadingPt_2match = (TH1D *)inputFile->Get("VBFJet_leadingPt_2match")->Clone();
    VBFJet_leadingPt_2match->Scale(scaleNum);
    VBFJet_leadingPt_1match = (TH1D *)inputFile->Get("VBFJet_leadingPt_1match")->Clone();
    VBFJet_leadingPt_1match->Scale(scaleNum);
    VBFJet_leadingPt_0match = (TH1D *)inputFile->Get("VBFJet_leadingPt_0match")->Clone();
    VBFJet_leadingPt_0match->Scale(scaleNum);
    VBFJet_subleadingPt_2match = (TH1D *)inputFile->Get("VBFJet_subleadingPt_2match")->Clone();
    VBFJet_subleadingPt_2match->Scale(scaleNum);
    VBFJet_subleadingPt_1match = (TH1D *)inputFile->Get("VBFJet_subleadingPt_1match")->Clone();
    VBFJet_subleadingPt_1match->Scale(scaleNum);
    VBFJet_subleadingPt_0match = (TH1D *)inputFile->Get("VBFJet_subleadingPt_0match")->Clone();
    VBFJet_subleadingPt_0match->Scale(scaleNum);

    VBFJet_Mjj_2match = (TH1D *)inputFile->Get("VBFJet_Mjj_2match")->Clone();
    VBFJet_Mjj_2match->Scale(scaleNum);
    VBFJet_Mjj_1match = (TH1D *)inputFile->Get("VBFJet_Mjj_1match")->Clone();
    VBFJet_Mjj_1match->Scale(scaleNum);
    VBFJet_Mjj_0match = (TH1D *)inputFile->Get("VBFJet_Mjj_0match")->Clone();
    VBFJet_Mjj_0match->Scale(scaleNum);

    VBFJet_DeltaEta_2match_total = (TH1D *)inputFile->Get("VBFJet_DeltaEta_2match_total")->Clone();
    VBFJet_DeltaEta_2match_total->Scale(scaleNum);
    VBFJet_DeltaEta_1match_total = (TH1D *)inputFile->Get("VBFJet_DeltaEta_1match_total")->Clone();
    VBFJet_DeltaEta_1match_total->Scale(scaleNum);
    VBFJet_DeltaEta_0match_total = (TH1D *)inputFile->Get("VBFJet_DeltaEta_0match_total")->Clone();
    VBFJet_DeltaEta_0match_total->Scale(scaleNum);

    matched_VBFJet_Pt_total = (TH1D *)inputFile->Get("matched_VBFJet_Pt_total")->Clone();
    matched_VBFJet_Pt_total->Scale(scaleNum);
    matched_VBFJet_Eta_total = (TH1D *)inputFile->Get("matched_VBFJet_Eta_total")->Clone();
    matched_VBFJet_Eta_total->Scale(scaleNum);
    matched_VBFJet_qgl_total = (TH1D *)inputFile->Get("matched_VBFJet_qgl_total")->Clone();
    matched_VBFJet_qgl_total->Scale(scaleNum);

    inputFile->Close();
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();
}

int main()
{
    TString signalNames[] = {"SSWWH_2", "OSWWH_2", "WZH_2", "ZZH_2"};
    for (int isize=0; isize<sizeof(signalNames)/sizeof(signalNames[0]); isize++)
    {
        scale(signalNames[isize]);
    }

    TString QCDNames[] = {"QCD_HT300to500_2", "QCD_HT500to700_2", "QCD_HT700to1000_2", "QCD_HT1000to1500_2", "QCD_HT1500to2000_2", "QCD_HT2000toInf_2"};
    for (int isize=0; isize<sizeof(QCDNames)/sizeof(QCDNames[0]); isize++)
    {
        scale(QCDNames[isize]);
    }

    TString ZJetNames[] = {"ZJetsToQQ_HT-200to400_2", "ZJetsToQQ_HT-400to600_2", "ZJetsToQQ_HT-600to800_2", "ZJetsToQQ_HT-800toInf_2"};
    for (int isize=0; isize<sizeof(ZJetNames)/sizeof(ZJetNames[0]); isize++)
    {
        scale(ZJetNames[isize]);
    }

    TString WJetNames[] = {"WJetsToQQ_HT-200to400_2", "WJetsToQQ_HT-400to600_2", "WJetsToQQ_HT-600to800_2", "WJetsToQQ_HT-800toInf_2"};
    for (int isize=0; isize<sizeof(WJetNames)/sizeof(WJetNames[0]); isize++)
    {
        scale(WJetNames[isize]);
    }

    TString TTToHadronicNames[] = {"TTToHadronic_2"};
    for (int isize=0; isize<sizeof(TTToHadronicNames)/sizeof(TTToHadronicNames[0]); isize++)
    {
        scale(TTToHadronicNames[isize]);
    }

    TString TTToSemiLeptonicNames[] = {"TTToSemiLeptonic_2"};
    for (int isize=0; isize<sizeof(TTToSemiLeptonicNames)/sizeof(TTToSemiLeptonicNames[0]); isize++)
    {
        scale(TTToSemiLeptonicNames[isize]);
    }
    
}