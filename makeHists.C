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

    char plotname[50];
    for (int icategory = 0; icategory < 14; icategory++)
    {
        sprintf(plotname, "cutflow%i", icategory);
        cutflow[icategory] = new TH1D(plotname, plotname, 20, 0, 20);

        sprintf(plotname, "number_of_fatjets%i", icategory);
        number_of_fatjets[icategory] = new TH1D(plotname, plotname, 10, 0, 10);

        sprintf(plotname, "number_of_jets%i", icategory);
        number_of_jets[icategory] = new TH1D(plotname, plotname, 10, 0, 10);

        sprintf(plotname, "number_of_central_jets%i", icategory);
        number_of_central_jets[icategory] = new TH1D(plotname, plotname, 10, 0, 10);

        sprintf(plotname, "VBFJet_DeltaEta_2match%i", icategory);
        VBFJet_DeltaEta_2match[icategory] = new TH1D(plotname, plotname, 100, 0, 10);
        VBFJet_DeltaEta_2match[icategory]->Sumw2();

        sprintf(plotname, "VBFJet_DeltaEta_1match%i", icategory);
        VBFJet_DeltaEta_1match[icategory] = new TH1D(plotname, plotname, 100, 0, 10);
        VBFJet_DeltaEta_1match[icategory]->Sumw2();

        sprintf(plotname, "VBFJet_DeltaEta_0match%i", icategory);
        VBFJet_DeltaEta_0match[icategory] = new TH1D(plotname, plotname, 100, 0, 10);
        VBFJet_DeltaEta_0match[icategory]->Sumw2();

        sprintf(plotname, "matched_VBFJet_Pt%i", icategory);
        matched_VBFJet_Pt[icategory] = new TH1D(plotname, plotname, 100, 0, 1000);
        matched_VBFJet_Pt[icategory]->Sumw2();

        sprintf(plotname, "matched_VBFJet_Eta%i", icategory);
        matched_VBFJet_Eta[icategory] = new TH1D(plotname, plotname, 100, -5, 5);
        matched_VBFJet_Eta[icategory]->Sumw2();

        sprintf(plotname, "matched_VBFJet_qgl%i", icategory);
        matched_VBFJet_qgl[icategory] = new TH1D(plotname, plotname, 100, 0, 1);
        matched_VBFJet_qgl[icategory]->Sumw2();

        for (int ifatjet=0; ifatjet<3; ifatjet++){
            sprintf(plotname, "fatjet_msoftdrop%i_%i", icategory, ifatjet);
            fatjet_msoftdrop[icategory][ifatjet] = new TH1D(plotname,plotname,100,40,1000);
            fatjet_msoftdrop[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_pt%i_%i", icategory, ifatjet);
            fatjet_pt[icategory][ifatjet] = new TH1D(plotname,plotname,100,250,3000);
            fatjet_pt[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_eta%i_%i", icategory, ifatjet);
            fatjet_eta[icategory][ifatjet] = new TH1D(plotname,plotname,100,-2.5,2.5);
            fatjet_eta[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_WvsQCD%i_%i", icategory, ifatjet);
            fatjet_WvsQCD[icategory][ifatjet] = new TH1D(plotname,plotname,100,0,1);
            fatjet_WvsQCD[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_mass%i_%i", icategory, ifatjet);
            fatjet_mass[icategory][ifatjet] = new TH1D(plotname,plotname,100,40,1000);
            fatjet_mass[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_Xbb_modified%i_%i", icategory, ifatjet);
            fatjet_Xbb_modified[icategory][ifatjet] = new TH1D(plotname,plotname,100,0,1);
            fatjet_Xbb_modified[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_Xcc%i_%i", icategory, ifatjet);
            fatjet_Xcc[icategory][ifatjet] = new TH1D(plotname,plotname,100,0,1);
            fatjet_Xcc[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_Xqq%i_%i", icategory, ifatjet);
            fatjet_Xqq[icategory][ifatjet] = new TH1D(plotname,plotname,100,0,1);
            fatjet_Xqq[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_QCD%i_%i", icategory, ifatjet);
            fatjet_QCD[icategory][ifatjet] = new TH1D(plotname,plotname,100,0,1);
            fatjet_QCD[icategory][ifatjet]->Sumw2();

            sprintf(plotname, "fatjet_Xccqq_modified%i_%i", icategory, ifatjet);
            fatjet_Xccqq_modified[icategory][ifatjet] = new TH1D(plotname,plotname,100,0,1);
            fatjet_Xccqq_modified[icategory][ifatjet]->Sumw2();
        }

        sprintf(plotname, "VBF_max_mass%i", icategory);
        VBF_max_mass[icategory] = new TH1D(plotname,plotname,100,500,2500);
        VBF_max_mass[icategory]->Sumw2();
    }

    fatjet_btag_score = new TH1D("fatjet_btag_score", "fatjet_btag_score", 100, 0, 1);
    fatjet_btag_score->Sumw2();
    first_fatjet_btag_score = new TH1D("first_fatjet_btag_score", "first_fatjet_btag_score", 100, 0, 1);
    first_fatjet_btag_score->Sumw2();
    second_fatjet_btag_score = new TH1D("second_fatjet_btag_score", "second_fatjet_btag_score", 100, 0, 1);
    second_fatjet_btag_score->Sumw2();
    third_fatjet_btag_score = new TH1D("third_fatjet_btag_score", "third_fatjet_btag_score", 100, 0, 1);
    third_fatjet_btag_score->Sumw2();

    first_fatjet_btag_score_2match = new TH1D("first_fatjet_btag_score_2match", "first_fatjet_btag_score_2match", 100, 0, 1);
    first_fatjet_btag_score_2match->Sumw2();
    second_fatjet_btag_score_2match = new TH1D("second_fatjet_btag_score_2match", "second_fatjet_btag_score_2match", 100, 0, 1);
    second_fatjet_btag_score_2match->Sumw2();
    third_fatjet_btag_score_2match = new TH1D("third_fatjet_btag_score_2match", "third_fatjet_btag_score_2match", 100, 0, 1);
    third_fatjet_btag_score_2match->Sumw2();

    first_fatjet_btag_score_1match = new TH1D("first_fatjet_btag_score_1match", "first_fatjet_btag_score_1match", 100, 0, 1);
    first_fatjet_btag_score_1match->Sumw2();
    second_fatjet_btag_score_1match = new TH1D("second_fatjet_btag_score_1match", "second_fatjet_btag_score_1match", 100, 0, 1);
    second_fatjet_btag_score_1match->Sumw2();
    third_fatjet_btag_score_1match = new TH1D("third_fatjet_btag_score_1match", "third_fatjet_btag_score_1match", 100, 0, 1);
    third_fatjet_btag_score_1match->Sumw2();

    first_fatjet_btag_score_0match = new TH1D("first_fatjet_btag_score_0match", "first_fatjet_btag_score_0match", 100, 0, 1);
    first_fatjet_btag_score_0match->Sumw2();
    second_fatjet_btag_score_0match = new TH1D("second_fatjet_btag_score_0match", "second_fatjet_btag_score_0match", 100, 0, 1);
    second_fatjet_btag_score_0match->Sumw2();
    third_fatjet_btag_score_0match = new TH1D("third_fatjet_btag_score_0match", "third_fatjet_btag_score_0match", 100, 0, 1);
    third_fatjet_btag_score_0match->Sumw2();

    jet_mass = new TH1D("jet_mass","jet_mass",100,0,200);
    jet_mass->Sumw2();

    VBFJet_leadingPt_2match = new TH1D("VBFJet_leadingPt_2match", "VBFJet_leadingPt_2match", 100, 0, 1000);
    VBFJet_leadingPt_2match->Sumw2();
    VBFJet_leadingPt_1match = new TH1D("VBFJet_leadingPt_1match", "VBFJet_leadingPt_1match", 100, 0, 1000);
    VBFJet_leadingPt_1match->Sumw2();
    VBFJet_leadingPt_0match = new TH1D("VBFJet_leadingPt_0match", "VBFJet_leadingPt_0match", 100, 0, 1000);
    VBFJet_leadingPt_0match->Sumw2();

    VBFJet_subleadingPt_2match = new TH1D("VBFJet_subleadingPt_2match", "VBFJet_subleadingPt_2match", 100, 0, 1000);
    VBFJet_subleadingPt_2match->Sumw2();
    VBFJet_subleadingPt_1match = new TH1D("VBFJet_subleadingPt_1match", "VBFJet_subleadingPt_1match", 100, 0, 1000);
    VBFJet_subleadingPt_1match->Sumw2();
    VBFJet_subleadingPt_0match = new TH1D("VBFJet_subleadingPt_0match", "VBFJet_subleadingPt_0match", 100, 0, 1000);
    VBFJet_subleadingPt_0match->Sumw2();

    VBFJet_Mjj_2match = new TH1D("VBFJet_Mjj_2match", "VBFJet_Mjj_2match", 100, 0, 6000);
    VBFJet_Mjj_2match->Sumw2();
    VBFJet_Mjj_1match = new TH1D("VBFJet_Mjj_1match", "VBFJet_Mjj_1match", 100, 0, 6000);
    VBFJet_Mjj_1match->Sumw2();
    VBFJet_Mjj_0match = new TH1D("VBFJet_Mjj_0match", "VBFJet_Mjj_0match", 100, 0, 6000);
    VBFJet_Mjj_0match->Sumw2();

    VBFJet_DeltaEta_2match_total = new TH1D("VBFJet_DeltaEta_2match_total", "VBFJet_DeltaEta_2match_total", 100, 0, 10);
    VBFJet_DeltaEta_2match_total->Sumw2();
    VBFJet_DeltaEta_1match_total = new TH1D("VBFJet_DeltaEta_1match_total", "VBFJet_DeltaEta_1match_total", 100, 0, 10);
    VBFJet_DeltaEta_1match_total->Sumw2();
    VBFJet_DeltaEta_0match_total = new TH1D("VBFJet_DeltaEta_0match_total", "VBFJet_DeltaEta_0match_total", 100, 0, 10);
    VBFJet_DeltaEta_0match_total->Sumw2();

    matched_VBFJet_Pt_total = new TH1D("matched_VBFJet_Pt_total", "matched_VBFJet_Pt_total", 100, 0, 1000);
    matched_VBFJet_Pt_total->Sumw2();
    matched_VBFJet_Eta_total = new TH1D("matched_VBFJet_Eta_total", "matched_VBFJet_Eta_total", 100, -5, 5);
    matched_VBFJet_Eta_total->Sumw2();
    matched_VBFJet_qgl_total = new TH1D("matched_VBFJet_qgl_total", "matched_VBFJet_qgl_total", 100, 0, 1);
    matched_VBFJet_qgl_total->Sumw2();

}

void makeHists::saveHists()
{
    hf->cd();
    hf->Write();
    hf->Close();
}
