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

        sprintf(plotname, "number_of_jets%i", icategory);
        number_of_jets[icategory] = new TH1D(plotname, plotname, 10, 0, 10);

        sprintf(plotname, "number_of_central_jets%i", icategory);
        number_of_central_jets[icategory] = new TH1D(plotname, plotname, 10, 0, 10);
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

    VBFJet_DeltaEta_2match = new TH1D("VBFJet_DeltaEta_2match", "VBFJet_DeltaEta_2match", 100, 0, 10);
    VBFJet_DeltaEta_2match->Sumw2();
    VBFJet_DeltaEta_1match = new TH1D("VBFJet_DeltaEta_1match", "VBFJet_DeltaEta_1match", 100, 0, 10);
    VBFJet_DeltaEta_1match->Sumw2();
    VBFJet_DeltaEta_0match = new TH1D("VBFJet_DeltaEta_0match", "VBFJet_DeltaEta_0match", 100, 0, 10);
    VBFJet_DeltaEta_0match->Sumw2();
}

void makeHists::saveHists()
{
    hf->cd();
    hf->Write();
    hf->Close();
}
