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
    char infileName[100] = "../../outfiles/C2V_3/OSWWH_0_selected.root";
    char outfileName[100] = "../../outfiles/C2V_3/OSWWH_scaled.root";

    TFile *inputFile = new TFile(infileName);
    TH1D *weight_Scale = (TH1D *)inputFile->Get("weight_Scale");
    double scaleNum = 1 / weight_Scale->GetBinContent(1);

    TFile *outputFile = new TFile(outfileName, "RECREATE");
    char plotname[50];
    TH1D *cutflow[14];
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

    for (int icategory = 0; icategory < 14; icategory++)
    {
        sprintf(plotname, "cutflow%i", icategory);
        cutflow[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "number_of_jets%i", icategory);
        number_of_jets[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        sprintf(plotname, "number_of_central_jets%i", icategory);
        number_of_central_jets[icategory] = (TH1D *)inputFile->Get(plotname)->Clone();
        cutflow[icategory]->Scale(scaleNum);
        number_of_jets[icategory]->Scale(scaleNum);
        number_of_central_jets[icategory]->Scale(scaleNum);
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

    inputFile->Close();
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();
}
