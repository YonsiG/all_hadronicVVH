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
 
    TH1D *cutflow = (TH1D *)inputFile->Get("cutflow")->Clone();;
    cutflow->Scale(scaleNum);

    inputFile->Close();
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();
}
