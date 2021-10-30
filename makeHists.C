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

    cutflow_1 = new TH1D("cutflow_1", "cutflow_1", 20, 0, 20);
    cutflow_2 = new TH1D("cutflow_2", "cutflow_2", 20, 0, 20);
    cutflow_3 = new TH1D("cutflow_3", "cutflow_3", 20, 0, 20);
    cutflow_4 = new TH1D("cutflow_4", "cutflow_4", 20, 0, 20);
    cutflow_5 = new TH1D("cutflow_5", "cutflow_5", 20, 0, 20);
    cutflow_6 = new TH1D("cutflow_6", "cutflow_6", 20, 0, 20);
    cutflow_7 = new TH1D("cutflow_7", "cutflow_7", 20, 0, 20);
    cutflow_8 = new TH1D("cutflow_8", "cutflow_8", 20, 0, 20);
    cutflow_9 = new TH1D("cutflow_9", "cutflow_9", 20, 0, 20);
    cutflow_10 = new TH1D("cutflow_10", "cutflow_10", 20, 0, 20);
    cutflow_11 = new TH1D("cutflow_11", "cutflow_11", 20, 0, 20);
    cutflow_12 = new TH1D("cutflow_12", "cutflow_12", 20, 0, 20);
    cutflow_13 = new TH1D("cutflow_13", "cutflow_13", 20, 0, 20);
    cutflow_14 = new TH1D("cutflow_14", "cutflow_14", 20, 0, 20);
    cutflow_15 = new TH1D("cutflow_15", "cutflow_15", 20, 0, 20);
    cutflow_16 = new TH1D("cutflow_16", "cutflow_16", 20, 0, 20);
    cutflow_17 = new TH1D("cutflow_17", "cutflow_17", 20, 0, 20);
}

void makeHists::saveHists()
{
    hf->cd();
    hf->Write();
    hf->Close();
}
