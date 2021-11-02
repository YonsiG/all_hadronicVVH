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
}

void makeHists::saveHists()
{
    hf->cd();
    hf->Write();
    hf->Close();
}
