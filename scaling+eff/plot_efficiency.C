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
   char infileName[100] = "../outfiles/C2V_4p5/C2V_4p5_selected.root";
   char outfileName[100] = "../outfiles/C2V_4p5/C2V_4p5_efficiency.root";
   TFile *inputFile = new TFile(infileName);

   TH1D *cutflow = (TH1D *)inputFile->Get("cutflow");
   TH1D *weight_Scale = (TH1D *)inputFile->Get("weight_Scale");

   double Nentries = cutflow->GetBinContent(1) / weight_Scale->GetBinContent(1);

   TH1D *total_Higgspt = (TH1D *)inputFile->Get("total_Higgspt");

   TH1D *cut_HiggsptLL = (TH1D *)inputFile->Get("cut_HiggsptLL");
   TH1D *cut_Higgspt_2passLL = (TH1D *)inputFile->Get("cut_Higgspt_2passLL");
   TH1D *cut_Higgspt_1passLL = (TH1D *)inputFile->Get("cut_Higgspt_1passLL");
   TH1D *cut_Higgspt_0passLL = (TH1D *)inputFile->Get("cut_Higgspt_0passLL");

   TH1D *cut_HiggsptMM = (TH1D *)inputFile->Get("cut_HiggsptMM");
   TH1D *cut_Higgspt_2passMM = (TH1D *)inputFile->Get("cut_Higgspt_2passMM");
   TH1D *cut_Higgspt_1passMM = (TH1D *)inputFile->Get("cut_Higgspt_1passMM");
   TH1D *cut_Higgspt_0passMM = (TH1D *)inputFile->Get("cut_Higgspt_0passMM");

   TH1D *cut_HiggsptTT = (TH1D *)inputFile->Get("cut_HiggsptTT");
   TH1D *cut_Higgspt_2passTT = (TH1D *)inputFile->Get("cut_Higgspt_2passTT");
   TH1D *cut_Higgspt_1passTT = (TH1D *)inputFile->Get("cut_Higgspt_1passTT");
   TH1D *cut_Higgspt_0passTT = (TH1D *)inputFile->Get("cut_Higgspt_0passTT");

   TH1D *total_quarkDeltaR = (TH1D *)inputFile->Get("total_quarkDeltaR");
   TH1D *cut_quarkDeltaRLL = (TH1D *)inputFile->Get("cut_quarkDeltaRLL");
   TH1D *cut_quarkDeltaRMM = (TH1D *)inputFile->Get("cut_quarkDeltaRMM");
   TH1D *cut_quarkDeltaRTT = (TH1D *)inputFile->Get("cut_quarkDeltaRTT");

   double xbin[14] = {0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
   double ybin[31];
   for (int ibin = 0; ibin < 20; ibin++)
   {
      ybin[ibin] = ibin * 50;
   }
   for (int ibin = 20; ibin < 31; ibin++)
   {
      ybin[ibin] = ibin * 100 - 1000;
   }

   TFile *outputFile = new TFile(outfileName, "RECREATE");

   TH1D *efficiency_vs_Higgs_PtLL = new TH1D("efficiency_vs_Higgs_PtLL", "efficiency_vs_Higgs_PtLL", 30, ybin);
   efficiency_vs_Higgs_PtLL->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_2passLL = new TH1D("efficiency_vs_Higgs_Pt_2passLL", "efficiency_vs_Higgs_Pt_2passLL", 30, ybin);
   efficiency_vs_Higgs_Pt_2passLL->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_1passLL = new TH1D("efficiency_vs_Higgs_Pt_1passLL", "efficiency_vs_Higgs_Pt_1passLL", 30, ybin);
   efficiency_vs_Higgs_Pt_1passLL->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_0passLL = new TH1D("efficiency_vs_Higgs_Pt_0passLL", "efficiency_vs_Higgs_Pt_0passLL", 30, ybin);
   efficiency_vs_Higgs_Pt_0passLL->Sumw2();

   TH1D *efficiency_vs_Higgs_PtMM = new TH1D("efficiency_vs_Higgs_PtMM", "efficiency_vs_Higgs_PtMM", 30, ybin);
   efficiency_vs_Higgs_PtMM->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_2passMM = new TH1D("efficiency_vs_Higgs_Pt_2passMM", "efficiency_vs_Higgs_Pt_2passMM", 30, ybin);
   efficiency_vs_Higgs_Pt_2passMM->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_1passMM = new TH1D("efficiency_vs_Higgs_Pt_1passMM", "efficiency_vs_Higgs_Pt_1passMM", 30, ybin);
   efficiency_vs_Higgs_Pt_1passMM->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_0passMM = new TH1D("efficiency_vs_Higgs_Pt_0passMM", "efficiency_vs_Higgs_Pt_0passMM", 30, ybin);
   efficiency_vs_Higgs_Pt_0passMM->Sumw2();

   TH1D *efficiency_vs_Higgs_PtTT = new TH1D("efficiency_vs_Higgs_PtTT", "efficiency_vs_Higgs_PtTT", 30, ybin);
   efficiency_vs_Higgs_PtTT->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_2passTT = new TH1D("efficiency_vs_Higgs_Pt_2passTT", "efficiency_vs_Higgs_Pt_2passTT", 30, ybin);
   efficiency_vs_Higgs_Pt_2passTT->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_1passTT = new TH1D("efficiency_vs_Higgs_Pt_1passTT", "efficiency_vs_Higgs_Pt_1passTT", 30, ybin);
   efficiency_vs_Higgs_Pt_1passTT->Sumw2();

   TH1D *efficiency_vs_Higgs_Pt_0passTT = new TH1D("efficiency_vs_Higgs_Pt_0passTT", "efficiency_vs_Higgs_Pt_0passTT", 30, ybin);
   efficiency_vs_Higgs_Pt_0passTT->Sumw2();

   TH1D *efficiency_vs_quark_DeltaRLL = new TH1D("efficiency_vs_quark_DeltaRLL", "efficiency_vs_quark_DeltaRLL", 13, xbin);
   efficiency_vs_quark_DeltaRLL->Sumw2();

   TH1D *efficiency_vs_quark_DeltaRMM = new TH1D("efficiency_vs_quark_DeltaRMM", "efficiency_vs_quark_DeltaRMM", 13, xbin);
   efficiency_vs_quark_DeltaRMM->Sumw2();

   TH1D *efficiency_vs_quark_DeltaRTT = new TH1D("efficiency_vs_quark_DeltaRTT", "efficiency_vs_quark_DeltaRTT", 13, xbin);
   efficiency_vs_quark_DeltaRTT->Sumw2();

   for (int ibin = 1; ibin < 31; ibin++)
   {
      double total = total_Higgspt->GetBinContent(ibin);

      double cut_LL = cut_HiggsptLL->GetBinContent(ibin);
      double cut_2passLL = cut_Higgspt_2passLL->GetBinContent(ibin);
      double cut_1passLL = cut_Higgspt_1passLL->GetBinContent(ibin);
      double cut_0passLL = cut_Higgspt_0passLL->GetBinContent(ibin);

      double cut_MM = cut_HiggsptMM->GetBinContent(ibin);
      double cut_2passMM = cut_Higgspt_2passMM->GetBinContent(ibin);
      double cut_1passMM = cut_Higgspt_1passMM->GetBinContent(ibin);
      double cut_0passMM = cut_Higgspt_0passMM->GetBinContent(ibin);

      double cut_TT = cut_HiggsptTT->GetBinContent(ibin);
      double cut_2passTT = cut_Higgspt_2passTT->GetBinContent(ibin);
      double cut_1passTT = cut_Higgspt_1passTT->GetBinContent(ibin);
      double cut_0passTT = cut_Higgspt_0passTT->GetBinContent(ibin);

      if (total > 0)
      {
         double eff = cut_LL / total;
         efficiency_vs_Higgs_PtLL->SetBinContent(ibin, eff);
         efficiency_vs_Higgs_PtLL->SetBinError(ibin, sqrt(eff * (1 - eff) / Nentries));

         double eff2 = cut_2passLL / total;
         efficiency_vs_Higgs_Pt_2passLL->SetBinContent(ibin, eff2);
         efficiency_vs_Higgs_Pt_2passLL->SetBinError(ibin, sqrt(eff2 * (1 - eff2) / Nentries));

         double eff1 = cut_1passLL / total;
         efficiency_vs_Higgs_Pt_1passLL->SetBinContent(ibin, eff1);
         efficiency_vs_Higgs_Pt_1passLL->SetBinError(ibin, sqrt(eff1 * (1 - eff1) / Nentries));

         double eff0 = cut_0passLL / total;
         if (eff0 < 1)
            efficiency_vs_Higgs_Pt_0passLL->SetBinContent(ibin, eff0);
         else
            efficiency_vs_Higgs_Pt_0passLL->SetBinContent(ibin, 0);
         //    myHists->efficiency_vs_Higgs_Pt_0passLL->SetBinError(ibin,sqrt(eff0*(1-eff0)/Nentries));
         efficiency_vs_Higgs_Pt_0passTT->SetBinError(ibin, 0); //manually set to 0, as it's not useful, only draw a stack plot.

         eff = cut_MM / total;
         efficiency_vs_Higgs_PtMM->SetBinContent(ibin, eff);
         efficiency_vs_Higgs_PtMM->SetBinError(ibin, sqrt(eff * (1 - eff) / Nentries));

         eff2 = cut_2passMM / total;
         efficiency_vs_Higgs_Pt_2passMM->SetBinContent(ibin, eff2);
         efficiency_vs_Higgs_Pt_2passMM->SetBinError(ibin, sqrt(eff2 * (1 - eff2) / Nentries));

         eff1 = cut_1passMM / total;
         efficiency_vs_Higgs_Pt_1passMM->SetBinContent(ibin, eff1);
         efficiency_vs_Higgs_Pt_1passMM->SetBinError(ibin, sqrt(eff1 * (1 - eff1) / Nentries));

         eff0 = cut_0passMM / total;
         if (eff0 < 1)
            efficiency_vs_Higgs_Pt_0passMM->SetBinContent(ibin, eff0);
         else
            efficiency_vs_Higgs_Pt_0passMM->SetBinContent(ibin, 0);
         //    myHists->efficiency_vs_Higgs_Pt_0passMM->SetBinError(ibin,sqrt(eff0*(1-eff0)/Nentries));
         efficiency_vs_Higgs_Pt_0passMM->SetBinError(ibin, 0); //manually set to 0, as it's not useful, only draw a stack plot.

         eff = cut_TT / total;
         efficiency_vs_Higgs_PtTT->SetBinContent(ibin, eff);
         efficiency_vs_Higgs_PtTT->SetBinError(ibin, sqrt(eff * (1 - eff) / Nentries));

         eff2 = cut_2passTT / total;
         efficiency_vs_Higgs_Pt_2passTT->SetBinContent(ibin, eff2);
         efficiency_vs_Higgs_Pt_2passTT->SetBinError(ibin, sqrt(eff2 * (1 - eff2) / Nentries));

         eff1 = cut_1passTT / total;
         efficiency_vs_Higgs_Pt_1passTT->SetBinContent(ibin, eff1);
         efficiency_vs_Higgs_Pt_1passTT->SetBinError(ibin, sqrt(eff1 * (1 - eff1) / Nentries));

         eff0 = cut_0passTT / total;
         if (eff0 < 1)
            efficiency_vs_Higgs_Pt_0passTT->SetBinContent(ibin, eff0);
         else
            efficiency_vs_Higgs_Pt_0passTT->SetBinContent(ibin, 0);
         //    myHists->efficiency_vs_Higgs_Pt_0passTT->SetBinError(ibin,sqrt(eff0*(1-eff0)/Nentries));
         efficiency_vs_Higgs_Pt_0passTT->SetBinError(ibin, 0); //manually set to 0, as it's not useful, only draw a stack plot.
      }
   }
   for (int ibin = 1; ibin < 14; ibin++)
   {
      double total_DeltaR = total_quarkDeltaR->GetBinContent(ibin);
      double cut_DeltaRLL = cut_quarkDeltaRLL->GetBinContent(ibin);
      double cut_DeltaRMM = cut_quarkDeltaRMM->GetBinContent(ibin);
      double cut_DeltaRTT = cut_quarkDeltaRTT->GetBinContent(ibin);

      if (total_DeltaR > 0)
      {
         double eff = cut_DeltaRLL / total_DeltaR;
         efficiency_vs_quark_DeltaRLL->SetBinContent(ibin, eff);
         efficiency_vs_quark_DeltaRLL->SetBinError(ibin, sqrt(eff * (1 - eff) / Nentries));

         eff = cut_DeltaRMM / total_DeltaR;
         efficiency_vs_quark_DeltaRMM->SetBinContent(ibin, eff);
         efficiency_vs_quark_DeltaRMM->SetBinError(ibin, sqrt(eff * (1 - eff) / Nentries));

         eff = cut_DeltaRTT / total_DeltaR;
         efficiency_vs_quark_DeltaRTT->SetBinContent(ibin, eff);
         efficiency_vs_quark_DeltaRTT->SetBinError(ibin, sqrt(eff * (1 - eff) / Nentries));
      }
   }

   inputFile->Close();
   outputFile->Write();
   outputFile->Close();
}
