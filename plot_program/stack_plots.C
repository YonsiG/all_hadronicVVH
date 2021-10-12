//Time:2021.6.29
//Last Modified:2021.6.29

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
#include "THStack.h"
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
   bool isLog = 0;

   TFile *signalfile = new TFile("../outfiles/C2V_4p5/C2V_4p5_scaled.root");
//   TFile *ttWfile = new TFile("../outfiles/ttW_scaled.root");
   TFile *TTJets_DiLeptfile = new TFile("../outfiles/TTJets_DiLept/TTJets_DiLept_scaled.root");
   TFile *TTJets_SingleLeptFromTfile = new TFile("../outfiles/TTJets_SingleLeptFromT/TTJets_SingleLeptFromT_scaled.root");
   TFile *TTJets_SingleLeptFromTbarfile = new TFile("../outfiles/TTJets_SingleLeptFromTbar/TTJets_SingleLeptFromTbar_scaled.root");

   string histName = "recoHiggs_msoftdrop_s+bTT";

   TCanvas *c = new TCanvas("plot1", "", 800, 600);
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1.0, 1.0);
   THStack *stack = new THStack("ts1", "");
   //   stack->GetYaxis()->SetRangeUser(1,3000);
   c->SetMargin(0.1, 0.05, 0.1, 0.05);

   TH1D *signalplot = (TH1D *)signalfile->Get("recoHiggs_msoftdropTT");
//   TH1D *ttWplot = (TH1D *)ttWfile->Get("recoHiggs_msoftdropTT");
   TH1D *TTJets_DiLeptplot = (TH1D *)TTJets_DiLeptfile->Get("recoHiggs_msoftdropTT");
   TH1D *TTJets_SingleLeptFromTplot = (TH1D *)TTJets_SingleLeptFromTfile->Get("recoHiggs_msoftdropTT");
   TH1D *TTJets_SingleLeptFromTbarplot = (TH1D *)TTJets_SingleLeptFromTbarfile->Get("recoHiggs_msoftdropTT");
   TH1D *dataplot = (TH1D *)TTJets_DiLeptplot->Clone();
// = (TH1D *)ttWplot->Clone();
//   dataplot->Add(TTJets_DiLeptplot);
   dataplot->Add(TTJets_SingleLeptFromTplot);
   dataplot->Add(TTJets_SingleLeptFromTbarplot);

   signalplot->Rebin(4);
//   ttWplot->Rebin(4);
   TTJets_DiLeptplot->Rebin(4);
   TTJets_SingleLeptFromTplot->Rebin(4);
   TTJets_SingleLeptFromTbarplot->Rebin(4);
   dataplot->Rebin(4);

   for (int ibin = 1; ibin < TTJets_DiLeptplot->GetNbinsX() + 1; ibin++)
   {
//      ttWplot->SetBinError(ibin, 0);
      TTJets_DiLeptplot->SetBinError(ibin, 0);
      TTJets_SingleLeptFromTplot->SetBinError(ibin, 0);
      TTJets_SingleLeptFromTbarplot->SetBinError(ibin, 0);
   }

//   ttWplot->SetFillColor(kOrange);
   TTJets_DiLeptplot->SetFillColor(kOrange + 1);
   TTJets_SingleLeptFromTplot->SetFillColor(kOrange + 2);
   TTJets_SingleLeptFromTbarplot->SetFillColor(kOrange + 3);

   /****************N1T*****************/
   pad1->Draw();
   if (isLog)
      pad1->SetLogy();
   pad1->cd();

   signalplot->GetYaxis()->SetRangeUser(10e-2, 3000);
   signalplot->GetXaxis()->SetRangeUser(0,1000);
   //   signalplot->GetYaxis()->SetMinimum(10-e2);
   signalplot->SetTitle("");
   signalplot->SetLineWidth(2.0);
   signalplot->SetStats(0);
   signalplot->SetLineColor(kRed + 1);
   signalplot->GetYaxis()->SetTitle("Entries");
   signalplot->GetYaxis()->SetTitleSize(0.04);
   signalplot->GetYaxis()->SetTitleOffset(1);
   signalplot->GetYaxis()->SetTitleFont(70);
   signalplot->GetYaxis()->SetLabelSize(0.03);
   signalplot->GetYaxis()->SetLabelFont(70);
   signalplot->GetXaxis()->SetTitle("Higgs msoftdrop / GeV");
   signalplot->GetXaxis()->SetTitleSize(0.04);
   signalplot->GetXaxis()->SetTitleOffset(1);
   signalplot->GetXaxis()->SetTitleFont(70);
   signalplot->GetXaxis()->SetLabelSize(0.03);
   signalplot->GetXaxis()->SetLabelFont(70);

   dataplot->SetLineColor(kBlack);
   dataplot->SetMarkerStyle(20);
   dataplot->GetXaxis()->SetRangeUser(0,1000);
   dataplot->GetYaxis()->SetRangeUser(10e-2, 3000);
//   ttWplot->GetYaxis()->SetRangeUser(10e-2, 3000);
   TTJets_DiLeptplot->GetXaxis()->SetRangeUser(0,1000);
   TTJets_SingleLeptFromTplot->GetXaxis()->SetRangeUser(0, 1000);
   TTJets_SingleLeptFromTplot->GetYaxis()->SetRangeUser(10e-2, 3000);
   TTJets_SingleLeptFromTbarplot->GetXaxis()->SetRangeUser(0, 1000);
   TTJets_SingleLeptFromTbarplot->GetYaxis()->SetRangeUser(10e-2, 3000);

   signalplot->Draw("HIST");
//   stack->Add(ttWplot);
   stack->Add(TTJets_DiLeptplot);
   stack->Add(TTJets_SingleLeptFromTplot);
   stack->Add(TTJets_SingleLeptFromTbarplot);
   stack->Draw("H");
   dataplot->Draw("E0 same");
   signalplot->Draw("same HIST");

   TLegend *legend = new TLegend(0.65, 0.75, 0.88, 0.88);
   legend->AddEntry(signalplot, "signal", "lpfe");
   legend->AddEntry(dataplot, "data", "lpfe");
//   legend->AddEntry(ttWplot, "ttW", "f");
   legend->AddEntry(TTJets_DiLeptplot, "TTJets_DiLept", "f");
   legend->AddEntry(TTJets_SingleLeptFromTplot, "TTJets_SingleLeptFromT", "f");
   legend->AddEntry(TTJets_SingleLeptFromTbarplot, "TTJets_SingleLeptFromTbar", "f");
   legend->SetTextFont(20);
   legend->SetTextSize(0.025);
   legend->SetBorderSize(0);
   legend->Draw("same");
   //   legend -> SetHeader("ATLAS Internal");

   if (isLog)
      c->SaveAs(("../plots/" + histName + "_log.pdf").c_str()); //save plots to disk
   if (!isLog)
      c->SaveAs(("../plots/" + histName + "_linear.pdf").c_str());
}
