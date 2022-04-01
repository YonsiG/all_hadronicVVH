// Time:2021.6.29
// Last Modified:2021.6.29

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
   bool isLog = 1;

   TFile *signalfile = new TFile("../../outfiles/XbbMD/selection4/signal_merged_2_scaled.root");
   TFile *QCDfile = new TFile("../../outfiles/XbbMD/selection4/QCD_merged_2_scaled.root");
   TFile *TTToHadronicfile = new TFile("../../outfiles/XbbMD/selection4/TTToHadronic_2_scaled.root");
   TFile *TTToSemiLeptonicfile = new TFile("../../outfiles/XbbMD/selection4/TTToSemiLeptonic_2_scaled.root");
   TFile *WJetfile = new TFile("../../outfiles/XbbMD/selection4/WJet_merged_2_scaled.root");
   TFile *ZJetfile = new TFile("../../outfiles/XbbMD/selection4/ZJet_merged_2_scaled.root");
   TFile *datafile = new TFile("../../outfiles/XbbMD/selection4/data_2_selected.root");

   string histName = "fatjet_pt_s+b";

   TCanvas *c = new TCanvas("plot1", "", 800, 600);
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1.0, 1.0);
   THStack *stack = new THStack("ts1", "");
   //   stack->GetYaxis()->SetRangeUser(1,3000);
   c->SetMargin(0.1, 0.05, 0.1, 0.05);

   gStyle->SetOptStat(0000);

   TH1D *signalplot = (TH1D *)signalfile->Get("ST0");
   TH1D *QCDplot = (TH1D *)QCDfile->Get("ST0");
   TH1D *TTToHadronicplot = (TH1D *)TTToHadronicfile->Get("ST0");
   TH1D *TTToSemiLeptonicplot = (TH1D *)TTToSemiLeptonicfile->Get("ST0");
   TH1D *WJetplot = (TH1D *)WJetfile->Get("ST0");
   TH1D *ZJetplot = (TH1D *)ZJetfile->Get("ST0");
   TH1D *dataplot = (TH1D *)datafile->Get("ST0");

   signalplot->Rebin(4);
   QCDplot->Rebin(4);
   TTToHadronicplot->Rebin(4);
   TTToSemiLeptonicplot->Rebin(4);
   WJetplot->Rebin(4);
   ZJetplot->Rebin(4);
   dataplot->Rebin(4);

   for (int ibin = 1; ibin < TTToHadronicplot->GetNbinsX() + 1; ibin++)
   {
      QCDplot->SetBinError(ibin, 0);
      TTToHadronicplot->SetBinError(ibin, 0);
      TTToSemiLeptonicplot->SetBinError(ibin, 0);
      WJetplot->SetBinError(ibin, 0);
      ZJetplot->SetBinError(ibin, 0);
   }

   QCDplot->SetFillColor(kOrange);
   TTToHadronicplot->SetFillColor(kOrange + 1);
   TTToSemiLeptonicplot->SetFillColor(kOrange + 2);
   WJetplot->SetFillColor(kOrange + 3);
   ZJetplot->SetFillColor(kOrange + 4);

   /****************N1T*****************/
   pad1->Draw();
   if (isLog)
      pad1->SetLogy();
   pad1->cd();

   signalplot->GetYaxis()->SetRangeUser(10e-2, 30);
   signalplot->GetXaxis()->SetRangeUser(-5, 1000);
   //   signalplot->GetYaxis()->SetMinimum(10e-2);
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
   signalplot->GetXaxis()->SetTitle("ST / GeV");
   signalplot->GetXaxis()->SetTitleSize(0.04);
   signalplot->GetXaxis()->SetTitleOffset(1);
   signalplot->GetXaxis()->SetTitleFont(70);
   signalplot->GetXaxis()->SetLabelSize(0.03);
   signalplot->GetXaxis()->SetLabelFont(70);

   dataplot->SetLineColor(kBlack);
   dataplot->SetMarkerStyle(20);
   //   dataplot->GetXaxis()->SetRangeUser(-5,1000);
   dataplot->GetYaxis()->SetRangeUser(10e-2, 30);
   //   QCDplot->GetXaxis()->SetRangeUser(0,1000);
   QCDplot->GetYaxis()->SetRangeUser(10e-2, 30);
   //   TTToHadronicplot->GetXaxis()->SetRangeUser(0,1000);
   TTToHadronicplot->GetYaxis()->SetRangeUser(10e-2, 30);
   //   TTToSemiLeptonicplot->GetXaxis()->SetRangeUser(0, 1000);
   TTToSemiLeptonicplot->GetYaxis()->SetRangeUser(10e-2, 30);
   //   WJetplot->GetXaxis()->SetRangeUser(0, 1000);
   WJetplot->GetYaxis()->SetRangeUser(10e-2, 30);
   //   ZJetplot->GetXaxis()->SetRangeUser(0, 1000);
   ZJetplot->GetYaxis()->SetRangeUser(10e-2, 30);

   stack->Add(QCDplot);
   stack->Add(TTToHadronicplot);
   stack->Add(TTToSemiLeptonicplot);
   stack->Add(WJetplot);
   stack->Add(ZJetplot);
   dataplot->Draw("E0");
   stack->Draw("H same");
   dataplot->Draw("E0 same");
   signalplot->Draw("same HIST");

   TLegend *legend = new TLegend(0.65, 0.75, 0.88, 0.88);
   legend->AddEntry(signalplot, "signal", "lpfe");
   legend->AddEntry(dataplot, "data", "lpfe");
   legend->AddEntry(QCDplot, "QCD", "f");
   legend->AddEntry(TTToHadronicplot, "TTToHadronic", "f");
   legend->AddEntry(TTToSemiLeptonicplot, "TTToSemiLeptonic", "f");
   legend->AddEntry(WJetplot, "WJet", "f");
   legend->AddEntry(ZJetplot, "ZJet", "f");
   legend->SetTextFont(20);
   legend->SetTextSize(0.025);
   legend->SetBorderSize(0);
   legend->Draw("same");
   //   legend -> SetHeader("ATLAS Internal");

   if (isLog)
      c->SaveAs(("../../plots/" + histName + "_log.pdf").c_str()); // save plots to disk
   if (!isLog)
      c->SaveAs(("../../plots/" + histName + "_linear.pdf").c_str());
}
