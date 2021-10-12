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

   TFile *signalfile = new TFile("../outfiles/C2V_4p5/C2V_4p5_selected.root");

   string histName = "recoHiggs_msoftdrop_stackTT";

   TCanvas *c = new TCanvas("plot1", "", 800, 600);
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1.0, 1.0);
   THStack *stack = new THStack("ts1", "");
   c->SetMargin(0.1, 0.05, 0.1, 0.05);

   TH1D *plot_2pass = (TH1D *)signalfile->Get("Higgs_msoftdrop_2passTT");
   TH1D *plot_1pass = (TH1D *)signalfile->Get("Higgs_msoftdrop_1passTT");
   TH1D *plot_0pass = (TH1D *)signalfile->Get("Higgs_msoftdrop_0passTT");

   for (int ibin = 1; ibin < plot_2pass->GetNbinsX() + 1; ibin++)
   {
      plot_2pass->SetBinError(ibin, 0);
      plot_1pass->SetBinError(ibin, 0);
      plot_0pass->SetBinError(ibin, 0);
   }

   plot_2pass->SetFillColor(kCyan);
   plot_1pass->SetFillColor(kPink);
   plot_0pass->SetFillColor(kYellow);

   /****************N1T*****************/
   pad1->Draw();
   if (isLog)
      pad1->SetLogy();
   pad1->cd();

   plot_2pass->SetTitle("");
   plot_2pass->SetLineWidth(2.0);
   plot_2pass->SetStats(0);
   plot_2pass->GetYaxis()->SetTitle("Entries");
   plot_2pass->GetYaxis()->SetTitleSize(0.04);
   plot_2pass->GetYaxis()->SetTitleOffset(1);
   plot_2pass->GetYaxis()->SetTitleFont(70);
   plot_2pass->GetYaxis()->SetLabelSize(0.03);
   plot_2pass->GetYaxis()->SetLabelFont(70);
   plot_2pass->GetXaxis()->SetTitle("Higgs msoftdrop / GeV");
   plot_2pass->GetXaxis()->SetTitleSize(0.04);
   plot_2pass->GetXaxis()->SetTitleOffset(1);
   plot_2pass->GetXaxis()->SetTitleFont(70);
   plot_2pass->GetXaxis()->SetLabelSize(0.03);
   plot_2pass->GetXaxis()->SetLabelFont(70);

   stack->Add(plot_2pass);
   stack->Add(plot_1pass);
   stack->Add(plot_0pass);
   stack->Draw("H");

   TLegend *legend = new TLegend(0.65, 0.75, 0.88, 0.88);
   legend->AddEntry(plot_2pass, "2pass", "f");
   legend->AddEntry(plot_1pass, "1pass", "f");
   legend->AddEntry(plot_0pass, "0pass", "f");
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
