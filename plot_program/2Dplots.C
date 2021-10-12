//Time:2021.4.15
//Last Modified:2021.4.15
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
  /*******Switches********/
  bool SetXRange = false;
  bool SetYRange = false;
  bool DoSave = true; //Save as pdf version
  bool DoRebin = false;

  char filename1[100] = "../outfiles/MC_selected.root";

  char plotname1[100] = "QJ_DeltaR_vs_HiggsMMM";

  char saveFileName[100] = "../plots/QJ_DeltaR_vs_HiggsMMM.pdf";

  char Xtitle[100] = "QJ DeltaR";
  char Ytitle[100] = "reco Higgs mass";

  char legend1[100] = "C2V_6";

  double minX = 25;  //if SetXRange
  double maxX = 120; //if SetXRange
  double minY = 0;   //if SetYRange
  double maxY = 1;   //if SetYRange

  int RebinNum = 2;

  /******Main functions*******/
  TFile *file1 = new TFile(filename1);
  TH1D *h1 = (TH1D *)file1->Get(plotname1);
  if (DoRebin)
    h1->Rebin(RebinNum);

  h1->SetBinContent(1, h1->GetBinContent(1) + h1->GetBinContent(0));
  h1->SetBinContent(h1->GetNbinsX(), h1->GetBinContent(h1->GetNbinsX() + 1) + h1->GetBinContent(h1->GetNbinsX()));

  /*********SetPad********/
  double canvasYsize, padYsize;
  canvasYsize = 500;
  padYsize = 0;
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, canvasYsize);
  TPad *pad1 = new TPad("pad1", "pad1", 0, padYsize, 1.0, 1.0); //Xbeginning,Ybeginning,Xending,Yending,100 per cent

  /*******SetCompare******/
  if (!SetXRange)
  {
    minX = h1->GetXaxis()->GetBinLowEdge(1);
    maxX = h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins());
  }
  h1->GetXaxis()->SetRangeUser(minX, maxX);
  pad1->Draw();
  pad1->cd();
  h1->Draw("colz");

  TLegend *legend = new TLegend(0.75, 0.75, 0.88, 0.86); //the coordination of the legend frame;
  legend->AddEntry(h1, legend1, "lpfe");
  legend->SetTextFont(70);
  legend->SetTextSize(0.04);
  legend->Draw("same");

  /*******h1 settings******/
  h1->SetTitle("");
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(2.0);
  h1->SetStats(0);
  h1->GetXaxis()->SetTitle(Xtitle);
  h1->GetYaxis()->SetTitle(Ytitle);
  h1->GetYaxis()->SetTitleSize(0.04);
  h1->GetYaxis()->SetTitleOffset(1);
  h1->GetYaxis()->SetTitleFont(70);
  h1->GetYaxis()->SetLabelSize(0.03);
  h1->GetYaxis()->SetLabelFont(70);
  h1->GetXaxis()->SetTitle(Xtitle);
  h1->GetXaxis()->SetTitleSize(0.04);
  h1->GetXaxis()->SetTitleOffset(1);
  h1->GetXaxis()->SetTitleFont(70); //the thickness of the title
  h1->GetXaxis()->SetLabelSize(0.03);
  h1->GetXaxis()->SetLabelFont(72);

  if (DoSave)
  {
    canvas->SaveAs(saveFileName);
  }
}
