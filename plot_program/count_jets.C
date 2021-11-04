//Time:2017.10.11
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
  bool NormalizeToOne = true; //Set Compare mode first
  bool DoSave = true;         //Save as pdf version

  char filename1[100] = "../../outfiles/C2V_3/OSWWH_scaled.root";

  char plotname1[100] = "number_of_central_jets0";
  char plotname2[100] = "number_of_central_jets1";
  char plotname3[100] = "number_of_central_jets2";

  char plotname4[100] = "number_of_central_jets3";
  char plotname5[100] = "number_of_central_jets4";
  char plotname6[100] = "number_of_central_jets5";
  char plotname7[100] = "number_of_central_jets6";

  char plotname8[100] = "number_of_central_jets7";
  char plotname9[100] = "number_of_central_jets8";
  char plotname10[100] = "number_of_central_jets9";
  char plotname11[100] = "number_of_central_jets10";
  char plotname12[100] = "number_of_central_jets11";
  char plotname13[100] = "number_of_central_jets12";

  char saveFileName1[100] = "../../plots/C2V_3/OSWWH_3fatjet_ncentral_jets.pdf";
  char saveFileName2[100] = "../../plots/C2V_3/OSWWH_2fatjet_ncentral_jets.pdf";
  char saveFileName3[100] = "../../plots/C2V_3/OSWWH_1fatjet_ncentral_jets.pdf";

  char Xtitle[100] = "number of central_jets";
  char Ytitle[100] = "entries";

  char legend1[100] = "OSWWH";

  double minX = 25;    //if SetXRange
  double maxX = 120;   //if SetXRange
  double minY = 10e-3; //if SetYRange
  double maxY = 2;     //if SetYRange

  /******Main functions*******/
  TFile *file1 = new TFile(filename1);
  TH1D *h1 = (TH1D *)file1->Get(plotname1);
  TH1D *h2 = (TH1D *)file1->Get(plotname2);
  TH1D *h3 = (TH1D *)file1->Get(plotname3);
  TH1D *h4 = (TH1D *)file1->Get(plotname4);
  TH1D *h5 = (TH1D *)file1->Get(plotname5);
  TH1D *h6 = (TH1D *)file1->Get(plotname6);
  TH1D *h7 = (TH1D *)file1->Get(plotname7);
  TH1D *h8 = (TH1D *)file1->Get(plotname8);
  TH1D *h9 = (TH1D *)file1->Get(plotname9);
  TH1D *h10 = (TH1D *)file1->Get(plotname10);
  TH1D *h11 = (TH1D *)file1->Get(plotname11);
  TH1D *h12 = (TH1D *)file1->Get(plotname12);
  TH1D *h13 = (TH1D *)file1->Get(plotname13);

  h1->SetBinContent(1, h1->GetBinContent(1) + h1->GetBinContent(0));
  h1->SetBinContent(h1->GetNbinsX(), h1->GetBinContent(h1->GetNbinsX() + 1) + h1->GetBinContent(h1->GetNbinsX()));
  h2->SetBinContent(1, h2->GetBinContent(1) + h2->GetBinContent(0));
  h2->SetBinContent(h2->GetNbinsX(), h2->GetBinContent(h2->GetNbinsX() + 1) + h2->GetBinContent(h2->GetNbinsX()));
  h3->SetBinContent(1, h3->GetBinContent(1) + h3->GetBinContent(0));
  h3->SetBinContent(h3->GetNbinsX(), h3->GetBinContent(h3->GetNbinsX() + 1) + h3->GetBinContent(h3->GetNbinsX()));

  h4->SetBinContent(1, h4->GetBinContent(1) + h4->GetBinContent(0));
  h4->SetBinContent(h4->GetNbinsX(), h4->GetBinContent(h4->GetNbinsX() + 1) + h4->GetBinContent(h4->GetNbinsX()));
  h5->SetBinContent(1, h5->GetBinContent(1) + h5->GetBinContent(0));
  h5->SetBinContent(h5->GetNbinsX(), h5->GetBinContent(h5->GetNbinsX() + 1) + h5->GetBinContent(h5->GetNbinsX()));
  h6->SetBinContent(1, h6->GetBinContent(1) + h6->GetBinContent(0));
  h6->SetBinContent(h6->GetNbinsX(), h6->GetBinContent(h6->GetNbinsX() + 1) + h6->GetBinContent(h6->GetNbinsX()));
  h7->SetBinContent(1, h7->GetBinContent(1) + h7->GetBinContent(0));
  h7->SetBinContent(h7->GetNbinsX(), h7->GetBinContent(h7->GetNbinsX() + 1) + h7->GetBinContent(h7->GetNbinsX()));

  h8->SetBinContent(1, h8->GetBinContent(1) + h8->GetBinContent(0));
  h8->SetBinContent(h8->GetNbinsX(), h8->GetBinContent(h8->GetNbinsX() + 1) + h8->GetBinContent(h8->GetNbinsX()));
  h9->SetBinContent(1, h9->GetBinContent(1) + h9->GetBinContent(0));
  h9->SetBinContent(h9->GetNbinsX(), h9->GetBinContent(h9->GetNbinsX() + 1) + h9->GetBinContent(h9->GetNbinsX()));
  h10->SetBinContent(1, h10->GetBinContent(1) + h10->GetBinContent(0));
  h10->SetBinContent(h10->GetNbinsX(), h10->GetBinContent(h10->GetNbinsX() + 1) + h10->GetBinContent(h10->GetNbinsX()));
  h11->SetBinContent(1, h11->GetBinContent(1) + h11->GetBinContent(0));
  h11->SetBinContent(h11->GetNbinsX(), h11->GetBinContent(h11->GetNbinsX() + 1) + h11->GetBinContent(h11->GetNbinsX()));
  h12->SetBinContent(1, h12->GetBinContent(1) + h12->GetBinContent(0));
  h12->SetBinContent(h12->GetNbinsX(), h12->GetBinContent(h12->GetNbinsX() + 1) + h12->GetBinContent(h12->GetNbinsX()));
  h13->SetBinContent(1, h13->GetBinContent(1) + h13->GetBinContent(0));
  h13->SetBinContent(h13->GetNbinsX(), h13->GetBinContent(h13->GetNbinsX() + 1) + h13->GetBinContent(h13->GetNbinsX()));

  cout << "Overflow: " << h1->GetBinContent(h1->GetNbinsX() + 1) << endl;

  h1->Add(h2, 1);
  h1->Add(h3, 1);

  h4->Add(h5, 1);
  h4->Add(h6, 1);
  h4->Add(h7, 1);

  h8->Add(h9, 1);
  h8->Add(h10, 1);
  h8->Add(h11, 1);
  h8->Add(h12, 1);
  h8->Add(h13, 1);

  /*********SetPad********/
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 500);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1.0, 1.0); //Xbeginning,Ybeginning,Xending,Yending,100 per cent

  /*******SetCompare******/
  if (!SetXRange)
  {
    minX = h1->GetXaxis()->GetBinLowEdge(1);
    maxX = h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetNbins());
  }
  h1->GetXaxis()->SetRangeUser(minX, maxX);
  pad1->Draw();
  pad1->cd();
  h1->Draw("TEXT45 HIST .2g");

  //   TLegend* legend = new TLegend(0.75,0.75,0.88,0.86); //the coordination of the legend frame;
  TLegend *legend = new TLegend(0.73, 0.73, 0.85, 0.8);
  legend->AddEntry(h1, legend1, "lpfe");
  legend->SetTextFont(70);
  legend->SetTextSize(0.04);
  legend->Draw("same");

  if (NormalizeToOne)
  {
    double inte1 = h1->Integral();
    h1->Scale(1 / inte1);
  }

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

  /*********SetPad********/
  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 800, 500);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1.0, 1.0); //Xbeginning,Ybeginning,Xending,Yending,100 per cent

  /*******SetCompare******/
  if (!SetXRange)
  {
    minX = h4->GetXaxis()->GetBinLowEdge(1);
    maxX = h4->GetXaxis()->GetBinUpEdge(h4->GetXaxis()->GetNbins());
  }
  h4->GetXaxis()->SetRangeUser(minX, maxX);
  pad2->Draw();
  pad2->cd();
  h4->Draw("TEXT45.2f HIST");

  //   TLegend* legend2 = new TLegend(0.75,0.75,0.88,0.86); //the coordination of the legend2 frame;
  TLegend *legend2 = new TLegend(0.73, 0.73, 0.85, 0.8);
  legend2->AddEntry(h4, legend1, "lpfe");
  legend2->SetTextFont(70);
  legend2->SetTextSize(0.04);
  legend2->Draw("same");

  if (NormalizeToOne)
  {
    double inte1 = h4->Integral();
    h4->Scale(1 / inte1);
  }

  /*******h4 settings******/
  h4->SetTitle("");
  h4->SetLineColor(kBlue);
  h4->SetLineWidth(2.0);
  h4->SetStats(0);
  h4->GetXaxis()->SetTitle(Xtitle);
  h4->GetYaxis()->SetTitle(Ytitle);
  h4->GetYaxis()->SetTitleSize(0.04);
  h4->GetYaxis()->SetTitleOffset(1);
  h4->GetYaxis()->SetTitleFont(70);
  h4->GetYaxis()->SetLabelSize(0.03);
  h4->GetYaxis()->SetLabelFont(70);
  h4->GetXaxis()->SetTitle(Xtitle);
  h4->GetXaxis()->SetTitleSize(0.04);
  h4->GetXaxis()->SetTitleOffset(1);
  h4->GetXaxis()->SetTitleFont(70); //the thickness of the title
  h4->GetXaxis()->SetLabelSize(0.03);
  h4->GetXaxis()->SetLabelFont(72);

  /*********SetPad********/
  TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 800, 500);
  TPad *pad3 = new TPad("pad3", "pad3", 0, 0, 1.0, 1.0); //Xbeginning,Ybeginning,Xending,Yending,100 per cent

  /*******SetCompare******/
  if (!SetXRange)
  {
    minX = h8->GetXaxis()->GetBinLowEdge(1);
    maxX = h8->GetXaxis()->GetBinUpEdge(h8->GetXaxis()->GetNbins());
  }
  h8->GetXaxis()->SetRangeUser(minX, maxX);
  pad3->Draw();
  pad3->cd();
  h8->Draw("TEXT45.2f HIST");

  //   TLegend* legend = new TLegend(0.75,0.75,0.88,0.86); //the coordination of the legend frame;
  TLegend *legend3 = new TLegend(0.73, 0.73, 0.85, 0.8);
  legend3->AddEntry(h8, legend1, "lpfe");
  legend3->SetTextFont(70);
  legend3->SetTextSize(0.04);
  legend3->Draw("same");

  if (NormalizeToOne)
  {
    double inte1 = h8->Integral();
    h8->Scale(1 / inte1);
  }

  /*******h8 settings******/
  h8->SetTitle("");
  h8->SetLineColor(kBlue);
  h8->SetLineWidth(2.0);
  h8->SetStats(0);
  h8->GetXaxis()->SetTitle(Xtitle);
  h8->GetYaxis()->SetTitle(Ytitle);
  h8->GetYaxis()->SetTitleSize(0.04);
  h8->GetYaxis()->SetTitleOffset(1);
  h8->GetYaxis()->SetTitleFont(70);
  h8->GetYaxis()->SetLabelSize(0.03);
  h8->GetYaxis()->SetLabelFont(70);
  h8->GetXaxis()->SetTitle(Xtitle);
  h8->GetXaxis()->SetTitleSize(0.04);
  h8->GetXaxis()->SetTitleOffset(1);
  h8->GetXaxis()->SetTitleFont(70); //the thickness of the title
  h8->GetXaxis()->SetLabelSize(0.03);
  h8->GetXaxis()->SetLabelFont(72);

  if (DoSave)
  {
    canvas->SaveAs(saveFileName1);
    canvas2->SaveAs(saveFileName2);
    canvas3->SaveAs(saveFileName3);
  }
}
