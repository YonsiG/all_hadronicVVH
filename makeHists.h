#ifndef _MAKEHIST_H_
#define _MAKEHIST_H_

#include <iostream>
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"
#include "TLegend.h"

class makeHists
{
public:
  TFile *hf;

  TH1D *weight_Scale;
  TH1D *cutflow_1;
  TH1D *cutflow_2;
  TH1D *cutflow_3;
  TH1D *cutflow_4;
  TH1D *cutflow_5;
  TH1D *cutflow_6;
  TH1D *cutflow_7;
  TH1D *cutflow_8;
  TH1D *cutflow_9;
  TH1D *cutflow_10;
  TH1D *cutflow_11;
  TH1D *cutflow_12;
  TH1D *cutflow_13;
  TH1D *cutflow_14;
  TH1D *cutflow_15;
  TH1D *cutflow_16;
  TH1D *cutflow_17;

  void createHists(const char *fileName);
  void saveHists();
};
#endif
