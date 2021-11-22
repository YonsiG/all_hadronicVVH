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
  TH1D *cutflow[14];

  TH1D *number_of_jets[14];
  TH1D *number_of_central_jets[14];

  TH1D *fatjet_btag_score;
  TH1D *first_fatjet_btag_score;
  TH1D *second_fatjet_btag_score;
  TH1D *third_fatjet_btag_score;

  TH1D *first_fatjet_btag_score_2match;
  TH1D *second_fatjet_btag_score_2match;
  TH1D *third_fatjet_btag_score_2match;

  TH1D *first_fatjet_btag_score_1match;
  TH1D *second_fatjet_btag_score_1match;
  TH1D *third_fatjet_btag_score_1match;

  TH1D *first_fatjet_btag_score_0match;
  TH1D *second_fatjet_btag_score_0match;
  TH1D *third_fatjet_btag_score_0match;

  TH1D *VBFJet_leadingPt_2match;
  TH1D *VBFJet_leadingPt_1match;
  TH1D *VBFJet_leadingPt_0match;
  TH1D *VBFJet_subleadingPt_2match;
  TH1D *VBFJet_subleadingPt_1match;
  TH1D *VBFJet_subleadingPt_0match;

  TH1D *VBFJet_leadingMass_2match;
  TH1D *VBFJet_leadingMass_1match;
  TH1D *VBFJet_leadingMass_0match;
  TH1D *VBFJet_subleadingMass_2match;
  TH1D *VBFJet_subleadingMass_1match;
  TH1D *VBFJet_subleadingMass_0match;

  TH1D *VBFJet_Mjj_2match;
  TH1D *VBFJet_Mjj_1match;
  TH1D *VBFJet_Mjj_0match;

  void createHists(const char *fileName);
  void saveHists();
};
#endif
