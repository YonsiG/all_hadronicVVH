#define _EFFICIENCY_C_
#include "efficiency.h"
#include "tools.h"
#include "EventSelection.h"

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
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "Math/VectorUtil.h"
#include "TH3F.h"
#include <TRandom3.h>
#include "TMinuit.h"
#include <TApplication.h>
#include <TEnv.h>
#include <TComplex.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TTree.h>

using namespace std;
void efficiency::Loop(const char *typeName)
{
   if (isRead == false || fChain == 0)
   {
      cout << "This file contains no events, skip to next file" << endl;
      return;
   }

   /********defining cross sections**********/
   double xsection;
   xsection = getXsection(string(typeName));

   /********defining total weights***********/
   int weightnum = runChain->GetEntries();
   for (int iweight = 0; iweight < weightnum; iweight++)
   {
      runChain->GetEntry(iweight);
      myHists->weight_Scale->Fill(0.5, genEventSumw);
   }

   Int_t Nentries = fChain->GetEntries();

   /*****************************************/
   /*********Main Looping Code Start*********/
   /*****************************************/

   for (int iLoop = 0; iLoop < Nentries; iLoop++)
   {
      fChain->GetEntry(iLoop);
      Sta_TotalNumber++;
      Sta_FileEventNumber++;

      double weight = double(genWeight * xsection * 137.0);
      for (int icutflow = 0; icutflow < 14; icutflow++)
         myHists->cutflow[icutflow]->Fill(0.5, weight);

      /*************gen particles****************/
      if (typeName == string("WZH") || typeName == string("ZZH") || typeName == string("OSWWH") || typeName == string("SSWWH"))
      {
         int Is_Hbb = 0;
         GenPart_Filter(Is_Hbb, GenHiggs, GenBquarkFromH, GenantiBquarkFromH,
                        GenVBFJets[0], GenVBFJets[1],
                        nGenPart, GenPart_status, GenPart_pdgId, GenPart_genPartIdxMother,
                        GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass);

         if (Is_Hbb != 1)
            continue; // only consider Generated Hbb events, other inclusive channels are not studied
      }

      for (int icutflow = 0; icutflow < 14; icutflow++)
         myHists->cutflow[icutflow]->Fill(1.5, weight);

      /*******************trigger********************/
      bool trig_status = 0;
      if (HLT_AK8PFHT800_TrimMass50)
         trig_status = 1;
      if (HLT_PFHT1050)
         trig_status = 1;
      if (HLT_PFJet500)
         trig_status = 1;
      if (HLT_AK8PFJet500)
         trig_status = 1;
      if (HLT_AK8PFJet400_TrimMass30)
         trig_status = 1;
      if (HLT_AK8PFJet420_TrimMass30)
         trig_status = 1;

      if (trig_status == 0)
         continue;
      for (int icutflow = 0; icutflow < 14; icutflow++)
         myHists->cutflow[icutflow]->Fill(2.5, weight);

      /**************lepton selection****************/
      int count_lepton = 0;

      for (int ielec = 0; ielec < nElectron; ielec++)
      {
         Electron[ielec].SetPtEtaPhiM(Electron_pt[ielec], Electron_eta[ielec], Electron_phi[ielec], Electron_mass[ielec]);
         if (Electron_LooseId_ttH(Electron[ielec], Electron_deltaEtaSC[ielec], Electron_dxy[ielec], Electron_dz[ielec], Electron_sip3d[ielec], Electron_miniPFRelIso_all[ielec], Electron_lostHits[ielec], Electron_mvaFall17V2noIso_WPL[ielec]))
            count_lepton++;
      }

      for (int imuon = 0; imuon < nMuon; imuon++)
      {
         Muon[imuon].SetPtEtaPhiM(Muon_pt[imuon], Muon_eta[imuon], Muon_phi[imuon], Muon_mass[imuon]);
         if (Muon_LooseId_ttH(Muon[imuon], Muon_dxy[imuon], Muon_dz[imuon], Muon_sip3d[imuon], Muon_miniPFRelIso_all[imuon], Muon_looseId[imuon]))
            count_lepton++;
      }

      if (count_lepton != 0)
         continue;

      for (int icutflow = 0; icutflow < 14; icutflow++)
         myHists->cutflow[icutflow]->Fill(3.5, weight);

      /****************fatjet selection**************/
      TLorentzVector tempFatJet;
      int sort_index = 0;
      double pass_cut_fatjet_msoftdrop[200];
      double pass_cut_fatjet_WvsQCD[200];
      double pass_cut_fatjet_mass[200];
      double pass_cut_fatjet_Xbb_modified[200];
      double pass_cut_fatjet_Xcc[200];
      double pass_cut_fatjet_Xqq[200];
      double pass_cut_fatjet_QCD[200];
      double pass_cut_fatjet_Xccqq_modified[200];
      int count_fatjet = 0;
      for (int ifatjet = 0; ifatjet < nFatJet; ifatjet++)
      {
         tempFatJet.SetPtEtaPhiM(FatJet_pt[ifatjet], FatJet_eta[ifatjet], FatJet_phi[ifatjet], FatJet_mass[ifatjet]);
         bool Fatjet_kinematic_filter = FatJet_kinematic_Select(tempFatJet, FatJet_jetId[ifatjet], FatJet_msoftdrop[ifatjet]);
         if (Fatjet_kinematic_filter == false)
            continue;
         count_fatjet++;
         FatJet[sort_index] = tempFatJet;
         pass_cut_fatjet_msoftdrop[sort_index] = FatJet_msoftdrop[ifatjet];
         pass_cut_fatjet_WvsQCD[sort_index] = FatJet_particleNet_WvsQCD[ifatjet];
         pass_cut_fatjet_mass[sort_index] = FatJet_particleNet_mass[ifatjet];
         pass_cut_fatjet_Xbb_modified[sort_index] = FatJet_particleNetMD_Xbb[ifatjet]/(FatJet_particleNetMD_QCD[ifatjet] + FatJet_particleNetMD_Xbb[ifatjet]);
         pass_cut_fatjet_Xcc[sort_index] = FatJet_particleNetMD_Xcc[ifatjet];
         pass_cut_fatjet_Xqq[sort_index] = FatJet_particleNetMD_Xqq[ifatjet];
         pass_cut_fatjet_QCD[sort_index] = FatJet_particleNetMD_QCD[ifatjet];
         pass_cut_fatjet_Xccqq_modified[sort_index] = (FatJet_particleNetMD_Xcc[ifatjet] + FatJet_particleNetMD_Xqq[ifatjet]) / (FatJet_particleNetMD_QCD[ifatjet] + FatJet_particleNetMD_Xcc[ifatjet] + FatJet_particleNetMD_Xqq[ifatjet]);
         
         sort_index++;
      }

      // sorting the fatjets according to their btag score
      for (int isort = 0; isort < sort_index; isort++)
      {
         int thisfatjet_rank = 0;
         for (int irankindex = 0; irankindex < sort_index; irankindex++)
         {
            if (isort == irankindex)
               continue;
            if (pass_cut_fatjet_Xbb_modified[isort] < pass_cut_fatjet_Xbb_modified[irankindex])
               thisfatjet_rank++;
         }
         FatJet_btagsort[thisfatjet_rank] = FatJet[isort];
         FatJet_msoftdrop_btagsort[thisfatjet_rank] = pass_cut_fatjet_msoftdrop[isort];
         FatJet_WvsQCD_btagsort[thisfatjet_rank] = pass_cut_fatjet_WvsQCD[isort];
         FatJet_mass_btagsort[thisfatjet_rank] = pass_cut_fatjet_mass[isort];
         FatJet_Xbb_modified_btagsort[thisfatjet_rank] = pass_cut_fatjet_Xbb_modified[isort];
         FatJet_Xcc_btagsort[thisfatjet_rank] = pass_cut_fatjet_Xcc[isort];
         FatJet_Xqq_btagsort[thisfatjet_rank] = pass_cut_fatjet_Xqq[isort];
         FatJet_QCD_btagsort[thisfatjet_rank] = pass_cut_fatjet_QCD[isort];
         FatJet_Xccqq_modified_btagsort[thisfatjet_rank] = pass_cut_fatjet_Xccqq_modified[isort];
      }

      /*****************jet selection****************/
      // pass the jet_kinematic_selection first
      int count_jet = 0;
      int count_central_jet = 0;
      TLorentzVector tempJet;
      for (int ijet = 0; ijet < nJet; ijet++)
      {
         tempJet.SetPtEtaPhiM(Jet_pt[ijet], Jet_eta[ijet], Jet_phi[ijet], Jet_mass[ijet]);
         if (!Jet_kinVBF_select(tempJet))
            continue;

         int count_overlap = 0;
         for (int ifatjet = 0; ifatjet < count_fatjet; ifatjet++)
            if (tempJet.DeltaR(FatJet[ifatjet]) < 0.8)
               count_overlap = 1;

         if (count_overlap == 1)
            continue;

         Jet[count_jet] = tempJet;
         count_jet++;
         if (Jet_central_Select(tempJet))
            count_central_jet++;
      }

      /*****VBF selection*****/
      bool VBF_selection = 0;
      int VBF_method = 2;
      int VBF_jet_index[2];
      float VBF_max_mass = 0;
      float VBF_max_DeltaEta = 0;

      if (count_jet >= 2)
      {
         VBF_selection = VBF_Selection(VBF_method, VBF_jet_index[0], VBF_jet_index[1],
                                       VBF_max_mass, VBF_max_DeltaEta, count_jet, Jet);
      }
      
      /*****************categorization****************/
      if (count_fatjet >= 3)
      {
         myHists->cutflow[0]->Fill(4.5, weight);
         myHists->cutflow[1]->Fill(4.5, weight);
         myHists->cutflow[2]->Fill(4.5, weight);
      }

      if (count_fatjet == 2)
      {
         myHists->cutflow[3]->Fill(4.5, weight);
         myHists->cutflow[4]->Fill(4.5, weight);
         myHists->cutflow[5]->Fill(4.5, weight);
         myHists->cutflow[6]->Fill(4.5, weight);
      }

      if (count_fatjet == 1)
      {
         myHists->cutflow[7]->Fill(4.5, weight);
         myHists->cutflow[8]->Fill(4.5, weight);
         myHists->cutflow[9]->Fill(4.5, weight);
         myHists->cutflow[10]->Fill(4.5, weight);
         myHists->cutflow[11]->Fill(4.5, weight);
         myHists->cutflow[12]->Fill(4.5, weight);
      }

      int category_number = categorize(count_fatjet, count_jet);
      myHists->cutflow[category_number]->Fill(5.5, weight);

      if (VBF_selection != 1)
         continue;
      myHists->cutflow[category_number]->Fill(6.5, weight);

      for (int ijet = 0; ijet < count_jet; ijet++)
      {
         bool isVBFJet = false;
         if (ijet == VBF_jet_index[0])
            isVBFJet = true;
         if (ijet == VBF_jet_index[1])
            isVBFJet = true;
      }

      /****************plot filling******************/
      myHists->number_of_fatjets[category_number]->Fill(count_fatjet, weight);
      myHists->number_of_jets[category_number]->Fill(count_jet, weight);
      myHists->number_of_central_jets[category_number]->Fill(count_central_jet, weight);

      myHists->fatjet_msoftdrop[category_number][0]->Fill(FatJet_msoftdrop_btagsort[0],weight);
      myHists->fatjet_pt[category_number][0]->Fill(FatJet_btagsort[0].Pt(),weight);
      myHists->fatjet_eta[category_number][0]->Fill(FatJet_btagsort[0].Eta(),weight);
      myHists->fatjet_WvsQCD[category_number][0]->Fill(FatJet_WvsQCD_btagsort[0],weight);
      myHists->fatjet_mass[category_number][0]->Fill(FatJet_mass_btagsort[0],weight);
      myHists->fatjet_Xbb_modified[category_number][0]->Fill(FatJet_Xbb_modified_btagsort[0],weight);
      myHists->fatjet_Xcc[category_number][0]->Fill(FatJet_Xcc_btagsort[0],weight);
      myHists->fatjet_Xqq[category_number][0]->Fill(FatJet_Xqq_btagsort[0],weight);
      myHists->fatjet_QCD[category_number][0]->Fill(FatJet_QCD_btagsort[0],weight);
      myHists->fatjet_Xccqq_modified[category_number][0]->Fill(FatJet_Xccqq_modified_btagsort[0],weight);

      if (sort_index>1){
         myHists->fatjet_msoftdrop[category_number][1]->Fill(FatJet_msoftdrop_btagsort[1],weight);
         myHists->fatjet_pt[category_number][1]->Fill(FatJet_btagsort[1].Pt(),weight);
         myHists->fatjet_eta[category_number][1]->Fill(FatJet_btagsort[1].Eta(),weight);
         myHists->fatjet_WvsQCD[category_number][1]->Fill(FatJet_WvsQCD_btagsort[1],weight);
         myHists->fatjet_mass[category_number][1]->Fill(FatJet_mass_btagsort[1],weight);
         myHists->fatjet_Xbb_modified[category_number][1]->Fill(FatJet_Xbb_modified_btagsort[1],weight);
         myHists->fatjet_Xcc[category_number][1]->Fill(FatJet_Xcc_btagsort[1],weight);
         myHists->fatjet_Xqq[category_number][1]->Fill(FatJet_Xqq_btagsort[1],weight);
         myHists->fatjet_QCD[category_number][1]->Fill(FatJet_QCD_btagsort[1],weight);
         myHists->fatjet_Xccqq_modified[category_number][1]->Fill(FatJet_Xccqq_modified_btagsort[1],weight);
      }

      if (sort_index>2){
         myHists->fatjet_msoftdrop[category_number][2]->Fill(FatJet_msoftdrop_btagsort[2],weight);
         myHists->fatjet_pt[category_number][2]->Fill(FatJet_btagsort[2].Pt(),weight);
         myHists->fatjet_eta[category_number][2]->Fill(FatJet_btagsort[2].Eta(),weight);
         myHists->fatjet_WvsQCD[category_number][2]->Fill(FatJet_WvsQCD_btagsort[2],weight);
         myHists->fatjet_mass[category_number][2]->Fill(FatJet_mass_btagsort[2],weight);
         myHists->fatjet_Xbb_modified[category_number][2]->Fill(FatJet_Xbb_modified_btagsort[2],weight);
         myHists->fatjet_Xcc[category_number][2]->Fill(FatJet_Xcc_btagsort[2],weight);
         myHists->fatjet_Xqq[category_number][2]->Fill(FatJet_Xqq_btagsort[2],weight);
         myHists->fatjet_QCD[category_number][2]->Fill(FatJet_QCD_btagsort[2],weight);
         myHists->fatjet_Xccqq_modified[category_number][2]->Fill(FatJet_Xccqq_modified_btagsort[2],weight);
      }

      myHists->VBF_max_mass[category_number]->Fill(VBF_max_mass,weight);

      for (int isort = 0; isort < sort_index; isort++)
         myHists->fatjet_btag_score->Fill(FatJet_Xbb_modified_btagsort[isort], weight);
      if (sort_index > 0)
         myHists->first_fatjet_btag_score->Fill(FatJet_Xbb_modified_btagsort[0], weight);
      if (sort_index > 1)
         myHists->second_fatjet_btag_score->Fill(FatJet_Xbb_modified_btagsort[1], weight);
      if (sort_index > 2)
         myHists->third_fatjet_btag_score->Fill(FatJet_Xbb_modified_btagsort[2], weight);

      if (sort_index > 0)
      {
         if (FatJet_btagsort[0].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[0].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->first_fatjet_btag_score_2match->Fill(FatJet_Xbb_modified_btagsort[0], weight);
         if (FatJet_btagsort[0].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[0].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->first_fatjet_btag_score_1match->Fill(FatJet_Xbb_modified_btagsort[0], weight);
         if (FatJet_btagsort[0].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[0].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->first_fatjet_btag_score_1match->Fill(FatJet_Xbb_modified_btagsort[0], weight);
         if (FatJet_btagsort[0].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[0].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->first_fatjet_btag_score_0match->Fill(FatJet_Xbb_modified_btagsort[0], weight);
      }

      if (sort_index > 1)
      {
         if (FatJet_btagsort[1].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[1].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->second_fatjet_btag_score_2match->Fill(FatJet_Xbb_modified_btagsort[1], weight);
         if (FatJet_btagsort[1].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[1].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->second_fatjet_btag_score_1match->Fill(FatJet_Xbb_modified_btagsort[1], weight);
         if (FatJet_btagsort[1].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[1].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->second_fatjet_btag_score_1match->Fill(FatJet_Xbb_modified_btagsort[1], weight);
         if (FatJet_btagsort[1].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[1].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->second_fatjet_btag_score_0match->Fill(FatJet_Xbb_modified_btagsort[1], weight);
      }

      if (sort_index > 2)
      {
         if (FatJet_btagsort[2].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[2].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->third_fatjet_btag_score_2match->Fill(FatJet_Xbb_modified_btagsort[2], weight);
         if (FatJet_btagsort[2].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[2].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->third_fatjet_btag_score_1match->Fill(FatJet_Xbb_modified_btagsort[2], weight);
         if (FatJet_btagsort[2].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[2].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->third_fatjet_btag_score_1match->Fill(FatJet_Xbb_modified_btagsort[2], weight);
         if (FatJet_btagsort[2].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[2].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->third_fatjet_btag_score_0match->Fill(FatJet_Xbb_modified_btagsort[2], weight);
      }

      if (VBF_selection == 1)
      {
         double distance00 = Jet[VBF_jet_index[0]].DeltaR(GenVBFJets[0]);
         double distance01 = Jet[VBF_jet_index[0]].DeltaR(GenVBFJets[1]);
         int condition;
         if (distance00 < distance01)
            condition = 1; // VBF0&Gen0 VBF1&Gen1
         if (distance00 > distance01)
            condition = 2; // VBF0&Gen1 VBF1&Gen0

         double distance10 = Jet[VBF_jet_index[1]].DeltaR(GenVBFJets[0]);
         double distance11 = Jet[VBF_jet_index[1]].DeltaR(GenVBFJets[1]);

         if (condition == 1 && distance00 < 0.4 && distance11 < 0.4)
         {
            myHists->VBFJet_leadingPt_2match->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->VBFJet_subleadingPt_2match->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->VBFJet_Mjj_2match->Fill(VBF_max_mass, weight);
            myHists->VBFJet_DeltaEta_2match_total->Fill(VBF_max_DeltaEta, weight);
            myHists->VBFJet_DeltaEta_2match[category_number]->Fill(VBF_max_DeltaEta, weight);

            myHists->matched_VBFJet_Pt[category_number]->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->matched_VBFJet_Pt[category_number]->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->matched_VBFJet_Pt_total->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->matched_VBFJet_Pt_total->Fill(Jet[VBF_jet_index[1]].Pt(), weight);

            myHists->matched_VBFJet_Eta[category_number]->Fill(Jet[VBF_jet_index[0]].Eta(), weight);
            myHists->matched_VBFJet_Eta[category_number]->Fill(Jet[VBF_jet_index[1]].Eta(), weight);
            myHists->matched_VBFJet_Eta_total->Fill(Jet[VBF_jet_index[0]].Eta(), weight);
            myHists->matched_VBFJet_Eta_total->Fill(Jet[VBF_jet_index[1]].Eta(), weight);

            myHists->matched_VBFJet_qgl[category_number]->Fill(Jet_qgl[VBF_jet_index[0]], weight);
            myHists->matched_VBFJet_qgl[category_number]->Fill(Jet_qgl[VBF_jet_index[1]], weight);
            myHists->matched_VBFJet_qgl_total->Fill(Jet_qgl[VBF_jet_index[0]], weight);
            myHists->matched_VBFJet_qgl_total->Fill(Jet_qgl[VBF_jet_index[1]], weight);
         }
         if (condition == 2 && distance01 < 0.4 && distance10 < 0.4)
         {
            myHists->VBFJet_leadingPt_2match->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->VBFJet_subleadingPt_2match->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->VBFJet_Mjj_2match->Fill(VBF_max_mass, weight);
            myHists->VBFJet_DeltaEta_2match_total->Fill(VBF_max_DeltaEta, weight);
            myHists->VBFJet_DeltaEta_2match[category_number]->Fill(VBF_max_DeltaEta, weight);

            myHists->matched_VBFJet_Pt[category_number]->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->matched_VBFJet_Pt[category_number]->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->matched_VBFJet_Pt_total->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->matched_VBFJet_Pt_total->Fill(Jet[VBF_jet_index[1]].Pt(), weight);

            myHists->matched_VBFJet_Eta[category_number]->Fill(Jet[VBF_jet_index[0]].Eta(), weight);
            myHists->matched_VBFJet_Eta[category_number]->Fill(Jet[VBF_jet_index[1]].Eta(), weight);
            myHists->matched_VBFJet_Eta_total->Fill(Jet[VBF_jet_index[0]].Eta(), weight);
            myHists->matched_VBFJet_Eta_total->Fill(Jet[VBF_jet_index[1]].Eta(), weight);

            myHists->matched_VBFJet_qgl[category_number]->Fill(Jet_qgl[VBF_jet_index[0]], weight);
            myHists->matched_VBFJet_qgl[category_number]->Fill(Jet_qgl[VBF_jet_index[1]], weight);
            myHists->matched_VBFJet_qgl_total->Fill(Jet_qgl[VBF_jet_index[0]], weight);
            myHists->matched_VBFJet_qgl_total->Fill(Jet_qgl[VBF_jet_index[1]], weight);
         }

         if (condition == 1 && distance00 < 0.4 && distance11 > 0.4)
         {
            myHists->VBFJet_leadingPt_1match->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->VBFJet_subleadingPt_1match->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->VBFJet_Mjj_1match->Fill(VBF_max_mass, weight);
            myHists->VBFJet_DeltaEta_1match_total->Fill(VBF_max_DeltaEta, weight);
            myHists->VBFJet_DeltaEta_1match[category_number]->Fill(VBF_max_DeltaEta, weight);

            myHists->matched_VBFJet_Pt[category_number]->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->matched_VBFJet_Pt_total->Fill(Jet[VBF_jet_index[0]].Pt(), weight);

            myHists->matched_VBFJet_Eta[category_number]->Fill(Jet[VBF_jet_index[0]].Eta(), weight);
            myHists->matched_VBFJet_Eta_total->Fill(Jet[VBF_jet_index[0]].Eta(), weight);

            myHists->matched_VBFJet_qgl[category_number]->Fill(Jet_qgl[VBF_jet_index[0]], weight);
            myHists->matched_VBFJet_qgl_total->Fill(Jet_qgl[VBF_jet_index[0]], weight);
         }

         if (condition == 1 && distance00 > 0.4 && distance11 < 0.4)
         {
            myHists->VBFJet_leadingPt_1match->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->VBFJet_subleadingPt_1match->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->VBFJet_Mjj_1match->Fill(VBF_max_mass, weight);
            myHists->VBFJet_DeltaEta_1match_total->Fill(VBF_max_DeltaEta, weight);
            myHists->VBFJet_DeltaEta_1match[category_number]->Fill(VBF_max_DeltaEta, weight);

            myHists->matched_VBFJet_Pt[category_number]->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->matched_VBFJet_Pt_total->Fill(Jet[VBF_jet_index[1]].Pt(), weight);

            myHists->matched_VBFJet_Eta[category_number]->Fill(Jet[VBF_jet_index[1]].Eta(), weight);
            myHists->matched_VBFJet_Eta_total->Fill(Jet[VBF_jet_index[1]].Eta(), weight);

            myHists->matched_VBFJet_qgl[category_number]->Fill(Jet_qgl[VBF_jet_index[1]], weight);
            myHists->matched_VBFJet_qgl_total->Fill(Jet_qgl[VBF_jet_index[1]], weight);
         }

         if (condition == 2 && distance01 < 0.4 && distance10 > 0.4)
         {
            myHists->VBFJet_leadingPt_1match->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->VBFJet_subleadingPt_1match->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->VBFJet_Mjj_1match->Fill(VBF_max_mass, weight);
            myHists->VBFJet_DeltaEta_1match_total->Fill(VBF_max_DeltaEta, weight);
            myHists->VBFJet_DeltaEta_1match[category_number]->Fill(VBF_max_DeltaEta, weight);

            myHists->matched_VBFJet_Pt[category_number]->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->matched_VBFJet_Pt_total->Fill(Jet[VBF_jet_index[0]].Pt(), weight);

            myHists->matched_VBFJet_Eta[category_number]->Fill(Jet[VBF_jet_index[0]].Eta(), weight);
            myHists->matched_VBFJet_Eta_total->Fill(Jet[VBF_jet_index[0]].Eta(), weight);

            myHists->matched_VBFJet_qgl[category_number]->Fill(Jet_qgl[VBF_jet_index[0]], weight);
            myHists->matched_VBFJet_qgl_total->Fill(Jet_qgl[VBF_jet_index[0]], weight);
         }
         if (condition == 2 && distance01 > 0.4 && distance10 < 0.4)
         {
            myHists->VBFJet_leadingPt_1match->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->VBFJet_subleadingPt_1match->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->VBFJet_Mjj_1match->Fill(VBF_max_mass, weight);
            myHists->VBFJet_DeltaEta_1match_total->Fill(VBF_max_DeltaEta, weight);
            myHists->VBFJet_DeltaEta_1match[category_number]->Fill(VBF_max_DeltaEta, weight);

            myHists->matched_VBFJet_Pt[category_number]->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->matched_VBFJet_Pt_total->Fill(Jet[VBF_jet_index[1]].Pt(), weight);

            myHists->matched_VBFJet_Eta[category_number]->Fill(Jet[VBF_jet_index[1]].Eta(), weight);
            myHists->matched_VBFJet_Eta_total->Fill(Jet[VBF_jet_index[1]].Eta(), weight);

            myHists->matched_VBFJet_qgl[category_number]->Fill(Jet_qgl[VBF_jet_index[1]], weight);
            myHists->matched_VBFJet_qgl_total->Fill(Jet_qgl[VBF_jet_index[1]], weight);
         }

         if (condition == 1 && distance00 > 0.4 && distance11 > 0.4)
         {
            myHists->VBFJet_leadingPt_0match->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->VBFJet_subleadingPt_0match->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->VBFJet_Mjj_0match->Fill(VBF_max_mass, weight);
            myHists->VBFJet_DeltaEta_0match_total->Fill(VBF_max_DeltaEta, weight);
            myHists->VBFJet_DeltaEta_0match[category_number]->Fill(VBF_max_DeltaEta, weight);
         }

         if (condition == 2 && distance01 > 0.4 && distance10 > 0.4)
         {
            myHists->VBFJet_leadingPt_0match->Fill(Jet[VBF_jet_index[0]].Pt(), weight);
            myHists->VBFJet_subleadingPt_0match->Fill(Jet[VBF_jet_index[1]].Pt(), weight);
            myHists->VBFJet_Mjj_0match->Fill(VBF_max_mass, weight);
            myHists->VBFJet_DeltaEta_0match_total->Fill(VBF_max_DeltaEta, weight);
            myHists->VBFJet_DeltaEta_0match[category_number]->Fill(VBF_max_DeltaEta, weight);
         }
      }

      /*
      if (FatJet_btagXbb_modified[maxindex] > 0.7)
      {

      }

      if (FatJet_btagXbb_modified[maxindex] > 0.89)
      {

      }

      if (FatJet_btagXbb_modified[maxindex] > 0.92)
      {

      }
*/
      /************clearing variables************/

   } // main looping

} // efficiency::loop ends
