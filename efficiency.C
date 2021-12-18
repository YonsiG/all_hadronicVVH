#define _EFFICIENCY_C_
#include "efficiency.h"
#include "tools.h"
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
   if (string(typeName) == string("OSWWH"))
      xsection = 5.652;
   if (string(typeName) == string("SSWWH"))
      xsection = 3.559;
   if (string(typeName) == string("WZH"))
      xsection = 3.742;
   if (string(typeName) == string("ZZH"))
      xsection = 2.994;

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

      int double_count_Hbb = 0;

      for (int ipart = 0; ipart < nGenPart; ipart++)
      {
         if (GenPart_status[ipart] == 62 && GenPart_pdgId[ipart] == 25)
            GenHiggs.SetPtEtaPhiM(GenPart_pt[ipart], GenPart_eta[ipart], GenPart_phi[ipart], GenPart_mass[ipart]);

         if (GenPart_pdgId[ipart] == -5 || GenPart_pdgId[ipart] == 5)
         {
            if (GenPart_pdgId[GenPart_genPartIdxMother[ipart]] != 25 || GenPart_status[GenPart_genPartIdxMother[ipart]] != 62)
               continue;

            if (GenPart_pdgId[ipart] == 5)
               GenBquarkFromH.SetPtEtaPhiM(GenPart_pt[ipart], GenPart_eta[ipart], GenPart_phi[ipart], GenPart_mass[ipart]);

            if (GenPart_pdgId[ipart] == -5)
               GenantiBquarkFromH.SetPtEtaPhiM(GenPart_pt[ipart], GenPart_eta[ipart], GenPart_phi[ipart], GenPart_mass[ipart]);

            double_count_Hbb++; // each Hbb event will count 2 in b, so called double count
         }
      }

      // keep track of the final W(pdgId=+-24)(status=62) and their daughters(2 quarks or 2 leptons/neutrino)

      // keep track of VBF quarks(pdgId=+-1,2,3,4,5,6), the mothers should be quarks from protons(search for pdgId), their daughters should be W and quark
      // easy:Instance 2,3 are always the VBF jets. GenPart[2] and GenPart[3]
      GenVBFJets[0].SetPtEtaPhiM(GenPart_pt[2], GenPart_eta[2], GenPart_phi[2], GenPart_mass[2]);
      GenVBFJets[1].SetPtEtaPhiM(GenPart_pt[3], GenPart_eta[3], GenPart_phi[3], GenPart_mass[3]);

      if (double_count_Hbb != 2)
         continue; // only consider Generated Hbb events, other inclusive channels are not studied

      for (int icutflow = 0; icutflow < 14; icutflow++)
         myHists->cutflow[icutflow]->Fill(1.5, weight);

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
         myHists->cutflow[icutflow]->Fill(2.5, weight);

      /****************fatjet selection**************/
      TLorentzVector tempFatJet;
      int sort_index = 0;
      double pass_cut_fatjet_btagDDBvL[200]; // temporarily not sorting
      int count_fatjet = 0;
      for (int ifatjet = 0; ifatjet < nFatJet; ifatjet++)
      {
         tempFatJet.SetPtEtaPhiM(FatJet_pt[ifatjet], FatJet_eta[ifatjet], FatJet_phi[ifatjet], FatJet_mass[ifatjet]);
         bool Fatjet_kinematic_filter = FatJet_kinematic_Select(tempFatJet, FatJet_jetId[ifatjet], FatJet_msoftdrop[ifatjet]);
         if (Fatjet_kinematic_filter == false)
            continue;
         count_fatjet++;
         FatJet[sort_index] = tempFatJet;
         pass_cut_fatjet_btagDDBvL[sort_index] = FatJet_btagDDBvL[ifatjet];
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
            if (pass_cut_fatjet_btagDDBvL[isort] < pass_cut_fatjet_btagDDBvL[irankindex])
               thisfatjet_rank++;
         }
         FatJet_btagsort[thisfatjet_rank] = FatJet[isort];
         FatJet_DDBvL_btagsort[thisfatjet_rank] = pass_cut_fatjet_btagDDBvL[isort];
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
      int VBF_selection = 0;
      int VBF_jet_index[2];
      float VBF_max_mass = 0;
      float VBF_max_DeltaEta = 0;
      int VBF_method = 1;

      if (count_jet >= 2)
      {
         VBF_max_mass = 0;
         if (VBF_method == 1)
         {
            for (int ivbf = 0; ivbf < count_jet; ivbf++)
            {
               for (int ivbf2 = ivbf; ivbf2 < count_jet; ivbf2++)
               {
                  if ((Jet[ivbf] + Jet[ivbf2]).M() > VBF_max_mass)
                  {
                     if (Jet[ivbf].Pt() > Jet[ivbf2].Pt())
                     {
                        VBF_jet_index[0] = ivbf;
                        VBF_jet_index[1] = ivbf2;
                     }
                     if (Jet[ivbf].Pt() < Jet[ivbf2].Pt())
                     {
                        VBF_jet_index[0] = ivbf2;
                        VBF_jet_index[1] = ivbf;
                     }
                     VBF_max_mass = (Jet[ivbf] + Jet[ivbf2]).M();
                  }
               }
            }
         }
         VBF_max_mass = 0;

         VBF_max_DeltaEta = 0;
         if (VBF_method == 2)
         {
            for (int ivbf = 0; ivbf < count_jet; ivbf++)
            {
               for (int ivbf2 = ivbf; ivbf2 < count_jet; ivbf2++)
               {
                  if (fabs(Jet[ivbf].Eta() - Jet[ivbf2].Eta()) > VBF_max_DeltaEta)
                  {
                     if (Jet[ivbf].Pt() > Jet[ivbf2].Pt())
                     {
                        VBF_jet_index[0] = ivbf;
                        VBF_jet_index[1] = ivbf2;
                     }
                     if (Jet[ivbf].Pt() < Jet[ivbf2].Pt())
                     {
                        VBF_jet_index[0] = ivbf2;
                        VBF_jet_index[1] = ivbf;
                     }
                     VBF_max_DeltaEta = fabs(Jet[ivbf].Eta() - Jet[ivbf2].Eta());
                  }
               }
            }
         }
         VBF_max_DeltaEta = 0;

         if (VBF_method == 3)
         {
            float VBF_max_Energy = 0;
            for (int ivbf = 0; ivbf < count_jet; ivbf++)
            {
               if (Jet[ivbf].E() > VBF_max_Energy)
               {
                  VBF_jet_index[0] = ivbf;
                  VBF_max_Energy = Jet[ivbf].E();
               }
            } // select max energy jet
            VBF_max_Energy = 0;
            int opposite_eta = 0;
            for (int ivbf = 0; ivbf < count_jet; ivbf++)
            {
               if (ivbf == VBF_jet_index[0])
                  continue;
               if (Jet[ivbf].Eta() * Jet[VBF_jet_index[0]].Eta() > 0)
                  continue;
               opposite_eta = 1;
               if (Jet[ivbf].E() > VBF_max_Energy)
               {
                  VBF_jet_index[1] = ivbf;
                  VBF_max_Energy = Jet[ivbf].E();
               }
            }

            VBF_max_DeltaEta = 0;
            if (opposite_eta == 0)
            {
               for (int ivbf = 0; ivbf < count_jet; ivbf++)
               {
                  if (ivbf == VBF_jet_index[0])
                     continue;
                  if (fabs(Jet[ivbf].Eta() - Jet[VBF_jet_index[0]].Eta()) > VBF_max_DeltaEta)
                  {
                     VBF_jet_index[1] = ivbf;
                     VBF_max_DeltaEta = fabs(Jet[ivbf].Eta() - Jet[VBF_jet_index[0]].Eta());
                  }
               }
            }
         }

         VBF_max_mass = (Jet[VBF_jet_index[0]] + Jet[VBF_jet_index[1]]).M();
         VBF_max_DeltaEta = fabs(Jet[VBF_jet_index[0]].Eta() - Jet[VBF_jet_index[1]].Eta());
         // if (VBF_max_mass > 500 && VBF_max_DeltaEta > 3)
         if (VBF_max_mass > 500)
            VBF_selection = 1;
      }

      /*****************categorization****************/
      int category_number = 13;
      if (count_fatjet >= 3)
      {
         myHists->cutflow[0]->Fill(3.5, weight);
         myHists->cutflow[1]->Fill(3.5, weight);
         myHists->cutflow[2]->Fill(3.5, weight);
      }
      if (count_fatjet >= 3 && count_jet >= 2)
         category_number = 0;
      if (count_fatjet >= 3 && count_jet == 1)
         category_number = 1;
      if (count_fatjet >= 3 && count_jet < 1)
         category_number = 2;

      if (count_fatjet == 2)
      {
         myHists->cutflow[3]->Fill(3.5, weight);
         myHists->cutflow[4]->Fill(3.5, weight);
         myHists->cutflow[5]->Fill(3.5, weight);
         myHists->cutflow[6]->Fill(3.5, weight);
      }
      if (count_fatjet == 2 && count_jet >= 4)
         category_number = 3;
      if (count_fatjet == 2 && count_jet == 3)
         category_number = 4;
      if (count_fatjet == 2 && count_jet == 2)
         category_number = 5;
      if (count_fatjet == 2 && count_jet < 2)
         category_number = 6;

      if (count_fatjet == 1)
      {
         myHists->cutflow[7]->Fill(3.5, weight);
         myHists->cutflow[8]->Fill(3.5, weight);
         myHists->cutflow[9]->Fill(3.5, weight);
         myHists->cutflow[10]->Fill(3.5, weight);
         myHists->cutflow[11]->Fill(3.5, weight);
         myHists->cutflow[12]->Fill(3.5, weight);
      }
      if (count_fatjet == 1 && count_jet >= 6)
         category_number = 7;
      if (count_fatjet == 1 && count_jet == 5)
         category_number = 8;
      if (count_fatjet == 1 && count_jet == 4)
         category_number = 9;
      if (count_fatjet == 1 && count_jet == 3)
         category_number = 10;
      if (count_fatjet == 1 && count_jet == 2)
         category_number = 11;
      if (count_fatjet == 1 && count_jet < 2)
         category_number = 12;

      myHists->cutflow[category_number]->Fill(4.5, weight);

      if (VBF_selection == 1)
         myHists->cutflow[category_number]->Fill(5.5, weight);

      /****************plot filling******************/
      myHists->number_of_jets[category_number]->Fill(count_jet, weight);
      myHists->number_of_central_jets[category_number]->Fill(count_central_jet, weight);

      for (int isort = 0; isort < sort_index; isort++)
         myHists->fatjet_btag_score->Fill(FatJet_DDBvL_btagsort[isort], weight);
      if (sort_index > 0)
         myHists->first_fatjet_btag_score->Fill(FatJet_DDBvL_btagsort[0], weight);
      if (sort_index > 1)
         myHists->second_fatjet_btag_score->Fill(FatJet_DDBvL_btagsort[1], weight);
      if (sort_index > 2)
         myHists->third_fatjet_btag_score->Fill(FatJet_DDBvL_btagsort[2], weight);

      if (sort_index > 0)
      {
         if (FatJet_btagsort[0].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[0].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->first_fatjet_btag_score_2match->Fill(FatJet_DDBvL_btagsort[0], weight);
         if (FatJet_btagsort[0].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[0].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->first_fatjet_btag_score_1match->Fill(FatJet_DDBvL_btagsort[0], weight);
         if (FatJet_btagsort[0].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[0].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->first_fatjet_btag_score_1match->Fill(FatJet_DDBvL_btagsort[0], weight);
         if (FatJet_btagsort[0].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[0].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->first_fatjet_btag_score_0match->Fill(FatJet_DDBvL_btagsort[0], weight);
      }

      if (sort_index > 1)
      {
         if (FatJet_btagsort[1].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[1].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->second_fatjet_btag_score_2match->Fill(FatJet_DDBvL_btagsort[1], weight);
         if (FatJet_btagsort[1].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[1].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->second_fatjet_btag_score_1match->Fill(FatJet_DDBvL_btagsort[1], weight);
         if (FatJet_btagsort[1].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[1].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->second_fatjet_btag_score_1match->Fill(FatJet_DDBvL_btagsort[1], weight);
         if (FatJet_btagsort[1].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[1].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->second_fatjet_btag_score_0match->Fill(FatJet_DDBvL_btagsort[1], weight);
      }

      if (sort_index > 2)
      {
         if (FatJet_btagsort[2].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[2].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->third_fatjet_btag_score_2match->Fill(FatJet_DDBvL_btagsort[2], weight);
         if (FatJet_btagsort[2].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[2].DeltaR(GenantiBquarkFromH) < 0.8)
            myHists->third_fatjet_btag_score_1match->Fill(FatJet_DDBvL_btagsort[2], weight);
         if (FatJet_btagsort[2].DeltaR(GenBquarkFromH) < 0.8 && FatJet_btagsort[2].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->third_fatjet_btag_score_1match->Fill(FatJet_DDBvL_btagsort[2], weight);
         if (FatJet_btagsort[2].DeltaR(GenBquarkFromH) > 0.8 && FatJet_btagsort[2].DeltaR(GenantiBquarkFromH) > 0.8)
            myHists->third_fatjet_btag_score_0match->Fill(FatJet_DDBvL_btagsort[2], weight);
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
      if (FatJet_btagDDBvL[maxindex] > 0.7)
      {

      }

      if (FatJet_btagDDBvL[maxindex] > 0.89)
      {

      }

      if (FatJet_btagDDBvL[maxindex] > 0.92)
      {

      }
*/
      /************clearing variables************/

   } // main looping

} // efficiency::loop ends
