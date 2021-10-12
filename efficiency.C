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

   // if (string(typeName) == string("TTJets_Dilept")) cout<<1<<endl;

   double xsection;
   if (string(typeName) == string("C2V_3"))
      xsection = 0.000203237376;
   if (string(typeName) == string("C2V_4p5"))
      xsection = 0.0005865984/0.3171/0.3171;
   if (string(typeName) == string("C2V_6"))
      xsection = 0.00116;
   if (string(typeName) == string("ttW"))
      xsection = 0.2043;
   if (string(typeName) == string("TTJets_DiLept"))
      xsection = 91.044;
   if (string(typeName) == string("TTJets_SingleLeptFromT"))
      xsection = 182.96;
   if (string(typeName) == string("TTJets_SingleLeptFromTbar"))
      xsection = 182.96;

   int weightnum = runChain->GetEntries();
   for (int iweight = 0; iweight < weightnum; iweight++)
   {
      runChain->GetEntry(iweight);
      myHists->weight_Scale->Fill(0.5, genEventSumw);
   }

   // if (string(typeName) == string("TTJets_DiLept")) cout<<1<<xsection<<endl;

   Int_t Nentries = fChain->GetEntries();

   /**********variables for looping**********/

   /*****************************************/
   /*********Main Looping Code Start*********/
   /*****************************************/

   for (int iLoop = 0; iLoop < Nentries; iLoop++)
   {
      fChain->GetEntry(iLoop);
      Sta_TotalNumber++;
      Sta_FileEventNumber++;
      //   double weight = double(genWeight*xsection*59600.0)/weightsum;
      //   double weight = genWeight;
      double weight = double(genWeight * xsection * 59600.0);

      /********pre-select validated events******/

      /************reconstructions***************/

      /*************event selections*************/

      /**************plot filling****************/
      double quark_pt, anti_quark_pt, nb = 0, Higgs_pt;
      double quark_eta, anti_quark_eta, quark_phi, anti_quark_phi;
      bool TestPass = false;
      for (int ipart = 0; ipart < nGenPart; ipart++)
      {
         if (GenPart_status[ipart] == 62 || GenPart_pdgId[ipart] == 25)
            Higgs_pt = GenPart_pt[ipart];

         if (GenPart_pdgId[ipart] == -5 || GenPart_pdgId[ipart] == 5)
         {
            if (GenPart_pdgId[GenPart_genPartIdxMother[ipart]] != 25)
               continue;
            if (GenPart_pdgId[ipart] == 5)
            {
               quark_pt = GenPart_pt[ipart];
               quark_eta = GenPart_eta[ipart];
               quark_phi = GenPart_phi[ipart];
            }
            if (GenPart_pdgId[ipart] == -5)
            {
               anti_quark_pt = GenPart_pt[ipart];
               anti_quark_eta = GenPart_eta[ipart];
               anti_quark_phi = GenPart_phi[ipart];
            }
            nb++;
         }
         if (GenPart_pdgId[ipart] == -11 || GenPart_pdgId[ipart] == 11)
         {
            if (GenPart_eta[ipart]<2.5 && GenPart_eta[ipart]>-2.5)
            {
                TestPass=1;
            }
         }

         if (GenPart_pdgId[ipart] == -13 || GenPart_pdgId[ipart] == 13)
         {
            if (GenPart_eta[ipart]<2.4 && GenPart_eta[ipart]>-2.4)
            {
                TestPass=1;
            }
         }
      }
         if (TestPass==1) Test_1MoreLeptPass++;
      double Bquark_deltaR;
      Bquark_deltaR = sqrt((quark_eta - anti_quark_eta) * (quark_eta - anti_quark_eta) + deltaPhi(quark_phi, anti_quark_phi) * deltaPhi(quark_phi, anti_quark_phi));

      myHists->cutflow->Fill(0.5, weight);
      myHists->total_Higgspt->Fill(Higgs_pt, weight);
      myHists->total_quarkDeltaR->Fill(Bquark_deltaR, weight);

      myHists->nBquark->Fill(nb);
      myHists->Bquark_Pt->Fill(quark_pt, weight);
      myHists->anti_Bquark_Pt->Fill(anti_quark_pt, weight);
      myHists->double_Bquark_Pt->Fill(quark_pt, weight);
      myHists->double_Bquark_Pt->Fill(anti_quark_pt, weight);
      myHists->double_Bquark_Pt_vs_Higgs_Pt->Fill(quark_pt, Higgs_pt, weight);
      myHists->double_Bquark_Pt_vs_Higgs_Pt->Fill(anti_quark_pt, Higgs_pt, weight);
      myHists->leading_Bquark_Pt->Fill(max(quark_pt, anti_quark_pt), weight);
      myHists->leading_Bquark_Pt_vs_Higgs_Pt->Fill(max(quark_pt, anti_quark_pt), Higgs_pt, weight);
      myHists->sub_Bquark_Pt->Fill(min(quark_pt, anti_quark_pt), weight);
      myHists->sub_Bquark_Pt_vs_Higgs_Pt->Fill(min(quark_pt, anti_quark_pt), Higgs_pt, weight);

      if (quark_pt > anti_quark_pt)
      {
         myHists->leading_Bquark_Eta->Fill(quark_eta, weight);
         myHists->sub_Bquark_Eta->Fill(anti_quark_eta, weight);
      }
      if (quark_pt < anti_quark_pt)
      {
         myHists->leading_Bquark_Eta->Fill(anti_quark_eta, weight);
         myHists->sub_Bquark_Eta->Fill(quark_eta, weight);
      }
      myHists->Bquark_DeltaR->Fill(Bquark_deltaR, weight);
      myHists->Bquark_DeltaR_vs_Higgs_Pt->Fill(Bquark_deltaR, Higgs_pt, weight);

      int sort_index = 0;
      int pass_cut_jet_index[200];
      double pass_cut_jet_btagDDBvL[200];

      double sum_electronLooseId = 0, sum_muonLooseId = 0, sum_electronTightId = 0, sum_muonTightId = 0;
      for (int ielectron = 0; ielectron < nElectron; ielectron++)
      {
         sum_electronLooseId += int(Electron_goodLooseId[ielectron]);
         sum_electronTightId += int(Electron_goodTightId[ielectron]);
      }
      for (int imuon = 0; imuon < nMuon; imuon++)
      {
         sum_muonLooseId += int(Muon_goodLooseId[imuon]);
         sum_muonTightId += int(Muon_goodTightId[imuon]);
      }

      double sum_leptonLooseId = sum_electronLooseId + sum_muonLooseId;
      double sum_leptonTightId = sum_electronTightId + sum_muonTightId;

      if (sum_leptonLooseId < 1)
         continue;
      myHists->cutflow->Fill(1.5, weight);

      if (sum_leptonTightId < 1)
         continue;
      myHists->cutflow->Fill(2.5, weight);

      for (int ijet = 0; ijet < nFatJet; ijet++)
      {
         if (FatJet_Filter[ijet] == false)
            continue;
         if (FatJet_pt[ijet] < 250)
            continue;
         if (FatJet_jetId[ijet] <= 0)
            continue;
         if (FatJet_eta[ijet] > 2.5)
            continue;
         if (FatJet_eta[ijet] < -2.5)
            continue;
         if (FatJet_mass[ijet] < 50)
            continue;
         if (FatJet_msoftdrop[ijet] < 40)
            continue;

         pass_cut_jet_index[sort_index] = ijet;
         pass_cut_jet_btagDDBvL[sort_index] = FatJet_btagDDBvL[ijet];
         sort_index++;
      }

      if (sort_index > 0)
      {
         int maxindex;
         for (int iindex = 0; iindex < sort_index; iindex++)
         {
            double smaller_than = 0;
            for (int iindex_compare = 0; iindex_compare < sort_index; iindex_compare++)
            {
               if (FatJet_pt[iindex_compare] > FatJet_pt[iindex])
                  smaller_than++;
            }
            if (smaller_than == 0)
               maxindex = pass_cut_jet_index[iindex];
         }

         TLorentzVector bjet1, recoHiggs;

         bjet1.SetPtEtaPhiM(FatJet_pt[maxindex], FatJet_eta[maxindex], FatJet_phi[maxindex], FatJet_mass[maxindex]);
         recoHiggs = bjet1;

         myHists->cutflow->Fill(3.5, weight);

         myHists->fatjet_btagger->Fill(FatJet_btagDDBvL[maxindex], weight);

         if (FatJet_btagDDBvL[maxindex] > 0.7)
         {
            myHists->recoHiggs_MassLL->Fill(recoHiggs.M(), weight);
            myHists->recoHiggs_msoftdropLL->Fill(FatJet_msoftdrop[maxindex], weight);
            myHists->recoHiggs_PtLL->Fill(recoHiggs.Pt(), weight);
            myHists->recoHiggs_EtaLL->Fill(recoHiggs.Eta(), weight);

            myHists->cutflow->Fill(4.5, weight);
            myHists->cut_HiggsptLL->Fill(Higgs_pt, weight);
            myHists->cut_quarkDeltaRLL->Fill(Bquark_deltaR, weight);

            //assign b quark to b jet
            double distance1 = (quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi());
            double distance2 = (anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi());
            double deltaR_qj1, deltaR_qj2;
            if (quark_pt > anti_quark_pt)
            {
               deltaR_qj1 = sqrt((quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi()));
               deltaR_qj2 = sqrt((anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi()));
            }
            if (quark_pt < anti_quark_pt)
            {
               deltaR_qj2 = sqrt((quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi()));
               deltaR_qj1 = sqrt((anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi()));
            }
            myHists->QJ_DeltaR1LL->Fill(deltaR_qj1, weight);
            myHists->QJ_DeltaR2LL->Fill(deltaR_qj2, weight);
            myHists->QJ_DeltaRLL->Fill(deltaR_qj1, weight);
            myHists->QJ_DeltaRLL->Fill(deltaR_qj2, weight);
            myHists->QJ_DeltaR1_vs_QJ_DeltaR2LL->Fill(deltaR_qj1, deltaR_qj2, weight);
            myHists->QJ_DeltaR1_vs_HiggsMLL->Fill(deltaR_qj1, recoHiggs.M(), weight);
            myHists->QJ_DeltaR2_vs_HiggsMLL->Fill(deltaR_qj2, recoHiggs.M(), weight);
            myHists->QJ_DeltaR_vs_HiggsMLL->Fill(deltaR_qj1, recoHiggs.M(), weight);
            myHists->QJ_DeltaR_vs_HiggsMLL->Fill(deltaR_qj2, recoHiggs.M(), weight);

            int region1, region2, region3;
            if (deltaR_qj1 < 0.1)
               region1 = 1;
            if (deltaR_qj1 > 0.1 && deltaR_qj1 < 0.25)
               region1 = 2;
            if (deltaR_qj1 > 0.25 && deltaR_qj1 < 0.5)
               region1 = 3;
            if (deltaR_qj2 < 0.1)
               region2 = 1;
            if (deltaR_qj2 > 0.1 && deltaR_qj2 < 0.25)
               region2 = 2;
            if (deltaR_qj2 > 0.25 && deltaR_qj2 < 0.5)
               region2 = 3;
            if (region1 == 1 && region2 == 1)
               myHists->HiggsM_QJ_DeltaR0_01LL->Fill(recoHiggs.M(), weight);
            if (region1 == 2 || region2 == 2)
               myHists->HiggsM_QJ_DeltaR01_025LL->Fill(recoHiggs.M(), weight);
            if (region1 == 3 || region2 == 3)
               myHists->HiggsM_QJ_DeltaR025_05LL->Fill(recoHiggs.M(), weight);

            if (deltaR_qj1 < 0.8 && deltaR_qj2 < 0.8)
            {
               myHists->cutflow->Fill(7.5,weight);
               myHists->cut_Higgspt_2passLL->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_2passLL->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_2passLL->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_2passLL->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 < 0.8 && deltaR_qj2 > 0.8)
            {
               myHists->cut_Higgspt_1passLL->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_1passLL->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_1passLL->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_1passLL->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 > 0.8 && deltaR_qj2 < 0.8)
            {
               myHists->cut_Higgspt_1passLL->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_1passLL->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_1passLL->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_1passLL->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 > 0.8 && deltaR_qj2 > 0.8)
            {
               myHists->cut_Higgspt_0passLL->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_0passLL->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_0passLL->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_0passLL->Fill(FatJet_msoftdrop[maxindex], weight);
            }
         }

         if (FatJet_btagDDBvL[maxindex] > 0.89)
         {
            myHists->recoHiggs_MassMM->Fill(recoHiggs.M(), weight);
            myHists->recoHiggs_msoftdropMM->Fill(FatJet_msoftdrop[maxindex], weight);
            myHists->recoHiggs_PtMM->Fill(recoHiggs.Pt(), weight);
            myHists->recoHiggs_EtaMM->Fill(recoHiggs.Eta(), weight);

            myHists->cutflow->Fill(5.5, weight);
            myHists->cut_HiggsptMM->Fill(Higgs_pt, weight);
            myHists->cut_quarkDeltaRMM->Fill(Bquark_deltaR, weight);

            //assign b quark to b jet
            double distance1 = (quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi());
            double distance2 = (anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi());

            double deltaR_qj1, deltaR_qj2;
            if (quark_pt > anti_quark_pt)
            {
               deltaR_qj1 = sqrt((quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi()));
               deltaR_qj2 = sqrt((anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi()));
            }
            if (quark_pt < anti_quark_pt)
            {
               deltaR_qj2 = sqrt((quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi()));
               deltaR_qj1 = sqrt((anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi()));
            }
            myHists->QJ_DeltaR1MM->Fill(deltaR_qj1, weight);
            myHists->QJ_DeltaR2MM->Fill(deltaR_qj2, weight);
            myHists->QJ_DeltaRMM->Fill(deltaR_qj1, weight);
            myHists->QJ_DeltaRMM->Fill(deltaR_qj2, weight);
            myHists->QJ_DeltaR1_vs_QJ_DeltaR2MM->Fill(deltaR_qj1, deltaR_qj2, weight);
            myHists->QJ_DeltaR1_vs_HiggsMMM->Fill(deltaR_qj1, recoHiggs.M(), weight);
            myHists->QJ_DeltaR2_vs_HiggsMMM->Fill(deltaR_qj2, recoHiggs.M(), weight);
            myHists->QJ_DeltaR_vs_HiggsMMM->Fill(deltaR_qj1, recoHiggs.M(), weight);
            myHists->QJ_DeltaR_vs_HiggsMMM->Fill(deltaR_qj2, recoHiggs.M(), weight);

            int region1, region2, region3;
            if (deltaR_qj1 < 0.1)
               region1 = 1;
            if (deltaR_qj1 > 0.1 && deltaR_qj1 < 0.25)
               region1 = 2;
            if (deltaR_qj1 > 0.25 && deltaR_qj1 < 0.5)
               region1 = 3;
            if (deltaR_qj2 < 0.1)
               region2 = 1;
            if (deltaR_qj2 > 0.1 && deltaR_qj2 < 0.25)
               region2 = 2;
            if (deltaR_qj2 > 0.25 && deltaR_qj2 < 0.5)
               region2 = 3;
            if (region1 == 1 && region2 == 1)
               myHists->HiggsM_QJ_DeltaR0_01MM->Fill(recoHiggs.M(), weight);
            if (region1 == 2 || region2 == 2)
               myHists->HiggsM_QJ_DeltaR01_025MM->Fill(recoHiggs.M(), weight);
            if (region1 == 3 || region2 == 3)
               myHists->HiggsM_QJ_DeltaR025_05MM->Fill(recoHiggs.M(), weight);

            if (deltaR_qj1 < 0.8 && deltaR_qj2 < 0.8)
            {
               myHists->cutflow->Fill(8.5, weight);
               myHists->cut_Higgspt_2passMM->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_2passMM->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_2passMM->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_2passMM->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 < 0.8 && deltaR_qj2 > 0.8)
            {
               myHists->cut_Higgspt_1passMM->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_1passMM->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_1passMM->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_1passMM->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 > 0.8 && deltaR_qj2 < 0.8)
            {
               myHists->cut_Higgspt_1passMM->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_1passMM->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_1passMM->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_1passMM->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 > 0.8 && deltaR_qj2 > 0.8)
            {
               myHists->cut_Higgspt_0passMM->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_0passMM->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_0passMM->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_0passMM->Fill(FatJet_msoftdrop[maxindex], weight);
            }
         }

         if (FatJet_btagDDBvL[maxindex] > 0.92)
         {
            myHists->recoHiggs_MassTT->Fill(recoHiggs.M(), weight);
            myHists->recoHiggs_msoftdropTT->Fill(FatJet_msoftdrop[maxindex], weight);
            myHists->recoHiggs_PtTT->Fill(recoHiggs.Pt(), weight);
            myHists->recoHiggs_EtaTT->Fill(recoHiggs.Eta(), weight);

            myHists->cutflow->Fill(6.5, weight);
            myHists->cut_HiggsptTT->Fill(Higgs_pt, weight);
            myHists->cut_quarkDeltaRTT->Fill(Bquark_deltaR, weight);

            //assign b quark to b jet
            double distance1 = (quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi());
            double distance2 = (anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi());
            double deltaR_qj1, deltaR_qj2;
            if (quark_pt > anti_quark_pt)
            {
               deltaR_qj1 = sqrt((quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi()));
               deltaR_qj2 = sqrt((anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi()));
            }
            if (quark_pt < anti_quark_pt)
            {
               deltaR_qj2 = sqrt((quark_eta - bjet1.Eta()) * (quark_eta - bjet1.Eta()) + deltaPhi(quark_phi, bjet1.Phi()) * deltaPhi(quark_phi, bjet1.Phi()));
               deltaR_qj1 = sqrt((anti_quark_eta - bjet1.Eta()) * (anti_quark_eta - bjet1.Eta()) + deltaPhi(anti_quark_phi, bjet1.Phi()) * deltaPhi(anti_quark_phi, bjet1.Phi()));
            }
            myHists->QJ_DeltaR1TT->Fill(deltaR_qj1, weight);
            myHists->QJ_DeltaR2TT->Fill(deltaR_qj2, weight);
            myHists->QJ_DeltaRTT->Fill(deltaR_qj1, weight);
            myHists->QJ_DeltaRTT->Fill(deltaR_qj2, weight);
            myHists->QJ_DeltaR1_vs_QJ_DeltaR2TT->Fill(deltaR_qj1, deltaR_qj2, weight);
            myHists->QJ_DeltaR1_vs_HiggsMTT->Fill(deltaR_qj1, recoHiggs.M(), weight);
            myHists->QJ_DeltaR2_vs_HiggsMTT->Fill(deltaR_qj2, recoHiggs.M(), weight);
            myHists->QJ_DeltaR_vs_HiggsMTT->Fill(deltaR_qj1, recoHiggs.M(), weight);
            myHists->QJ_DeltaR_vs_HiggsMTT->Fill(deltaR_qj2, recoHiggs.M(), weight);

            int region1, region2, region3;
            if (deltaR_qj1 < 0.1)
               region1 = 1;
            if (deltaR_qj1 > 0.1 && deltaR_qj1 < 0.25)
               region1 = 2;
            if (deltaR_qj1 > 0.25 && deltaR_qj1 < 0.5)
               region1 = 3;
            if (deltaR_qj2 < 0.1)
               region2 = 1;
            if (deltaR_qj2 > 0.1 && deltaR_qj2 < 0.25)
               region2 = 2;
            if (deltaR_qj2 > 0.25 && deltaR_qj2 < 0.5)
               region2 = 3;
            if (region1 == 1 && region2 == 1)
               myHists->HiggsM_QJ_DeltaR0_01TT->Fill(recoHiggs.M(), weight);
            if (region1 == 2 || region2 == 2)
               myHists->HiggsM_QJ_DeltaR01_025TT->Fill(recoHiggs.M(), weight);
            if (region1 == 3 || region2 == 3)
               myHists->HiggsM_QJ_DeltaR025_05TT->Fill(recoHiggs.M(), weight);

            if (deltaR_qj1 < 0.8 && deltaR_qj2 < 0.8)
            {
               myHists->cutflow->Fill(9.5, weight);
               myHists->cut_Higgspt_2passTT->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_2passTT->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_2passTT->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_2passTT->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 < 0.8 && deltaR_qj2 > 0.8)
            {
               myHists->cut_Higgspt_1passTT->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_1passTT->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_1passTT->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_1passTT->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 > 0.8 && deltaR_qj2 < 0.8)
            {
               myHists->cut_Higgspt_1passTT->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_1passTT->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_1passTT->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_1passTT->Fill(FatJet_msoftdrop[maxindex], weight);
            }
            if (deltaR_qj1 > 0.8 && deltaR_qj2 > 0.8)
            {
               myHists->cut_Higgspt_0passTT->Fill(Higgs_pt, weight);
               myHists->recoHiggs_Eta_0passTT->Fill(recoHiggs.Eta(), weight);
               myHists->recoHiggs_Pt_0passTT->Fill(recoHiggs.Pt(), weight);
               myHists->recoHiggs_msoftdrop_0passTT->Fill(FatJet_msoftdrop[maxindex], weight);
            }
         }
      }

      /************clearing variables************/

   } //main looping

} //efficiency::loop ends
