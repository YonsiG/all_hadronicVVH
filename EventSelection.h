#ifndef _EventSelection_H_
#define _EventSelection_H_

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TH3F.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TApplication.h>
#include <TEnv.h>
#include <TComplex.h>
#include <TH2D.h>
#include "Math/LorentzVector.h"

using namespace std;

void GenPart_Filter(int &Is_Hbb, TLorentzVector &GenHiggs, TLorentzVector &GenBquarkFromH, TLorentzVector &GenantiBquarkFromH,
                    TLorentzVector &GenVBFJets0, TLorentzVector &GenVBFJets1,
                    int nGenPart, int *GenPart_status, int *GenPart_pdgId, int *GenPart_genPartIdxMother,
                    float *GenPart_pt, float *GenPart_eta, float *GenPart_phi, float *GenPart_mass)
{
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

    // keep track of VBF quarks(pdgId=+-1,2,3,4,5,6), the mothers should be quarks from protons(search for pdgId), their daughters should be W and quark. Instance 2,3 are always the VBF jets. GenPart[2] and GenPart[3]
    GenVBFJets0.SetPtEtaPhiM(GenPart_pt[2], GenPart_eta[2], GenPart_phi[2], GenPart_mass[2]);
    GenVBFJets1.SetPtEtaPhiM(GenPart_pt[3], GenPart_eta[3], GenPart_phi[3], GenPart_mass[3]);

    if (double_count_Hbb == 2)
        Is_Hbb = 1;
    else
        Is_Hbb = 0;
        
}

int VBF_Selection(int VBF_method, int &VBF_jet_index0, int &VBF_jet_index1,
                  float &VBF_max_mass, float &VBF_max_DeltaEta, int count_jet, TLorentzVector *Jet)
{
    int VBF_selection = 0;

    /********* VBF max mass selection method ***********/
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
                        VBF_jet_index0 = ivbf;
                        VBF_jet_index1 = ivbf2;
                    }
                    if (Jet[ivbf].Pt() < Jet[ivbf2].Pt())
                    {
                        VBF_jet_index0 = ivbf2;
                        VBF_jet_index1 = ivbf;
                    }
                    VBF_max_mass = (Jet[ivbf] + Jet[ivbf2]).M();
                }
            }
        }
    }

    /************ VBF max DeltaEta selection method ***************/
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
                        VBF_jet_index0 = ivbf;
                        VBF_jet_index1 = ivbf2;
                    }
                    if (Jet[ivbf].Pt() < Jet[ivbf2].Pt())
                    {
                        VBF_jet_index0 = ivbf2;
                        VBF_jet_index1 = ivbf;
                    }
                    VBF_max_DeltaEta = fabs(Jet[ivbf].Eta() - Jet[ivbf2].Eta());
                }
            }
        }
    }

    /************ VBF max energy method ************/
    VBF_max_DeltaEta = 0;
    if (VBF_method == 3)
    {
        float VBF_max_Energy = 0;
        for (int ivbf = 0; ivbf < count_jet; ivbf++)
        {
            if (Jet[ivbf].E() > VBF_max_Energy)
            {
                VBF_jet_index0 = ivbf;
                VBF_max_Energy = Jet[ivbf].E();
            }
        } // select max energy jet
        VBF_max_Energy = 0;
        int opposite_eta = 0;
        for (int ivbf = 0; ivbf < count_jet; ivbf++)
        {
            if (ivbf == VBF_jet_index0)
                continue;
            if (Jet[ivbf].Eta() * Jet[VBF_jet_index0].Eta() > 0)
                continue;
            opposite_eta = 1;
            if (Jet[ivbf].E() > VBF_max_Energy)
            {
                VBF_jet_index1 = ivbf;
                VBF_max_Energy = Jet[ivbf].E();
            }
        }

        VBF_max_Energy = 0;
        if (opposite_eta == 0)
        {
            for (int ivbf = 0; ivbf < count_jet; ivbf++)
            {
                if (ivbf == VBF_jet_index0)
                    continue;
                if (Jet[ivbf].E() > VBF_max_Energy)
                {
                    VBF_jet_index1 = ivbf;
                    VBF_max_Energy = Jet[ivbf].E();
                }
            }
        }
    }

    VBF_max_mass = (Jet[VBF_jet_index0] + Jet[VBF_jet_index1]).M();
    VBF_max_DeltaEta = fabs(Jet[VBF_jet_index0].Eta() - Jet[VBF_jet_index1].Eta());
    if (VBF_max_mass > 500 && VBF_max_DeltaEta > 3)
        // (VBF_max_mass > 500)
        VBF_selection = 1;
    return VBF_selection;
}

#endif // end define _EventSelection_H_