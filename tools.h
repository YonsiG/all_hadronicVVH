#ifndef _TOOLS_H_
#define _TOOLS_H_

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

double deltaPhi(double phi1, double phi2)
{

    double PI = 3.14159265;
    double result = phi1 - phi2;
    while (result > PI)
        result -= 2 * PI;
    while (result <= -PI)
        result += 2 * PI;
    return result;
}

bool FatJet_kinematic_Select(TLorentzVector FatJet, int FatJet_jetId, float FatJet_msoftdrop)
{
    if (FatJet.Pt() < 250)
        return false;
    if (FatJet_jetId <= 0) //seperate something clustered together but not jets
        return false;
    if (fabs(FatJet.Eta()) > 2.5)
        return false;
    if (FatJet_msoftdrop < 40)
        return false;
    return true;
}

bool Jet_kinematic_Select(TLorentzVector Jet)
{
    if (Jet.Pt() < 30)
        return false;
    if (fabs(Jet.Eta() > 2.5))
        return false;

    /*      if (Jet_cleanmask[ijet] == 0)
        continue;
      if (Jet_puId[ijet] != 7)
        continue;
      if (Jet_jetId[ijet] != 6)
        continue;
    */

    return true;
}

bool Electron_LooseId_ttH(TLorentzVector Electron, float Electron_deltaEtaSC, float Electron_dxy, float Electron_dz, float Electron_sip3d, float Electron_miniPFRelIso_all, int Electron_lostHits, bool Electron_mvaFall17V2noIso_WPL)
{
    if (Electron.Pt() <= 7.)
    {
        return false;
    }
    if (fabs(Electron.Eta() + Electron_deltaEtaSC) >= 2.5) // track Eta + cluster Eta (not perfectly aligned)
    {
        return false;
    }
    if (fabs(Electron_dxy) >= 0.05)
    {
        return false;
    }
    if (fabs(Electron_dz) >= 0.1)
    {
        return false;
    }
    if (fabs(Electron_sip3d) >= 8)
    {
        return false;
    }
    if (Electron_miniPFRelIso_all >= 0.4)
    {
        return false;
    }
    if (int(Electron_lostHits) > 1)
    {
        return false;
    }
    if (!Electron_mvaFall17V2noIso_WPL)
    {
        return false;
    }
    return true;
}

bool Muon_LooseId_ttH(TLorentzVector Muon, float Muon_dxy, float Muon_dz, float Muon_sip3d, float Muon_miniPFRelIso_all, bool Muon_looseId)
{
    if (Muon.Pt() <= 5.)
    {
        return false;
    }
    if (fabs(Muon.Eta()) >= 2.4)
    {
        return false;
    }
    if (fabs(Muon_dxy) >= 0.05)
    {
        return false;
    } // W bosons will have prompt muons, won't have displacement muons away from primary vertex
    if (fabs(Muon_dz) >= 0.1)
    {
        return false;
    }
    if (fabs(Muon_sip3d) >= 8)
    {
        return false;
    }
    if (Muon_miniPFRelIso_all >= 0.4)
    {
        return false;
    }
    if (!Muon_looseId)
    {
        return false;
    } // loose POG ID
    return true;
}

#endif //end define _TOOLS_H_
