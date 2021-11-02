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

bool Jet_central_Select(TLorentzVector Jet)
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

bool Jet_kinVBF_select(TLorentzVector Jet)
{
   if (Jet.Pt() < 25)
       return false;
   if (fabs(Jet.Eta() > 4.7))
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


bool Electron_fakableId_ttH(bool looseId, float pt, bool convVeto, int tightCharge, int lostHits, float JetDeepFlav, bool triggerSafeNoIso, float mvaTTH, bool mvaFall17V2noIso_WP80, float jetRelIso, float WP)
{
   if (looseId)return false;
   if (pt <= 10.)  return false; 
   if (!convVeto ) return false;
   if (tightCharge!=2) return false;
   if (int(lostHits > 0)) return false;
   if (JetDeepFlav>WP) return false;
   if (!triggerSafeNoIso) return false;
   if (mvaTTH<= 0.8) 
   {
      if (!mvaFall17V2noIso_WP80) return false; 
      if (jetRelIso >= 0.7) return false; 
   }
   return true;
}

bool Muon_fakableId_ttH(bool looseId, float pt, float JetDeepFlav, float ptErr, float mvaTTH, float jetRelIso, bool WP)
{
  if (looseId)return false;
  if (pt<=10.) return false;
  if (JetDeepFlav>WP) return false;
  if (ptErr / pt >= 0.2) return false;
  if (mvaTTH <= 0.85) 
  { 
     if (jetRelIso >= 0.5) return false;
  }
  return true;
}

bool Electron_tightId_ttH(bool looseId, float pt, bool convVeto, int tightCharge, int lostHits, float JetDeepFlav, bool triggerSafeNoIso, float mvaTTH, float WP)
{
   if (looseId)return false;
   if (pt <= 10.)  return false; 
   if (!convVeto ) return false;
   if (tightCharge!=2) return false;
   if (int(lostHits > 0)) return false;
   if (JetDeepFlav>WP) return false;
   if (!triggerSafeNoIso) return false;
   if (mvaTTH<= 0.8) return false;
   return true;
}

bool Muon_tightId_ttH(bool looseId, float pt, float JetDeepFlav, float ptErr, float mvaTTH, bool mediumId, float WP)
{
  if (looseId)return false;
  if (pt<=10.) return false;
  if (JetDeepFlav>WP) return false;
  if (ptErr / pt >= 0.2) return false;
  if (mvaTTH <= 0.85) return false;
  if (!mediumId) return false;
  return true;
}

bool Electron_triggerSafeNoIso(float eta, float deltaEtaSC, float hoe, float eInvMinusPInv, float sieie)
{
    // Calculate absolute value of supercluster eta
    float SC_absEta = fabs(eta) + deltaEtaSC;
    if (hoe >= 0.1) return false; 
    if (eInvMinusPInv <= -0.04) return false; 
    if (SC_absEta <= 1.479) 
    {  // Barrel
        if (sieie >= 0.011) return false; 
    }
    else if (SC_absEta > 1.479 && SC_absEta < 2.5)
    {  // Endcaps
       if (sieie >= 0.030) { return false; }
    }
    return true;
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
float mediumWP(int year)
{
 if (year == 2016) return 0.2489;
 else if (year == 2017) return 0.3040;
 else return 0.2783;
} 
#endif //end define _TOOLS_H_
