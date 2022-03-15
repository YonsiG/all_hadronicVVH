#ifndef _EFFICIENCY_H_
#define _EFFICIENCY_H_

#include "makeHists.h"
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
class efficiency
{
public:
    /***************test variables*****************/
    bool isRead;

    /*****************Tree input*******************/
    TTree *fChain;
    TTree *runChain;

    /***********Tree branch vairables**************/
    Float_t genWeight;
    Double_t genEventSumw;

    Bool_t HLT_AK8PFHT800_TrimMass50;
    Bool_t HLT_PFHT1050;
    Bool_t HLT_PFJet500;
    Bool_t HLT_AK8PFJet500;
    Bool_t HLT_AK8PFJet400_TrimMass30;
    Bool_t HLT_AK8PFJet420_TrimMass30;

    UInt_t nGenPart;
    UInt_t nFatJet;
    UInt_t nJet;
    UInt_t nElectron;
    UInt_t nMuon;

    Int_t GenPart_genPartIdxMother[200];
    Int_t GenPart_pdgId[200];
    Int_t GenPart_status[200];
    Float_t GenPart_mass[200];
    Float_t GenPart_pt[200];
    Float_t GenPart_eta[200];
    Float_t GenPart_phi[200];

    Float_t FatJet_eta[200];
    Float_t FatJet_pt[200];
    Float_t FatJet_phi[200];
    Float_t FatJet_mass[200];
    Int_t FatJet_jetId[200];
    Float_t FatJet_msoftdrop[200];
    Float_t FatJet_particleNet_WvsQCD[200];
    Float_t FatJet_particleNet_ZvsQCD[200];
    Float_t FatJet_particleNet_HbbvsQCD[200];
    Float_t FatJet_particleNet_mass[200];
    Float_t FatJet_particleNetMD_Xbb[200];
    Float_t FatJet_particleNetMD_Xcc[200];
    Float_t FatJet_particleNetMD_Xqq[200];
    Float_t FatJet_particleNetMD_QCD[200];

    Float_t Jet_eta[200];
    Float_t Jet_phi[200];
    Float_t Jet_pt[200];
    Float_t Jet_mass[200];
    Float_t Jet_qgl[200];

    Float_t Electron_mass[200];
    Float_t Electron_pt[200];
    Float_t Electron_eta[200];
    Float_t Electron_phi[200];
    Float_t Electron_deltaEtaSC[200];
    Float_t Electron_dxy[200];
    Float_t Electron_dz[200];
    Float_t Electron_sip3d[200];
    Float_t Electron_miniPFRelIso_all[200];
    UChar_t Electron_lostHits[200];
    Bool_t Electron_mvaFall17V2noIso_WPL[200];

    Float_t Muon_mass[200];
    Float_t Muon_pt[200];
    Float_t Muon_eta[200];
    Float_t Muon_phi[200];
    Float_t Muon_dxy[200];
    Float_t Muon_dz[200];
    Float_t Muon_sip3d[200];
    Float_t Muon_miniPFRelIso_all[200];
    Bool_t Muon_looseId[200];

    /***************particle variables***************/
    TLorentzVector GenBquarkFromH;
    TLorentzVector GenantiBquarkFromH;
    TLorentzVector GenHiggs;
    TLorentzVector GenVBFJets[2];

    TLorentzVector FatJet[200];
    TLorentzVector Jet[200];
    TLorentzVector Electron[200];
    TLorentzVector Muon[200];

    /************statistical variables*************/
    Int_t Sta_TotalNumber = 0;
    Int_t Sta_FileEventNumber;

    /*************sorting variables****************/
    TLorentzVector FatJet_btagsort[200];
    Float_t FatJet_msoftdrop_btagsort[200]; 
    Float_t FatJet_WvsQCD_btagsort[200];
    Float_t FatJet_ZvsQCD_btagsort[200];
    Float_t FatJet_HbbvsQCD_btagsort[200];
    Float_t FatJet_mass_btagsort[200];
    Float_t FatJet_Xbb_modified_btagsort[200];
    Float_t FatJet_Xcc_btagsort[200];
    Float_t FatJet_Xqq_btagsort[200];
    Float_t FatJet_QCD_btagsort[200];
    Float_t FatJet_Xccqq_modified_btagsort[200];

    TLorentzVector FatJet_allsort[200];
    Float_t FatJet_msoftdrop_allsort[200]; 
    Float_t FatJet_WvsQCD_allsort[200];
    Float_t FatJet_ZvsQCD_allsort[200];
    Float_t FatJet_HbbvsQCD_allsort[200];
    Float_t FatJet_mass_allsort[200];
    Float_t FatJet_Xbb_modified_allsort[200];
    Float_t FatJet_Xcc_allsort[200];
    Float_t FatJet_Xqq_allsort[200];
    Float_t FatJet_QCD_allsort[200];
    Float_t FatJet_Xccqq_modified_allsort[200];
    /******************functions*******************/
    efficiency(const char *infileName, const char *typeName, const char *fileNumber);
    virtual ~efficiency();
    virtual void Initial(const char *rootName, int rootNumber, const char *typeName);
    virtual void Loop(const char *typeName);
    virtual void End(int rootNumber);
    virtual void Save(int rootNumber);

    /********************plots*********************/
    makeHists *myHists;

}; // class efficiency definition ends

#endif // end define _EFFICIENCY_H_

#ifdef _EFFICIENCY_C_ // if already defined, then add functions
efficiency::efficiency(const char *infileName, const char *typeName, const char *fileNumber)
{
    myHists = new makeHists();
    TString histName = "../outfiles/" + (TString)typeName + "_" + (TString)fileNumber + "_selected.root";
    myHists->createHists(histName);
}

efficiency::~efficiency() {}

void efficiency::Initial(const char *rootName, int rootNumber, const char *typeName)
{
    Sta_FileEventNumber = 0;
    cout << "**Running starting Rootfile: " << rootNumber << endl;

    TTree *tree;
    TFile *file = (TFile *)gROOT->GetListOfFiles()->FindObject(rootName);
    if (!file)
        file = new TFile(rootName);
    tree = (TTree *)gDirectory->Get("Events");

    TTree *runtree;
    runtree = (TTree *)gDirectory->Get("Runs");

    fChain = tree;
    runChain = runtree;
    isRead = true;

    double fileEntries = fChain->GetEntries();
    if (fileEntries == 0)
        isRead = false;

    /***************set branches**************/
    runChain->SetBranchAddress("genEventSumw", &genEventSumw);
    fChain->SetBranchAddress("genWeight", &genWeight);

    fChain->SetBranchAddress("HLT_AK8PFHT800_TrimMass50",&HLT_AK8PFHT800_TrimMass50);
    fChain->SetBranchAddress("HLT_PFHT1050",&HLT_PFHT1050);
    fChain->SetBranchAddress("HLT_PFJet500",&HLT_PFJet500);
    fChain->SetBranchAddress("HLT_AK8PFJet500",&HLT_AK8PFJet500);
    fChain->SetBranchAddress("HLT_AK8PFJet400_TrimMass30",&HLT_AK8PFJet400_TrimMass30);
    fChain->SetBranchAddress("HLT_AK8PFJet420_TrimMass30",&HLT_AK8PFJet420_TrimMass30);

    fChain->SetBranchAddress("nGenPart", &nGenPart);
    fChain->SetBranchAddress("nFatJet", &nFatJet);
    fChain->SetBranchAddress("nJet", &nJet);
    fChain->SetBranchAddress("nElectron", &nElectron);
    fChain->SetBranchAddress("nMuon", &nMuon);

    fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    fChain->SetBranchAddress("GenPart_status", GenPart_status);
    fChain->SetBranchAddress("GenPart_mass", GenPart_mass);
    fChain->SetBranchAddress("GenPart_pt", GenPart_pt);
    fChain->SetBranchAddress("GenPart_eta", GenPart_eta);
    fChain->SetBranchAddress("GenPart_phi", GenPart_phi);

    fChain->SetBranchAddress("FatJet_eta", FatJet_eta);
    fChain->SetBranchAddress("FatJet_phi", FatJet_phi);
    fChain->SetBranchAddress("FatJet_pt", FatJet_pt);
    fChain->SetBranchAddress("FatJet_mass", FatJet_mass);
    fChain->SetBranchAddress("FatJet_jetId", FatJet_jetId);
    fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop);
    fChain->SetBranchAddress("FatJet_particleNet_WvsQCD", FatJet_particleNet_WvsQCD);
    fChain->SetBranchAddress("FatJet_particleNet_ZvsQCD", FatJet_particleNet_ZvsQCD);
    fChain->SetBranchAddress("FatJet_particleNet_HbbvsQCD", FatJet_particleNet_HbbvsQCD);
    fChain->SetBranchAddress("FatJet_particleNet_mass", FatJet_particleNet_mass);
    fChain->SetBranchAddress("FatJet_particleNetMD_Xbb", FatJet_particleNetMD_Xbb);
    fChain->SetBranchAddress("FatJet_particleNetMD_Xcc", FatJet_particleNetMD_Xcc);
    fChain->SetBranchAddress("FatJet_particleNetMD_Xqq", FatJet_particleNetMD_Xqq);
    fChain->SetBranchAddress("FatJet_particleNetMD_QCD", FatJet_particleNetMD_QCD);

    fChain->SetBranchAddress("Jet_eta", Jet_eta);
    fChain->SetBranchAddress("Jet_phi", Jet_phi);
    fChain->SetBranchAddress("Jet_pt", Jet_pt);
    fChain->SetBranchAddress("Jet_mass", Jet_mass);
    fChain->SetBranchAddress("Jet_qgl", Jet_qgl);

    fChain->SetBranchAddress("Electron_mass", Electron_mass);
    fChain->SetBranchAddress("Electron_pt", Electron_pt);
    fChain->SetBranchAddress("Electron_eta", Electron_eta);
    fChain->SetBranchAddress("Electron_phi", Electron_phi);
    fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC);
    fChain->SetBranchAddress("Electron_dxy", Electron_dxy);
    fChain->SetBranchAddress("Electron_dz", Electron_dz);
    fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d);
    fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all);
    fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits);
    fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", Electron_mvaFall17V2noIso_WPL);

    fChain->SetBranchAddress("Muon_mass", Muon_mass);
    fChain->SetBranchAddress("Muon_pt", Muon_pt);
    fChain->SetBranchAddress("Muon_eta", Muon_eta);
    fChain->SetBranchAddress("Muon_phi", Muon_phi);
    fChain->SetBranchAddress("Muon_dxy", Muon_dxy);
    fChain->SetBranchAddress("Muon_dz", Muon_dz);
    fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d);
    fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all);
    fChain->SetBranchAddress("Muon_looseId", Muon_looseId);
}

void efficiency::End(int rootNumber)
{
    cout << "**Running: Free Rootfile: " << rootNumber << endl;
    cout << "*************************" << endl;
    cout << "*********outputs*********" << endl;
    cout << "Total event number in this file:      " << Sta_FileEventNumber << endl;
    delete fChain->GetCurrentFile();
}

void efficiency::Save(int rootNumber)
{
    cout << "**Running: " << rootNumber << "  rootfiles finished" << endl;
    cout << "**The files contain " << Sta_TotalNumber << " events" << endl;

    if (myHists)
    {
        myHists->saveHists();
        delete myHists;
    }
    cout << "**Running histograms saved" << endl;
}
#endif // end ifdef _EFFICIENCY_C_
