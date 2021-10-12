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
    // std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>>*       fatjets_p4=0;
    Float_t GenPart_mass[200];
    Float_t GenPart_pt[200];
    Float_t GenPart_eta[200];
    Float_t GenPart_phi[200];
    Int_t GenPart_genPartIdxMother[200];
    Int_t GenPart_status[200];
    Int_t GenPart_pdgId[200];
    UInt_t nGenPart;
    Float_t genWeight;
    Double_t genEventSumw;

    Float_t FatJet_btagDDBvL[200];
    Float_t FatJet_eta[200];
    Float_t FatJet_pt[200];
    Float_t FatJet_phi[200];
    Float_t FatJet_mass[200];
    Int_t FatJet_jetId[200];
    // Int_t                              Jet_puId[200];
    // UChar_t                            Jet_cleanmask[200];
    Float_t FatJet_msoftdrop[200];
    Bool_t FatJet_Filter[200];

    UInt_t nFatJet;
    UInt_t nElectron;
    UInt_t nMuon;

    Bool_t Electron_goodLooseId[200];
    Bool_t Muon_goodLooseId[200];
    Bool_t Electron_goodTightId[200];
    Bool_t Muon_goodTightId[200];

    /************statistical variables*************/
    Int_t Sta_TotalNumber = 0;
    Int_t Sta_FileEventNumber;
    Int_t Test_1MoreLeptPass = 0;

    /******************functions*******************/
    efficiency(const char *infileName, const char *typeName, const char *fileNumber);
    virtual ~efficiency();
    virtual void Initial(const char *rootName, int rootNumber, const char *typeName);
    virtual void Loop(const char *typeName);
    virtual void End(int rootNumber);
    virtual void Save(int rootNumber);

    /********************plots*********************/
    makeHists *myHists;

}; //class efficiency definition ends

#endif //end define _EFFICIENCY_H_

#ifdef _EFFICIENCY_C_ //if already defined, then add functions
efficiency::efficiency(const char *infileName, const char *typeName, const char *fileNumber)
{
    myHists = new makeHists();
    TString histName = "../outfiles/" + (TString)typeName + "/" + (TString)typeName + "_" + (TString)fileNumber + "_selected.root";
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
    // if(string(typeName)==string("C2V_4p5")) runChain->SetBranchAddress("genEventSumw",&genEventSumw);
    // else runChain->SetBranchAddress("genEventSumw_",&genEventSumw);
    runChain->SetBranchAddress("genEventSumw", &genEventSumw);
    fChain->SetBranchAddress("nFatJet", &nFatJet);
    fChain->SetBranchAddress("nElectron", &nElectron);
    fChain->SetBranchAddress("nMuon", &nMuon);
    fChain->SetBranchAddress("GenPart_mass", GenPart_mass);
    fChain->SetBranchAddress("GenPart_pt", GenPart_pt);
    fChain->SetBranchAddress("GenPart_eta", GenPart_eta);
    fChain->SetBranchAddress("GenPart_phi", GenPart_phi);
    fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId);
    fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
    fChain->SetBranchAddress("GenPart_status", GenPart_status);
    fChain->SetBranchAddress("nGenPart", &nGenPart);
    fChain->SetBranchAddress("genWeight", &genWeight);
    fChain->SetBranchAddress("FatJet_btagDDBvL", FatJet_btagDDBvL);
    fChain->SetBranchAddress("FatJet_eta", FatJet_eta);
    fChain->SetBranchAddress("FatJet_phi", FatJet_phi);
    fChain->SetBranchAddress("FatJet_pt", FatJet_pt);
    fChain->SetBranchAddress("FatJet_mass", FatJet_mass);
    fChain->SetBranchAddress("FatJet_jetId", FatJet_jetId);
    // fChain->SetBranchAddress("Jet_puId",Jet_puId);
    // fChain->SetBranchAddress("Jet_cleanmask",Jet_cleanmask);
    fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop);
    fChain->SetBranchAddress("FatJet_Filter", &FatJet_Filter);

    fChain->SetBranchAddress("Electron_goodLooseId", Electron_goodLooseId);
    fChain->SetBranchAddress("Muon_goodLooseId", Muon_goodLooseId);
    fChain->SetBranchAddress("Electron_goodTightId", Electron_goodTightId);
    fChain->SetBranchAddress("Muon_goodTightId", Muon_goodTightId);

    /*
 int weightnum = runtree->GetEntries();
 for (int iweight=0;iweight<weightnum;iweight++)
 {
  runtree->GetEntry(iweight);
  weightsum=weightsum+genEventSumw;
 }
 cout<<weightsum<<endl;
*/
    //  weightsum = 6442.43; //C2V_4p5
    // weightsum = 2.71361e+07; //TTJets_DiLept
}

void efficiency::End(int rootNumber)
{
    cout << "**Running: Free Rootfile: " << rootNumber << endl;
    cout << "*************************" << endl;
    cout << "*********outputs*********" << endl;
    cout << "Total event number in this file:      " << Sta_FileEventNumber << endl;
    cout << "**The files contain " << Test_1MoreLeptPass << " events have at least one pass eta cut lepton" <<endl;
    delete fChain->GetCurrentFile();
}

void efficiency::Save(int rootNumber)
{
    cout << "**Running: " << rootNumber << "  rootfiles finished" << endl;
    cout << "**The files contain " << Sta_TotalNumber << " events" << endl;
    cout << "**The files contain " << Test_1MoreLeptPass << " events have at least one pass eta cut lepton" <<endl;

    if (myHists)
    {
        myHists->saveHists();
        delete myHists;
    }
    cout << "**Running histograms saved" << endl;
    // TFile *outfile = new TFile("evts.txt","W");
    // outfile->Write();
    // outfile->Close();
}
#endif //end ifdef _EFFICIENCY_C_
