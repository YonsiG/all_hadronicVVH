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

#endif //end define _TOOLS_H_
