#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TObjString.h"
#include "TGraphErrors.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooPolyVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataSet.h"

#include "RooGaussian.h"
#include "RooFFTConvPdf.h"

#include <ROOT/RDataFrame.hxx>
#include "ROOT/RVec.hxx"
#include <ROOT/RLogger.hxx>
#include <ROOT/RSnapshotOptions.hxx>

#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace std;
using namespace RooFit;
using namespace ROOT;
using namespace ROOT::VecOps;

#define RVecD ROOT::RVec<double>

//ECAL[2][5][6] cell
int ecalID(int iplane, int ix, int iy) {
  int ret = 0;
  ret = ix * 6 + iy;
  return ret;
}

//HCAL[4][3][3]
int hcalID(int iplane, int ix, int iy) {
  int ret = 0;
  ret = iplane * 9 + ix * 3 + iy;
  return ret;
}

double dimuAngle(double &DIMU_momX, double &DIMU_momZ){
  return std::atan2(DIMU_momX, DIMU_momZ);
}

double dimuDeflection(double &DIMU_momX, double &DIMU_momZ, double &DIMU_vertexX, double &DIMU_vertexZ){
  double ECALendZ = 3430.96;
  double HCAL2startZ = 7264.2;
  double totalZ = (HCAL2startZ - ECALendZ) + (ECALendZ - DIMU_vertexZ);
  double angle = std::atan2(DIMU_momX, DIMU_momZ);
  
  return angle*totalZ + DIMU_vertexX;
}

double HCAL12rvalue(RVecD &hcal) {
  double total = 0.;
  double per = 0.;
  double ret = -1.;
  for (int imod = 1; imod < 3; imod++) {
    for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 3; iy++) {
        auto id = hcalID(imod, ix, iy);
        total += hcal[id];
        if (ix == 1 && iy == 1) continue;
        per += hcal[id];
      }
    }
  }
  ret = (total == 0) ? 1 : per/total;
  return ret;
}
/*the following function is used to take r-value of HCAL modules
 */
double HCAL0rvalue(RVecD &hcal) {
  double total = 0.;
  double per = 0.;
  double ret = -1.;
    for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 3; iy++) {
        auto id = hcalID(0, ix, iy);
        total += hcal[id];
        if (ix == 1 && iy == 1) continue;
        per += hcal[id];
      }
    }
  ret = (total == 0) ? 1 : per/total;
  return ret;
}

double HCAL1rvalue(RVecD &hcal) {
  double total = 0.;
  double per = 0.;
  double ret = -1.;
    for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 3; iy++) {
        auto id = hcalID(1, ix, iy);
        total += hcal[id];
        if (ix == 1 && iy == 1) continue;
        per += hcal[id];
      }
    }
  ret = (total == 0) ? 1 : per/total;
  return ret;
}

double HCAL2rvalue(RVecD &hcal) {
  double total = 0.;
  double per = 0.;
  double ret = -1.;
    for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 3; iy++) {
        auto id = hcalID(2, ix, iy);
        total += hcal[id];
        if (ix == 1 && iy == 1) continue;
        per += hcal[id];
      }
    }
  ret = (total == 0) ? 1 : per/total;
  return ret;
}

/* HC0ratio
 * (C0 - C2)/(C0 + C2)
 1:left side, 0:center, -1:right side
 */
double Hcal0Ratio (RVecD &hcal) {
  double SumLeft(0), SumRight(0);
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(0, ix, iy);
      if (ix == 0) SumLeft += hcal[id];
      if (ix == 2) SumRight += hcal[id];
    }
  }
  if ((SumLeft+SumRight) <= 0) return 0.;
  return ((SumLeft - SumRight) / (SumLeft + SumRight));
}
/* HC0ratio
 * (C0 - C2)/(C0 + C2)
 1:left side, 0:center, -1:right side
 */
double Hcal1Ratio (RVecD &hcal) {
  double SumLeft(0), SumRight(0);
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(1, ix, iy);
      if (ix == 0) SumLeft += hcal[id];
      if (ix == 2) SumRight += hcal[id];
    }
  }
  if ((SumLeft+SumRight) <= 0) return 0.;
  return ((SumLeft - SumRight) / (SumLeft + SumRight));
}
/* HC0ratio
 * (C0 - C2)/(C0 + C2)
 1:left side, 0:center, -1:right side
 */
double Hcal2Ratio (RVecD &hcal) {
  double SumLeft(0), SumRight(0);
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(2, ix, iy);
      if (ix == 0) SumLeft += hcal[id];
      if (ix == 2) SumRight += hcal[id];
    }
  }
  if ((SumLeft+SumRight) <= 0) return 0.;
  return ((SumLeft - SumRight) / (SumLeft + SumRight));
}
/* ECratio
 * (5x6 - 3x3)/(5x6)

 */
double EcalRatio (RVecD &ECALenergy) {
  double SumAll(0), Sum9(0);
  for (int is = 0; is < 2; is++) {
    for (int ix = 0; ix < 5; ix++) {
      for (int iy = 0; iy < 6; iy++) {
        auto id = ecalID(is, ix, iy);
        SumAll += ECALenergy[id];
        if ((ix >= 1) && (ix <= 3) && (iy >= 2) && (iy <= 4))
          Sum9 += ECALenergy[id];
      }
    }
  }
  if (SumAll <= 0) return -1;
  return ((SumAll - Sum9) / SumAll);
}

//Straw multiplicity cut
bool filterFunSTmult(int &straw3X, int &straw3Y) {
  bool ret;
  ret = ((straw3X > 0) && (straw3X <= 5) && (straw3Y > 0) && (straw3Y <= 5));
  return ret;
}

//Centermost cell cut: maximum in-time energy of ECAL only is in center cell
bool filterFunECAL23(RVecD &ECAL_ene) {
  double eMax = 0;
  int xMax, yMax;
  for (int ix = 0; ix < 5; ix++) {
    for (int iy = 0; iy < 6; iy++) {
      auto id = ecalID(1, ix, iy);
      if (ECAL_ene[id] > eMax) {
        eMax = ECAL_ene[id];
        xMax = ix;
        yMax = iy;
      }
    }
  }
  return ((xMax == 2) && (yMax == 2));
}

bool filterFunSRD(RVecD &SRD_ene) {
  bool srdtok[3];
  for (int ii = 0; ii < 3; ii++) {
    srdtok[ii] = ((SRD_ene[ii] > 0.001) && (SRD_ene[ii] < 0.08));
  }
  return (srdtok[0] && srdtok[1] && srdtok[2]);
  //return (srdtok[0] && srdtok[1] && SRD_ene[2] > 0.002) || (srdtok[0] && srdtok[2] && SRD_ene[1] > 0.002) || (srdtok[1] && srdtok[2] && SRD_ene[0] > 0.002);
}

bool filterFunVETO(RVecD &VETO_ene){
  const double vetoCut0 = 0.032;
  const double vetoCut1 = 0.032;
  const double vetoCut2 = 0.032;
  return ((VETO_ene[0] < vetoCut0) && (VETO_ene[1] < vetoCut1) && (VETO_ene[2] < vetoCut2));
} 
bool filterFunVETO2(RVecD &VETO_ene){
  return (VETO_ene[1] > 0.015);
} 

/*ECcuts: EcalTotal<150 and PreShower>0.3
 * EcalTotal: ECAL+PRS
 */
bool filterFunEC(RVecD &PRS_ene, RVecD &ECAL_ene) {
  double sECAL_ene = 0;
  double sPS_ene = 0;
  for (int ix = 0; ix < 5; ix++) {
    for (int iy = 0; iy < 6; iy++) {
      auto id = ecalID(0, ix, iy);
      sPS_ene += PRS_ene[id];
      sECAL_ene += (ECAL_ene[id]+PRS_ene[id]);
    }
  }
  return ((sECAL_ene < 150.)  && (sPS_ene > 0.7));
}
//track
bool filterFunTrack(double &mom_genfit_upMM, double &chi2_genfit_upMM) {
  return (mom_genfit_upMM > 0. && mom_genfit_upMM < 9999. && chi2_genfit_upMM > 0.01);
}

/* ECratio cut
 * (5x6 - 3x3)/(5x6) < 0.07
*/
bool filterFunECratio(RVecD &PRS_ene, RVecD &ECAL_ene) {
  double SumAll(0), Sum9(0);
  for (int ix = 0; ix < 5; ix++) {
    for (int iy = 0; iy < 6; iy++) {
      auto id = ecalID(0, ix, iy);
      SumAll += PRS_ene[id];
      SumAll += ECAL_ene[id];
      if ((ix >= 1) && (ix <= 3) && (iy >= 2) && (iy <= 4)){
        Sum9 += PRS_ene[id];
        Sum9 += ECAL_ene[id];
      }
    }
  }
  if (SumAll <= 0)
    return false;
  return ((SumAll - Sum9) / SumAll < .07);
}

/* HCmuon cut
 * HC2 > 2 GeV || (HC0(1,1) > 1.5 && HC1(1,1) > 1.5 && HC2 > 1.3 GeV)
 *
 */
bool filterFun12(RVecD &HCAL_ene) {
  double HC2(0), HC011(0), HC111(0);
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(2, ix, iy);
      HC2 += HCAL_ene[id];
      auto id2 = hcalID(0, 1, 1);
      HC011 = HCAL_ene[id2];
      auto id3 = hcalID(1, 1, 1);
      HC111 = HCAL_ene[id3];
    }
  }

  bool muon = (HC2 > 2.) || (HC011 > 1.5 && HC111 > 1.5 && HC2 > 1.3);
  return !muon;
}

/* HCtopo
 * Reject event if (NumHCcells ≥ 1 (Ecell > 0.5) && EnergyHC – E(1,1) > E(1,1)), (in HC0), E(1,1) – central cell.
 *
 */
bool filterFunHCALtopo(RVecD &HCAL_ene) {
  int nCells = 0;
  double thr = 0.5;
  double HC0 = 0;
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(0, ix, iy);
      HC0 += HCAL_ene[id];
      if (HCAL_ene[id] > thr) //toropin 40 adc
        nCells++;
    }
  }
  if (nCells >= 1) {
    return (HC0 < 2 * HCAL_ene[hcalID(0, 1, 1)]);
  }
  return true;
}

/*SRD cut total in-time energy cut < 120 MeV
*/
bool filterFunSRDtot(RVecD &SRD_ene) {
  double eSRDtot(0);
  for (int ii = 0; ii < 3; ii++) {
    eSRDtot += SRD_ene[ii];
  }
  //cout<<eSRDtot<<endl;
  return (eSRDtot < 0.120);
}

/*ECAL shower chi2 cut based on total ECAL energy
 * As implemented by L. Marsicano based on dimuons and setting the efficiency to be 95%
 */
const int nsteps = 16;
const double ChiCuts[nsteps] = {2.5, 3.86, 3.25, 2.915, 2.83, 2.85, 2.93, 3.035, 3.16, 3.285, 3.385, 3.51, 3.6, 3.685, 3.74, 3.695};
const double Eranges[nsteps+1] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80};

bool filterFunECALchi2(double &ECAL_chi2, double &ECAL_ene_tot) {
  bool ret = false;
  double chicut = 0;
  for (int ii = 0; ii < nsteps; ii++) {
    if (ECAL_ene_tot >= Eranges[ii] && ECAL_ene_tot < Eranges[ii+1]) {
      chicut = ChiCuts[ii];
    }
  }
  if (ECAL_ene_tot >= Eranges[nsteps]) {
    chicut = 4;
  }
  if (ECAL_chi2 < chicut) {
    ret = true;
  }

  return ret;
}

bool physTrigger(RVecD &PRS_ene, RVecD &ECAL_ene) {
  double sPS_ene = 0;
  double sECAL_ene = 0;
  for (int ix = 1; ix < 4; ix++) { //Threshold from 3x4 central cells
    for (int iy = 1; iy < 5; iy++) {
      auto id = ecalID(0, ix, iy);
      sPS_ene += PRS_ene[id];
      sECAL_ene += ECAL_ene[id]; //Threshold on ECAL1 only
    }
  }
  return ((sECAL_ene < 76)  && (sPS_ene > 0.4)); // < 76, > 0.4
}

/* HC2muon cut
 * HC2 > 4 GeV && HC2 < 8 GeV
 *
 */
bool filterFunMuonHC2(RVecD &HCALenergy) {
  double HC2(0);
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(2, ix, iy);
      HC2 += HCALenergy[id];
    }
  }

  bool muon = (HC2 > 4.) && (HC2 < 8.);
  return muon;
}

/* HC1muon cut
 * HC1 > 4 GeV && HC1 < 8 GeV
 *
 */
bool filterFunMuonHC1(RVecD &HCALenergy) {
  double HC1(0);
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(1, ix, iy);
      HC1 += HCALenergy[id];
    }
  }

  bool muon = (HC1 > 4.) && (HC1 < 8.);
  return muon;
}

/* HC0muon cut
 * HC0 > 4 GeV && HC0 < 8 GeV
 *
 */
bool filterFunMuonHC0(RVecD &HCALenergy) {
  double HC0(0);
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(0, ix, iy);
      HC0 += HCALenergy[id];
    }
  }

  bool muon = (HC0 > 4.) && (HC0 < 8.);
  return muon;
}

bool filterFunMuonHC0_center(RVecD &HCALenergy){
  double HC0(0); 
  
  auto id = hcalID(0, 1, 1);
  HC0 += HCALenergy[id];  

  bool center_muon = (HC0 > 4.) && (HC0 < 8.);
  return center_muon;
}

bool filterFunMuonHC1_center(RVecD &HCALenergy){
  double HC1(0); 

  auto id = hcalID(1, 1, 1);
  HC1 += HCALenergy[id];

  bool center_muon = (HC1 > 4.) && (HC1 < 8.);
  return center_muon;
}

bool filterFunMuonHC2_center(RVecD &HCALenergy){
  double HC2(0); 

  auto id = hcalID(2, 1, 1);
  HC2 += HCALenergy[id];

  bool center_muon = (HC2 > 4.) && (HC2 < 8.);
  return center_muon;
}

ROOT::RVec<double> HcalPedestral(RVecD &ene) {

  const double pedestral = 0.2;

  ROOT::RVec<double> HCAL_ene_out;
  for (int is = 0; is < 4; is++) {
    for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 3; iy++) {
        auto id = hcalID(is, ix, iy);
        if (ene[id] > pedestral)
          HCAL_ene_out.push_back(ene[id]);
        else
          HCAL_ene_out.push_back(0.);
      }
    }
  }
  return HCAL_ene_out;
}

  //Straw multiplicity cut
bool filterStrawMult(int &straw3X, int &straw3Y) {
  bool ret;
  ret = ((straw3X > 0) && (straw3X <= 5) && (straw3Y > 0) && (straw3Y <= 5));
  return ret;
}

int StrawMultHist(vector<double> &straw) {
  return straw.size();
}

//Straw Acceptence cut
bool filterAcceptenceST(vector<int> &straw11_PDGID, vector<int> &straw12_PDGID) { 
  // Number of mu+ and mu- in ST11 and ST12
  int N_mup11 = 0; // mu+
  int N_mun11 = 0; // mu-

  int N_mup12 = 0;
  int N_mun12 = 0;

  for(const auto &pdgid : straw11_PDGID) {
    // If any other particle other than muons hit the Straws, return 0.
    if (!(pdgid == 13 || pdgid == -13)) return 0;
    if (pdgid == -13) ++N_mup11; // PGDID mu+ = -13 
    if (pdgid == 13) ++N_mun11; // PGDID mu- = 13
  }

  // As long as we see hits of mu+ AND mu- in BOTH straws, return 1.
  if ((N_mup11 > 0 && N_mun11 > 0)) return 1;

  return 0;
}

//Straw multiplicity cut
bool filterMPMM(vector<double> &MM6x, vector<double> &MM7x) {
  if(MM6x.size() == MM7x.size() && MM6x.size() == 1) return 1;
  else return 0;
}

// Filter of straw11 so that we have only one hit in the + and - region.
bool filterMMEachSide(vector<double> &MM6x, vector<double> &MM7x) {
  int leftSide = 0;
  int rightSide = 0;
  // Check X position of straws
  for (size_t i = 0; i < MM6x.size(); ++i) {
    // If more than 1 in each side, continue
    if(leftSide > 1 || rightSide > 1) continue;
    if(MM6x[i] > 0.) leftSide++; // mu+
    else if(MM7x[i] < 0.) rightSide++; // mu-
  }
  // Check if y position is ok
  //for(auto y : straw11Y){
  //  if(!(y > -93. && y < 106.)) return 0;
  //}
  // Return true only if one hit in each side
  if(leftSide == 1 && rightSide == 1) return 1;

  return 0;
}

//Straw multiplicity cut
bool filterLargeStrawMult(vector<double> &strawX, vector<double> &strawY) {
  if(strawX.size() == strawY.size() && strawX.size() == 2) return 1;
  else return 0;
}

//Straw multiplicity cut
bool filterSingleMuonST(vector<double> &strawX, vector<double> &strawY) {
  if(strawX.size() == strawY.size() && strawX.size() == 1) return 1;
  else return 0;
}


bool filterStrawFalseMuons(vector<int> &straw11_PDGID){
  int mu_p = 0;
  int mu_n = 0;
  for(const auto &pdgid : straw11_PDGID){
    if (pdgid == 13) mu_n = -1;
    else if(pdgid == -13) mu_p = 1;
  }
  if(!(mu_n == -1 && mu_p == 1) && straw11_PDGID.size() > 3 && straw11_PDGID.size() < 9) return 1;
  return 0;
}

bool filterSingleMuons(vector<int> &straw11_PDGID){
  int mu_p = 0;
  int mu_n = 0;
  for(const auto &pdgid : straw11_PDGID){
    if (pdgid == 13) mu_n = 1;
    else if(pdgid == -13) mu_p = 1;
  }
  if((mu_n == 1 && mu_p == 0 || mu_n == 0 && mu_p == 1) && straw11_PDGID.size() > 0 && straw11_PDGID.size() < 5) return 1;
  //if(straw11_PDGID.size() > 0 && straw11_PDGID.size() < 5) return 1;
  return 0;
}

bool filterDimuons(vector<int> &straw11_PDGID){
  int mu_p = 0;
  int mu_n = 0;
  //cout << straw11_PDGID.size() << endl;
  for(const auto &pdgid : straw11_PDGID){
    //cout << pdgid << endl;
    if (pdgid == 13) mu_n = 1;
    else if(pdgid == -13) mu_p = 1;
  }
  if((mu_n == 1 && mu_p == 1) && straw11_PDGID.size() > 4 && straw11_PDGID.size() < 9) return 1;
  //if(straw11_PDGID.size() > 4 && straw11_PDGID.size() < 9) return 1;
  return 0;
}

bool filterStrawBackground(vector<int> &straw11_PDGID){
  int mu_p = 0;
  int mu_n = 0;
  int background = 0;
  for(const auto &pdgid : straw11_PDGID){
    if (pdgid == 13) mu_n = -1;
    else if(pdgid == -13) mu_p = 1;
    else background = 1;
  }
  if(background == 1) return 1;
  return 0;
}

// Filter of straw11 so that we have only one hit in the + and - region.
bool filterStraw11EachSide(vector<double> &straw11X, vector<double> &straw11Y) {
  int leftSide = 0;
  int rightSide = 0;
  // Check X position of straws
  for (size_t i = 0; i < straw11X.size(); ++i) {
    // If more than 1 in each side, continue
    if(leftSide > 1 || rightSide > 1) continue;
    if(straw11X[i] > 0.) leftSide++; // mu+
    else if(straw11X[i] < 0.) rightSide++; // mu-
  }
  // Check if y position is ok
  //for(auto y : straw11Y){
  //  if(!(y > -93. && y < 106.)) return 0;
  //}
  // Return true only if one hit in each side
  if(leftSide == 1 && rightSide == 1) return 1;

  return 0;
}

// Filter of straw11 so that we have only one hit in the + and - region.
bool filterStraw11EachSide2(vector<double> &straw11X, vector<double> &straw11Y, vector<int> &straw11_PDGID) {
  int leftSide = 0;
  int rightSide = 0;
  // Check X position of straws
  for (size_t i = 0; i < straw11X.size(); ++i) {
    // If more than 1 in each side, continue
    //if(leftSide > 1 || rightSide > 1) continue;
    if(straw11X[i] > -749. && straw11_PDGID[i] == -13) leftSide++; // mu+
    else if(straw11X[i] < -749. && straw11_PDGID[i] == 13) rightSide++; // mu-
  }
  // Check if y position is ok
  //for(auto y : straw11Y){
  //  if(!(y > -93. && y < 106.)) return 0;
  //}
  // Return true only if one hit in each side
  if(leftSide > 1 && rightSide > 1) return 1;

  return 0;
}

// Filter of straw11 so that we have only one hit in the + and - region.
bool filterStraw11MuonSeperation(vector<double> &straw11X, vector<int> &straw11_PDGID) {
  // Check X position of straws
  for (size_t i = 0; i < straw11X.size(); ++i) {
    // -749 is the center in x position of straw in experiment, not straw cooordinate
    if(straw11X[i] > (-749.) && straw11_PDGID[i] == 13) return 0; // Left side of straw. Remove events that are mu-
    if(straw11X[i] < (-749.) && straw11_PDGID[i] == -13) return 0; // Right side of straw. Remove events that are mu+
  }
  return 1;
}

// Filter of straw11 with a two regions, being 2 sigma from the center of the two peaks.
bool filterStraw11Region(vector<double> &straw11X, vector<double> &straw11Y) {
  int leftSide = 0;
  int rightSide = 0;
  // Check X position of straws
  for(auto x : straw11X){
    // If more than 1 in each side, continue
    if(leftSide > 1 || rightSide > 1) continue;
    if(x > -950. && x < -780.) rightSide++;
    else if(x < -420. && x > -640.) leftSide++;
    else return 0;
  }
  // Check if y position is ok
  //for(auto y : straw11Y){
  //  if(!(y > -93. && y < 106.)) return 0;
  //}
  // Return true only if one hit in each side
  if(leftSide == 1 && rightSide == 1) return 1;

  return 0;
}

bool filterMagnetIronBlock(vector<double> &straw11Y, vector<double> &straw12X, vector<double> &straw12Y) {
  for(auto y : straw11Y){
    if(!(y > -100. && y < 105.)) return 0;
  }
  for(auto y : straw12Y){
    if(!(y > -110. && y < 130.)) return 0;
  }
  //for(auto x : straw12X){
  //  for(auto y : straw12Y){
  //    if(((x > -130. && x < 60.) && (y > -150. && y < -100.))) return 0;
  //  }
  //}
  return 1;
}

// Filter of straw12 w.r.t. position in x and y
bool filterStraw12(vector<double> &straw12X, vector<double> &straw12Y) {
  int leftSide = 0;
  int rightSide = 0;
  // Check X position of straws
  for(auto x : straw12X){
    // If more than 1 in each side, continue
    if(leftSide > 1 || rightSide > 1) continue;
    if(x > 80.) rightSide++;
    else if(x < -150.) leftSide++;
  }
  // Check if y position is ok
  for(auto y : straw12Y){
    if(!(y > -150. && y < 150.)) return 0;
  }
  // Return true only if one hit in each side
  if(leftSide == 1 && rightSide == 1) return 1;

  return 0;
}

/* ECratio cut
 * (5x6 - 3x3)/(5x6) < 0.06
 * Center is in cell (2,2)
 */
bool filterECALratio(RVecD &ECALenergy) {
  double SumAll(0), Sum9(0);
  for (int is = 0; is < 2; is++) {
    for (int ix = 0; ix < 5; ix++) {
      for (int iy = 0; iy < 6; iy++) {
        auto id = ecalID(is, ix, iy);
        SumAll += ECALenergy[id];
        if ((ix >= 1) && (ix <= 3) && (iy >= 1) && (iy <= 3))
          Sum9 += ECALenergy[id];
      }
    }
  }
  if (SumAll <= 0)
    return false;
  return ((SumAll - Sum9) / SumAll < .06);
}

double filterDimuAngles(double &DIMU_momX, double &DIMU_momZ){
  double angle = std::atan2(DIMU_momX, DIMU_momZ);
  return ((angle > -0.0225) && (angle < -0.0215));  
}

//Function to produce histograms, to be called on the result of the "Filter" call
std::vector<ROOT::RDF::RResultPtr<TH1>> doHisto(ROOT::RDF::RNode df, int icut, string cutName) {

  std::vector<ROOT::RDF::RResultPtr<TH1>> ret;

  // ECAL
  df = df.Define("eecal", "Sum(ECALenergy)+Sum(ECALpenergy)");
  df = df.Define("eprs", "Sum(ECALpenergy)"); //take PRS energies: first 5x6 in ECAL
  df = df.Define("eecal_ratio", EcalRatio, { "ECALenergy" });
  df = df.Define("eprs_central", "ECALenergy[14]"); // ECAL[0][2][2]
  df = df.Define("eecal_central", "ECALenergy[44]"); // ECAL[1][2][2]
  

  //df = df.Define("eEcentral", "ECALenergy[2][2]");
  //HCAL
  df = df.Define("ehcal", "Sum(eH_PDS)"); //total energy
  df = df.Define("eH0", "Sum(Take(eH_PDS,9))"); //take first 9 HC0 distribution
  df = df.Define("eH1", "Sum(Take(eH_PDS,18)) - Sum(Take(eH_PDS,9))"); //take first 9 HC1 distribution
  df = df.Define("eH2", "Sum(Take(eH_PDS,27)) - Sum(Take(eH_PDS,18))"); //take first 9 HC2 distribution
  df = df.Define("eH3", "Sum(Take(eH_PDS,-9))"); //take first 9 HC3 distribution
  df = df.Define("ehcal01", "Sum(Take(eH_PDS,18))"); //take first 9+9 energies: HC0, HC1
  df = df.Define("ehcal012", "Sum(Take(eH_PDS,27))"); //take first 9+9+9 energies: HC0, HC1, HC2

  // Periferal HCAL cells
  df = df.Define("eH0_00", "eH_PDS[0]"); 
  df = df.Define("eH0_01", "eH_PDS[1]"); 
  df = df.Define("eH0_02", "eH_PDS[2]");
  df = df.Define("eH0_10", "eH_PDS[3]");
  df = df.Define("eH0_11", "eH_PDS[4]");
  df = df.Define("eH0_12", "eH_PDS[5]");
  df = df.Define("eH0_20", "eH_PDS[6]");
  df = df.Define("eH0_21", "eH_PDS[7]");
  df = df.Define("eH0_22", "eH_PDS[8]");

  df = df.Define("eH1_00", "eH_PDS[9]");  
  df = df.Define("eH1_01", "eH_PDS[10]");  
  df = df.Define("eH1_02", "eH_PDS[11]");
  df = df.Define("eH1_10", "eH_PDS[12]");
  df = df.Define("eH1_11", "eH_PDS[13]");
  df = df.Define("eH1_12", "eH_PDS[14]");
  df = df.Define("eH1_20", "eH_PDS[15]");
  df = df.Define("eH1_21", "eH_PDS[16]");
  df = df.Define("eH1_22", "eH_PDS[17]");

  df = df.Define("eH2_00", "eH_PDS[18]");  
  df = df.Define("eH2_01", "eH_PDS[19]");  
  df = df.Define("eH2_02", "eH_PDS[20]");
  df = df.Define("eH2_10", "eH_PDS[21]");
  df = df.Define("eH2_11", "eH_PDS[22]");
  df = df.Define("eH2_12", "eH_PDS[23]");
  df = df.Define("eH2_20", "eH_PDS[24]");
  df = df.Define("eH2_21", "eH_PDS[25]");
  df = df.Define("eH2_22", "eH_PDS[26]");
  
  df = df.Define("eH0rvalue", HCAL0rvalue, { "eH_PDS" });// HC0 r-value
  df = df.Define("eH1rvalue", HCAL1rvalue, { "eH_PDS" });// HC1 r-value
  df = df.Define("eH2rvalue", HCAL2rvalue, { "eH_PDS" });// HC2 r-value
  df = df.Define("eH12rvalue", HCAL12rvalue, { "eH_PDS" });
  df = df.Define("eH0_ratio", Hcal0Ratio, { "eH_PDS" });// left-right asymmetry in HC0
  df = df.Define("eH1_ratio", Hcal1Ratio, { "eH_PDS" });// left-right asymmetry in HC1
  df = df.Define("eH2_ratio", Hcal2Ratio, { "eH_PDS" });// left-right asymmetry in HC2

  df = df.Define("dimuonAngleP", dimuAngle, { "DIMUp_momX", "DIMUp_momZ" });
  df = df.Define("dimuonAngleN", dimuAngle, { "DIMUn_momX", "DIMUn_momZ" });

  df = df.Define("dimuonDeflectionP", dimuDeflection, { "DIMUp_momX", "DIMUp_momZ", "DIMU_vertexX", "DIMU_vertexZ" });
  df = df.Define("dimuonDeflectionN", dimuDeflection, { "DIMUn_momX", "DIMUn_momZ", "DIMU_vertexX", "DIMU_vertexZ" });

  df = df.Define("mpStrawRecon11X", StrawMultHist, { "Strawx11" });
  df = df.Define("mpStrawRecon11Y", StrawMultHist, { "Strawy11" });
  df = df.Define("mpStrawRecon12X", StrawMultHist, { "Strawx12" });
  df = df.Define("mpStrawRecon12Y", StrawMultHist, { "Strawy12" });

  // Dimuon information
 ret.push_back(df.Histo1D( { Form("DIMU_Etot_%i", icut), "DIMU_Etot", 100u, 0., 100.}, "DIMU_Etot"));
 ret.push_back(df.Histo1D( { Form("DIMUp_E_%i", icut), "DIMUp_E", 100u, 0., 100.}, "DIMUp_E"));
 ret.push_back(df.Histo1D( { Form("DIMUn_E_%i", icut), "DIMUn_E", 100u, 0., 100.}, "DIMUn_E"));

  //Straw multiplicity
  ret.push_back(df.Histo1D( { Form("mpStraw3X_%i", icut), "mpStraw3X", 20u, -0.5,19.5}, "mpStraw3X"));
  ret.push_back(df.Histo1D( { Form("mpStraw3Y_%i", icut), "mpStraw3Y", 20u, -0.5,19.5}, "mpStraw3Y"));
  ret.push_back(df.Histo1D( { Form("mpStraw4X_%i", icut), "mpStraw4X", 20u, -0.5,19.5}, "mpStraw4X"));
  ret.push_back(df.Histo1D( { Form("mpStraw4Y_%i", icut), "mpStraw4Y", 20u, -0.5,19.5}, "mpStraw4Y"));

  ret.push_back(df.Histo1D( { Form("mpStrawRecon11X_%i", icut), "mpStrawRecon11X", 20u, -0.5,19.5}, "mpStrawRecon11X"));
  ret.push_back(df.Histo1D( { Form("mpStrawRecon11Y_%i", icut), "mpStrawRecon11Y", 20u, -0.5,19.5}, "mpStrawRecon11Y"));
  ret.push_back(df.Histo1D( { Form("mpStrawRecon12X_%i", icut), "mpStrawRecon12X", 20u, -0.5,19.5}, "mpStrawRecon12X"));
  ret.push_back(df.Histo1D( { Form("mpStrawRecon12Y_%i", icut), "mpStrawRecon12Y", 20u, -0.5,19.5}, "mpStrawRecon12Y"));

  ret.push_back(df.Histo1D( { Form("Straw5X_%i", icut), "Straw5X", 200u, -100., 100.}, "Strawx5"));
  ret.push_back(df.Histo1D( { Form("Straw5Y_%i", icut), "Straw5Y", 200u, -100., 100.}, "Strawy5"));
  ret.push_back(df.Histo1D( { Form("Straw6X_%i", icut), "Straw6X", 200u, -100., 100.}, "Strawx6"));
  ret.push_back(df.Histo1D( { Form("Straw6Y_%i", icut), "Straw6Y", 200u, -100., 100.}, "Strawy6"));
  ret.push_back(df.Histo1D( { Form("Straw11X_%i", icut), "Straw11X", 300u, -600., 600.}, "Strawx11"));
  ret.push_back(df.Histo1D( { Form("Straw11Y_%i", icut), "Straw11Y", 300u, -600., 600.}, "Strawy11"));
  ret.push_back(df.Histo1D( { Form("Straw12X_%i", icut), "Straw12X", 300u, -600., 600.}, "Strawx12"));
  ret.push_back(df.Histo1D( { Form("Straw12Y_%i", icut), "Straw12Y", 300u, -600., 600.}, "Strawy12"));
  ret.push_back(df.Histo1D( { Form("Straw11X_true_%i", icut), "Straw11X true", 200u, -1500., 0.}, "Strawx11_truth"));
  //ret.push_back(df.Histo1D( { Form("Straw11Y_true_%i", icut), "Straw11Y true", 200u, -600., 600.}, "Strawy11_truth"));
  ret.push_back(df.Histo1D( { Form("Straw12X_true_%i", icut), "Straw12X true", 200u, -1500., 0.}, "Strawx12_truth"));
  //ret.push_back(df.Histo1D( { Form("Straw12Y_true_%i", icut), "Straw12Y true", 200u, -600., 600.}, "Strawy12_truth"));
  
  //ret.push_back(df.Histo1D( { Form("Straw11_TrackEnergy_%i", icut), "Straw11 Energy", 100u, 0., 100.}, "Straw11_TrackEnergy"));
  //ret.push_back(df.Histo1D( { Form("Straw12_TrackEnergy_%i", icut), "Straw12 Energy", 100u, 0., 100.}, "Straw12_TrackEnergy"));

  //ret.push_back(df.Histo1D( { Form("mpMuonsHCAL_%i", icut), "mpMuonsHCAL", 20u, -0.5,19.5}, "HCALTrue_muons"));
  //ret.push_back(df.Histo1D( { Form("mpMuonsST12_%i", icut), "mpMuonsST12", 20u, -0.5,19.5}, "ST12True_muons"));

  ret.push_back(df.Histo2D( { Form("Straw5map_%i", icut), Form("%s ; Straw 5 X [mm]; Straw 5 Y [mm]", cutName.c_str()), 400u, -100., 100., 400u, -100., 100. }, "Strawx5", "Strawy5"));
  ret.push_back(df.Histo2D( { Form("Straw6map_%i", icut), Form("%s ; Straw 6 X [mm]; Straw 6 Y [mm]", cutName.c_str()), 400u, -100., 100., 400u, -100., 100. }, "Strawx6", "Strawy6"));
  ret.push_back(df.Histo2D( { Form("Straw11map_%i", icut), Form("%s ; Straw 11 X [mm]; Straw 11 Y [mm]", cutName.c_str()), 400u, -1200., 1200., 400u, -600., 600. }, "Strawx11", "Strawy11"));
  ret.push_back(df.Histo2D( { Form("Straw12map_%i", icut), Form("%s ; Straw 12 X [mm]; Straw 12 Y [mm]", cutName.c_str()), 400u, -1200., 1200., 400u, -600., 600. }, "Strawx12", "Strawy12"));
  ret.push_back(df.Histo2D( { Form("Straw11map_true_%i", icut), Form("%s ; Straw 11 X [mm]; Straw 11 Y [mm]", cutName.c_str()), 400u, -1400., 100., 400u, -400., 400. }, "Strawx11_truth", "Strawy11_truth"));
  ret.push_back(df.Histo2D( { Form("Straw12map_true_%i", icut), Form("%s ; Straw 12 X [mm]; Straw 12 Y [mm]", cutName.c_str()), 400u, -1500., 0., 400u, -400., 400. }, "Strawx12_truth", "Strawy12_truth"));
  ret.push_back(df.Histo2D( { Form("HCAL_truedetector_%i", icut), Form("%s ; ST HCAL X [mm]; ST HCAL Y [mm]", cutName.c_str()), 400u, -700., 700., 400u, -400., 400. }, "HCALTrueX", "HCALTrueY"));
  //ret.push_back(df.Histo2D( { Form("ST12_truedetector_%i", icut), Form("%s ; Straw 12 X [mm]; Straw 12 Y [mm]", cutName.c_str()), 400u, -700., 700., 400u, -400., 400. }, "ST12TrueX", "ST12TrueY"));

  // MM
  ret.push_back(df.Histo1D( { Form("MM3X_%i", icut), "MM3X", 100u, -125.,125.}, "MM3x"));
  ret.push_back(df.Histo1D( { Form("MM3Y_%i", icut), "MM3Y", 100u, -45.,45.}, "MM3y"));
  ret.push_back(df.Histo1D( { Form("MM4X_%i", icut), "MM4X", 100u, -125.,125.}, "MM4x"));
  ret.push_back(df.Histo1D( { Form("MM4Y_%i", icut), "MM4Y", 100u, -45.,45.}, "MM4y"));
  ret.push_back(df.Histo1D( { Form("MM5X_%i", icut), "MM5X", 100u, -125.,125.}, "MM5x"));
  ret.push_back(df.Histo1D( { Form("MM5Y_%i", icut), "MM5Y", 100u, -45.,45.}, "MM5y"));
  ret.push_back(df.Histo1D( { Form("MM6X_%i", icut), "MM6X", 100u, -125.,125.}, "MM6x"));
  ret.push_back(df.Histo1D( { Form("MM6Y_%i", icut), "MM6Y", 100u, -45.,45.}, "MM6y"));
  ret.push_back(df.Histo1D( { Form("MM7X_%i", icut), "MM7X", 100u, -125.,125.}, "MM7x"));
  ret.push_back(df.Histo1D( { Form("MM7Y_%i", icut), "MM7Y", 100u, -45.,45.}, "MM7y"));

  ret.push_back(df.Histo2D( { Form("MM5map_%i", icut), Form("%s ; MM 5 X [mm]; MM 5 Y [mm]", cutName.c_str()), 200u, -125.,125., 200u, -45.,45. }, "MM5x", "MM5y"));
  ret.push_back(df.Histo2D( { Form("MM6map_%i", icut), Form("%s ; MM 6 X [mm]; MM 6 Y [mm]", cutName.c_str()), 200u, -125.,125., 200u, -45.,45. }, "MM6x", "MM6y"));
  ret.push_back(df.Histo2D( { Form("MM7map_%i", icut), Form("%s ; MM 7 X [mm]; MM 7 Y [mm]", cutName.c_str()), 200u, -125.,125., 200u, -45.,45. }, "MM7x", "MM7y"));

  //Veto energy
  ret.push_back(df.Define("eveto", "VETO[0]").Histo1D( { Form("Veto01_E_%i", icut), "Veto01_E", 400u, 0., 0.1 }, "eveto"));
  ret.push_back(df.Define("eveto", "VETO[1]").Histo1D( { Form("Veto23_E_%i", icut), "Veto23_E", 400u, 0., 0.1 }, "eveto"));
  ret.push_back(df.Define("eveto", "VETO[2]").Histo1D( { Form("Veto45_E_%i", icut), "Veto45_E", 400u, 0., 0.1 }, "eveto"));

  //ECAL
  ret.push_back(df.Histo1D( { Form("ECAL_%i", icut), "ECAL", 400u, 0., 200. }, "eecal"));
  ret.push_back(df.Histo1D( { Form("PRS_%i", icut), "PRS", 400u, 0., 200. }, "eprs"));
  ret.push_back(df.Histo1D( { Form("ECAL_ratio_%i", icut), "ECAL total ratio (all - 3x3) / all", 201u, -0.0025, 1.0025 }, "eecal_ratio"));
  ret.push_back(df.Histo1D( { Form("ECALcentral_%i", icut), "ECAL center", 400u, 0., 200. }, "eecal_central"));
  ret.push_back(df.Histo1D( { Form("PRScentral_%i", icut), "PRS center", 400u, 0., 200. }, "eprs_central"));

  //HCAL
  ret.push_back(df.Histo1D( { Form("HCAL_%i", icut), "HCAL", 400u, 0., 100. }, "ehcal"));
  ret.push_back(df.Histo1D( { Form("HCAL0_%i", icut), "HCAL0", 400u, 0., 120. }, "eH0"));
  ret.push_back(df.Histo1D( { Form("HCAL1_%i", icut), "HCAL1", 400u, 0., 120. }, "eH1"));
  ret.push_back(df.Histo1D( { Form("HCAL2_%i", icut), "HCAL2", 400u, 0., 120. }, "eH2"));
  ret.push_back(df.Histo1D( { Form("HCAL3_%i", icut), "HCAL3", 400u, 0., 120. }, "eH3"));
  ret.push_back(df.Histo1D( { Form("HCAL01_%i", icut), "HCAL01", 400u, 0., 200. }, "ehcal01"));
  ret.push_back(df.Histo1D( { Form("HCAL012_%i", icut), "HCAL012", 400u, 0., 200. }, "ehcal012"));
  ret.push_back(df.Histo1D( { Form("HCAL_PDS_%i", icut), "HCAL Pedestral", 400u, 0., 200. }, "eH_PDS"));
  ret.push_back(df.Histo1D( { Form("HCAL0_R_%i", icut), "HCAL0 rvalue", 101u, -0.005, 1.005 }, "eH0rvalue"));
  ret.push_back(df.Histo1D( { Form("HCAL1_R_%i", icut), "HCAL1 rvalue", 101u, -0.005, 1.005 }, "eH1rvalue"));
  ret.push_back(df.Histo1D( { Form("HCAL2_R_%i", icut), "HCAL2 rvalue", 101u, -0.005, 1.005 }, "eH2rvalue"));
  ret.push_back(df.Histo1D( { Form("HCAL12_R_%i", icut), "HCAL12 rvalue", 101u, -0.005, 1.005 }, "eH12rvalue"));
  ret.push_back(df.Histo1D( { Form("HCAL0_ratio_%i", icut), "HCAL0 ratio", 400u, -1., 1. }, "eH0_ratio"));
  ret.push_back(df.Histo1D( { Form("HCAL1_ratio_%i", icut), "HCAL1 ratio", 400u, -1., 1. }, "eH1_ratio"));
  ret.push_back(df.Histo1D( { Form("HCAL2_ratio_%i", icut), "HCAL2 ratio", 400u, -1., 1. }, "eH2_ratio"));

  ret.push_back(df.Histo1D( { Form("HCAL0_00_%i", icut), "HCAL0_00", 100u, 0., 10. }, "eH0_00"));
  ret.push_back(df.Histo1D( { Form("HCAL0_01_%i", icut), "HCAL0_01", 100u, 0., 10. }, "eH0_01"));
  ret.push_back(df.Histo1D( { Form("HCAL0_02_%i", icut), "HCAL0_02", 100u, 0., 10. }, "eH0_02"));
  ret.push_back(df.Histo1D( { Form("HCAL0_10_%i", icut), "HCAL0_10", 100u, 0., 10. }, "eH0_10"));
  ret.push_back(df.Histo1D( { Form("HCAL0_11_%i", icut), "HCAL0_11", 100u, 0., 30. }, "eH0_11"));
  ret.push_back(df.Histo1D( { Form("HCAL0_12_%i", icut), "HCAL0_12", 100u, 0., 10. }, "eH0_12"));
  ret.push_back(df.Histo1D( { Form("HCAL0_20_%i", icut), "HCAL0_20", 100u, 0., 10. }, "eH0_20"));
  ret.push_back(df.Histo1D( { Form("HCAL0_21_%i", icut), "HCAL0_21", 100u, 0., 10. }, "eH0_21"));
  ret.push_back(df.Histo1D( { Form("HCAL0_22_%i", icut), "HCAL0_22", 100u, 0., 10. }, "eH0_22"));

  ret.push_back(df.Histo1D( { Form("HCAL1_00_%i", icut), "HCAL1_00", 100u, 0., 10. }, "eH1_00"));
  ret.push_back(df.Histo1D( { Form("HCAL1_01_%i", icut), "HCAL1_01", 100u, 0., 10. }, "eH1_01"));
  ret.push_back(df.Histo1D( { Form("HCAL1_02_%i", icut), "HCAL1_02", 100u, 0., 10. }, "eH1_02"));
  ret.push_back(df.Histo1D( { Form("HCAL1_10_%i", icut), "HCAL1_10", 100u, 0., 10. }, "eH1_10"));
  ret.push_back(df.Histo1D( { Form("HCAL1_11_%i", icut), "HCAL1_11", 100u, 0., 30. }, "eH1_11"));
  ret.push_back(df.Histo1D( { Form("HCAL1_12_%i", icut), "HCAL1_12", 100u, 0., 10. }, "eH1_12"));
  ret.push_back(df.Histo1D( { Form("HCAL1_20_%i", icut), "HCAL1_20", 100u, 0., 10. }, "eH1_20"));
  ret.push_back(df.Histo1D( { Form("HCAL1_21_%i", icut), "HCAL1_21", 100u, 0., 10. }, "eH1_21"));
  ret.push_back(df.Histo1D( { Form("HCAL1_22_%i", icut), "HCAL1_22", 100u, 0., 10. }, "eH1_22"));

  ret.push_back(df.Histo1D( { Form("HCAL2_00_%i", icut), "HCAL2_00", 100u, 0., 10. }, "eH2_00"));
  ret.push_back(df.Histo1D( { Form("HCAL2_01_%i", icut), "HCAL2_01", 100u, 0., 10. }, "eH2_01"));
  ret.push_back(df.Histo1D( { Form("HCAL2_02_%i", icut), "HCAL2_02", 100u, 0., 10. }, "eH2_02"));
  ret.push_back(df.Histo1D( { Form("HCAL2_10_%i", icut), "HCAL2_10", 100u, 0., 10. }, "eH2_10"));
  ret.push_back(df.Histo1D( { Form("HCAL2_11_%i", icut), "HCAL2_11", 100u, 0., 30. }, "eH2_11"));
  ret.push_back(df.Histo1D( { Form("HCAL2_12_%i", icut), "HCAL2_12", 100u, 0., 10. }, "eH2_12"));
  ret.push_back(df.Histo1D( { Form("HCAL2_20_%i", icut), "HCAL2_20", 100u, 0., 10. }, "eH2_20"));
  ret.push_back(df.Histo1D( { Form("HCAL2_21_%i", icut), "HCAL2_21", 100u, 0., 10. }, "eH2_21"));
  ret.push_back(df.Histo1D( { Form("HCAL2_22_%i", icut), "HCAL2_22", 100u, 0., 10. }, "eH2_22"));

  //ECAL vs HCAL
  ret.push_back(df.Histo2D( { Form("hHC01vsEC_cuts_%i", icut), Form("%s ; ECAL0+ECAL1 [GeV]; HCAL0+HCAL1 [GeV]", cutName.c_str()), 200u, 0., 120., 200u, 0.,
        120. }, "eecal", "ehcal01"));
  ret.push_back(df.Histo2D( { Form("hHC012vsEC_cuts_%i", icut), Form("%s ; ECAL0+ECAL1 [GeV]; HCAL0+HCAL1+HCAL2 [GeV]", cutName.c_str()), 200u, 0., 120., 200u, 0.,
        120. }, "eecal", "ehcal012"));

  // Dimuon angle
  //ret.push_back(df.Histo1D( { Form("dimuonAngleP_%i", icut), "dimuonAngleP", 400u, -0.1,0.1}, "dimuonAngleP"));
  //ret.push_back(df.Histo1D( { Form("dimuonAngleN_%i", icut), "dimuonAngleN", 400u, -0.1,0.1}, "dimuonAngleN"));

  // Dimuon deflection
  //ret.push_back(df.Histo1D( { Form("dimuonDeflectionP_%i", icut), "dimuonDeflectionP", 400u, -1000., 1000.}, "dimuonDeflectionP"));
  //ret.push_back(df.Histo1D( { Form("dimuonDeflectionN_%i", icut), "dimuonDeflectionN", 400u, -1000., 1000.}, "dimuonDeflectionN"));
 
  return ret;
}

//Get all outputs
std::vector<ROOT::RDF::RResultPtr<TH1>> allH;
std::vector<ROOT::RDF::RResultPtr<ULong64_t>> allN;
std::vector<std::string> allCutName;

std::vector<ROOT::RDF::RResultPtr<TH1>> allH_ord;
std::vector<ROOT::RDF::RResultPtr<ULong64_t>> allN_ord;
std::vector<std::string> allCutName_ord;

void postCutProcess(int idx, string cutName, ROOT::RDF::RNode df, bool order) {
  if(!order){
    allN.push_back(df.Count());
    allCutName.push_back(cutName);
    auto ret = doHisto(df, idx, cutName);
    for (auto h : ret) {
      allH.push_back(h);
    }
  }
  else{
    allN_ord.push_back(df.Count());
    allCutName_ord.push_back(cutName);
    auto ret = doHisto(df, idx, cutName);
    for (auto h : ret) {
      allH_ord.push_back(h);
    }
  }
}

void sim_ana_dimuons(string fIn = "", string fOut = "") {
  string fnameIn, fnameOut;
  if (fIn == "") {
    fnameIn = "./19032024_ST11fix.root";
  } else {
    fnameIn = fIn;
  }

  if (fOut == "") {
    fnameOut = "simulation_out.root";
  } else {
    fnameOut = fOut;
  }
  
  // Load full dataframe
  ROOT::RDataFrame d0("ana_tree", fnameIn.c_str());

  auto d1 = d0.Define("eEC", "Sum(ECALenergy)+Sum(ECALpenergy)");
  d1 = d1.Define("eH_PDS", HcalPedestral , { "HCALenergy" });
  d1 = d1.Define("eH0_blind", "Sum(Take(eH_PDS,9))"); //take first 9 HC0 distribution
  d1 = d1.Define("eH1_blind", "Sum(Take(eH_PDS,18)) - Sum(Take(eH_PDS,9))"); //take first 9 HC1 distribution
  d1 = d1.Define("eH2_blind", "Sum(Take(eH_PDS,27)) - Sum(Take(eH_PDS,18))"); //take first 9 HC2 distribution
  d1 = d1.Define("eH3_blind", "Sum(Take(eH_PDS,-9))"); //take first 9 HC3 distribution

  auto dMagnet = d1.Filter(filterMagnetIronBlock, { "Strawy11", "Strawx12", "Strawy12" });
  //auto dMagnet_truth = d1.Filter(filterMagnetIronBlock, { "Strawy11_truth", "Strawx12_truth", "Strawy12_truth" });

  auto dBlind = dMagnet.Filter("eH0_blind > 4 && eH1_blind > 4 && eH2_blind > 4");
  postCutProcess(0, "AllEvents", dMagnet, 0);
  auto dF1 = dMagnet.Filter(physTrigger, { "ECALpenergy", "ECALenergy" });
  postCutProcess(1, "EC < 76 && PS > 0.4", dF1, 0);
  auto dF3 = dF1.Filter(filterFunSRD, { "SRD" });
  postCutProcess(2, "SRD ", dF3, 0);
  auto dF4 = dF3.Filter(filterFunTrack, { "BestMomupGenfit", "BestMomupChisq" });
  postCutProcess(3, "Track Quality", dF4, 0);
  auto dF5 = dF4.Filter("abs(BestMomupGenfit-100)<10");
  postCutProcess(4, "Momentum", dF5, 0);
  auto dF6 = dF5.Filter("abs(InAngle)<3.");
  postCutProcess(5, "Angle12", dF6, 0); 
  auto dF7 = dF6.Filter(filterStrawMult, { "mpStraw3X", "mpStraw3Y" });
  postCutProcess(6, "Straw3 Multiplicity Cut", dF7, 0);
  auto dF72 = dF7.Filter("Sum(VHCALenergy) < 1.5");
  postCutProcess(7, "VHCAL cut", dF72, 0);
  auto dF8 = dF72.Filter(filterStrawMult, { "mpStraw4X", "mpStraw4Y" });
  postCutProcess(8, "Straw4 Multiplicity Cut", dF8, 0);
  auto dF9 = dF8.Filter(filterFunECAL23, { "ECALenergy" }); // Check if correct!
  postCutProcess(9, "cutEC(2,2) ", dF9, 0); // CENTER ECAL CELL IS INDEED (2,2)
  auto dF10 = dF9.Filter(filterFunEC, { "ECALpenergy", "ECALenergy" });
  postCutProcess(10, "EC ", dF10, 0);
  auto dF12 = dF10.Filter(filterECALratio, { "ECALenergy" });
  postCutProcess(11, "ECratio ", dF12, 0);
  auto dF15 = dF12.Filter("Sum(Take(eH_PDS,-9))<2.");
  postCutProcess(12, "HC3", dF15, 0);
  auto dF14 = dF15.Filter(filterFunHCALtopo, { "eH_PDS" });
  postCutProcess(13, "HC0 topo", dF14, 0);
  auto dF16 = dF14.Filter(filterFunMuonHC2, { "eH_PDS" });
  postCutProcess(14, "HC2 muon cut", dF16, 0);
  auto dF17 = dF16.Filter(filterFunMuonHC1, { "eH_PDS" });
  postCutProcess(15, "HC1 muon cut", dF17, 0);
  auto dF18 = dF17.Filter(filterFunMuonHC0, { "eH_PDS" });
  postCutProcess(16, "HC0 muon cut", dF18, 0);
  auto dF19 = dF18.Filter("DIMU_NDimuonsTrue == 1");
  postCutProcess(17, "NDimuonsCut == 1", dF19, 0);
  auto dF20 = dF19.Filter("DIMU_vertexZ > 2980.96");
  postCutProcess(18, "DIMU_BeforeECAL cut ", dF20, 0);
  auto dF21 = dF20.Filter("DIMU_vertexZ < 3430.96");
  postCutProcess(19, "DIMU_AfterECAL cut ", dF21, 0);
  auto dF22 = dF21.Filter(filterFunVETO, { "VETO" });
  postCutProcess(20, "VETO < 0.05", dF22, 0);
  auto dF23 = dF22.Filter(filterFunVETO2, { "VETO" });
  postCutProcess(21, "Veto23 > 0.01", dF23, 0);
  auto dF24 = dF23.Filter(filterMagnetIronBlock, { "Strawy11", "Strawx12", "Strawy12" });
  postCutProcess(22, "Magnet Acceptence (Iron Cut)", dF24, 0); 
  auto dF25 = dF24.Filter(filterLargeStrawMult, { "Strawx11", "Strawy11" });
  postCutProcess(23, "Straw Acceptence Cut", dF25, 0);
  auto dF26 = dF25.Filter(filterStraw11EachSide, { "Strawx11", "Strawy11" });
  postCutProcess(24, "Straw11 +- cut", dF26, 0);
  auto dF27 = dF24.Filter(filterSingleMuonST, { "Strawx11", "Strawy11" });
  postCutProcess(25, "Straw Single Muon Cut", dF27, 0);
  auto dF28 = dF24.Filter(filterSingleMuons, { "Straw12_PDGID" });
  postCutProcess(26, "Straw11 single muons cut [TRUE]", dF28, 0);
  auto dF29 = dF24.Filter(filterDimuons, { "Straw12_PDGID"});
  postCutProcess(27, "Straw11 dimuon seperation cut [TRUE]", dF29, 0);

/*
  auto dF24 = dF23.Filter(filterMagnetIronBlock, { "Strawy11_truth", "Strawx12_truth", "Strawy12_truth" });
  postCutProcess(22, "Magnet Iron Cut", dF24, 0);
  auto dF25 = dF24.Filter(filterAcceptenceST, { "Straw11_PDGID", "Straw12_PDGID" });
  postCutProcess(23, "Straw Acceptence Cut", dF25, 0);
  auto dF26 = dF24.Filter(filterStrawFalseMuons, { "Straw11_PDGID" });
  postCutProcess(24, "Straw11 false muons cut", dF26, 0);
  auto dF27 = dF24.Filter(filterSingleMuons, { "Straw11_PDGID" });
  postCutProcess(25, "Straw11 single muons cut", dF27, 0);
  auto dF28 = dF24.Filter(filterStrawBackground, { "Straw11_PDGID" });
  postCutProcess(26, "Straw11 background cut", dF28, 0);

  auto dF27 = dF25.Filter(filterStraw11EachSide, { "Strawx11", "Strawy11" });
  postCutProcess(25, "Straw11 +- cut", dF26, 0);
  auto dF28 = dF26.Filter(filterLargeStrawMult, { "Strawx12", "Strawy12" });
  postCutProcess(26, "Straw12 mult == 2 cut", dF27, 0);
  auto dF29 = dF27.Filter(filterStraw11EachSide, { "Strawx12", "Strawy12" });
  postCutProcess(27, "Straw12 +- cut", dF28, 0);
  auto dF30 = dF25.Filter(filterStraw11MuonSeperation, { "Strawx11_truth", "Straw11_PDGID"});
  postCutProcess(27, "Straw11 dimuon seperation cut", dF30, 0);
  auto dF31 = dF25.Filter(filterStraw11Region, { "Strawx11_truth", "Strawx12_truth"});
  postCutProcess(28, "Straw11 dimuon region cut", dF31, 0);
  postCutProcess(29, "ALL dimuons", d1, 0);
  postCutProcess(30, "After Blind", dBlind, 0);
  //postCutProcess(31, "Without iron magnet", dMagnet_truth, 0);
  
  //auto dF27 = dF26.Filter(filterLargeStrawMult, { "Strawx11", "Strawy11" });
  //postCutProcess(26, "Straw11 region cut", dF27, 0);
  //auto dF26 = dF25.Filter(filterStraw12, { "Strawx12", "Strawy12" });
  //postCutProcess(24, "Straw12 cut", dF26, 0);
*/
// ------------------------------------------------------------------------------------------------------------ //
  
  auto dBlind_ord = d1.Filter("eH0_blind > 2 && eH1_blind > 2 && eH2_blind > 2");
  postCutProcess(0, "AllEvents",dBlind_ord, 1);
  auto dF1_ord = dBlind_ord.Filter(physTrigger, { "ECALpenergy", "ECALenergy" });
  postCutProcess(1, "EC < 76 && PS > 0.4", dF1_ord, 1);
  auto dF3_ord = dF1_ord.Filter(filterFunSRD, { "SRD" });
  postCutProcess(2, "SRD ", dF3_ord, 1);
  auto dF4_ord = dF3_ord.Filter(filterFunTrack, { "BestMomupGenfit", "BestMomupChisq" });
  postCutProcess(3, "Track Quality", dF4_ord, 1);
  auto dF5_ord = dF4_ord.Filter("abs(BestMomupGenfit-100)<10");
  postCutProcess(4, "Momentum", dF5_ord, 1);
  auto dF6_ord = dF5_ord.Filter("abs(InAngle)<3.");
  postCutProcess(5, "Angle12", dF6_ord, 1);
  auto dF7_ord = dF6_ord.Filter(filterStrawMult, { "mpStraw3X", "mpStraw3Y" });
  postCutProcess(6, "Straw3 Multiplicity Cut", dF7_ord, 1);
  auto dF72_ord = dF7_ord.Filter("Sum(VHCALenergy) < 1.5");
  postCutProcess(7, "VHCAL cut", dF72_ord, 1);
  auto dF8_ord = dF72_ord.Filter(filterStrawMult, { "mpStraw4X", "mpStraw4Y" });
  postCutProcess(8, "Straw4 Multiplicity Cut", dF8_ord, 1);
  auto dF9_ord = dF8_ord.Filter(filterFunECAL23, { "ECALenergy" }); // Check if correct!
  postCutProcess(9, "cutEC(2,2) ", dF9_ord, 1); // CENTER ECAL CELL IS INDEED (2,2)
  auto dF10_ord = dF9_ord.Filter(filterFunEC, { "ECALpenergy", "ECALenergy" });
  postCutProcess(10, "EC ", dF10_ord, 1);
  auto dF11_ord = dF10_ord.Filter(filterECALratio, { "ECALenergy" });
  postCutProcess(11, "ECratio ", dF11_ord, 1);
  //auto dF12_ord = dF11_ord.Filter(filterLargeStrawMult, { "Strawx11", "Strawy11" });
  //postCutProcess(12, "Straw11 mult == 2 cut", dF12_ord, 1);
  //auto dF13_ord = dF12_ord.Filter(filterStraw11EachSide, { "Strawx11", "Strawy11" });
  //postCutProcess(13, "Straw11 cut", dF13_ord, 1);
  //auto dF13_ord = dF12_ord.Filter(filterLargeStrawMult, { "Strawx12", "Strawy12" });
  //postCutProcess(13, "Straw12 mult == 2 cut", dF13_ord, 1);
  //auto dF14_ord = dF13_ord.Filter(filterStraw12, { "Strawx12", "Strawy12" });
  //postCutProcess(14, "Straw12 cut", dF14_ord, 1);
  auto dF15_ord = dF11_ord.Filter("Sum(Take(eH_PDS,-9))<2.");
  postCutProcess(12, "HC3", dF15_ord, 1);
  auto dF16_ord = dF15_ord.Filter(filterFunHCALtopo, { "eH_PDS" });
  postCutProcess(13, "HC0 topo", dF16_ord, 1);
  auto dF17_ord = dF16_ord.Filter(filterFunMuonHC2_center, { "eH_PDS" });
  postCutProcess(14, "HC2 muon cut", dF17_ord, 1);
  auto dF18_ord = dF17_ord.Filter(filterFunMuonHC1_center, { "eH_PDS" });
  postCutProcess(15, "HC1 muon cut", dF18_ord, 1);
  auto dF19_ord = dF18_ord.Filter(filterFunMuonHC0_center, { "eH_PDS" });
  postCutProcess(16, "HC0 muon cut", dF19_ord, 1);
  auto dF20_ord = dF19_ord.Filter("DIMU_NDimuonsTrue == 1");
  postCutProcess(17, "NDimuonsCut == 1", dF20_ord, 1);
  auto dF21_ord = dF20_ord.Filter("DIMU_vertexZ > 2980.96");
  postCutProcess(18, "DIMU_BeforeECAL cut ", dF21_ord, 1);
  auto dF22_ord = dF21_ord.Filter("DIMU_vertexZ < 3430.96");
  postCutProcess(19, "DIMU_AfterECAL cut ", dF22_ord, 1);
  auto dF23_ord = dF22_ord.Filter(filterFunVETO, { "VETO" });
  postCutProcess(20, "VETO < 0.05", dF23_ord, 1);
  auto dF24_ord = dF23_ord.Filter(filterFunVETO2, { "VETO" });
  postCutProcess(21, "Veto23 > 0.01", dF24_ord, 1);
  auto dF25_ord = dF24_ord.Filter(filterMagnetIronBlock, { "Strawy11", "Strawx12", "Strawy12" });
  postCutProcess(22, "Magnet Acceptence (Iron Cut)", dF25_ord, 1);
  auto dF26_ord = dF25_ord.Filter(filterLargeStrawMult, { "Strawx11", "Strawy11" });
  postCutProcess(23, "Straw Acceptence Cut", dF26_ord, 1);
  auto dF27_ord = dF26_ord.Filter(filterStraw11EachSide, { "Strawx11", "Strawy11" });
  postCutProcess(24, "Straw11 +- cut", dF27_ord, 1);
  //auto dF24 = dF23.Filter(filterStraw11, { "Strawx11", "Strawx11" });
  //postCutProcess(22, "Straw11 mult == 2 cut", dF24);

  double bias = 15.;
  double EOT = 884630574.;
  std::size_t size = allN.size();

  cout << "TotalEvents: " << (allN[0].GetValue()) << endl;

  cout << endl;
  cout << "Ordered cutflow" << endl;
  printf("Cut:|CutName\t\t\t   |NumEvents\t\t|Eff in %%\t\t|EffRel in %%\t\t\n");
  for (int ii = 1; ii < allN.size(); ii++) {
    double eff = (allN[ii - 1].GetValue() - allN[ii].GetValue()) / (EOT * bias);
    double effRel = (allN[ii - 1].GetValue() - allN[ii].GetValue()) / (1. * allN[ii - 1].GetValue());
    eff = eff * 100.;
    effRel = effRel * 100.;

    printf("%.2i: |%-30s|%llu\t\t|%.7f\t\t|%.7f\n", ii, allCutName[ii].c_str(), allN[ii].GetValue(), eff, effRel);
  }
  cout << endl;
  cout << "Simulated EOTs (EOT): " << EOT << std::endl;
  cout << "Total number of EOTs (EOT * BIAS): " << EOT * bias << std::endl;
  cout << "Total efficency: " << *allN[size-1] / (EOT * bias)  << " +- " << sqrtf(*allN[size-1]) / (EOT * bias) <<  std::endl; 
  //cout << "----------------------  ACCEPTENCE  ------------------" << std::endl;
  //cout << "Acceptence cut: " << *allN[size-6] / (EOT * bias)  << " +- " << sqrtf(*allN[size-6]) / (EOT * bias) <<  std::endl;
  //cout << "----------------------  DIMUON CUT  ------------------" << std::endl; 
  //cout << "Efficiency Straw11 mult: " << *allN[size-1]*1. / (*allN[size-6]*1.) << " +- " << sqrtf(*allN[size-1])*1. / (*allN[size-6]*1.) <<  std::endl;
  //cout << "----------------------  TOTAL EFFICIENCIES  ------------------" << std::endl;
  //cout << "Total efficency: " << *allN[size-1] / (EOT * bias)  << " +- " << sqrtf(*allN[size-1]) / (EOT * bias) <<  std::endl; 

  /*cout << "----------------------  CUT EFFICIENCIES  ------------------" << std::endl;
  cout << "Efficiency Straw11 mult: " << *allN[size-4]*1. / (*allN[size-5]*1.) << " +- " << sqrtf(*allN[size-4])*1. / (*allN[size-5]*1.) <<  std::endl;
  cout << "Efficiency Straw11 +- cut: " << *allN[size-3]*1. / (*allN[size-5]*1.)  << " +- " << sqrtf(*allN[size-3])*1. / (*allN[size-5]*1.) <<  std::endl;
  cout << "Efficiency Straw12 mult cut: " << *allN[size-2]*1. / (*allN[size-5]*1.)  << " +- " << sqrtf(*allN[size-2])*1. / (*allN[size-5]*1.) <<  std::endl;
  cout << "Efficiency Straw12 +- cut: " << allN.back().GetValue()*1. / (*allN[size-5]*1.) << " +- " << sqrtf(allN.back().GetValue())*1. / (*allN[size-5]*1.) <<  std::endl;

  cout << "----------------------  TOTAL EFFICIENCIES  ------------------" << std::endl;
  cout << "Final eff (cut 1): " << *allN[size-4] / (EOT * bias)  << " +- " << sqrtf(*allN[size-4]) / (EOT * bias) <<  std::endl;
  cout << "Final eff (cut 2): " << *allN[size-3] / (EOT * bias)  << " +- " << sqrtf(*allN[size-3]) / (EOT * bias) <<  std::endl;
  cout << "Final eff (cut 3): " << *allN[size-2] / (EOT * bias)  << " +- " << sqrtf(*allN[size-2]) / (EOT * bias) <<  std::endl;
  cout << "Final eff (cut 4): " << allN.back().GetValue() / (EOT * bias)  << " +- " << sqrtf(allN.back().GetValue()) / (EOT * bias) <<  std::endl;*/

  cout << endl;
  cout << "Ordered cutflow" << endl;
  printf("Cut:|CutName\t\t\t   |NumEvents\t\t|Eff in %%\t\t|EffRel in %%\t\t\n");
  for (int ii = 1; ii < allN_ord.size(); ii++) {
    double eff = (allN_ord[ii - 1].GetValue() - allN_ord[ii].GetValue()) / (EOT * bias);
    double effRel = (allN_ord[ii - 1].GetValue() - allN_ord[ii].GetValue()) / (1. * allN_ord[ii - 1].GetValue());
    eff = eff * 100.;
    effRel = effRel * 100.;

    printf("%.2i: |%-30s|%llu\t\t|%.7f\t\t|%.7f\n", ii, allCutName_ord[ii].c_str(), allN_ord[ii].GetValue(), eff, effRel);
  }

  cout << endl;
  cout << "Simulated EOTs (EOT): " << EOT << std::endl;
  cout << "Total number of EOTs (EOT * BIAS): " << EOT * bias << std::endl;
  cout << "  ------------------------------------------  " << std::endl;
  cout << "Final eff: " << allN_ord.back().GetValue() / (EOT * bias)  << std::endl;

  //Write to the output
  TFile *fout = new TFile(fnameOut.c_str(), "recreate");
  fout -> mkdir("Dimuons");
  fout -> cd("Dimuons");
  for (auto h : allH) {
    h->Write();
  }
  fout -> mkdir("Dimuons_ord");
  fout -> cd("Dimuons_ord");
  for (auto h : allH_ord) {
    h->Write();
  }
  fout->Close();
}
