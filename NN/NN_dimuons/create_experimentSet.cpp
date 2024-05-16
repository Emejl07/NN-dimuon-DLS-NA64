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

// Simple event struct
struct Event {
  UInt_t number;
  bool isEqual(UInt_t other) const { return (number == other); }
};

// Structure to place bad events in spills
struct EventsInSpill {
  UInt_t run;
  UInt_t spill;
  std::vector<Event> events;

  // method to check if a given event is in this spill interval
  bool contains(UInt_t runN, UInt_t spillN, UInt_t spilleventN) const
  {
    // Check if run matches
    if (runN != run) return false;

    // Check if spill matches
    if (spillN != spill) return false;

    // Loop over all events in struct and check if matches
    for (auto event : events) {
      if (event.isEqual(spilleventN)) return true;
    }
    return false;
  }
};

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

//VHCAL[1][4][4]
int vhcalID(int iplane, int ix, int iy) {
  int ret = 0;
  ret = iplane * 16 + ix * 4 + iy;
  return ret;
}

ROOT::RVec<double> HcalPedestal(RVecD &ene) {

  const double pedestal = 0.2;

  ROOT::RVec<double> HCAL_ene_out;
  for (int is = 0; is < 4; is++) {
    for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 3; iy++) {
        auto id = hcalID(is, ix, iy);
        if (ene[id] > pedestal)
          HCAL_ene_out.push_back(ene[id]);
        else
          HCAL_ene_out.push_back(0.);
      }
    }
  }
  return HCAL_ene_out;
}

/*The following function is used to define in-time energy for ECAL cells.
 */
ROOT::RVec<double> EcalEneInTime(RVecD &ene, RVecD &t0, RVecD &tE, RVecD &tSigma) {

  const int nsigma = 3.;

  ROOT::RVec<double> ECAL_eneT;
  for (int is = 0; is < 2; is++) {
    for (int ix = 0; ix < 5; ix++) {
      for (int iy = 0; iy < 6; iy++) {
        auto id = ecalID(is, ix, iy);
        auto deltaT = t0[id] - tE[id];
        if (fabs(deltaT) < nsigma * tSigma[id])
          ECAL_eneT.push_back(ene[id]);
        else
          ECAL_eneT.push_back(0.);
      }
    }
  }
  return ECAL_eneT;
}

/*The following function is used to define in-time energy for HCAL cells.
 */
ROOT::RVec<double> HcalEneInTime(RVecD &ene, RVecD &t0, RVecD &tE, RVecD &tSigma) {

  const int nsigma = 3.;

  ROOT::RVec<double> HCAL_eneT;
  for (int is = 0; is < 4; is++) {
    for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 3; iy++) {
        auto id = hcalID(is, ix, iy);
        auto deltaT = t0[id] - tE[id];
        if (fabs(deltaT) < nsigma * 5.5)
          HCAL_eneT.push_back(ene[id]);
        else
          HCAL_eneT.push_back(0.);
      }
    }
  }
  return HCAL_eneT;
}

/*The following function is used to define in-time energy for VHCAL cells.
 */
ROOT::RVec<double> VHcalEneInTime(RVecD &ene, RVecD &t0, RVecD &tE, RVecD &tSigma) {

  const int nsigma = 3.;

  ROOT::RVec<double> VHCAL_eneT;
  for (int is = 0; is < 1; is++) {
    for (int ix = 0; ix < 4; ix++) {
      for (int iy = 0; iy < 4; iy++) {
        auto id = vhcalID(is, ix, iy);
        auto deltaT = t0[id] - tE[id];
        if (fabs(deltaT) < nsigma * tSigma[id])
          VHCAL_eneT.push_back(ene[id]);
        else
          VHCAL_eneT.push_back(0.);
      }
    }
  }
  return VHCAL_eneT;
}


/*the following function is used to define in-time energy for srd cells.
 */
ROOT::RVec<double> SRDEneInTime(RVecD &ene, RVecD &t0, RVecD &tE, RVecD &tSigma) {

  const int nsigma = 3.;
  //const double corrFactor[3] = {0.5394, 0.5512, 0.7489};

  ROOT::RVec<double> SRD_eneT;
  for (int ix = 0; ix < 3; ix++) {
    auto deltaT = t0[ix] - tE[ix];
    if (fabs(deltaT) < nsigma * tSigma[ix])
      //SRD_eneT.push_back(ene[ix]/corrFactor[ix]);
      SRD_eneT.push_back(ene[ix]);
    else
      SRD_eneT.push_back(0.);
  }
  return SRD_eneT;
}

/*The following function is used to define in-time energy for VETO cells.
 */
ROOT::RVec<double> VetoEneInTime(RVecD &ene, RVecD &t0, RVecD &tE, RVecD &tSigma) {

  const int nsigma = 3.;

  ROOT::RVec<double> VETO_eneT;
  for (int ix = 0; ix < 6; ix++) {
    auto deltaT = t0[ix] - tE[ix];
    if (fabs(deltaT) < nsigma * tSigma[ix/2])
      VETO_eneT.push_back(ene[ix]);
    else
      VETO_eneT.push_back(0.);
  }
  return VETO_eneT;
}

std::vector<EventsInSpill> badSpills;

/*
 * Filter bad events in spills
 */
bool filterBadSpills(UInt_t runN, UInt_t spillN, UInt_t spilleventN) {
  for (auto spill : badSpills) {
    if (spill.contains(runN, spillN, spilleventN)) return false;
  }

  return true;
}

// FILTER FUNCTIONS --------------------------------------------------
//Trigger time filter
bool filterTriggerTime(double &masterT) {
  bool ret = ((masterT >= 86.5) && (masterT <= 113.5));
  return ret;
}

//Straw multiplicity cut
bool filterStrawMult(int &straw3X, int &straw3Y) {
  bool ret;
  ret = ((straw3X > 0) && (straw3X <= 5) && (straw3Y > 0) && (straw3Y <= 5));
  return ret;
}

//Centermost cell cut: maximum in-time energy of ECAL only is in center cell
bool filterECALcenter(RVecD &ECAL_eneT) {
  double eMax = 0;
  int xMax, yMax;
  for (int ix = 0; ix < 5; ix++) {
    for (int iy = 0; iy < 6; iy++) {
      auto id = ecalID(1, ix, iy);
      if (ECAL_eneT[id] > eMax) {
        eMax = ECAL_eneT[id];
        xMax = ix;
        yMax = iy;
      }
    }
  }
  return ((xMax == 2) && (yMax == 2));
}

/*SRD cut. From Toropin:
 Srdtok[ix] = (SRDsignal[ix] > 1.0) && (SRDsignal[ix] < 80) && intime, ix = 0,1,2.
 SRD_OK = (srdtok[0] && srdtok[1] &&  srdtok[2])
 */
bool filterSRD(RVecD &SRD_eneT) {
  bool srdtok[3];

  for (int ii = 0; ii < 3; ii++) {
    srdtok[ii] = (SRD_eneT[ii] > 0.001) && (SRD_eneT[ii] < 0.08);
  }

  return (srdtok[0] && srdtok[1] && srdtok[2]);
}

/*VETO cut. Toropin: sum veto pmts < 4.0 MeV*/
/*MIP MPV ~ 10 MeV*/
/*2*MIP in sum ~ 20 MeV*/
bool filterVETO(RVecD &VETO_eneT){
  // TODO: Cut Value to be checked!
  const double vetoCut01 = 0.0076;
  const double vetoCut23 = 0.008;
  const double vetoCut45 = 0.0076;

  double eVETO01(0), eVETO23(0), eVETO45(0);
  // VETO 0
  eVETO01 += VETO_eneT[0]; eVETO01 += VETO_eneT[1];
  // VETO 1
  eVETO23 += VETO_eneT[2]; eVETO23 += VETO_eneT[3];
  // VETO 2
  eVETO45 += VETO_eneT[4]; eVETO45 += VETO_eneT[5];

  return ((eVETO01 < vetoCut01) && (eVETO23 < vetoCut23) && (eVETO45 < vetoCut45));
}

/*ECcuts: EcalBadEnergy<30 GeV and EcalTotal<150 and PreShower>0.7
 * EcalTotal: all energies (in-time,out-of-time)
 * EcalBadEnergy: total - in_time
 */
bool filterECALenergy(RVecD &ECAL_ene, RVecD &ECAL_eneT) {
  double sECAL_ene = 0;
  double sECAL_eneT = 0;
  double sPS_ene = 0;
  double sPS_eneT = 0;
  int xMax, yMax;
  for (int is = 0; is < 2; is++) {
    for (int ix = 0; ix < 5; ix++) {
      for (int iy = 0; iy < 6; iy++) {
        auto id = ecalID(is, ix, iy);
        if (is == 0) { // PS energy
          sPS_ene += ECAL_ene[id];
          sPS_eneT += ECAL_eneT[id];
        }
        sECAL_ene += ECAL_ene[id];
        sECAL_eneT += ECAL_eneT[id];
      }
    } 
  } 
  return ((sECAL_ene < 120.) && (sECAL_ene - sECAL_eneT < 30.) && (sPS_ene > 0.7));
} 

/* HCcuts:
 * SumHC012-SumHC2 < 50.; where:
 SumHC012 is the total energy in HCAL0+HCAL1+HCAL2 (no time cut)
 SumHC2 is the total energy in HCAL0+HCAL1+HCAL2 with time cut on individual cells
 below is done this->SumHC2 is the total energy in HCAL0+HCAL1+HCAL2 with time cut on individual cells TO BE CHECKED
 *
 */
bool filterHCALOOTenergy(RVecD &HCAL_ene, RVecD &HCAL_eneT) {
  double SumHC = 0;
  double SumHC_T = 0;
  for (int is = 0; is < 3; is++) {
    for (int ix = 0; ix < 3; ix++) {
      for (int iy = 0; iy < 3; iy++) {
        auto id = hcalID(is, ix, iy);
        //if (is < 3) {//I also sum HCAL2, otherwise the cut above has no meaning in case of large energy deposition in HCAL2 (in time)-> modified
        if (is < 2) { //HCAL2 not included as Alexander
          SumHC_T += HCAL_eneT[id];
        }
        SumHC += HCAL_ene[id];
      }
    }
  }
  return (SumHC - SumHC_T < 50.);
}


//track
bool filterTrackQuality(double &mom_genfit_upMM, double &pvalue_genfit_upMM) {
  return (mom_genfit_upMM > 0. && mom_genfit_upMM < 9999. && pvalue_genfit_upMM > 0.01);
}

/* ECratio cut
 * (5x6 - 3x3)/(5x6) < 0.06
 * Center is in cell (2,2)
 */
bool filterECALratio(RVecD &ECAL_eneT) {
  double SumAll(0), Sum9(0);
  for (int is = 0; is < 2; is++) {
    for (int ix = 0; ix < 5; ix++) {
      for (int iy = 0; iy < 6; iy++) {
        auto id = ecalID(is, ix, iy);
        SumAll += ECAL_eneT[id];
        if ((ix >= 1) && (ix <= 3) && (iy >= 1) && (iy <= 3))
          Sum9 += ECAL_eneT[id];
      }
    }
  }
  if (SumAll <= 0)
    return false;
  return ((SumAll - Sum9) / SumAll < .06);
}

/* HCtopo
 * Reject event if (NumHCcells ≥ 1 (Ecell > 0.5) && EnergyHC – E(1,1) > E(1,1)), (in HC0), E(1,1) – central cell.
 *
 */
bool filterHCtopo(RVecD &HCAL_eneT) {
  int nCells = 0;
  double thr = 0.5;
  double HC0 = 0;
  for (int ix = 0; ix < 3; ix++) {
    for (int iy = 0; iy < 3; iy++) {
      auto id = hcalID(0, ix, iy);
      HC0 += HCAL_eneT[id];
      if (HCAL_eneT[id] > thr) //toropin 40 adc
        nCells++;
    }
  } 
  if (nCells >= 1) {
    return (HC0 < 2 * HCAL_eneT[hcalID(0, 1, 1)]);  
  } 
  return true;
}

int mpStrawHits(vector<double> &strawX){
  return strawX.size();
}

int strawSelection(vector<double> &straw11X){
  int mpST11 = straw11X.size();
 
  if(mpST11 == 0) return -30; // No hit
  
  if(mpST11 == 1) { // Single hit
    if (straw11X[0] > -147) return -10; // mu-
    else return 10; // mu+
  }
  if(mpST11 == 2){ // two hit
    if (straw11X[0] > -147 && straw11X[1] < -147 || straw11X[0] < -147 && straw11X[1] > -147) return 20;
    else return -20;
  }
  else return 30; // More than 2 hits

}

void create_experimentSet(string fIn = "", string fOut = "") {
   
    //
    // DIMUON EVENTS DATA
    //

    std::string fnameIn = (fIn == "") ? "/afs/cern.ch/work/e/ezaya/public/data_partblinded_pass5_small/job008747.root" : fIn;  // Data set of Dimuons
    std::string fnameOut = (fOut == "") ? "experimentSet_out.root" : fOut;

    // Define the RDataFrame with the input ROOT file
    ROOT::RDataFrame d0("tout", fnameIn.c_str());

    // PHYSICS TRIGGER
    auto dTrigger = d0.Filter("triggerSources & 4");
    // Skip bad events
    auto dRun = dTrigger.Filter(filterBadSpills, {"runN", "spillN", "spilleventtN"});
    //Define columns
    auto d2 = dRun.Define("eecal", "eecal0+eecal1");
    d2 = d2.Define("ehcal", "ehcal0+ehcal1");
    d2 = d2.Define("chi2_genfit", "chi2[0]");
    d2 = d2.Define("chi2_mm34", "chi2[1]");
    d2 = d2.Define("HCAL_ene_PDS", HcalPedestal , { "HCAL_ene" });
    d2 = d2.Define("ECAL_eneT",     EcalEneInTime,  { "ECAL_ene", "ECAL_t0", "ECAL_tE", "ECAL_tSigma" });
    d2 = d2.Define("HCAL_eneT",     HcalEneInTime,  { "HCAL_ene_PDS", "HCAL_t0", "HCAL_tE", "HCAL_tSigma" });
    d2 = d2.Define("VHCAL_eneT",    VHcalEneInTime, { "VHCAL_ene", "VHCAL_t0", "VHCAL_tE", "VHCAL_tSigma" });
    d2 = d2.Define("SRD_eneT",      SRDEneInTime,   { "SRD_ene", "SRD_t0", "SRD_tE", "SRD_tSigma" });
    d2 = d2.Define("VETO_enePMT_T", VetoEneInTime,  { "VETO_enePMT", "VETO_t0", "VETO_tE", "VETO_tSigma" });
    d2 = d2.Define("ECAL_OOTratio", "(Sum(ECAL_ene) - Sum(ECAL_eneT))/(Sum(ECAL_eneT))");
    d2 = d2.Define("energy", "Sum(ECAL_eneT) + Sum(Take(HCAL_eneT,18))");
    d2 = d2.Define("energyDiff", "abs(mom_genfit_upMM - energy)");
    d2 = d2.Define("HCAL0", "Sum(Take(HCAL_eneT,9))"); //take first 9 HC0 distribution
    d2 = d2.Define("HCAL1", "Sum(Take(HCAL_eneT,18)) - Sum(Take(HCAL_eneT,9))"); //take first 9 HC1 distribution
    d2 = d2.Define("HCAL2", "Sum(Take(HCAL_eneT,27)) - Sum(Take(HCAL_eneT,18))"); //take first 9 HC2 distribution
    d2 = d2.Define("HCAL012", "Sum(Take(HCAL_eneT,27))");
    d2 = d2.Define("eH0_11", "HCAL_eneT[4]");
    d2 = d2.Define("eH1_11", "HCAL_eneT[13]");
    d2 = d2.Define("eH2_11", "HCAL_eneT[22]");
    d2 = d2.Define("ECAL", "Sum(ECAL_eneT)");
    d2 = d2.Define("mpST11", mpStrawHits, { "ST12x" });
    d2 = d2.Define("mpST12", mpStrawHits, { "ST11x" });
    d2 = d2.Define("strawSelection", strawSelection, { "ST12x" });

    auto dF1 = d2.Filter(filterTriggerTime, { "masterT" });
    auto dF2 = dF1.Filter(filterSRD, { "SRD_eneT" });
    auto dF3 = dF2.Filter(filterTrackQuality, { "mom_genfit_upMM", "pvalue_genfit_upMM" });
    auto dF4 = dF3.Filter("abs(mom_genfit_upMM-100)<10");
    auto dF5 = dF4.Filter("abs(inangle_genfit_upMM)<3.");
    auto dF6 = dF5.Filter(filterStrawMult, { "mpStraw3X", "mpStraw3Y" });
    auto dF7 = dF6.Filter("Sum(VHCAL_ene) < 1.5");
    auto dF8 = dF7.Filter(filterStrawMult, { "mpStraw4X", "mpStraw4Y" });
    auto dF9 = dF8.Filter(filterECALcenter, { "ECAL_eneT" });
    auto dF10 = dF9.Filter(filterECALenergy, { "ECAL_ene", "ECAL_eneT" });
    auto dF11 = dF10.Filter(filterECALratio, { "ECAL_eneT" });
    auto dF12 = dF11.Filter("chi2_genfit<100");  // Set very large chi2 cut to study dimuon response
    auto dF14 = dF12.Filter(filterHCALOOTenergy, { "HCAL_ene_PDS", "HCAL_eneT" });
    auto dF15 = dF14.Filter("Sum(Take(HCAL_ene_PDS,-9))<2.");
    auto dF16 = dF15.Filter(filterHCtopo, { "HCAL_eneT"} );




    // Take a snapshot of the processed dataframe
    ROOT::RDF::RSnapshotOptions opts;
    opts.fMode = "update";
    opts.fOverwriteIfExists = true;
    dF16.Snapshot("training_set", fnameOut.c_str(), {"HCAL012", "HCAL0", "HCAL1", "HCAL2", "ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection" }, opts);

}
