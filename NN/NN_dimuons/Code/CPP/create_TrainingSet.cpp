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

bool filterMagnetIronBlock(vector<double> &straw11Y, vector<double> &straw12X, vector<double> &straw12Y) {
  for(auto y : straw11Y){
    if(!(y > -120. && y < 120.)) return 0;
  }
  for(auto y : straw12Y){
    if(!(y > -150. && y < 150.)) return 0;
  }
  for(auto x : straw12X){
    for(auto y : straw12Y){
      if(((x > -130. && x < 60.) && (y > -150. && y < -100.))) return 0;
    }
  }
  return 1;
}

int mpStrawHits(vector<double> &strawX){
  return strawX.size();
}

int strawSelection(vector<double> &straw11X){
  int mpST11 = straw11X.size();
 
  if(mpST11 == 0) return -30; // No hit
  
  if(mpST11 == 1) { // Single hit
    if (straw11X[0] > -150) return -10; // mu-
    else return 10; // mu+
  }
  if(mpST11 == 2){ // two hit
    if (straw11X[0] > -150 && straw11X[1] < -150 || straw11X[0] < -150 && straw11X[1] > -150) return 20;
    else return -20;
  }
  else return 30; // More than 2 hits

}

void create_TrainingSet(string fIn = "", string fOut = "") {
   
    //
    // DIMUON EVENTS DATA
    //

    //std::string fnameIn = (fIn == "") ? "/eos/user/e/ezaya/simulation_output/DIMUONS/data/09042024_TrueDetector_beforeMagnet.root" : fIn;  // Data set of Dimuons 01052024_Mag185Mag215
    std::string fnameIn = (fIn == "") ? "/eos/user/e/ezaya/simulation_output/DIMUONS/data/01052024_Mag185Mag215.root" : fIn;  // Data set of Dimuons 01052024_Mag185Mag215
    std::string fnameOut = (fOut == "") ? "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_dimuon_out.root" : fOut;

    // Define the RDataFrame with the input ROOT file
    ROOT::RDataFrame d0("ana_tree", fnameIn);

    // Make all the necessary filters
    auto dFiltered = d0.Filter(("DIMU_NDimuonsTrue < 2")); // Save single dimuon event or no dimuon event

    auto dMagnet = dFiltered.Filter(filterMagnetIronBlock, { "Strawy11", "Strawx12", "Strawy12" });

    // Define the dataframe transformations
    auto dOut_dimuon = dMagnet.Define("IsDimuon", "DIMU_NDimuonsTrue")
                 .Define("eH_PDS", HcalPedestral, {"HCALenergy"})
                 .Define("HCAL0", "Sum(Take(eH_PDS, 9))")  // take first 9 HC0 distribution
                 .Define("HCAL1", "Sum(Take(eH_PDS, 18)) - Sum(Take(eH_PDS, 9))")  // take first 9 HC1 distribution
                 .Define("HCAL2", "Sum(Take(eH_PDS, 27)) - Sum(Take(eH_PDS, 18))")  // take first 9 HC2 distribution
                 .Define("HCAL012", "Sum(Take(eH_PDS, 27))")
                 .Define("ECAL", "Sum(ECALenergy)")
                 .Define("eH0_11", "eH_PDS[4]")
                 .Define("eH1_11", "eH_PDS[13]")
                 .Define("eH2_11", "eH_PDS[22]")
                 .Define("mpST11", mpStrawHits, { "Strawx11" })
                 .Define("mpST12", mpStrawHits, { "Strawx12" })
                 .Define("strawSelection", strawSelection, { "Strawx11" })
                 .Define("veto01", "VETO[0]")
                 .Define("veto23", "VETO[1]")
                 .Define("veto45", "VETO[2]");

    // Take a snapshot of the processed dataframe
    ROOT::RDF::RSnapshotOptions opts;
    opts.fMode = "update";
    opts.fOverwriteIfExists = true;
    dOut_dimuon.Snapshot("training_set", fnameOut.c_str(), {"IsDimuon", "HCAL012", "HCAL0", "HCAL1", "HCAL2", "ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection", "veto01", "veto23", "veto45"  }, opts);

    //
    // PION EVENTS DATA
    //

    std::string fnameIn_pion = (fIn == "") ? "/eos/user/e/ezaya/simulation_output/PION/data/pion-2023A.root" : fIn; // "./../old_simulations/Dimuons250224_15.root" : fIn; // Data set of Dimuons
    std::string fnameOut_pion = (fOut == "") ? "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_pion_out.root" : fOut;

    // Define the RDataFrame with the input ROOT file
    ROOT::RDataFrame d0_pion("ana_tree", fnameIn_pion);

    // Make all the necessary filters
    //auto dMuonHit = d0.Filter("CDIMU_PDGID == 13") // Remove 

    auto dMagnet_pion = d0_pion.Filter(filterMagnetIronBlock, { "Strawy11", "Strawx12", "Strawy12" });

    // Define the dataframe transformations
    auto dOut_pion = dMagnet_pion.Define("IsDimuon", "DIMU_NDimuonsTrue")
                 .Define("eH_PDS", HcalPedestral, {"HCALenergy"})
                 .Define("HCAL0", "Sum(Take(eH_PDS, 9))")  // take first 9 HC0 distribution
                 .Define("HCAL1", "Sum(Take(eH_PDS, 18)) - Sum(Take(eH_PDS, 9))")  // take first 9 HC1 distribution
                 .Define("HCAL2", "Sum(Take(eH_PDS, 27)) - Sum(Take(eH_PDS, 18))")  // take first 9 HC2 distribution
                 .Define("HCAL012", "Sum(Take(eH_PDS, 27))")
                 .Define("ECAL", "Sum(ECALenergy)")
                 .Define("eH0_11", "eH_PDS[4]")
                 .Define("eH1_11", "eH_PDS[13]")
                 .Define("eH2_11", "eH_PDS[22]")
                 .Define("mpST11", mpStrawHits, { "Strawx11" })
                 .Define("mpST12", mpStrawHits, { "Strawx12" })
                 .Define("strawSelection", strawSelection, { "Strawx11" })
                 .Define("Weight", "AWeight")
                 .Define("veto01", "VETO[0]")
                 .Define("veto23", "VETO[1]")
                 .Define("veto45", "VETO[2]");

    // Take a snapshot of the processed dataframe
    dOut_pion.Snapshot("training_set", fnameOut_pion.c_str(), { "IsDimuon", "HCAL012", "HCAL0", "HCAL1", "HCAL2", "ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection", "Weight", "veto01", "veto23", "veto45"  }, opts);

    //
    // KAON EVENTS DATA
    //

    std::string fnameIn_kaon = (fIn == "") ? "/eos/user/e/ezaya/simulation_output/KAON/data/kaon-2023A.root" : fIn; // Dataset of Kaons
    std::string fnameOut_kaon = (fOut == "") ? "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_kaon_out.root" : fOut;

    // Define the RDataFrame with the input ROOT file
    ROOT::RDataFrame d0_kaon("ana_tree", fnameIn_kaon);

    // Make all the necessary filters
    //auto dFiltered = d0.Filter(("DIMU_NDimuonsTrue < 2")); // Save single dimuon event or no dimuon event
    //auto dMuonHit = dFiltered.Filter("DIMU_PDGID == 13")

    auto dMagnet_kaon = d0_kaon.Filter(filterMagnetIronBlock, { "Strawy11", "Strawx12", "Strawy12" });

    // Define the dataframe transformations
    auto dOut_kaon = dMagnet_kaon.Define("IsDimuon", "DIMU_NDimuonsTrue")
                 .Define("eH_PDS", HcalPedestral, {"HCALenergy"})
                 .Define("HCAL0", "Sum(Take(eH_PDS, 9))")  // take first 9 HC0 distribution
                 .Define("HCAL1", "Sum(Take(eH_PDS, 18)) - Sum(Take(eH_PDS, 9))")  // take first 9 HC1 distribution
                 .Define("HCAL2", "Sum(Take(eH_PDS, 27)) - Sum(Take(eH_PDS, 18))")  // take first 9 HC2 distribution
                 .Define("HCAL012", "Sum(Take(eH_PDS, 27))")
                 .Define("ECAL", "Sum(ECALenergy)")
                 .Define("eH0_11", "eH_PDS[4]")
                 .Define("eH1_11", "eH_PDS[13]")
                 .Define("eH2_11", "eH_PDS[22]")
                 .Define("mpST11", mpStrawHits, { "Strawx11" })
                 .Define("mpST12", mpStrawHits, { "Strawx12" })
                 .Define("strawSelection", strawSelection, { "Strawx11" })
                 .Define("Weight", "AWeight")
                 .Define("veto01", "VETO[0]")
                 .Define("veto23", "VETO[1]")
                 .Define("veto45", "VETO[2]");

    // Take a snapshot of the processed dataframe
    dOut_kaon.Snapshot("training_set", fnameOut_kaon.c_str(), { "IsDimuon", "HCAL012", "HCAL0", "HCAL1", "HCAL2", "ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection", "Weight", "veto01", "veto23", "veto45" }, opts);

    //
    // COMMON EVENTS: NO HITS IN HCAL
    //
    
    std::string fnameIn_all = (fIn == "") ? "/eos/user/e/ezaya/simulation_output/DIMUONS/data/25032024_standardRun2.root" : fIn; // Dataset of all events, not only dimuons
    std::string fnameOut_all = (fOut == "") ? "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_common_out.root" : fOut; // Create a seperate out-file, seperating the dimuon out-file.

    // Define the RDataFrame with the input ROOT file
    ROOT::RDataFrame d0_all("ana_tree", fnameIn_all);

    // Make all the necessary filters
    //auto dMagnet = dFiltered.Filter(filterMagnetIronBlock, { "Strawy11", "Strawx12", "Strawy12" });

    // Define the dataframe transformations
    auto dOut_all = d0_all.Define("IsDimuon", "DIMU_NDimuonsTrue")
                 .Define("eH_PDS", HcalPedestral, {"HCALenergy"})
                 .Define("HCAL0", "Sum(Take(eH_PDS, 9))")  // take first 9 HC0 distribution
                 .Define("HCAL1", "Sum(Take(eH_PDS, 18)) - Sum(Take(eH_PDS, 9))")  // take first 9 HC1 distribution
                 .Define("HCAL2", "Sum(Take(eH_PDS, 27)) - Sum(Take(eH_PDS, 18))")  // take first 9 HC2 distribution
                 .Define("HCAL012", "Sum(Take(eH_PDS, 27))")
                 .Define("ECAL", "Sum(ECALenergy)")
                 .Define("eH0_11", "eH_PDS[4]")
                 .Define("eH1_11", "eH_PDS[13]")
                 .Define("eH2_11", "eH_PDS[22]")
                 .Define("mpST11", mpStrawHits, { "Strawx11" })
                 .Define("mpST12", mpStrawHits, { "Strawx12" })
                 .Define("strawSelection", strawSelection, { "Strawx11" })
                 .Define("veto01", "VETO[0]")
                 .Define("veto23", "VETO[1]")
                 .Define("veto45", "VETO[2]");

    // Take a snapshot of the processed dataframe
    dOut_all.Snapshot("training_set", fnameOut_all.c_str(), { "IsDimuon", "HCAL012", "HCAL0", "HCAL1", "HCAL2", "ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection", "veto01", "veto23", "veto45" }, opts);
}
