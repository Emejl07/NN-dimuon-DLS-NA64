#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"

void create_ResultSet(std::string fIn = "", std::string fOut = "") {
    std::string fnameIn = (fIn == "") ? "/afs/cern.ch/work/e/ezaya/public/data_partblinded_pass5_small/job008747.root" : fIn;
    std::string fnameOut = (fOut == "") ? "NN_Selected_experimentSet.root" : fOut;

    // Load selected indices from CSV file
    std::ifstream file("Selected_Indices.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Vector to store the selected indices
    std::vector<int> selected_indices;
    std::string line;
    // Skip the header
    std::getline(file, line);

    // Read the indices
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int index;
        ss >> index;
        selected_indices.push_back(index);
    }
    file.close();

    // Open the original ROOT file
    TFile *inputFile = TFile::Open(fnameIn.c_str());
    if (!inputFile) {
        std::cerr << "Error opening ROOT file!" << std::endl;
        return;
    }

    // Get the tree from the ROOT file
    TTree *tree = (TTree*)inputFile->Get("tout");  // Replace with actual tree name

    // Create a new ROOT file to store the selected events
    TFile *outputFile = new TFile(fnameOut.c_str(), "RECREATE");
    TTree *newTree = tree->CloneTree(0);

    // Loop over the selected indices and fill the new tree
    for (int index : selected_indices) {
        tree->GetEntry(index);
        newTree->Fill();
    }

    // Write the new tree to the output file
    newTree->Write();
    outputFile->Close();
    inputFile->Close();
}
