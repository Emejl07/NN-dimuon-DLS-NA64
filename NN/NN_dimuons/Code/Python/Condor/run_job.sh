#!/bin/bash

# Set up the environment
source /afs/cern.ch/user/e/ezaya/setup_el9.sh

# Stage in the script from EOS
cp /eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Code/Python/DNN_dimuons_pyTorch.py .

# Run the Python script
python3 DNN_dimuons_pyTorch.py

# Stage out the output model to EOS
cp dimuon_selection_model_test.pt /eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Models/


