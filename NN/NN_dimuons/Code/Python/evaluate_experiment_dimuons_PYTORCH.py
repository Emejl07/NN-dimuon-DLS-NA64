import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import seaborn as sns
import torch
from NN_dimuons_pyTorch import NeuralNetwork
from DNN_dimuons_pyTorch import DeepNeuralNetwork

def load_data(file_path, treeName):
    # Load data from ROOT file using uproot
    file = uproot.open(file_path)
    tree = file[treeName]

    # Convert ROOT tree to pandas DataFrame
    df = tree.arrays(library="pd")
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)  
    #print(df["eH0_11"])
    return df

def preprocess_data(df):
    # Convert DataFrame to numpy arrays
    features = df[["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection", "veto01", "veto23", "veto45"]].values
    indices = df.index.values
    return features, indices

def load_model(filename):
    # Load the model state dictionary
    model = torch.load(filename)
    return model

def main():
    # Load dimuon data
    file_dimuon = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/experimentSet_out.root"
    df_dimuon = load_data(file_dimuon, "training_set")
    
    # Load trained model
    model = load_model("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Models/dimuon_selection_model_test.pt")

    # Preprocess data
    features, indices = preprocess_data(df_dimuon)

    # Convert features to PyTorch tensor
    X_tensor = torch.Tensor(features)

    # Perform inference
    with torch.no_grad():
        model.eval()  # Set the model to evaluation mode
        predictions = model(X_tensor).numpy()

    # Filter predicted dimuon events (label = 1)
    treshold = 0.5
    dimuon_indices = predictions[:, 0] > treshold

    dimuon_features = features[dimuon_indices]
    selected_indices = indices[dimuon_indices]
    print(r"Number of $\mu\mu$: %2i", len(dimuon_features))

    pd.DataFrame(selected_indices, columns=["index"]).to_csv("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Input/Selected_Indices.csv", index=False)

    selected_data = pd.DataFrame(dimuon_features, columns=["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection", "veto01", "veto23", "veto45"])
    selected_data["index"] = selected_indices
    selected_data.to_csv("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Input/selected_data.csv", index=False)
    np.save('/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Input/NN_Selected_experiment.npy', dimuon_features)

    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot histograms of features
    sns.histplot(dimuon_features[:, 1], bins=60, color='blue', edgecolor='black', ax=axes[0])
    axes[0].set_ylabel('Frequency')
    axes[0].set_title(r'Energy Center HCAL$_0$')
    axes[0].set_xlim(0, 20)

    sns.histplot(dimuon_features[:, 2], bins=60, color='red', edgecolor='black', ax=axes[1])
    axes[1].set_xlabel('Energy [GeV]')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Energy Center HCAL$_1$')
    axes[1].set_xlim(0, 20)

    sns.histplot(dimuon_features[:, 3], bins=60, color='green', edgecolor='black', ax=axes[2])
    axes[2].set_ylabel('Frequency')
    axes[2].set_title('Energy Center HCAL$_2$')
    axes[2].set_xlim(0, 20)

    # Adjust layout
    straw_indices = dimuon_features[:, 6] == 20
    print(len(dimuon_features[straw_indices]))

    plt.tight_layout()
    plt.savefig("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Plots/NN_results.pdf")
    plt.show()

if __name__ == "__main__":
    main()
