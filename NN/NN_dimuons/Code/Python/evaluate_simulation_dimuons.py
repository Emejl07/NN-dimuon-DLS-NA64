import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot 
import seaborn as sns
import torch
from sklearn.model_selection import train_test_split
from NN_dimuons_pyTorch import NeuralNetwork
from DNN_dimuons_pyTorch import DeepNeuralNetwork


def load_data(file_path, num_samples, treeName):
    file = uproot.open(file_path)
    tree = file[treeName]
    df = tree.arrays(library="pd").head(num_samples)

    df = tree.arrays(library="pd")
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    #print(df)

    return df

def preprocess_data(df):
    # Convert DataFrame to numpy arrays
    features = df[["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]].values
    labels = df["IsDimuon"].values  # Assuming the label column is "IsDimuon"
    indices = df.index.values
    return features, labels, indices

#def preprocess_data(df):
    # Convert DataFrame to numpy arrays
#    features = df[["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]].values
#    return features

def load_model(filename):
    # Load the model state dictionary
    model = torch.load(filename)
    return model

def calculate_tpr_fpr(y_true, y_scores, threshold):
    # Calculate True Positive Rate (TPR) and False Positive Rate (FPR)
    y_pred = y_scores > threshold
    tp = np.sum((y_pred == 1) & (y_true == 1))
    fp = np.sum((y_pred == 1) & (y_true == 0))
    fn = np.sum((y_pred == 0) & (y_true == 1))
    tn = np.sum((y_pred == 0) & (y_true == 0))

    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0
    
    return tpr, fpr

def plot_signal_background_histogram(y_train, y_scores_train, y_test, y_scores_test):
    fig, ax = plt.subplots(figsize=(9, 7))

    # Training histograms
    plt.hist(y_scores_train[y_train == 1], bins=50, alpha=0.5, label='Dimuons (Training)', color='blue', histtype='stepfilled', density=True)
    plt.hist(y_scores_train[y_train == 1], bins=50, color='blue', histtype='step', linewidth=1.5, density=True)  # Outline
    plt.hist(y_scores_train[y_train == 0], bins=50, alpha=0.5, label='Non-Dimuons (Training)', color='red', histtype='step', hatch='///', edgecolor='red', linewidth=1.5, density=True)
    plt.hist(y_scores_train[y_train == 0], bins=50, color='red', histtype='step', linewidth=1.5, density=True)  # Outline

    # Test histograms and error calculations
    counts_signal, bin_edges_signal = np.histogram(y_scores_test[y_test == 1], bins=50, density=True)
    bin_centers_signal = 0.5 * (bin_edges_signal[1:] + bin_edges_signal[:-1])  # Calculate bin centers
    errors_signal = np.sqrt(counts_signal) / counts_signal.sum()

    counts_bkg, bin_edges_bkg = np.histogram(y_scores_test[y_test == 0], bins=50, density=True)
    bin_centers_bkg = 0.5 * (bin_edges_bkg[1:] + bin_edges_bkg[:-1])  # Calculate bin centers
    errors_bkg = np.sqrt(counts_bkg) / counts_bkg.sum()

    # Plot test data with errors
    plt.errorbar(bin_centers_signal, counts_signal, yerr=errors_signal, fmt='o', color='blue', label='Dimuons (Test)', zorder=3, capsize=5)
    plt.errorbar(bin_centers_bkg, counts_bkg, yerr=errors_bkg, fmt='o', color='red', label='Non-Dimuons (Test)', zorder=3, capsize=5)

    plt.xlabel('Predicted Probability', fontsize=20)
    plt.ylabel('Normalized Frequency', fontsize=20)
    leg = ax.legend(fontsize=18, loc="upper center")
    ax.tick_params(axis='y', labelsize=18)
    ax.tick_params(axis='x', labelsize=18)
    ax.set_yscale('log')
    leg.get_frame().set_edgecolor('none')
    leg.get_frame().set_facecolor('none')
    plt.tight_layout()
    plt.savefig("signal_background_hist_simulationSet.pdf")
    plt.show()

# Assuming you have data to use, call the function as needed:
# plot_signal_background_histogram(y_train, y_scores_train, y_test, y_scores_test)

def main():
    # Load dimuon data
    #file_dimuon = "SimulationSet_dimuon_out.root"
   # file_dimuon = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/SimulationSet_dimuon_out.root"
   # df_dimuon = load_data(file_dimuon, "training_set")

    file_dimuon = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_dimuon_out.root"
    file_all = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_common_out.root"
    file_pion = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_pion_out.root"
    file_kaon = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_kaon_out.root"
    num_samples = 400000
    df_dimuon = load_data(file_dimuon, num_samples, "training_set")
    df_all = load_data(file_all, num_samples, "training_set")
    df_pion = load_data(file_pion, num_samples, "training_set")
    df_kaon = load_data(file_kaon, num_samples, "training_set")
    df = pd.concat([df_dimuon, df_pion, df_kaon])

    # Load trained model
    model = load_model("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Models/dimuon_selection_model.pt")
    
    # Preprocess data
    features, labels, _ = preprocess_data(df)

    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)

    X_train_tensor = torch.Tensor(X_train)
    X_test_tensor = torch.Tensor(X_test)

    # Convert features to PyTorch tensor
   # X_tensor = torch.Tensor(features)

    # Perform inference
    with torch.no_grad():
        model.eval()  # Set the model to evaluation mode
        y_scores_train = model(X_train_tensor).numpy()[:, 0]
        y_scores_test = model(X_test_tensor).numpy()[:, 0]
        #predictions = model(X_tensor).numpy()

    # Filter predicted dimuon events (label = 1)
    threshold = 0.5
    dimuon_indices_train = y_scores_train > threshold
    dimuon_indices_test = y_scores_test > threshold
    plot_signal_background_histogram(y_train, y_scores_train, y_test, y_scores_test)

   # tpr, fpr = calculate_tpr_fpr(labels, predictions[:, 0], threshold)
    #print(f"TPR: {tpr}, FPR: {fpr}") 

    #dimuon_features = features[dimuon_indices]
    # Save dimuon_features to a .npy file
    #np.save("dimuon_features.npy", dimuon_features)
    #print(r"Number of $\mu\mu$: %2i", len(dimuon_features))

'''
    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot histograms of features
    sns.histplot(dimuon_features[:, 1], bins='auto', color='blue', edgecolor='black', ax=axes[0])
    #axes[0].set_xlabel('Feature 1')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title(r'Energy Center HCAL$_0$')
    axes[0].set_xlim(0, 20)
        
    sns.histplot(dimuon_features[:, 2], bins='auto', color='red', edgecolor='black', ax=axes[1])
    axes[1].set_xlabel('Energy [GeV]')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Energy Center HCAL$_1$')
    axes[1].set_xlim(0, 20)
    
    sns.histplot(dimuon_features[:, 3], bins='auto', color='green', edgecolor='black', ax=axes[2])
    #axes[2].set_xlabel('Feature 3')
    axes[2].set_ylabel('Frequency')
    axes[2].set_title('Energy Center HCAL$_2$')
    axes[2].set_xlim(0, 20)

    # Adjust layout
    plt.tight_layout()
    plt.savefig("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Plots/NN_results_simulation.pdf")
    plt.show()
'''
#    plot_signal_background_histogram(labels, predictions[:, 0])
   # plot_signal_background_histogram(y_train, y_scores_train, y_test, y_scores_test)

if __name__ == "__main__":
    main()
