import numpy as np
from sklearn.metrics import roc_curve, auc, accuracy_score
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import joblib
import seaborn as sns

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
    # features = ["eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]
    features = df[["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]].values
    #labels = df["IsDimuon"].values
    return features

def plot_roc_curve(y_test, y_scores):
    # Calculate ROC curve
    fpr, tpr, thresholds = roc_curve(y_test, y_scores)
    roc_auc = auc(fpr, tpr)

    # Calculate accuracy
    y_pred = (y_scores > 0.5).astype(int)  # Convert scores to binary predictions
    accuracy = accuracy_score(y_test, y_pred)

    # Print accuracy
    print("Accuracy:", roc_auc)
    
    # Plot ROC curve
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (AUC = %0.2f)' % roc_auc)
    
    # Highlight area under the curve with striped lines
    plt.fill_between(fpr, tpr, color='darkorange', alpha=0.2, hatch='/')
    
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC)\nAccuracy = %0.2f' % accuracy)
    plt.legend(loc="lower right")
    plt.savefig("ROC_curve.pdf")
    plt.show()

def load_model(filename):
    # Load the trained model
    return joblib.load(filename)

def evaluate_model(model, X_test, y_test):
    # Evaluate the model on the testing set
    test_acc = model.score(X_test, y_test)
    print('Test accuracy:', test_acc)
    return test_acc

def main():
    # Load dimuon data
    file_dimuon = "experimentSet_out.root"
    file_sim = "TrainingSet_dimuon_out.root"
    df_dimuon = load_data(file_dimuon, "training_set")
    df_sim = load_data(file_sim, "training_set")
    
    # Load trained model
    model = load_model("models/dimuon_selection_model.pkl")

    # Preprocess data
    features = preprocess_data(df_dimuon)

    predictions = model.predict(features)

    print(predictions)

    # Filter predicted dimuon events (label = 1)
    dimuon_indices = predictions == 1
    dimuon_features = features[dimuon_indices]
    print(dimuon_features[:, 3].size)

    # Create subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Plot histograms of features
    sns.histplot(dimuon_features[:, 1], bins=50, color='blue', edgecolor='black', ax=axes[0])
    axes[0].set_xlabel('Feature 1')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Histogram of Feature 1')
    axes[0].set_xlim(0, 20)

    sns.histplot(dimuon_features[:, 2], bins=50, color='red', edgecolor='black', ax=axes[1])
    axes[1].set_xlabel('Feature 2')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Histogram of Feature 2')
    axes[1].set_xlim(0, 20)

    sns.histplot(dimuon_features[:, 3], bins=50, color='green', edgecolor='black', ax=axes[2])
    axes[2].set_xlabel('Feature 3')
    axes[2].set_ylabel('Frequency')
    axes[2].set_title('Histogram of Feature 3')
    axes[2].set_xlim(0, 20)

    # Adjust layout
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
