import uproot
import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score
import matplotlib.pyplot as plt
import pandas as pd
import joblib
import torch

def load_data(file_path, num_samples, treeName):
    # Load data from ROOT file using uproot
    file = uproot.open(file_path)
    tree = file[treeName]
    
    # Convert ROOT tree to pandas DataFrame
    df = tree.arrays(library="pd").head(num_samples)
    print(df.size)
    return df

def preprocess_data(df):
    # Convert DataFrame to numpy arrays
    # features = ["eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]
    features = df[["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]].values
    labels = df["IsDimuon"].values
    return features, labels

def train_model(X_train, y_train):
    # Define the neural network model (MLPClassifier)
    model = MLPClassifier(hidden_layer_sizes=(64, 32), activation='relu', solver='adam', max_iter=300, verbose=True)
    
    # Train the model
    model.fit(X_train, y_train)
    return model

def evaluate_model(model, X_test, y_test):
    # Evaluate the model on the testing set
    test_acc = model.score(X_test, y_test)
    print('Test accuracy:', test_acc)
    return test_acc

def plot_loss_curve(loss_curve):
    # Plot the loss curve
    plt.plot(loss_curve, 'b', label='Training loss')
    plt.title('Training loss curve')
    plt.xlabel('Iterations')
    plt.ylabel('Loss')
    plt.yscale('log') 
    plt.legend()
    plt.savefig("loss_curve.pdf")
    plt.show()

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

def plot_feature_label_relationships(df, features, label):
    # Create a 3x3 grid of scatter plots
    fig, axes = plt.subplots(nrows=7, ncols=7, figsize=(15, 15))
    
    # Flatten the axes array for easier iteration
    axes = axes.flatten()
    
    # Define colors for the labels
    colors = {0: 'blue', 1: 'red'}
    
    for i, feature1 in enumerate(features):
        for j, feature2 in enumerate(features):
            if i != j:  # Exclude diagonal plots
                ax = axes[i * len(features) + j]
                for lbl, color in colors.items():
                    ax.scatter(df[df[label] == lbl][feature1], df[df[label] == lbl][feature2], c=color, label=f'IsDimuon={lbl}', alpha=0.5)
                ax.set_xlabel(feature1)
                ax.set_ylabel(feature2)
                ax.legend()
                ax.set_title(f"{feature1} vs {feature2}")
    
    # Hide empty subplot(s)
    for k in range(len(features) ** 2):
        if k >= len(features) * len(features) - (len(features)):
            axes[k].axis('off')
    
    fig.tight_layout()
    plt.show()

def plot_feature_relationship(df, feature1, feature2, label, filename):
    # Plot scatter plot of selected features
    plt.figure(figsize=(8, 6))
    
    # Define colors for the labels
    colors = {0: 'blue', 1: 'red'}
    
    for lbl, color in colors.items():
        plt.scatter(df[df[label] == lbl][feature1], df[df[label] == lbl][feature2], c=color, label=f'IsDimuon={lbl}', alpha=0.5, s=20)  # Adjust the size (s) as needed
    
    plt.xlabel(feature1)
    plt.ylabel(feature2)
    plt.legend()
    plt.title(f"{feature1} vs {feature2}")
    plt.savefig(filename)  # Save the plot as PDF with specified filename
    plt.show()

def main():
    # Load data
    file_dimuon = "TrainingSet_dimuon_out.root"
    file_all = "TrainingSet_common_out.root"
    file_pion = "TrainingSet_pion_out.root"
    file_kaon = "TrainingSet_kaon_out.root"
    num_samples = 24000
    df_dimuon = load_data(file_dimuon, num_samples, "training_set")
    df_all = load_data(file_all, num_samples, "training_set")
    df_pion = load_data(file_pion, num_samples, "training_set")
    df_kaon = load_data(file_kaon, num_samples, "training_set")
    df = pd.concat([df_dimuon, df_pion, df_kaon, df_all], ignore_index=True)

    print(df)
    
    # Plot feature-label relationships
    # features = ["eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]
    features = ["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]
    label = "IsDimuon"
    #plot_feature_label_relationships(df, features, label)
    plot_feature_relationship(df, "eH1_11", "eH1_11", "IsDimuon", "plot1.pdf")
    plot_feature_relationship(df, "ECAL", "eH1_11", "IsDimuon", "plot1.pdf")
    plot_feature_relationship(df, "eH2_11", "mpST12", "IsDimuon", "plot1.pdf")
 
    # Preprocess data
    features, labels = preprocess_data(df)
    
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)
    
    # Train the model
    model = train_model(X_train, y_train)
    loss_curve = model.loss_curve_
    plot_loss_curve(loss_curve)
    
    # Evaluate the model
    test_acc = evaluate_model(model, X_test, y_test)
    
    # Calculate ROC curve
    y_scores = model.predict_proba(X_test)[:, 1]
    plot_roc_curve(y_test, y_scores)

    # Save model
    joblib.dump(model, "models/dimuon_selection_model.pkl")

if __name__ == "__main__":
    main()
