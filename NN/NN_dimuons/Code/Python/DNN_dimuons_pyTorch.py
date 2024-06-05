import uproot
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

class DeepNeuralNetwork(nn.Module):
    def __init__(self, input_size):
        super(DeepNeuralNetwork, self).__init__()
        self.fc1 = nn.Linear(input_size, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, 32)
        self.fc4 = nn.Linear(32, 16)
        self.fc5 = nn.Linear(16, 1)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(p=0.5)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.dropout(out)
        out = self.fc2(out)
        out = self.relu(out)
        out = self.dropout(out)
        out = self.fc3(out)
        out = self.relu(out)
        out = self.dropout(out)
        out = self.fc4(out)
        out = self.relu(out)
        out = self.dropout(out)
        out = self.fc5(out)
        out = self.sigmoid(out)
        return out

def load_data(file_path, num_samples, treeName):
    file = uproot.open(file_path)
    tree = file[treeName]
    df = tree.arrays(library="pd").head(num_samples)
    return df

def preprocess_data(df):
    features = df[["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection", "veto01", "veto23", "veto45"]].values
    #features = df[["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]].values
    labels = df["IsDimuon"].values
    return features, labels

def train_model(model, X_train, y_train, X_test, y_test, criterion, optimizer, num_epochs=200, stopping_threshold=1e-7):
    X_train_tensor = torch.Tensor(X_train)
    y_train_tensor = torch.Tensor(y_train).unsqueeze(1)
    X_test_tensor = torch.Tensor(X_test)
    y_test_tensor = torch.Tensor(y_test).unsqueeze(1)

    loss_curve = []
    train_acc_curve = []
    test_acc_curve = []

    prev_loss = float('inf')
    for epoch in range(num_epochs):
        model.train()
        optimizer.zero_grad()
        outputs = model(X_train_tensor)
        loss = criterion(outputs, y_train_tensor)
        loss.backward()
        optimizer.step()
        loss_curve.append(loss.item())

        # Calculate training accuracy
        with torch.no_grad():
            train_preds = torch.round(model(X_train_tensor))
            train_acc = (train_preds.eq(y_train_tensor).sum().item()) / y_train_tensor.size(0)
            train_acc_curve.append(train_acc)
            
            # Calculate test accuracy
            test_preds = torch.round(model(X_test_tensor))
            test_acc = (test_preds.eq(y_test_tensor).sum().item()) / y_test_tensor.size(0)
            test_acc_curve.append(test_acc)

        print(f'Epoch [{epoch + 1}/{num_epochs}], Loss: {loss.item()}, Train Acc: {train_acc}, Test Acc: {test_acc}')

        if abs(prev_loss - loss.item()) < stopping_threshold:
            print(f'Loss change below threshold ({stopping_threshold}), stopping training.')
            break

        prev_loss = loss.item()

    return model, loss_curve, train_acc_curve, test_acc_curve

def evaluate_model(model, X_test, y_test):
    X_test_tensor = torch.Tensor(X_test)
    outputs = model(X_test_tensor)
    y_pred = torch.round(outputs).detach().numpy().flatten()
    accuracy = np.mean(y_pred == y_test)
    print('Test accuracy:', accuracy)
    return accuracy

def plot_loss_and_accuracy_curve(loss_curve, train_acc_curve, test_acc_curve):
    fig, ax1 = plt.subplots(figsize=(7, 7))

    ax1.plot(loss_curve, 'b', label='Loss', lw=3)
    ax1.set_xlabel('Epochs', fontsize=20)
    ax1.set_ylabel('Binary Cross-Entropy (BCE) Loss', fontsize=20)
    ax1.set_yscale('log')
    ax1.tick_params(axis='y', labelsize=18)
    ax1.tick_params(axis='x', labelsize=18)

    ax2 = ax1.twinx()
    ax2.plot(train_acc_curve, 'g', label='Training accuracy', lw=3)
    ax2.plot(test_acc_curve, 'r--', label='Test accuracy', lw=3)
    ax2.set_ylabel('Accuracy', fontsize=20)
    ax2.tick_params(axis='y', labelsize=18)
    ax2.tick_params(axis='x', labelsize=18)
    ax2.set_ylim([0, 1])

    fig.tight_layout()
    leg1 = ax1.legend(fontsize=18, loc="lower left")
    leg1.get_frame().set_edgecolor('none')
    leg1.get_frame().set_facecolor('none') 

    leg2 = ax2.legend(fontsize=18, loc="center right")
    leg2.get_frame().set_edgecolor('none')
    leg2.get_frame().set_facecolor('none')
    plt.savefig("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Plots/loss_and_accuracy_curve_test.pdf")
    plt.show()

def plot_roc_curve(y_test, y_scores, accuracy):
    fpr, tpr, thresholds = roc_curve(y_test, y_scores)
    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(7, 7))
    plt.plot(fpr, tpr, color='darkorange', lw=3, label='ROC curve (AUC = 0.985)')
    plt.fill_between(fpr, tpr, color='darkorange', alpha=0.2, hatch='/')
    plt.plot([0, 1], [0, 1], color='navy', lw=3, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=20)
    plt.ylabel('True Positive Rate', fontsize=20)
    plt.title('Receiver Operating Characteristic (ROC)', fontsize=22)
    ax.tick_params(axis='y', labelsize=18)
    ax.tick_params(axis='x', labelsize=18)
    #plt.legend(loc="lower right")
    leg = ax.legend(fontsize=18, loc="lower right")
    leg.get_frame().set_edgecolor('none')
    leg.get_frame().set_facecolor('none')
    plt.tight_layout()
    plt.savefig("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Plots/ROC_curve_test.pdf")
    plt.show()

def plot_signal_background_histogram(y_test, y_scores):
    fig, ax = plt.subplots(figsize=(7, 7))
    plt.hist(y_scores[y_test == 1], bins=150, alpha=0.5, label='Dimuons (Signal)', color='blue', histtype='stepfilled')
    plt.hist(y_scores[y_test == 0], bins=150, alpha=0.5, label='Non-Dimuons (Background)', color='red', histtype='stepfilled')
    plt.xlabel('Predicted Probability', fontsize=20)
    plt.ylabel('Frequency', fontsize=20)
    plt.title('Predicted Probabilities', fontsize=22)
    leg = ax.legend(fontsize=18, loc="upper right")
    ax.tick_params(axis='y', labelsize=18)
    ax.tick_params(axis='x', labelsize=18)
    ax.set_yscale('log')
    leg.get_frame().set_edgecolor('none')
    leg.get_frame().set_facecolor('none')
    plt.tight_layout()
    plt.savefig("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Plots/signal_background_histogram_test.pdf")
    plt.show()

def main():
    # Load data
    file_dimuon = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_dimuon_out.root"
    file_common = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_common_out.root"
    file_pion = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_pion_out.root"
    file_kaon = "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/TrainingSet_kaon_out.root"
    num_samples = 400 # Max value: 556279
    df_dimuon = load_data(file_dimuon, num_samples, "training_set")
    df_all = load_data(file_common, num_samples, "training_set")
    df_pion = load_data(file_pion, num_samples, "training_set")
    df_kaon = load_data(file_kaon, num_samples, "training_set")
    df = pd.concat([df_dimuon, df_pion, df_kaon])

    # Preprocess data
    features, labels = preprocess_data(df)
    
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)
    
    # Define the neural network model
    input_size = X_train.shape[1]
    model = DeepNeuralNetwork(input_size)

    # Define loss function and optimizer
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # Train the model
    trained_model, loss_curve, train_acc_curve, test_acc_curve = train_model(model, X_train, y_train, X_test, y_test, criterion, optimizer)

    # Evaluate the model
    test_acc = evaluate_model(trained_model, X_test, y_test)
    
    # Calculate ROC curve
    X_test_tensor = torch.Tensor(X_test)
    outputs = trained_model(X_test_tensor)
    y_scores = outputs.detach().numpy()
    plot_roc_curve(y_test, y_scores, test_acc)
    
    # Plot loss and accuracy curves
    plot_loss_and_accuracy_curve(loss_curve, train_acc_curve, test_acc_curve)
    
    # Plot signal vs background histogram
    plot_signal_background_histogram(y_test, y_scores)
    
    # Save model
    torch.save(trained_model, "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Models/dimuon_selection_model_test.pt")

if __name__ == "__main__":
    main()

