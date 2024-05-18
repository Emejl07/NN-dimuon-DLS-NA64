import uproot
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score

class NeuralNetwork(nn.Module):
    def __init__(self, input_size, hidden_size1, hidden_size2, output_size):
        super(NeuralNetwork, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size1)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size1, hidden_size2)
        self.fc3 = nn.Linear(hidden_size2, output_size)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        out = self.relu(out)
        out = self.fc3(out)
        out = self.sigmoid(out)
        return out

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
    features = df[["ECAL", "eH0_11", "eH1_11", "eH2_11", "mpST11", "mpST12", "strawSelection"]].values
    labels = df["IsDimuon"].values
    return features, labels

def train_model(model, X_train, y_train, criterion, optimizer, num_epochs=10000, stopping_threshold=1e-7):
    X_train_tensor = torch.Tensor(X_train)
    y_train_tensor = torch.Tensor(y_train).unsqueeze(1)

    loss_curve = []

    prev_loss = float('inf')
    for epoch in range(num_epochs):
        optimizer.zero_grad()
        outputs = model(X_train_tensor)
        loss = criterion(outputs, y_train_tensor)
        loss.backward()
        optimizer.step()
        loss_curve.append(loss.item())
        print(f'Epoch [{epoch + 1}/{num_epochs}], Loss: {loss.item()}')

        # Check if loss change is below the stopping threshold
        if abs(prev_loss - loss.item()) < stopping_threshold:
            print(f'Loss change below threshold ({stopping_threshold}), stopping training.')
            break

        prev_loss = loss.item()

    return model, loss_curve

def evaluate_model(model, X_test, y_test):
    X_test_tensor = torch.Tensor(X_test)
    outputs = model(X_test_tensor)
    y_pred = torch.round(outputs).detach().numpy().flatten()  # Round the output probabilities to obtain binary predictions
    accuracy = np.mean(y_pred == y_test)
    print('Test accuracy:', accuracy)
    return accuracy

def plot_loss_curve(loss_curve):
    plt.plot(loss_curve, 'b', label='Training loss')
    plt.title('Training loss curve')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.yscale('log') 
    plt.legend()
    plt.tight_layout()
    plt.savefig("loss_curve.pdf")
    plt.show()

def plot_roc_curve(y_test, y_scores, accuracy):
    fpr, tpr, thresholds = roc_curve(y_test, y_scores)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (AUC = %0.2f)' % roc_auc)
    plt.fill_between(fpr, tpr, color='darkorange', alpha=0.2, hatch='/')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC)\nAccuracy = %0.2f' % accuracy)
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig("ROC_curve.pdf")
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
    df = pd.concat([df_dimuon, df_pion, df_kaon])

    # Preprocess data
    features, labels = preprocess_data(df)
    
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)
    
    # Define the neural network model
    input_size = X_train.shape[1]
    hidden_size1 = 64
    hidden_size2 = 32
    output_size = 1
    model = NeuralNetwork(input_size, hidden_size1, hidden_size2, output_size)

    # Define loss function and optimizer
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # Train the model
    trained_model, loss_curve = train_model(model, X_train, y_train, criterion, optimizer)

    # Evaluate the model
    test_acc = evaluate_model(trained_model, X_test, y_test)
    
    # Calculate ROC curve
    X_test_tensor = torch.Tensor(X_test)
    outputs = trained_model(X_test_tensor)
    y_scores = outputs.detach().numpy()
    plot_roc_curve(y_test, y_scores, test_acc)
    plot_loss_curve(loss_curve)

    # Save model
    torch.save(trained_model, "models/dimuon_selection_model.pt")

if __name__ == "__main__":
    main()

