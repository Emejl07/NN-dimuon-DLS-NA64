import uproot
import numpy as np
import tensorflow as tf
from tensorflow import keras
from sklearn.model_selection import train_test_split
#from sklearn.metrics import f1_score, roc_curve, auc
import matplotlib.pyplot as plt

# Load data from ROOT file using uproot
file = uproot.open("TrainingSet_out.root")
tree = file["training_set"]

# Convert ROOT tree to pandas DataFrame
df = tree.arrays(library="pd")
print(df)

# Convert DataFrame to numpy arrays
features = df[["HCAL012", "HCAL0", "HCAL1", "HCAL2", "ECAL"]].values
labels = df["IsDimuon"].values

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)

# Define the neural network model
model = keras.Sequential([
    keras.layers.Dense(64, activation='relu', input_shape=(5,)),
    keras.layers.Dense(32, activation='relu'),
    keras.layers.Dense(1, activation='sigmoid')
])

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# Train the model
model.fit(X_train, y_train, epochs=10, batch_size=32, validation_data=(X_test, y_test))

# Evaluate the model on the testing set
test_loss, test_acc = model.evaluate(X_test, y_test)
print('Test accuracy:', test_acc)

# Calculate F1-score
y_pred = model.predict(X_test)
y_pred_binary = np.round(y_pred).astype(int)

'''
f1 = f1_score(y_test, y_pred_binary)
print('F1-score:', f1)

# Plot ROC curve
fpr, tpr, _ = roc_curve(y_test, y_pred)
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()
'''
