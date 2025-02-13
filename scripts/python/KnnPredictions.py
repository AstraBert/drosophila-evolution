from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np

# Load datasets
csv = pd.read_csv("reordered_param_swd.tsv", sep="\t")
df = pd.read_csv("median_values.csv")
csv = csv.drop(columns=["scenario"])

keys = list(csv.keys())  # Column names

preds = []
truths = []

for k in keys:
    knn = KNeighborsClassifier(n_neighbors=20)
    
    X = csv.drop(columns=[k])  # Features for training
    y = csv[k]  # Target variable
    
    knn.fit(X, y)  # Train KNN model
    
    X_test = df.drop(columns=[k])  # Test features
    y_test = df[k].values  # Convert to NumPy array for compatibility
    
    truths.append(y_test)  # Store true values
    y_pred = knn.predict(X_test)  # Predictions
    preds.append(y_pred)

print(preds)
print(truths)