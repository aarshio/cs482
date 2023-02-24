from joblib import dump, load
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn import svm
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# Load real and fake peptide datasets
with open('real_peptides.txt') as f:
    real_data = f.readlines()
real_data = [x.strip() for x in real_data]

with open('fake_peptides.txt') as f:
    fake_data = f.readlines()
fake_data = [x.strip() for x in fake_data]

# Create labels for real and fake peptides
y = np.concatenate((np.ones(len(real_data)), np.zeros(len(fake_data))))


# Combine real and fake peptides into a single array
X = np.concatenate((real_data, fake_data))


# Transform peptides into numerical feature vectors
# You will need to decide on an appropriate feature representation for your peptides
# One example is to use the amino acid composition of the peptides
# You can use a library such as BioPython to calculate the amino acid composition
# Alternatively, you can use other feature representations such as physicochemical properties or k-mer counts
# Here is an example code snippet to calculate amino acid composition
X_features = []
for peptide in X:
    protein = ProteinAnalysis(peptide)
    aa_comp = protein.get_amino_acids_percent()
    X_features.append(list(aa_comp.values()))


# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_features, y, test_size=0.2, random_state=42)
# Train an SVM on the training data
clf = svm.SVC(kernel='linear', C=1, probability=True)
clf.fit(X_train, y_train)

# Test the SVM on the testing data
y_pred = clf.predict(X_test)

# Evaluate performance using metrics such as accuracy, precision, recall, F1-score, and AUC-ROC
print('Accuracy:', accuracy_score(y_test, y_pred))
print('Precision:', precision_score(y_test, y_pred))
print('Recall:', recall_score(y_test, y_pred))
print('F1-score:', f1_score(y_test, y_pred))
print('AUC-ROC:', roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1]))

# Save the trained model using joblib
dump(clf, 'svm_model.joblib')
