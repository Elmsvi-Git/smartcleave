________________________________________

CleaveSmart: Deciphering the Puzzle of Rational 10–23 DNAzyme Selection through Interpretable AI Insights 
________________________________________

CleaveSmart is a transparent, data-driven machine learning framework to replace trial-and-error with a transparent, data-driven DNAzyme 10-23 selection process. It provides a robust pipeline to predict whether the sustained cleavage efficiency (measured after 60 minutes) of DNAzyme 10–23 on RNA substrates exceeds the 40% threshold.



## ✔ Project Overview
The CleaveSmart framework addresses the practical challenges of DNAzyme selection by:
- Integrating an expanded set of thermodynamic and structural descriptors.
- Utilizing an optimized Random Forest architecture (the final CleaveSmart model).
- Providing Interpretable AI insights through permutation feature importance analysis to reveal the mechanistic drivers of cleavage.

---

## 🛠 Feature Engineering
CleaveSmart utilizes five high-impact molecular features:
1. RNA–DNA Duplex Binding Energy
2. DNAzyme Internal Energy
3. DNAzyme Homodimer Energy
4. RNA-Cost Energy
5. Total Energy (Delta G_binding - Delta G_internal + Delta G_RNA-cost)

---

## 🤖 Machine Learning Workflow
The following workflow was applied to both Sustained (60-min) and Rapid (10-min) activity datasets:

1. Feature dataset loading and preprocessing (StandardScaler).
2. Training and testing data partitioning.
3. Model training using four classifiers (Logistic Regression, Random Forest, SVC, and Gradient Boosting).
4. Hyperparameter optimization via GridSearchCV with cross-validation.
5. Performance evaluation using F1-score, Balanced Accuracy, Precision, Recall, and AUC.
6. Selection of the best model based on F1-score (CleaveSmart).
7. Permutation feature importance analysis.
8. External validation using an independent dataset (achieving Accuracy = 0.91, F1-score = 0.95 for the sustained model).

---

## 📊 Framework Outputs
The CleaveSmart pipeline is designed to provide actionable insights for experimental design, generating the following outputs:

- **Top Candidate Identification**: Automatically identifies and highlights the most promising DNAzyme sequence with the highest predicted efficiency for the target RNA.
- **Comprehensive Enzyme Library**: Generates a complete list of all potential DNAzymes that can target the provided substrate.
- **Detailed Feature Profiling**: Exports a full dataset of all calculated thermodynamic and structural features for every candidate in the library.
- **Probability-Based Scoring**: Provides the specific probability (ranging from 0 to 1) for each candidate to achieve Class 1 status (>40% cleavage efficiency), allowing for fine-grained ranking.
---

## 📂 Supplementary Folder Content
1. **Feature_Calculation_Script.py**: Functions for calculating thermodynamic features.
2. **Sustained_Activity_Workflow.py**: Script for CleaveSmart training, optimization, and testing.
3. **Rapid_Activity_Workflow.py**: Parallel analysis for rapid (10-min) cleavage activity.
4. **Datasets**: Includes Sustained, Rapid, and External Validation datasets.

---

## 💻 Requirements
- Python 3.11
- NumPy 2.2.1
- Pandas 2.2.2
- scikit-learn 1.5.1
- ViennaRNA 2.7.0
- Matplotlib 3.9.2 
- Seaborn 0.13.2
- Streamlit 1.54.0
- Biopython 

---

## 🚀 Access
The CleaveSmart platform is publicly available at:  
🔗 [https://cleavesmart.app](https://cleavesmart.app)

