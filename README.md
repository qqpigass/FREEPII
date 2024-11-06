## FREEPII <br />

### Basic data
The human CF-MS datasets (PXD002892, PXD014820, and PXD015406) were downloaded from Zenodo (doi: 10.5281/zenodo.4106578) [1]. Among the various files corresponding to different protein quantification strategies provided by the authors, we selected those containing iBAQ intensity of chromatograms for further analysis. <br /> 

The human protein complex dataset was downloaded from the CORUM database, the human protein sequences and subcellular locations of human proteins were retrieved from UniProt database. <br />

The semantics and relationships between GO terms were retrieved from the GO Consortium released on November 4, 2022. <br /><br />

### Preprocess data and generate input data for the model
Please follow the steps below to generate the data required for model training and performance evaluation. <br />
1. The protein complexes set is used to generate labels for PPIs in our task. In this preprocessing step, we will generate a new set of protein complexes containing only members from the CF-MS data and filter out protein complexes consisting of fewer than three genes. Please follow the script: code/preprocess/preprocess_complexes.R,  to generate filtered protein complexes set for further analysis.
2. To converted the protein sequences into frequency matrix chaos game representation (FCGR) [2], please follow the script: code/preprocess/Generate_seq_FCGR_16x.R. The scaling factor is set to 0.863271 to prevent the overlap of attractors.
3. Preprocessing of CF-MS data consists of removing samples with missing or all-zero values ​​and then normalizing the values ​​to a range between 0-1. Please follow the script: code/preprocess/Preprocess_EPF.R, to generate needed data.
4. To fairly compare with other models, we generated fixed splits for FREEPII and other models. Please follow the scripts: code/preprocess/Generate_split.R, and code/preprocess/Generate_cv_split_csv.R, to generate needed split data.
5. Finally, use the codes in code/preprocess/Prepare_name-idx-dict_and_cv_input.ipynb, and code/preprocess/Prepare_seq_FCGR_16x_input.ipynb to convert the pre-processed data into input data suitable for FREEPII.
<br />

### Run FREEPII
The source code of FREEPII is in the code/FREEPII/FREEPII.ipynb. This code will use SEC2-heavy as an example. <br />
This model is trained under a 5-fold cross-validation scheme. For each fold, the ratio between train:test = 70:30 and the ratio between positive labels:negative labels = 1:1. <br />
The training acc, testing acc, training loss and testing loss of each fold will be recorded and used to confirm whether the model is overfitted. <br />
The best model in each fold is determined by the lowest testing loss in that fold and will be evaluated by true positive and true negative scores. <br />
The final best model is determined by the lowest testing loss across all folds and will be used to generate predictions for all PPIs in the experimental data for cluster analysis. <br />
The output of FREEPII will include the PPI prediction results of training set, heldout set, del-fold set as well as the prediction results of all PPIs, the protein feature representation learned by the model, and the predicted protein complexes.
<br />
<br />

### Evaluate model performance
Note: The go-basic in GO-data needs to be decompressed after downloading. <br />
1. PPI classification performance: The average predictions of FREEPII across all folds will be used to evaluate model overall performance on PPI classification. Evaluation indicators include sensitivity, specificity, Matthews correlation coefficient (MCC), area under the receiver-operator characteristic curve (AUC of ROC), etc. The heldout set is constructed under imbalanced conditions, while the del-fold set consists entirely of negative labeled data. Evaluating the model's performances on these sets can reflect how well the model handles imbalanced data under balanced training. Follow the scripts: code/performance/Performance_PPI.R to get the performance scores. 
2. Evaluation of the quality of predicted protein complexes (with gold standard): To evaluate the structural compositional similarity between predicted protein complexes and a reference protein complex dataset (gold standard), we use the composite score [3]. The composite score is the sum of three components: Overlap, Accuracy, and Maximum Matching Ratio (MMR) [4]. Follow the scripts: code/performance/Composite_score_Complex.R to calculate the composite scores of the predicted protein complexes.
4. Evaluation of the quality of predicted protein complexes (without gold-standard): To assess the similarity of protein localization within a protein complex, we measured the colocalization score (defined from [5]) of each predicted protein complex. To assess the functional similarity between proteins within the same protein complex, we used the GOGO method described in [6] to calculate the GOGO score for each predicted protein complex. Follow the code: code/performance/CoLocalization_score_Complex.R, and code/performance/GO_score_Complex.ipynb, to compute these scores.
<br />
<br />

### References
1. Skinnider,M.A. and Foster,L.J. (2021) Meta-analysis defines principles for the design and analysis of co-fractionation mass spectrometry experiments. Nat. Methods, 18, 806–815.
2. Almeida,J.S., Carriç,J.A., Ant´,A., Maretzek,A., Noble,P.A. and Fletcher,M. (2001) Analysis of genomic sequences by Chaos Game Representation. Bioinformatics, 17, 429–437.
3. Hu LZM, Goebels F, Tan JH, et al. (2019) EPIC: software toolkit for elution profile-based inference of protein complexes. Nat Method, 16, 737–742.
4. Nepusz T, Yu H, Paccanaro A. (2012) Detecting overlapping protein complexes in protein-protein interaction networks. Nat Methods, 9, 471–472
5. Krumsiek J, Zimmer R, Friedel CC. (2009) Bootstrapping the interactome: unsupervised identification of protein complexes in yeast. J. Comput. Biol, 16, 971–987.
6. Zhao C, Wang Z. (2018) GOGO: An improved algorithm to measure the semantic similarity between gene ontology terms. Sci Rep, 8, 15107.
<br />
