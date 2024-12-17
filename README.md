## FREEPII <br />
> FREEPII (Feature Representation Enhancement End-to-end Protein Interaction Inference) is a end-to-end learning method encompassing autonomous feature extraction and feature representation enhancement for PPIs and protein complexes inference. <br />
> FREEPII establishes feature maps that are consistent with model training, and emphasizes learning the feature representation of proteins rather than the feature representation of PPI, reducing the computational complexity from N^2 to N (N: number of proteins). <br />
> FREEPII uses protein sequences to expand the information available for calculating data similarity between proteins. Furthermore, it introduces network-level information into the final feature representation to rescale the strength of interactions present in CF-MS data.

### Basic data
The human CF-MS datasets (PXD002892, PXD014820, and PXD015406) were downloaded from Zenodo (doi: 10.5281/zenodo.4106578) [1]. Among the various files corresponding to different protein quantification strategies provided by the authors, we selected those containing iBAQ intensity of chromatograms for further analysis. <br /> 

The human protein complex dataset was downloaded from the CORUM database, the human protein sequences and subcellular locations of human proteins were retrieved from UniProt database. <br />

The semantics and relationships between GO terms were retrieved from the GO Consortium released on November 4, 2022. <br /><br />

### Getting started
The go-basic in GO-data needs to be decompressed after downloading. Execute the following command
```
sudo apt update
sudo apt install p7zip-full
7z x /GO-data/go-basic.7z
```

To start running FREEPII, we create an environment for FREEPII and acticate it by executing the following command
```
conda env create –f environment.yml -n FREEPII
conda activate FREEPII
```
<br />

### Preprocess data and generate input data for the model
The protein complexes set is used to generate labels for PPIs. <br />
In this preprocessing step, we first generate a new set of protein complexes containing only members from the CF-MS data and filter out protein complexes consisting of fewer than three genes.
<br />
Next, we preprocessed the CF-MS data, including removing samples with missing or all-zero values, and then normalizing the values ​​to a range between 0-1.<br />
After processing the data, we generate the input required by FREEPII.

Execute the following command to run the preprocessing process
```
python -W ignore ./Code/preprocess_step.py
```
<br />

The default setting for the model is PXD002892, SEC2-heavy. <br />
If you want to analyze other data, you can specify the experiment name (-e_name) and experimental conditions (-e_cond)
```
python -W ignore ./Code/preprocess_step.py -e_name PXD014820 -e_cond Ctrl
```
<br />

### Run FREEPII
FREEPII is trained under a 5-fold cross-validation scheme. For each fold, the ratio between train:test = 70:30 and the ratio between positive labels:negative labels = 1:1. <br />
The training acc, testing acc, training loss and testing loss of each fold will be recorded and used to confirm whether the model is overfitted. <br />
The best model in each fold is determined by the lowest testing loss in that fold and will be evaluated by true positive and true negative scores. <br />
The final best model is determined by the lowest testing loss across all folds and will be used to generate predictions for all PPIs in the experimental data for cluster analysis. <br />
The output of FREEPII will include the PPI prediction results of training set, heldout set, del-fold set as well as the prediction results of all PPIs, the protein feature representation learned by the model, and the predicted protein complexes.
<br />

Execute the following command to train FREEPII and generate predicted PPIs and clusters
```
python -W ignore ./Code/run_FREEPII.py
```
<br />
<br />

### Evaluate model performance
Note: The go-basic in GO-data needs to be decompressed after downloading. The current GO association is established based on the go.obo in the folder. If you download a new go.obo, please re-establish the GO association, otherwise there may be an error when running the GOGO operation. <br />
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
