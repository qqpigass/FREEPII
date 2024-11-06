## Source code for FREEPII <br />
>> Currently only human data are provided as examples <br />

### Basic data
The human CF-MS datasets (PXD002892, PXD014820, and PXD015406) were downloaded from Zenodo (doi: 10.5281/zenodo.4106578). <br />
Among the various files corresponding to different protein quantification strategies provided by the authors, we selected those containing iBAQ intensity of chromatograms for further analysis. <br /> 

The human protein complex dataset was downloaded from the CORUM database, and the human protein sequences were retrieved from UniProt database. <br />

We retrieved the semantics and relationships between GO terms from the GO Consortium released on November 4, 2022. <br /><br />

### Preprocess data and generate input data for the model
Please follow the steps below to generate the data required for model training and performance evaluation. <br />
1. The protein complexes set is used to generate labels for PPIs in our task. In this preprocessing step, we will generate a new set of protein complexes containing only members from the CF-MS data and filter out protein complexes consisting of fewer than three genes. Please follow the script <<< code/preprocess/preprocess_complexes.R >>> ,  to generate filtered protein complexes set for further analysis.
2. To converted the protein sequences into frequency matrix chaos game representation (FCGR), please follow the script <<< code/preprocess/Generate_seq_FCGR_16x.R >>> . The scaling factor is set to 0.863271 to prevent the overlap of attractors.
3. Preprocessing of CF-MS data consists of removing samples with missing or all-zero values ​​and then normalizing the values ​​to a range between 0-1. Please follow the script <<< code/preprocess/Preprocess_EPF.R >>> , to generate needed data.
4. To fairly compare with other models, we generated fixed splits for FREEPII and other models. Please follow the scripts <<< code/preprocess/Generate_split.R >>> , and <<< code/preprocess/Generate_cv_split_csv.R >>> , to generate needed split data.
5. Finally, use the codes in <<< code/preprocess/Prepare_name-idx-dict_and_cv_input.ipynb >>> , and <<< Prepare_seq_FCGR_16x_input.ipynb >>> to convert the pre-processed data into input data suitable for FREEPII.
<br />

### Run FREEPII
The source code of FREEPII is in the code/FREEPII/FREEPII.ipynb. This code will use SEC2-heavy as an example. <br />
The model is trained under the scheme of training: testing = 70:30, positive labels: negative labels = 1:1, and performs 5-fold cross-validation. <br />
The training set is further divided into training set and validation set in a ratio of 80:20. <br />
The training acc, validation acc, training loss and validation loss of each fold will be recorded and used to confirm whether the model is overfitted. <br />
The best model in each fold will be evaluated by the true positive and true negative scores. <br />
The best model among all folds will be used to generate predictions for all PPIs in the experimental data for clustering analysis. <br />
The output of FREEPII will include the PPI prediction results of training split, held-out split, del-fold split as well as the prediction results of all PPIs, the protein feature representation learned by the model, and the predicted protein complexes. 
<br />
<br />

### Evaluation model performance
The average prediction across all folds will be used to evaluate model overall performance on PPI classification. <br />


### References

<br />
1. Zhao C, Wang Z. (2018) GOGO: An improved algorithm to measure the semantic similarity between gene ontology terms. Sci Rep, 8, 15107.  (GOGO scores)
<br />
<br />
2. Krumsiek J, Zimmer R, Friedel CC. (2009) Bootstrapping the interactome: unsupervised identification of protein complexes in yeast. J. Comput. Biol, 16, 971–987.  (Co-localization scores)
<br />
<br />
3. Hu LZM, Goebels F, Tan JH, et al. (2019) EPIC: software toolkit for elution profile-based inference of protein complexes. Nat Method, 16, 737–742.  (Composite scores)
<br />
<br />
4. Nepusz T, Yu H, Paccanaro A. (2012) Detecting overlapping protein complexes in protein-protein interaction networks. Nat Methods, 9, 471–472.  (MMR scores used in composite scores calculation)
<br />
