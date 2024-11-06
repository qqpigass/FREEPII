## Source code for FREEPII <br />
>> Currently only human data are provided as examples <br />

## Basic data
The human CF-MS datasets (PXD002892, PXD014820, and PXD015406) were downloaded from Zenodo (doi: 10.5281/zenodo.4106578). <br />
Among the various files corresponding to different protein quantification strategies provided by the authors, we selected those containing iBAQ intensity of chromatograms for further analysis. <br /> <br />

The human protein complex dataset was downloaded from the CORUM database, and the human protein sequences were retrieved from UniProt database. <br /><br />

We retrieved the semantics and relationships between GO terms from the GO Consortium released on November 4, 2022. <br /><br />

To generate predicted PPIs, follow FREEPII_SEC2-heavy_cv1_example.ipynb in the code folder, and use MCL_Clusters.ipynb to generate predicted clusters.
<br />
<br />
To evaluate model performance on PPI classification, follow code in Performance_PPI.R.
<br />
<br />
To measure the quality of predicted clusters, use GO_score.ipynb, CoLocalization_score_Complex.R, Composite_score_Complex.R in the code folder to calculate GOGO scores, co-localization scores, and composite scores. <br />
<br />
<br />

## References
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
