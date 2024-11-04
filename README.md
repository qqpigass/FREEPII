# FREEPII
## Source code for FREEPII <br />
>> Currently only human data are provided as examples
 <br />
To generate predicted PPIs, follow FREEPII_SEC2-heavy_cv1_example.ipynb in the code folder, and use MCL_Clusters.ipynb to generate predicted clusters.
<br />
<br />
To evaluate model performance on PPI classification, follow code in Performance_PPI.R.
<br />
<br />
To measure the quality of predicted clusters, use GO_score.ipynb, CoLocalization_score_Complex.R, Composite_score_Complex.R in the code folder to calculate GOGO scores, co-localization scores, and composite scores. <br />
<br />
<br />

## References:
<br />
Zhao C, Wang Z. (2018) GOGO: An improved algorithm to measure the semantic similarity between gene ontology terms. Sci Rep, 8, 15107. (GOGO scores)
<br />
<br />
Krumsiek J, Zimmer R, Friedel CC. (2009) Bootstrapping the interactome: unsupervised identification of protein complexes in yeast. J. Comput. Biol, 16, 971–987.  (Co-localization scores)
<br />
<br />
Hu LZM, Goebels F, Tan JH, et al. (2019) EPIC: software toolkit for elution profile-based inference of protein complexes. Nat Method, 16, 737–742. (Composite scores)
<br />
Nepusz T, Yu H, Paccanaro A. (2012) Detecting overlapping protein complexes in protein-protein interaction networks. Nat Methods, 9, 471–472. (MMR scores used in composite scores calculation)
<br />
