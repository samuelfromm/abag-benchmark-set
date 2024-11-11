# Introduction



# METHODS

## Dataset creation

In order to create the dataset, we use a slightly modified version of the Antibody Antigen Dataset Maker (AADaM) [1] available under https://github.com/samuelfromm/AADaM-fork.
As a starting point, we use all structures deposited in SAbDab on **2024-10-23 11:47:11**. The run our version of the Antibody Antigen Dataset Maker using the following parameters:

```
--inputDb=<All structures from SAbdAB downloaded --2024-10-23 11:47:11-->
--date=2021/09/30 #  training cutoff for AF 2.3 (see https://github.com/google-deepmind/alphafold/blob/main/docs/technical_note_v2.3.0.md)
--abCompSeqCut=80
--withinDatasetCut=80
--resCut=3.5
--cutoffStrict=True
--nx=True
--globalSeqID=True
--minAtomSeqresFraction=0.9
--methodsAllowed='X-RAY DIFFRACTION,ELECTRON MICROSCOPY'
```

For details on what this exactly means see the description for https://github.com/samuelfromm/AADaM-fork. 


The resulting dataset consists of 116 antibody-antigen structures. Each complex consists of one or several antigen chains, a heavy chain and a light chain (if applicable).

# Model generation

We run Alphafold 2.3.0.

### MSA generation

We run the default Alphafold pipeline for MSA generation. Out of the 116 complexes, one complex did not run due to a failure of HHblits.

### Model generation



# References

[1] McCoy KM, Ackerman ME, Grigoryan G. A comparison of antibody-antigen complex sequence-to-structure prediction methods and their systematic biases. Protein Sci. 2024 Sep;33(9):e5127. doi: 10.1002/pro.5127. PMID: 39167052; PMCID: PMC11337930.

[2] Constantin Schneider, Matthew I J Raybould, Charlotte M Deane. SAbDab in the age of biotherapeutics: updates including SAbDab-nano, the nanobody structure tracker. Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D1368–D1372, https://doi.org/10.1093/nar/gkab1050.




 # QUESTIONS

 - WHich parameters to use for global aligner
 - which antigen-type to include in after cut off set / before cutoff set? (see filterbycutoffdate)
 ['Hapten' nan 'protein | protein' 'protein' 'protein | peptide'
 'nucleic-acid' 'protein | protein | protein' 'peptide' 'carbohydrate'
 'protein | protein | protein | protein' 'peptide | peptide'
 'peptide | protein' 'protein | protein | peptide'
 'protein | protein | protein | protein | protein'
 'protein | peptide | protein' 'protein | protein | protein | peptide'
 'protein | nucleic-acid' 'unknown' 'carbohydrate | protein | protein'
 'peptide | peptide | peptide' 'peptide | protein | protein'
 'carbohydrate | protein' 'nucleic-acid | nucleic-acid'
 'nucleic-acid | nucleic-acid | nucleic-acid']

 - "Single-Chain Fragment Variable (scFv) anti-
bodies were excluded, as many methods predict the
heavy and light chains separately, but DOCKQ score can-
not pair multiple chains to a single chain. "
Should scFv be excluded?
- allow different methods in the trainign dataset? (SOLUTION, etc.)