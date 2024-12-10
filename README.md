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

We run Alphafold 2.3.2.

### MSA generation

We run the default Alphafold pipeline for MSA generation. Out of the 116 complexes, one complex did not run due to a failure of HHblits.

### Model generation

#### AF default

We run 40 predictions per model with the AlphaFold 2 standard parameters, resulting in 40x5=200 models.

# References

[1] McCoy KM, Ackerman ME, Grigoryan G. A comparison of antibody-antigen complex sequence-to-structure prediction methods and their systematic biases. Protein Sci. 2024 Sep;33(9):e5127. doi: 10.1002/pro.5127. PMID: 39167052; PMCID: PMC11337930.

[2] Constantin Schneider, Matthew I J Raybould, Charlotte M Deane. SAbDab in the age of biotherapeutics: updates including SAbDab-nano, the nanobody structure tracker. Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D1368–D1372, https://doi.org/10.1093/nar/gkab1050.




 # QUESTIONS

- For 1 ID, the MSA generation did not run
- For 4 IDs, AADaM does not run properly (i.e. 8r4q,8ezl,8cz8,7tzh there is an with the L chain identifier being the same as the H chain identifier, but the complex fasta lists different identifiers)
- [default, norecycles] For 2 IDs, dockq returns an empty dictionary in some cases (i.e. 8hbi,8f8w) (NOTE: This can possibly be fixed)