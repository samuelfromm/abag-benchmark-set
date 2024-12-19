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

For details on what the respective parameters exactly mean see the description for https://github.com/samuelfromm/AADaM-fork. 

The resulting dataset consists of 115 antibody-antigen structures. Each complex consists of one or several antigen chains, a heavy chain and a light chain (if applicable). 

## Dateset structure


data/
    db/
        antibodyFastas                  
        antibodyFastasColonSep
        antigenFastas
        antigenFastasColonSep
        complexFastas               # contains the fasta files for each (interface) complex
        complexFastasColonSep
        fullDb.json                 # contains some metadata about each (interface) complex
        fullDb.pkl
        IDs.csv                     # contains the complex IDs
        IDs_msas.csv
        IDs_working.csv
        lightDb.txt                 # a csv file with columns "pdb ID,A chain(s),H chain(s),L chain(s)"
        structures/                 # pdb files for the ground truth complexes
        structures_filtered/
    MSA/
        (contents of MSA directory here)
    MSA.zip

## MSA generation

The MSAs were generated using the standard We run Alphafold 2.3.2 pipeline (specifically commit `f251de6`). Out of the 115 complexes, one complex did not run due to a failure of HHblits.

# MODEL GENERATION

## AF2 based 

For the AF2 based methods, we use Alphafold 2.3.2 (specifically commit `f251de6`) unless stated otherwise.


*default*: We run 40 predictions per model with the AlphaFold 2 standard parameters, resulting in 40x5=200 models.



# References

[1] McCoy KM, Ackerman ME, Grigoryan G. A comparison of antibody-antigen complex sequence-to-structure prediction methods and their systematic biases. Protein Sci. 2024 Sep;33(9):e5127. doi: 10.1002/pro.5127. PMID: 39167052; PMCID: PMC11337930.

[2] Constantin Schneider, Matthew I J Raybould, Charlotte M Deane. SAbDab in the age of biotherapeutics: updates including SAbDab-nano, the nanobody structure tracker. Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D1368–D1372, https://doi.org/10.1093/nar/gkab1050.




 # ISSUES

- For 1 ID, 8w90, the MSA generation did not run due to a failure of HHblits
- For 4 IDs, AADaM does not run properly (i.e. 8r4q,8ezl,8cz8,7tzh: the issue is that the L chain identifier output by AADaM is the same as the H chain identifier. One could probably fix these manually by looking up the correct L chain identifier in the ground truth PDB file but I decided to ignore this and move on.)