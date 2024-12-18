
Link to overview file:

https://docs.google.com/document/d/1mWBy-rSmfUf1kbHbQ-7C6SYpg5zamuJTPoh0yb8Vusk/edit?tab=t.0

# Directories

- `/data/db/` contains the benchmarking dataset
    - `data/db/IDs_msas.csv` are all the IDs for which there are MSAs for AF2. When running AF3, we can also restrict ourselves to these IDs. Note that the IDs are in 4 letter form but in some of the paths they might appear as 4 letter + "complex", such as *7ox2complex*.
    - `data/db/lightDb.txt` contains a csv with columns *pdb ID,A chain(s),H chain(s),L chain(s)*. Note that in some cases L chain(s) is NA or A chain(s) contains more than one entry.
    - `data/db/complexFastas` contains the fasta files for each ID (resp. complex).
 - `/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/data/MSAs/` contains the MSAs for all IDs in `IDs_msas.csv`. Note that the individual directories contain a "complex" suffix, e.g. `/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/data/MSAs/7ox2complex`

# Output

## AF 2 based output

For the AF 2 based models, you output should ideally follow the following pattern:

- pkl: `${OUTPUTDIR}/${RUNNAME}/${ID}/result_model_${i}_multimer_v3_pred_${k}_${RUNNAME}.pkl`
- pdb: `${OUTPUTDIR}/${RUNNAME}/${ID}/unrelaxed_model_${i}_multimer_v3_pred_${k}_${RUNNAME}.pdb`

where `OUTPUTDIR` is the output directory, `RUNNAME` is the name of the run (could be something like `afsample2` or `subsampling-32-64`), `i` is the number of the model, `k` is the model number (starting at 0) and `ID` is the ID in 4 letter + "complex" shape, e.g. *7ox2complex*. 
You can use the `/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/scripts-berzelius/check_run.sh` script to check that your output is complete.

## AF 3 based output

For AF 3, the output should ideally follow a similar pattern. Something along the lines of `..._${i}_..._${k}_..._${RUNNAME}` where `k` would be the prediction number (corresponding to the model seed), and `i` is the diffusion sample number. 
I have not looked at other AF 3 clones but ideally the output should follow the same pattern. If possible (I have not checked), I would prefer the output to be .pdb and .pkl files instead of .cif and .json (if they have a flag for that).



In both cases, if you would rather use a different notation or the program requires it, that is fine too (for instance having the runname in the base directory name or something like that). 
What matters is that I can run a script (see `/proj/berzelius-2021-29/users/x_safro/git/abag-benchmark-set/scripts-tetralith/create_samples_per_id.sh`), where I can collect for each model from your output the pdb file (or cif file in the case of AF3), pkl file (or json file in the case of AF3), features.pkl file (in the case of AF2), the prediction number `k`, the model number `i`.


# What to run

## AF2 based models

- afsample2default: 5 x 20 = 200 models (the 5 different model weights and 20 seeds per id. Note that for AF2 the seeds are randomized so you dont have to specify the seed.) (use method 'afsample2') (https://github.com/iamysk/AFsample2)
- afsample2speachaf: 5 x 20 = 200 models (use method 'speachaf')
- afsample: 5 x 20 = 200 models (alphafold 2.1 version available here, https://github.com/bjornwallner/alphafoldv2.2.0. One has to update afsample to use alphafold 2.3 first though.)

## AF3 based models

- af3default: 5 x 20 = 200 models (5 diffusion samples and 20 seeds per id)
- boltzdefault: 5 x 20 = 200 models (5 diffusion samples and 20 seeds per id)
- chai1default: 5 x 20 = 200 models (5 diffusion samples and 20 seeds per id)

Please make sure to document your runs (e.g. which settings and parameters you use for each of the different runs)

# Timeline

It would be great if we could finish the benchmarking until say middle of January or third week of January.

Depending on how AF3 performs, it might be interesting to run AF3 with msa subsampling (subsampling increases performance of AF2 by a fair bit). However, I am not sure if we are allowed to modify AF3 so that is why we would need Boltz too (since it is opensource).

Depending on how much time we have we might want to add other benchmarking models too but we see how it goes.


# Raum für Gedanken (space for thoughts)