# 01_subSample

## Purposes
* Subsample (all lineages) alignment using the global country-specific case counts (from the day of most recent sample) by month
    * Optional: Subsample Canadian sequences proportionally to provinces' burden
    * Export alignments and corresponding metadata for b boostraps

## Dependencies
See parent directory

## Usage in terminal
* Change directory to 01_subSample

* subsample sequences, specifying the number of global sequences to sample
    $Rscript scripts/subsample_n.R <"seq_in.fasta"> <"meta_in.csv"> <number of bootstraps> <number of total samples to take, must be greater than n Canada>
    $Rscript scripts/subsample_n.R "../00_cleanData/masked/mask_clean_fake.fasta" "../00_cleanData/cleaned/clean_fake_meta.csv" 10 150
    
    * Outputs:
        * for each bootstrap: 
            * bootsamples/subsamp_meta_*.csv 
            * bootsamples/subsamp_align_*.fasta
