# 01_subSample

## Purposes
* Subsample (all lineages) alignment using the global country-specific case counts (from the day of most recent sample) by month
    * Optional: Subsample Canadian sequences proportionally to provinces' burden
    * Export alignments and corresponding metadata for b boostraps

## Dependencies
See parent directory

## Usage in terminal
* Change directory to 01_subSample
```console    
$cd 01_subSample
```
* Subsample sequences, specifying the number of global sequences to sample
   * Outputs for each bootstrap: bootsamples/subsamp_meta_\*.csv, bootsamples/subsamp_align_\*.fasta
```console    
#generally,
$Rscript scripts/subsample_n.R <"seq_in.fasta"> <"meta_in.csv"> <number of bootstraps> <number of total samples to take, must be greater than n Canada>
#for example,
$Rscript scripts/subsample_n.R "../00_cleanData/masked/mask_clean_fake.fasta" "../00_cleanData/cleaned/clean_fake_meta.csv" 10 150
```    

## References
**Provincial COVID-19 cases**  
Public Health Agency of Canada. (2021). Coronavirus disease 2019 (COVID-19): Epidemiology update. Government of Canada. https://health-infobase.canada.ca/covid-19/epidemiological-summary-covid-19-cases.html?stat=num&measure=active#a2 (accessed 1 April 2021).

**Population counts**  
Statistics Canada. (2021). Canada’s population clock (real-time model). https://www150.statcan.gc.ca/n1/pub/71-607-x/71-607-x2018005-eng.htm (accessed 7 April 2021).
