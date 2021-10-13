# 00_cleanData

## Purposes
* Clean GISAID data
    * remove sequences with >10% gaps, N's, ambigs; non-human animals and environmental samples; temporal outliers; exclude.txt list from Nextstrain; incomplete dates (must be yyyy-mm-dd, except for Canadian sequences)
    
* Clean GISAID metadata in same manner

* Align and mask clean GISAID data

## Dependencies
See parent directory

## Usage
* Add sequences and metadata to the GISAID folder or run the example dataset, 'fake.fasta' as shown for scripts below

* Change directory to 00_cleanData
```console
foo@bar:~$cd 00_cleanData
```

* Align sequences 
   * Output: aligned/fake.fasta.aln 
```console
foo@bar:~$ViralMSA.py -s GISAID/fake.fasta -r scripts/wuhan-hu-1.fasta -o aligned -a minimap2 -t 4 -e <your-email>
```

* Clean data
   * Outputs: cleaned/clean_fake.fasta, cleaned/clean_fake_meta.csv
```console
foo@bar:~$Rscript scripts/cleanAlignedGISAID.R aligned/fake.fasta.aln GISAID/fake_meta.tsv
```

* Run pangolin on cleaned sequences pre-masking
    * Output: cleaned/pangolin_results.csv
    * Note that fake example mostly falls in basal lineages B and A.1 because the simulated sequences have mutations added randomly across the genome (starting with wuhan-hu-1), so we don't expect them to have any set of lineage-characteristic mutations

* Mask problematic sites in the alignment, as per de Maio et al.
    * Note that need to make sure the reference_id matches wuhan-hu-1 sequence
```console
foo@bar:~$python3 scripts/mask_alignment_using_vcf.py -i cleaned/clean_fake.fasta -o masked/mask_clean_fake.fasta -v scripts/problematic_sites_vcfv5.vcf -n "n" --reference_id "China/fake_1/2020/EPI_ISL_FAKE_1/2019-12-26"
```
