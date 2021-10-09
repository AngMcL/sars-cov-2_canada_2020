# sars-cov-2_canada_2020
Angela McLaughlin et al. 

## Description
This repo contains the scripts needed to reproduce our analysis as described in the manuscript, 'Early introductions of SARS-CoV-2 sublineages into Canada drove the 2020 epidemic', in submission at eLife. Pre-print available at medRxiv doi:10.1101/2021.04.09.21255131v3

## Dependencies
Note this pipeline was developed and tested on mac OX Catalina 10.15.7 (check work desktop also)

### Programs
* R version 4.0.3
* Python 3.8.5
* FastTree v2.2.1 
* IQTREE-2.1.2 (includes LSD2)
* TempEst v1.5 

### R packages:
ape 5.4-1 49, Biostrings 3.1.3 50, phytools 0.7-70 51, phangorn 2.5.5 52, forcats 0.5.0 53, coronavirus 0.3.0.9000 8, dplyr 1.0.2 54, tidyr 1.1.2 55, plyr 1.8.6 56, lubridate 1.7.9.2 57, stringr 1.4.0 58, stringi 1.5.3 59, zoo 1.8-8 60, ggplot2 3.3.3 61, RColorBrewer 1.1-2 62, ochRe 1.0.0 63, cowplot 1.1.164, ggstance 0.3.5 65, ggalluvial 0.12.3 66, ggmosaic 0.3_3 67, ggtree 2.2.4 68, ggplotify 0.0.5 69, ggrepel 0.9.1 70, and MASS 7.3-53 46. Additional R packages used to generate the maps included rgeos 0.5-5 71, maptools 1.0-2 72, ggsn 0.5.0 73, broom 0.7.6 74, and rgdal 1.5-18 75.

### Python packages:
* viralMSA.py https://github.com/niemasd/ViralMSA
* minimap2

### Additional files from cloned repos
* nextstrain: https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt
* problematic sites: https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf

## Usage
The readme files in the individual subfolders 00 - 05 have details on running the scripts

## Thank you to.. 
* GISAID.org
* Originating and contributing laboratories who contributed to GISAID database (4.1 million and counting on 2021-10-06)
* Co-authors, including memebers of Canadian COVID-19 Genomics Network (CanCOGen) Consortium
* The many software and package developers; see manuscript for full citations

