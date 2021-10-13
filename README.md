# sars-cov-2_canada_2020
Angela McLaughlin et al. 

## Description
Scripts needed to reproduce our analysis as described in the manuscript, 'Early introductions of SARS-CoV-2 sublineages into Canada drove the 2020 epidemic', in submission at eLife. Pre-print available at medRxiv doi:10.1101/2021.04.09.21255131v3

Since we are not allowed to share GISAID data publicly, a fake dataset of SARS-CoV-2 sequences simulated from a Wuhan-hu-1 origin was run through the pipeline in lieu.

## Dependencies
Note this pipeline was developed and tested on mac OX Catalina 10.15.7 

### Programs
* R version 4.0.3
* Python 3.8.5
* pangolin v3.1.14 https://github.com/cov-lineages/pangolin
* FastTree v2.2.1 
* IQTREE-2.1.2 (includes LSD2)
* TempEst v1.5 
* minimap2 

### R packages
ape 5.4-1, Biostrings 3.1.3, phytools 0.7-70, phangorn 2.5.5, forcats 0.5.0, coronavirus 0.3.0.9000, dplyr 1.0.2, tidyr 1.1.2, plyr 1.8.6, lubridate 1.7.9.2, stringr 1.4.8, stringi 1.5.3, zoo 1.8-8, ggplot2 3.3.3, RColorBrewer 1.1-2, ochRe 1.0.0, cowplot 1.1.164, ggstance 0.3.5, ggalluvial 0.12.3, ggmosaic 0.3_3, ggtree 2.2.4, ggplotify 0.0.5, ggrepel 0.9.1, and MASS 7.3-53. Additional R packages used to generate the maps included rgeos 0.5-5, maptools 1.0-2, ggsn 0.5.0, broom 0.7.6, and rgdal 1.5-18.

### Python packages/scripts
* viralMSA.py https://github.com/niemasd/ViralMSA
* masking script https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/src/mask_alignment_using_vcf.py

### Additional files cloned from github
* nextstrain exclude list https://github.com/nextstrain/ncov/blob/master/defaults/exclude.txt
* problematic sites https://github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf

## Usage
The readme files in the individual subfolders have details on running scripts

## Thank you
* GISAID.org
* Originating and contributing laboratories who contributed to GISAID database (4.1 million and counting on 2021-10-06)
* Co-authors: Vincent Montoya, Rachel Miller, Gideon Mordecai, Canadian COVID-19 Genomics Network (CanCOGeN) Consortium, Michael Worobey, Art F. Y. Poon, Jeffrey B. Joy
* The many software and package developers
* See manuscript for full citations

## Funding sources
AM was supported by a Canadian Institutes for Health Research (CIHR) Doctoral grant and a Natural Sciences and Engineering Research Council of Canada (NSERC) CREATE scholarship. RLM was supported by an NSERC CREATE scholarship. GM was supported by the Liber Ero Fellowship Programme. MW was supported by the David and Lucile Packard Foundation. AFYP was supported by a CIHR Project Grant PJT-156178. JBJ was supported by Genome Canada BCB 287PHY grant, an operating grant from the CIHR Coronavirus Rapid Response Programme number 440371, and a CIHR variant of concern supplement. The British Columbia Centre for Excellence in HIV/AIDS also provided support.

