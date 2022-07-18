# Genomic epidemiology of the first two waves of SARS-CoV-2 in Canada

## Description
Reproducible example of the analyses described in the manuscripts, 

McLaughlin, A., Montoya, V., Miller, R. L., Mordecai, G. J., Canadian COVID-19 Genome Network Consortium (CanCOGen), Worobey, M., Poon, A. F. Y., & Joy, J. B. (2022). Genomic epidemiology of the first two waves of SARS-CoV-2 in Canada. eLife. https://doi.org/10.7554/eLife.73896.

McLaughlin, A., Montoya, V., Miller, R. L., Mordecai, G. J., Worobey, M., Poon, A. F. Y., & Joy, J. B. (2021). Early and ongoing importations of SARS-CoV-2 in Canada [pre-print]. medRxiv. https://www.medrxiv.org/content/10.1101/2021.04.09.21255131v3.

Since we are not allowed to share GISAID data publicly, a fake dataset of SARS-CoV-2 sequences simulated from a Wuhan-hu-1 origin was run through the pipeline as a reproducible example.

## Dependencies
This pipeline has been tested on macOS Catalina 10.15.7 and Monterey 12.0.1.

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

