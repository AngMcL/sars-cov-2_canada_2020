# 02_trees

## Description
For all boostrap subsamples, infer approx max likelihood trees using FastTree. Then root the trees on Wuhan-hu-1 (or earliest A sample), calculate residuals in root-to-tip regression in Tempest, remove temporal outliers (3 or more stdev from the mean residual OR pendant edge longer than 12 mutations). Export fasta, tree, and datefile in format for LSD2. Timescale trees using LSD2.

## Usage in terminal
* change directory to 02_trees

* Generate fasttree for each bootstrap
    * output: ft/\*.tre
```console
$FastTree -gtr -nt ../01_subSample/bootsamples/subsamp_align_fake150_1.fasta > ft/subsamp_align_fake150_1.tre
```
    
* Root the trees
    * output: ft_root/\* .tre
    * note EPI_ISL_FAKE_1 represents Wuhan-hu-1
```console
#generally,
$Rscript scripts/Root.R <"folder with trees to root"> <"rooting accession ID">

#for example,
$Rscript scripts/Root.R "ft" "EPI_ISL_FAKE_1"
```

* Calculate residuals of divergence over time regression; read each rooted tree into TempEst v1.5 (see dependencies in parental folder), parse dates, export data as .tsv. Do not re-root
    * output: ft_root/\*.tsv

* Remove temporal outliers from the tree and alignment, and generate date file for LSD2
    * output: ft_root_res/\*
```console
#generally,
$Rscript scripts/RemoveOutliers.R <"folder with rooted trees"> <"folder with corresponding fasta">

#for example,
$Rscript scripts/RemoveOutliers.R "ft_root" "../01_subSample/bootsamples"
```

* If running tree files >1000 tips, zip and run on computing cluster
```console
$tar -czf ft_root_res.tar.gz ft_root_res/
```

* Run IQTREE2 and LSD2: infer incomplete dates and time-scale tree
    * notes:
        * the nt and mem may need to be tweaked for individual machine
            * for 50000 tips on the computing cluster, -nt 64 -mem 64G (seems to be close to upper limit)
            * for 150 tips run locally, -nt AUTO
        * normal outlier (-o) for wuhan-hu-1 without fake data is "China/2019-12-26/EPI_ISL_402125"
        * may have to change to suit the location of iqtree2 on your computer
        * GTR+I+R3 model selection based Rob Lanfear SARS-CoV-2 tree comparisons
            * https://github.com/roblanf/sarscov2phylo/blob/master/tree_estimation.md
```console
$cd ft_root_res/

#generally,
$/usr/bin/iqtree-2.1.2-Linux/bin/iqtree2 -s <final.fasta> --date <dateFile.txt> -te <starting_fasttree.tre> --dating LSD --clock-sd 0.2 -nt AUTO -m GTR+I+R3 --date-ci 100 -o "China/fake_1/2020/EPI_ISL_FAKE_1/2019-12-26"

#for example,
$/usr/bin/iqtree-2.1.2-Linux/bin/iqtree2 -s res_rooted_subsamp_align_fake150_10.fasta --date dates_res_rooted_subsamp_align_fake150_10.txt -te res_rooted_subsamp_align_fake150_10.tre --dating LSD --clock-sd 0.2 -nt AUTO -m GTR+I+R3 --date-ci 100 -o "China/fake_1/2020/EPI_ISL_FAKE_1/2019-12-26"
```

* Move all LSD2 output to ft_root_res_time
