# 03_stateReconstruction

## Description
Max likelihood inference of the ancestral geographic state of internal nodes in the phylogeny, quantification of intros of SARS-CoV-2 to Canada, as well as downstream analysis and visualization of this output

## Usage
* Conduct ancestral reconstruction of geographic state for internal nodes of the time-scaled phylogenies
    * Ouputs: DF/\*, dataframes of tips' ancestral state, meta restricted to boot, internal node states, sublineages
```console
$Rscript AncestralReconstruction.R <tree folder in> <clean meta in>
```

* Pull the inferred dates of incomplete Canadian collection dates from the timetree and merge into the metadata for each bootstrap
    * Ouputs: DF/InfDates_metab\*.csv 
```console
$Rscript InferredDates.R <tree folder in> 
```

* Use meta to generate lineage and geography color schemes. Also, cross-references lineages' aliases and groups lineages into groups: A\*, B\*, B.1\*, B.1.1\*
    * Ouputs: DF/lineagecolors.tsv, DF/globalcolors.tsv, DF/lineageGroups.csv
```console
$Rscript scripts/MakeColors.R 
```

* Analyze the output of the ancestral reconstruction, generate figures,etc
    * Outputs: Results/\*
```console
$Rscript scripts/AnalyzeStates.R
```

* Make phylogeny figure
    * note that the new ggtree version is not playing nicely with dplyr and other pkgs
```console
$Rscript scripts/PlotTreeFigures.R
```
