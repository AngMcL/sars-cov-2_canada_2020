# 03_stateReconstruction

## Description
Max likelihood inference of the ancestral geographic state of internal nodes in the phylogeny, quantification of intros of SARS-CoV-2 to Canada, as well as downstream analysis and visualization of this output

## Usage
* First conduct ancestral reconstruction of geographic state for internal nodes of the time-scaled phylogenies
    $Rscript AncestralReconstruction.R <tree folder in> <clean meta in>
    * Ouputs: DF/* 
        * dfs of tips' ancestral state, meta restricted to boot, internal node states, sublineages

* Pull the inferred dates of incomplete Canadian collection dates from the timetree and merge into the metadata for each bootstrap
    $Rscript InferredDates.R <tree folder in> 
    * Ouputs: DF/InfDates_metab*.csv 

* Use meta to generate lineage and geography color schemes
    $Rscript scripts/MakeColors.R 

* Analyze the output of the ancestral reconstruction, generate figures,etc
    $Rscript scripts/AnalyzeStates.R
    * Outputs: Results/*
    
* Make phylogeny figure
    $Rscript scripts/PlotTreeFigures.R
    * note that the new ggtree version is not playing nicely
