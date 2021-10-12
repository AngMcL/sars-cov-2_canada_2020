#!/usr/bin/env Rscript
# pull the inferred dates for tips with incomplete dates

#usage:
# $Rscript InferredDates.R <"LSD tree folder in">
# $Rscript InferredDates.R "../02_trees/ft_root_res_time"

#lib setup 
library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(plyr)
library(ggplot2)
# BiocManager::install("ggtree")
library(ggtree)
# BiocManager::install("ggimage")
library(ggimage)
library(phytools)
library(phangorn)
library(forcats)
library(lubridate)
library(ggstance)
library(treeio)
library(gtools)

#manual input
trees.in.f<-"../02_trees/ft_root_res_time"

trees.in<-list.files(trees.in.f, pattern="fasta.timetree.nex$",full.names = T)
trees.in<-mixedsort(trees.in)

#add'l inputs
b<-length(trees.in)
BOOTS<-1:b
#sublineage summary 
sum.b.in<-paste("DF/sublin.b",BOOTS,".csv",sep="")
#meta from AncestralReconstruction
meta.b.in<-paste("DF/meta.b",BOOTS,".csv",sep="")
  
#output setup
meta.out<-paste("DF/InfDates_meta.b",BOOTS,".csv",sep="")

for (i in 1:b){ #loop across all boots
  #sublineage df
  sum.b<-read.csv(sum.b.in[i])
  meta.b.inf<-read.csv(meta.b.in[i])

  ## extract the dates and uncertainty
  t<-treeio::read.beast(trees.in[i])
  td<-get.data(t) #data on inferred internal node dates and incomplete dates on tips
  tp<-get.tree(t) #non-binary nexus time tree
  
  #### Pull the information about the inferred tip dates for incomplete dates #####
  tips.with.CI<-which(!is.na(td$CI_height) & td$node<Ntip(tp)) 
  tips.incomp<-tp$tip.label[td$node[tips.with.CI]]
  date.inf<-td$date[tips.with.CI]
  date.CI.inf<-td$CI_date[tips.with.CI]
  
  #parse the low and high CI
  date.CI.low<-c()
  date.CI.high<-c()
  for (j in 1:length(date.CI.inf)){
    date.CI.low<-c(date.CI.low,date.CI.inf[[j]][1])
    date.CI.high<-c(date.CI.high,date.CI.inf[[j]][2])
  }
  
  #make into a little df that you can then use to replace values in the meta
  incomp.df<-data.frame(tip.label=tips.incomp, date.lsd=date.inf, date.lsd.low=date.CI.low, date.lsd.high=date.CI.high)
  
  #issues with some being inferred as only "2021" (also upper, but not lower)
  #for these instances, add a manual CI
  for (j in 1:nrow(incomp.df)){
    if (incomp.df$date.lsd[j]==2021){
      incomp.df$date.lsd[j]<-"2021-01-01"
      incomp.df$date.lsd.high[j]<-"2021-01-15"
    }
  }
  
  #join the dates onto meta
  # meta.b.inf$date.inf[meta.b.inf$tip.label %in% incomp.df$tip.label]
  meta.b.2<-left_join(meta.b.inf, incomp.df,by="tip.label")
  meta.b.2[meta.b.inf$tip.label %in% incomp.df$tip.label,
           c("date","date.lsd","date.lsd.low","date.lsd.high")]
  
  #Make another column that either has the full date already or has the lsd date, "date.lsd.full"
  #need a column that combines the dates if complete and lsd-inferred dates
  meta.b.2$date.lsd.full<-meta.b.2$date
  sum(table(meta.b.2$date.lsd.full))
  
  for (j in 1:nrow(meta.b.2)){
    dt<-unlist(strsplit(meta.b.2$date[j],split="-"))
    #if initial date was incomplete, use the lsd inferred date
    if(length(dt<3) & !is.na(meta.b.2$date.lsd[j])){
      meta.b.2$date.lsd.full[j]<-meta.b.2$date.lsd[j]}
  }
  # sum(table(meta.b.2$date.lsd.full))
  # table(meta.b.2$date.lsd.full) #sweeeet
  
  #export a new meta version that has inferred dates in them (use this for TMRCA)
  write.csv(meta.b.2,meta.out[i])
  
}

