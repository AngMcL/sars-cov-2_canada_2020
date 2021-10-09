#!/usr/bin/env Rscript
# Root tree and remove temporal outliers

#usage:
# $Rscript Root.R "ft" "EPI_ISL_FAKE_1"

#setup libraries
library(stringr)
library(ape)
library(tidyverse)
library(gtools)

#### Inputs and outputs here ####
args<-commandArgs(trailingOnly=TRUE)
treebase<-args[1] #folder with trees
root.accession<-args[2]

## Or, change inputs manually
# treebase<-"ft"
# root.accession<-"EPI_ISL_FAKE_1" #"EPI_ISL_402125" #use the latter if using real data

#setup treefiles
trees.in<-list.files(treebase,pattern = ".tre",full.names = T) 
#sort these in numeric order
trees.in<-mixedsort (trees.in)

b<-length(trees.in) #number of bootstraps

#filebase
fb<-unlist(str_split(last(unlist(str_split(trees.in[1],"/"))),"_"))
fb<-fb[-length(fb)]
filebase<-paste0(fb,collapse="_")

#set up outputs
tree.out.base<-paste("ft_root/rooted_",filebase,"_",sep="") #add boot and .tre at the end

# On rooting 
# In the manuscript, we followed recommendation (Rambaut 2020 Nature) to root on earliest lineage A viruses:Wuhan/WH04/2020 (EPI_ISL_406801), sampled on 5 January 2020, 
# however here, we root on Wuhan-Hu-1 (GenBank accession no. MN908947) sampled on 26 December 2019, because this fake dataset doesn't have EPI_ISL_406801
#read in trees into a list
trees<-replicate(n=b, vector)
for (i in 1:b){
  trees[[i]]<-read.tree(file=trees.in[i])

  #force into binary, ie resolve polytomies randomly
  if (!is.binary(trees[[i]])){
      trees[[i]]<-multi2di(trees[[i]])
  }

  #earliest lineage A virus:Wuhan/WH04/2020 sampled 2020-01-05,
  ##### According to Rambaut, MRCA of whole outbreak as it shares two muts not found in wuhan-hu-1
  # linA<-"EPI_ISL_406801"
  # linA.tip<-str_which(trees[[i]]$tip.label,linA)
  
  #earliest lineage B virus: "Wuhan-Hu-1" (Genbank accession MN908947) sampled on 2019-12-26
  linB<-root.accession
  linB.tip<- grep(paste("\\b",linB,"\\b",sep=""), trees[[i]]$tip.label)
  #find the MRCA node of lineages A and B
  # AB.mrca<-getMRCA(trees[[i]],tip=c(linB.tip,linA.tip))
  
  #re-root it on outgroup
  trees[[i]]<-ape::unroot(trees[[i]])
  trees[[i]]<-ape::root(trees[[i]],outgroup=linB.tip,resolve.root=TRUE)
  #is.rooted(trees[[i]])
  
  # export trees
  write.tree(trees[[i]],file=paste(tree.out.base,i,".tre",sep=""))
}

