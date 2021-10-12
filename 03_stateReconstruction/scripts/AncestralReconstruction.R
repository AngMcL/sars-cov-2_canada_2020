#!/usr/bin/env Rscript
# Ancestral state reconstruction of viral geography

#usage:
# $Rscript AncestralReconstruction.R <tree folder in> <clean meta in>

##objectives
# use the subsampled data/trees that have had outliers removed based on tempest analysis
# conduct ancestral reconstruction on the overall tree (not split by lineage)
# Then, identify introductions (can nodes preceded by non-can nodes)
# split each intro into a sublineage
# parse the lineage from each and then add a .1, .2, .3
# for each sublineage intro, quantify number of desc, whether provincial mixing, from where

#library setup
library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
# BiocManager::install("ggtree")
library(ggtree)
# BiocManager::install("ggimage")
library(ggimage)
library(phytools)
library(phangorn)
library(forcats)
library(gtools)

#### Setup inputs and outputs ####
args<-commandArgs(trailingOnly=TRUE)
trees.in.f<-args[1] 
clean.meta.in<-args[2]

#manual input
trees.in.f<-"../02_trees/ft_root_res_time"
clean.meta.in<-"../00_cleanData/cleaned/clean_fake_meta.csv"

#input cont'd
pango.meta<-"../00_cleanData/cleaned/pangolin_results.csv"
trees.in<-list.files(path=trees.in.f,pattern=".timetree.nex",full.names = T)
trees.in<-mixedsort (trees.in)
b<-length(trees.in) #number of bootstraps
BOOTS<-1:b

#output DFs
can.anc.out<-paste("DF/can.anc.b",BOOTS,".csv",sep="")
sum.out<-paste("DF/sublin.b",BOOTS,".csv",sep="")
meta.out<-paste("DF/meta.b",BOOTS,".csv",sep="")
can.states.out<-paste("DF/tip.anc.b",BOOTS,".csv",sep="")

#make directory if it doesn't exist
if(!dir.exists("DF/")){dir.create("DF/")}

#### read and clean metadata (pre-subsample) common to entire build ####
#all boots share the same meta
meta.all<-read.csv(clean.meta.in,header=T)

#change some colnames and reduce size
meta.all<-meta.all[,c("new.names", "gisaid_epi_isl","strain","date","region","country","division", "pangolin_lineage", "GISAID_clade")]
colnames(meta.all)<-c("tip.label", "GISAID_ID","strain","date","region","country","division", "lineageOld", "GISAID_clade")

## update pangolin lineages in meta
lookup<-read.csv(pango.meta)
lookup<-lookup[,c(1,2)]
colnames(lookup)<-c("tip.label","Lineage")
meta.all<-left_join(meta.all,lookup,"tip.label")

#some anomalies that are supposed to be in QC
if (any(meta.all$division=="Canada")){
  meta.all$division[which(meta.all$division=="Canada")]<-"Quebec" 
}

#remove duplicates if needed (shouldn't be any)
if (any(duplicated(meta.all$GISAID_ID))){
  meta.all<-meta.all[-dups,]
}

#### Loop workflow through each boot ####

for (i in 1:b){
  #read output from iqtree
  tree<-read.nexus(trees.in[[i]])
  
  #force binary (necessary for recontstr.)
  #note that this is a random resolution and will result in different trees run to run 
  set.seed(111)
  tree<-multi2di(tree)
  
  #need to eliminate branch length==0 with v v small values
  tree$edge.length[tree$edge.length==0]<-max(nodeHeights(tree))*1e-8
  tree$edge.length[tree$edge.length<0]<-max(nodeHeights(tree))*1e-8

  #write for future use
  if(!file.exists(paste(trees.in[[i]],".binary.tre",sep=""))){
    write.tree(tree,paste(trees.in[[i]],".binary.tre",sep="")) 
  }
  
  if(!file.exists(paste(trees.in[[i]],".binary.nex",sep=""))){
    writeNexus(tree,paste(trees.in[[i]],".binary.nex",sep="")) 
  #for some reason, this file doesn't appear to be binary although the newick above is
  }
  
  ## Exclude meta not represented in bootstrap tree
  # tree$tip.label[which(!tree$tip.label %in% meta.b$tip.label)]
  these<-which(!meta.all$tip.label %in% tree$tip.label)
  if (length(these)>0){
    meta.b<-meta.all[-these,]
  }

  ## Add a state column, which will be country for most but for canadian rows add prov also
  for (j in 1:nrow(meta.b)){
    meta.b$state[j]<-meta.b$country[j]
    if (meta.b$state[j]=="Canada") {meta.b$state[j]<-paste(meta.b$state[j], meta.b$division[j], sep="_")}
  }
  
  #make sure meta.b in same order as the tree
  meta.b<-meta.b[order(match(meta.b[,"tip.label"], tree$tip.label)),]
  
  #### Use trees to reconstruct ancestral state ####
  ## run ace, ML method, equilib
  anc.b<-ace(x=meta.b$state,tree,type="discrete",method="ML",model="ER")
  
  #pull the likelihood of each state at each node into a df
  anc.df.b<-as.data.frame(anc.b$lik.anc)
  
  #### Identify Canadian sublineages ####
  #ie internal nodes with non-Can parental node
  #use the anc object to identify which nodes have majority Can support
  #add node number in order, see http://www.phytools.org/Cordoba2017/ex/8/Anc-states-discrete.html
  n.tip.b<-Ntip(tree)
  anc.df.b$node <- n.tip.b + 1:(n.tip.b-1) 
  
  #pull the majority (Max Likeli state) for each node 
  anc.df.b$node.state<-NA
  anc.df.b$node.lik<-NA
  
  for (j in 1:nrow(anc.df.b)){
    max<-max(anc.df.b[j,1:(ncol(anc.df.b)-3)])
    maxcol<-which(anc.df.b[j,]==max)[1]
    anc.df.b$node.state[j]<-colnames(anc.df.b)[maxcol]
    anc.df.b$node.lik[j]<-anc.df.b[j,maxcol]
  }
  
  #ID the Canadian internal nodes and quantify
  can<-str_which(anc.df.b$node.state,"Canada")
  # length(can); length(tree$tip.label) 


  ##Query states for canadian nodes parents and desc
  #prime col number for locations
  colz<-ncol(anc.df.b)-3 #exclude node info
  
  #for each of these, look at parent internal node and descendant tips
  anc.df.b$node.par<-NA
  anc.df.b$node.desc<-NA
  
  anc.df.b$par.state<-NA
  anc.df.b$desc.state<-NA
  
  anc.df.b$desc.n<-NA
  anc.df.b$desc.accession<-NA
  
  anc.df.b$par.lik<-NA
  anc.df.b$desc.lik<-NA
  
  #loop through canadian nodes and query parent, descendants, likelihoods
  ####Note could change this to loop through all, but will result in a bigger object
  for (j in 1:length(can)){
    par<-Ancestors(tree, node=anc.df.b$node[can[j]], type="parent") #parent node ID
    desc<-as.numeric(c(Descendants(tree, node=anc.df.b$node[can[j]]))[[1]]) #desc nodes
    desc.c<-paste0(as.character(desc),collapse=", ")
    
    #nodes of desc and par
    anc.df.b$node.desc[can[j]]<-desc.c #descendents as a character string
    anc.df.b$node.par[can[j]]<-par #parent
    
    #parent's geography state
    anc.df.b$par.state[can[j]]<- anc.df.b$node.state [which(anc.df.b$node==par)]
    #descendant geography state
    anc.df.b$desc.state[can[j]]<-paste0(meta.b$state[desc],collapse=", ") #meta sorted the same as tree
    
    #number of descendants
    anc.df.b$desc.n[can[j]]<-length(desc)
    #Accession ID of descendants
    anc.df.b$desc.accession[can[j]]<-paste0(meta.b$GISAID_ID[desc],collapse=", ")
    
    #likelihood of desc and par
    anc.df.b$par.lik[can[j]]<-anc.df.b$node.lik [which(anc.df.b$node==par)]
    anc.df.b$desc.lik[can[j]]<-paste0(anc.df.b$node.lik[desc],collapse=", ")
  }

  ## Subset state object to Can states only
  can.anc.b<-anc.df.b[can,]
    
  #check
  # table(can.anc.b$node.state)
  
  #exclude rows if parent node was also Canadian (not an intro to Canada)
  #quantify interprovincial transm overall in subseq. script
  # nrow(can.anc.b)
  if (any(str_detect(can.anc.b$par.state, "Canada"))) {
    can.anc.b<-can.anc.b[-str_which(can.anc.b$par.state, "Canada"),]
  }
  # nrow(can.anc.b)
  
  # ADD a new column for whether nested or not, then can exclude on this basis later if desired
  can.anc.b$nested<-NA
  for (j in 1:nrow(can.anc.b)){
    dd<-can.anc.b$node.desc[j]
    if (j>1){
      #if all node descendants are found in any other rows, then nested
      if (mean(str_detect(can.anc.b$node.desc[-j], dd))>0) can.anc.b$nested[j]<-"nested"
    }
  }
  
  #make this into a df with a row for each introduction 
  sum.b<-can.anc.b[,c("node","node.state","node.lik","node.par","par.state","par.lik","desc.n","nested","node.desc","desc.state", "desc.accession","desc.lik")]
  
  #summarize desc state into one character
  for (j in 1:nrow(sum.b)){
    t<-table(unlist(strsplit(sum.b$desc.state[j],split=", ")))
    des<-c()
    for (k in 1:length(t)){
      n<-paste(names(t)[k], ": ",t[k],sep="")
      des<-c(des,n)
    }
    sum.b$desc.state[j]<-paste0(des,collapse="; ")
  }
  
  colnames(sum.b)<-c("Node","Node.Location","Node.Likelihood","Parent.Node","Parent.Location","Parent.Likelihood","Number.Descendants","Nested","Descendant.Node","Descendant.Location","Descendant.AccessionID","Descendant.Likelihood")
  
  #order rows in descending order of n.desc and name them!
  sum.b<-sum.b[with(sum.b, order(Number.Descendants, decreasing = TRUE)),]
 
  ## Assign lineages to sum.b sublineages using pango updated Lineage
  # Note: moved sublineage naming into analyze states because we want to know the first sample date (which requires pulling the LSD-inferred dates)
  
  #Give each one a majority lineage designation
  sum.b$Lineage<-NA
  sum.b$Lin.Support<-NA #percent of descendants with same lineage #note could have derivatives inside
  
  #Pull the majority lineage and the support for each type within
  for (j in 1:nrow(sum.b)){
    epis<-unlist(strsplit(sum.b$Descendant.AccessionID[j],split=", "))
    lins<-meta.b$Lineage[meta.b$GISAID_ID %in% epis]
    t<-table(lins)
    all<-sum(t)
    t<-rev(sort(t)) #sort by most frequent
    des<-c()
    for (k in 1:length(t)){
      n<-paste(names(t)[k], ": ",round(t[k]/all*100,digits=2),"%",sep="")
      des<-c(des,n)
    }
    sum.b$Lin.Support[j]<-paste0(des,collapse="; ")
    sum.b$Lineage[j]<-names(rev(sort(t)))[1] 
  }
  
  #### Identify Canadian tips' direct ancestral (parental) state ####
  #identify Can tips, labs, and states
  can.tips<-str_which(tree$tip.label,"Canada")
  can.tiplabs<-tree$tip.label[can.tips]
  can.tip.state<-meta.b$state[match(can.tiplabs,meta.b$tip.label)]
  
  #identify parental nodes
  can.par.nodes<-Ancestors(tree,node=can.tips,type="parent")
  
  #start a df for all this
  can.states.df<-data.frame("tip.id"=can.tips, "tip.label"=can.tiplabs, "tip.state"=can.tip.state,"par.id"=can.par.nodes,"par.state"=NA, "par.lik"=NA)
  
  #go through each row, pull the parent state and its likelihood
  for (j in 1:nrow(can.states.df)){
    par<-can.states.df$par.id[j]
    #lookup the node state in anc.df
    can.states.df$par.state[j]<-anc.df.b$node.state[anc.df.b$node==par]
    can.states.df$par.lik[j]<-anc.df.b$node.lik[anc.df.b$node==par]
  }
  #### write these objects #### 
  
  #NO OVERWRITES
  if (!file.exists(anc.df.out[i])){
    write.csv(can.anc.b,can.anc.out[i],row.names = F)
    write.csv(sum.b, sum.out[i],row.names = F)
    write.csv(meta.b,meta.out[i],row.names = F)
    write.csv(can.states.df, can.states.out[i],row.names = F)
  }

}# final loop closure over boot i
