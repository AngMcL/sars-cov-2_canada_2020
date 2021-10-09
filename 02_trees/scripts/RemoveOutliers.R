#!/usr/bin/env Rscript
# Root tree and remove temporal outliers
# export fasta and tree with no outliers (input for iqtree2/lsd2)
# export a dateFile corresponding to tips for lsd2

#usage:
# $Rscript RemoveOutliers.R "ft_root" "../01_subSample/bootsamples"

# setup libraries
library(stringr)
library(ape)
library(treeio)
library(adephylo)
library(tidyverse)
library(gtools)

#### Inputs and outputs ####
args<-commandArgs(trailingOnly=TRUE)
treebase<-args[1] #folder with trees
fastabase<-args[2] #folder with alignments corresponding to trees

# ## Or, change inputs manually
# treebase<-"ft_root"
# fastabase<-"../01_subSample/bootsamples"

#setup treefiles, alignment files using folders
trees.in<-list.files(treebase,pattern = ".tre",full.names = T)
#sort these in numeric order
trees.in<-mixedsort (trees.in)

#number of bootstraps
b<-length(trees.in)

#setup tempest data
tempest.in<-list.files(treebase,pattern = ".tsv",full.names = T)
tempest.in<-mixedsort (tempest.in)

#alignments in (to remove outliers from)
aligns.in<-list.files(fastabase,pattern=".fasta",full.names=T)
aligns.in<-mixedsort (aligns.in)

#filebase
fb<-unlist(str_split(last(unlist(str_split(trees.in[1],"/"))),"_"))
fb<-fb[-length(fb)]
filebase<-paste0(fb,collapse="_")

#set up outputs
tree.out.base<-paste("ft_root_res/res_",filebase,"_",sep="") #add boot and .tre at the end
align.out.base<-paste("ft_root_res/res_",filebase,"_",sep="") #add boot and .fasta at the end
date.out.base<-paste("ft_root_res/dates_res_",filebase,"_",sep="") 

fold.out<-"ft_root_res/results"
if(!dir.exists(fold.out)){dir.create(fold.out)}


#### READ IN DATA #####

## Read in the trees into a list
trees<-replicate(n=b, vector)
for (i in 1:b){
  trees[[i]]<-read.tree(file=trees.in[i])
}

## Read in tempest data with residuals here
resid<-replicate(n=b,vector())
for (i in 1:b){
  resid[[i]]<-read.table(file=tempest.in[i] ,sep="\t",header=T)
  resid[[i]]$boot<-i
}

##REMOVE residual datapoint if date=0 (incomplete dates)
#can't use for analysis, so will carry forward blindly... 
for (i in 1:b){
  resid[[i]]<-resid[[i]][-which(resid[[i]]$date==0),]
}

#### PLOT RESIDUALS ####
##Look at distribution of residuals and clock in each
##  plot residuals by BOOTS
plot.residuals<-function(res){
  #root-to-tip regression
  lm<-summary(lm(res$distance ~ res$date))
  slope<-formatC(coef(lm)[[2]],format="e", digits=2)
  slope.n<-coef(lm)[[2]] #numeric
  int<-coef(lm)[[1]]
  eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
  R2<-signif(lm$adj.r.squared,digits = 3)
  
  #set positions
  max.x<-max(res$date)
  max.y<-max(res$distance)
  
  BOOT<-res$boot[1]
  
  #overlay line on scatter plot, colored by residual
  ggplot(res)+
    geom_point(aes(x=date,y=distance, color=abs(residual)))+
    scale_color_continuous(type = "viridis")+
    labs(x="Year",y="Root-to-tip distance",color="Abs.Value\nResidual",title=BOOT)+
    theme_bw()+
    geom_abline(slope=slope.n,intercept = int,color="black")+
    annotate("text",label=paste("Adj.R^2=",R2,"\n",eq,sep=""),x=max.x,y=max.y*0.9,hjust=1,color="black",size=4)
  ggsave(paste(fold.out,"/TemporalSignalResidualsIncludingOutliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
}

#run it over the list and save plots into results
lapply(resid,plot.residuals)


##Look at distribution of residuals in each
##  plot density residuals with cutoff}
plot.noabs.residuals<-function(res){
  BOOT<-res$boot[1]
  
  #calculate mean + 2 stdev for each for residual
  potential.cutoff<-mean((res$residual))+(3*sd((res$residual)))
  
  #calculate n cutoff (on negative or positive side)
  n.cutoff<-length(which(abs(res$residual)>potential.cutoff))
  
  max.x<-max((res$residual))
  max.y<-.95
  
  #density plot, colored by residual
  ggplot(res)+
    geom_density(aes(x=(residual), y=..scaled..))+
    labs(x="Residual",y="Scaled density",title=BOOT)+
    theme_bw()+
    geom_vline(xintercept=potential.cutoff,color="red",linetype=2)+
    geom_vline(xintercept=-potential.cutoff,color="red",linetype=2)+
    annotate("text",label=paste("cutoff = ",signif(potential.cutoff,digits=3),"\n", "n excluded = ",n.cutoff,sep=""),
             x=max.x,y=max.y*0.8,hjust=1,color="black",size=3)
  ggsave(paste(fold.out,"/densityResiduals_N-outliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
}

lapply(resid,plot.noabs.residuals)

## Remove residuals past a cutoff
#apply cutoff of mean(residual)+3*sd(residual)
for (i in 1:b){
  cut<-mean((resid[[i]]$residual))+(3*sd((resid[[i]]$residual)))
  outliers<-which(abs(resid[[i]]$residual)>cut) #abs() here cuts off too low and too high
  #make sure wuhan-hu-1 is not removed,need it to root the timedtree
  wu.out<-str_which(resid[[i]]$tip[outliers],"2019-12-26")
  if(length(wu.out)>0){outliers<-outliers[-wu.out]}
  #remove the outliers
  if(length(outliers)>0){
    trees[[i]]<-ape::drop.tip(phy=trees[[i]],tip=resid[[i]]$tip[outliers])
    resid[[i]]<-resid[[i]][-outliers,]
  }
}

#check
# Ntip(trees[[1]])==nrow(resid[[1]])
#not the same if excluded incomplete dates (0s) from resid

# note that the fake data was simulated under a strict clock so shouldn't really have temp outliers, unless made more stringent

# plot the residuals again
plot.new.residuals<-function(res){
  
  #root-to-tip regression
  lm<-summary(lm(res$distance ~ res$date))
  res$residual<-residuals(lm(res$distance ~ res$date))
  slope<-formatC(coef(lm)[[2]],format="e", digits=2)
  slope.n<-coef(lm)[[2]] #numeric
  int<-coef(lm)[[1]]
  eq<-paste("y = ",slope," x + ",signif(int,digits=3),sep="")
  R2<-signif(lm$adj.r.squared,digits = 3)
  
  #set positions
  max.x<-max(res$date)
  max.y<-max(res$distance)
  min.y<-0.0002
  
  BOOT<-res$boot[1]
  
  #overlay line on scatter plot, colored by residual
  ggplot(res)+
    geom_point(aes(x=date,y=distance, color=abs(residual)))+
    scale_color_continuous(type = "viridis")+
    labs(x="Year",y="Root-to-tip distance",color="Abs.Value\nResidual",title=BOOT)+
    theme_bw()+
    geom_abline(slope=slope.n,intercept = int,color="black")+
    annotate("text",label=paste("Adj.R^2=",R2,"\n",eq,sep=""),x=max.x,y=min.y*0.9,hjust=1,color="black",size=3)
  ggsave(paste(fold.out,"/TemporalSignalResidualsExcludingOutliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
  
}

lapply(resid,plot.new.residuals)

## calculate the pendant edge length (last branch leading to tip) for all tips 
#PE Function (from Jeffrey Joy)
pendant.edge <-function(tree, scale=F){
  
  if(is.rooted(tree)==FALSE)
    warning("A rooted phylogeny is required for meaningful output of this function", call.=FALSE)
  if(scale==TRUE){
    
    #Scale tree to have unit depth (for an ultrametric tree) or scale all branches to unit length (for an additive tree)
    if(is.ultrametric(tree)==TRUE)
      tree<- rescaleTree(tree, 1) else
        tree$edge.length<- tree$edge.length/sum(tree$edge.length)
  }
  
  edge.index<-which(tree$edge[,2] %in% 1:length(tree$tip.label))
  w<- as.data.frame(tree$edge.length[edge.index])
  results<- cbind(tree$tip.label, w)
  
  names(results)<- c("Species", "w")
  results
}

## See if long branches still an issue
##  pendant edge with cutoffs}
####THIS CHUNK NOT ADAPTED FOR MULTIPLE BOOTS YET####
for (i in 1:b){
  pend<-pendant.edge(trees[[i]])
   #max number of muts expected on one branch if 6 months under strict clock 8e-4
  # (29903*8e-4)/12 #2 mut/month
  #if 12 muts max, what length branch in subs/site?
  pend.cutoff<-12/29903
  n.cutoff<-length(which(pend$w>pend.cutoff))
  
  max.x<-max(pend$w)
  max.y<-.95
  BOOT<-i
  
  #density plot of # mutations
  ggplot(pend)+
    geom_density(aes(x=round(w*29903), y=..scaled..))+
    labs(x="Number of mutations on pendant edge",y="Scaled density")+
    theme_bw()+
    geom_vline(xintercept=12,color="red",linetype=2)+
    annotate("text",label=paste("cutoff = ",signif(pend.cutoff,digits=3),"\n", "n excluded = ",n.cutoff,sep=""),
             x=15,y=0.95,hjust=1,color="black",size=3)
  ggsave(paste(fold.out,"/PendantEdges_N-outliers_",BOOT,".png",sep=""),height=3,width=4,units="in")
}


## Exclude long pendant edges (more than 12 muts)

#apply cutoff of 12 muts
for (i in 1:b){
  cut<-12/29903  
  pend<-pendant.edge(trees[[i]])
  
  outliers<-which(pend$w > cut) 
  if(length(outliers)>0){
    trees[[i]]<-ape::drop.tip(phy=trees[[i]],tip=pend$Species[outliers])
    pend2<-pend[-outliers,]
    resid[[i]]<-resid[[i]][-which(resid[[i]]$tip %in% pend$Species[outliers] ),]
  }
}

#check for one boot
# Ntip(trees[[length(trees)]])==nrow(pend2)
# Ntip(trees[[i]])

## Re-plot and save the residuals, excluding seqs with high residuals and long pendant edges (final strict clock)
lapply(resid,plot.new.residuals)

## exclude these outliers from the alignment
align<-replicate(n=b,vector)

for (i in 1:b){
  fast<-read.FASTA(aligns.in[i])
  
  #check
  length(which(names(fast) %in% trees[[i]]$tip.label))==length(trees[[i]]$tip.label)
  names(fast)[which(!names(fast) %in% trees[[i]]$tip.label)]
  trees[[i]]$tip.label [which(!trees[[i]]$tip.label %in% names(fast))]
  
  #now can match on name
  #only keep the keepers
  keep<-which(names(fast)%in% trees[[i]]$tip.label)
  align[[i]]<-fast[keep]
}

# test, these should be the same
for (i in 1:b){
  print(length(align[[i]])==Ntip(trees[[i]]))
}

## make sure that naming order is the same in fasta and tree (for lsd2)
for (i in 1:b){
  if(!all(names(align[[i]])==trees[[i]]$tip.label))
    #reorder the fasta to match tree
    align[[i]]<-align[[i]][order(match(names(align[[i]]), trees[[i]]$tip.label))]
    print(all(names(align[[i]])==trees[[i]]$tip.label)) #should all be T after
}

## Extract the dates for LSDin datefile

## export a datefile from this tree
dates.list<-replicate(b,vector())
for (i in 1:b){
  #extract tipnames
  nm<-trees[[i]]$tip.label
  #parse the dates
  dates.list[[i]]<-data.frame(node_name=nm,date=NA)
  for (j in 1:length(nm)){
    dates.list[[i]]$date[j]<-unlist(strsplit(nm[j],split="/"))[[5]]
  }
}
  
# head(dates.list[[1]])

##export the trees, alignments, datefiles
for (i in 1:b){
  write.tree(trees[[i]], paste(tree.out.base, i, ".tre",sep=""))
  write.FASTA(align[[i]],paste(align.out.base, i, ".fasta",sep=""))
  write.table(dates.list[[i]], paste(date.out.base, i, ".txt",sep=""), row.names=F, sep="\t", quote = F, col.names=F)
}
