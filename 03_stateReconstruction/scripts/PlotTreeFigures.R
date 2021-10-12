#!/usr/bin/env Rscript

# Make timetree figure for the highest likelihood tree of the bootstrap set

#usage:
# $Rscript PlotTreeFigure.R 

#### Inputs and outputs ####
# args<-commandArgs(trailingOnly=TRUE)

#setup libs
library(ape)
library(ggplot2)
library(ggtree)
library(RColorBrewer)
library(lubridate)
library(stringr)
library(phangorn)
library(stringr)
library(phylobase)
library(dplyr)

## inputs and outputs
#color schemes
colors.in<-"../03_stateReconstruction/DF/globalcolors.tsv" 
lin.colors.in<-"../03_stateReconstruction/DF/lineagecolors.tsv"
lin.grp.colors.in<-"../03_stateReconstruction/DF/lineageGroupColors.tsv"

#lineage lookup table
lin.lookup.in<-"../03_stateReconstruction/DF/lineageGroups.csv"

#bootstrap setup
boot<-"b9"
tree.in<-"../02_trees/ft_root_res_time/res_rooted_subsamp_align_fake150_9.fasta.timetree.nex.binary.tre"
meta.in<-"../03_stateReconstruction/DF/meta.b9.csv"

#overall across bootstraps
#most recent sampling date
mrsd="2021-02-04"

#tmrca 
tmrca<-"2019-12-08" 
tmrca.dt<-format(as.Date(tmrca),"%d %b %Y")
tmrca.dt.up<-format(as.Date("2019-10-28"),"%d %b %Y")
tmrca.dt.low<-format(as.Date("2019-12-09"),"%d %b %Y")

#set up directory if doesn't exist
fold.out<-"results/"
if(!dir.exists(fold.out)){dir.create(fold.out)}

## setup province colors
#import a tsv of name and hex color 
globalPalette<-read.table(colors.in,sep="\t")

## make color scheme for Lineages
glob.colz<-row.names(globalPalette)
globalPalette.ch<-as.character(globalPalette$globalPalette)
names(globalPalette.ch)<-glob.colz

provPalette<-globalPalette.ch[str_which(names(globalPalette.ch),"Alberta|British Columbia|Manitoba|Maritimes|Ontario|Quebec") ] 
#Prince Edward Island|Yukon|Northwest Territories|Nunavut

colScale<-scale_colour_manual(name = "Canadian\nprovince",values = provPalette)
fillScale<-scale_fill_manual(name = "Canadian\nprovince",values = provPalette)

## Read in tree and data
#binarized newick tree (from multifurcating nexus)
tn<-read.tree(tree.in)
# tn<-ladderize(tn,right=F)
#import appropriate metadata
meta<-read.csv(meta.in,header=T)

#remove any tips not found in the tree
meta<-meta[which(meta$tip.label %in% tn$tip.label), ]
nrow(meta)==Ntip(tn)

#sort to match tree
meta<-meta[match(tn$tip.label,meta$tip.label),]

#convert martimes provinces
meta$division<-str_replace_all(meta$division, "Nova Scotia|New Brunswick|Newfoundland and Labrador|Newfoundland\nand Labrador","Maritimes")

#only keep cols you want to use as features in tree and tip labels
t.df<-meta[,c("tip.label","country","division")]
colnames(t.df)<-c("tip.label","country","division")

#make a column for canada yes or no
t.df$canada<-0
t.df$province<-NA

for (i in 1:nrow(t.df)){
  if (t.df$country[i]=="Canada") {t.df$canada[i]<-1; t.df$province[i]<-t.df$division[i] }
}
# table(t.df$canada)
table(as.factor(t.df$province))

## Basic trees
##  baseplot}
#Plot the timetree, plain

#issue with new ggtree version
p1<-ggtree(tn,mrsd=mrsd, as.Date=TRUE, aes(na.rm=T,color="state"),lwd=0.05,color="grey70")+
  theme_tree2() +#add scale bar on x axis
  theme(axis.text=element_text(size=18))
# p1
# ggsave("timeTreePlots/0211_timetree.png",width=9, height=10, units = "in")

## build up detail for time tree
#make dataframe for single point in annotation
point<-data.frame(x=as.Date("2019-12-16"), y=50) #y=7050
## timetree root date annotate
p2<-p1+
  geom_rootpoint(color="darkolivegreen", size=4, shape=18,alpha=0.9)+
  annotate("text", as.Date("2019-11-15"), y=40, hjust=0, size=3, label=paste("tMRCA\n",tmrca.dt,"\n(",tmrca.dt.up," -\n",tmrca.dt.low,")", sep=""),fontface="bold")+ #y=5700
  geom_point(data=point, aes(x=x,y=y),color="darkolivegreen", size=4, shape=18,alpha=0.9)

#merge metadata with tree 
p3<- p2 %<+% t.df [,c("tip.label","canada","province")]

## color provinces
colScale<-scale_colour_manual(name = "Canadian\nprovince",values = provPalette, labels=names(provPalette), breaks=names(provPalette))

p5<-p3+
  geom_tippoint(aes(color=as.factor(province),size=as.factor(canada)), shape=16, na.rm = TRUE,alpha=0.9)+
  colScale+
  scale_size_manual(values = c(-1,1),guide=FALSE)+
  scale_y_continuous(expand=c(0.01,0)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y",
               limits=as.Date(c("2019-11-10","2021-02-11")), expand=c(0.01,0))+
  theme(axis.text.x=element_text(angle=40,hjust = 1,size=10),
        axis.title.x = element_text(size=12),
        plot.margin=unit(c(1.5,0,0,0),"cm"),
        legend.title = element_text(face="bold"),
        legend.text=element_text(size=10), 
        legend.key.height=unit(.2, "in"),
        legend.key.width=unit(.1, "in"), 
        legend.position=c(0.5,0.99), 
        legend.direction = "horizontal", 
        legend.justification = "bottom")+
  guides(color=guide_legend(override.aes = list(size=4),ncol=6,title.position="top", title="Province"))+
  coord_cartesian(clip="off")
p5
ggsave(paste(fold.out,boot,"_timetree.png",sep=""),width=8, height=10, units = "in")
# 
# 
# # Lineage Group
# ## Setup for lineage group
# 
# #import a tsv of name and hex color 
# linGrpPalette<-read.table(lin.grp.colors.in,sep="\t")
# 
# ## make color scheme for Lineage groups
# linGrp.colz<-row.names(linGrpPalette)
# linGrpPalette.ch<-as.character(linGrpPalette$mygrpz.col)
# names(linGrpPalette.ch)<-linGrp.colz
# 
# # collinGrpScale<-scale_colour_manual(name = "Canadian\nprovince",values = linGrpPalette.ch)
# filllinGrpScale<-scale_fill_manual(name = "Lineage\nGroup",values = linGrpPalette.ch,na.translate = F)
# collinGrpScale<-scale_color_manual(name = "Lineage\nGroup",values = linGrpPalette.ch,na.translate = F)
# 
# #make lineage group col in meta (not in this version)
# lookup.lin<-read.csv(lin.lookup.in,header = T) #lineage group, alias (longform), lineagre
# meta$Lineage.grp<-NA
# for (i in 1:nrow(meta)){
#   match<-which(lookup.lin$lineage==meta$Lineage[i])
#   if(length(match)>0){
#       meta$Lineage.grp[i]<-lookup.lin$lineagegroup[which(lookup.lin$lineage==meta$Lineage[i])]
#   }
# }
# length(which(is.na(meta$Lineage.grp)))
# meta[which(is.na(meta$Lineage.grp)),]



## Add in the lineage GROUP as an annotation on the right side
# t.df2<-meta[,c("tip.label","country","division","Lineage.grp")]
# colnames(t.df2)<-c("tip.label","country","division","Lineage Group")
# 
# #make a column for canada yes or no
# t.df2$canada<-0
# t.df2$province<-NA
# 
# for (i in 1:nrow(t.df2)){
#   if (t.df2$country[i]=="Canada") {t.df2$canada[i]<-1; t.df2$province[i]<-t.df2$division[i] }
# }
# 
# #add lineage annotation layer
# t.df3<-t.df2
# rownames(t.df3)<-t.df3$tip.label
# t.df3<-t.df3[,'Lineage Group',drop=FALSE]
# colnames(t.df3)<-"Lineage\nGroup"
# 
# #lineage group layer
# plin<-gheatmap(p5, t.df3, offset = 0.01, color=NULL, width=0.03,
#          colnames_position="top",colnames_offset_y=1410,family = "bold") +
#   scale_fill_manual(values=linGrpPalette.ch, guide=FALSE,na.translate = F)
# 
# #tree layer
# plin2<- plin + 
#   scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", limits=as.Date(c("2019-11-10","2021-03-07")), expand=c(0.01,0))+
#   theme(axis.text.x=element_text(angle=40,hjust = 1,size=10),
#         legend.text=element_text(size=10), 
#         legend.key.height=unit(.2, "in"),
#         legend.key.width=unit(.1, "in"), 
#         axis.title.x = element_text(size=12))+
#   guides(color=guide_legend(override.aes = list(size=4),ncol=6,title.position="top", title="Province"))+
#   theme(plot.margin=unit(c(0.6,0,0,0),"cm"),
#         legend.position=c(0.5,0.96), 
#         legend.direction = "horizontal", 
#         legend.justification = "bottom")+
#   coord_cartesian(clip="off")
# 
# plin2
# ggsave(paste(fold.out,boot,"_timetree_linGrps.png",sep=""),width=8,height=10,units="in")

# Color branches by tip state
# t.p2<-t.p %<+% t.df [,c("tip.label","state")]
# ggtree(t,aes(color="state"))+
#     geom_tippoint(aes(color=state),size=2, shape=17)

# d = data.frame(node=1:59, color=sample(c('red', 'blue', 'green'), 59, replace=T))
# ggtree(tr) %<+% d + aes(color=I(color))

# set.seed(123)
# tr = rtree(30)
# g1 = as(tr, 'phylo4')
# d = data.frame(color=sample(c('red', 'blue', 'green'), 30, replace=T))
# rownames(d) = tr$tip.label
# g2 = phylo4d(g1, d)
# 
# rNodeData <- data.frame(randomTrait = rnorm(nNodes(g1)),
#                         color = sample(c('purple', 'yellow', 'black'), nNodes(g1), replace=T),
#                         row.names = nodeId(g1, "internal"))
# nodeData(g2) <- rNodeData
# ggtree(g2, aes(color=I(color)))

# g1 = as(tn, 'phylo4')
# collz<-as.data.frame(globalPalette.ch)
# collz$state<-rownames(collz)
# d = left_join(meta[,c("tip.label","state")], collz,by="state")
# rownames(d) = tn$tip.label
# g2 = phylo4d(g1, d)
# rNodeData <- data.frame(randomTrait = rnorm(nNodes(g1)),
#                         color = "grey70",
#                         row.names = nodeId(g1, "internal"))
# nodeData(g2) <- rNodeData
# tre1<-ggtree(g2, aes(color=I(color)))
# ggsave(tre1,"coloredtree.png")
# tre1
