#!/usr/bin/env Rscript
# Analyze output from ancestral state reconstruction

#usage:
# $Rscript AnalyzeStates.R 

# analyze output of the ancestral reconstruction 
# summarize the mean estimates across boots
# Look at origins, intros, descendants over time 
# sublineage, singleton modeling

# note that the fake example leads to some issues in the script b/c of low n (zero for some locations), some boots have no singletons, etc
# most parts work, sections that don't work with fake example are mostly commented out.

## setup libs
library(ape)
library(stringr)
library(dplyr)
library(tidyr)
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(ggtree)
library(ggimage)
library(phytools)
library(phangorn)
library(forcats)
library(lubridate)
library(gridExtra)
library(cowplot)
library(ggstance)
library(ggalluvial)
library(ggmosaic)
library(MASS)

#### In n outs ####
## setup bootstraps and inputs, outputs
BOOTS<-1:10
b<-length(BOOTS)

#input clean metadata for all and pangolin updated lineages (common to all boots)
clean.meta.in<-"../00_cleanData/cleaned/clean_fake_meta.csv"
pango.in<-"../00_cleanData/cleaned/pangolin_results.csv"

#input trees
trees.in.f<-"../02_trees/ft_root_res_time/"
trees.in<-list.files(trees.in.f, pattern="fasta.timetree.nex.binary.tre$", full.names = T)

#input DFs
df.in<-"DF"
can.anc.in<-list.files(df.in,"can.anc.b", full.names = T) %>% mixedsort()
sum.in<-list.files(df.in, "sublin.b", full.names = T) %>% mixedsort()
meta.in<-list.files(df.in, "InfDates", full.names = T) %>% mixedsort()
can.states.in<-list.files(df.in, "tip.anc", full.names = T) %>% mixedsort()

#make output directory if it doesn't exist
if(!dir.exists("results/")){dir.create("results/")}

#output of modified sublineage df with first Canadian date, etc
sum.out<-paste("DF/firstCanDate_sum.b",BOOTS,".csv",sep="")

#### read in color schemes for Lineages and Global ####
#import a tsv of name and hex color for LINEAGES 
#should be able to just use the tsv from the first bootstrap... we'll see
mycolz<-read.table("DF/lineagecolors.tsv",sep="\t")
lin.colz<-row.names(mycolz)
mycolz.ch<-as.character(mycolz$mycolz)
names(mycolz.ch)<-lin.colz
colScale<-scale_colour_manual(name = "Lineage",values = mycolz.ch,na.value="grey60")
fillScale<-scale_fill_manual(name = "Lineage",values = mycolz.ch,na.value="grey60")

#read in lineage groups
lookup.lin<-read.csv("DF/lineageGroups.csv",header = T) #lineage group, alias (longform), lineagre
lin.group.col<-read.table("DF/lineageGroupColors.tsv",sep="\t") #4 colors
lin.grp.colz<-row.names(lin.group.col)
lin.grp.ch<-as.character(lin.group.col[,1])
names(lin.grp.ch)<-lin.grp.colz
LinGrpColScale<-scale_colour_manual(name = "Lineage",values = lin.grp.ch,na.value="grey60")
LinGrpFillScale<-scale_fill_manual(name = "Lineage",values = lin.grp.ch,na.value="grey60")

### GEO
#import a tsv of name and hex color for LOCATIONS
globalPalette<-read.table("DF/globalcolors.tsv",sep="\t")
## make color scheme
glob.colz<-row.names(globalPalette)
globalPalette.ch<-as.character(globalPalette$globalPalette)
names(globalPalette.ch)<-glob.colz
GlobColScale<-scale_colour_manual(name = "Location",values = globalPalette.ch,na.value="grey60")
GlobFillScale<-scale_fill_manual(name = "Location",values = globalPalette.ch,na.value="grey60")

#read in lookup table to reduce cateogires to apply this color scheme (if haven't run MakeGlobalColors_boots.Rmd in this project already)
lookup.geo<-read.csv("DF/lookup.geo.csv",header = T)

## Figure publication themes
pubTheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background=element_rect("grey95"), axis.line = element_line(colour = "black"),
            legend.key.size = unit(0.5,"line"),
            text=element_text(size=10,face="bold"),
            legend.text=element_text(size=6))

pubThemeDate<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background=element_rect("grey95"), axis.line = element_line(colour = "black"),
                    text=element_text(size=10,face="bold"),
                    legend.key.size = unit(0.4,"line"),
                    legend.text=element_text(size=6),
                    axis.text.x=element_text(angle = 45,hjust=1,size=rel(1)))

# summary functions
#function to spit out mean and 95% CI (for small n) for a vector containing 10 observations from the bootstraps
mean.95ci<-function(vector){ 
  m<-mean(vector)
  sd<-sd(vector)
  up<-m+(sd/sqrt(10)*2.26)
  low<-m-(sd/sqrt(10)*2.26)
  print(paste(round(m,digits=2), " (",round(low,digits=2)," - ",round(up,digits=2),")",sep=""))
}

#1 digit
mean.95ci.1<-function(vector){ 
  m<-mean(vector)
  sd<-sd(vector)
  up<-m+(sd/sqrt(10)*2.26)
  low<-m-(sd/sqrt(10)*2.26)
  print(paste(round(m,digits=1), " (",round(low,digits=1)," - ",round(up,digits=1),")",sep=""))
}

#0 digit
mean.95ci.0<-function(vector){ 
  m<-mean(vector)
  sd<-sd(vector)
  up<-m+(sd/sqrt(10)*2.26)
  low<-m-(sd/sqrt(10)*2.26)
  print(paste(round(m,digits=0), " (",round(low,digits=0)," - ",round(up,digits=0),")",sep=""))
}


#X digit
mean.95ci.X<-function(vector, digits){ 
  m<-mean(vector)
  sd<-sd(vector)
  up<-m+(sd/sqrt(10)*2.26)
  low<-m-(sd/sqrt(10)*2.26)
  print(paste(round(m,digits=digits), " (",round(low,digits=digits)," - ",round(up,digits=digits),")",sep=""))
}


#scientific notation, 2 digits
mean.95ci.SC<-function(vector){ 
  m<-mean(vector)
  sd<-sd(vector)
  up<-m+(sd/sqrt(10)*2.26)
  low<-m-(sd/sqrt(10)*2.26)
  print(paste(formatC(m,format = "e", digits = 2), " (",formatC(low,format = "e", digits = 2)," - ",formatC(up,format = "e", digits = 2),")",sep=""))
}

#as above, but with input for n also
mean.95ci.X.n<-function(vector, digits, n){ 
  m<-mean(vector)
  sd<-sd(vector)
  t<-qt(0.025,(n-1),lower.tail=F)
  up<-m+(sd/sqrt(n)*t)
  low<-m-(sd/sqrt(n)*t)
  print(paste(round(m,digits=digits), " (",round(low,digits=digits)," - ",round(up,digits=digits),")",sep=""))
}

##### Summarize results for each bootstrap ####

#setup lists for dfs of interest
# can.anc.boots<-replicate(n=b,vector())
meta.boots<-replicate(n=b,vector())
sum.boots<-replicate(n=b,vector())
tree.boots<-replicate(n=b,vector())

#mid script lists for rolling means
Parent.Location.sum.boots<-replicate(n=b,vector())
Parent.Location.sum.boots.full<-replicate(n=b,vector())
Node.Location.sum.boots<-replicate(n=b,vector())
Node.Location.sum.boots.full<-replicate(n=b,vector())

#df of descendants
sublin.long<-replicate(n=b,vector())

# Loop the boots to summarize objects
for (k in BOOTS){
  #### Read in objects from ancestralreconstruction.rmd ####
  # can.anc.boots[[k]]<-read.csv(can.anc.in[k])
  sum.boots[[k]]<-read.csv(sum.in[k])
  meta.boots[[k]]<-read.csv(meta.in[k])
  tree.boots[[k]]<-read.tree(trees.in[k])
  
  #Rename provinces so Canada not in label
  sum.boots[[k]]$Node.Location<-sapply(sum.boots[[k]]$Node.Location, function (x) str_replace_all(x, "Canada_", ""))
  
  #make sure no "." instead of space 
  sum.boots[[k]]$Node.Location<-str_replace_all(sum.boots[[k]]$Node.Location, "\\.", " ")
  sum.boots[[k]]$Parent.Location<-str_replace_all(sum.boots[[k]]$Parent.Location, "\\.", " ")
  sum.boots[[k]]$Descendant.Location<-str_replace_all(sum.boots[[k]]$Descendant.Location, "\\.", " ")
  
  ## APPLY the lookup table to change names
  ## GROUP MARITIMES TOGETHER and small contributors (as per lookup table made in GlobalColors_boots.Rmd)
  for (i in 1:nrow(lookup.geo)){
    if(lookup.geo$new.loc[i] != lookup.geo$og.loc[i]){ #if they don't match, go replace instances
      #find matches in node loc, and replace with new loc
      pr.ma<-which(sum.boots[[k]]$Node.Location==lookup.geo$og.loc[i])
      if (length(pr.ma)>0) {sum.boots[[k]]$Node.Location[pr.ma]<-lookup.geo$new.loc[i]; next}
      par.ma<-which(sum.boots[[k]]$Parent.Location==lookup.geo$og.loc[i])
      if (length(par.ma)>0) {sum.boots[[k]]$Parent.Location[par.ma]<-lookup.geo$new.loc[i]; next}
    }
  }

  # make an object of provs
  provs<-c("Ontario","Quebec","Maritimes","Alberta","British Columbia","Manitoba")

  ## ADD A NEW COLUMN for lineage group using lookup table generated in LineageColors.Rmd (in the _timetree not _boots project)
    sum.boots[[k]]$Lineage.grp<-NA
  for (i in 1:nrow(sum.boots[[k]])){
    sum.boots[[k]]$Lineage.grp[i]<-lookup.lin$lineagegroup[which(lookup.lin$lineage==sum.boots[[k]]$Lineage[i])]
  }
  # head(sum.boots[[k]][,c('Lineage','Lineage.grp')],n=25)

  ## Make list of sublineages, where each one is a df in order to map descendants
  n.sl<-nrow(sum.boots[[k]])
  desc.df.list<-replicate(n=n.sl,vector())
  for (i in 1:n.sl){
    #name it by defining internal node
    names(desc.df.list)[[i]]<-sum.boots[[k]]$Node[i] 
    
    #extract accession ID (to lookup date in next chunk)
    accessions<-unlist(strsplit(sum.boots[[k]]$Descendant.AccessionID[i], ", "))
  
    #add each desc as a row in the df
    desc.df.list[[i]]<-data.frame(Lineage=sum.boots[[k]]$Lineage[i],
                                  Node=sum.boots[[k]]$Node[i], #use node to match back later
                                  Descendant.AccessionIDs=accessions,
                                  State=NA,
                                  Date=NA) 
  } #
  
  ## Look up each accession number for its sampling location and LSD-inferred date
  for (i in 1:n.sl){
    # i=10
    metaRows<-match (desc.df.list[[i]]$Descendant.AccessionIDs, meta.boots[[k]]$GISAID_ID) # in order
    #extract location and dates
    desc.df.list[[i]]$State<-meta.boots[[k]]$state [metaRows]
    desc.df.list[[i]]$Date<-as.Date(meta.boots[[k]]$date.lsd.full[metaRows])
  }

  ## Add the LSD dates onto the sumb1 object
  ## for each introduction, find the earliest sample collection date and where sampled
  sum.boots[[k]]$FirstCanDateLSD<-as.Date(NA)
  sum.boots[[k]]$FirstCanLocLSD<-NA
  sum.boots[[k]]$FirstCanGISAID<-NA 
  #repeat for lasts
  sum.boots[[k]]$LastCanDateLSD<-as.Date(NA)
  sum.boots[[k]]$LastCanLocLSD<-NA
  sum.boots[[k]]$LastCanGISAID<-NA 

    
  for (i in 1:nrow(sum.boots[[k]])){
    #list of desc is ordered same as sum.boot
    #subset to the descendants in Canada
    mtchCan<-desc.df.list[[i]] [str_detect(desc.df.list[[i]]$State, "Canada") ,]
    
    #order the subsetted dataframe by date
    mtchCan<-mtchCan[order(mtchCan$Date),]
    
    sum.boots[[k]]$FirstCanDateLSD[i]<-as.Date(first(mtchCan$Date))
    sum.boots[[k]]$FirstCanLocLSD[i]<-first(mtchCan$State)
    sum.boots[[k]]$FirstCanGISAID[i]<-first(mtchCan$Descendant.AccessionIDs)
    sum.boots[[k]]$LastCanDateLSD[i]<-as.Date(last(mtchCan$Date))
    sum.boots[[k]]$LastCanLocLSD[i]<-last(mtchCan$State)
    sum.boots[[k]]$LastCanGISAID[i]<-last(mtchCan$Descendant.AccessionIDs)
  }

  ## Name the sublineages as Lineage.canX, where X is the order of first sample detection for a given lineage
  #sort by first date
  sum.boots[[k]]<-with(sum.boots[[k]],sum.boots[[k]][(order(sum.boots[[k]]$FirstCanDateLSD)),])
  
  #make a new column 
  sum.boots[[k]]$Sublineage<-NA 
  sublins<-c()
  for (i in 1:nrow(sum.boots[[k]])){
    subl<-paste(sum.boots[[k]]$Lineage[i],".can",sep="")
    ss<-0
    subl2<-str_replace_all(subl,"\\.","\\\\.") #weird R shit
    ss<-length(str_which(sublins,subl2)) #how many sublins already designated?
    sum.boots[[k]]$Sublineage[i]<-paste(subl,(ss+1),sep="")
    sublins<-c(sublins,subl)
  }

  ## ANOTHER sublinaege name that incorps the first gisaid id
  sum.boots[[k]]$Sublineage2<-NA 
  for (i in 1:nrow(sum.boots[[k]])){
    subl<-paste(sum.boots[[k]]$Lineage[i],".can.",sum.boots[[k]]$FirstCanGISAID[i],sep="")
    sum.boots[[k]]$Sublineage2[i]<-subl
  }
  #check unique
  length(sum.boots[[k]]$Sublineage2)==unique(length(sum.boots[[k]]$Sublineage2))

  ## Add sublineage back into the desc.df.list object
  #use the node column to matchy
  for (i in 1:length(desc.df.list)){
    match<-which(sum.boots[[k]]$Node== desc.df.list[[i]]$Node[1])
    desc.df.list[[i]]$Sublineage<-sum.boots[[k]]$Sublineage[match]
    desc.df.list[[i]]$Sublineage2<-sum.boots[[k]]$Sublineage2[match]
  }
  #check
  # all(subs %in% sum.boots[[k]]$Sublineage)
  # all(sum.boots[[k]]$Sublineage %in% subs )
    
  #rbind them all into one df to refer to later
  sublin.long[[k]]<-rbind(desc.df.list[[1]],desc.df.list[[2]])
  for (i in 3:length(desc.df.list)){
    sublin.long[[k]]<-rbind(sublin.long[[k]], desc.df.list[[i]])
  }
  
  #### Pull tmrca ####
  ## Find MRCA node in time tree linking the descendent accession IDs in each sublineage (identified in divergence tree)
  tmrca.df<-data.frame(Node=sum.boots[[k]]$Node, mrca_timetree=NA, #mrca node in subs/site and time tree
                       tmrca=NA, tmrca.upper=NA, tmrca.lower=NA, #decimal format, rel. to root date
                       tmrca.dt=NA, tmrca.upper.dt=NA, tmrca.lower.dt=NA) #date format
  
  #from the LSD2 log file:
  #tMRCA 2019-12-09 [2019-10-27; 2019-12-10]
  root.date<-as.Date("2019-12-09") #these vary by one day
  #note this is from actual analysis in paper, not fake data date
  root.date.dec<-decimal_date(root.date)
  
  for (i in 1:nrow(sum.boots[[k]])){
    #find tip labels matching any ID
    tip.match<-str_which(tree.boots[[k]]$tip.label, str_replace_all(sum.boots[[k]]$Descendant.AccessionID[i],", ","\\|"))
    # length(tip.match)==sum.boots[[k]]$Number.Descendants[i]
    
    #find the MRCA of these tip matches
    mrca.i<-getMRCA(tree.boots[[k]],tip.match) #change this to t.p?
    
    if(is.null(mrca.i)){next}
    tmrca.df$mrca_timetree[i]<-mrca.i #add to df
  
    #check all the descendants are in the sublineage
    mrca.desc<-unlist(Descendants(tree.boots[[k]],node=mrca.i,"tips"))
    
    #pull the node height for the mrca
    tmrca.i<-nodeheight(tree.boots[[k]],mrca.i)
    tmrca.df$tmrca[i]<-tmrca.i
    
    #convert to calendar date
    tmrca.df$tmrca.dt[i]<-format(date_decimal(root.date.dec+tmrca.i),"%Y-%m-%d")
  }  
  
  ## connect TMRCA back to sum.boots[[k]] object and EXPORT
  # all(sum.boots[[k]]$Node %in% tmrca.df$Node)
  sum.boots[[k]]<-left_join(sum.boots[[k]],tmrca.df,by="Node")
  
  ## calculate the detection lag (FirstCanDate - tMRCA)
  ## TO DOOO also incorporate the tMRCA uncertainty here
  sum.boots[[k]]<-sum.boots[[k]] %>% mutate(detection.lag = as.Date(FirstCanDateLSD) -as.Date(tmrca.dt))
  sum.boots[[k]]$detection.lag<-as.numeric(sum.boots[[k]]$detection.lag)
  
  #these are issues related to incomplete dates...and lsd assigned their median and high estimates as 2021 despite the prior being 2020-12... so I think it is fine to manually change these to a month earlier for this analysis, but it points to the need for a better way to infer dates
  # sum.boots[[k]][sum.boots[[k]]$detection.lag<0,]
  sum.boots[[k]]$tmrca.dt<-as.Date(sum.boots[[k]]$tmrca.dt)
  # sum.boots[[k]]$tmrca.dt[sum.boots[[k]]$detection.lag<0]<-sum.boots[[k]]$tmrca.dt[sum.boots[[k]]$detection.lag<0] %m-% months(1)
  #recalculate detection lag 
  sum.boots[[k]]<-sum.boots[[k]] %>% mutate(detection.lag = as.Date(FirstCanDateLSD) -as.Date(tmrca.dt))
  sum.boots[[k]]$detection.lag<-as.numeric(sum.boots[[k]]$detection.lag)
  
  #write this for future use (in TMRCA, now included below, but just in case)
  write.csv(sum.boots[[k]], sum.out[k])
  
######
} #end of k loop across bootstraps
#####


# to skip chunk above, run

# for (k in 1:b){
#   sum.boots[[k]]<-read.csv(sum.out[k],header=T)
# }


#overall number of sublins across bootstraps
sublin.tot<-c()
for (k in 1:b){
  sublin.tot<-c(sublin.tot,nrow(sum.boots[[k]]))
}
mean.95ci(sublin.tot)


# Figure 2 plots: Sankey, mosaic, sizes
## summarize the total number and percent of intros by 1) par. loc, 2) prov of intro, maybe also 3) lineage - need to run this chunk for below 

#prep list item
sum.Par.Prov.l<-replicate(b,vector())

#go through each list item/subsample
## tabulate instances of location pairs
for (k in 1:b){
  sum.Par.Prov.l[[k]]<-sum.boots[[k]] %>% dplyr::group_by (Parent.Location, Node.Location) %>% 
    dplyr::summarize (.groups="rowwise", n.Par= n()) %>% 
    as.data.frame() %>% 
    mutate(perc.Par=(n.Par/sum(n.Par))*100)

  ## add a column for subsample
  sum.Par.Prov.l[[k]]$SampleSet<-paste("Sample",BOOTS[k])

}

## rbind the summaries

sum.Par.Prov<-rbind(sum.Par.Prov.l[[1]],sum.Par.Prov.l[[2]])
for (k in 3:b){
  sum.Par.Prov<-rbind(sum.Par.Prov,sum.Par.Prov.l[[k]])
}

#summarize the mean and range for each of 1), 2) and 3)
sum.Par.Prov.summary<-sum.Par.Prov %>% dplyr::group_by(Parent.Location, Node.Location) %>%
  dplyr::summarize(.groups="rowwise",
                   mean.n=round(mean(n.Par)),
                   sd.n=round(sd(n.Par,na.rm=T),digits=2),
                   mean.perc=round(mean(perc.Par),digits=2),
                   sd.perc=round(sd(perc.Par,na.rm=T),digits=2)) %>%
  as.data.frame()

tot.Par<-sum.Par.Prov.summary %>% dplyr::group_by(Parent.Location) %>%
  dplyr::summarise(.groups = "rowwise", totalImports=sum(mean.n)) %>%
  as.data.frame()

## Alluvial plot for Figure 2
#make a "subject column"
sum.Par.Prov.summary$subject<-1:nrow(sum.Par.Prov.summary)
# sum.Par.Prov.summary<-sum.Par.Prov.summary[-which(sum.Par.Prov.summary$mean.n<1),]

#order the geos
ord.count<-sum.Par.Prov.summary%>% dplyr::group_by(Parent.Location) %>%
  dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
ord.count<-ord.count[rev(order(ord.count$n)),'Parent.Location']
sum.Par.Prov.summary$Parent.Location<-factor(sum.Par.Prov.summary$Parent.Location,levels=ord.count)

ord.prov<-sum.Par.Prov.summary%>% dplyr::group_by(Node.Location) %>%
  dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
ord.prov<-ord.prov[rev(order(ord.prov$n)),'Node.Location']
sum.Par.Prov.summary$Node.Location<-factor(sum.Par.Prov.summary$Node.Location,levels=ord.prov)


#make it long
sum.Par.Prov.summary.long<-sum.Par.Prov.summary %>% pivot_longer(1:2, names_to = "geo.type", values_to = "geo")
#order the geo types
sum.Par.Prov.summary.long$geo.type<-factor(sum.Par.Prov.summary.long$geo.type,levels=c("Parent.Location","Node.Location"),labels=c("Global origin","Canadian province"))

# sum.Par.Prov.summary.long$geo[which(!sum.Par.Prov.summary.long$geo %in% names(globalPalette.ch))]


## Alluvial plot

P1<- ggplot(sum.Par.Prov.summary.long,
       aes(x = geo.type, stratum = geo, alluvium = subject,
           y = mean.n,
           fill = geo, label = geo)) +
  scale_x_discrete(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_flow(alpha = .6,width=0.45) +
  geom_stratum(alpha = .8,width=0.45) +
  geom_text(stat = "stratum", size = 3.4,min.y=0,fontface="bold") + #note to change min.y if many cats and intros
  pubTheme+
  theme(legend.position = "none", axis.line = element_blank(), text =element_text(size=14),
        axis.ticks.x = element_blank(),axis.text.x = element_text(hjust=c(0.5,0.6)))+
  labs(x=NULL,y="# sublineages")+
  GlobFillScale
P1
ggsave("results/sankeyPlot.Parent.Node.png",height=10,width=5,units="in")


## Repeat but stratified by lineage group for flows - DONT NEED TO RUN THIS

# #summarize the total number and percent of intros by 1) par. loc, 2) prov of intro, maybe also 3) lineage
# #prep list item
# sum.Par.Prov.L.l<-replicate(b,vector())
# 
# #go through each list item/subsample
# ## tabulate instances of location pairs
# for (k in 1:b){
#   sum.Par.Prov.L.l[[k]]<-sum.boots[[k]] %>% dplyr::group_by (Parent.Location, Node.Location, Lineage.grp) %>% dplyr::summarize (.groups="rowwise", n.Par= n()) %>% as.data.frame() %>% mutate(perc.Par=(n.Par/sum(n.Par))*100)
# 
#   ## add a column for subsample
#   sum.Par.Prov.L.l[[k]]$SampleSet<-paste("Sample",BOOTS[k])
# 
# }
# 
# sum.Par.Prov.L<-rbind(sum.Par.Prov.L.l[[1]],sum.Par.Prov.L.l[[2]])
# for (k in 3:b){
#   sum.Par.Prov.L<-rbind(sum.Par.Prov.L,sum.Par.Prov.L.l[[k]])
# }
# 
# 
# #summarize the mean and range for each of 1), 2) and 3)
# sum.Par.Prov.L.summary<-sum.Par.Prov.L %>% dplyr::group_by(Parent.Location, Node.Location, Lineage.grp) %>%
#   dplyr::summarize(.groups="rowwise",
#                    mean.n=round(mean(n.Par)),
#                    sd.n=round(sd(n.Par,na.rm=T),digits=2),
#                    mean.perc=round(mean(perc.Par),digits=2),
#                    sd.perc=round(sd(perc.Par,na.rm=T),digits=2)) %>%
#   as.data.frame()
# 
# ### Alluvial plot for Figure 2
# 
# #make a "subject column"
# sum.Par.Prov.L.summary$subject<-1:nrow(sum.Par.Prov.L.summary)
# # sum.Par.Prov.L.summary<-sum.Par.Prov.L.summary[-which(sum.Par.Prov.L.summary$mean.n<1),]
# 
# #order the geos by frequency
# ord.count<-sum.Par.Prov.L.summary%>% dplyr::group_by(Parent.Location) %>%
#   dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
# ord.count<-ord.count[rev(order(ord.count$n)),'Parent.Location']
# sum.Par.Prov.L.summary$Parent.Location<-factor(sum.Par.Prov.L.summary$Parent.Location,levels=ord.count)
# 
# ord.prov<-sum.Par.Prov.L.summary%>% dplyr::group_by(Node.Location) %>%
#   dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
# ord.prov<-ord.prov[rev(order(ord.prov$n)),'Node.Location']
# sum.Par.Prov.L.summary$Node.Location<-factor(sum.Par.Prov.L.summary$Node.Location,levels=ord.prov)
# 
# #make it long
# sum.Par.Prov.L.summary.long<-sum.Par.Prov.L.summary %>% pivot_longer(1:2, names_to = "geo.type", values_to = "geo")
# #order the geo types
# sum.Par.Prov.L.summary.long$geo.type<-factor(sum.Par.Prov.L.summary.long$geo.type,levels=c("Parent.Location","Node.Location"),labels=c("Origin location","Province of introduction"))
# 
# #plot it
# ggplot(sum.Par.Prov.L.summary.long,
#        aes(x = geo.type, stratum = geo, alluvium = subject,
#            y = mean.n, label = geo)) +
#   scale_x_discrete(expand = c(0.01,0.01)) +
#   scale_y_continuous(expand = c(0,0)) +
#   geom_flow(aes(fill = Lineage.grp),alpha = .8,width=0.4) +
#   scale_fill_manual(name = "Lineage.grp",values = lin.grp.ch,na.value="grey60")+
#   # new_scale("fill")+ #sneaky
#   geom_stratum(alpha = .4,width=0.4) +
#   scale_color_manual(name = "geo",values = globalPalette.ch,na.value="grey60")+
#   geom_text(stat = "stratum", size = 4,min.y=4) +
#   pubTheme+
#   theme(legend.position = "none",text =element_text(size=14),
#         axis.ticks.x = element_blank(),axis.text.x = element_text(hjust=c(0.55,0.65)))+
#   labs(x=NULL,y="# Sublineages seeded")
# #ggsave("results/sankeyPlot.Parent.Node.Lin.png",height=10,width=5,units="in")


## Mosaic plot of lineage by origin location

#similar to above, need to generate an object summarizing the mean of each insatance (origin, lineage) across boots
#then use this frequency table to make a mosaic plot, colored by lineage or lineage group

#prep list item
sum.Par.L.l<-replicate(b,vector())

#go through each list item/subsample
## tabulate instances of location pairs
for (k in 1:b){
  sum.Par.L.l[[k]]<-sum.boots[[k]] %>% dplyr::group_by (Parent.Location, Lineage.grp) %>% dplyr::summarize (.groups="rowwise", n.Par= n()) %>% as.data.frame() %>% mutate(perc.Par=(n.Par/sum(n.Par))*100)
  
  ## add a column for subsample
  sum.Par.L.l[[k]]$SampleSet<-paste("Sample",BOOTS[k])

}

sum.Par.L<-rbind(sum.Par.L.l[[1]],sum.Par.L.l[[2]])
for (k in 3:b){
  sum.Par.L<-rbind(sum.Par.L,sum.Par.L.l[[k]])
}


#summarize the mean and range for each of 1), 2) and 3)
sum.Par.L.summary<-sum.Par.L %>% dplyr::group_by(Parent.Location,  Lineage.grp) %>%
  dplyr::summarize(.groups="rowwise",
                   mean.n=round(mean(n.Par)),
                   sd.n=round(sd(n.Par,na.rm=T),digits=2),
                   mean.perc=round(mean(perc.Par),digits=2),
                   sd.perc=round(sd(perc.Par,na.rm=T),digits=2)) %>%
  as.data.frame()

#go from having frequencies to having indiviudal occurences

ord.count<-sum.Par.Prov.L.summary%>% dplyr::group_by(Parent.Location) %>%
  dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
ord.count<-ord.count[,'Parent.Location']

#order as above
sum.Par.L.summary$Parent.Location<-factor(sum.Par.L.summary$Parent.Location,levels=rev(ord.count))
sum.Par.L.summary$Lineage.grp<-factor(sum.Par.L.summary$Lineage,levels=(row.names(lin.group.col)))

P2<-sum.Par.L.summary %>% 
  ggplot()+
  geom_mosaic(aes(weight=mean.n,x=product(Lineage.grp,Parent.Location),
                  fill=Lineage.grp), offset = .002,alpha=0.9)+
  pubTheme+
  theme(axis.text.y = element_text(size=12),
        # axis.text.y = element_text(size=rev(c(12,rep(9,times=6),rep(8,times=4),6,6,6,4,4))),
        axis.line=element_blank(),
        legend.position = "top", 
        axis.ticks=element_blank(),
        legend.text=element_text(size=rel(1)),
        legend.margin=margin(unit(c(c(0,0,-10,-5)), units="line")),
        legend.title = element_text(size=rel(1)),
        axis.text.x = element_blank(),
        panel.background = element_blank(), 
        axis.title.y=element_text(vjust=-10,size=rel(1.2)), 
        plot.margin=unit(c(0.1,0.5,0,0),"cm"))+
    guides(fill=guide_legend(title.position = "top", keywidth = 1.1,keyheight = 1.1))+
    scale_fill_manual(name = "Lineage group",values = lin.grp.ch,na.value="grey60")+
  labs(y="",x="Global origin")+
  scale_x_productlist(expand=c(0.02,0))+
  scale_y_productlist(expand=c(0,0))+
  # annotate(geom="text",y=0.07,x=1.015,label="A*")+
  # annotate(geom="text",y=0.15,x=1.015,label="B*")+
  # annotate(geom="text",y=0.55,x=1.015,label="B.1*")+
  # annotate(geom="text",y=0.94,x=1.015,label="B.1.1*")+
  coord_flip()
P2
ggsave("results/MosaicPlot_LinGrpParentLoc2.png",height=6,width=4,units="in")

#### Summary of sublineage size ####
## okay, so here, need to summarize the mean sublineage size for each sublineage 
#need to be able to cross reference the sublineages across bootstraps
#want to convey that a small number of sublineages accounted for majority of sampled canadian cases

#first need a way of comparing sublineages... cross reference to see how common
sublin.df<-data.frame(Sublineage=NA,Sublineage2=NA,Lineage=NA,Number.Descendants=NA,
                      Lineage.grp=NA,SampleSet=NA)
for (k in 1:b){
  sum.boots[[k]]$SampleSet<-k
  #add sublineages to the dataframe
  sublin.df<-rbind(sublin.df,sum.boots[[k]][,c("Sublineage","Sublineage2","Lineage","Number.Descendants", "Lineage.grp","SampleSet")])
}
sublin.df<-sublin.df[-1,] #NA row from rbind



## Binned hist of sublineage size
# Plots of sublineage size
#not using this object outside here, so don't need list for dfs
#Bin sublineage size by frequency
sum.boots.binSize<-sum.boots[[2]] %>% group_by(Number.Descendants,Lineage.grp) %>% dplyr::summarize(.groups="rowwise",n=n())

#bin the larger sizes
# table(sum.boots.binSize$Number.Descendants)
#1:20 individual; 21-30; 31-40; 41-50;51-100;100-1000; 1000+
sum.boots.binSize$Number.Desc.Factor<-as.character("NA")
for (i in 1:nrow(sum.boots.binSize)){
  if (sum.boots.binSize$Number.Descendants[i]<10){
    sum.boots.binSize$Number.Desc.Factor[i]<-as.character(sum.boots.binSize$Number.Descendants[i]);next
  }
  if (sum.boots.binSize$Number.Descendants[i]>=10 & 
      sum.boots.binSize$Number.Descendants[i]<20){
    sum.boots.binSize$Number.Desc.Factor[i]<-as.character("10-19");next}
  if (sum.boots.binSize$Number.Descendants[i]>=20 & 
      sum.boots.binSize$Number.Descendants[i]<30){
    sum.boots.binSize$Number.Desc.Factor[i]<-as.character("20-29");next}
  if (sum.boots.binSize$Number.Descendants[i]>=30 & 
      sum.boots.binSize$Number.Descendants[i]<40){
    sum.boots.binSize$Number.Desc.Factor[i]<-as.character("30-39");next}
  if (sum.boots.binSize$Number.Descendants[i]>=40 & 
      sum.boots.binSize$Number.Descendants[i]<50){
    sum.boots.binSize$Number.Desc.Factor[i]<-as.character("40-49");next}
  if (sum.boots.binSize$Number.Descendants[i]>=50 & 
      sum.boots.binSize$Number.Descendants[i]<100){
    sum.boots.binSize$Number.Desc.Factor[i]<-as.character("50-99");next}
  if (sum.boots.binSize$Number.Descendants[i]>=100 & 
      sum.boots.binSize$Number.Descendants[i]<1000){
    sum.boots.binSize$Number.Desc.Factor[i]<-as.character("100-999");next}
  if (sum.boots.binSize$Number.Descendants[i]>=1000){
    sum.boots.binSize$Number.Desc.Factor[i]<-as.character("1000+");next}
}

# order numerically by Number.Desc.Factor then sublineage
sum.boots.binSize$Number.Desc.Factor<-factor(sum.boots.binSize$Number.Desc.Factor, levels=c(as.character(2:9),'10-19','20-29','30-39','40-49','50-99','100-999','1000+'))

sum.boots.binSize$Lineage.grp<-factor(sum.boots.binSize$Lineage.grp,levels=rev(row.names(lin.group.col)))

#make a summary df of n in each category for text and pos
sum.boots.binSize.summ<-sum.boots.binSize%>% dplyr::group_by(Number.Desc.Factor) %>%
  dplyr::summarise(.groups="rowwise",total.subs=sum(n)) %>% as.data.frame()

#PLOT IT 
plot.lim.y<-max(sum.boots.binSize.summ$total.subs)+5
P3<-ggplot(data=sum.boots.binSize)+
  geom_bar(aes(x=as.factor(Number.Desc.Factor),y=n, group=Lineage.grp,fill=Lineage.grp),
               position="stack",stat="identity",alpha=0.9)+
  pubTheme+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=0,hjust=1,size=rel(1.2)),
        # legend.position = c(0.7,0.8),
        legend.position = "none",
        legend.text = element_text(size=9),
        # axis.title.y=element_text(vjust=4),       
        axis.title.y=element_text(vjust=6,size=rel(1.2)), 
        plot.margin=unit(c(0,0.5,0,0),"cm"))+
  LinGrpFillScale+
  labs(x="Sublineage size",y="# sublineages")+
  guides(fill=guide_legend(ncol= 1,title="Lineage group",title.position="top",reverse = T, keywidth = 1.2,keyheight = 1.2))+
  geom_text(data=sum.boots.binSize.summ,aes(y=total.subs,x=Number.Desc.Factor,label=paste("n=",total.subs,sep=""),hjust=-0.1),size=2.7)+  
  scale_y_continuous(limits=c(0,plot.lim.y),expand=c(0.01,0))+
  coord_flip()

P3
ggsave("results/Sublineagesizes.bylin.Frequencies.png",width=3,height=5,units="in")



#### Export composite figure for manuscript ####
right_col <- plot_grid(P2,P3, labels = c('B', 'C'), label_size = 14,ncol=1, align="hv", axis="lr")
plot_grid(P1, right_col, labels = c('A', ''), label_size = 14, ncol = 2, align="h",axis="b",rel_widths = c(0.55,0.45))
ggsave(file="results/Composite_AlluvialMosaicSize.png",width=8.5,height=10,units = "in")




#### Figure: rolling sublineage importation rates over time ####
## Calculate rolling mean of importations by 1) origin, 2) destination using tmrca
for (k in 1:b){
  #make sure this a date
  sum.boots[[k]]$tmrca.dt<-as.Date(sum.boots[[k]]$tmrca.dt)
  
  #### Calculate a rolling 7-day mean for origins ####
  
  ## count the importations by origin location over time
  Parent.Location.sum.boots[[k]]<-sum.boots[[k]] %>%
    dplyr::select(Parent.Location, Lineage,tmrca.dt) %>%
    dplyr::group_by(tmrca.dt, Parent.Location) %>%
    dplyr::summarize(.groups="rowwise", total=n()) %>%
    dplyr::arrange(desc(Parent.Location)) %>% 
    dplyr::group_by(Parent.Location) 
  
  #need to add rows for missing dates
  alldays<-seq(ymd(first(sort(Parent.Location.sum.boots[[k]]$tmrca.dt))),
      ymd(last(sort(Parent.Location.sum.boots[[k]]$tmrca.dt))),
      by='1 day')
  
  #make a empty df in same structure as above then populate it
  nL<-length(unique(Parent.Location.sum.boots[[k]]$Parent.Location))
  nD<-length(alldays)
  Parent.Location.sum.boots.full[[k]]<-data.frame(tmrca.dt=rep(alldays,times=nL),
Parent.Location = sort(rep(unique(Parent.Location.sum.boots[[k]]$Parent.Location),times=nD)),
                                 total=0)
  # nrow(Parent.Location.sum.boots[[k]].empty)==nD*nL      
  
  #populate it
  for (i in 1:nrow(Parent.Location.sum.boots.full[[k]])){
    #look for a match
    match<-which(Parent.Location.sum.boots[[k]]$Parent.Location==Parent.Location.sum.boots.full[[k]]$Parent.Location[i] &
            Parent.Location.sum.boots[[k]]$tmrca.dt==Parent.Location.sum.boots.full[[k]]$tmrca.dt[i])
    if(length(match)==0) next #no match, no change
    #else, replace:
    Parent.Location.sum.boots.full[[k]]$total[i]<-Parent.Location.sum.boots[[k]]$total[match]
  }
  
  # sum(Parent.Location.sum.boots.full[[k]]$total[Parent.Location.sum.boots.full[[k]]$Parent.Location=="USA"])==sum(Parent.Location.sum.boots[[k]]$total[Parent.Location.sum.boots[[k]]$Parent.Location=="USA"])
  
  Parent.Location.sum.boots.full[[k]]<-Parent.Location.sum.boots.full[[k]] %>% 
    dplyr::mutate(intros_mean7d = zoo::rollmean(total, k = 7, fill = NA),
                  intros_mean14d = zoo::rollmean(total, k = 14, fill = NA),
                  intros_median7d = zoo::rollmedian(total, k = 7, fill = NA),
                  intros_sum7d = zoo::rollsum(total, k = 7, fill = NA,align="right")) %>% #right align to sum all prev
    #rolling mean and median weekly importation rate
    dplyr::mutate(intros_meansum7d = zoo::rollmean(intros_sum7d, k = 7, fill = NA), 
                  intros_mediansum7d = zoo::rollmedian(intros_sum7d, k = 7, fill = NA),) %>%
    dplyr::ungroup()
   
  ## count the importations by node.location over time
  Node.Location.sum.boots[[k]]<-sum.boots[[k]] %>%
    dplyr::select(Node.Location, Lineage,tmrca.dt) %>%
    group_by(tmrca.dt, Node.Location) %>%
    dplyr::summarize(.groups="rowwise", total=n()) %>%
    dplyr::arrange(desc(Node.Location)) %>% 
    dplyr::group_by(Node.Location) 
  
  #need to add rows for missing dates
  alldays<-seq(ymd(first(sort(Node.Location.sum.boots[[k]]$tmrca.dt))),
      ymd(last(sort(Node.Location.sum.boots[[k]]$tmrca.dt))),
      by='1 day')
  
  #make a empty df in same structure as above then populate it
  nL<-length(unique(Node.Location.sum.boots[[k]]$Node.Location))
  nD<-length(alldays)
  Node.Location.sum.boots.full[[k]]<-data.frame(tmrca.dt=rep(alldays,times=nL),
                                 Node.Location=sort(rep(unique(Node.Location.sum.boots[[k]]$Node.Location),times=nD)),
                                 total=0)
  # nrow(Node.Location.sum.boots[[k]].empty)==nD*nL      
  
  #populate it
  for (i in 1:nrow(Node.Location.sum.boots.full[[k]])){
    #look for a match
    match<-which(Node.Location.sum.boots[[k]]$Node.Location==Node.Location.sum.boots.full[[k]]$Node.Location[i] &
            Node.Location.sum.boots[[k]]$tmrca.dt==Node.Location.sum.boots.full[[k]]$tmrca.dt[i])
    if(length(match)==0) next #no match, no change
    #else, replace:
    Node.Location.sum.boots.full[[k]]$total[i]<-Node.Location.sum.boots[[k]]$total[match]
  }
  
  # sum(Node.Location.sum.boots.full[[k]]$total[Node.Location.sum.boots.full[[k]]$Node.Location=="USA"])==sum(Node.Location.sum.boots[[k]]$total[Node.Location.sum.boots[[k]]$Node.Location=="USA"])
  
  Node.Location.sum.boots.full[[k]]<-Node.Location.sum.boots.full[[k]] %>% 
    dplyr::mutate(intros_mean7d = zoo::rollmean(total, k = 7, fill = NA),
                  intros_mean14d = zoo::rollmean(total, k = 14, fill = NA),
                  intros_median7d = zoo::rollmedian(total, k = 7, fill = NA),
                  intros_sum7d = zoo::rollsum(total, k = 7, fill = NA,align="right")) %>% #right align to sum all prev
    #rolling mean and median weekly importation rate
    dplyr::mutate(intros_meansum7d = zoo::rollmean(intros_sum7d, k = 7, fill = NA), 
                  intros_mediansum7d = zoo::rollmedian(intros_sum7d, k = 7, fill = NA),) %>%
    dplyr::ungroup()

}



## Summarize the rolling means overall

## BY ORIGINS

#go through each list item/subsample
for (k in 1:b){
  ## add a column for subsample
  Parent.Location.sum.boots.full[[k]]$SampleSet<-paste("Sample",BOOTS[k])
}

## rbind the summaries
sum.Par.Roll<-rbind(Parent.Location.sum.boots.full[[1]],Parent.Location.sum.boots.full[[2]])
for (k in 3:b){
  sum.Par.Roll<-rbind(sum.Par.Roll,Parent.Location.sum.boots.full[[k]])
}

#summarize the mean and confint
sum.Par.Roll.summary<-sum.Par.Roll %>% dplyr::group_by(tmrca.dt, Parent.Location) %>%
  dplyr::summarize(.groups="rowwise",
                   intros_meansum7d.mean=(mean(intros_meansum7d,na.rm=T)),
                   intros_meansum7d.sd=sd(intros_meansum7d,na.rm=T)) %>%
  as.data.frame()

#order the geos
sum.Par.Roll.summary$Parent.Location<-factor(sum.Par.Roll.summary$Parent.Location,levels=ord.count)

## BY PROVINCE

#go through each list item/subsample
for (k in 1:b){
  ## add a column for subsample
  Node.Location.sum.boots.full[[k]]$SampleSet<-paste("Sample",BOOTS[k])
}

## rbind the summaries
sum.Prov.Roll<-rbind(Node.Location.sum.boots.full[[1]],Node.Location.sum.boots.full[[2]])
for (k in 3:b){
  sum.Prov.Roll<-rbind(sum.Prov.Roll,Node.Location.sum.boots.full[[k]])
}

#summarize the mean and confint
sum.Prov.Roll.summary<-sum.Prov.Roll %>% dplyr::group_by(tmrca.dt, Node.Location) %>%
  dplyr::summarize(.groups="rowwise",
                   intros_meansum7d.mean=(mean(intros_meansum7d,na.rm=T)),
                   intros_meansum7d.sd=(sd(intros_meansum7d,na.rm=T))) %>%
  as.data.frame()

#order the geos
sum.Prov.Roll.summary$Node.Location<-factor(sum.Prov.Roll.summary$Node.Location,levels=ord.prov)

## Plot sublineages over time using rolling rates, by tmrca, by origin or dest

#dataframe for travel restrictions
trav.df<-data.frame(x1=as.Date("2020-03-21"),y1=-0.5, x2=as.Date("2020-03-16"),y2=-1.5)

#rolling importation rate, by origin
GlobFillScale2<-scale_fill_manual(name = "Location",values = globalPalette.ch[7:(length(globalPalette.ch)-1)],na.value="grey60")
plot.lim.y<-3 #for fake example only, modify
p1<-sum.Par.Roll.summary %>%
  ggplot(aes(x=tmrca.dt,y=intros_meansum7d.mean,group=Parent.Location,fill=Parent.Location))+
  geom_density(stat="identity", position="stack",lwd=0,alpha=0.9)+
  GlobFillScale2+
  pubThemeDate+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text=element_text(size=10,face="bold"),
        legend.position = "none")+ #change to none if separate legend
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits = c(as.Date("2020-01-01"),as.Date("2021-01-01")) ) +
  labs(x="Date of most recent common ancestor",y="# sublineages introduced per week", fill="Origin Location")+
  #travel restrictions
  geom_vline(xintercept=as.Date("2020-03-21"),color="grey20",linetype=2)+
  annotate("text",label="Maximum\nstringency",
           x=as.Date("2020-03-28"),y=53,size=4,hjust=0,vjust=1,color="grey20")+
  #minimal stringency
  geom_vline(xintercept=as.Date("2020-10-10"),color="grey20",linetype=2)+
  annotate("text",label="Reduced\nstringency",
           x=as.Date("2020-10-14"),y=53,size=4,hjust=0,vjust=1,color="grey20")+
  scale_y_continuous(expand=c(0,0),limits=c(-0.1,plot.lim.y),breaks=seq(0,plot.lim.y,1)) #modify breaks for other ex
# p1

#fake plots to take the legends
# rolling importation rate, by origin
p1.f<-sum.Par.Roll.summary %>%
  ggplot(aes(x=tmrca.dt,y=intros_meansum7d.mean,group=Parent.Location,fill=Parent.Location))+
  geom_density(stat="identity", position="stack",lwd=0,alpha=0.9)+
  GlobFillScale2+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=9.5,face="bold"),
        legend.position = c(0.1,1),
        legend.margin=margin(0,0,0,0),
        legend.background = element_blank(),
        legend.justification = c(0,1))+
  guides(fill = guide_legend(keywidth = 0.7,keyheight=0.7,title.position = "top",title="Global origin",legend.spacing=0,ncol=4))

p1.guide<-get_legend(p1.f)

#rolling importation rate, by province destination
plot.lim.y<-5
provFillScale<-scale_fill_manual(name = "Location",values = globalPalette.ch[1:6],na.value="grey60")

p2<-sum.Prov.Roll.summary %>%
  ggplot(aes(x=tmrca.dt,y=intros_meansum7d.mean,group=Node.Location,fill=Node.Location))+
  geom_density(stat="identity", position="stack",lwd=0,alpha=0.9)+
  provFillScale+
  pubThemeDate+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")+ #change to none if running guide sep
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits = c(as.Date("2020-01-01"),as.Date("2021-01-01")) ) +
  labs(x="Date of most recent common ancestor",y="# sublineages introduced per week")+
  #travel restrictions
  geom_vline(xintercept=as.Date("2020-03-21"),color="grey20",linetype=2)+
  annotate("text",label="Maximum\nstringency",
           x=as.Date("2020-03-28"),y=53,size=4,hjust=0,vjust=1,color="grey20")+
  #minimal stringency
  geom_vline(xintercept=as.Date("2020-10-10"),color="grey20",linetype=2)+
  annotate("text",label="Reduced\nstringency",
           x=as.Date("2020-10-14"),y=53,size=4,hjust=0,vjust=1,color="grey20")+
  scale_y_continuous(expand=c(0,0),limits=c(-0.1,plot.lim.y),breaks=seq(0,plot.lim.y,1)) #modify for more
p2

#fake plot to take legend
p2.f<-sum.Prov.Roll.summary %>%
  ggplot(aes(x=tmrca.dt,y=intros_meansum7d.mean,group=Node.Location,fill=Node.Location))+
  geom_density(stat="identity", position="stack",lwd=0,alpha=0.9)+
  provFillScale+
  theme(axis.text.x=element_text(angle=45,hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background=element_rect("grey95"),
        text=element_text(size=9.5,face="bold"),
        legend.position = c(0.1,1),
        legend.margin=margin(0,0,0,0),
        legend.background = element_blank(),
        legend.justification = c(0,1))+
    guides(fill = guide_legend(keywidth = 0.7,keyheight=0.7,title.position = "top",title="Province",legend.spacing=0,ncol=3))
p2.guide<-get_legend(p2.f)

#make a list of grobs
## add in the Canadian province representation plots
plot.row<-plot_grid(p1,p2,ncol=2,labels=c("A" ,"B"))
plot.leg<-plot_grid(p1.guide,p2.guide,ncol=2,align="h")
plot.all<-plot_grid(plot.leg,plot.row,nrow=2,rel_heights = c(0.25,1),align = "v")
ggsave(plot.all,file="results/OriginsDestinationsInclRates0verTime.png",width=7,height=4.5,units = "in")

## summarize Maximum rolling rates

#Maximum importation rate
head(na.omit(sum.Par.Roll.summary[rev(order(sum.Par.Roll.summary$intros_meansum7d.mean)),]))
29.64286+(2.26* (6.230366/sqrt(9)))
29.64286-(2.26* (6.230366/sqrt(9)))
#maximum of 29.6 (25.0-34.3) new sublineages per week from USA on 2020-03-22 

#max by prov
head(na.omit(sum.Prov.Roll.summary[rev(order(sum.Prov.Roll.summary$intros_meansum7d.mean)),]))
#588 2020-03-22        Quebec              30.57143            5.826673
30.57143+(2.26* ( 5.826673/sqrt(9)))
30.57143-(2.26* ( 5.826673/sqrt(9)))
#maximum of 30.6 (26.2-35.0) new sublineages per week into Quebec on 2020-03-22 

#max into Ontario
temp<-sum.Prov.Roll.summary[sum.Prov.Roll.summary$Node.Location=="Ontario",]
head(na.omit(temp[rev(order(temp$intros_meansum7d.mean)),]))
# 635 2020-03-30       Ontario              9.742857            2.698618
9.742857+(2.26* (2.698618/sqrt(9)))
9.742857-(2.26* (2.698618/sqrt(9)))
#maximum of 9.7 (7.7-11.8) new sublineages per week into Ontario on 2020-03-30

#max into Ontario in Fall
temp<-sum.Prov.Roll.summary[sum.Prov.Roll.summary$Node.Location=="Ontario" & sum.Prov.Roll.summary$tmrca.dt>as.Date("2020-08-30"),]
head(na.omit(temp[rev(order(temp$intros_meansum7d.mean)),]))
#2117 2020-12-02       Ontario              6.814286            1.863907
6.814286+(2.26* (1.863907/sqrt(9)))
6.814286-(2.26* (1.863907/sqrt(9)))
#In the Fall, Ontario reached 6.8 (5.4-8.2) new sublineages per week on 2020-12-02

#overall sum of new sublineages/week in Canada max
sum.Roll.total<-sum.Prov.Roll.summary %>% dplyr::group_by(tmrca.dt) %>%
  dplyr::summarize(.groups="rowwise",
                   total_meansum7d.mean=(sum(intros_meansum7d.mean,na.rm=T)),
                   total_meansum7d.sd=(sum(intros_meansum7d.sd,na.rm=T))) %>%
  as.data.frame()
head(na.omit(sum.Roll.total[rev(order(sum.Roll.total$total_meansum7d.mean)),]),n=20)
#98  2020-03-22             53.24286           12.52869
 53.24286 +(2.26* (12.52869/sqrt(9)))
 53.24286 -(2.26* (12.52869/sqrt(9)))
#The total sublineage imporation rate peaked on 22 March 2020 with 53.2 (43.8-62.7) new sublineages imported weekly 
 
#what was the rate two weeks after the max?
sum.Roll.total[sum.Roll.total$tmrca.dt==as.Date("2020-04-05"),]
#112 2020-04-05             15.68571           5.814042
15.68571 +(2.26* (5.814042/sqrt(9)))
15.68571 -(2.26* (5.814042/sqrt(9)))
53.24286/15.68571 #3.4 fold reduction in 2 weeks
(53.24286 +(2.26* (12.52869/sqrt(9))))/ (15.68571 +(2.26* (5.814042/sqrt(9)))) #3.1
(53.24286 -(2.26* (12.52869/sqrt(9))))/ (15.68571 -(2.26* (5.814042/sqrt(9)))) #3.9

#fold decrease in prov-spec rolling importation rates between day X and Y
sum.Prov.Roll.summary$Node.Location[which(sum.Prov.Roll.summary$tmrca.dt == as.Date("2020-03-21"))]
sum.Prov.Roll.summary$intros_meansum7d.mean[which(sum.Prov.Roll.summary$tmrca.dt == as.Date("2020-03-22"))]/sum.Prov.Roll.summary$intros_meansum7d.mean[which(sum.Prov.Roll.summary$tmrca.dt == as.Date("2020-04-05"))]

#for stdev, need to calculate from each bootstrap fold change separately 
all.folds<-data.frame(Node.Location=NA,fold=NA)
for (k in 1:b){
  ss<-paste("Sample ",k,sep="")
  df<-sum.Prov.Roll[which(sum.Prov.Roll$SampleSet==ss),]
  folds<-df %>% dplyr::group_by(Node.Location) %>% dplyr::summarise(.groups="rowwise",
      fold=intros_meansum7d[tmrca.dt == as.Date("2020-03-22")]/
           intros_meansum7d[tmrca.dt == as.Date("2020-04-05")])
  all.folds<-rbind(all.folds,folds)
}
all.folds %>% dplyr::group_by(Node.Location) %>% 
  dplyr::summarise(.groups="rowwise",mean=mean(fold),
                   sd=sd(fold)) %>%
  mutate( upper=mean + 2.26* (sd/sqrt(9)),
          lower=mean - 2.26* (sd/sqrt(9))) %>%
  dplyr::ungroup() %>% as.data.frame()

#Range of total importations from summer to fall
range(sum.Roll.total$total_meansum7d.mean[which(sum.Roll.total$tmrca.dt>as.Date("2020-05-01") & sum.Roll.total$tmrca.dt<as.Date("2020-10-01"))])
#0.8142857 5.5714286
range(sum.Roll.total$total_meansum7d.sd[which(sum.Roll.total$tmrca.dt>as.Date("2020-05-01") & sum.Roll.total$tmrca.dt<as.Date("2020-10-01"))])
0.8142857 -(2.26* (0.8142857/sqrt(9))) #lowest of low
5.5714286 + (2.26* (3.3801159/sqrt(9))) #highest of high

#overall high in early nov and dec
sum.Roll.total[which(sum.Roll.total$tmrca.dt>as.Date("2020-11-01") & sum.Roll.total$tmrca.dt<as.Date("2020-12-31")),]
# 323 2020-11-02            7.7714286          2.3377286
#354 2020-12-03            7.4428571          2.5734602
7.7714286 -(2.26* (2.3377286/sqrt(9))) 
7.7714286 +(2.26* (2.3377286/sqrt(9))) 
# 7.8 (6.0-9.5)

7.4428571 - (2.26* (2.5734602/sqrt(9))) #
7.4428571 + (2.26* (2.5734602/sqrt(9))) #
#7.4 (5.5 - 9.4)



## Sublineage longevity and recency analysis AND ACTIVE
for (k in 1:b){
    
  #calculate sublienage longevity (first to last canadian sample date)
  #change: had previously been going to last sample overall, now focusing on Ca
  sum.boots[[k]]<- sum.boots[[k]] %>% mutate(longevitySample=LastCanDateLSD-FirstCanDateLSD,
                             longevityTMRCA=LastCanDateLSD-tmrca.dt)
  
  #make another column for active (in Canada) or not (case in last 3 months) 
  sum.boots[[k]]$active<-"no"
  for (i in 1:nrow(sum.boots[[k]])){
    if(sum.boots[[k]]$LastCanDateLSD[i]>as.Date("2020-11-11")){sum.boots[[k]]$active[i]<-"yes"}
  }

}

head(sum.boots[[k]]$longevityTMRCA)
table(sum.boots[[k]]$active)

#distrib of longevity, active and inactive separately
k=4
P3<-ggplot(sum.boots[[k]])+
  geom_density(aes(x=longevityTMRCA,fill=active),alpha=0.5)+
  scale_fill_manual(values=c("deepskyblue3","coral2"), labels=c("Inactive","Active"))+
  pubTheme+
  labs(x="Sublineage lifespan (days)",y="Density")
  # annotate(geom="text",x=50,y=0.01,label="Inactive",color="deepskyblue3",hjust=0)+
  # annotate(geom="text",x=180,y=.004,label="Active",color="coral2",hjust=0)+
  # annotate(geom="text",x=50,y=0.009,label=paste("Median = ", median(sum.boots[[k]]$longevityTMRCA[sum.boots[[k]]$active=="no"],na.rm=T)," days",sep=""),hjust=0,color="deepskyblue3")+
  # annotate(geom="text",x=180,y=0.003,label=paste("Median = ", median(sum.boots[[k]]$longevityTMRCA[sum.boots[[k]]$active=="yes"],na.rm=T)," days",sep=""),hjust=0,color="coral2")+
  # theme(legend.position="none")
ggsave("results/SublineageLifespanDensityBYACTIVE.png",height=4,width=4, units="in")


#is it associated with tmrca (did they reduce over time?)
ggplot(sum.boots[[k]])+
  geom_point(aes(y=longevityTMRCA, x=tmrca.dt))

#looks like yes
#for only inactive sublins
# note that all inactive in fake dataset
# mod.long<-lm(as.numeric(longevityTMRCA) ~ tmrca.dt + active,data=sum.boots[[k]])
# summary(mod.long)
# newdata2 <- data.frame(
  # tmrca.dt = rep(seq(from = min(sum.boots[[k]]$tmrca.dt), to = max(sum.boots[[k]]$tmrca.dt), length.out = 100), 2),
  # active = factor(rep(1:2, each = 100), levels = 1:2, labels =
  # c("no","yes")))

# newdata2 <- cbind(newdata2, predict(mod.long, newdata2, type = "response", se.fit=TRUE))
# newdata2 <- within(newdata2, {
#   longevityTMRCA <- (fit)
#   LL <- (fit - 1.96 * se.fit)
#   UL <- (fit + 1.96 * se.fit)
# })

ggplot(sum.boots[[k]])+
  geom_point(aes(y=longevityTMRCA, x=tmrca.dt,color=active),alpha=0.6)+
  # geom_ribbon(data=newdata2,aes(x=tmrca.dt, ymin=(LL),ymax=(UL),fill=active),alpha = 0.25) +
  # geom_line(data=newdata2,aes(x=tmrca.dt, y=longevityTMRCA, color=active), size = 0.5,alpha=0.5) +
  scale_color_manual(values=c("deepskyblue3","coral2"))+
  scale_fill_manual(values=c("deepskyblue3","coral2"))+
  pubThemeDate+
  theme(legend.position="none")+
  labs(x="Date of most recent common ancestor",y="Sublineage lifespan (days)")+
  # annotate(geom="text",x=as.Date("2020-01-01"),y=160,label="Inactive",color="deepskyblue3",hjust=0)+
  # annotate(geom="text",x=as.Date("2020-05-01"),y=220,label="Active sublineages\nin past 3 months",color="coral2",hjust=0)+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits = c(as.Date("2020-01-01"),as.Date("2021-01-01")),expand=c(0.01,0))+
  scale_y_continuous(expand=c(0.01,0), limits=c(0,390))
ggsave("results/LongevityOverTimebyActive.png",height=4,width=4, units="in")

#repeat this with separate models (stratified instead of adjusted)
# #note, no active in fake data

# dat.act<-sum.boots[[k]] %>% filter( active=="yes")
# dat.inact<-sum.boots[[k]] %>% filter( active=="no")
# 
# #active model and predictions
# mod.long.act<-lm(as.numeric(longevityTMRCA) ~ tmrca.dt,data=dat.act)
# summary(mod.long.act)
# # plot(mod.long.act)
# newdata.act <- data.frame(
#   tmrca.dt = rep(seq(from = min(sum.boots[[k]]$tmrca.dt), to = max(sum.boots[[k]]$tmrca.dt), length.out = 100), 1),
#   active = factor(rep(1:1, each = 100), levels = 1:1, labels =
#   c("yes")))
# newdata.act <- cbind(newdata.act, predict(mod.long.act, newdata.act, type = "response", se.fit=TRUE))
# newdata.act <- within(newdata.act, {
#   longevityTMRCA <- (fit)
#   LL <- (fit - 1.96 * se.fit)
#   UL <- (fit + 1.96 * se.fit)
# })
# 
# #INactive model and predictions
# mod.long.inact<-lm(as.numeric(longevityTMRCA) ~ tmrca.dt,data=dat.inact)
# summary(mod.long.inact)
# # plot(mod.long.inact)
# newdata.inact <- data.frame(
#   tmrca.dt = rep(seq(from = min(sum.boots[[k]]$tmrca.dt), to = max(sum.boots[[k]]$tmrca.dt), length.out = 100), 1),
#   active = factor(rep(1:1, each = 100), levels = 1:1, labels =
#   c("no")))
# newdata.inact <- cbind(newdata.inact, predict(mod.long.inact, newdata.inact, type = "response", se.fit=TRUE))
# newdata.inact <- within(newdata.inact, {
#   longevityTMRCA <- (fit)
#   LL <- (fit - 1.96 * se.fit)
#   UL <- (fit + 1.96 * se.fit)
# })
# 
# #plot again
# P5<-ggplot(sum.boots[[k]])+
#   geom_point(aes(y=longevityTMRCA, x=tmrca.dt,color=active),alpha=0.6)+
#   scale_color_manual(values=c("deepskyblue3","coral2"))+
#   geom_ribbon(data=newdata.act,aes(x=tmrca.dt, ymin=(LL),ymax=(UL)),fill="coral2",alpha = 0.25) +
#   geom_line(data=newdata.act,aes(x=tmrca.dt, y=longevityTMRCA), color="coral2",size = 0.5,alpha=0.5) +
#   geom_ribbon(data=newdata.inact,aes(x=tmrca.dt, ymin=(LL),ymax=(UL)),fill="deepskyblue3",alpha = 0.25,) +
#   geom_line(data=newdata.inact,aes(x=tmrca.dt, y=longevityTMRCA), color="deepskyblue3",size = 0.5,alpha=0.5,) +
#   pubThemeDate+
#   theme(legend.position="none")+
#   labs(x="Date of most recent common ancestor",y="Sublineage lifespan (days)")+
#   annotate(geom="text",x=as.Date("2020-01-01"),y=160,label="Inactive",color="deepskyblue3",hjust=0)+
#   annotate(geom="text",x=as.Date("2020-03-20"),y=310,label="Active sublineages\nin past 3 months",color="coral2",hjust=0)+
#   scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits = c(as.Date("2020-01-01"),as.Date("2021-01-01")),expand=c(0.01,0))+
#   scale_y_continuous(expand=c(0.01,0), limits=c(-50,390))+
#     coord_cartesian(clip="on",ylim=c(0,390))
# P5

#ggsave("results/LifespanOverTimebyActive.png",height=4,width=4, units="in")

P5<-ggplot(sum.boots[[k]])+
  geom_point(aes(y=longevityTMRCA, x=tmrca.dt,color=active),alpha=0.6)+
  scale_color_manual(values=c("deepskyblue3","coral2"))+
  # geom_ribbon(data=newdata.act,aes(x=tmrca.dt, ymin=(LL),ymax=(UL)),fill="coral2",alpha = 0.25) +
  # geom_line(data=newdata.act,aes(x=tmrca.dt, y=longevityTMRCA), color="coral2",size = 0.5,alpha=0.5) +
  # geom_ribbon(data=newdata.inact,aes(x=tmrca.dt, ymin=(LL),ymax=(UL)),fill="deepskyblue3",alpha = 0.25,) +
  # geom_line(data=newdata.inact,aes(x=tmrca.dt, y=longevityTMRCA), color="deepskyblue3",size = 0.5,alpha=0.5,) +
  pubThemeDate+
  theme(legend.position="none")+
  labs(x="Date of most recent common ancestor",y="Sublineage lifespan (days)")+
  # annotate(geom="text",x=as.Date("2020-01-01"),y=160,label="Inactive",color="deepskyblue3",hjust=0)+
  # annotate(geom="text",x=as.Date("2020-03-20"),y=310,label="Active sublineages\nin past 3 months",color="coral2",hjust=0)+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits = c(as.Date("2020-01-01"),as.Date("2021-01-01")),expand=c(0.01,0))
  # coord_cartesian(clip="on",ylim=c(0,390))
P5


## Composite figure of transmission lifespan
plot_grid(P3,P5,ncol=2,labels=c("A" ,"B"))
ggsave("results/SublineageLifespanbyActive.png",height=4,width=8, units="in")

## Number of descendent vs tmrca, negative binomial
#### # descendant vs tmrca.dt ######
k=4 #the highest likelihood boot in the empirical analysis

#basic Y ~ X
mod.desc<-glm.nb((Number.Descendants) ~ as.Date(tmrca.dt), data=sum.boots[[k]])
summary(mod.desc)

#adjust for node location
mod.desc2<-glm.nb((Number.Descendants) ~ as.Date(tmrca.dt) + Node.Location, data=sum.boots[[k]])
summary(mod.desc2)

#compare to model 1
anova(mod.desc,mod.desc2,test="lrt")
BIC(mod.desc);BIC(mod.desc2)

#instead of node location, active vs inactive in 2021
mod.desc3<-glm.nb((Number.Descendants) ~ as.Date(tmrca.dt) + active, data=sum.boots[[k]])
summary(mod.desc3)

#compare to model 1
anova(mod.desc,mod.desc3,test="lrt")
BIC(mod.desc);BIC(mod.desc3)

#now add node location on top of active
mod.desc4<-glm.nb((Number.Descendants) ~ as.Date(tmrca.dt) + active + Node.Location, data=sum.boots[[k]])
summary(mod.desc4)

#compare to mod3 (nested)
anova(mod.desc3,mod.desc4,test="lrt")
BIC(mod.desc3);BIC(mod.desc4)
AIC(mod.desc3);AIC(mod.desc4)
# neither AIC or BIC support also adding node.location...

plot(mod.desc3) #residuals and QQ plot don't look great...
#try logging
# mod.desc3b<-glm.nb((Number.Descendants) ~ as.Date(tmrca.dt) + active, data=sum.boots[[k]])


# go with model 3
mod.desc<-mod.desc3

#summarize coeff, exponentiated because negbin
nb.summ<-as.data.frame(cbind(exp(mod.desc$coefficients),exp(confint(mod.desc))))
pvalz<-as.data.frame(summary(mod.desc)[12])[,4]
nb.summ$pvalue<-pvalz
colnames(nb.summ)<-c("Estimate","Lower 95% bound","Upper 95% bound","p-value")
nb.summ[,1:3]<-format(nb.summ[,1:3],scientific=T,digits=3)
nb.summ[,4]<-format(nb.summ[,4], scientific = TRUE,digits = 3)
rownames(nb.summ)<-str_replace_all(rownames(nb.summ),"Node.Location","")
# rownames(nb.summ)<-str_replace_all(rownames(nb.summ),"Newfoundland\nand Labrador","Newfoundland and Labrador")
rownames(nb.summ)[1:2]<-c("Intercept","Detection lag")
write.csv(nb.summ,"results/tmrca_vs_logNumberDescendants_byLocation_NEGBIN.csv")
nb.summ

#obtain the mean predicted number of events for values of detection lag across its entire range for each level of ACTIVE and plot 
# newdata2 <- data.frame(
#   tmrca.dt = rep(seq(from = min(sum.boots[[k]]$tmrca.dt), to = max(sum.boots[[k]]$tmrca.dt), length.out = 100), 1))
# 
# newdata2 <- cbind(newdata2, predict.(mod.desc, newdata2, type = "link", se.fit=TRUE))
# newdata2 <- within(newdata2, {
#   NumberDescendent <- exp(fit)
#   LL <- exp(fit - 1.96 * se.fit)
#   UL <- exp(fit + 1.96 * se.fit)
# })
sum.boots[[k]]$active<-factor(sum.boots[[k]]$active)
# newdata1 <- data.frame(tmrca.dt = mean(sum.boots[[k]]$tmrca.dt), active = factor(1:2, levels = 1:2, labels = levels(sum.boots[[k]]$active)))
# 
# newdata1$phat<-predict(mod.desc,newdata1, type="response")
# newdata3<-cbind(sum.boots[[k]],newdata2)

newdata2 <- data.frame(
  tmrca.dt = rep(seq(from = min(sum.boots[[k]]$tmrca.dt), to = max(sum.boots[[k]]$tmrca.dt), length.out = 100), 2),
  active = factor(rep(1:2, each = 100), levels = 1:2, labels =
  levels(sum.boots[[k]]$active)))

newdata2 <- cbind(newdata2, predict(mod.desc, newdata2, type = "link", se.fit=TRUE))
newdata2 <- within(newdata2, {
  Number.Descendants <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

#PUBLICATION: MAKE A PLOT WITH NEGATIVE BINOMIAL FITS FOR ACTIVE vs INACTIVE
p3<-ggplot(data=sum.boots[[k]],aes(x=tmrca.dt,y=log(Number.Descendants)))+
  geom_point(aes(color=active),alpha=0.8,color="deepskyblue3")+
  geom_ribbon(data=newdata2,aes(ymin=log(LL),ymax=log(UL),fill=active),fill="deepskyblue3",alpha = 0.15) +
  geom_line(data=newdata2,aes(color=active), color="deepskyblue3",size = 0.5,alpha=0.4) +
  coord_cartesian(ylim=c(0.3,8.2))+
  pubThemeDate+
  theme(legend.position="none")+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",expand=c(0.01,0))+
  scale_y_continuous(expand=c(0,0))+
  labs(x="Date of most recent common ancestor",y="log(Number sampled descendants)",color="Province\nof introduction")
  # scale_color_manual(values=c("deepskyblue3","coral2"))+
  # scale_fill_manual(values=c("deepskyblue3","coral2"))+
  # annotate(geom="text",x=as.Date("2020-03-28"),y=6.2,label="Active sublineages\nin past 3 months",color="coral2",
           # hjust=0,size=4)+
  # annotate(geom="text",x=as.Date("2020-04-26"),y=2.3,label="Inactive",color="deepskyblue3",
           # hjust=0,size=4)
p3
ggsave("results/tmrca_vs_logNumberDescendants_byProv_NEGBIN.png",height=4,width=4,units="in")

## Time since importation
sublin.long2<-replicate(b,vector())
k=4
#join select columns from sum.boots onto sublin.long by sublineage2
sublin.long2[[k]]<-left_join(sublin.long[[k]],sum.boots[[k]][,c('Sublineage2','Parent.Location','tmrca.dt')],
                             by="Sublineage2")

# what is the average time since importation for sampled cases over time?
# for each sampled descendant of a sublineage, how long since the tmrca?
sublin.long2[[k]]$sublinAge<-NA
sublin.long2[[k]]<- sublin.long2[[k]] %>% mutate(sublinAge=Date - tmrca.dt)

#only keep canadian descendants
# table(sublin.long2[[k]]$State)
sublin.long2[[k]]$State<-str_replace_all(sublin.long2[[k]]$State,"Canada_","")
keepers<-which(sublin.long2[[k]]$State %in% provs)
sublin.long.Can<-sublin.long2[[k]] [keepers,] 

#calculate this as a rolling mean over time
summAge<-sublin.long.Can %>% dplyr::group_by(Date) %>%
  dplyr::summarize(.groups="rowwise",
                   meanAge=mean(sublinAge,na.rm=T),
                   sdAge=as.numeric(sd(sublinAge,na.rm=T)),
                   n=n())
summAge[is.na(summAge)]<-0 #get rid of NAs
summAge2<- summAge %>% as.data.frame() %>%
  dplyr::summarize(.groups="keep",
                   Date=Date,
                   meanAge2=meanAge,
                   #note that qt(.025, n, lower.tail=F) looks up the critical value for a given n
                   upperAge=as.numeric(meanAge)+( (sdAge/sqrt(n-1))*qt(.025, n, lower.tail=F) ), 
                   lowerAge=as.numeric(meanAge)-( (sdAge/sqrt(n-1))*qt(.025, n, lower.tail=F) ) ) %>% 
  dplyr::ungroup() %>% as.data.frame()
head(summAge2,n=20)
#if we can't calculate a sd, just make the upper/lower the mean
for (i in 1:nrow(summAge2)){
  if (is.na(summAge2$upperAge[i])){
    summAge2$lowerAge[i]<-summAge2$meanAge2[i]
    summAge2$upperAge[i]<-summAge2$meanAge2[i]
  }
  #if negative, make it zero
  if (summAge2$lowerAge[i]<0){summAge2$lowerAge[i]<-0}
}
# summAge2[is.na(summAge2)] #no more NAs

summAge3<- summAge2 %>% as.data.frame() %>%
  dplyr::mutate(meanAge7d = zoo::rollmean(as.numeric(meanAge2), k = 7, fill=NA,align="center"),
                upperAge7d = zoo::rollmean(as.numeric(upperAge),k=7, fill=NA,align="center"),
                lowerAge7d = zoo::rollmean(as.numeric(lowerAge),k=7, fill=NA,align="center"),
                meanAge14d = zoo::rollmean(as.numeric(meanAge2), k = 14, fill=NA,align="center"),
                upperAge14d = zoo::rollmean(as.numeric(upperAge),k=14, fill=NA,align="center"),
                lowerAge14d = zoo::rollmean(as.numeric(lowerAge),k=14, fill=NA,align="center")) %>% 
  dplyr::ungroup() %>% as.data.frame()
head(summAge3)

#Add to the model
p4<-sublin.long2[[k]] %>% filter(State %in% provs) %>%
  ggplot()+
  geom_point(aes(x=Date,y=sublinAge,color=State),alpha=0.3,size=0.6)+
  geom_line(data=summAge3,aes(x=Date,y=meanAge14d),color="grey10",alpha=0.7,size=1)+
  geom_ribbon(data=summAge3,aes(x=Date,ymin=lowerAge14d,ymax=upperAge14d),fill="grey10",alpha=0.5)+
  # geom_line(data=pred.val,aes(x=Date,y=fit),color="grey10",alpha=0.8,size=1.5)+
  pubThemeDate+
  GlobColScale+
  theme(legend.position="none")+
  labs(x="Canadian descendant sampling date",y="Days since importation")+
  # annotate(geom="text",x=as.Date("2020-09-08"),y=210,label="Rolling mean",hjust=1,size=3)+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits = c(as.Date("2020-01-01"),as.Date("2021-01-01")) )+
  scale_y_continuous(limits=c(0,320),expand=c(0.01,0))
p4
# #ggsave("results/TimeSinceImportationVsSampleDate.png",height=4,width=5,units="in")

# do we see that the cases in 2021 were due to recent importations or old?



## assemble figure 3 part 1

plot.modz<-plot_grid(p3,p4, labels=c("C","D"))
plot_grid(plot.all, plot.modz,rel_heights = c(1.2,1),align="hv",nrow = 2)
ggsave("results/Fig3_SublineagesOvertime_mod.png",height=8,width=8,units="in")


# SUMMARY STATS FOR ALL BOOTSTRAPS
#tmrca of the b1.1.7 lineage
#this won't work for fake data, no b.1.1.7
# datez<-c()
# for (k in 1:b){
#   dt<-as.character(sum.boots[[k]][which(sum.boots[[k]]$Lineage=="B.1.1.7"),"tmrca.dt"])
#   datez<-c(datez,dt)
# }
# mean.95ci(as.Date(datez))
# 
# #origin of the b1.1.17 lineage
# parz<-c()
# for (k in 1:b){
#   par<-as.character(sum.boots[[k]][which(sum.boots[[k]]$Lineage=="B.1.1.7"),"Parent.Location"])
#   parz<-c(parz,par)
# }
# parz 
# #40% supported Turkey as the origin: 2,3,8,10
# 
# #investigate Turkey origin further... 
# #table of B.1.1.7 sequence locations with collection date before Jan 2020 for each boot
# for (k in 1:b){
#   print(rev(sort( table( meta.boots[[k]]$country [meta.boots[[k]]$Lineage=="B.1.1.7" & meta.boots[[k]]$date.lsd.full < as.Date("2020-12-31") ] )))[1:5])
# }
# #hmmm doesn't resolve this
# 
# #what are the descendant ids? investigate the tree
# sum.boots[[k]]$Descendant.AccessionID[which(sum.boots[[k]]$Lineage=="B.1.1.7")]
# #"EPI_ISL_751796, EPI_ISL_751797"

# how many sublineagse had >100 descendants?
hund<-c()
hund.perc<-c()
hund.des<-c()
for (k in 1:b){ 
  big<-which(sum.boots[[k]]$Number.Descendants>99)
  hund<-c(hund, length(big))
  hund.perc<-c(hund.perc, length(big)/nrow(sum.boots[[k]])*100)
  # What percent of sampled cases in Canada were downstream of these sublineages?
  can.des<-length(which(sublin.long[[k]]$State [ sublin.long[[k]]$Sublineage %in% sum.boots[[k]]$Sublineage[big] ] %in% provs))
  hund.des<-c(hund.des, can.des/nrow(meta.boots[[k]][meta.boots[[k]]$country=="Canada",])*100)
}
mean.95ci.1(hund)
mean.95ci.1(hund.perc)
mean.95ci.1(hund.des)

#how many total sampled descendants
samp.dec<-c()
for (k in 1:b){
  dec<-length(unique(sublin.long[[k]]$Descendant.AccessionIDs))
  samp.dec<-c(samp.dec,dec)
}
mean.95ci.1(samp.dec) #"10589.3 (9662.9-11515.7)"

#mean and conf.int of sampled descendants
samp.dec.mean<-c()
samp.dec.lower<-c()
samp.dec.upper<-c()

for (k in 1:b){
  dec<-sum.boots[[k]]$Number.Descendants
  n<-nrow(sum.boots[[k]])
  m<-mean(dec)
  sd<-sd(dec)
  t<-qt(0.025,(n-1),lower.tail=F)
  up<-m+(sd/sqrt(n)*t)
  low<-m-(sd/sqrt(n)*t)
  
  samp.dec.mean<-c(samp.dec.mean, m)
  samp.dec.lower<-c(samp.dec.lower, low)
  samp.dec.upper<-c(samp.dec.upper, up)
}
  
mean.95ci.X.n(samp.dec.mean, 1, 10) 
mean.95ci.X.n(samp.dec.lower, 1, 10) 
mean.95ci.X.n(samp.dec.upper, 1, 10) 
#mean 26.6 (13.0-40.1)

#how many unique locations (cross ref to meta)
samp.countr<-c()
for (k in 1:b){
  dec<-unique(sublin.long[[k]]$Descendant.AccessionIDs)
  countr<-length(unique(meta.boots[[k]]$country [which( meta.boots[[k]]$GISAID_ID %in% dec) ]))
  samp.countr<-c(samp.countr,countr)
}
mean.95ci.1(samp.countr) #"71.3 (68.5-74.1)"

#how many sampled desc in Canada only
samp.dec<-c()
for (k in 1:b){
  dec<-length(unique(sublin.long[[k]]$Descendant.AccessionIDs[sublin.long[[k]]$State %in% provs]))
  samp.dec<-c(samp.dec,dec)
}
mean.95ci.1(samp.dec) #""8281.3 (8242.9-8319.7)"

#how many importations before 25 January?
imports<-c()
for (k in 1:b){
  imports<-c(imports, length(which(sum.boots[[k]]$tmrca.dt<as.Date("2020-01-25"))))
}
mean.95ci.1(imports)

#form where?
#how many importations before 25 January?
imports.loc<-replicate(b,vector())
for (k in 1:b){
  imports.loc[[k]]<-table(sum.boots[[k]] [which(sum.boots[[k]]$tmrca.dt<as.Date("2020-01-25")),
                                          c('Parent.Location',"Lineage")])
}
imports.loc

#how many tmrca before first canadian covid case
first<-"2020-01-25"
early<-which(sum.boots[[k]]$tmrca.dt<as.Date("2020-01-25"))
length(early) #11
table(sum.boots[[k]][early,'Node.Location'])
table(sum.boots[[k]][early,'Parent.Location'])
sum.boots[[k]][which(sum.boots[[k]]$tmrca.dt<as.Date("2020-01-25")),]

#how many importations before 21 March?
imports<-c()
for (k in 1:b){
  imports<-c(imports, length(which(sum.boots[[k]]$tmrca.dt<as.Date("2020-03-21"))))
}
mean.95ci.1(imports)

#number of final canadian sequences after temporal outliers removed
can.count<-c()
for (k in 1:b){ 
  can.count<-c(can.count,nrow(cana.boots[[k]]))
}
mean.95ci(can.count)
range(can.count)

#total number of migration events into Canada 
for (k in 1:b){  
  print(nrow(sum.boots[[k]]))
}  

#what percent of sublineages came from the top sources
top.sources<-names(rev(sort(table(sum.boots[[1]]$Parent.Location))))
for (i in 1:length(top.sources)){
  source<-c()
  for (k in 1:b){
    perc<- signif(length(which(sum.boots[[k]]$Parent.Location==top.sources[i]))/
      nrow(sum.boots[[k]])*100, digits = 3)
    source<-c(source, perc)
  }
  print(top.sources[i])
  mean.95ci.1(source)
}

#top provinces of intro
top.provs<-names(rev(sort(table(sum.boots[[1]]$Node.Location))))
for (i in 1:length(top.provs)){
  prov<-c()
  for (k in 1:b){
    perc<- signif(length(which(sum.boots[[k]]$Node.Location==top.provs[i]))/
      nrow(sum.boots[[k]])*100, digits = 3)
    prov<-c(prov, perc)
  }
  print(top.provs[i])
  mean.95ci.1(prov)
}

#strict clock rate
str<-c(7.27E-04,8.53E-04,6.77E-04,7.77E-04,7.74E-04,7.52E-04,7.56E-04,7.50E-04,7.32E-04,7.10E-04)
mean.95ci.SC(str)

# # unique pango lineages
un.lins<-c()
for (k in 1:b){
    un.lins<-c(un.lins, length(unique(sum.boots[[k]]$Lineage)))
  }
mean.95ci.1(un.lins)

#top lineages (percentage)
top.lins<-names(rev(sort(table(sum.boots[[1]]$Lineage))))[1:15]
for (i in 1:length(top.lins)){
  lin<-c()
  for (k in 1:b){
    perc<- signif(length(which(sum.boots[[k]]$Lineage==top.lins[i]))/
      nrow(sum.boots[[k]])*100, digits = 3)
    lin<-c(lin, perc)
  }
  print(top.lins[i])
  mean.95ci.1(lin)
}

##top lineages (number)
top.lins<-names(rev(sort(table(sum.boots[[1]]$Lineage))))[1:100]
for (i in 1:length(top.lins)){
  lin<-c()
  for (k in 1:b){
    n <- length(which(sum.boots[[k]]$Lineage==top.lins[i]))
    lin<-c(lin, n)
  }
  print(top.lins[i])
  mean.95ci.1(lin)
}

#B.1.1.7
##top lineages (number)
lin<-c()
for (k in 1:b){
  n <- length(which(sum.boots[[k]]$Lineage=="B.1.1.7"))
  lin<-c(lin, n)
}
mean.95ci.1(lin)

#other "canadian lineages"
#B.1.1.7
##top lineages (number)
lin<-c()
for (k in 1:b){
  n <- length(which(sum.boots[[k]]$Lineage=="C.8"))
  lin<-c(lin, n)
}
mean.95ci.1(lin)

lin<-c()
for (k in 1:b){
  n <- length(which(sum.boots[[k]]$Lineage=="AM.2"))
  lin<-c(lin, n)
}
mean.95ci.1(lin)



# supplementary figure of descendant locations

#summarize the total descendants by location

#go through each list item/subsample
## tabulate instances of location pairs
sum.Desc.l<-replicate(b,vector())

for (k in 1:b){
  
  sublin.long[[k]]$State<-str_replace_all(sublin.long[[k]]$State, "Canada_","") 
  
  #convert using lookup table
   for (i in 1:nrow(lookup.geo)){
      if(lookup.geo$new.loc[i] != lookup.geo$og.loc[i]){ #if they don't match, go replace instances
        #find matches in node loc, and replace with new loc
        pr.ma<-which(sublin.long[[k]]$State==lookup.geo$og.loc[i])
        if (length(pr.ma)>0) {sublin.long[[k]]$State[pr.ma]<-lookup.geo$new.loc[i]; next}
      }
    }
    
  sum.Desc.l[[k]]<-sublin.long[[k]] %>% dplyr::group_by (State) %>% 
    dplyr::summarize (.groups="rowwise", n.Par= n()) %>% as.data.frame() 

  ## add a column for subsample
  sum.Desc.l[[k]]$SampleSet<-paste("Sample",BOOTS[k])

}

## rbind the summaries
sum.Desc<-rbind(sum.Desc.l[[1]],sum.Desc.l[[2]])
for (k in 3:b){
  sum.Desc<-rbind(sum.Desc,sum.Desc.l[[k]])
}

#summarize the mean and range for each 
sum.Desc.summary<-sum.Desc %>% dplyr::group_by(State) %>%
  dplyr::summarize(.groups="rowwise",
                   mean.n=round(mean(n.Par)),
                   sd.n=round(sd(n.Par,na.rm=T),digits=2)) %>%
  as.data.frame()

## plot it
sum.Desc.summary<-sum.Desc.summary[rev(order(sum.Desc.summary$mean.n)),]
sum.Desc.summary$State<-factor(sum.Desc.summary$State,
                               levels=sum.Desc.summary$State)
#big ones
sum.Desc.summary %>%
      ggplot(aes(x=State,y=mean.n,group=State,fill=State))+
        geom_bar(stat="identity")+
        theme_bw()+
        pubThemeDate+
        theme(legend.position = "none")+
        labs(x="Sublineage descendant location", y="# sampled descendants")+
        GlobFillScale
# ggsave("results/DescendantsBarPlot.png",height=4, width=6,units="in")


## Plot LSD-inferred dates: median and range
k=4
#note that incomp.df is generated from the InferredDatesCI script
#reorder by LSDlow
#for this build, use incomp.df.2 (which also has the 0115 inferred dates) from inferredDatesCI
incomp.df.2<-meta.boots[[k]][!is.na(meta.boots[[k]]$date.lsd),]

# plot a comparison of the median inferred date and LSD-inferred dates with CIs
meta.boots[[k]] %>% 
  dplyr::filter(tip.label %in% incomp.df.2$tip.label) %>%
  mutate(tip.label = fct_reorder(tip.label, desc(as.Date(date.lsd.low)))) %>%
  ggplot(aes(y=tip.label,group=tip.label))+
  #make line and point smaller for more sublins, but large here
  geom_linerangeh(aes(xmin=as.Date(date.lsd.low), xmax=as.Date(date.lsd.high)),linetype=1,color="blue",alpha=0.2,width=4)+
  geom_point(aes(x=as.Date(date.lsd)),color="blue",size=2, alpha=0.8)+
  # geom_point(aes(x=as.Date(date.median)),color="red",size=0.2, alpha=0.8)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_rect("grey95"), axis.line = element_line(colour = "black"),
  axis.text.x=element_text(angle=45,hjust = 1))+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y")+
  labs(x="Inferred sample collection date",y=paste("Samples with incomplete dates: n=",(nrow(incomp.df.2)),sep=""))

ggsave("results/InferredDatesLSDonly.png",height=8, width=6,unit="in")

# Detection lag modelling
## Detection lag over time
k=4
#overall distribution of values
px1<-ggplot(sum.boots[[k]],aes(x=as.numeric(detection.lag)))+
  geom_density(fill="grey20")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background=element_rect("grey95"),
        legend.key.size = unit(0.5,"line"),text=element_text(size=10,face="bold"),
        legend.text=element_text(size=8),axis.text.x=element_text(angle=45,hjust = 1))+
  labs(x="Detection lag (days)",y="Density")+
  theme(legend.position = "none")

#write a generic lm function
make.lm.eqn<-function(lm){
  slope<-coef(lm)[2]
  yint<-coef(lm)[1]
  pval<-anova(lm)$'Pr(>F)'[1]
  radj<-summary(lm)$adj.r.squared
  eqn<-paste("y = ",signif(slope,digits=3)," x + ",signif(yint,digits=3),"\n",
             "adj. R^2 = ",signif(radj,digits = 2),"\n",
             "p-value = ", signif(pval,digits = 2), sep="")
  return(eqn)
}

#make a model and add the lm fit to the plot
dl.mod<-lm((detection.lag) ~ as.Date(tmrca.dt), data=sum.boots[[k]])
dl.mod.2<-lm((detection.lag) ~ as.Date(tmrca.dt)+Node.Location, data=sum.boots[[k]])
anova(dl.mod,dl.mod.2)
summary(dl.mod.2)
#support for including node location
sum.boots[[k]]$Node.Location<-factor(sum.boots[[k]]$Node.Location)
newdata2 <- data.frame(
  tmrca.dt = rep(seq(from = min(sum.boots[[k]]$tmrca.dt,na.rm=T), 
                     to = max(sum.boots[[k]]$tmrca.dt,na.rm=T), length.out = 100), 6),
  Node.Location = factor(rep(1:6, each = 100), levels = 1:6, labels =
  levels(sum.boots[[k]]$Node.Location)))

newdata2 <- cbind(newdata2, predict(dl.mod.2, newdata2, type = "response", se.fit=TRUE))
newdata2 <- within(newdata2, {
  detection.lag<- (fit)
  LL <- (fit - 1.96 * se.fit)
  UL <- (fit + 1.96 * se.fit)
})

#Color by province and add linear fit overall
ggplot(sum.boots[[k]],aes(x=as.Date(tmrca.dt),y=detection.lag))+
  geom_point(aes(color=Node.Location))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background=element_rect("grey95"),
        legend.key.size = unit(0.5,"line"),text=element_text(size=10,face="bold"),
        legend.text=element_text(size=8),axis.text.x=element_text(angle=45,hjust = 1))+
  pubThemeDate+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y")+
  labs(x="Date of most recent common ancestor",y="Detection lag (days)",color="Province\nof introduction")+
  theme(legend.position = "none")+
  GlobColScale+
  GlobFillScale
  # annotate("text",label=dl.mod.text,x=as.Date(max.x),y=max.y,hjust=1,color="black",size=3,vjust=1)+
  # geom_line(data=newdata2,aes(y=detection.lag,color=Node.Location),lwd=0.4)+
  # geom_ribbon(data=newdata2,aes(ymin = UL, ymax = LL,fill=Node.Location),alpha = .5) 
ggsave("results/DetectionLag_SampleDateLSD-tMRCA_AdjbyProv.png",height=4,width=4,units="in")


#repeat with one fit
dl.mod

newdata2 <- data.frame(
  tmrca.dt = rep(seq(from = min(sum.boots[[k]]$tmrca.dt), to = max(sum.boots[[k]]$tmrca.dt), length.out = 100), 6))

newdata2 <- cbind(newdata2, predict(dl.mod, newdata2, type = "response", se.fit=TRUE))
newdata2 <- within(newdata2, {
  detection.lag<- (fit)
  LL <- (fit - 1.96 * se.fit)
  UL <- (fit + 1.96 * se.fit)
})

#Color by province and add linear fit overall
px3<-ggplot(sum.boots[[k]],aes(x=as.Date(tmrca.dt),y=detection.lag))+
  geom_point(aes(color=Node.Location))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background=element_rect("grey95"),
        legend.key.size = unit(0.5,"line"),text=element_text(size=10,face="bold"),
        legend.text=element_text(size=8),axis.text.x=element_text(angle=45,hjust = 1))+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y")+
  labs(x="Date of most recent common ancestor",y="Detection lag (days)",color="Province\nof introduction")+
  theme(legend.position = "none")+
  GlobColScale+
  geom_line(data=newdata2,aes(y=detection.lag),lwd=0.4,color="grey30")+
  geom_ribbon(data=newdata2,aes(ymin = UL, ymax = LL),fill="grey30",alpha = .5) 
px3
ggsave("results/DetectionLag_SampleDateLSD-tMRCA.png",height=4,width=4,units="in")

## plot distrib of detection lag by prov
dl.summ<-sum.boots[[k]] %>% group_by(Node.Location) %>% 
  dplyr::summarize(.groups="rowwise",
                   mean.dl=mean(detection.lag), 
                   median.dl=median(detection.lag), 
                   total=n())
# dl.summ
# are detection lags sign'ly diff b/w provinces
kruskal.test(sum.boots[[k]]$detection.lag, sum.boots[[k]]$Node.Location)
pairwise.wilcox.test(sum.boots[[k]]$detection.lag, sum.boots[[k]]$Node.Location,
                 p.adjust.method = "bonferroni")
library(dunn.test)
dunn.test::dunn.test(sum.boots[[k]]$detection.lag, sum.boots[[k]]$Node.Location, method="bonferroni")
# p-value = 0.002* b/w ontario and BC in real data
#none in fake data

#which ones are 
ord1<-dl.summ[rev(order(dl.summ$median.dl)),]$Node.Location
sum.boots[[k]]$Node.Location<-factor(sum.boots[[k]]$Node.Location,levels=ord1)
px2<-ggplot(sum.boots[[k]])+
  geom_boxplot(aes(x=Node.Location,group=Node.Location,y=detection.lag))+
  geom_point(aes(x=Node.Location,group=Node.Location,y=detection.lag,color=Node.Location),alpha=0.8)+
  GlobColScale+
  theme_bw()+
  pubThemeDate+
  labs(x=NULL,y="Detection lag (days)")+
  theme(legend.position="none")+
  geom_text(data=dl.summ,aes(label=paste("n=",total,sep=""),x=Node.Location,y=-5),size=3,vjust=1)
  # geom_segment(aes(x=2, xend=5, y=310, yend=310),size=0.1)
  # geom_text(aes(label="p = 0.0002",x=3.5,  y=314 ),size=2.5,fontface="plain",vjust=0)
px2
ggsave("results/DetectionLag_by_province_boxplotpoints.png",width=4,height=4,units="in")

## Plot number of descendants by detection lag, colored by prov, linear model
k=4
mod.9<-lm(log(Number.Descendants) ~ detection.lag,data=sum.boots[[k]])
# summary(mod.9) #SIGNIFICANT

#SETUP FOR PLOT
mod.9.text<-make.lm.eqn(mod.9)
#set positions
max.x<-max(sum.boots[[k]]$detection.lag)
max.y<-max(log(sum.boots[[k]]$Number.Descendants))
pred.mod.9<-predict(mod.9,interval="confidence")
dat.mod.9 <-cbind(sum.boots[[k]], pred.mod.9)

#### Plot number of descendants by detection lag, colored by prov, 
###with eqn for line and line
# ggplot(dat.mod.9,aes(x=detection.lag))+
#   geom_point(aes(y=log(Number.Descendants), color=Node.Location))+
#   theme_bw()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black"), panel.background=element_rect("grey95"),
#         legend.key.size = unit(0.5,"line"),text=element_text(size=10,face="bold"),
#         legend.text=element_text(size=8),axis.text.x=element_text(angle=45,hjust = 1))+
#   # scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y")+
#   labs(x="Detection lag (days)",y="log(Number Sampled Descendants)",color="Province\nof introduction")+
#   GlobColScale+
#   # annotate("text",label=mod.9.text,x=max.x,y=max.y,hjust=1,color="black",size=3,vjust=1)+
#   geom_line(aes(y=fit),color="darkblue",lwd=0.4)+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), fill="lightblue",alpha = .5) 

#ggsave("results/DetectionLag_vs_logNumberDescendants_byProv_LM.png",height=4,width=6,units="in")

#look at adjusting for tmrca.dt here
mod.9b<-lm(log(Number.Descendants) ~ detection.lag+as.Date(tmrca.dt),data=sum.boots[[k]])
summary(mod.9)
summary(mod.9b)
AIC(mod.9b);AIC(mod.9) #lower AIC for first model
BIC(mod.9b);BIC(mod.9) #lower BIC for first model
anova(mod.9,mod.9b,test="LRT") #not significanlty better


## detection lag vs number descendants, adjusted by prov (in supp publication)


mod.9.pois<-glm(Number.Descendants ~ detection.lag,data=sum.boots[[k]], family="poisson")
mean(sum.boots[[k]]$Number.Descendants)
sd(sum.boots[[k]]$Number.Descendants) #overdispersed, need negative binomial
summary(mod.9.pois) #SIGNIFICANT


mod.9.nb<-glm.nb((Number.Descendants) ~ detection.lag,data=sum.boots[[k]])
mod.9.nb.2<-glm.nb((Number.Descendants) ~ detection.lag + Node.Location,data=sum.boots[[k]])


summary(mod.9.nb) #SIGNIFICANT
summary(mod.9.nb.2) #SIGNIFICANT
anova(mod.9.nb, mod.9.nb.2)

#summarize coeff
nb.summ<-as.data.frame(cbind(exp(mod.9.nb.2$coefficients),exp(confint(mod.9.nb.2))))
pvalz<-as.data.frame(summary(mod.9.nb.2)[12])[,4]
nb.summ$pvalue<-pvalz
colnames(nb.summ)<-c("Estimate","Lower 95% bound","Upper 95% bound","p-value")
nb.summ[,1:3]<-round(nb.summ[,1:3],digits=3)
nb.summ[,4]<-round(nb.summ[,1:3],digits=3)
nb.summ[,4]<-format(nb.summ[,4], scientific = TRUE,digits = 3)
rownames(nb.summ)<-str_replace_all(rownames(nb.summ),"Node.Location","")
rownames(nb.summ)[1:2]<-c("Intercept","Detection lag")

#generate predictions for each province
sum.boots[[k]]$Node.Location<-factor(sum.boots[[k]]$Node.Location, levels=provs)
newdata1 <- data.frame(detection.lag = mean(sum.boots[[k]]$detection.lag), Node.Location = factor(1:6, levels = 1:6, 
    labels = levels(sum.boots[[k]]$Node.Location)))
newdata1$phat <- predict(mod.9.nb.2, newdata1, type = "response")

#obtain the mean predicted number of events for values of detection lag across its entire range for each level of node.location and plot 
newdata2 <- data.frame(
  detection.lag = rep(seq(from = min(sum.boots[[k]]$detection.lag), to = max(sum.boots[[k]]$detection.lag), length.out = 100), 6),
  Node.Location = factor(rep(1:6, each = 100), levels = 1:6, labels =
  levels(sum.boots[[k]]$Node.Location)))

newdata2 <- cbind(newdata2, predict(mod.9.nb.2, newdata2, type = "link", se.fit=TRUE))
newdata2 <- within(newdata2, {
  NumberDescendent <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

#PUBLICATION: MAKE A PLOT WITH NEGATIVE BINOMIAL FITS FOR EACH PROV
px4<-ggplot(dat.mod.9,aes(x=detection.lag))+
  geom_point(aes(y=log(Number.Descendants), color=Node.Location),alpha=0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), panel.background=element_rect("grey95"),
        legend.key.size = unit(0.5,"line"),text=element_text(size=10,face="bold"),
        legend.text=element_text(size=8),axis.text.x=element_text(angle=45,hjust = 1))+
  # scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y")+
  labs(x="Detection lag (days)",y="log(Number sampled descendants)",color="Province\nof introduction")+
  # geom_ribbon(data=newdata2, aes(detection.lag, log(NumberDescendent),ymin = log(LL), ymax = log(UL), fill = Node.Location), alpha = .1) +
  # geom_line(data=newdata2, aes(detection.lag, log(NumberDescendent), colour = Node.Location), size = 0.5,alpha=0.7) +
  GlobFillScale+
  GlobColScale+
  coord_cartesian(ylim=c(0,10))+
  theme(legend.position = "none")
px4
ggsave("results/DetectionLag_vs_logNumberDescendants_byProv_NEGBIN.png",height=4,width=4,units="in")

## Merge the plots of detection lag for supp fig
plot_grid(px1,px2,px3,px4,labels=c("A","B","C","D"))
ggsave("results/SuppFig_DetectionLagModels.png",width=8,height=8,units="in")


# Publication supplementary figure of sublineages' dates over time

for (k in 1:b){
  #find min and max dates, then leftjoin to select sumboots cols by sublineage
  sublin.long[[k]]$DateMin<-as.Date(NA)
  sublin.long[[k]]$DateMax<-as.Date(NA)
  
  for(i in 1:nrow(sum.boots[[k]])){
    IDz<-str_replace_all(sum.boots[[k]]$Descendant.AccessionID[i],", ","\\|")
    rowz<-str_which(sublin.long[[k]]$Descendant.AccessionIDs, IDz)
    datez<-sort(sublin.long[[k]]$Date[rowz])
    sublin.long[[k]]$DateMin[rowz]<-as.Date(first(datez))
    sublin.long[[k]]$DateMax[rowz]<-as.Date(last(datez))
  }
  
  #add on by sublineage
  add.on<-sum.boots[[k]][,c("Sublineage","Node.Location","tmrca.dt","FirstCanDateLSD")]
  add.on$Node.Location<-factor(add.on$Node.Location)
  sublin.long[[k]]<-left_join(sublin.long[[k]],add.on,by="Sublineage")
  
}

# head(sublin.long[[k]])

k=4
# plot these dates over time for each sublinage
sublin.long[[k]]$Sublineage <- reorder(sublin.long[[k]]$Sublineage, sublin.long[[k]]$DateMin)
provColScale<-scale_color_manual(name = "Location",values = globalPalette.ch[1:6],na.value="grey60")

ggplot(sublin.long[[k]],aes(y=Sublineage,color=Node.Location))+
  #make larger for this example
  geom_linerangeh(aes(xmin=DateMin,xmax=DateMax),lwd=0.5)+
  geom_linerangeh(aes(xmin=tmrca.dt,xmax=FirstCanDateLSD),linetype="dotted",lwd=0.5)+ #detection lag  
  geom_point(aes(x=Date),alpha=0.7,size=1)+
  geom_point(aes(x=FirstCanDateLSD),shape=22,size=1.5)+
  geom_point(aes(x=tmrca.dt),shape=18,size=2)+
  theme_bw()+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background=element_rect("grey95"),
          legend.key.size = unit(0.5,"line"),
          text=element_text(size=10,face="bold"),
          legend.text=element_text(size=8),
        legend.position="top")+
    scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y")+
  labs(x=NULL)+
  provColScale+
  scale_y_discrete(expand = c(0.01,0.01))+
  guides(color=guide_legend(title="Province\nof Introduction",title.position = "left"))
  # ylim(c(-1,nrow(sum.boots[[k]])+1))
ggsave("results/SamplesOverTimeForSublineages.png",height=5,width=7,units="in")

# Transmission source of all Canadian tips, not just intros: singeltons +

## Introductions over time by MONTH and province

#setup lists
can.states.boots<-replicate(b,vector())
cana.boots<-replicate(b,vector())
prov.summ.mo<-replicate(b,vector())
prov.summ.mo.long<-replicate(b,vector())
prov.summ.mo.long.2<-replicate(b,vector())

#color scheme from colorbrewer, n=4, qual
prop.colz<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c")

for (k in 1:b){
  
  ##Read in the Canadian tip state (and parent state) object 
  can.states.boots[[k]]<-read.csv(can.states.in[k])
  
  ## Calculate relative contribution of within-province, between province, USA, and other international   transmission sources
  #canadian metadata
  
  #make year-month a new column
  meta.boots[[k]]$month<-as.character(NA)
  for (i in 1:nrow(meta.boots[[k]])){
    meta.boots[[k]]$month[i]<-format(as.Date(meta.boots[[k]]$date.lsd.full[i]), "%Y-%m")
  }
  cana.boots[[k]]<-meta.boots[[k]] [str_which(meta.boots[[k]]$tip.label, "Canada"),]

  #make a new column in sumb1 to match (useful in another chunk)
  sum.boots[[k]]$month<-as.character(NA)
  for (i in 1:nrow(sum.boots[[k]])){
    sum.boots[[k]]$month[i]<-format(as.Date(sum.boots[[k]]$tmrca.dt[i]), "%Y-%m")
  }
    
  #merge maritime provinces
  cana.boots[[k]]$division<-str_replace_all(cana.boots[[k]]$division, "Nova Scotia|New Brunswick|Newfoundland and Labrador", "Maritimes")
  can.states.boots[[k]]$par.state<-str_replace_all(can.states.boots[[k]]$par.state, "Nova Scotia|New Brunswick|Newfoundland and Labrador", "Maritimes")
  
  #prime new df of canadian samples, grouped by province and month of sample
  prov.summ.mo[[k]]<-data.frame(table(cana.boots[[k]]$division,cana.boots[[k]]$month))
  colnames(prov.summ.mo[[k]])<-c("Province","month","N.Sample")
  prov.summ.mo[[k]]<-prov.summ.mo[[k]] %>% arrange(-N.Sample)

  #set up empty columns
  prov.summ.mo[[k]]$N.tips<-NA
  prov.summ.mo[[k]]$N.provincial<-NA
  prov.summ.mo[[k]]$N.domestic<-NA
  prov.summ.mo[[k]]$N.USA<-NA
  prov.summ.mo[[k]]$N.international<-NA
  prov.summ.mo[[k]]$N.totalinternational<-NA
  prov.summ.mo[[k]]$Prop.provincial<-NA
  prov.summ.mo[[k]]$Prop.domestic<-NA
  prov.summ.mo[[k]]$Prop.USA<-NA
  prov.summ.mo[[k]]$Prop.international<-NA
  prov.summ.mo[[k]]$Prop.totalinternational<-NA
  
  i=1
  #For each unique province and month
  for (i in 1:nrow(prov.summ.mo[[k]])) {  
    pr<-as.character(prov.summ.mo[[k]]$Province[i])
    seas<-as.character(prov.summ.mo[[k]]$month[i])
    
    #identify all can tip label in that prov and month
    pr.labs<-cana.boots[[k]]$tip.label[which(cana.boots[[k]]$division==pr & cana.boots[[k]]$month==seas)]
    
    #identify all can tip label in that prov
    pr.idz<-match(pr.labs,can.states.boots[[k]]$tip.label)
    pr.state<-can.states.boots[[k]]$par.state[pr.idz]
    
    #replace the "." with " "
    pr.state<-str_replace_all(pr.state,"\\."," ")
    
    prov.summ.mo[[k]]$N.tips[i]<-length(pr.state)
    #what number of them are the same province?
    prov.summ.mo[[k]]$N.provincial[i]<-length(str_which(str_replace_all(pr.state, "Canada_", "") ,pr))
    #a different province?
    prov.summ.mo[[k]]$N.domestic[i]<-length(str_which(pr.state,"Canada")) - prov.summ.mo[[k]]$N.provincial[i]
    #USA?
    prov.summ.mo[[k]]$N.USA[i]<-length(str_which(pr.state,"USA"))
    #Other international?
    prov.summ.mo[[k]]$N.international[i]<-length(pr.state) - length(str_which(pr.state,"Canada") ) - length(str_which(pr.state,"USA"))
    #total international?
    prov.summ.mo[[k]]$N.totalinternational[i]<-length(pr.state) - length(str_which(pr.state,"Canada") )
    
    #calculate a proportion for each of these
    prov.summ.mo[[k]]$Prop.provincial[i]<-prov.summ.mo[[k]]$N.provincial[i]/length(pr.state)
    prov.summ.mo[[k]]$Prop.domestic[i]<-prov.summ.mo[[k]]$N.domestic[i]/length(pr.state)
    prov.summ.mo[[k]]$Prop.USA[i]<-prov.summ.mo[[k]]$N.USA[i]/length(pr.state)
    prov.summ.mo[[k]]$Prop.international[i]<-prov.summ.mo[[k]]$N.international[i]/length(pr.state)
    prov.summ.mo[[k]]$Prop.totalinternational[i]<-prov.summ.mo[[k]]$N.totalinternational[i]/length(pr.state)
  
  }
  
  #Go from wide to long (gather proportions)
  prov.summ.mo.long[[k]]<-prov.summ.mo[[k]] %>% gather(Type,Proportion,Prop.provincial:Prop.totalinternational)

  #Plot the number of each type of transmission in a line plot
  prov.summ.mo.long[[k]]$Type<-factor(prov.summ.mo.long[[k]]$Type, levels=c("Prop.USA","Prop.international","Prop.totalinternational","Prop.domestic","Prop.provincial"))

  #repeat with total numbers instead of proportions
  #Go from wide to long (gather proportions)
  prov.summ.mo.long.2[[k]]<-prov.summ.mo[[k]] %>% gather(Type,N,N.provincial:N.totalinternational)
  
  #Plot the number of each type of transmission in a stacked bar plot
  prov.summ.mo.long.2[[k]]$Type<-factor(prov.summ.mo.long.2[[k]]$Type, levels=c("N.USA","N.international","N.totalinternational","N.domestic","N.provincial"))

######
} #end of k loop for BOOTS
#####

# Merge and summarize across bootstraps
prov.summ.mo.long.all<-rbind(prov.summ.mo.long[[1]],prov.summ.mo.long[[2]])
for (k in 3:b){
  prov.summ.mo.long.all<-rbind(prov.summ.mo.long.all,prov.summ.mo.long[[2]])
}
SUM.prov.summ.mo.long<-prov.summ.mo.long.all %>% dplyr::group_by (Province, month, Type) %>%
  dplyr::summarise(.groups="rowwise",mean.n=mean(Proportion),sd.n=sd(Proportion))

# Merge and summarize across bootstraps
prov.summ.mo.long.2.all<-rbind(prov.summ.mo.long.2[[1]],prov.summ.mo.long.2[[2]])
for (k in 3:b){
  prov.summ.mo.long.2.all<-rbind(prov.summ.mo.long.2.all,prov.summ.mo.long.2[[2]])
}
SUM.prov.summ.mo.long.2<-prov.summ.mo.long.2.all %>% dplyr::group_by (Province, month, Type) %>%
  dplyr::summarise(.groups="rowwise",mean.n=mean(N),sd.n=sd(N))

SUM.prov.summ.mo.long$Province<-factor(SUM.prov.summ.mo.long$Province,levels=ord.prov)
SUM.prov.summ.mo.long.2$Province<-factor(SUM.prov.summ.mo.long.2$Province,levels=ord.prov)

# summarize mean 

SUM.prov.summ.mo.long %>% dplyr::group_by(month,Type) %>%
  dplyr::summarise(.groups="rowwise",
                   mean.overall=mean(mean.n),
                   sd.overall=mean(sd.n)) %>%
  as.data.frame() %>%
  mutate(upper=mean.overall+(sd.overall/sqrt(6)*qt(.025, 6, lower.tail=F)),
         lower=mean.overall-(sd.overall/sqrt(6)*qt(.025, 6, lower.tail=F))) %>%
  filter(month %in% c("2020-03","2020-04","2020-05"))

SUM.prov.summ.mo.long %>% filter(month=="2020-04") %>% filter(Type=="totalinternational") %>% mutate(
  upper=mean.n+(sd.n/sqrt(10)*qt(.025, 10, lower.tail=F)),
  lower=mean.n-(sd.n/sqrt(10)*qt(.025, 10, lower.tail=F)))



## Plot proportion and total events separately
PZ1<-SUM.prov.summ.mo.long %>%
  filter(Type!="Prop.totalinternational") %>%
  # filter(month %in% c("2020-03","2020-04","2020-05","2020-06","2020-07")) %>%
  ggplot()+
  geom_bar(aes(x=month, y=mean.n, fill=Type, group=Type), stat="identity", position=position_stack(),alpha=0.8)+
  pubThemeDate+
  theme(strip.text = element_text(size=7),legend.text=element_text(size=7),
        legend.title = element_text(size=7,face="bold"),
        legend.position="top",
        plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"))+
  guides(fill=guide_legend(title.position="top"))+
  labs(x=NULL, y="Proportion of transmission events")+
  scale_fill_manual(values=prop.colz,name = "Transmission source", labels = c("USA","Other International","Between-province","Within-province"))+
  facet_wrap(.~(Province),nrow=1,scales="free_y")+
  # scale_x_discrete(breaks=c("2020-03","2020-04","2020-05","2020-06","2020-07"),
  #                  labels=c("Mar 2020","Apr 2020","May 2020","June 2020","July 2020"),expand=c(0.01,0))+
  scale_x_discrete(expand=c(0.01,0))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),expand=c(0.01,0))
PZ1
ggsave("results/BOOTS_TransmissionSourceByProvince_Proportions_monthly.png",height=3,width=8)
  

PZ2<-SUM.prov.summ.mo.long.2  %>%
  filter(Type!="N.totalinternational") %>%
  # filter(month %in% c("2020-03","2020-04","2020-05")) %>%
  ggplot()+
  geom_bar(aes(x=month, y=mean.n, fill=Type, group=Type), stat="identity", position=position_stack(),alpha=0.8)+
  pubThemeDate+
  theme(strip.text = element_text(size=7),
        legend.position="none",
        legend.spacing=unit(c(-1,-1,-1,-1),"cm"),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))+
  labs(x=NULL, y="Total transmssion events")+
  scale_fill_manual(values=prop.colz,name = "Transmission source", labels = c("USA","Other International","Between-province","Within-province"))+
  # scale_x_discrete(breaks=c("2020-03","2020-04","2020-05"),
  #                  labels=c("Mar 2020","Apr 2020","May 2020"),expand=c(0.01,0))+
  scale_x_discrete(expand=c(0.01,0))+
  scale_y_continuous(expand=c(0.01,0))+
  facet_wrap(.~(Province),ncol=6,scales="free_y")
PZ2
ggsave("results/BOOTS_TransmissionSourceByProvince_Counts_monthly.png",height=3,width=7)

plot.trans<-plot_grid(PZ1,PZ2,nrow=2, align="hv",axis="tl",rel_heights = c(1,1.1))
ggsave(plot.trans,file="results/BOOTS_monthal_TransmissionSourceByProvince.png",width=8,height=5,units = "in")

#export for mapping
# write.csv(SUM.prov.summ.mo.long,"../04_Map/BOOTS_Provinces_monthly_transmission.csv")

## Make both these plots using one object

SUM.prov.summ.mo.long$Type<-str_replace_all(SUM.prov.summ.mo.long$Type,"Prop.","")
SUM.prov.summ.mo.long.2$Type<-str_replace_all(SUM.prov.summ.mo.long.2$Type,"N.","")

SUM.prov.summ.mo.long$Param<-"Proportion"
SUM.prov.summ.mo.long.2$Param<-"Total"
SUM.prov.summ.mo.long.all<-rbind(SUM.prov.summ.mo.long,SUM.prov.summ.mo.long.2)
SUM.prov.summ.mo.long.all <-SUM.prov.summ.mo.long.all  [-which(SUM.prov.summ.mo.long.all  $Type=="totalinternational"),]

SUM.prov.summ.mo.long.all$Type<-factor(SUM.prov.summ.mo.long.all$Type, c("USA","international","domestic","provincial"))

PZ4<-SUM.prov.summ.mo.long.all  %>%
  filter(month %in% c("2020-03","2020-04","2020-05")) %>%
  ggplot()+
  geom_bar(aes(x=month, y=mean.n, fill=Type, group=Type), stat="identity", position=position_stack(),alpha=0.8)+
  pubThemeDate+
  theme(strip.text = element_text(size=9),
        legend.position="top",
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
        legend.text = element_text(size=9))+
  guides(fill=guide_legend(keywidth = 1,keyheight = 1))+
  labs(x=NULL, y="Sampled transmssion events")+
  scale_fill_manual(values=prop.colz,name = "Transmission source", labels = c("USA","Other International","Between-province","Within-province"))+
  scale_x_discrete(breaks=c("2020-03","2020-04","2020-05"),
                   labels=c("Mar 2020","Apr 2020","May 2020"),expand=c(.01,0))+
  scale_y_continuous(expand=c(0.01,0))+
  facet_grid(rows = vars(Param), cols=vars(Province),scales="free")
PZ4
#ggsave(file="results/BOOTS_monthly_TransmissionSourceByProvince.png",width=8,height=4,units = "in")

#export the dataframe so we can recreate the plot in the map project
# write.csv(SUM.prov.summ.mo.long.all  %>%
  # filter(month %in% c("2020-03","2020-04","2020-05")), file="DF/SUM.prov.summ.mo.long.all.csv")




## Identify singletons and link the parent state to canadian metadata
can.sing.boots<-replicate(b,vector())

for (k in 1:b){
  
  #join canadian metadata by tiplabel to can.states.boots[[k]]
  if (!"par.lik" %in% colnames(cana.boots[[k]])){ #this is a check to make sure don't do it twice
    cana.boots[[k]]<-left_join(cana.boots[[k]],can.states.boots[[k]],by="tip.label")
  }
  
  #make a vector of all unique descendants within all sublins
  all.des.un<-unique(sublin.long[[k]]$Descendant.AccessionIDs)
  
  #subset these to Canadian
  can.des.un<-all.des.un[which(all.des.un %in% cana.boots[[k]]$GISAID_ID)]

  #Can representation
  table(cana.boots[[k]]$state[cana.boots[[k]]$GISAID_ID%in%can.des.un])
  
  #change names
  cana.boots[[k]]$par.state<-str_replace_all(cana.boots[[k]]$par.state,"Canada_","")
  cana.boots[[k]]$tip.state<-str_replace_all(cana.boots[[k]]$tip.state,"Canada_","")
  
  #find Canadian tip labels in can.states.boots[[k]] that aren't a descendant in any of the sublienages
  #subset the singletons
  can.sing.boots[[k]]<-cana.boots[[k]][which(!cana.boots[[k]]$GISAID_ID %in% can.des.un),]

  #apply lookup table
  for (i in 1:nrow(lookup.geo)){
    if(lookup.geo$new.loc[i] != lookup.geo$og.loc[i]){ #if they don't match, go replace instances
      #find matches in node loc, and replace with new loc
      pr.ma<-which(can.sing.boots[[k]]$tip.state==lookup.geo$og.loc[i])
      if (length(pr.ma)>0) {can.sing.boots[[k]]$tip.state[pr.ma]<-lookup.geo$new.loc[i]}
      par.ma<-which(can.sing.boots[[k]]$par.state==lookup.geo$og.loc[i])
      if (length(par.ma)>0) {can.sing.boots[[k]]$par.state[par.ma]<-lookup.geo$new.loc[i]; next}
    }
  }
  
  ## different def'n of singleton (can be inside a sublineage)
  singlez2<-cana.boots[[k]][-which(cana.boots[[k]]$par.state %in% provs),]
  nrow(singlez2) #1407
  nrow(singlez2)-nrow(can.sing.boots[[k]]) #=re-intro singletons
}



## summarize singletons for results text

#bootstraps summary of singletons
sing.tot<-c()
for (k in 1:b){
   sing.tot<-c(sing.tot,nrow(can.sing.boots[[k]] ))
}
mean.95ci.X(sing.tot,0) #1342 (1311 - 1374)"

#what was the TOTAL number of intros (singles + sublins)
tot.tot<-c()
for (k in 1:b){
   tot.tot<-c(tot.tot,(nrow(can.sing.boots[[k]] ) + nrow(sum.boots[[k]])))
}
mean.95ci.X(tot.tot,0) # "1754 (1722 - 1786)"

#what proportion were singeltons
mean.95ci.X((sing.tot/tot.tot*100),1) #"76.5 (75.9 - 77.1)"

#what percent did the few large intros account for?
mean.95ci.X( hund/tot.tot*100, 1)


# PROPORTION OF importations leading to sublineage analysis
## proportion of importations from each country of origin leading to singleton vs an introduction BY COUNTRY

allimportz.boots<-replicate(b,vector())
for (k in 1:b){
  #what proportion of importations from each country resulted in a singleton vs an introdctuion leading to ongoing transmission, by month?
  singz<-as.data.frame.matrix(table(can.sing.boots[[k]]$par.state,can.sing.boots[[k]]$month))
  subz<-as.data.frame.matrix(table(sum.boots[[k]]$Parent.Location,sum.boots[[k]]$month))
  
  #go from wide to long
  singz$country<-rownames(singz)
  singz.long<-singz %>% pivot_longer(1:(ncol(singz)-1),names_to = "month")
  colnames(singz.long)[3]<-"singz"
    
  subz$country<-rownames(subz)
  subz.long<-subz %>% pivot_longer(1:(ncol(subz)-1),names_to = "month")
  colnames(subz.long)[3]<-"subz"
  
  if(nrow(singz.long)>0){
    allimportz<-full_join(singz.long,subz.long,by=c("country","month"))
    allimportz<-allimportz %>% mutate(total=subz+singz, prop.subz=subz/total)
  }
  if(nrow(singz.long)==0){
    allimportz<-subz.long
    allimportz<-allimportz %>% mutate(singz=0, total=subz, prop.subz=1)
  }
    
  allimportz.boots[[k]]<-allimportz[,c("country","month","singz","subz","total","prop.subz")]
}

#summarize across boots
summ.allimportz.boots<-rbind(allimportz.boots[[1]],allimportz.boots[[2]])
for (k in 3:b){
  summ.allimportz.boots<-rbind(summ.allimportz.boots, allimportz.boots[[k]])
}
summary.allimportz.boots<-summ.allimportz.boots %>% 
  dplyr::group_by (country, month) %>%
  dplyr::summarise(.groups="rowwise",
                   mean.propsubz=mean(prop.subz)) %>%
  dplyr::ungroup() %>% as.data.frame()

  # #statistical test (different over time)
  kruskal.test(summary.allimportz.boots$month, summary.allimportz.boots$mean.propsubz) #not signf
  kruskal.test(summary.allimportz.boots$country, summary.allimportz.boots$mean.propsubz) #not signif
  
## consider showing as boxplot with points behind
#USE IN PUBLlCATION GROB BELOW
# PX2<-summary.allimportz.boots  %>%
#   ggplot()+
#     geom_boxplot(aes(x=month,y=mean.propsubz))+
#     geom_point(aes(x=month,y=mean.propsubz,group=country,color=country),position=position_dodge2(width=0.4),alpha=0.9)+
#     pubThemeDate+
#     theme(legend.position="top")+
#   labs(x=NULL,y="Prop. of importations resulting in sublineage")+
#   guides(color=guide_legend(title="Origin\nlocation"))+
#   GlobColScale+
#   coord_cartesian(clip="off")
# PX2
# #ggsave("results/Boxplot with POINT Prop. of importations resulting in sublineage BY ORIGIN.png",width=5,height=4,units="in")

## proportion of importations into each province leading to singleton vs an introduction BY PROVINCE
allimportz.boots<-replicate(b,vector())
for (k in 1:b){
  #what proportion of importations from each country resulted in a singleton vs an introdctuion leading to ongoing transmission, by month?
  singz<-as.data.frame.matrix(table(can.sing.boots[[k]]$tip.state,can.sing.boots[[k]]$month))
  subz<-as.data.frame.matrix(table(sum.boots[[k]]$Node.Location,sum.boots[[k]]$month))
  
  #go from wide to long
  singz$country<-rownames(singz)
  singz.long<-singz %>% pivot_longer(1:(ncol(singz)-1),names_to = "month")
  colnames(singz.long)[3]<-"singz"
    
  subz$country<-rownames(subz)
  subz.long<-subz %>% pivot_longer(1:(ncol(subz)-1),names_to = "month")
  colnames(subz.long)[3]<-"subz"
  
  allimportz<-left_join(singz.long,subz.long,by=c("country","month"))
  allimportz.boots[[k]]<-allimportz %>% mutate(total=subz+singz, prop.subz=subz/total)
}

#summarize across boots
summ.allimportz.boots<-rbind(allimportz.boots[[1]],allimportz.boots[[2]])
for (k in 3:b){
  summ.allimportz.boots<-rbind(summ.allimportz.boots, allimportz.boots[[k]])
}
summary.allimportz.boots<-summ.allimportz.boots %>% 
  dplyr::group_by (country, month) %>%
  dplyr::summarise(.groups="rowwise",
                   mean.propsubz=mean(prop.subz)) %>%
  dplyr::ungroup() %>% as.data.frame()

  # #statistical test (different over time)
  kruskal.test(summary.allimportz.boots$month, summary.allimportz.boots$mean.propsubz) #0.4528
  kruskal.test(summary.allimportz.boots$country, summary.allimportz.boots$mean.propsubz) #not signif = 0.5896
  
#boxplot with points
# PX1<-summary.allimportz.boots %>%
#   ggplot()+
#     geom_boxplot(aes(x=month,y=mean.propsubz))+
#     geom_point(aes(x=month,y=mean.propsubz,group=country,color=country),position=position_dodge2(width=0.4),alpha=0.9)+
#   pubThemeDate+
#   theme(legend.position="top")+
#   labs(x=NULL,y="Prop. importations resulting in sublineage")+
#   guides(color=guide_legend(title="Province"))+
#   GlobColScale+
#   coord_cartesian(clip="off")
# PX1
# #ggsave("results/Boxplot with POINT Prop. of importations resulting in sublineage BY PROV.png",width=5,height=4,units="in")


##figure of proportion of importations resulting in sublineage
# plot_grid(PX1, PX2, labels=c("A","B"))
#ggsave(file="results/PUBL_boxplot with POINT Prop. of importations resulting in dom trans BY PROV and ORIGIN.png",width=8,height=4,units = "in")

# Plot singletons as in Main figures
## Alluvial for singletons

#summarize the total number and percent of intros by 1) par. loc, 2) prov of intro, maybe also 3) lineage 
#prep list item
sing.Par.Prov.l<-replicate(b,vector())

#go through each list item/subsample
## tabulate instances of location pairs
rem<-c()
for (k in 1:b){
  if (nrow(can.sing.boots[[k]])<1){rem<-c(rem,k);next}
  sing.Par.Prov.l[[k]]<-can.sing.boots[[k]] %>% dplyr::group_by (par.state, tip.state) %>% dplyr::summarize (.groups="rowwise", n.Par= n()) %>% as.data.frame() %>% mutate(perc.Par=(n.Par/sum(n.Par))*100)

  ## add a column for subsample
  sing.Par.Prov.l[[k]]$SampleSet<-paste("Sample",BOOTS[k])
}

## rbind the summaries except rem
sum.Par.Prov<-rbind(sing.Par.Prov.l[[1]],sing.Par.Prov.l[[2]])
for (k in BOOTS[-c(1,2,rem)]){
  sum.Par.Prov<-rbind(sum.Par.Prov,sing.Par.Prov.l[[k]])
}

#summarize the mean and range for each of 1), 2) and 3)
sing.Par.Prov.summary<-sum.Par.Prov %>% dplyr::group_by(par.state, tip.state) %>%
  dplyr::summarize(.groups="rowwise",
                   mean.n=round(mean(n.Par)),
                   sd.n=round(sd(n.Par,na.rm=T),digits=2),
                   mean.perc=round(mean(perc.Par),digits=2),
                   sd.perc=round(sd(perc.Par,na.rm=T),digits=2)) %>%
  as.data.frame()

tot.Par<-sing.Par.Prov.summary %>% dplyr::group_by(par.state) %>%
  dplyr::summarise(.groups = "rowwise", totalImports=sum(mean.n)) %>%
  as.data.frame()

## Alluvial plot for Figure 2
#make a "subject column"
sing.Par.Prov.summary$subject<-1:nrow(sing.Par.Prov.summary)
# sing.Par.Prov.summary<-sing.Par.Prov.summary[-which(sing.Par.Prov.summary$mean.n<1),]

#order the geos
ord.count<-sing.Par.Prov.summary%>% dplyr::group_by(par.state) %>%
  dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
ord.count<-ord.count[rev(order(ord.count$n)),'par.state']
sing.Par.Prov.summary$par.state<-factor(sing.Par.Prov.summary$par.state,levels=ord.count)

ord.prov<-sing.Par.Prov.summary%>% dplyr::group_by(tip.state) %>%
  dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
ord.prov<-ord.prov[rev(order(ord.prov$n)),'tip.state']
sing.Par.Prov.summary$tip.state<-factor(sing.Par.Prov.summary$tip.state,levels=ord.prov)


#make it long
sing.Par.Prov.summary.long<-sing.Par.Prov.summary %>% pivot_longer(1:2, names_to = "geo.type", values_to = "geo")
#order the geo types
sing.Par.Prov.summary.long$geo.type<-factor(sing.Par.Prov.summary.long$geo.type,levels=c("par.state","tip.state"),labels=c("Global origin","Canadian province"))

## Alluvial plot
P1<- ggplot(sing.Par.Prov.summary.long,
       aes(x = geo.type, stratum = geo, alluvium = subject,
           y = mean.n,
           fill = geo, label = geo)) +
  scale_x_discrete(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_flow(alpha = .6,width=0.45) +
  geom_stratum(alpha = .8,width=0.45) +
  geom_text(stat = "stratum", size = 3.4,min.y=0,fontface="bold") +
  pubTheme+
  theme(legend.position = "none", axis.line = element_blank(), text =element_text(size=14),
        axis.ticks.x = element_blank(),axis.text.x = element_text(hjust=c(0.5,0.6)))+
  labs(x=NULL,y="# singletons")+
  GlobFillScale
P1
ggsave("results/sankeyPlot.singletons.Parent.Node.png",height=10,width=5,units="in")


## singelton mosaic of lineages
#add lineage group to singeltons object
for (k in 1:b){
  for (i in 1:nrow(can.sing.boots[[k]])){
    can.sing.boots[[k]]$Lineage.grp[i]<-lookup.lin$lineagegroup[which(lookup.lin$lineage==can.sing.boots[[k]]$Lineage[i])]
  
  }
}
# table(can.sing.boots[[k]]$Lineage.grp)

#prep list item
sing.Par.L.l<-replicate(b,vector())

#go through each list item/subsample
## tabulate instances of location pairs
for (k in 1:b){
  sing.Par.L.l[[k]]<-can.sing.boots[[k]] %>% dplyr::group_by (par.state, Lineage.grp) %>% dplyr::summarize (.groups="rowwise", n.Par= n()) %>% as.data.frame() %>% mutate(perc.Par=(n.Par/sum(n.Par))*100)
  
  ## add a column for subsample
  sing.Par.L.l[[k]]$SampleSet<-paste("Sample",BOOTS[k])

}

sing.Par.L<-rbind(sing.Par.L.l[[1]],sing.Par.L.l[[2]])
for (k in 3:b){
  sing.Par.L<-rbind(sing.Par.L,sing.Par.L.l[[k]])
}


#summarize the mean and range for each of 1), 2) and 3)
sing.Par.L.summary<-sing.Par.L %>% dplyr::group_by(par.state,  Lineage.grp) %>%
  dplyr::summarize(.groups="rowwise",
                   mean.n=round(mean(n.Par)),
                   sd.n=round(sd(n.Par,na.rm=T),digits=2),
                   mean.perc=round(mean(perc.Par),digits=2),
                   sd.perc=round(sd(perc.Par,na.rm=T),digits=2)) %>%
  as.data.frame()

#go from having frequencies to having indiviudal occurences

#order as above
sing.Par.L.summary$par.state<-factor(sing.Par.L.summary$par.state,levels=rev(ord.count))
sing.Par.L.summary$Lineage.grp<-factor(sing.Par.L.summary$Lineage,levels=(row.names(lin.group.col)))

P2<-sing.Par.L.summary %>% 
  ggplot()+
  geom_mosaic(aes(weight=mean.n,x=product(Lineage.grp,par.state),
                  fill=Lineage.grp), offset = .002,alpha=0.9)+
  pubTheme+
  theme(axis.text.y = element_text(size=12),
        axis.line=element_blank(),
        legend.position = "top", axis.ticks=element_blank(),
        legend.text=element_text(size=rel(1.1)),
        legend.title = element_text(size=rel(1.2)),
        axis.text.x = element_blank(),
        panel.background = element_blank(), 
        axis.title.y=element_text(vjust=-11,hjust=0.6,size=rel(1.4)), 
        plot.margin=unit(c(0.5,0.5,6,0.1),"cm"))+
    scale_fill_manual(name = "Lineage group",values = lin.grp.ch,na.value="grey60")+
  labs(y="",x="Global origin")+
  scale_x_productlist(expand=c(0.02,0))+
  scale_y_productlist(expand=c(0,0))+
  guides(fill=guide_legend(title.position = "top", keywidth = 1.1,keyheight = 1.1))+
  # annotate(geom="text",y=0.01,x=1.015,label="A*")+
  # annotate(geom="text",y=0.1,x=1.015,label="B*")+
  # annotate(geom="text",y=0.55,x=1.015,label="B.1*")+
  # annotate(geom="text",y=0.92,x=1.015,label="B.1.1*")+
  coord_flip()
P2
ggsave("results/SING_MosaicPlot_LinGrpParentLoc2.png",height=4,width=4,units="in")


## Composite plot of singletons alluvial and mosaic, exactly as in sublineages 
## Add a blank grob to mimic the same layout with sublineage size
# 
# right_col <- plot_grid(P2, labels = c('B'), label_size = 14,ncol=1, align="hv", axis="lr")
# 
# # right_col <- plot_grid(P2,NULL, labels = c('B', NULL), label_size = 14,ncol=1, align="hv", axis="lr")
# plot_grid(P1, right_col, labels = c('A', ''), label_size = 14, ncol = 2, align="h",axis="b",rel_widths = c(0.55,0.45))
# 
# ggsave(file="results/Singletons_sankeyMosaicSize_MOD.png",width=8.5,height=10,units = "in")



## singletons: rolling rates

#mid script lists for rolling means
par.state.sing.boots<-replicate(n=b,vector())
par.state.sing.boots.full<-replicate(n=b,vector())
tip.state.sing.boots<-replicate(n=b,vector())
tip.state.sing.boots.full<-replicate(n=b,vector())


for (k in 1:b){
  #make sure this a date
  can.sing.boots[[k]]$date.lsd.full<-as.Date(can.sing.boots[[k]]$date.lsd.full)
  
  #### Calculate a rolling 7-day mean for origins ####
  
  ## count the importations by origin location over time
  par.state.sing.boots[[k]]<-can.sing.boots[[k]] %>%
    dplyr::select(par.state, Lineage,date.lsd.full) %>%
    dplyr::group_by(date.lsd.full, par.state) %>%
    dplyr::summarize(.groups="rowwise", total=n()) %>%
    dplyr::arrange(desc(par.state)) %>% 
    dplyr::group_by(par.state) 
  
  #need to add rows for missing dates
  alldays<-seq(ymd(first(sort(par.state.sing.boots[[k]]$date.lsd.full))),
      ymd(last(sort(par.state.sing.boots[[k]]$date.lsd.full))),
      by='1 day')
  
  #make a empty df in same structure as above then populate it
  nL<-length(unique(par.state.sing.boots[[k]]$par.state))
  nD<-length(alldays)
  par.state.sing.boots.full[[k]]<-data.frame(date.lsd.full=rep(alldays,times=nL),
par.state = sort(rep(unique(par.state.sing.boots[[k]]$par.state),times=nD)),
                                 total=0)
  # nrow(par.state.sing.boots[[k]].empty)==nD*nL      
  
  #populate it
  for (i in 1:nrow(par.state.sing.boots.full[[k]])){
    #look for a match
    match<-which(par.state.sing.boots[[k]]$par.state==par.state.sing.boots.full[[k]]$par.state[i] &
            par.state.sing.boots[[k]]$date.lsd.full==par.state.sing.boots.full[[k]]$date.lsd.full[i])
    if(length(match)==0) next #no match, no change
    #else, replace:
    par.state.sing.boots.full[[k]]$total[i]<-par.state.sing.boots[[k]]$total[match]
  }
  
  # sum(par.state.sing.boots.full[[k]]$total[par.state.sing.boots.full[[k]]$par.state=="USA"])==sum(par.state.sing.boots[[k]]$total[par.state.sing.boots[[k]]$par.state=="USA"])
  
  par.state.sing.boots.full[[k]]<-par.state.sing.boots.full[[k]] %>% 
    dplyr::mutate(intros_mean7d = zoo::rollmean(total, k = 7, fill = NA),
                  intros_mean14d = zoo::rollmean(total, k = 14, fill = NA),
                  intros_median7d = zoo::rollmedian(total, k = 7, fill = NA),
                  intros_sum7d = zoo::rollsum(total, k = 7, fill = NA,align="right")) %>% #right align to sum all prev
    #rolling mean and median weekly importation rate
    dplyr::mutate(intros_meansum7d = zoo::rollmean(intros_sum7d, k = 7, fill = NA), 
                  intros_mediansum7d = zoo::rollmedian(intros_sum7d, k = 7, fill = NA),) %>%
    dplyr::ungroup()
   
  ## count the importations by tip.state over time
  tip.state.sing.boots[[k]]<-can.sing.boots[[k]] %>%
    dplyr::select(tip.state, Lineage,date.lsd.full) %>%
    group_by(date.lsd.full, tip.state) %>%
    dplyr::summarize(.groups="rowwise", total=n()) %>%
    dplyr::arrange(desc(tip.state)) %>% 
    dplyr::group_by(tip.state) 
  
  #need to add rows for missing dates
  alldays<-seq(ymd(first(sort(tip.state.sing.boots[[k]]$date.lsd.full))),
      ymd(last(sort(tip.state.sing.boots[[k]]$date.lsd.full))),
      by='1 day')
  
  #make a empty df in same structure as above then populate it
  nL<-length(unique(tip.state.sing.boots[[k]]$tip.state))
  nD<-length(alldays)
  tip.state.sing.boots.full[[k]]<-data.frame(date.lsd.full=rep(alldays,times=nL),
                                 tip.state=sort(rep(unique(tip.state.sing.boots[[k]]$tip.state),times=nD)),
                                 total=0)
  # nrow(tip.state.sing.boots[[k]].empty)==nD*nL      
  
  #populate it
  for (i in 1:nrow(tip.state.sing.boots.full[[k]])){
    #look for a match
    match<-which(tip.state.sing.boots[[k]]$tip.state==tip.state.sing.boots.full[[k]]$tip.state[i] &
            tip.state.sing.boots[[k]]$date.lsd.full==tip.state.sing.boots.full[[k]]$date.lsd.full[i])
    if(length(match)==0) next #no match, no change
    #else, replace:
    tip.state.sing.boots.full[[k]]$total[i]<-tip.state.sing.boots[[k]]$total[match]
  }
  
  # sum(tip.state.sing.boots.full[[k]]$total[tip.state.sing.boots.full[[k]]$tip.state=="USA"])==sum(tip.state.sing.boots[[k]]$total[tip.state.sing.boots[[k]]$tip.state=="USA"])
  
  tip.state.sing.boots.full[[k]]<-tip.state.sing.boots.full[[k]] %>% 
    dplyr::mutate(intros_mean7d = zoo::rollmean(total, k = 7, fill = NA),
                  intros_mean14d = zoo::rollmean(total, k = 14, fill = NA),
                  intros_median7d = zoo::rollmedian(total, k = 7, fill = NA),
                  intros_sum7d = zoo::rollsum(total, k = 7, fill = NA,align="right")) %>% #right align to sum all prev
    #rolling mean and median weekly importation rate
    dplyr::mutate(intros_meansum7d = zoo::rollmean(intros_sum7d, k = 7, fill = NA), 
                  intros_mediansum7d = zoo::rollmedian(intros_sum7d, k = 7, fill = NA),) %>%
    dplyr::ungroup()

}

## Summarize the rolling means overall

## BY ORIGINS

# #go through each list item/subsample
# for (k in 1:b){
#   ## add a column for subsample
#   par.state.sing.boots.full[[k]]$SampleSet<-paste("Sample",BOOTS[k])
# }
# 
# ## rbind the summaries
# sing.Par.Roll<-rbind(par.state.sing.boots.full[[1]],par.state.sing.boots.full[[2]])
# for (k in 3:b){
#   sing.Par.Roll<-rbind(sing.Par.Roll,par.state.sing.boots.full[[k]])
# }
# 
# #summarize the mean and confint
# sing.Par.Roll.summary<-sing.Par.Roll %>% dplyr::group_by(date.lsd.full, par.state) %>%
#   dplyr::summarize(.groups="rowwise",
#                    intros_meansum7d.mean=(mean(intros_meansum7d,na.rm=T)),
#                    intros_meansum7d.sd=sd(intros_meansum7d,na.rm=T)) %>%
#   as.data.frame()
# 
# #order the geos
# sing.Par.Roll.summary$par.state<-factor(sing.Par.Roll.summary$par.state,levels=ord.count)
# 
# ## BY PROVINCE
# 
# #go through each list item/subsample
# for (k in 1:b){
#   ## add a column for subsample
#   tip.state.sing.boots.full[[k]]$SampleSet<-paste("Sample",BOOTS[k])
# }
# 
# ## rbind the summaries
# sing.Prov.Roll<-rbind(tip.state.sing.boots.full[[1]],tip.state.sing.boots.full[[2]])
# for (k in 3:b){
#   sing.Prov.Roll<-rbind(sing.Prov.Roll,tip.state.sing.boots.full[[k]])
# }
# 
# #summarize the mean and confint
# sing.Prov.Roll.summary<-sing.Prov.Roll %>% dplyr::group_by(date.lsd.full, tip.state) %>%
#   dplyr::summarize(.groups="rowwise",
#                    intros_meansum7d.mean=(mean(intros_meansum7d,na.rm=T)),
#                    intros_meansum7d.sd=(sd(intros_meansum7d,na.rm=T))) %>%
#   as.data.frame()
# 
# #order the geos
# sing.Prov.Roll.summary$tip.state<-factor(sing.Prov.Roll.summary$tip.state,levels=ord.prov)
# 
# ## Plot singletons over time using rolling rates, by tmrca, by origin or dest
# #dataframe for travel restrictions
# trav.df<-data.frame(x1=as.Date("2020-03-21"),y1=-0.5, x2=as.Date("2020-03-16"),y2=-1.5)
# 
# #rolling importation rate, by origin
# plot_lim_y<-max(sing.Par.Roll.summary$intros_meansum7d.mean,na.rm=T)*10
# p1<-sing.Par.Roll.summary %>%
#   ggplot(aes(x=date.lsd.full,y=intros_meansum7d.mean,group=par.state,fill=par.state))+
#   geom_density(stat="identity", position="stack",lwd=0,alpha=0.9)+
#   GlobFillScale+
#   pubThemeDate+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         text=element_text(size=10,face="bold"),
#         legend.position = "none")+
#   scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits = c(as.Date("2020-01-01"),as.Date("2021-01-01")) ) +
#   labs(x="Sampling date",y="# singletons introduced per week", fill="Origin Location")+
#   #travel restrictions
#   geom_vline(xintercept=as.Date("2020-03-21"),color="grey20",linetype=2)+
#   annotate("text",label="Maximum\nstringency",
#            x=as.Date("2020-03-30"),y=plot_lim_y,size=4,hjust=0,vjust=1,color="grey20")+
#   #minimal stringency
#   geom_vline(xintercept=as.Date("2020-10-10"),color="grey20",linetype=2)+
#   annotate("text",label="Reduced\nstringency",
#            x=as.Date("2020-10-12"),y=plot_lim_y,size=4,hjust=0,vjust=1,color="grey20")+
#   scale_y_continuous(expand=c(0,0),limits=c(-0.1,plot_lim_y),breaks=seq(0,plot_lim_y,0.1))
# p1
# 
# ##fake plots to take the legends
# #rolling importation rate, by origin
# p1.f<-sing.Par.Roll.summary %>%
#   ggplot(aes(x=date.lsd.full,y=intros_meansum7d.mean,group=par.state,fill=par.state))+
#   geom_density(stat="identity", position="stack",lwd=0,alpha=0.9)+
#   GlobFillScale+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         text=element_text(size=9,face="bold"),
#         legend.position = c(0.1,1),
#         legend.margin=margin(0,0,0,0),
#         legend.background = element_blank(),
#         legend.justification = c(0,1))+
#   guides(fill = guide_legend(keywidth = 0.8,keyheight=0.8,title.position = "top",title="Global origin",legend.spacing=0,ncol=4))
# 
# p1.guide<-get_legend(p1.f)
# 
# 
# #rolling importation rate, by province destination
# p2<-sing.Prov.Roll.summary %>%
#   ggplot(aes(x=date.lsd.full,y=intros_meansum7d.mean,group=tip.state,fill=tip.state))+
#   geom_density(stat="identity", position="stack",lwd=0,alpha=0.9)+
#   GlobFillScale+
#   pubThemeDate+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         legend.position = "none")+
#   scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits = c(as.Date("2020-01-01"),as.Date("2021-01-01")) ) +
#   labs(x="Sampling date",y="# singletons introduced per week")+
#   #travel restrictions
#   geom_vline(xintercept=as.Date("2020-03-21"),color="grey20",linetype=2)+
#   annotate("text",label="Maximum\nstringency",
#            x=as.Date("2020-03-30"),y=400,size=4,hjust=0,vjust=1,color="grey20")+
#   #minimal stringency
#   geom_vline(xintercept=as.Date("2020-10-10"),color="grey20",linetype=2)+
#   annotate("text",label="Reduced\nstringency",
#            x=as.Date("2020-10-12"),y=400,size=4,hjust=0,vjust=1,color="grey20")+
#   scale_y_continuous(expand=c(0,0),limits=c(-1,410),breaks=seq(0,400,100))
# p2
# #fake plot to take legend
# p2.f<-sing.Prov.Roll.summary %>%
#   ggplot(aes(x=date.lsd.full,y=intros_meansum7d.mean,group=tip.state,fill=tip.state))+
#   geom_density(stat="identity", position="stack",lwd=0,alpha=0.9)+  
#   GlobFillScale+
#   theme(axis.text.x=element_text(angle=45,hjust = 1),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background=element_rect("grey95"),
#         text=element_text(size=9,face="bold"),
#         legend.position = c(0.1,1),
#         legend.margin=margin(0,0,0,0),
#         legend.background = element_blank(),
#         legend.justification = c(0,1))+
#     guides(fill = guide_legend(keywidth = 0.8,keyheight=0.8,title.position = "top",title="Province",legend.spacing=0,ncol=3))
# p2.guide<-get_legend(p2.f)
# 
# #make a list of grobs
# ## add in the Canadian province representation plots
# plot.row<-plot_grid(p1,p2,ncol=2,labels=c("A" ,"B"))
# plot.leg<-plot_grid(p1.guide,p2.guide,ncol=2,align="h")
# plot.all<-plot_grid(plot.leg,plot.row,nrow=2,rel_heights = c(0.25,1),align = "v")
# #ggsave(plot.all,file="results/Singletons_OriginsDestinationsInclRatesOVERTIME.png",width=7,height=4.5,units = "in")


## Likelihood sensitivity analysis

likes<-replicate(b,vector())
for (k in 1:b){
  likes[[k]]<-sum.boots[[k]][,c('Node.Likelihood','Parent.Likelihood')]
}
sum.likes<-rbind(likes[[1]],likes[[2]])
for (k in 3:b){
  sum.likes<-rbind(sum.likes,likes[[k]])
}

#summary stats
summary(sum.likes$Node.Likelihood)
summary(sum.likes$Parent.Likelihood)

Nlik<-nrow(sum.likes)

#text for plot
p.5.n.5<-round(nrow(sum.likes[sum.likes$Parent.Likelihood<0.5 & sum.likes$Node.Likelihood<0.5,]) / Nlik * 100 ,digits=2)
p.5.n.9<-round(nrow(sum.likes[sum.likes$Parent.Likelihood<0.5 & sum.likes$Node.Likelihood<0.9 & sum.likes$Node.Likelihood>=0.5,]) / Nlik * 100,digits=2)
p.9.n.5<-round(nrow(sum.likes[sum.likes$Parent.Likelihood<0.9 & sum.likes$Parent.Likelihood>=0.5 & sum.likes$Node.Likelihood<0.5,]) / Nlik * 100,digits=2)
p.9.n.9<-round(nrow(sum.likes[sum.likes$Parent.Likelihood<0.9 & sum.likes$Parent.Likelihood>=0.5 & sum.likes$Node.Likelihood>=0.5 & sum.likes$Node.Likelihood<0.9,]) / Nlik * 100,digits=2)

p.5<-round(nrow(sum.likes[sum.likes$Parent.Likelihood<0.5 & sum.likes$Node.Likelihood>=0.9,]) / Nlik * 100,digits=2)
n.5<-round(nrow(sum.likes[sum.likes$Parent.Likelihood>=0.9 & sum.likes$Node.Likelihood<0.5,]) / Nlik * 100,digits=2)
p.9<-round(nrow(sum.likes[sum.likes$Parent.Likelihood<0.9 & sum.likes$Parent.Likelihood>=0.5 & sum.likes$Node.Likelihood>=0.9,]) / Nlik * 100,digits=2)
n.9<-round(nrow(sum.likes[sum.likes$Parent.Likelihood>=0.9 & sum.likes$Node.Likelihood>=0.5 & sum.likes$Node.Likelihood<0.9,]) / Nlik * 100,digits=2)

p.n<-round(nrow(sum.likes[sum.likes$Parent.Likelihood>=0.9 & sum.likes$Node.Likelihood>=0.9,]) / Nlik * 100,digits=2)



#Plot this with 0.5 and 0.9 cutoffs shown
ggplot(sum.likes)+
  geom_point(aes(x=Node.Likelihood, y=Parent.Likelihood),color="cyan4",alpha=0.5)+
  geom_vline(xintercept=c(0.5,0.9),linetype="dashed", color="red",size=0.5)+
  geom_hline(yintercept=c(0.5,0.9),linetype="dashed", color="blue",size=0.5)+
  pubTheme+
  annotate("text",x=0.5,y=0.5,label=paste(p.5.n.5,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  annotate("text",x=0.9,y=0.5,label=paste(p.5.n.9,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  annotate("text",x=0.5,y=0.9,label=paste(p.9.n.5,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  annotate("text",x=0.9,y=0.9,label=paste(p.9.n.9,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  annotate("text",x=1,y=0.5,label=paste(p.5,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  annotate("text",x=1,y=0.9,label=paste(p.9,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  annotate("text",x=0.5,y=1,label=paste(n.5,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  annotate("text",x=0.9,y=1,label=paste(n.9,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  annotate("text",x=1,y=1,label=paste(p.n,"%",sep=""),size=3.5,hjust=1,vjust=1, fontface="bold")+
  labs(x="Introduction node likelihood",y="Parent node likelihood")+
  theme(panel.background = element_rect(fill="grey95"))

ggsave("results/boot.Likelihoods.Nodes.Parents.png",width=5,height=5)

#### note: should we apply a cutoff here? #####

#### Between-province dynamics ####

## looking just at cana (ie sampled tips from Canada): 
## what were the relative flows in/out of provinces?
k=1

sum.bw.Prov.l<-replicate(b,vector())
#go through each list item/subsample
## tabulate instances of location pairs
for (k in 1:b){
  can.can<-cana.boots[[k]][which(cana.boots[[k]]$par.state %in% provs),] #temp object
  #apply lookup to maritimes
  can.can$tip.state<-str_replace_all(can.can$tip.state, "Nova Scotia|New Brunswick|Newfoundland and Labrador", "Maritimes")
  sum.bw.Prov.l[[k]]<-can.can %>% dplyr::group_by (par.state, tip.state) %>% dplyr::summarize (.groups="rowwise", n.Par= n()) %>% as.data.frame() %>% mutate(perc.Par=(n.Par/sum(n.Par))*100)

  ## add a column for subsample
  sum.bw.Prov.l[[k]]$SampleSet<-paste("Sample",BOOTS[k])

}

# bind the summaries
sum.bw.Prov<-bind_rows(sum.bw.Prov.l, .id = "SampleSet")

#summarize the mean and range for each of 1), 2) and 3)
sum.bw.Prov.summary<-sum.bw.Prov %>% dplyr::group_by(par.state, tip.state) %>%
  dplyr::summarize(.groups="rowwise",
                   mean.n=round(mean(n.Par)),
                   sd.n=round(sd(n.Par,na.rm=T),digits=2),
                   mean.perc=round(mean(perc.Par),digits=2),
                   sd.perc=round(sd(perc.Par,na.rm=T),digits=2)) %>%
  as.data.frame() 

tot.Par<-sum.bw.Prov.summary %>% dplyr::group_by(par.state) %>%
  dplyr::summarise(.groups = "rowwise", 
                   totalImports=sum(mean.n)) %>%
  as.data.frame()


#REMOVE WITHIN PROVINCE
# head(sum.bw.Prov.summary)
rem<-c()
for (i in 1:nrow(sum.bw.Prov.summary)){
  if (sum.bw.Prov.summary$par.state[i]==sum.bw.Prov.summary$tip.state[i]) {rem<-c(rem,i)}
}
sum.bw.Prov.summary<-sum.bw.Prov.summary[-rem,]

## Alluvial plot for Figure 2
#make a "subject column"
sum.bw.Prov.summary$subject<-1:nrow(sum.bw.Prov.summary)
# sum.bw.Prov.summary<-sum.bw.Prov.summary[-which(sum.bw.Prov.summary$mean.n<1),]

#order the geos
ord.prov.2<-sum.bw.Prov.summary%>% dplyr::group_by(par.state) %>%
  dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
ord.prov.2<-ord.prov.2[rev(order(ord.prov.2$n)),'par.state']
sum.bw.Prov.summary$par.state<-factor(sum.bw.Prov.summary$par.state,levels=ord.prov.2)
 
ord.prov.3<-sum.bw.Prov.summary%>% dplyr::group_by(tip.state) %>%
  dplyr::summarise(.groups="rowwise",n=sum(mean.n)) %>% as.data.frame()
ord.prov.3<-ord.prov.3[rev(order(ord.prov.3$n)),'tip.state']
sum.bw.Prov.summary$tip.state<-factor(sum.bw.Prov.summary$tip.state,levels=ord.prov.3)

#make it long
sum.bw.Prov.summary.long<-sum.bw.Prov.summary %>% pivot_longer(1:2, names_to = "geo.type", values_to = "geo")
#order the geo types
sum.bw.Prov.summary.long$geo.type<-factor(sum.bw.Prov.summary.long$geo.type,levels=c("par.state","tip.state"),labels=c("Source","Recipient"))

#write this dataframe to recreate in 
write.csv(sum.bw.Prov.summary.long,"DF/sum.bw.Prov.summary.long.csv")

# sum.bw.Prov.summary.long$geo[which(!sum.bw.Prov.summary.long$geo %in% names(globalPalette.ch))]

PZ3<- ggplot(sum.bw.Prov.summary.long,
       aes(x = geo.type, stratum = geo, alluvium = subject,
           y = mean.n,
           fill = geo, label = geo)) +
  scale_x_discrete(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_flow(alpha = .6,width=0.6) +
  geom_stratum(alpha = .8,width=0.6) +
  geom_text(stat = "stratum", size = 3.2,min.y=0,fontface="bold") +
  pubTheme+
  theme(legend.position = "none", axis.line = element_blank(), text =element_text(size=14),
        axis.ticks.x = element_blank(),axis.text.x = element_text(hjust=c(0.5,0.6)))+
  labs(x=NULL,y="Between-province transmission events")+
  GlobFillScale

PZ3

ggsave( "results/Sankey_Betweenprovince.png",height=4,width=5,units="in")


##numerical summary of b/w prov
tot.Par<-sum.bw.Prov.summary %>% dplyr::group_by(par.state) %>%
  dplyr::summarise(.groups = "rowwise", 
                   totalImports=sum(mean.n,na.rm=T),
                   lowerImports=sum(lower,na.rm=T),
                   upperImports=sum(upper,na.rm=T)) %>%
  as.data.frame()
tot.Par

tot.tip<-sum.bw.Prov.summary %>% dplyr::group_by(tip.state) %>%
  dplyr::summarise(.groups = "rowwise", 
                   totalImports=sum(mean.n,na.rm=T),
                   lowerImports=sum(lower,na.rm=T),
                   upperImports=sum(upper,na.rm=T)) %>%
  as.data.frame()
# tot.tip
# sum.bw.Prov.summary

# ## diagonal net transmisison b/w provs
# mat<-table(can.can$tip.state,can.can$par.state)
# 
# #melt long
# mat2<-mat %>% melt()
# colnames(mat2)<-c("Recipient","Source","value")
# 
# #add empty rows for Newf as source
# empty<-data.frame(mat2[1:8,])
# empty$Source<-"Newfoundland\nand Labrador"
# empty$value<-0
# mat3<-rbind(mat2,empty)
# 
# #remove within-provinces
# rem<-c()
# for (i in 1:nrow(mat3)){
#   if(mat3$Source[i]==mat3$Recipient[i]){rem<-c(rem,i)}
# }
# length(rem)
# mat4<-mat3[-rem,]
# 
# #order Var1 and Var2
# mat4$Recipient<-factor(mat4$Recipient, levels = provs)
# mat4$Source<-factor(mat4$Source, levels = provs)
# 
# #Calculate sums
# mat4
# 
# ## plot as a heatmap 
# 
# #example
# # mat<-table(cana.boots[[k]]$state, cana.boots[[k]]$state)
# 
# P1<- mat4 %>% ggplot()+
#   geom_tile(aes(x=Source,y=Recipient,fill=sqrt(value)))+
#   labs(x="Source",y="Recipient",fill="Between-\nprovince\ntransmission\nevents")+
#   theme(axis.text.x=element_text(angle=45,hjust=1))+
#   scale_fill_viridis(discrete="F", breaks = seq(1:7), 
#                         labels = (seq(1:7))^2, limits=c(0.0001,7.5))+
#   # guides(fill=guide_legend(values=seq(0,1.5,0.5),labels=2^(seq(0,1.5,0.5))))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#             panel.background=element_rect("grey95"),
#             text=element_text(size=10,face="bold"),
#             legend.text=element_text(size=6))
# P1
#ggsave("results/Between province transmission events heatmap.png",height=4,width=5,units="in")
# 
# totalSource<-as.data.frame( mat4 %>% group_by(Source) %>%
#                               dplyr::summarise(.groups="rowwise",totalSource=sum(value)))
# totalRecip<-as.data.frame( mat4 %>% group_by(Recipient) %>%
#                               dplyr::summarise(.groups="rowwise",totalRecipient=sum(value)))
# 
# totalz<-cbind(totalSource,totalRecip)
# totalz <- totalz %>% mutate(net=totalSource-totalRecipient)
# P2<-ggplot(totalz)+
#   geom_bar(aes(x=Source, y=net,group=Source,fill=Source),stat="identity")+
#   GlobFillScale+
#   labs(x=NULL, y="Net transmission = total sources - total recipients")+
#   pubThemeDate+
#   theme(legend.position = "none")
# P2
#ggsave("results/NET Between province transmission events.png",height=5,width=4,units="in")

#### SUMMARY FOR SINGLETONS ACROSS BOOTSTRAPS ####
#what percent of singletons came from the top sources
top.sources<-names(rev(sort(table(can.sing.boots[[1]]$par.state))))
for (i in 1:length(top.sources)){
  source<-c()
  for (k in 1:b){
    perc<- signif(length(which(can.sing.boots[[k]]$par.state==top.sources[i]))/
      nrow(can.sing.boots[[k]])*100, digits = 3)
    source<-c(source, perc)
  }
  print(top.sources[i])
  mean.95ci.1(source)
}


#what percent of singletons came from the top sources
top.sources<-names(rev(sort(table(sum.boots[[1]]$Parent.Location))))
for (i in 1:length(top.sources)){
  source<-c()
  for (k in 1:b){
    perc<- signif(length(which(sum.boots[[k]]$Parent.Location==top.sources[i]))/
      nrow(sum.boots[[k]])*100, digits = 3)
    source<-c(source, perc)
  }
  print(top.sources[i])
  mean.95ci.1(source)
}




# Sublineges with 20+ Desc
## Find sublineage to focus on with more than 20 descendants 

sublinLast<-sublin.long2 %>% group_by(Sublineage) %>% summarize(.groups="rowwise",lastDate=first(rev(sort(DatesLSD))), n=n()) %>% as.data.frame()

focus.sublins<-sublinLast[sublinLast$n>20,]


## Look at samps over time just for the focus sublins

sublin.long2 %>% filter(Sublineage %in% focus.sublins$Sublineage) %>%
ggplot(aes(y=Sublineage,color=Lineage))+
  geom_linerangeh(aes(xmin=DateMin,xmax=DateMax),lwd=0.3)+
  geom_linerangeh(aes(xmin=tmrca,xmax=FirstCan),linetype="dotted",lwd=0.3)+ #detection lag  
  geom_point(aes(x=DatesLSD),alpha=0.7,size=0.5)+
  geom_point(aes(x=FirstCan),shape=22,size=0.7)+
  geom_point(aes(x=tmrca),shape=18,size=1)+
  theme_bw()+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background=element_rect("grey95"),
          legend.key.size = unit(0.5,"line"),
          text=element_text(size=10,face="bold"),
          legend.text=element_text(size=8),
        legend.position="none")+
    scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",
                 limits=c(as.Date("2020-01-01"),as.Date("2021-04-01")))+
  labs(x=NULL,y="Sublineages with >=20 sampled desc.")+
  colScale+
  scale_y_discrete(expand = c(0.05,0.05))+
  guides(color=guide_legend(title="Province of Introduction",title.position = "top",ncol=4))+
  geom_text(aes(label=Sublineage,x=as.Date("2021-02-10",y=Sublineage)),vjust=0.5,hjust=0,size=2,color="grey40")
#ggsave("results/SamplesOverTimeFor_FOCUS_Sublineages_lincolor.png",height=7,width=5,units="in")

## Repeat colored by province of intro
sublin.long2 %>% filter(Sublineage %in% focus.sublins$Sublineage) %>%
ggplot(aes(y=Sublineage,color=Node.Location))+
  geom_linerangeh(aes(xmin=DateMin,xmax=DateMax),lwd=0.3)+
  geom_linerangeh(aes(xmin=tmrca,xmax=FirstCan),linetype="dotted",lwd=0.3)+ #detection lag  
  geom_point(aes(x=DatesLSD),alpha=0.7,size=0.5)+
  geom_point(aes(x=FirstCan),shape=22,size=0.7)+
  geom_point(aes(x=tmrca),shape=18,size=1)+
  theme_bw()+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background=element_rect("grey95"),
          legend.key.size = unit(0.5,"line"),
          text=element_text(size=10,face="bold"),
          legend.text=element_text(size=8),
        legend.position="none")+
    scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",
                 limits=c(as.Date("2020-01-01"),as.Date("2021-04-01")))+
  labs(x=NULL,y="Sublineages with >=20 sampled desc.")+
  GlobColScale+
  scale_y_discrete(expand = c(0.05,0.05))+
  guides(color=guide_legend(title="Province of Introduction",title.position = "top",ncol=1))+
  geom_text(aes(label=Sublineage,x=as.Date("2021-02-10",y=Sublineage)),vjust=0.5,hjust=0,size=2,color="grey40")
#ggsave("results/SamplesOverTimeFor_FOCUS_Sublineages_Provcolor.png",height=7,width=5,units="in")



## upper limits of introductions
seq.case<-9657/866487

sing.tot<-c()
for (k in 1:b){
   sing.tot<-c(sing.tot,nrow(can.sing.boots[[k]] ))
}
mean.95ci.X(sing.tot,0) #1342 (1311 - 1374)"

sublin.tot<-c()
for (k in 1:b){
  sublin.tot<-c(sublin.tot,nrow(sum.boots[[k]]))
}
mean.95ci.X(sublin.tot,0) #"411 (402 - 420)"

(1342 + 411 )* 1/seq.case
(1311 + 402 )* 1/seq.case
(1374 + 420 )* 1/seq.case

#upper estimate of introductions
# 157,290 (153,701 - 160,969)


## Look at whether the country of exposure is consistent with the parental origins of sublineage

k=8

allmet<-read.csv(clean.meta)
canmet<-allmet[which(allmet$country=="Canada"),]
nrow(canmet)
table(canmet$country_exposure) 
cantravel<-canmet[which(canmet$country_exposure != "Canada"),]

#which canadians with travel history were a descendant of a sublineage
travdesc<-cantravel[which(cantravel$gisaid_epi_isl %in% sublin.long[[k]]$Descendant.AccessionIDs),]
travdesc$inferredParent<-NA
travdesc$sublineage2<-NA
for (i in 1:nrow(travdesc)){
  matchy<-which(sublin.long[[k]]$Descendant.AccessionIDs== travdesc$gisaid_epi_isl[i])
  travdesc$sublineage2[i]<-sublin.long[[k]]$Sublineage2[matchy]
  travdesc$inferredParent[i]<-sum.boots[[k]]$Parent.Location[which(sum.boots[[k]]$Sublineage2==sublin.long[[k]]$Sublineage2[matchy])]
}

travdesc$inferredParent
travdesc$country_exposure
travdesc$sublineage2
#100% agreement

#what about the singletons?



