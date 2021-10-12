#!/usr/bin/env Rscript
# Make geography colors and lineage colors

#usage:
# $Rscript 

#setup libs
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(stringi)
# library(devtools)
# devtools::install_github("ropenscilabs/ochRe")
library(ochRe)
# library(cowplot)
library(grid)
library(gridExtra) 
# library(ggplotify)
library(gtools)

## Setup inputs ##
#sublineages input
sum.in<-list.files("DF",pattern="sublin",full.names = T)
sum.in<-mixedsort(sum.in)
##boots
b<-length(sum.in)
BOOTS<-1:b
#meta in
meta.in<-"../00_cleanData/cleaned/clean_fake_meta.csv"
meta.b.in<-paste("DF/meta.b",BOOTS,".csv",sep="")

#### Make global and provincial color scheme ####
## read in clean meta (no subsamp)
meta.b1<-read.csv(meta.in)

### find locations that are the highest contributors to importations####
#group the smaller ones in analysis/viz
#summarize the mean number of imports from each location in each sum.b* file
sumz<-replicate(b,vector())
sumz.t<-replicate(b,vector())
sum.b1<-read.csv(sum.in[1])

for (k in 1:b){
  sumz[[k]]<-read.csv(sum.in[k])
  sumz.t[[k]]<-as.data.frame(table(sumz[[k]]$Parent.Location))
  sumz.t[[k]]$SampleSet<-k
  if (k>1){sum.b1<-rbind(sum.b1,sumz[[k]])} #useful below
}

#rbind the summaries of number of parental imports
sumz.all<-rbind(sumz.t[[1]],sumz.t[[2]])
for (k in 1:b){
  sumz.all<-rbind(sumz.all,sumz.t[[k]])
}

#now take mean for each location
sumzImports<-sumz.all %>% dplyr::group_by (Var1) %>% 
  dplyr::summarise(.groups="rowwise",meanImport=round(mean(Freq),digits=0)) %>%
  as.data.frame()
sumzImports<-sumzImports[rev(order(sumzImports$meanImport)),]

#consider geo with 5 or more importations their own entity, otherwise group to continent
big.origins<-sumzImports$Var1[sumzImports$meanImport>=5] 

#### Name lookups and manips ####
#Rename provinces so Canada not in label
sum.b1$Node.Location<-sapply(sum.b1$Node.Location, function (x) str_replace_all(x, "Canada_", ""))

#make sure no "." instead of space 
sum.b1$Node.Location<-str_replace_all(sum.b1$Node.Location, "\\.", " ")
sum.b1$Parent.Location<-str_replace_all(sum.b1$Parent.Location, "\\.", " ")
sum.b1$Descendant.Location<-str_replace_all(sum.b1$Descendant.Location, "\\.", " ")

## Make a lookup table of location used in ASR, location used for colors (smaller groups)
#build lookup table
lookup.geo<-data.frame(og.loc=c(unique(sum.b1$Node.Location),unique(meta.b1$country)), new.loc=NA)

#for all provs, keep same, except maritimes grouped
provs<-names(table(sum.b1$Node.Location))
marit<-c("Newfoundland and Labrador","Nova Scotia","New Brunswick")
for (i in 1:nrow(lookup.geo)){
  if(lookup.geo$og.loc[i] %in% provs){
    ifelse(lookup.geo$og.loc[i] %in% marit, 
           lookup.geo$new.loc[i]<-"Maritimes",
           lookup.geo$new.loc[i]<-lookup.geo$og.loc[i])
    next
  }
  ifelse (lookup.geo$og.loc[i] %in% big.origins,
          lookup.geo$new.loc[i]<-lookup.geo$og.loc[i],
          lookup.geo$new.loc[i]<-meta.b1$region[which(meta.b1$country==lookup.geo$og.loc[i])[1]]
  )
}

lookup.geo
#note in the fake example, so few importations that most countries are not "big contribs"

#Now replace the continent names as needed:
lookup.geo$new.loc<-stri_replace_all_fixed(lookup.geo$new.loc, pattern = c("South America","North America"), replacement = c("Latin America","Latin America"), vectorize_all = FALSE)

#take Canada out 
lookup.geo$new.loc[lookup.geo$og.loc=="Canada"]<-"Canada"

#write the lookuptable for future use
# lookup.geo
write.csv(lookup.geo,"DF/lookup.geo.csv",row.names=F)

## APPLY the lookup table to change names
for (i in 1:nrow(lookup.geo)){
  if(lookup.geo$new.loc[i] != lookup.geo$og.loc[i]){ #if they don't match, go replace instances
    #find matches in node loc, and replace with new loc
    pr.ma<-which(sum.b1$Node.Location==lookup.geo$og.loc[i])
    if (length(pr.ma)>0) {sum.b1$Node.Location[pr.ma]<-lookup.geo$new.loc[i]; next}
    par.ma<-which(sum.b1$Parent.Location==lookup.geo$og.loc[i])
    if (length(par.ma)>0) {sum.b1$Parent.Location[par.ma]<-lookup.geo$new.loc[i]; next}
  }
}

#first pick out divergent color scheme for provinces
provcols<-ochre_palettes$namatjira_div[c(1,2,3,5,7,8)]
names(provcols)<-c("Ontario","Quebec","Maritimes","Alberta","British Columbia","Manitoba")
plot(1:length(provcols),1:length(provcols),col=provcols,pch=16)

#assign continents grouped colors
asia.cols<-c("#d9f0a3","#78c679","#41ab5d","#238443","#005a32") #greens
names(asia.cols)<-c("Asia","Iran","China","India","Russia")
# plot(1:length(asia.cols),col=asia.cols,pch=16)
eur.cols<-c("#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58")
names(eur.cols)<-c("France","Germany","Italy","Spain","United Kingdom","Europe")
# plot(1:length(eur.cols),col=eur.cols,pch=16)
amer.cols<-c("#8c96c6","#8c6bb1","#6e016b") #purples
names(amer.cols)<-c("Latin America","Brazil","USA")
afr.cols<-c("#f768a1") #pink
names(afr.cols)<-"Africa"
oce.cols<-c("#fbb4b9") #pink
names(oce.cols)<-"Oceania"
#no antartica

#single color to represent canada 
can.cols<-provcols[3]
names(can.cols)<-"Canada"
#Add all together
globalPalette<-c(provcols,asia.cols,eur.cols,amer.cols,afr.cols,oce.cols,can.cols)
# nm<-names(globalPalette)
# globalPalette<-as.data.frame(globalPalette)
# rownames(globalPalette)<-nm

#export a tsv of name and hex color for others to use
write.table(as.data.frame(globalPalette),"DF/globalcolors.tsv",sep="\t",row.names = T)

#example plot
plot(1:length(globalPalette),1:length(globalPalette),col=globalPalette,pch=20)


#### Make lineage colors #####
## for subsampled bootstraps, need to make sure all boots' lineages are included
## read in data
meta.b<-replicate(n=b,vector())
for(k in 1:b){
  meta.b[[k]]<-read.csv(meta.b.in[k])
}

#merge them all and then identify unique lins
meta.b1<-rbind(meta.b[[1]],meta.b[[2]])
for(k in 3:b){
  meta.b1<-rbind(meta.b1,meta.b[[k]])
}

#different colors with meaningful shades
mylinz<-sort(as.character(unique(meta.b1$Lineage )))

## sort out aliases in new table
aliases<-data.frame(lineage=mylinz, alias=NA)
for (i in 1:nrow(aliases)){
  if (aliases$lineage[i]=="A"|aliases$lineage[i]=="B") {aliases$alias[i]<-aliases$lineage[i];next}
  if(str_detect(aliases$lineage[i],"AA\\.")){next}
  if(str_detect(aliases$lineage[i],"AB\\.")){next}
  if(str_detect(aliases$lineage[i],"B\\.")){aliases$alias[i]<-aliases$lineage[i];next}
  if(str_detect(aliases$lineage[i],"A\\.")){aliases$alias[i]<-aliases$lineage[i];next}
}
#how many remaining
sh<-which(is.na(aliases$alias))

#go through manually I guess....can't find this table on pango repo
if(length(sh)>0){
  al.sh<-aliases[sh,]
  al.sh$alias<-as.vector(al.sh$lineage) %>% str_replace_all(c("AA.1" = "B.1.177.15.1", 
                                                              "AA.4" = "B.1.177.15.4", 
                                                              "AC.1" = "B.1.1.405.1",
                                                              "AD.2" = "B.1.1.315.2", 
                                                              "AD.2.1" = "B.1.1.315.2.1", 
                                                              
                                                              "AE.1" = "B.1.1.306.1",
                                                              "AE.2" = "B.1.1.306.2", 
                                                              "AE.4" = "B.1.1.306.4", 
                                                              "AE.6" = "B.1.1.306.6",
                                                              "AE.7" = "B.1.1.306.7", 
                                                              "AE.8" = "B.1.1.306.8", 
                                                              "AG.1" = "B.1.1.297.1",
                                                              "AH.2" = "B.1.1.241.2", 
                                                              "AH.3" = "B.1.1.241.3", 
                                                              "AK.1" = "B.1.1.232.1",
                                                              "AK.2" = "B.1.1.232.2", 
                                                              "AL.1" = "B.1.1.231.1", 
                                                              "AM.1" = "B.1.1.216.1",
                                                              "AM.2" = "B.1.1.216.2", 
                                                              "AM.4" = "B.1.1.216.4", 
                                                              "AP.1" = "B.1.1.70.1",
                                                              "AS.2" = "B.1.1.317.2", 
                                                              
                                                              "C.11" = "B.1.1.1.11", 
                                                              "C.13" = "B.1.1.1.13",
                                                              "C.14" = "B.1.1.1.14",
                                                              "C.16" = "B.1.1.1.16",
                                                              "C.17" = "B.1.1.1.17", 
                                                              "C.18" = "B.1.1.1.18",
                                                              "C.21" = "B.1.1.1.21",
                                                              "C.22" = "B.1.1.1.22",
                                                              "C.23" = "B.1.1.1.23",
                                                              "C.25" = "B.1.1.1.25", 
                                                              "C.26" = "B.1.1.1.26",
                                                              "C.28" = "B.1.1.1.28",
                                                              "C.29" = "B.1.1.1.29", 
                                                              "C.30" = "B.1.1.1.30", 
                                                              "C.31" = "B.1.1.1.31",
                                                              "C.32" = "B.1.1.1.32",
                                                              "C.33" = "B.1.1.1.33",
                                                              "C.34" = "B.1.1.1.34", 
                                                              "C.35" = "B.1.1.1.35",
                                                              "C.36.1" = "B.1.1.1.36.1",
                                                              "C.36" = "B.1.1.1.36",
                                                              "C.3" = "B.1.1.1.3", #ORDER IMPORTANT 
                                                              "C.2" = "B.1.1.1.2", 
                                                              "C.1" = "B.1.1.1.1",
                                                              "C.4" = "B.1.1.1.4",
                                                              "C.5" = "B.1.1.1.5", 
                                                              "C.6" = "B.1.1.1.6",
                                                              "C.8" = "B.1.1.1.8",
                                                              "C.9" = "B.1.1.1.9",
                                                              
                                                              "D.2" = "B.1.1.25.2",
                                                              "D.5" = "B.1.1.25.5",
                                                              
                                                              "K.1" = "B.1.1.277.1",
                                                              "K.2" = "B.1.1.277.2",    
                                                              
                                                              "L.1" = "B.1.1.10.1",
                                                              "L.2" = "B.1.1.10.2",
                                                              "L.3" = "B.1.1.10.3",
                                                              "L.4" = "B.1.1.10.4",     
                                                              
                                                              "M.1" = "B.1.1.294.1",     
                                                              
                                                              "N.1" = "B.1.1.33.1",
                                                              "N.2" = "B.1.1.33.2",
                                                              "N.3" = "B.1.1.33.3",
                                                              "N.4" = "B.1.1.33.4",
                                                              "N.5" = "B.1.1.33.5",               
                                                              "N.6" = "B.1.1.33.6",
                                                              "N.8" = "B.1.1.33.8",
                                                              "N.9" = "B.1.1.33.9",   
                                                              
                                                              "P.1" = "B.1.1.28.1",
                                                              "P.2" = "B.1.1.28.2",    
                                                              
                                                              "R.1" = "B.1.1.316.1",                                             
                                                              "W.1" = "B.1.177.53.1",
                                                              "W.2" = "B.1.177.53.2",
                                                              "W.3" = "B.1.177.53.3",   
                                                              "W.4" = "B.1.177.53.4",      
                                                              
                                                              "Y.1" = "B.1.177.52.1",    
                                                              "Z.1" = "B.1.177.50.1"
  ))                                                  
  
  #join back onto aliases
  aliases[sh,]<-al.sh  
}

## test to see if all lineages have an alias
if(length(which(is.na(aliases$alias)))>0)  warning('not all lineages have been assigned an alias')
any(!is.na(aliases$alias)) #should be tr

## make color scheme for all lineages present in the meta (all 50000)
## CHANGE: make scheme based on the alias names, but then repalce with the original names
mylinz<-aliases$alias

#Red colors for A lineages
id.A<-str_which(mylinz,"A")
# mylinz[id.A]
n.A<-length(id.A)
A<-rev(colorRampPalette(brewer.pal(n=9, "YlOrRd")[3:9])(n.A)) #exclude weak colors

#Purples for B.1.1.X and descendents 
#pretty sure all C.:Z. are B.1.1.X.X)
#C.18=B.1.1.1.18; C.36=B.1.1.1.36; C.36.1=B.1.1.1.36.1; C.8=B.1.1.1.8; L.1=B.1.1.10.1; P2=B.1.1.28.2; R.1=B.1.1.316.1
# id.B11<-which(grepl("B.1.1\\.",mylinz) | grepl("C.",mylinz) | grepl("D.",mylinz)| grepl("E.",mylinz)| grepl("H.",mylinz)| grepl("I.",mylinz) |grepl("K.",mylinz)|grepl("L.",mylinz)|grepl("M.",mylinz) | grepl("N.",mylinz))
id.B11<-which(grepl("B.1.1\\.",mylinz))
id.B11<-sort(c(id.B11, which(mylinz=="B.1.1"))) #add on the grandparent of all :)
# mylinz[id.B11]
b11<-length(id.B11)
B11<-rev(colorRampPalette(brewer.pal(n=9, "Purples")[3:9])(b11)) #exclude weak colors

#Blues for B.1.X and descendents
# id.B1<-str_which(mylinz,"B.1")
id.B1<-which(!grepl("B.1.1\\.",mylinz) & grepl("B.1",mylinz))
id.B1<-id.B1[-which(mylinz[id.B1]=="B.1.1")]
mylinz[id.B1]
remov<-which(mylinz[id.B1] %in% as.character(paste("B.",10:19,sep="")))
if (length(remov)>0){
  id.B1<-sort(id.B1[-remov])
}
# mylinz[id.B1]
b1<-length(id.B1)
B1<-rev(colorRampPalette(brewer.pal(n=9, "PuBu")[3:9])(b1)) #exclude weak colors

#Greens for remaining B descendents (B3,B4,B6,B7,...)
id.B<-str_which(paste0(mylinz[c(id.B1,id.B11,id.A)],collapse="|"),mylinz,negate=T)
# mylinz[id.B]
id.B<-sort(id.B)
b<-length(id.B)
# sort(mylinz[id.B])
B<-rev(colorRampPalette(brewer.pal(n=9, "Greens")[3:9])(b)) #exclude weak colors

## Make an object to lookup lineage group
df.A<-data.frame(lineage=mylinz[id.A],lineagegroup="A*")
df.B<-data.frame(lineage=mylinz[id.B],lineagegroup="B*")
df.B1<-data.frame(lineage=mylinz[id.B1],lineagegroup="B.1*")
df.B11<-data.frame(lineage=mylinz[id.B11],lineagegroup="B.1.1*")
#bandaid because my fake data doesn't have B.1* or B.1.1* lineages
if(length(id.B1)==0){
  df.B1<-data.frame(lineage="NA",lineagegroup="B.1*")
}
if(length(id.B11)==0){
  df.B11<-data.frame(lineage="NA",lineagegroup="B.1.1*")
}
df.lingrp<-rbind(df.A,df.B,df.B1,df.B11)
# head(df.lingrp)

#Reorder for matching to palettes
mylinz<-mylinz[c(id.A,id.B1,id.B11,id.B)]
mycolz<-c(A,B1,B11,B)

#replace my linz with non-alias names
#re-order the aliases dataframe to match my linz and then sub out the columns
mylinz<-aliases[match(mylinz,aliases$alias),'lineage']
names(mycolz)<-as.character(mylinz)

#export a tsv of name and hex color for others to use
write.table(as.data.frame(mycolz),"DF/lineagecolors.tsv",sep="\t",row.names=T)

## add back lineage names via alias to lineage group object
#this is actually alias (long form)
df.lingrp$alias<-df.lingrp$lineage
#leftjoin via alias to add the official lineage names
df.lingrp<-df.lingrp[,-1]
df.lingrp<-left_join(df.lingrp,aliases,"alias")
# tail(df.lingrp)
write.csv(df.lingrp,"DF/lineageGroups.csv",row.names = F)

#write a color object for lineage groups based on above
zA<-"#d73027"
zB<-"#f46d43"
zB1<-"#fdae61"
zB11<-"#fee08b"

mygrpz<-unique(df.lingrp$lineagegroup)
mygrpz.col<-c(zA,zB,zB1,zB11)
names(mygrpz.col)<-as.character(mygrpz)

#export a tsv of name and hex color for others to use
write.table(as.data.frame(mygrpz.col),"DF/lineageGroupColors.tsv",sep="\t",row.names=T)

#Example plots
plot(1:length(mygrpz),1:length(mygrpz),col=mygrpz.col,pch=15)

