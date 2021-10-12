#!/usr/bin/env Rscript

# Summarize the metadata
## Objectives 
# Summarize new cases and cumulative cases over time for global and Canada (with and without rolling 7 day average)
# Summarize for 1) clean, no subsampling 2) all_subsample:
#   number of (cumulative) sequences by collection date over time for global and Canada
#   relationship between number of sequences and number of cases by geography and by month, plot + regression
#   lineages over time by geography

# setup libraries
library(knitr)
library(tidyverse)
library(stringr)
library(coronavirus)
library(ggplot2)
library(RColorBrewer)
library(zoo)
library(ggstance)
library(cowplot)
library(lubridate)
library(MASS)
library(ggrepel)

## inputs 
# inputs from commandline
#Rscript summarize_meta.R <full_clean_meta.csv> <subsampled_meta_boot.csv>
args<-commandArgs(trailingOnly=TRUE)
meta.in<-args[1]
sub.meta.in<-args[2]

# manual
#full metadata  cleaned in 00_CleanData
meta.in<-"../00_cleanData/cleaned/clean_fake_meta.csv"
#subsampled metadata, single bootstrap 
sub.meta.in<-"../01_subSample/bootsamples/subsamp_meta_fake.fasta150_1.csv"
global.colors.in<-"DF/globalcolors.tsv"

#outputs
f.out<-"results"
if(!dir.exists(f.out)){dir.create(f.out)}

## objects that are pulled from further down the pipeline
## Geo lookup tables for grouping
#subsampled with inferred dates
meta.inferred.dates<-"DF/InfDates_meta.b1.csv"
# Apply this at the very last df before plotting
lookup.geo<-read.csv("DF/lookup.geo.csv",header=T)

## Setup themes
#publication themes
pubTheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background=element_rect("grey95"),
            axis.line = element_line(colour = "black"),
            legend.key.size = unit(0.5,"line"),
            text=element_text(size=10,face="bold"),
            legend.text=element_text(size=8))

pubThemeDate<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_rect("grey95"), 
                    axis.line = element_line(colour = "black"),
                    text=element_text(size=10,face="bold"),
                    legend.text=element_text(size=8),
                    axis.text.x=element_text(angle = 45,hjust=1,size=rel(1))) #this is the only difference b/w them

## Read in clean metadata pre-sub-sample

#read in the data
# meta<-read.csv(meta.in)
meta<-read.csv(file = meta.in, header = TRUE)
meta<-meta[,-1] #first column is empty
data.date<-last(sort(meta$date,decreasing = F)) # the day of most recent sample

#remove any with only year
# rem<-c()
# for (i in 1:nrow(meta)){
#   if(length(unlist(strsplit(meta$date[i],"-")))==1) {rem<-c(rem,i)}
# }
# table(meta$country[rem]);length(rem)

#change a column name
# colnames(meta)[which(colnames(meta)=="lineage")]<-"Lineage"

#Remove any with only year...come on!
# rem<-c()
# for (i in 1:nrow(meta)){
#   if(length(unlist(strsplit(meta$date[i],"-")))==1) {rem<-c(rem,i)}
# }
# table(meta$date[rem]);length(rem)

#K bye
# if (length(rem)>0){
#   meta<-meta[-rem,]
# }

#make year-month a new column
meta$month<-as.character(NA)
for (i in 1:nrow(meta)){
  if(length(unlist(strsplit(meta$date[i],"-")))==2) {
    meta$month[i]<-format(as.Date(paste(meta$date[i],"-01",sep="")), "%Y-%m");next
    }
  if(length(unlist(strsplit(meta$date[i],"-")))==3) {meta$month[i]<-format(as.Date(meta$date[i]), "%Y-%m")}
}

#fix issue
meta$division[meta$division=="Canada"]<-"Quebec"

## Separate Canadian metadata
meta.can<-meta[meta$country=="Canada",]
l.can<-nrow(meta.can) 
meta.all<-meta #backup that has everything
# length(meta$new.names)==length(unique(meta$new.names))
provs<-unique(meta.can$division)

## Prepare the coronavirus case data over time for global and canada
#update the dataset
coronavirus <- refresh_coronavirus_jhu()

#add a column in the coronavirus df for month
coronavirus$month<-format(as.Date(coronavirus$date), "%Y-%m")

#interested in country-based here 
coronavirus$country<-NA
for (i in 1:nrow(coronavirus)){
  if(coronavirus$location_type[i]=="country") {coronavirus$country[i]<-coronavirus$location[i]; next}
  #otherwise stringsplit it
  coronavirus$country[i]<-unlist(str_split(coronavirus$location[i],", "))[[2]]
}

coronavirus$country[coronavirus$country=="US"]<-"USA"

#want province where country==canada
coronavirus$prov<-NA
for (i in 1:nrow(coronavirus)){
  if (coronavirus$country[i]=="Canada"){
    coronavirus$prov[i]<-unlist(strsplit(coronavirus$location[i],split=", "))[1]
  }
}


## Ensure case and meta can be matched by country, then remove countries with no seqs

#Make exception for those with country_exposure regarding travel history to italy or iran (as these are underrepresented in dataset)
#also diamond princess 
for (i in 1:nrow(meta)){
  if (str_detect(meta$country_exposure[i], "Iran")) meta$country[i]<-"Iran"
  if (str_detect(meta$country_exposure[i], "Diamond Princess")) meta$country[i]<-"Diamond Princess"
  if (str_detect(meta$country_exposure[i], "Italy")) meta$country[i]<-"Italy"
}
#make some changes to meta country
meta$country[which(meta$country=="Republic of the Congo")]<-"Democratic Republic of the Congo"
meta$country[which(meta$country=="Hong Kong")]<-"China"
#check
# table(meta$country) #looks good
meta.country<-unique(meta$country)
case.country<-unique(as.character(coronavirus$country))

# make changes to align with meta, as these are represented by sequences
coronavirus$country<-as.character(coronavirus$country) #change to character for all
coronavirus$country[which(coronavirus$country=="US")]<-"USA" #US to USA
coronavirus$country[str_which(coronavirus$country,"Congo")]<-"Democratic Republic of the Congo"
coronavirus$country[str_which(coronavirus$country,"Taiwan*")]<-"Taiwan"
coronavirus$country[str_which(coronavirus$country,"Czechia")]<-"Czech Republic"

#for this purpose, remove any countries that have no sequences available
#remove sequences from countries with no case contrib
unique(coronavirus[which(!coronavirus$country %in% meta$country),'country'])
unique(meta$country[which(!meta$country %in% coronavirus$country)])

coronavirus<-coronavirus[coronavirus$country %in% meta$country,]
meta<-meta[which(meta$country %in% coronavirus$country),]



## Make a color scheme for all countries with seqs

#Use new color scheme with fewer groupings
#import a tsv of name and hex color for LOCATIONS
globalPalette<-read.table(global.colors.in,sep="\t")

## make color scheme
glob.colz<-row.names(globalPalette)
globalPalette.ch<-as.character(globalPalette$globalPalette)
names(globalPalette.ch)<-glob.colz
GlobColScale<-scale_colour_manual(name = "Location",values = globalPalette.ch,na.value="grey60")
GlobFillScale<-scale_fill_manual(name = "Location",values = globalPalette.ch,na.value="grey60")

# #make a second color scale that is keeping the original names, but assigning them the new colors using the lookup table
# globalPalette.big.nm<-as.vector(c(unique(meta$country), unique(meta$division[meta$country=="Canada"])))
# globalPalette.big.ch<-c()
# for (i in 1:length(globalPalette.big.nm)){
#   if (globalPalette.big.nm[i] %in% lookup.geo$new.loc){ #if it matches new names
#     col<-globalPalette.ch[which(names(globalPalette.ch)==globalPalette.big.nm[i])]
#   }
#   if (!globalPalette.big.nm[i] %in% lookup.geo$new.loc){ #if it doesn't match the new locations
#     col<-globalPalette.ch[which(names(globalPalette.ch)== lookup.geo$new.loc [which(lookup.geo$og.loc == globalPalette.big.nm[i] )] )]
#   }
#   globalPalette.big.ch<-c(globalPalette.big.ch,as.character(col))
# }
# # length(globalPalette.big.ch)==length(globalPalette.big.nm)
# names( globalPalette.big.ch)<-globalPalette.big.nm
# GlobBigColScale<-scale_colour_manual(name = "Location",values = globalPalette.big.ch,na.value="grey60")
# GlobBigFillScale<-scale_fill_manual(name = "Location",values = globalPalette.big.ch,na.value="grey60")


## Can cases over time
#sum values by date and geography
can.cases<-coronavirus %>%
  filter(country=="Canada") %>%
  filter(data_type=="cases_new") %>%
  group_by(date, prov) %>%
  dplyr::summarize(.groups="rowwise", total=sum(value))

#calculate rolling 7-day average of new cases
can.cases <- can.cases %>%
    dplyr::arrange(desc(prov)) %>% 
    dplyr::group_by(prov) %>% 
    dplyr::mutate(cases_7d = zoo::rollmean(total, k = 7, fill = NA)) %>% 
  dplyr::ungroup()

#remove cruise ships from Can provinces
if (length(str_which(can.cases$prov,"Princess"))>0){
  can.cases<-can.cases[-which(can.cases$prov=="Diamond Princess" | can.cases$prov=="Grand Princess"),]
}



## Global cases over time 
#sum values by date and geography
global.cases<-coronavirus %>% 
  filter(data_type=="cases_new") %>%
  group_by(date, country) %>%
  dplyr::summarise(.groups="rowwise",total = sum(value)) %>%
  arrange(date)

#calculate rolling 7-day average of new cases
global.cases <- global.cases %>%
    dplyr::arrange(desc(country)) %>% 
    dplyr::group_by(country) %>% 
    dplyr::mutate(cases_7d = zoo::rollmean(total, k = 7, fill = NA)) %>% 
  dplyr::ungroup()


## Calculate and visualize sequences available over time by country/province

#sum values by date and geography
can.seqs<-meta.can %>%
  group_by(date, division) %>%
  dplyr::summarize(.groups="rowwise", total=n())

global.seqs<-meta %>%
  group_by(date, country) %>%
  dplyr::summarize(.groups="rowwise", total=n())

#calculate rolling 7-day average of new seqs
global.seqs <- global.seqs %>%
    dplyr::arrange(desc(country)) %>% 
    dplyr::group_by(country) %>% 
    dplyr::mutate(seqs_7d = zoo::rollmean(total, k = 7, fill = NA)) %>%
  dplyr::ungroup()

can.seqs <- can.seqs %>%
    dplyr::arrange(desc(division)) %>% 
    dplyr::group_by(division) %>% 
    dplyr::mutate(seqs_7d = zoo::rollmean(total, k = 7, fill = NA)) %>%
  dplyr::ungroup()


## Calculate countries' contribution to global new cases by month

#up to what year-month?
last.mon<-format(as.Date(data.date), "%m")
last.mon.l<-12+as.numeric(last.mon) #starts at 0, so actually leads to 13+..

##MADE CHANGES TO ACCOUNT FOR repeat months (now that we're in 2021...)
mons <- format(ymd(as.Date("2019-12-01")) %m+% months(0:last.mon.l),"%Y-%m")
l.mons<-length(mons)

##Make a list of dataframes where each df has the cases by country in a given month
monthly.cases<-replicate(data.frame(),n=l.mons)
names(monthly.cases)<-mons


#for each month, filter the coronavirus cases to that month only, sum as above and sort in descending order
for (i in 1:l.mons){
  monthly.cases[[i]]<-as.data.frame(coronavirus %>%
                          # filter(country != "Canada") %>% #uncomment if want to exclude canada from these proportions
                          filter(month == mons[i]) %>% #limit to within each month
                          filter(data_type == "cases_new") %>% # cases
                          group_by(country) %>%
                          dplyr::summarise(.groups="rowwise",total = sum(value)) %>%
                          arrange(-total))
  monthly.cases[[i]]$country<-reorder(monthly.cases[[i]]$country, -monthly.cases[[i]]$total)
}

#check for negatives and turn into zero to avoid issues
for (i in 1:l.mons){
  zeros<-which(monthly.cases[[i]]$total<0)
  if(length(zeros)>0) {monthly.cases[[i]]$total[zeros]<-0}
}

#Now for each of these dataframes, make a new column for the proportion of cases in a given month in a given country
all.counts<-c()
for (i in 1:l.mons){
  all.d<-sum(monthly.cases[[i]]$total)
  monthly.cases[[i]]$proportion<-monthly.cases[[i]]$total / all.d
  all.counts<-c(all.counts,all.d)
}

#no data for dec, this R package doesn't include december counts
monthly.cases[[1]]<-monthly.cases[[2]] 
monthly.cases[[1]]$total<-0
monthly.cases[[1]]$proportion<-0
monthly.cases[[1]]$total[monthly.cases[[1]]$country=="China"]<-20
monthly.cases[[1]]$proportion[monthly.cases[[1]]$country=="China"]<-1
#what is each month's contribution to the overall cases
# all.counts/sum(all.counts)

# if we wanted to sample based on this, we would end up massively favouring the present...
# not great for ancestral reconstruction

#add a column for month to pull in the function
for (i in 1:length(monthly.cases)){
  monthly.cases[[i]]<-cbind( monthly.cases[[i]], data.frame(month=rep(mons[i],times=nrow(monthly.cases[[i]]))))
}


## Merge list of monthly cases into a dataframe

#change case counts into dataframe to link it up by month, country
for (i in 1:length(monthly.cases)){
  monthly.cases[[i]]$month<-names(monthly.cases)[[i]] # add a column for month
}

total.cases<-rbind(monthly.cases[[1]],monthly.cases[[2]])
for (i in 3:length(monthly.cases)){
  total.cases<-rbind(total.cases, monthly.cases[[i]])
}


## APPLY THE LOOKUP TABLE TO CHANGE COUNTRY TO REFLECT LARGER GROUPINGS

#for total cases
total.cases$new.geo<-NA
total.cases$country<-as.character(total.cases$country)
for (i in 1:nrow(total.cases)){
    #if it is found in the new locs, then use it
    if(total.cases$country[i] %in% lookup.geo$new.loc){ 
      total.cases$new.geo[i]<-total.cases$country[i]; next
    }
    #if it is NOT found in the new locs, then find it in old locs and use new loc for that
    if(! total.cases$country[i] %in% lookup.geo$new.loc){ 
      match<-which(lookup.geo$og.loc == total.cases$country[i])
      total.cases$new.geo[i]<- lookup.geo$new.loc [match]
    }
}

#check all have a value
which(is.na( total.cases$new.geo))

#replace the columns for ease below
total.cases$country<-total.cases$new.geo
total.cases<-total.cases[,-ncol(total.cases)] #remove the last col

#make total.cases so that it only has one row per location
total.cases<-total.cases %>% dplyr::group_by(country,month) %>%
  dplyr::summarise(.groups="rowwise",total=sum(total),proportion=sum(proportion))

colnames(total.cases)<-c("country","month","n.cases","proportion.cases")


## apply lookup to meta
meta$new.geo<-NA

for (i in 1:nrow(meta)){
    #if it is found in the new locs, then use it
    if(meta$country[i] %in% lookup.geo$new.loc){ 
      meta$new.geo[i]<-meta$country[i]; next
    }
    #if it is NOT found in the new locs, then find it in old locs and use new loc for that
    if(! meta$country[i] %in% lookup.geo$new.loc){ 
      match<-which(lookup.geo$og.loc == meta$country[i])
      meta$new.geo[i]<- lookup.geo$new.loc [match]
    }
}

#check all have a value
which(is.na( meta$new.geo))
meta$country<-meta$new.geo



## Calculate number and proportion of sequences contributed to whole by country and by month. Join to cases by country, month. Join to metaseq

#Need to calculate number of seqs in each country in each month
global.seqs.mo.sub<-meta %>% 
  dplyr::group_by(month,country) %>% 
  dplyr::summarize(.groups="rowwise",n = n())
global.seqs.mo.sub<-as.data.frame(global.seqs.mo.sub)

# global.seqs.mo.sub$country [which(!global.seqs.mo.sub$country %in% total.cases$country)]

#add a column for proportion of cases that were in that country compared to globe
global.seqs.mo.sub$prop<-NA
for (i in 1:nrow(global.seqs.mo.sub)){
  global.seqs.mo.sub$prop[i]<-global.seqs.mo.sub$n[i]/sum(global.seqs.mo.sub$n [which(global.seqs.mo.sub$month==global.seqs.mo.sub$month[i])])
}

#change some column names to facilitate merge
colnames(global.seqs.mo.sub)[3:4]<-c("n.sequence","proportion.sequence")

#do the join
class(total.cases$month)==class(global.seqs.mo.sub$month)

global.cases.seqs<-left_join(total.cases,global.seqs.mo.sub,by=c("country","month") )
# head(global.cases.seqs)

#calculate the net.difference and fold.difference between cases and sequences by month and country
#and overrep=proportion seqs/proportion cases (ie how over-represeneted (>1) or under-represented (<1)
global.cases.seqs <- global.cases.seqs %>% mutate(net.difference=n.cases-n.sequence, fold.difference=n.sequence/n.cases, overrep=proportion.sequence/proportion.cases)

#put zeroes where nec
# table(global.cases.seqs$proportion.sequence)

#make a version with no zeroes (to log)
global.cases.seqs.pl<-global.cases.seqs %>% filter(n.sequence>0 & n.cases>0)


## sequence representation by country

#fit a line
ft<-lm(log(n.sequence) ~ 0 + log(n.cases), data=global.cases.seqs.pl) #force y-intercept=0
R2<-round(as.numeric(summary(ft)[8]),digits = 2) # adjusted R2
ave.sp<-round(mean(global.cases.seqs.pl$n.sequence / global.cases.seqs.pl$n.cases,na.rm=T),digits=2)
Pears1<-cor((global.cases.seqs.pl$n.sequence), (global.cases.seqs.pl$n.cases), method="pearson")
Spear1<-cor((global.cases.seqs.pl$n.sequence), (global.cases.seqs.pl$n.cases), method="spearman")
Pears1;Spear1

#repeat with color by country and multiple months shown
global.cases.seqs.pl %>%
  ggplot()+
  geom_point(aes(x=log(n.cases), y=log(n.sequence), color=country))+
  geom_smooth(aes(x=log(n.cases), y=log(n.sequence)),method="lm",se=T,formula=y~x-1,lwd=0.5,alpha=0.6)+
  annotate(geom="text",label=paste("Adj. R^2 = ",R2,sep=""), x=0,y=6.5,hjust=0)+
  annotate(geom="text",label=paste("Ave. sequence representation = ",ave.sp,sep=""), x=0,y=6.0,hjust=0)+
  annotate(geom="text",label=paste("Spearman rank corr. coef. = ",round(Spear1,digits=2),sep=""), x=0,y=5.5,hjust=0)+
  labs(x="log(# monthly new diagnoses)",y="log(# monthly new sequences)")+
  theme_bw()+
  pubTheme+
  theme(legend.position="none")+
  geom_text(data=subset(global.cases.seqs.pl, log(n.sequence)/log(n.cases)>0.4),
            aes(x=log(n.cases), y=log(n.sequence), label=country),size=1,
            position=position_dodgev(height=1,preserve="total"), hjust=1,vjust=-1)+
  GlobColScale
  # lims(x=c(0,15.5),y=c(0,10.5))
ggsave("results/Global sequence representation before subsampling with labels.png",width=5,height=4,units="in")



## summarize sequence representation overall for all countries before subsampling
global.cases.seqs.over<-as.data.frame(global.cases.seqs.pl %>% dplyr::group_by(country) %>% dplyr::summarize(.groups="rowwise",n.cases=sum(n.cases),n.sequence=sum(n.sequence)))

global.cases.seqs.over<-global.cases.seqs.over%>% mutate(sampling.proportion=n.sequence/n.cases)


ft<-lm(log(n.sequence) ~ log(n.cases), data=global.cases.seqs.over) #force y-intercept=0
R2<-round(as.numeric(summary(ft)[8]),digits = 2) # adjusted R2
ave.sp.pp<-round(mean(global.cases.seqs.over$n.sequence / global.cases.seqs.over$n.cases,na.rm=T),digits=3)
pears2<-cor(global.cases.seqs.over$n.sequence, global.cases.seqs.over$n.cases, method="pearson")
spear2<-cor(global.cases.seqs.over$n.sequence, global.cases.seqs.over$n.cases, method="spearman")
pears2;spear2

#repeat with color by country
#Medrxiv
pp1<-global.cases.seqs.over %>% 
  ggplot()+
    geom_smooth(aes(x=log(n.cases), y=log(n.sequence)),method="lm",se=T,
                formula=y~x,lwd=0.2,alpha=0.4, fullrange=T)+
  geom_point(aes(x=log(n.cases), y=log(n.sequence), color=country))+
  annotate(geom="text",label=paste("Ave. # seq./# case = ",
                                   ave.sp.pp,sep=""), x=10,y=12.6,hjust=0,size=3)+
  annotate(geom="text",label=paste("Pearson's corr. coef = ",
                                   round(pears2,digits=2),sep=""), x=10,y=12,hjust=0,size=3)+
  labs(x="log(total # cases)",y="log(total # sequences before subsampling)")+
  pubTheme+
  theme(legend.position="none")+
  geom_text_repel(aes(x=log(n.cases), y=log(n.sequence), label=country),size=3)+
  GlobColScale+
  coord_cartesian(clip="off")
  # scale_x_continuous(expand=c(0.01,0), limits=c(10,18),breaks=seq(10,18,1))+
  # scale_y_continuous(expand=c(0.01,0), limits=c(3,13),breaks=seq(3,13,1))
pp1
ggsave("results/Global sequence representation before subsampling all months with labels.png",width=5,height=4,units="in")



## stacked density plot of country contributions

## Repeat but with bins for months
global.cases.seqs.df<-global.cases.seqs

#change month into a numeric date variable to smoothen
global.cases.seqs.df <-global.cases.seqs.df %>% mutate(date.smooth=as.Date(paste(month,"-01",sep="")))
#change NA to zeros ### IMPORTANT#
global.cases.seqs.df$proportion.sequence[which(is.na(global.cases.seqs.df$proportion.sequence))]<-0
global.cases.seqs.df$proportion.cases[which(is.na(global.cases.seqs.df$proportion.cases))]<-0
global.cases.seqs.df$n.sequence[which(is.na(global.cases.seqs.df$n.sequence))]<-0
global.cases.seqs.df$n.cases[which(is.na(global.cases.seqs.df$n.cases))]<-0
## Ordering of new geographies
country.ord<-c("Asia","Iran","China","India","Russia","France","Germany","Italy","Spain","United Kingdom","Europe","Latin America","Brazil","USA","Africa","Oceania")
global.cases.seqs.df$country<-factor(global.cases.seqs.df$country,levels = country.ord)



#PUBLICATION MEDRXIV FIGURE 1 part 1 and 2
P1<-global.cases.seqs.df %>% filter (country != "Canada") %>%
  ggplot( aes(x=date.smooth,y=(n.cases),group=country, fill=country))+
  geom_density(position="stack",stat="identity",lwd=0,alpha=0.9)+
  guides(fill=guide_legend(ncol=2,title.position = "top",keywidth = 1,keyheight = 1, title="Region"))+
  GlobFillScale+
  labs(y="Total monthly cases",x=NULL)+
  pubThemeDate+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits=as.Date(c("2020-01-01","2021-01-01")),expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,2e7,5e6),labels=c("0","5M","10M","15M","20M"),expand=c(0,0))+
  theme(panel.background=element_rect("grey95"),
        legend.position = c(0.05,0.95),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5),
        legend.justification=c(0,1),
        legend.background = element_blank())

P1
ggsave("results/TotalGlobalCASEContributionsByCountryAndMonth_densityStack.png",width=5,height=5,units="in")


#check proportions
# global.cases.seqs %>% dplyr::group_by(month) %>% dplyr::summarize(.groups="rowwise",sumseqprop=sum(proportion.sequence,na.rm=T),sumcaseprop=sum(proportion.cases))

##Repeat with sequence representation before sampling
#Plot it
P2<-global.cases.seqs.df %>% filter (country != "Canada") %>%
  ggplot(aes(x=date.smooth,y=n.sequence,group=country, fill=country))+
  geom_density(position="stack",stat="identity",lwd=0,alpha=0.9)+
  GlobFillScale+
  labs(y="Total sequences available",x=NULL)+
  pubThemeDate+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits=as.Date(c("2020-01-01","2021-01-01")),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), breaks=seq(0,70000,10000),labels=c("0","10K","20K","30K","40K","50K","60K","70K"))+
    theme(legend.position = "none",panel.background=element_rect("grey95"))
P2
#ggsave("results/TotalGlobalSEQUENCEContributionsByCountryAndMonth_BEFORE_densityStack.png",width=5,height=5,units="in")


## Calcuate probability of selecting Canadian sequences based on provinces' proportional contribution to new cases in Canada
## Canadian case proportions
## how many are we working with
# nrow(meta.can)
# table(meta.can$division, meta.can$Lineage)
# table(meta.can$Location, meta.can$month)

#fix meta.can and create prov.seqs

# LOOKUP: for meta.can
meta.can$new.geo<-NA
meta.can$division<-str_replace_all(meta.can$division,"Newfoundland\nand Labrador","Newfoundland and Labrador")
for (i in 1:nrow(meta.can)){
    #if it is found in the new locs, then use it
    if(meta.can$division[i] %in% lookup.geo$new.loc){ 
      meta.can$new.geo[i]<-meta.can$division[i]; next
    }
    #if it is NOT found in the new locs, then find it in old locs and use new loc for that
    if(! meta.can$division[i] %in% lookup.geo$new.loc){ 
      match<-which(lookup.geo$og.loc == meta.can$division[i])
      meta.can$new.geo[i]<- lookup.geo$new.loc [match]
    }
}

#check all have a value
which(is.na( meta.can$new.geo))
meta.can$division<-meta.can$new.geo

#Count total and proportion of sequences
#Need to calculate number of seqs in each prov in each month
prov.seqs<-as.data.frame(meta.can %>% #Use canadian meta
  dplyr::group_by(month,division) %>% 
  dplyr::summarize(.groups = "rowwise", n.sequence = n()))

#add a column for proportion of seqs that were in that prov compared to country
prov.seqs$proportion.sequence<-NA
for (i in 1:nrow(prov.seqs)){
  prov.seqs$proportion.sequence[i]<-prov.seqs$n[i]/sum(prov.seqs$n [which(prov.seqs$month==prov.seqs$month[i])])
}

## Read in the provincial summary data
prov<-read.csv("../01_subSample/Canada_cases/20210211_covid19_casesbyprov.csv", header=T)
prov<-prov[-1,]
prov<-prov[-which(prov$prname=="Canada" |prov$prname=="Repatriated travellers"),] 
#unfortunately we don't know where people repatriated to...
prov$date<-as.Date(prov$date, "%Y-%m-%d") 

# prov$prname<-str_replace_all(prov$prname,"Newfoundland and Labrador","Newfoundland\nand Labrador")

#add a column for month with year
prov$month<-format(as.Date(prov$date), "%Y-%m")

#group by province and month
monthly.cases.prov<-as.data.frame(prov %>% dplyr::group_by(prname, month) %>% 
                                    dplyr::summarize(.groups="rowwise",
                                                     n.cases=sum(numtoday)))
monthly.cases.prov2<-monthly.cases.prov
#use the dup later

#remove provinces with no sequences
### COULD COME BACK AND CHANGE THIS BUT COMPLICATES colors etc ####
# meta.can$division<-str_replace_all(meta.can$division,"Newfoundland and Labrador","Newfoundland\nand Labrador")
monthly.cases.prov<-monthly.cases.prov[which(monthly.cases.prov$prname %in% unique(can.seqs$division)),]

#calculate the proportion of cases in a given month from each province  
monthly.cases.prov$proportion.cases<-NA                                                            
for (i in 1:nrow(monthly.cases.prov)){
  monthly.cases.prov$proportion.cases[i]<-monthly.cases.prov$n.cases[i]/sum(monthly.cases.prov$n.cases [which(monthly.cases.prov$month==monthly.cases.prov$month[i])])
}
            
##LOOKUP TABLE
colnames(monthly.cases.prov)[1]<-"division"
# monthly.cases.prov$division<-str_replace_all(monthly.cases.prov$division,"Newfoundland\nand Labrador","Newfoundland and Labrador")

monthly.cases.prov$new.geo<-NA
monthly.cases.prov$division<-as.character(monthly.cases.prov$division)
for (i in 1:nrow(monthly.cases.prov)){
    #if it is found in the new locs, then use it
    if(monthly.cases.prov$division[i] %in% lookup.geo$new.loc){ 
      monthly.cases.prov$new.geo[i]<-monthly.cases.prov$division[i]; next
    }
    #if it is NOT found in the new locs, then find it in old locs and use new loc for that
    if(! monthly.cases.prov$division[i] %in% lookup.geo$new.loc){ 
      match<-which(lookup.geo$og.loc == monthly.cases.prov$division[i])
      monthly.cases.prov$new.geo[i]<- lookup.geo$new.loc [match]
    }
}

#check all have a value
# which(is.na( monthly.cases.prov$new.geo))
# table(monthly.cases.prov$new.geo)

#replace the columns for ease below
monthly.cases.prov$division<-monthly.cases.prov$new.geo
monthly.cases.prov<-monthly.cases.prov[,-ncol(monthly.cases.prov)] #remove the last col

#make monthly.cases.prov so that it only has one row per location/month
monthly.cases.prov<-monthly.cases.prov %>% dplyr::group_by(division,month) %>%
  dplyr::summarise(.groups="rowwise",total=sum(n.cases),proportion=sum(proportion.cases))

colnames(monthly.cases.prov)<-c("division","month","n.cases","proportion.cases")

                        
# need to put in zero values for early months for all provs
monz<-sort(unique(monthly.cases.prov$month))
provz<-sort(unique(monthly.cases.prov$division))
#make an empty df with zero contribs for provs with no cases
emptydf<-as.data.frame(monthly.cases.prov[1,])
emptydf[1,]<-c("Alberta","2020-01",0,0)
emptydf[2,]<-c("Manitoba","2020-01",0,0)
emptydf[3,]<-c("Maritimes","2020-01",0,0)
emptydf[4,]<-c("Ontario","2020-01",0,0)
emptydf[5,]<-c("Quebec","2020-01",0,0)
emptydf[6,]<-c("Alberta","2020-02",0,0)
emptydf[7,]<-c("Manitoba","2020-02",0,0)
emptydf[8,]<-c("Maritimes","2020-02",0,0)
emptydf[9,]<-c("Quebec","2020-02",0,0)
monthly.cases.prov<-rbind(monthly.cases.prov, emptydf)
monthly.cases.prov<-monthly.cases.prov[with(monthly.cases.prov,order(division)),]
monthly.cases.prov$n.cases<-as.numeric(monthly.cases.prov$n.cases)
monthly.cases.prov$proportion.cases<-as.numeric(monthly.cases.prov$proportion.cases)

#do the join to sequence counts
monthly.cases.prov<-left_join(monthly.cases.prov, prov.seqs, by=c("division","month") )

#turn NAs into zeroes
monthly.cases.prov[is.na(monthly.cases.prov)] <- 0

# #calculate the net.difference and fold.difference between cases and sequences by month and prov
# #how about the proportion seqs/proportion cases (ie how over-represeneted (>1) or under-represented (<1)
monthly.cases.prov <- monthly.cases.prov %>% mutate(net.difference=n.cases-n.sequence,
                                                    fold.difference=n.sequence/n.cases,
                                                    overrep=proportion.sequence/proportion.cases)

#simple n.cases and n.sequences, points
monthly.cases.prov %>% filter(n.sequence>0) %>%
ggplot()+
  geom_smooth(aes(x=log(n.cases), y=log(n.sequence)),method="lm",se=T,formula=y~x,lwd=0.5,alpha=0.6)+
  geom_point(aes(x=log(n.cases), y=log(n.sequence),group=division, color=as.factor(division)))+
  pubTheme+
  GlobColScale+
  labs(x="log(# monthly new diagnoses)",y="log(# monthly new sequences)",color="Province")
# #ggsave("results/CanadaProvinces_NewDiagvsNewSequences.png",width=5,height=3,units="in")



# proportional canadian representation plots}
# Make an overall plot (bin totals across months)
# Add in R2, pval, pearson
ft<-lm(log(n.sequence) ~ 0 + log(n.cases), data=monthly.cases.prov) #force y-intercept=0
ft2<-glm.nb(n.sequence ~ n.cases, data=monthly.cases.prov) 
summary(ft2)
R2<-round(as.numeric(summary(ft)[8]),digits = 2) # adjusted R2
###CHANGE: new way of calucaulting avesp to remove zeroes leading to isseus
ave.sp<-round(mean(as.data.frame(monthly.cases.prov %>% filter(n.sequence>0) %>% filter(n.cases>0) %>% dplyr::group_by(division, month) %>% dplyr::summarize(.groups="rowwise",sp=n.sequence/n.cases))$sp),digits=2)

#overall average (not mean across months)

ave.sp.df<-monthly.cases.prov %>% 
  dplyr::group_by(division) %>%
  dplyr::summarize(.groups="rowwise",tot.seq=sum(as.numeric(n.sequence)),tot.case=sum(as.numeric(n.cases))) %>% 
  as.data.frame() %>% 
  mutate(ave.sp=tot.seq/tot.case)
ave.sp.df
ave.sp<-round(mean(ave.sp.df$ave.sp),digits=2) #proport seq rep
# 
# pears3<-cor(log(ave.sp.df$tot.seq), log(ave.sp.df$tot.case), method="pearson")
# spear3<-cor(log(ave.sp.df$tot.seq), log(ave.sp.df$tot.case), method="spearman")
# pears3;spear3
# 
# P6<-ave.sp.df %>%
#   ggplot()+
# geom_smooth(aes(x=log(tot.case), y=log(tot.seq)),method="glm",se=T,formula=y~x,lwd=0.5,alpha=0.3)+
#   geom_point(aes(x=log(tot.case), y=log(tot.seq),group=division, color=as.factor(division)),size=4)+
#   annotate(geom="text",label=paste("Mean sequence representation = ",ave.sp,sep=""), x=6,y=9,hjust=0)+
#   annotate(geom="text",label=paste("Pearson's correlation coefficient = ",round(pears3,digits=2),sep=""),x=6,y=8.4,hjust=0)+
#   pubTheme+
#   GlobColScale+
#   labs(x="log(total # new diagnoses)",y="log(total # sequences)",color="Province")+
#   theme(legend.position="bottom",legend.text = element_text(size=8))+
#   scale_x_continuous(breaks=c(8,10,12,14),limits=c(7.8,13),expand=c(0,0))+
#   scale_y_continuous(expand=c(0,0))+
#   guides(color=guide_legend(ncol=4,title.position="top",title="Province"))
# P6
# #ggsave("results/CanadaProvinces_NewDiagvsNewSequences_OVERALL.png",width=5,height=4,units="in")






#overall average sequence representation of provinces
overall<-monthly.cases.prov %>% dplyr::group_by(division) %>% 
  dplyr::summarize(total.cases=sum(n.cases),
                   total.sequence=sum(n.sequence),
                   total.samp.prop=sum(n.sequence)/sum(n.cases))
summary(overall$total.samp.prop) #.052
sd(overall$total.samp.prop)
# $proportion.cases[which(is.na(monthly.cases.prov$proportion.cases))]<-0

#DENSITY PROPRTIONAL CONTRIBUTIONS
#publication
#change month into a numeric date variable to smoothen
monthly.cases.prov<-monthly.cases.prov %>% mutate(date.smooth=as.Date(paste(month,"-01",sep="")))
#change NA to zeros ### IMPORTANT#
monthly.cases.prov$proportion.cases[which(is.na(monthly.cases.prov$proportion.cases))]<-0
monthly.cases.prov$proportion.sequence[which(is.na(monthly.cases.prov$proportion.sequence))]<-0

P4<-ggplot(monthly.cases.prov, aes(x=date.smooth,y=n.cases,group=division, fill=division))+
  geom_density(position="stack",stat="identity",lwd=0,alpha=0.9)+
  theme(legend.position = c(0.05,0.95),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5),
        legend.justification=c(0,1),
        legend.background = element_blank())+
  guides(fill=guide_legend(title.position = "top",keywidth = 1,keyheight = 1, title="Province"))+  GlobFillScale+
  labs(y="Total monthly cases",x=NULL)+
  pubThemeDate+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits=as.Date(c("2020-01-01","2021-01-01")),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),breaks=seq(0,200000,50000),labels=c("0","50K","100K","150K","200K"))+
  theme(panel.background=element_rect("grey95"))

    
P4
#ggsave("results/TotalCanadianCASEContributionsByProvinceAndMonth_densityStack.png",width=5,height=5,units="in")

#check proportions
# global.cases.seqs %>% dplyr::group_by(month) %>% dplyr::summarize(.groups="rowwise",sumseqprop=sum(proportion.sequence,na.rm=T),sumcaseprop=sum(proportion.cases))

##Repeat with sequence representation before sampling
P5<-ggplot(monthly.cases.prov, aes(x=date.smooth,y=n.sequence,group=division, fill=division))+
  geom_density(position="stack",stat="identity",lwd=0,alpha=0.9)+
  theme(legend.position = "none")+
  GlobFillScale+
  labs(y="Total sequences available",x=NULL)+
  pubThemeDate+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits=as.Date(c("2020-01-01","2021-01-01")),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),breaks=seq(0,3000,1000),labels=c("0","1K","2K","3K"))+
  theme(legend.position = "none",panel.background=element_rect("grey95"))

P5
#ggsave("results/TotalCanadianSEQUENCEContributionsByProvinceAndMonth_BEFORE_densityStack.png",width=5,height=5,units="in")



# AFTER SUBSAMPLING
## Repeat global rep plots with subsampled data
# Post subsampling sequence representation}
meta.sub<-read.csv(meta.inferred.dates,header = T)

#change a column name
colnames(meta.sub)[which(colnames(meta.sub)=="pangolin_lineage")]<-"Lineage"
#make month a new column
meta.sub$month<-format(as.Date(meta.sub$date), "%Y-%m")

## Separate Canadian meta.subdata
meta.sub.can<-meta.sub[meta.sub$country=="Canada",]

### CHANGES AS ABOVE ###
#fix issue
meta.sub$division[meta.sub$division=="Canada"]<-"Quebec"
#change names to merge
for (i in 1:nrow(meta.sub)){
  if (str_detect(meta.sub$country_exposure[i], "Iran")) meta.sub$country[i]<-"Iran"
  if (str_detect(meta.sub$country_exposure[i], "Diamond Princess")) meta.sub$country[i]<-"Diamond Princess"
  if (str_detect(meta.sub$country_exposure[i], "Italy")) meta.sub$country[i]<-"Italy"
}
#make some changes to meta.sub country
meta.sub$country[which(meta.sub$country=="Republic of the Congo")]<-"Democratic Republic of the Congo"
meta.sub$country[which(meta.sub$country=="Hong Kong")]<-"China"
#check
meta.sub<-meta.sub[which(meta.sub$country %in% coronavirus$country),]

## APPLY LOOKUP to meta.sub
# LOOKUP: for meta.sub
meta.sub$new.geo<-NA
for (i in 1:nrow(meta.sub)){
    #if it is found in the new locs, then use it
    if(meta.sub$country[i] %in% lookup.geo$new.loc){ 
      meta.sub$new.geo[i]<-meta.sub$country[i]; next
    }
    #if it is NOT found in the new locs, then find it in old locs and use new loc for that
    if(! meta.sub$country[i] %in% lookup.geo$new.loc){ 
      match<-which(lookup.geo$og.loc == meta.sub$country[i])
      meta.sub$new.geo[i]<- lookup.geo$new.loc [match]
    }
}

#check all have a value
which(is.na( meta.sub$new.geo))
#replace
meta.sub$country<-meta.sub$new.geo

#####GLOBAL######
global.seqs.mo.sub<-meta.sub %>% 
  dplyr::group_by(month,country) %>% 
  dplyr::summarize(.groups="rowwise",n = n()) %>%
  as.data.frame()

# global.seqs.mo.sub$country [which(!global.seqs.mo.sub$country %in% total.cases$country)]

#add a column for proportion of cases that were in that country compared to globe
global.seqs.mo.sub$prop<-NA
for (i in 1:nrow(global.seqs.mo.sub)){
  global.seqs.mo.sub$prop[i]<-global.seqs.mo.sub$n[i]/sum(global.seqs.mo.sub$n [which(global.seqs.mo.sub$month==global.seqs.mo.sub$month[i])])
}

#change some column names to facilitate merge
colnames(global.seqs.mo.sub)[3:4]<-c("n.sequence","proportion.sequence")

#do the join
#note total cases already had lookup applied
global.cases.seqs.sub<-left_join(total.cases,global.seqs.mo.sub,by=c("country","month") )
# head(global.cases.seqs.sub)

#calculate the net.difference and fold.difference between cases and sequences by month and country
#and overrep=proportion seqs/proportion cases (ie how over-represeneted (>1) or under-represented (<1)
global.cases.seqs.sub <- global.cases.seqs.sub %>% mutate(net.difference=n.cases-n.sequence, fold.difference=n.sequence/n.cases, overrep=proportion.sequence/proportion.cases)

#make a version with no zeroes (to log)
global.cases.seqs.sub.pl<-global.cases.seqs.sub %>% filter(n.sequence>0 & n.cases>0)


## Plot overall sequence representation after subsampling

global.cases.seqs.over.sub<-as.data.frame(global.cases.seqs.sub.pl %>% dplyr::group_by(country) %>% dplyr::summarize(.groups="rowwise",n.cases=sum(n.cases),n.sequence=sum(n.sequence)))

global.cases.seqs.over.sub<-global.cases.seqs.over.sub %>% mutate(sampling.proportion=n.sequence/n.cases)

# global.cases.seqs.over.sub<-global.cases.seqs.over.sub[-which(global.cases.seqs.over.sub$n.sequence<1),]
# global.cases.seqs.over.sub<-global.cases.seqs.over.sub[-which(global.cases.seqs.over.sub$n.cases<1),]


ft<-lm(log(n.sequence) ~ log(n.cases), data=global.cases.seqs.over.sub) #force y-intercept=0
R2<-round(as.numeric(summary(ft)[8]),digits = 2) # adjusted R2
ave.sp<-round(mean(as.data.frame(global.cases.seqs.over.sub%>% filter(n.sequence>0) %>% filter(n.cases>0) %>% dplyr::group_by(country) %>% dplyr::summarize(.groups="rowwise",sp=n.sequence/n.cases))$sp),digits=4)
ave.sp

Pears.over<-cor((global.cases.seqs.over.sub$n.sequence), (global.cases.seqs.over.sub$n.cases), method="pearson")
Spear.over<-cor((global.cases.seqs.over.sub$n.sequence), (global.cases.seqs.over.sub$n.cases), method="spearman")

#repeat with color by country
pp2<- global.cases.seqs.over.sub %>%
  ggplot()+
    geom_smooth(aes(x=log(n.cases), y=log(n.sequence)),method="lm",se=T,
                formula=y~x,lwd=0.2,alpha=0.4, fullrange=T)+
  geom_point(aes(x=log(n.cases), y=log(n.sequence), color=country))+
  annotate(geom="text",label=paste("Ave. # seq./# case = ",ave.sp,sep=""), x=10,y=12.6,hjust=0,size=3)+
  annotate(geom="text",label=paste("Pearson's corr. coef = ",round(Pears.over,digits=2),sep=""), x=10,y=12,hjust=0,size=3)+
  labs(x="log(total # cases)",y="log(total # sequences after subsampling)")+
  theme_bw()+
  pubTheme+
  theme(legend.position="none")+
  geom_text_repel(aes(x=log(n.cases), y=log(n.sequence), label=country),size=3, max.overlaps=20)+
  GlobColScale+
  coord_cartesian(clip="off")+
  scale_x_continuous(expand=c(0.01,0), limits=c(10,18),breaks=seq(10,18,1))+
  scale_y_continuous(expand=c(0.01,0), limits=c(3,13),breaks=seq(3,13,1))

pp2
#ggsave("results/Global sequence representation after subsampling all months with labels.png",width=5,height=4,units="in")



## Join the two plots of beffore and after subsampling

plotz<-plot_grid(pp1,pp2,labels=c("A","B"))
titleplotz <- ggdraw() + 
  draw_label("Before                                                                 After",
             fontface = 'bold',x = 0,hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 140))
#combine the big titels and vertplot
plot_grid(titleplotz, plotz,ncol = 1,rel_heights = c(0.08,  1))
#ggsave("results/BeforeAfterSubsamplingSequenceRep_FigS3.png",width=8,height=4,units="in")




#fit a line
# ft<-lm(log(n.sequence) ~ 0 + log(n.cases), data=global.cases.seqs.sub.pl)
# R2<-round(as.numeric(summary(ft)[8]),digits = 3) # adjusted R2 
# ave.sp<-round(mean(global.cases.seqs$n.sequence / global.cases.seqs$n.cases,na.rm=T),digits=4)
#   
#repeat with color by country
# global.cases.seqs.sub.pl %>%
#   ggplot()+
#   geom_point(aes(x=log(n.cases), y=log(n.sequence), color=country))+
#   geom_smooth(aes(x=log(n.cases), y=log(n.sequence)),method="lm",se=T,formula=y~x-1,lwd=0.5,alpha=0.6)+
#   annotate(geom="text",label=paste("Adj. R^2 = ",R2,sep=""), x=1,y=10,hjust=0)+
#   annotate(geom="text",label=paste("Ave. sequence representation = ",ave.sp,sep=""), x=1,y=9.3,hjust=0)+
#   labs(x="log(# monthly new diagnoses)",y="log(# monthly new sequences)")+
#   theme_bw()+
#   pubTheme+
#   theme(legend.position="none")+
#   geom_text(data=subset(global.cases.seqs.sub.pl, log(n.sequence)/log(n.cases)>0.4),
#             aes(x=log(n.cases), y=log(n.sequence), label=country),size=1,
#             position=position_dodgev(height=1,preserve="total"), hjust=1,vjust=-1)+ 
#   GlobColScale+
#   lims(x=c(0,15.5),y=c(0,10))
# #ggsave("results/Global sequence representation after subsampling with labels.png",width=5,height=4,units="in")
#could add confidence and predictive intervals here, or label if residuals >x

# # summarize this by month in facetted plot
# global.cases.seqs.sub.pl %>%
#   ggplot()+
#   geom_point(aes(x=log(n.cases), y=log(n.sequence), color=country),size=0.8)+
#   geom_smooth(aes(x=log(n.cases), y=log(n.sequence)),method="lm",se=T,alpha=0.6,lwd=0.2)+
#   # annotate(geom="text",label=paste("Adjusted R^2 = ",R2,sep=""), x=4,y=10)+
#   labs(x="log(# monthly new diagnoses)",y="log(# monthly new sequences)")+
#   theme_bw()+
#   pubTheme+
#   theme(legend.position="none")+
#   GlobColScale+
#   geom_text(data=subset(global.cases.seqs.sub.pl, log(n.sequence)/log(n.cases)>0.7 ),
#             aes(x=log(n.cases), y=log(n.sequence), label=country),size=1,
#             position=position_dodgev(height = 3,preserve="total"), hjust=1,vjust=-1)+ 
#   facet_wrap(.~month, ncol=5,scales="free")
##ggsave("results/Global sequence representation after subsampling_facetMonth.png",width=10,height=6,units="in")



###### Plot proportional contribution of sequences after #####
## Repeat but with bins for months

global.cases.seqs.sub.df<-global.cases.seqs.sub

# global.cases.seqs.sub.df[which(global.cases.seqs.sub.df$proportion.cases<0),'proportion.cases']<-0
# global.cases.seqs.sub.df<-global.cases.seqs.sub.df[which(global.cases.seqs.sub.df$country %in% meta$country),]

#change month into a numeric date variable to smoothen
global.cases.seqs.sub.df <-global.cases.seqs.sub.df %>% mutate(date.smooth=as.Date(paste(month,"-01",sep="")))
#change NA to zeros ### IMPORTANT#
global.cases.seqs.sub.df$proportion.sequence[which(is.na(global.cases.seqs.sub.df$proportion.sequence))]<-0
global.cases.seqs.sub.df$n.sequence[which(is.na(global.cases.seqs.sub.df$n.sequence))]<-0

#PUBLICATION
##Repeat with sequence representation after sampling
global.cases.seqs.sub.df$country<-factor(global.cases.seqs.sub.df$country,levels=country.ord)

P3<-global.cases.seqs.sub.df %>% filter (country != "Canada") %>%
  ggplot(aes(x=date.smooth,y=n.sequence,group=country, fill=country))+
  geom_density(position="stack",stat="identity",lwd=0,alpha=0.9)+
  theme(legend.position = "none")+
  GlobFillScale+
  labs(y="Total sequences subsampled",x=NULL)+
  pubThemeDate+
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits=as.Date(c("2020-01-01","2021-01-01")),expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),breaks=seq(0,3000,1000),labels=c("0","1K","2K","3K"))+    theme(legend.position = "none",panel.background=element_rect("grey95"))

P3
# #ggsave("results/TotalGlobalSEQUENCEContributionsByCountryAndMonth_AFTER_densityStack.png",width=5,height=5,units="in")

# ggplot(data=subset(global.cases.seqs.sub.df,proportion.sequence>0.1),
#        aes(x=date.smooth,y=proportion.sequence,group=country, fill=country))+
#   geom_bar(position="dodge",stat="identity")+
#   geom_text(aes(label=country),position="dodge")+
#   theme(legend.position = "none")+
#   GlobFillScale+
#   labs(y="Prop. global sequences before sub-sampling",x=NULL)+
#   pubThemeDate+
#   scale_x_date(date_breaks = "1 month", date_minor_breaks = "2 weeks", date_labels = "%b %Y",limits=as.Date(c("2019-12-01","2021-03-01")))





###Make another version of P6 including all the global regions and the provinces
ave.sp.df.glob<-global.cases.seqs.sub.df %>% dplyr::group_by(country) %>% 
  filter (country != "Canada") %>% 
  dplyr::summarize(.groups="rowwise",tot.seq=sum(n.sequence,na.rm=T),tot.case=sum(n.cases)) %>% 
  as.data.frame() %>% mutate(ave.sp=tot.seq/tot.case)

colnames(ave.sp.df)[1]<-"country"
ave.sp.df.all<-rbind(ave.sp.df.glob, ave.sp.df)
ave.sp.all<-round(mean(ave.sp.df.all$ave.sp)*100,digits=2) #proport seq rep

all.all<-sum(ave.sp.df.all$tot.seq)/sum(ave.sp.df.all$tot.case)
all.all

pears3<-cor(log(ave.sp.df.all$tot.seq), log(ave.sp.df.all$tot.case), method="pearson")
spear3<-cor(log(ave.sp.df.all$tot.seq), log(ave.sp.df.all$tot.case), method="spearman")
pears3;spear3

P6<-ave.sp.df.all %>%
  ggplot()+
geom_smooth(aes(x=log10(tot.case), y=log10(tot.seq)),method="glm",se=T,formula=y~x,lwd=0.1,alpha=0.3,fullrange=T)+
  geom_point(aes(x=log10(tot.case), y=log10(tot.seq),group=country, color=as.factor(country)),size=2,alpha=0.5)+
  geom_text_repel(aes(x=log10(tot.case), y=log10(tot.seq),group=country, color=as.factor(country),label=country),size=3.5,fontface="bold", force=35,max.overlaps=30)+
  annotate(geom="text",label=paste("Mean ",ave.sp.all," sequences / 100 cases",sep=""), x=3.1,y=4.4,hjust=0)+
  annotate(geom="text",label=paste("Pearson's corr. coef. = ",round(pears3,digits=2),sep=""),x=3.1,y=4.2,hjust=0)+
  pubTheme+
  GlobColScale+
  labs(x="Cumulative cases",y="Cumulative sequences subsampled")+
  theme(legend.position="none")+
  scale_x_continuous(expand=c(0,0),limits=c(3,7.6),breaks=seq(3,7,1),labels=c("1K","10K","100K","1M","10M"))+
  scale_y_continuous(expand=c(0,0),limits=c(1.8,4.5),breaks=seq(2,4,1),labels=c("100","1K","10K"))+
  coord_cartesian()

# position=position_jitter(width=0.1,height=0.1)
P6
# ggsave("results/NAMESNewDiagvsNewSequences_OVERALL.png",width=5,height=4,units="in")




###Make a cowplot for figure 1###
# PUBLICATION FIG 1}
# plot_grid(P1,P2,P3,ncol=3,labels=c("a","b","c"))
# ##ggsave("results/Figure1_contributions.png",width=15,height=5,units="in")

## add in the Canadian province representation plots
# plot_grid(P1,P2,P3,P4,P5,P6,ncol=3,labels=c("A" ,"B", "C", "D", "E", "F"))
##ggsave("results/Figure1_contributions.png",width=15,height=10,units="in")

## add in the Canadian province representation plots
vertplot<-plot_grid(P1,P4,P2,P5,P3,P6,ncol=2,labels=c("A" ,"D", "B", "E", "C", "F"))
vertplot
# ##ggsave("results/MEDRXIV_Figure1_contributions_vert.png",width=10,height=15,units="in")
# ##ggsave("results/MEDRXIV_SMALL_Figure1_contributions_vert.png",width=7,height=10.5,units="in")
##ggsave("results/MEDRXIV_med_Figure1_contributions_vert.png",width=8.5,height=12,units="in")
#P6<-
title1 <- ggdraw() + 
  draw_label("Global                                                                    Canada",
             fontface = 'bold',x = 0,hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 140))
#combine the big titels and vertplot
plot_grid(title1, vertplot,ncol = 1,rel_heights = c(0.02,  1))
ggsave("results/Figure1_Sampling_log10.png",width=8.5,height=10,units="in")



##How many incomplete dates in the subsampled dataset?

s.meta<-read.csv(sub.meta.in)
s.meta<-s.meta[,-1]

s.meta$date<-str_replace_all(s.meta$date,"-XX","")

incomp<-c()
for (i in 1:nrow(s.meta)){
  if(length(unlist(strsplit(s.meta$date[i],split="-")))<3) 
    incomp<-c(incomp,s.meta$gisaid_epi_isl[i])
}
length(incomp) #number of subsampled seqs with incomplete dates to be inferred

#table of number of sequences with incomplete dates included in subsample, by country
#mostly Canada (n=2183), but 240 Japan, 95 form Poland, 79 India,...
rev(sort(table(s.meta[which(s.meta$gisaid_epi_isl %in% incomp),c('division')])))

s.meta0115<-read.csv("../../../20210115_v2_phylogeo/01_subSample/subsamp_all/20210115_subsamp_meta_50000_1.csv")
s.meta0115$date<-str_replace_all(s.meta0115$date,"-XX","")
incomp0115<-c()
for (i in 1:nrow(s.meta0115)){
  if(length(unlist(strsplit(s.meta0115$date[i],split="-")))<3) 
    incomp0115<-c(incomp0115,s.meta0115$gisaid_epi_isl[i])
}
length(incomp0115) #number of subsampled seqs with incomplete dates to be inferred in previous biuld
rev(sort(table(s.meta0115[which(s.meta0115$gisaid_epi_isl %in% incomp0115),c('country')])))
#all Canada, 1359




## print a supplementary table of sequences available and sequences sampled by month
globclean<-data.frame(table(meta$month[meta$country!="Canada"]))
globsamp<-data.frame(table(s.meta$month[s.meta$country!="Canada"]))

seqs.monthly<-data.frame('Var1'=globclean$Var1,'GlobalClean'=globclean$Freq,'GlobalSample'=globsamp$Freq)
cansamp<-data.frame(table(s.meta$month[s.meta$country=="Canada"]))
canclean<-data.frame(table(meta$month[meta$country=="Canada"]))
seqs.monthly<-left_join(seqs.monthly,canclean,by="Var1")
seqs.monthly<-left_join(seqs.monthly,cansamp,by="Var1")
colnames(seqs.monthly)<-c("Month","Clean Global","Sampled Global","Clean Canadian","Sampled Canadian")
seqs.monthly

#why the issue with 2020-02 and 2020-03, where there are fewer Canadian sampled sequences than expected, yet, for 02, there are 5 extra global sequences from Canada...
table(meta$country[meta$month=="2020-02"])
table(s.meta$country[s.meta$month=="2020-02"])
#write.csv(seqs.monthly,"results/SequencesPrePostSample.csv")


#make a supplementary table for canadian cases by month and province

head(monthly.cases.prov2)
#go from long to wide (column for each province)
monthly.cases.prov3<-monthly.cases.prov2 %>% spread(key=prname,value=n.cases)
#write.csv(monthly.cases.prov3,"results/CasesByMonthAndProvince.csv")


## calculate rolling sequence representation by province

can.cases[1500:1800,]
head(can.seqs)
can.seqs2<-can.seqs 
colnames(can.seqs2)[2]<-"prov"
can.seqs2$date<-as.Date(can.seqs2$date)

#join by day and prov
seq.rep<-left_join(can.cases,can.seqs2,by=c("date","prov"))

#rolling y day ave of seqs/rolling 7 day ave of samples
# seq.rep<-seq.rep %>% mutate(seq.rep_7d=seqs_7d/cases_7d)
# seq.rep<-seq.rep %>% mutate(seq.rep_daily=total.y/total.x)

seq.rep<-seq.rep%>%
  dplyr::group_by(prov)%>%
  dplyr::mutate(cases_7d2=zoo::rollmean(total.y,k=7,fill=NA,align="left"),
                seqs_7d2=zoo::rollmean(total.x,k=7,fill=NA,align="left"),
                seq.rep_7d2=seqs_7d2/cases_7d2) %>%
  dplyr::ungroup()

head(seq.rep)

#change NAs to zeroes
seq.rep$seq.rep_7d2[which(is.na(seq.rep$seq.rep_7d2))]<-0

range(seq.rep$seq.rep_7d)

#plot this out,might need to deal with zeroes
seq.rep %>% filter (prov %in% provs) %>%
  ggplot()+
  geom_line(aes(x=date,y=seq.rep_7d2,group=prov,color=prov))+
  labs(y="")


## Model the submission delay

nrow(meta.all)
meta.all$date_submitted<-as.Date(meta.all$date_submitted, format="%Y/%m/%d")

#fill in date with 15th day of month for incomplete for this 
rem<-c()
meta.all$date.inf<-as.Date(NA)
for (i in 1:nrow(meta.all)){
  dt<-unlist(strsplit(meta.all$date[i],split="-"))
  if(length(dt)==3) {meta.all$date.inf[i]<-meta.all$date[i];next}
  if(length(dt)==1) {rem<-c(rem,i); next} 
  if(length(dt)==2) {dt[3]<-15; 
                     meta.all$date.inf[i]<-paste0(dt,collapse="-")}
}
length(rem)
meta.all<-meta.all[-rem,]
meta.all$date.inf<-as.Date(meta.all$date.inf, format="%Y/%m/%d")
class(meta.all$date_submitted)
# meta.all$submission.lag = meta.all$date_submitted - meta.all$date
meta.all<-meta.all %>% mutate(submission.lag = date_submitted - date.inf)

head(meta.all)

#Make a density plot about it
library(ggridges)
# meta.all$submission.lag<-as.numeric(meta.all$submission.lag)

#filter to countries with most contributions
keepers<-names(rev(sort(table(meta.all$country)))[1:10])

#summary stats
lag.summ<-meta.all  %>% filter (country %in% keepers) %>% dplyr::group_by(country) %>%
  dplyr::summarize(.groups="rowwise", 
            med.delay=median(submission.lag), mean.delay=mean(submission.lag)) %>% as.data.frame()
lag.summ<-lag.summ[with(lag.summ, order(med.delay)),]

meta.keeper<-meta.all %>%filter (country %in% keepers) 
meta.keeper$country<-factor(meta.keeper$country, levels=lag.summ$country)
meta.keeper %>%
ggplot(aes(x=submission.lag,y=country,fill=country)) +
  geom_density_ridges(alpha=0.5,lwd=0.3)+
  theme_classic() +
  labs(x="Submission delay (days)",y=NULL)+
  theme(legend.position = "none", text=element_text(size=12,face="bold",color="black"))+
  scale_x_continuous(limits = c(-30, 325))+
  scale_y_discrete(expand = c(0.2,0))+
  geom_text(data=lag.summ, aes(label=med.delay,x=320,y=country),hjust=1,vjust=1,size=4)+
  annotate(geom="text",label="Median submission\ndelay (days)",x=320,y=11,hjust=1,size=4)+
  GlobFillScale
##ggsave("results/ggridge_submissionDelay.png",width=4,height=3.5,units="in")

#plot Canada's submission lag over time

meta.keeper %>% filter (country=="Canada") %>%
  ggplot(aes(x=date_submitted,y=submission.lag,group=country,color=country))+
  geom_point(alpha=0.5,size=0.5)+
  geom_smooth(method="gam",se=F)+
  theme_classic() +
  GlobColScale+
  labs(y="Submission delay (days)",x=NULL)+
  theme(legend.position = "top", text=element_text(size=12,face="bold",color="black"))
##ggsave("results/Canada_submissionDelay_overtime.png",width=4,height=3.5,units="in")

#by province
meta.keeper$division<-str_replace_all(meta.keeper$division,"Newfoundland and Labrador","Newfoundland/n and Labrador")
meta.keeper %>% filter (country=="Canada") %>%
  ggplot(aes(x=date_submitted,y=submission.lag,group=division,color=division))+
  geom_point(alpha=0.5,size=0.5)+
  theme_classic() +
  GlobColScale+
  labs(y="Submission delay (days)",x=NULL)+
  theme(legend.position = "right", text=element_text(size=8,face="bold",color="black"))
##ggsave("results/CanadaProvs_submissionDelay_overtime.png",width=6,height=3.5,units="in")

## Rank countries in terms of sequences, sequences/cases, submission delay

case.tots<-global.cases %>% group_by(country) %>% dplyr::summarize(.groups="rowwise",total.cases=sum(total)) %>% as.data.frame()

meta.all %>% dplyr::group_by(country) %>% dplyr::summarise(.groups="rowwise",
                                                           seq.tots=n(),
                                                           med.delay=median(submission.lag))
