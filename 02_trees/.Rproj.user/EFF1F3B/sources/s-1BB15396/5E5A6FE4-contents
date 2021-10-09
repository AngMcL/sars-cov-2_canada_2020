#!/usr/bin/env Rscript

# this script subsamples the cleaned GISAID data based on specified number of total samples
# sampling probability proportional to each country's contribution to global cases in each month
# outputs a subsampled fasta and metadata for each of b bootstraps

#usage:
# $Rscript subsample_n.R <seq.in> <meta.in> <number of boots> <number of samples to take total>

## load libraries
library(knitr)
library(tidyverse)
library(stringr)
library(ape)
library(coronavirus)
library(Biostrings)
library(lubridate)

#### Inputs and outputs here ####
# inputs from commandline
args<-commandArgs(trailingOnly=TRUE)
# test if all arguments specified
if (length(args)<4) {stop("Four arguments must be supplied", call.=FALSE)} 

align.in<-args[1]
meta.in<-args[2]
b<-as.numeric(args[3])
n.samp<-as.numeric(args[4])

## Or, change inputs manually
# align.in<-"../00_cleanData/masked/mask_clean_fake.fasta"
# meta.in<-"../00_cleanData/cleaned/clean_fake_meta.csv"
# b<-10
# n.samp<-150 #note this is the total (Canadian + global), default is take all Canadian

#outputs
align.base<-unlist(str_split(last(unlist(str_split(last(unlist(str_split(align.in,"/"))),"_"))),"\\."))[[1]]
if(!dir.exists("bootsamples")){dir.create("bootsamples")}
align.out<-paste("bootsamples/subsamp_align_",align.base,sep="")
meta.out<-paste("bootsamples/subsamp_meta_",align.base,sep="")

## Read in and link the alignment and metadata
#read in the data
align<-readDNAStringSet(align.in)
l.align<-length(align)
nmz<-names(align)

#META
meta<-read.csv(meta.in)
meta<-meta[,-1] #first column is empty
data.date<-last(sort(meta$date)) # the day of most recent sample

#change a column name
colnames(meta)[which(colnames(meta)=="pangolin_lineage")]<-"Lineage"

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

#fix issue with improper division
meta$division[meta$division=="Canada"]<-"Quebec"
   
#check that all names in alignment match meta and vice versa
all(nmz %in% meta$new.names)
all(meta$new.names %in% nmz)

#check no duplicate names
length(meta$new.names)==length(unique(meta$new.names))
length(names(align))==length(unique(names(align)))
# length(names(align))


#### Calculate countries' contribution to global new cases by month ####


## Pull country-specific case counts from coronavirus lib
#update the dataset (can be slow)
coronavirus <- refresh_coronavirus_jhu()

#add a column in the coronavirus df for month
coronavirus$month<-format(as.Date(coronavirus$date), "%Y-%m")
# table(coronavirus$month)

#only interested in country-based here 
coronavirus$country<-NA
for (i in 1:nrow(coronavirus)){
  if(coronavirus$location_type[i]=="country") {coronavirus$country[i]<-coronavirus$location[i]; next}
  #otherwise stringsplit it
  coronavirus$country[i]<-unlist(str_split(coronavirus$location[i],", "))[[2]]
}

## Set up list of df for each month 
last.mon<-format(as.Date(data.date), "%m")
last.mon.l<-12+as.numeric(last.mon) 

# mons<-c(0:last.mon) #changed Dec 2019 to 0
mons <- format(ymd(as.Date("2019-12-01")) %m+% months(0:last.mon.l),"%Y-%m")
l.mons<-length(mons)

##Make a list of dataframes where each df has the cases by country in a given month
monthly.cases<-replicate(data.frame(),n=l.mons)
names(monthly.cases)<-mons

#for each month, filter the coronavirus cases to that month only, sum as above and sort in descending order
for (i in 1:l.mons){
  monthly.cases[[i]]<-as.data.frame(coronavirus %>%
                          filter(country != "Canada") %>% #exclude canada from these proportions
                          filter(month == mons[i]) %>% #limit to within each month
                          filter(data_type == "cases_new") %>% # cases
                          group_by(country) %>%
                          dplyr::summarise(.groups="rowwise",total = sum(value)) %>%
                          arrange(-total))
  monthly.cases[[i]]$country<-reorder(monthly.cases[[i]]$country, -monthly.cases[[i]]$total)
}

#for each df, make a new column for the proportion of cases in a given month in a given country, relative to global
#also pull a cumulative number of cases from each month into allcounts
all.counts<-c()
for (i in 1:l.mons){
  all.d<-sum(monthly.cases[[i]]$total)
  monthly.cases[[i]]$proportion<-monthly.cases[[i]]$total / all.d
  all.counts<-c(all.counts,all.d)
}

#no data for dec: we presume was all (almost?) in china
monthly.cases[[1]]<-monthly.cases[[2]] #use its structr
monthly.cases[[1]]$total<-0
monthly.cases[[1]]$proportion<-0
monthly.cases[[1]]$total[monthly.cases[[1]]$country=="China"]<-20 #conservative total
monthly.cases[[1]]$proportion[monthly.cases[[1]]$country=="China"]<-1

## Ensure case and sequence data can be matched by country
#Make exception for those with country_exposure regarding travel history to italy, iran, or diamond princess 
#(as these are underrepresented in early dataset, see Worobey et al 2020)
for (i in 1:nrow(meta)){
  if (str_detect(meta$month[i],"2020")){
    if (str_detect(meta$country_exposure[i], "Iran")) meta$country[i]<-"Iran"
    # if (str_detect(meta$country_exposure[i], "Diamond Princess")) meta$country[i]<-"Diamond Princess"  #no longer in the data
    if (str_detect(meta$country_exposure[i], "Italy")) meta$country[i]<-"Italy"
  }
}

#check
# table(meta$country) #looks good
meta.country<-unique(meta$country)
case.country<-unique(as.character(monthly.cases[[2]]$country)) #same in all df in the list

#check how well they match to the countries in monthly cases 
# meta.country[which(!meta.country %in% case.country)] #mostly places that need to be matched
# case.country[which(!case.country %in% meta.country)] #places we don't have seqs for 

# make changes to align with meta, as these are represented by sequences
for (i in 1:length(monthly.cases)){
  monthly.cases[[i]]$country<-as.character(monthly.cases[[i]]$country) #change to character for all
  monthly.cases[[i]]$country[which(monthly.cases[[i]]$country=="US")]<-"USA" #US to USA
  monthly.cases[[i]]$country[str_which(monthly.cases[[i]]$country,"Congo")]<-"Democratic Republic of the Congo"
  monthly.cases[[i]]$country[str_which(monthly.cases[[i]]$country,"Taiwan*")]<-"Taiwan"
  monthly.cases[[i]]$country[str_which(monthly.cases[[i]]$country,"Czechia")]<-"Czech Republic"
}

#make some changes to meta country
meta$country[which(meta$country=="Republic of the Congo")]<-"Democratic Republic of the Congo"

#check again
# meta.country<-unique(meta$country)
# case.country<-unique(as.character(monthly.cases[[1]]$country))
# meta.country[which(!meta.country %in% case.country)]
# case.country[which(!case.country %in% meta.country)]


#change case counts into dataframe to link it up by month, country
for (i in 1:length(monthly.cases)){
  monthly.cases[[i]]$month<-names(monthly.cases)[[i]] # add a column for month
}

total.cases<-rbind(monthly.cases[[1]],monthly.cases[[2]])
for (i in 3:length(monthly.cases)){
  total.cases<-rbind(total.cases, monthly.cases[[i]])
}

#Remove zeroes
rem<-c()
for (i in 1:nrow(total.cases)){
  if (total.cases$total[i]==0) rem<-c(rem,i)
}

if (length(rem)>0){
  total.cases<-total.cases[-rem,]
}


## Separate Canadian metadata
meta.all<-meta #backup for downstream

#subset the metadata
meta.can<-meta[which(meta$country=="Canada"),]
l.can<-nrow(meta.can) 

#remove canada from meta, but see backup
meta<-meta[meta$country!="Canada",]
length(meta$new.names)==length(unique(meta$new.names))

#change any NA month issues
for (i in 1:nrow(meta.can)){
  if(length(unlist(strsplit(meta.can$date[i],split="-")))==2){
    meta.can$month[i]<-format(as.Date(paste(meta.can$date[i],"01",sep="-")),"%Y-%m")
  }
}

#some wrongly put canada as prov (montreal submitting lab)
if (length(which(meta.can$division=="Canada"))>1){
  meta.can$division[which(meta.can$division=="Canada")]<-"Quebec"
}

can.seq.rep<-table(meta.can$month,meta.can$division)

# if(!file.exists("CanSequenceRepresentationPreSample.csv")){
#   write.csv(can.seq.rep,"CanSequenceRepresentationPreSample.csv")
# }

## Calculate probability of selecting Canadian sequences based on provinces' proportional contribution to new cases in Canada
## In this BUILD keep ALL CANADIAN SEQUENCES
# table(meta.can$division, meta.can$Lineage)
# table(meta.can$Location, meta.can$month)

#Read in the provincial summary data
prov<-read.csv("Canada_cases/20210211_covid19_casesbyprov.csv", header=T)
prov<-prov[-1,]
prov<-prov[-which(prov$prname=="Canada" |prov$prname=="Repatriated travellers"),] 
#unfortunately we don't know where people repatriated to/from...
prov$date<-as.Date(prov$date, "%Y-%m-%d") 

#add a column for month
prov$month<-format(as.Date(prov$date), "%Y-%m")

#group by province and month
monthly.cases.prov<-as.data.frame(prov %>% dplyr::group_by(prname, month) %>% 
                                    summarize(.groups="rowwise", n.cases=sum(numtoday)))

#calculate the proportion of cases in a given month from each province 
monthly.cases.prov$proportion.cases<-NA                                                                      
for (i in 1:nrow(monthly.cases.prov)){
  monthly.cases.prov$proportion.cases[i]<-monthly.cases.prov$n.cases[i]/sum(monthly.cases.prov$n.cases [monthly.cases.prov$month==monthly.cases.prov$month[i]])
}

#Count total and proportion of sequences
#Need to calculate number of seqs in each prov in each month
prov.seqs<-as.data.frame(meta.can %>% #Use canadian meta
  dplyr::group_by(month,division) %>% 
  dplyr::summarize(.groups = "rowwise", n.sequence = n()))

#add a column for proportion of seqs that were in that prov compared to country
prov.seqs$proportion.sequence<-NA
for (i in 1:nrow(prov.seqs)){
  prov.seqs$proportion.sequence[i]<-prov.seqs$n[i]/sum(prov.seqs$n [prov.seqs$month==prov.seqs$month[i]])
}

#Remove zeroes, don't need to run
# rem<-c()
# for (i in 1:nrow(total.cases.prov)){
#   if (monthly.cases.prov$n.cases[i]==0) rem<-c(rem,i)
# }
# monthly.cases.prov<-monthly.cases.prov[-rem,]

#check congruency
colnames(monthly.cases.prov)[1]<-"division"
# monthly.cases.prov$division[which(!monthly.cases.prov$division %in% prov.seqs$division )] #these are provs with cases and no sequences in a given month

#ensure classes match,full match if zero
which(!prov.seqs$division %in% monthly.cases.prov$division) 
which (!prov.seqs$month %in% monthly.cases.prov$month)
prov.seqs$month[which (!prov.seqs$month %in% monthly.cases.prov$month)]
prov.seqs[which (!prov.seqs$month %in% monthly.cases.prov$month),]
# table(meta.can$date) 
#because I didn't remove incomplete dates, there are instances of '2020'... as well as just y-m, ignore for now

#do the join
monthly.cases.prov<-left_join(monthly.cases.prov, prov.seqs, by=c("division","month") )

#turn NAs into zeroes, dont need to run
# monthly.cases.prov[is.na(monthly.cases.prov)] <- 0
#actually, NAs probably better to exclude zeroes from plotting

#Calculate a probability of selecting an individual sequence in downsampled dataset
monthly.cases.prov<-monthly.cases.prov %>% mutate ("probability"=proportion.cases/n.sequence)

#Change any NA probabilities to 0
monthly.cases.prov$probability[which(is.na(monthly.cases.prov$probability))]<-0

#calculate the net.difference and fold.difference between cases and sequences by month and prov
monthly.cases.prov <- monthly.cases.prov %>% mutate(net.difference=n.cases-n.sequence, fold.difference=n.sequence/n.cases)

#how about the proportion seqs/proportion cases (ie how over-represeneted (>1) or under-represented (<1)
monthly.cases.prov<- monthly.cases.prov %>% mutate(overrep=proportion.sequence/proportion.cases)

#we join these probabilities onto meta.can 

# meta.can$month<-format(as.Date(meta.can$date), "%Y-%m")
meta.can<-left_join(meta.can,monthly.cases.prov, by=c("division","month"))
meta.can$probability[which(is.na(meta.can$probability))]<-0

#change this for plotting
meta.can$division<-str_replace_all(meta.can$division, "Newfoundland and Labrador", "Newfoundland\nand Labrador")

## Calculate number and proportion of sequences contributed to whole by country and by month. Join to cases by country, month. Join to metaseq
#Need to calculate number of seqs in each country in each month
country.seqs<-meta %>% 
  dplyr::group_by(month,country) %>% 
  dplyr::summarize(
    .groups="rowwise",
    n = n())
country.seqs<-as.data.frame(country.seqs)

#add a column for proportion of cases that were in that country compared to globe
country.seqs$prop<-NA
for (i in 1:nrow(country.seqs)){
  country.seqs$prop[i]<-country.seqs$n[i]/sum(country.seqs$n [country.seqs$month==country.seqs$month[i]])
}

#change some column names to facilitate merge
colnames(country.seqs)[3:4]<-c("n.sequence","proportion.sequence")
colnames(total.cases)<-c("country","n.cases","proportion.cases","month")

#do the join
total.cases<-left_join(total.cases,country.seqs,by=c("country","month") )

#Calculate a probability of selecting an individual sequence in downsampled dataset
total.cases<-total.cases %>% mutate ("probability"=proportion.cases/n.sequence)

#Change any NA probabilities to 0
if (length(is.na(total.cases$probability))>0){
  total.cases$probability[which(is.na(total.cases$probability))]<-0
}

#calculate the net.difference and fold.difference between cases and sequences by month and country
# proportion seqs/proportion cases (ie how over-represeneted (>1) or under-represented (<1)
total.cases <- total.cases %>% mutate(net.difference=n.cases-n.sequence, 
                                      fold.difference=n.sequence/n.cases,
                                      overrep=proportion.sequence/proportion.cases)

#join these probabilities onto meta 
meta<-left_join(meta, total.cases, by=c("country","month"),keep=F)
if (length(which(is.na(meta$probability)))>0){
  meta$probability[which(is.na(meta$probability))]<-0
}
#not sure why, but there were dups genererated that have to be removed
dupys<-which(duplicated(meta$new.names))
if (length(dupys)>0){
  meta<-meta[-dupys,]
}
length(meta$new.names)==length(unique(meta$new.names)) #should be true

# Subsample ALL Canadian sequences + remainder up to specified total nsamp of global sequences 
## Set up subsampling conditions

#What proportion of total global cases are canadian?
# sum(prov$numtoday)/sum(total.cases$n.cases) #0.007 # knowingly oversampling canada for focus 

#as this is a canada focused analyis, take ALL Canadian, remainder global
n.can<-nrow(meta.can) ### CHANGE here if you don't want to take all Canadian
n.glob<-n.samp-n.can
if (n.glob<1){print("Warning, specify more samples in n.samp or fewer n.can")}

#Setup a list for downsampled meta.seq x b bootstraps (specified above)
#will restrict alignments to IDs in metaseq after subsampling
samp.meta<-replicate(n=b, vector())
names(samp.meta)<-as.character(1:b)

####subsample Canadian sequences#####

##USE THIS IF KEEP ALL SEQs
for (i in 1:b){
  samp.meta[[i]]<-meta.can}

### UNCOMMENT/RUN below if subsampling Canada 
#Setup up number of seqs per month table
#Take an equal number of seqs from each month if possible
#If there aren't enough in a given month, take them all and then re-distribute remaining seqs between other months

# #table of available seqs by month
# t<-table(meta.can$month)
# 
# #distribute the desired sequences between the number of months with sequences
# mm<-trunc(n.can/length(t))
# #check which are "short months" where all seqs will be taken
# short<-names(which(t<mm))
# sh<-0
# while(sh==0){ #until sh=1, optimize number of seqs for months to be evenly distributed where possible
#     t.sh<-t[!names(t)%in%short]
#     #calculate month-specific mm.2, the number of seqs to take each month for non-short months
#     mm<-trunc((n.can- sum(t[short])) / length(t.sh) )
#     #which are still short?
#     st.short<-which(t.sh<mm)
#     #if all remaining months have that many seqs, great, sample that many
#     ifelse(length(st.short)==0, sh<-1, #if none are short, break the loop
#            short<-c(short,names(st.short)) #if at least one still short, add it to short list and repeat
#            )
# }
# t[!names(t)%in%short]<-mm
# 
# #if still don't have enough seqs (because of truncating above), add the diff
# if (sum(t)<n.can) {
#   diff<-as.numeric(n.can)-sum(t)
#   t[!names(t)%in%short][1:diff]<-t[!names(t)%in%short][1:diff] + 1
# }

#subsample according to provinces' monthly case proportions to the Canadian total
# for (i in 1:b){
  
  # #make an empty object for final IDs
  # can.samp<-vector()
  # #make sure no probabilities are zero
  # meta.can$probability[which(meta.can$probability==0)]<-0.00000000000001
  # #for each month in the table of seqs for a lineage, extract # seqs identified above
  # for (j in 1:length(t)){
  #   row.mm<-which(meta.can$month==names(t[j]))
  #   #take sample of size specified in table t, but if short, will just take that many
  #   can.samp<-c(can.samp, sample(meta.can$new.names[row.mm],
  #                     size=t[j], replace=FALSE,
  #                     prob=meta.can$probability[row.mm]) )
  # }
  # #check
  # # length(can.samp)==n.can
  # #add these to the meta.can (already populated with canada)
  # samp.meta[[i]]<-meta.can[which(meta.can$new.names %in% can.samp),]
# }

#check number of seqs in each boot
for (i in 1:b){
  print(nrow(samp.meta[[i]]))
  length(which(names(align) %in% samp.meta[[i]]$new.names))
}

####subsample global sequences#####
# Setup up number of seqs per month table
#Take an equal number of seqs from each month if possible
#If there aren't enough in a given month, take them all and then re-distribute remaining seqs between other months
#table of available seqs by month
t<-table(meta$month)
#distribute the desired sequences between the number of months with sequences
mm<-trunc(n.glob/length(t))
#check which are "short months" where all seqs will be taken
short<-names(which(t<mm))
sh<-0
while(sh==0){ #until sh=1, optimize number of seqs for months to be evenly distributed where possible
    t.sh<-t[!names(t)%in%short]
    #calculate month-specific mm.2, the number of seqs to take each month for non-short months
    mm<-trunc((n.glob- sum(t[short])) / length(t.sh) )
    #which are still short?
    st.short<-which(t.sh<mm)
    #if all remaining months have that many seqs, great, sample that many
    ifelse(length(st.short)==0, sh<-1, #if none are short, break the loop
           short<-c(short,names(st.short)) #if at least one still short, add it to short list and repeat
           )
}
t[!names(t)%in%short]<-mm

#if still don't have enough seqs (because of truncating above), add the diff
if (sum(t)<n.glob) {
  diff<-as.numeric(n.glob)-sum(t)
  t[!names(t)%in%short][1:diff]<-t[!names(t)%in%short][1:diff] + 1
}

sum(t)==n.glob #should be true

#make sure no probabilities in meta are zero
if (any(meta$probability==0)){
  meta$probability[which(meta$probability==0)]<-0.00000000000001
}
#and none are negative
if (any(meta$probability<0)){
  meta$probability[which(meta$probability<0)]<-0.00000000000001
}

#subsample according to monthly case proportions to the global total
set.seed(111)
for (i in 1:b) {
  #make an empty object for final IDs
  glob.samp<-vector()
  #for each month in the table of seqs for a lineage, extract # seqs identified above
  for (j in 1:length(t)){
    row.mm<-which(meta$month==names(t[j]))
    #take sample of size specified in table t, but if short, will just take that many
    glob.samp<-c(glob.samp, sample(meta$new.names[row.mm],
                      size=t[j], replace=FALSE,
                      prob=meta$probability[row.mm]) )
  }
  #check
  # length(glob.samp)==n.glob
  #add these to the meta (already populated with canada)
  samp.meta[[i]]<-rbind(samp.meta[[i]],meta[which(meta$new.names %in% glob.samp),])
}

#check number of seqs in each lineage
for (i in 1:b){
  print(nrow(samp.meta[[i]])==n.samp)
}

#check for wuhan-hu-1 (needed to root)
#should always be selected because we aren't subsampling from months with low number of seqs (ie December)
for (i in 1:b){
  print(any(samp.meta[[i]]$date=="2019-12-26")) #this is the only sample on this day
}
  

#look at representation
# table(samp.meta[[i]]$country)

## Generate bootstrap alignments from bootstrap subsampled metadata
#make alignments using sequence names in the bootstrap dataframes from both the global and canada samples
samp.align<-replicate(vector(),n=b)
length(which(names(align) %in% samp.meta[[i]]$new.names))
nrow(samp.meta[[1]])

for (i in 1:b){
  samp.align[[i]]<-align [names(align) %in% samp.meta[[i]]$new.names]
}

# check lengths are as expected
for (i in 1:b){
  print(length(samp.align[[i]]))
}


##Export full subsampled metadata and alignments for bootstraps
#export metadata
for (i in 1:b){
  write.csv(samp.meta[[i]],file = paste(meta.out,n.samp ,"_",i,".csv",sep=""))
}

#export fasta 
for (i in 1:b){
  writeXStringSet(samp.align[[i]], paste(align.out,n.samp ,"_",i, ".fasta",sep=""),format="fasta")
}

