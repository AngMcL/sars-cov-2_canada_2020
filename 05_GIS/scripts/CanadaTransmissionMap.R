#!/usr/bin/env Rscript
# Canadian covid transmission heatmap with US context

#usage:
# $Rscript CanadaTransmissionMap.R <>

#setup libraries
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(ggplot2)
library(rgdal)
library(rgeos)
library(maptools)
library(ggsn)
library(broom)
library(RColorBrewer)

#### Set up transmission input ####

#### Inputs and outputs here ####
args<-commandArgs(trailingOnly=TRUE)
transmission.in<-args[1]
# manual setup
transmission.in<-"DF/BOOTS_Provinces_monthly_transmission.csv"

#### Read in shapefiles ####
can<-readOGR("shapefiles/Can_provinces/lpr_000b16a_e.shp")
can.fort<-fortify(can)
can$id <- row.names(can)
can.fort1<-dplyr::left_join(can.fort, can@data, by="id")

us<-readOGR("shapefiles/USA_states/US_states_nad83.shp")
# us@proj4string
# us@bbox
us.fort<-fortify(us)
us$id <- row.names(us)
us.fort1<-dplyr::left_join(us.fort, us@data, by="id")

####  Canada map setup with partial US ####
#plot lims
can.lim.x<-range(can.fort1$long)
can.lim.y<<-range(can.fort1$lat)

#make lims narrower all around
#adj factors
x.adj<-(can.lim.x[2]-can.lim.x[1])/100
y.adj<-(can.lim.y[2]-can.lim.y[1])/100
#apply adj
can.lim.x[1]<-can.lim.x[1]+3*x.adj  #West
can.lim.x[2]<-can.lim.x[2]-3*x.adj #East
can.lim.y[1]<-can.lim.y[1]+3*y.adj  #South
can.lim.y[2]<-can.lim.y[2]-3*y.adj #North

#for ggsave
xx<-diff(can.lim.x)/100000/2
yy<-diff(can.lim.y)/100000/2

# Basemap
# for north and scalebar
# bb<-data.frame(long = can.lim.x*0.959, lat = can.lim.y*0.959) #want the north arrow slightly up and right
# bb2<-data.frame(long = can.lim.x, lat = can.lim.y)
# note these are quite slow
# p<-ggplot(data=can.fort1,aes(x=long,y=lat,group=group))+
#   geom_polygon(data=can.fort1,fill="burlywood",colour="black",lwd=0.5)+
#   geom_polygon(data=us.fort1,fill="grey65",colour="black",lwd=0.5)+
#   labs(fill=NULL)+
#   theme(axis.text = element_blank(),
#   axis.title = element_blank(),
#   axis.ticks = element_blank(),
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(),
#   panel.background = element_rect(fill="lightcyan1"))+
#   coord_fixed()+
#   coord_cartesian(xlim =c(can.lim.x),ylim=c(can.lim.y))
# p2<-p+north(bb,symbol=10,scale=0.06, location="topright")
# p3<-p2+scalebar(data=bb2,dist = 500, dist_unit = "km",
#                 st.size=6, st.dist=0.01,
#                 transform = FALSE, model = "NAD83",
#                 location="topright", inherit.aes =FALSE)
# ggsave(p3,filename = "maps_out/CanadaBaseMap.png",width=xx,height=yy,units = "in")


####  Prepare for plotting heatmap ####
#remove unncessary cols
can@data<-can@data[,-which(colnames(can@data) %in% c("PRFNAME","PREABBR","PRFABBR","PRNAME"))]

#read in transmission data from ancestral reconstruction analysis
prov.trans<-read.csv(transmission.in)
prov.trans<-prov.trans[,-1]

#Remove columns that aren't going to be plotted
prov.trans<-prov.trans %>% dplyr::filter (!Type %in%  c("N.tips","N.Sample","N.provincial","Prop.provincial","N.domestic","Prop.domestic","N.international","Prop.international", "Prop.USA","N.USA"))

#repalce the name of all total international to avoid issues with separate 
prov.trans$Type<-str_replace_all(prov.trans$Type, "N.total.international","N.totalinternational")
prov.trans$Type<-str_replace_all(prov.trans$Type, "Prop.total.international","Prop.totalinternational")

#remove sd
prov.trans<-prov.trans[,-ncol(prov.trans)]

#get rid of type (only have prop international in here)
prov.trans<-prov.trans[,-3]

#make it wide by separating each month 
prov.trans.df<-prov.trans %>% pivot_wider(names_from = 2, values_from = 3)
colnames(prov.trans.df)[-1]<-paste("Prop.totalinternational.",colnames(prov.trans.df)[-1],sep="")

#use english province name
can$PRNAME<-as.character(can$PRENAME)
can@data<-can@data[,-which(colnames(can@data)=="PRENAME")]

#change the colname for prov.trans
colnames(prov.trans.df)[which(colnames(prov.trans.df)=="Province")]<-"PRNAME"

#workaround for Maritimes: triplicate row and name by prov
prov.trans.df<-rbind(prov.trans.df, prov.trans.df[which(prov.trans.df$PRNAME=="Maritimes"),],
                     prov.trans.df[which(prov.trans.df$PRNAME=="Maritimes"),],
                     prov.trans.df[which(prov.trans.df$PRNAME=="Maritimes"),])
prov.trans.df$PRNAME[7:9]<-c("Nova Scotia","New Brunswick","Newfoundland and Labrador")
prov.trans.df<-prov.trans.df[-which(prov.trans.df$PRNAME=="Maritimes"),]

#check full match
all(prov.trans.df$PRNAME %in% can$PRNAME)
# can$PRNAME[which(!can$PRNAME%in% prov.trans.df$PRNAME)] #provinces we have no seqs from
#some provinces with no data - will just have NAs

#join the provincial transmission dataframe by PRNAME
can@data<-can@data %>% left_join(prov.trans.df, by="PRNAME")
# head(can@data)

#change the spdf into a df using tidy
# can.fort<-fortify(can)
can.fort<-broom::tidy(can) #replaces fortify

#lost region names, but the ordering of "id" is the same as that stored in can@data$region
# Recover row name 
can_df<-data.frame(can@data$PRNAME)
names(can_df)<-c("PRNAME")

# Create and append "id"
can_df$id<-seq(0,nrow(can_df)-1)
#Merge row names with new dataframe using "id": Finally, let us put the names back into the new dataframe:
can.fort.df <- plyr::join(can.fort, can_df, by="id")
# head(can.fort.df, n = 2) # peak at the fortified data
# head(can@data, n = 2) # final check before join 
can.fort.df2 <- left_join(can.fort.df, can@data, by="PRNAME")
can.fort.df2<-can.fort.df2[,-which(colnames(can.fort.df2)=="id.y")]

#change from wide to long format for hiv diagnoses using tidyr::gather
#data %>% gather(key, value, ..., na.rm = FALSE, convert = FALSE)
#SLOWWWWW
#get rid of big objects we're done with.
rm(can.fort.df)
rm(can.fort)

#could replace with pivot_longer or melt
#watch out for the as.numeric(10:(ncol(can.fort.df2))) - might need to change based on cols included
can.fort.df3 <-can.fort.df2 %>%
  gather(key="Measurement",value="Value",as.numeric(10:(ncol(can.fort.df2))),na.rm=TRUE)
  #bc.fort1<-bc.fort%>%gather(key="Measurement",value="Value",as.numeric(9:20),na.rm=TRUE) 
rm(can.fort.df2)

#separate measurement into interval and measurement
can.fort.df4<-can.fort.df3 %>% 
  separate(col=`Measurement`, into=c("Calc","Param", "Interval"),sep = "\\.")
rm(can.fort.df3)

#### Setup df with prov names and centroids for labels ####
idList <- unique(can.fort1$id)
prList<-unique(can.fort1$PRENAME)
centroids.df <- as.data.frame(coordinates(can))
names(centroids.df) <- c("long", "lat")
cent.df <- data.frame(id = idList,PRNAME=prList, centroids.df)

#join the Incidences into cent.df
# cent.df<-left_join(cent.df,canPopInc,by="PRNAME")
#make a lower lat for the Inc
cent.df$latlow<-cent.df$lat
cent.df$longlow<-cent.df$long
cent.df$verbose<-paste(cent.df$Inc," cases\nper 100,000",sep="")
#custom changes for names
cent.df$PRNAME[which(cent.df$PRNAME=="British Columbia")]<-"British\nColumbia"
cent.df$PRNAME[which(cent.df$PRNAME=="Northwest Territories")]<-"Northwest\nTerritories"
# only label as Maritimes
cent.df$PRNAME[which(cent.df$PRNAME=="New Brunswick")]<-""
cent.df$PRNAME[which(cent.df$PRNAME=="Nova Scotia")]<-""
cent.df$PRNAME[which(cent.df$PRNAME=="Newfoundland and Labrador")]<-"Maritimes"
cent.df$PRNAME[which(cent.df$PRNAME=="Prince Edward Island")]<-"PEI"

#custom changes for text locations
# cent.df$longlow[which(cent.df$PRNAME=="Newfoundland\nand Labrador")]<-cent.df$longlow[which(cent.df$PRNAME=="Newfoundland\nand Labrador")]+300000
# cent.df$latlow[which(cent.df$PRNAME=="Newfoundland\nand Labrador")]<-cent.df$latlow[which(cent.df$PRNAME=="Newfoundland\nand Labrador")]+100000
# cent.df$longlow[which(cent.df$PRNAME=="Nova\nScotia")] <-cent.df$longlow[which(cent.df$PRNAME=="Nova\nScotia")]+300000
# cent.df$longlow[which(cent.df$PRNAME=="New\nBrunswick")]<- cent.df$longlow[which(cent.df$PRNAME=="New\nBrunswick")]-10000
# cent.df$latlow[which(cent.df$PRNAME=="New\nBrunswick")]<- cent.df$latlow[which(cent.df$PRNAME=="New\nBrunswick")]-10000
cent.df$latlow[which(cent.df$PRNAME=="PEI")]<-cent.df$latlow[which(cent.df$PRNAME=="PEI")]+100000
cent.df$longlow[which(cent.df$PRNAME=="Maritimes")]<-cent.df$longlow[which(cent.df$PRNAME=="Maritimes")]-20000
cent.df$latlow[which(cent.df$PRNAME=="Maritimes")]<-cent.df$latlow[which(cent.df$PRNAME=="Maritimes")]-80000
cent.df$longlow[which(cent.df$PRNAME=="Yukon")] <-cent.df$longlow[which(cent.df$PRNAME=="Yukon")]-50000
cent.df$latlow[which(cent.df$PRNAME=="Saskatchewan")]<- cent.df$latlow[which(cent.df$PRNAME=="Saskatchewan")]-80000
cent.df$latlow[which(cent.df$PRNAME=="Manitoba")]<- cent.df$latlow[which(cent.df$PRNAME=="Manitoba")]+50000
cent.df$longlow[which(cent.df$PRNAME=="Manitoba")]<- cent.df$longlow[which(cent.df$PRNAME=="Manitoba")]+30000
cent.df$longlow[which(cent.df$PRNAME=="British\nColumbia")]<- cent.df$longlow[which(cent.df$PRNAME=="British\nColumbia")]-60000
cent.df$longlow[which(cent.df$PRNAME=="Northwest\nTerritories")]<-cent.df$longlow[which(cent.df$PRNAME=="Northwest\nTerritories")]-100000
cent.df$latlow[which(cent.df$PRNAME=="Quebec")]<- cent.df$latlow[which(cent.df$PRNAME=="Quebec")]-100000

#### FIGURE: hotspot/chloropleth maps using transmission data ####
#for north arrow
bb<-data.frame(long = c(can.lim.x[1],can.lim.x[2]*0.966), lat = c(can.lim.y[1]*1.08,can.lim.y[2])) #north arrow
bb2<-data.frame(long = can.lim.x, lat = c(can.lim.y[1]*0.94,can.lim.y[2])) #scalebar

# PLOT IT UP!
label<-'Proportion\ninternational\ntransmission\nApril 2020' #for legend
p<-can.fort.df4  %>%
  filter(Interval=="2020-04") %>%
  filter(Calc=="Prop") %>%
  filter(Param=="totalinternational") %>%
    ggplot(aes(x=long, y=lat, group=group))+
      geom_polygon(data=can.fort1,fill="grey65",colour="grey35",lwd=0.2)+
      geom_polygon(data=us.fort1,fill="grey85",colour="grey35",lwd=0.2)+
      geom_polygon(aes(fill=Value,group = group), lwd=0.2,  colour="grey35",stat="identity") +
      labs(x = NULL, y=NULL,fill = label) +
      scale_fill_gradientn(colours=brewer.pal(n=8,"YlOrRd"),na.value="grey65")+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = 'lightcyan1'),
            plot.margin = unit(c(0, 0, 0, 0), "null"),
            panel.margin = unit(c(0, 0, 0, 0), "null"),            
            axis.line = element_blank(),
            axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            legend.title = element_text(face="bold",size=10),
            legend.text = element_text(face="bold",size=10,hjust=0.5),
            legend.position = c(0.97, 0.97),legend.justification=c(1,1),
            legend.spacing=unit(c(0.05,0.05,0.05,0.05),"cm")) +
      guides(fill=guide_colorbar(barheight = unit(2,"cm")))+
      coord_fixed()+
      coord_cartesian(xlim =c(can.lim.x),ylim=c(can.lim.y))+
      #province labels
      geom_text(data=cent.df,aes(x=longlow, y=latlow,label=PRNAME,group=NULL),size=2.8,fontface="bold",lineheight = .7)+ 
      #alaska
      annotate(geom="text",label="USA",y=(bbox(can)[2,2]-800000), x=(bbox(can)[1,1]+460000),size=4)+
      #mainland US
      annotate(geom="text",label="USA",y=(bbox(can)[2,1])+100000, x=mean(bbox(can)[1,])-805000,size=4)+
      #Canadia
      annotate(geom="text",label="CANADA",y=mean(bbox(can)[2,])-300000,
               x=mean(bbox(can)[1,])-800000,size=4)

p2<-p+north(bb,symbol=1,scale=0.08, location="bottomright")
p3<-p2+scalebar(data=bb2,dist = 500, dist_unit = "km",
                st.size=2.5, st.dist=0.02,
                transform = FALSE, model = "NAD83",
                location="bottomright", inherit.aes =FALSE, border.size=0.3)
ggsave(p3,file=paste("maps_out/","Canada_ProportionIntlTransm_April2020",".png",sep=""),width=xx/5,height=yy/5,units="in")


#### re-create the bar plots and alluvial plots to make a composite image for the manuscript ####
## Redundant for regular folk: don't run 
# #read in the dataframes from the other R project
# SUM.prov.summ.mo.long.all <-read.csv(file="../03c_MLstate_boots/DF/SUM.prov.summ.mo.long.all.csv")
# sum.bw.Prov.summary.long <-read.csv(file="../03c_MLstate_boots/DF/sum.bw.Prov.summary.long.csv")
# 
# #setup nec libraries
# library(ggalluvial)
# library(cowplot)
# 
# #setup themes
# pubThemeDate<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                     panel.background=element_rect("grey95"), 
#                     axis.line = element_line(colour = "black"),
#                     text=element_text(size=10,face="bold"),
#                     legend.key.size = unit(0.4,"line"),
#                     legend.text=element_text(size=6),
#                     axis.text.x=element_text(angle = 45,hjust=1,size=rel(1)))
# pubTheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#             panel.background=element_rect("grey95"), axis.line = element_line(colour = "black"),
#             legend.key.size = unit(0.5,"line"),
#             text=element_text(size=10,face="bold"),
#             legend.text=element_text(size=6))
# prop.colz<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c")
# 
# ## read in color schemes geographies
# #import a tsv of name and hex color for LOCATIONS
# globalPalette<-read.table("../03c_MLstate_boots/DF/globalcolorsNew.tsv",sep="\t")
# ## make color scheme
# glob.colz<-row.names(globalPalette)
# globalPalette.ch<-as.character(globalPalette$globalPalette)
# names(globalPalette.ch)<-glob.colz
# 
# GlobColScale<-scale_colour_manual(name = "Location",values = globalPalette.ch,na.value="grey60")
# GlobFillScale<-scale_fill_manual(name = "Location",values = globalPalette.ch,na.value="grey60")
# 
# #order the provinces according to greatest source
# prov.ord<-c("Quebec","Ontario","British Columbia","Alberta","Manitoba","Maritimes")
# 
# SUM.prov.summ.mo.long.all$Province<-factor(SUM.prov.summ.mo.long.all$Province,levels=prov.ord)
# 
# #make the plots
# PZ4<-SUM.prov.summ.mo.long.all  %>%
#   # filter(month %in% c("2020-03","2020-04","2020-05")) %>%
#   ggplot()+
#   geom_bar(aes(x=month, y=mean.n, fill=Type, group=Type), stat="identity", position=position_stack(),alpha=0.8)+
#   pubThemeDate+
#   theme(strip.text = element_text(size=9),
#         legend.position="top",
#         plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
#         text =element_text(size=12),
#         legend.text = element_text(size=9))+
#   guides(fill=guide_legend(keywidth = 1,keyheight = 1))+
#   labs(x=NULL, y="Sampled transmssion events")+
#   scale_fill_manual(values=prop.colz,name = "Transmission source", labels = c("USA","Other International","Between-province","Within-province"))+
#   scale_x_discrete(breaks=c("2020-03","2020-04","2020-05"),
#                    labels=c("Mar 2020","Apr 2020","May 2020"),expand=c(.01,0))+
#   scale_y_continuous(expand=c(0.01,0))+
#   facet_grid(rows = vars(Param), cols=vars(Province),scales="free")
# # PZ4
# 
# #order the provinces
# sum.bw.Prov.summary.long$geo.type<-factor(sum.bw.Prov.summary.long$geo.type,levels=c("Source","Recipient"))
# sum.bw.Prov.summary.long$geo <- factor(sum.bw.Prov.summary.long$geo, levels=prov.ord)
# 
# PZ3<- ggplot(sum.bw.Prov.summary.long,
#        aes(x = geo.type, stratum = geo, alluvium = subject,
#            y = mean.n,
#            fill = geo, label = geo)) +
#   scale_x_discrete(expand = c(0.01,0.01)) +
#   scale_y_continuous(expand = c(0,0)) +
#   geom_flow(alpha = .6,width=0.6) +
#   geom_stratum(alpha = .8,width=0.6) +
#   geom_text(stat = "stratum", size = 3.2,min.y=4,fontface="bold") +
#   pubTheme+
#   theme(legend.position = "none", 
#         axis.line = element_blank(), 
#         text =element_text(size=12),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_text(hjust=c(0.5,0.6)))+
#   labs(x=NULL,y="Between-province transmission events")+
#   GlobFillScale
# PZ3
# 
# #Make fake plots to work the gridding out
# # ppdata<-data.frame(x=1:10,y=1:10)
# # pp1<-ggplot(ppdata)+geom_point(aes(x=x,y=y))
# # bottom_row <- plot_grid(pp1,pp1, labels = c('B', 'C'), label_size = 14,ncol=2, rel_widths = c(0.55,0.5))
# # finalPZ<-plot_grid(pp1, bottom_row, labels = c('A', ''), label_size = 14, nrow=2, greedy=T,rel_heights = c(0.4,0.6))
# # finalPZ
# 
# #make the composite plot
# bottom_row <- plot_grid(p3,PZ3, labels = c('B', 'C'), label_size = 14,ncol=2,rel_widths = c(0.55,0.5))
# 
# finalPZ<-plot_grid(PZ4, bottom_row, labels = c('A', ''), label_size = 14, nrow=2,rel_heights = c(0.5,0.5))
# 
# ggsave(finalPZ, file="OutputMaps/Figure5_BW-province_propIntl.png",width=8.34,height=8.6,units = "in")


#### Map of Canada with labeled cumulative incidence ####
## Read and clean cumulative incidence
## read in and calculate cum.inc.
canPop<-read.csv("../01_subSample/Canada_cases/2016CanCensusPopProv.csv")
canInc<-read.csv("../01_subSample/Canada_cases/20210211_covid19_casesbyprov.csv")
# head(canPop)
canPop<-canPop[,c(2,4)]
colnames(canPop)<-c("province","population")
# head(canInc)
canInc<-canInc[,c(2,4,9)]
colnames(canInc)<-c("province","date","totalcases")
canInc2<-canInc[canInc$date=="2021-02-11",]
canInc2<-canInc2%>%filter(province!="Canada") %>%filter(province!="Repatriated travellers")
# head(canInc2)

#merge em
canPopInc<-left_join(canInc2,canPop,by="province")

#calculate cumulative incidence = total number of confirmed covid cases/population*100,000
canPopInc<-canPopInc %>% mutate(Inc=totalcases/population*100000)
canPopInc$Inc<-round(canPopInc$Inc) #round for clean plotty
# head(canPopInc)

#use english province name
can<-readOGR("shapefiles/Can_provinces/lpr_000b16a_e.shp")
can$PRNAME<-as.character(can$PRENAME)
can@data<-can@data[,-which(colnames(can@data)=="PRENAME")]

#change the colname
colnames(canPopInc)[which(colnames(canPopInc)=="province")]<-"PRNAME"

#check full match
all(canPopInc$PRNAME %in% can$PRNAME)

#join the provincial transmission dataframe by PRNAME
can@data<-can@data %>% left_join(canPopInc, by="PRNAME")

#change the spdf into a df using tidy
can.fort<-broom::tidy(can) #replaces fortify

#lost region names, but the ordering of "id" is the same as that stored in can@data$region
# Recover row name 
can_df<-data.frame(can@data$PRNAME)
names(can_df)<-c("PRNAME")
# Create and append "id"
can_df$id<-seq(0,nrow(can_df)-1)
#Merge row names with new dataframe using "id": Finally, let us put the names back into the new dataframe:
can.fort.df <- plyr::join(can.fort, can_df, by="id")
# head(can.fort.df, n = 2) # peak at the fortified data
# head(can@data, n = 2) # final check before join 
can.fort.df2 <- left_join(can.fort.df, can@data, by="PRNAME")
can.fort.df2<-can.fort.df2[,-which(colnames(can.fort.df2)=="id.y")]



#### Make incidence map ####
#read in global palette from reconstruction output
globalPalette<-read.table("../03_stateReconstruction/DF/globalcolors.tsv",sep="\t")

## make color scheme
glob.colz<-row.names(globalPalette)
globalPalette.ch<-as.character(globalPalette$globalPalette)
names(globalPalette.ch)<-glob.colz

#add the separate maritimes provs
mar.colz<-rep("#a86030",times=3)
names(mar.colz)<-c("Nova Scotia","Newfoundland and Labrador","New Brunswick")
globalPalette.ch<-c(globalPalette.ch,mar.colz)

GlobColScale<-scale_colour_manual(name = "Location",values = globalPalette.ch,na.value="grey60")
GlobFillScale<-scale_fill_manual(name = "Location",values = globalPalette.ch,na.value="grey60")

idList <- unique(can.fort1$id)
prList<-unique(can.fort1$PRENAME)
centroids.df <- as.data.frame(coordinates(can))
names(centroids.df) <- c("long", "lat")
cent.df <- data.frame(id = idList,PRNAME=prList, centroids.df)

#join the incidences into cent.df
cent.df<-left_join(cent.df,canPopInc,by="PRNAME")
#Adjust some positions of prov labels
cent.df$lat[which(cent.df$PRNAME=="Prince Edward Island")]<-cent.df$lat[which(cent.df$PRNAME=="Prince Edward Island")]+70000
cent.df$long[which(cent.df$PRNAME=="Newfoundland and Labrador")]<-cent.df$long[which(cent.df$PRNAME=="Newfoundland and Labrador")]+110000
cent.df$long[which(cent.df$PRNAME=="Nova Scotia")]<-cent.df$long[which(cent.df$PRNAME=="Nova Scotia")]+60000
cent.df$lat[which(cent.df$PRNAME=="Nova Scotia")]<-cent.df$lat[which(cent.df$PRNAME=="Nova Scotia")]-20000
cent.df$long[which(cent.df$PRNAME=="New Brunswick")]<-cent.df$long[which(cent.df$PRNAME=="New Brunswick")]-40000
cent.df$lat[which(cent.df$PRNAME=="Alberta")]<-cent.df$lat[which(cent.df$PRNAME=="Alberta")]+70000
cent.df$lat[which(cent.df$PRNAME=="Saskatchewan")]<-cent.df$lat[which(cent.df$PRNAME=="Saskatchewan")]-100000
cent.df$long[which(cent.df$PRNAME=="Saskatchewan")]<-cent.df$long[which(cent.df$PRNAME=="Saskatchewan")]-50000
cent.df$lat[which(cent.df$PRNAME=="Manitoba")]<-cent.df$lat[which(cent.df$PRNAME=="Manitoba")]+50000

#custom changes for names
cent.df$PRNAME[which(cent.df$PRNAME=="Newfoundland and Labrador")]<-"Newfoundland\nand Labrador"
cent.df$PRNAME[which(cent.df$PRNAME=="Northwest Territories")]<-"Northwest\nTerritories"
cent.df$PRNAME[which(cent.df$PRNAME=="British Columbia")]<-"British\nColumbia"
cent.df$PRNAME[which(cent.df$PRNAME=="Nova Scotia")]<-"Nova\nScotia"
cent.df$PRNAME[which(cent.df$PRNAME=="New Brunswick")]<-"New\nBrunswick"
cent.df$PRNAME[which(cent.df$PRNAME=="Prince Edward Island")]<-"PEI"

#make a lower lat for the Inc
cent.df$latlow<-cent.df$lat-200000
cent.df$longlow<-cent.df$long
cent.df$verbose<-paste(cent.df$Inc," cases\nper 100,000",sep="")

#change position of cum.incidence
cent.df$latlow[which(cent.df$PRNAME=="British\nColumbia")]<-cent.df$latlow[which(cent.df$PRNAME=="British\nColumbia")]-30000
cent.df$latlow[which(cent.df$PRNAME=="Northwest\nTerritories")]<-cent.df$latlow[which(cent.df$PRNAME=="Northwest\nTerritories")]-30000
cent.df$latlow[which(cent.df$PRNAME=="PEI")]<-cent.df$latlow[which(cent.df$PRNAME=="PEI")]+200000
cent.df$longlow[which(cent.df$PRNAME=="PEI")]<-cent.df$longlow[which(cent.df$PRNAME=="PEI")]+320000
cent.df$latlow[which(cent.df$PRNAME=="Newfoundland\nand Labrador")]<-cent.df$latlow[which(cent.df$PRNAME=="Newfoundland\nand Labrador")]-30000
cent.df$longlow[which(cent.df$PRNAME=="Nova\nScotia")] <-cent.df$longlow[which(cent.df$PRNAME=="Nova\nScotia")]+150000
cent.df$latlow[which(cent.df$PRNAME=="Nova\nScotia")]<-cent.df$latlow[which(cent.df$PRNAME=="Nova\nScotia")]-30000
cent.df$longlow[which(cent.df$PRNAME=="New\nBrunswick")]<- cent.df$longlow[which(cent.df$PRNAME=="New\nBrunswick")]-30000
cent.df$latlow[which(cent.df$PRNAME=="New\nBrunswick")]<-cent.df$latlow[which(cent.df$PRNAME=="New\nBrunswick")]-30000

#boxes 
#for north arrow
bb<-data.frame(long = c(can.lim.x[1],can.lim.x[2]*0.965), lat = c(can.lim.y[1],can.lim.y[2]*0.96)) #north arrow
bb2<-data.frame(long = can.lim.x, lat = c(can.lim.y[1]*0.94,can.lim.y[2])) #scalebar=

#nowplotit: Inc heatmap
p<-ggplot(data=can.fort1,aes(x=long, y=lat, group=group))+
        geom_polygon(data=can.fort1,fill="grey65",colour="grey35",lwd=0.1)+
        geom_polygon(data=us.fort1,fill="grey85",colour="grey35",lwd=0.1)+
        geom_polygon(data=can.fort1, aes(fill=PRENAME,group = group),colour="grey35", lwd=0.1) +
        GlobFillScale+
        geom_text(data=cent.df,aes(x=long, y=lat,
                                   label=PRNAME,group=NULL),size=2.5,fontface="bold",lineheight = .7)+
        geom_label(data=cent.df,aes(label=verbose,x=longlow,y=latlow,group=NULL),size=2)+
        labs(x = NULL, y=NULL,fill = NULL) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = 'lightcyan1'), axis.line = element_blank(),
              axis.text.x=element_blank(), axis.ticks.x=element_blank(),
              axis.text.y=element_blank(), axis.ticks.y=element_blank(),
              legend.position = "none") +
        coord_fixed()+
        coord_cartesian(xlim =c(can.lim.x),ylim=c(can.lim.y))+
  #alaska
  annotate(geom="text",label="USA",y=(bbox(can)[2,2]-800000), x=(bbox(can)[1,1]+400000),size=4)+
  #mainland US
  annotate(geom="text",label="USA",y=(bbox(can)[2,1])+100000, x=mean(bbox(can)[1,])-800000,size=4)+
  #Canadia
  annotate(geom="text",label="CANADA",y=mean(bbox(can)[2,])-300000, x=mean(bbox(can)[1,])-800000,size=4)
  p2<-p+north(bb,symbol=1,scale=0.08, location="topright")
  p3<-p2+scalebar(data=bb2,dist = 500, dist_unit = "km",
                  st.size=2, st.dist=0.02,
                  transform = FALSE, model = "NAD83",
                  location="topright", inherit.aes =FALSE, border.size=0.3)
ggsave(p3,file=paste("maps_out/","Canada_CumulativeIncidence",".png",sep=""),width=xx/4,height=yy/4,units="in")
