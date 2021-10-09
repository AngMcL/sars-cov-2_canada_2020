#!/usr/bin/env Rscript

# Canada timeline

#usage:
# $Rscript CanadaTimeline.R

# lib setup
library(lubridate)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ggplotify)
library(gridExtra)
library(grid)
library(lattice)
library(reshape2)

# read in the data
cor<-read.csv("resources/covid-interventions.csv")
#read in stringency for oxford
ox<-read.csv("resources/oxford_stringency_index.csv")
ox<-ox[,-1]
ox.can<-ox[which(ox$country_name=="Canada"),]
colnames(ox.can)<-str_replace_all(colnames(ox.can),"X","")

#melt 
ox.can.melt<-melt(ox.can)
colnames(ox.can.melt)[3]<-"date"
colnames(ox.can.melt)[4]<-"stringency"

#convert dates
ox.can.melt$date<-as.Date(ox.can.melt$date,format="%d%B%Y")

# ggplot(ox.can.melt)+
#   geom_line(aes(x=date,y=stringency))

#subset to national
can<-cor[which(cor$Jurisdiction=="Can."),]
can$Date<-as.Date(can$Date)

#order by date
can<-can[with(can,order(Date)),]
rownames(can)<-1:nrow(can)

#keep important rows
can.2<-can[c(1,3:5,8:10,12,14,17,22:23,25,32,34,35,38,42,43,47,50,51,54),]

#empty timeline for 2020
time<-data.frame(date=c(as.Date("2020-01-01"),as.Date("2020-12-31")),y=0.5)

#shorten descriptions
can.2$Description<-str_replace_all(can.2$Description," issued","")
can.2$Description<-str_replace_all(can.2$Description," announced","")
can.2$Description[which(rownames(can.2)==54)]<-"International travellers permitted to leave quarantine with negative test upon arrival and retest 1 wk later"
can.2$Description[which(rownames(can.2)==43)]<-"Canada’s Flight Plan for Navigating COVID-19, aimed at reducing risks of air travel"
can.2$Description[which(rownames(can.2)==51)]<-"Travel restrictions for international students eased"
can.2$Description[which(rownames(can.2)==50)]<-"Travel restrictions for extended family members of Canadian citizens and residents eased"

#shift one day back for viz purposes
# can.2$Date[which(can.2$Description=="Cruise ship season postponed")]<-as.Date("2020-03-12")

#populate a df with important events in covid epidemic
keydates<-c("2020-01-25","2020-02-20","2020-03-05","2020-03-11","2020-03-23","2020-05-04","2020-12-09")
keydescr<-c("First case of novel coronavirus in Canada","First COVID-19 case in Canada from outside mainland China","First Canadian case from community transmission","WHO declares COVID-19 a global pandemic","Federal gov. announces repatriation flights for Canadians stranded abroad","Restrictions relaxed in several provinces","Health Canada approves first vaccine")
#7 events
n.ev<-7
inter.type<-rep("Key Canadian COVID-19 event",n.ev)

#make an empty df the same structure as can.2
key.df<-can.2[1:n.ev,]
key.df[1:n.ev,]<-NA
key.df$Date<-keydates
key.df$Intervention.Category<-inter.type
key.df$Description<-keydescr

#merge key.df onto can.2
can.3<-rbind(can.2,key.df)

#re-sort by date
can.3<-can.3[with(can.3,order(Date)),]

#make a column for daily change in stringency
ox.can.melt$stringency.change<-NA
for (i in 1:nrow(ox.can.melt)){
  if (i>1){
      ox.can.melt$stringency.change[i]<-ox.can.melt$stringency[i]-ox.can.melt$stringency[i-1]
  }
}
# ox.can.melt$date [which(ox.can.melt$stringency.change==max(ox.can.melt$stringency.change,na.rm = T))]
# march 16 was the day that the stringency increased the most

#some pointers from https://benalexkeen.com/creating-a-timeline-graphic-using-r-and-ggplot2/
## Make an improved timeline plot using y positions for text
# resort by date and generate positions and direction for plotting
can.3<-can.3[with(can.3,order(Date)),]
positions <- 1:nrow(can.3)
# positions <- 0.5+0.5*(0:(length(can.4$Date)-1))
directions <- c(1, 1)

line_pos <- data.frame(
    "Date"=unique(can.3$Date),
    "position"=rep(positions, length.out=length(unique(can.3$Date))),
    "direction"=rep(directions, length.out=length(unique(can.3$Date))))

can.4 <- merge(x=can.3, y=line_pos, by="Date", all = TRUE)

#adjust descriptions for conciseness
can.4$Description [24]<-"Entry conditions for travellers transiting through Canada to Alaska"
can.4$Description [29]<-"Quarantine shortened if neg. test upon arrival and retest"
can.4$Description [27]<-"Travel restrictions eased for extended family members"
can.4$Description [30]<- "Health Canada approves first vaccine"
can.4$Description [25]<- "Canada releases Flight Plan for Navigating COVID-19"

# add modified type column with distinction for restriction up /down + added/eased
# use to color text
can.4$Intervention.Category<-str_replace_all(can.4$Intervention.Category,"Key Canadian COVID-19 event","COVID-19 event")
can.4$Intervention.Category<-str_replace_all(can.4$Intervention.Category,"Travel","Restriction added")
can.4$Intervention.Category<-str_replace_all(can.4$Intervention.Category,"Case management","Restriction added")
can.4$Intervention.Category<-str_replace_all(can.4$Intervention.Category,"Public information","Restriction added")
#specify if restrictions eased
can.4$Intervention.Category[c(19,22,24,27,28,29)]<-"Restriction eased"

## rescale positions based on the oxford stringency index range (want a nested figure)
max.ox<-max(ox.can.melt$stringency,na.rm=T)
jump.ox<-max.ox/nrow(can.4)
can.4$position.y<- -16+jump.ox*(1:nrow(can.4)) 

#better colors
ox.color<-'grey20' #dark orange
type.colors<- c('firebrick','darkgreen','dodgerblue4')#red green blue

p1<-ggplot()+
  theme_classic()+
  theme(axis.title.x=element_text(NULL),
        panel.background = element_rect(fill="grey95"),
        axis.text.x=element_text(angle=45,hjust=1),
        axis.title.y=element_text(size=rel(1.25),face="bold"),
        legend.position="none")+
  labs(x=NULL,y=NULL,color="Intervention category")+
  #description and lines to axis
  geom_segment(data=can.4,aes(y=position.y,yend=-16, 
                              x=Date,xend=Date,color=Intervention.Category),
               size=0.2,alpha=0.4)+
  geom_text(data=can.4,aes(x=Date,label=Description,y=position.y,
                           color=Intervention.Category),  size=2.7,hjust=0)+
  
  #Add a manual legend of colored text in a box
    annotate(geom="rect",xmin = as.Date("2020-12-13"), xmax = as.Date("2021-03-01"),
             ymin = -2, ymax =19, color="black",size=0.3,  fill="white")+
    annotate(geom="text",x=as.Date("2020-12-15"),y=15,label="Event type",size=4,fontface="bold",hjust=0,vjust=0,alpha=0.9, color="black")+
  annotate(geom="text",x=as.Date("2020-12-15"),y=10,label="Restriction added",size=4,hjust=0,vjust=0,alpha=0.9, color=type.colors[1])+
    annotate(geom="text",x=as.Date("2020-12-15"),y=5,label="Restriction eased",size=4,hjust=0,vjust=0,alpha=0.9, color=type.colors[2])+
    annotate(geom="text",x=as.Date("2020-12-15"),y=0,label="COVID-19 event",size=4,hjust=0,vjust=0,alpha=0.9, color=type.colors[3])+

  # add on the ox stringency
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y", 
               expand=expansion(add = c(0, 50)),
               limits=c(as.Date("2020-01-01"),as.Date("2021-03-05")),
               labels=c(as.Date("2020-01-01"),as.Date("2021-03-05")),
               breaks=c(as.Date("2020-01-01"),as.Date("2021-03-05")))+
  scale_color_manual(values=type.colors, breaks = c("Restriction added","Restriction eased","COVID-19 event"))+
  scale_y_continuous(expand = expansion(add=c(0, 0)))+
  coord_cartesian(ylim=c(-16,78))+
  labs(x=NULL,y="Oxford Stringency Index")+
  geom_line(data=ox.can.melt,aes(x=date,y=stringency),color=ox.color,alpha=0.5,size=1)
p1
ggsave("results/Canada-COVID-timeline.png",height=4.2,width=8.5,units="in")





