###############################
###Transect variation figure###
###############################

#Run appropriate datpath from data_pathway script in Data Cleaning Folder

library(tidyverse)
library(ggpubr)

#SET UP DATA:

#load updated master data, "alldat", in grazing_recovery.R
alldat<-read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep="")) %>%
  select(-1)%>%
  filter(transect%in%c("THBUGM1", "THBUGM2", "THM1", "THM2", "THM3", "THM4", "THUBUGM1", "THUBUGM2"))%>%
  filter(!quadratNew%in%c("THM1-1", "THM3-3", "THM1-10"))%>%
  filter(thermal=="moderate")%>%
  group_by(year, spname, spcode, quadratNew, status, type, transect, burn, graze)%>%
  summarize(cover=sum(cover))

##########
#richness#
##########
#load "rich" from grazing_recovery.R
rich <- alldat %>%
  filter(type!="NA", status!="NA")%>%
  filter(spname !="Unknown", spname!="Moss") %>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  mutate(present=ifelse(cover>0, 1, 0))%>%
  group_by(year, quadratNew, func, trt, type, status, burn, graze, transect, present)%>%
  summarize(richness = sum(present))%>%
  group_by(year, quadratNew, func, trt, type, status, burn, graze, transect)%>%
  summarize(richness=sum(richness))
#separate transect and quadrat
richrich<-rich%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 
#set the intercept to ungrazed burned for the third pair comparison
richrich2<-richrich%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))%>%
  filter(!is.na(richness))
richrich2$trt <- factor(richrich2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))
#unique quadrat name
richrich2$quadUnique <- paste(richrich2$transect, richrich2$quadrat, sep = "_")
richrich2$year<- factor(richrich2$year, ordered=TRUE)

#######
#cover#
#######
#load "cov" from grazing_recovery.R
cov<-alldat%>%
  filter(type!="NA", status!="NA")%>%
  mutate(func=paste(type, status))%>%
  mutate(trt=paste(graze, burn))%>%
  filter(spname !="Unknown", spname!="Moss") %>%
  group_by(year, quadratNew, trt, func, burn, graze, transect)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(trt), !is.na(func), !is.na(sumcov)) %>%
  group_by(year, quadratNew, trt)%>%
  mutate(totcov=sum(sumcov))%>%
  mutate(relcov=sumcov/totcov)
#separate transect and quadrat
covcov<-cov%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 
#set the intercept to ungrazed burned for the third pair comparison
covcov2<-covcov%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))
covcov2$trt <- factor(covcov2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))
#unique quadrat name
covcov2$quadUnique <- paste(covcov2$transect, covcov2$quadrat, sep = "_")
covcov2$year<- factor(covcov2$year, ordered=TRUE)

theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 10),
              strip.text= element_text(size = 10),
              axis.text = element_text(size = 10))

p1 <- ggplot(data = subset(richrich, func=="grass non-native"&year>2004&year<2013), aes(x=transect, y=richness, fill=trt))+
  geom_boxplot()+
  geom_jitter(aes(color=year))+
  ylab("Richness")+
  theme(axis.text.x = element_text(angle = 90, hjust=1), axis.title.x=element_blank())+
  scale_fill_manual(values=c("grey0","grey85","grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                    labels=c("Burned & Ungrazed", "Unburned & Ungrazed","Burned & Grazed")) #change labels in the legend)+
#facet_grid(~time)

p2 <- ggplot(data = subset(richrich, func=="forb native"&year>2004&year<2013), aes(x=transect, y=richness, fill=trt))+
  geom_boxplot()+
  geom_jitter(aes(color=year))+
  ylab("Richness")+
  theme(axis.text.x = element_text(angle = 90, hjust=1), axis.title.x=element_blank())+
  scale_fill_manual(values=c("grey0","grey85","grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                    labels=c("Burned & Ungrazed", "Unburned & Ungrazed","Burned & Grazed")) #change labels in the legend)+
#facet_grid(~time)

p3 <- ggplot(data = subset(covcov, func=="grass non-native"&year>2004&year<2013), aes(x=transect, y=relcov, fill=trt))+
  geom_boxplot()+
  geom_jitter(aes(color=year))+
  ylab("Cover (%)")+
  theme(axis.text.x = element_text(angle = 90, hjust=1), axis.title.x=element_blank())+
  scale_fill_manual(values=c("grey0","grey85","grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                    labels=c("Burned & Ungrazed", "Unburned & Ungrazed","Burned & Grazed"))


p4<-ggplot(data = subset(covcov, func=="forb native"&year>2004&year<2013), aes(x=transect, y=relcov, fill=trt))+
  geom_boxplot()+
  geom_jitter(aes(color=year))+
  ylab("Cover (%)")+
  theme(axis.text.x = element_text(angle = 90, hjust=1), axis.title.x=element_blank())+
  scale_fill_manual(values=c("grey0","grey85","grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                    labels=c("Burned & Ungrazed", "Unburned & Ungrazed","Burned & Grazed")) #change labels in the legend)+
#facet_grid(~time)

ggarrange(ggparagraph(text=".", size = 16, color = "white"), ggparagraph(text=".", size = 16, color = "white"),p4, p2, p3, p1,  ncol = 2, nrow = 3,
          labels = c("","","a) Native forb cover", "b) Native forb richness",
                     "c) Non-native grass cover", "d) Non-native grass richness"),
          common.legend = TRUE, legend = "right",
          font.label = list(size = 12),
          heights=c(0.15,1.5,1.5),
          hjust = c(0,0,-0.45, -0.35, -0.35, -0.25),
          vjust = c(0.4,0.4,0.4,0.4))