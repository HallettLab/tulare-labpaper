## Set your datpath!! (file in Data_cleaning)
library(tidyverse)
library(readr)
library(vegan)
library(MASS)
library(dplyr)
library(RVAideMemoire) #for posthoc tests on permanova
library(indicspecies)
library(gridExtra)
library(goeveg)#for scree plot of NMDS to test number of dimensions
library(ggpubr)

## Set ggplot2 theme
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 14),
              strip.text= element_text(size = 12), 
              axis.text = element_text(size = 10))

# Import csv file, transform to wide data
all.dat <- read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep=""))
all.dat<-as.data.frame(all.dat)
all.dat<-all.dat %>% dplyr::select(-status, -type)

#import env data (litter, cowpies, rocks, etc)
env <- read_csv(paste(datpath_clean, "/envdat.csv", sep=""))
all.dat<- right_join(all.dat,env, by=c("quadratNew", "year"))
all.dat<- all.dat %>% group_by(quadratNew) %>% filter(max(rock)<79.9) #group by plot, then remove plot based on if max rock is over 80%

#create wide data, first filter so year is 2005 to 2012
all.data<- all.dat %>% dplyr::select(-X1.x, -spcode, -bare,-rock,-litter,-cowpie,-gopher) %>% filter(year>2003 & year < 2013) %>% spread(spname, cover)
all.data[is.na(all.data)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)
levels(as.factor(all.dat$quadratNew))
levels(as.factor(all.dat$thermal))
levels(as.factor(all.dat$spname)) #any to remove?
levels(as.factor(all.dat$burn))
levels(as.factor(all.dat$graze))

#############################
#Test effects of burning#####
#########MODERATE############
#######2005-2008#############

#create wide data, first filter so year is 2005 to 2008
mod.data.early<- all.dat %>% dplyr::select(-X1.x, -spcode, -bare, -rock, -litter, -cowpie, -gopher) %>% filter(year>2003 & year < 2009) %>% filter(thermal=="moderate") %>% spread(spname, cover)
mod.data.early[is.na(mod.data.early)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

#check count of thermal factor
mod.data.early %>% 
  group_by(graze, burn, transect) %>%
  summarise(no_rows = length(thermal)) #note PGEM transects have very little data, they are only in 2006

#remove transects that do not span all years
mod.data.early <- mod.data.early %>% filter(transect!="PGEM1", transect!="PGEM2", transect!="PGEM3",transect!="PGEM4")

#combine Festuca myuros and Festuca bromoides
mod.data.early <- mod.data.early %>% group_by(year, quadratNew) %>% mutate('Festuca sp.'=sum(`Festuca bromoides`, `Festuca myuros`)) %>%
  dplyr::select(-`Festuca bromoides`,-`Festuca myuros`) %>% ungroup()

#see which plots were removed due to high rock cover:
mod.check<-mod.data.early %>% 
  group_by(quadratNew) %>%
  summarise(no_years = length(year))
#TMH1-1
#TMH1-10
#THM3-3

#check count of graze and burn factor
mod.data.early %>% 
  group_by(graze, burn) %>%
  summarise(no_rows = length(graze))

#see which species are not present using colmax
cover.colmax<-sapply(mod.data.early ,max)
i<-cover.colmax != 0
mod.data.early.nozero<-mod.data.early[, i]

plotnames<-mod.data.early.nozero[,1]
cover.mod.early<- mod.data.early.nozero %>% dplyr::select(-c(1:7), -Unknown) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.mod.early)<-plotnames

#check for empty rows
cover.Biodrop.mod.early<-cover.mod.early[rowSums(cover.mod.early[, (1:66)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

#if needed, relativize by row or column or calculate presence/absence
cover.rowsums.me <- rowSums(cover.mod.early [1:66])
cover.relrow.me <- data.frame(cover.mod.early /cover.rowsums.me)
#cover.colmax<-sapply(cover.mod.early ,max)
#cover.relcolmax <- data.frame(sweep(cover.Bio ,2,cover.colmax,'/'))
#cover.pa <- cover.Bio %>% mutate_each(funs(ifelse(.>0,1,0)), 1:57)
mod.data.early$covercheck<-cover.rowsums.me#remove plots with very low cover?

######################
#NMDS FOR 
######################
#make bray-curtis dissimilarity matrix
mod.bcd.early <- vegdist(cover.relrow.me)
dimcheckMDS(cover.relrow.me, distance = "bray", k = 8, trymax = 20, autotransform = F) #check for optimal dimensions - choose when starts to flatten out
mod.mds.early<-metaMDS(cover.relrow.me, distance="bray", trace = TRUE, noshare=0.02, autotransform=F, trymax=1000, k=5) #runs several with different starting configurations
mod.mds.early #solution did not converge after 100 tries, try 1000 more runs
#mod.mds.early<-metaMDS(cover.relrow.me, distance="bray", previous.best = mod.mds.early, noshare=0.02, trace = TRUE, autotransform=T, trymax=1000, k=8)
#mod.mds.early<-metaMDS(cover.relrow.me, distance="bray", previous.best = mod.mds.early, noshare=0.02, trace = TRUE, autotransform=T, trymax=5000, k=8)
#mod.mds.early<-metaMDS(cover.relrow.me, distance="bray", previous.best = mod.mds.early, noshare=0.02, trace = TRUE, autotransform=T, trymax=10000, k=8)
summary(mod.mds.early)

#quick plot of results
stressplot(mod.mds.early, mod.bcd.early) #stressplot to show fit
ordiplot(mod.mds.early)

#store scores in new dataframe
mod.data.early <- mod.data.early %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
spscores1.mod.e<-scores(mod.mds.early,display="sites",choices=1)*-1 #note multiplication by -1 to make this ordination match orientation of grazing ordination
spscores2.mod.e<-scores(mod.mds.early,display="sites",choices=2)*-1 #note multiplication by -1 to make this ordination match orientation of grazing ordination
year<-mod.data.early$year
treatment<-mod.data.early$treatment
spscoresall.mod.e<-data.frame(year,treatment,spscores1.mod.e,spscores2.mod.e)

############################
#Indicator species for early MODERATE#
############################
mod.trt.early <- mod.data.early %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
mod_trt_isa_early = multipatt(cover.relrow.me, mod.trt.early$treatment, control=how(nperm=999))
summary(mod_trt_isa_early)

#create plot in ggplot 
fig1b<-ggplot(subset(spscoresall.mod.e, year==2005), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("b) 2005: 1 year post fire")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="none")
fig1b

fig1c<-ggplot(subset(spscoresall.mod.e, year==2006), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("c) 2006: 2 years post fire")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="none")
fig1c

fig1d<-ggplot(subset(spscoresall.mod.e, year==2007), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("d) 2007: 3 years post fire")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
fig1d

fig1e<-ggplot(subset(spscoresall.mod.e, year==2008), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("e) 2008: 4 years post fire")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
#theme(legend.position="bottom", legend.title=element_text(size=11), legend.text=element_text(size=10), axis.text=element_text(size=8), axis.title=element_text(size=11))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig1e

fig1f<-ggplot(spscoresall.mod.e, aes(x=NMDS1, y=NMDS2, fill=mod.data.early$treatment))+
  geom_point(cex=4.5, pch=21)+
  ggtitle("f) 2005 - 2008")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  #theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="bottom", legend.title=element_text(size=20), legend.text=element_text(size=17))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig1f

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(fig1f)

##test for differences in treatment by year#######
mod.2005<-subset(mod.data.early, year==2005)
cover.2005<-mod.2005 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.05 <- rowSums(cover.2005 [1:156])
cover.relrow.05 <- data.frame(cover.2005/cover.rowsums.05)
mod.bcd.05 <- vegdist(cover.relrow.05)
permanova05<-adonis(cover.relrow.05~mod.2005$treatment,perm=1000, method="bray")
permanova05
pairwise.perm.manova(mod.bcd.05,mod.2005$treatment, nperm=1000) #all three treatments differ in 2005
mod_isa_05 = multipatt(cover.relrow.05, mod.2005$treatment, control=how(nperm=999))
summary(mod_isa_05)

mod.2006<-subset(mod.data.early, year==2006)
cover.2006<-mod.2006 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.06 <- rowSums(cover.2006 [1:156])
cover.relrow.06 <- data.frame(cover.2006/cover.rowsums.06)
mod.bcd.06 <- vegdist(cover.relrow.06)
permanova06<-adonis(cover.relrow.06~mod.2006$treatment, perm=1000, method="bray")
permanova06
pairwise.perm.manova(mod.bcd.06,mod.2006$treatment, nperm=1000) #all three treatments differ in 2006
mod_isa_06 = multipatt(cover.relrow.06, mod.2006$treatment, control=how(nperm=999))
summary(mod_isa_06)

mod.2007<-subset(mod.data.early, year==2007)
cover.2007<-mod.2007 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.07 <- rowSums(cover.2007 [1:156])
cover.relrow.07 <- data.frame(cover.2007/cover.rowsums.07)
mod.bcd.07 <- vegdist(cover.relrow.07)
permanova07<-adonis(cover.relrow.07~mod.2007$treatment, perm=1000, method="bray")
permanova07
pairwise.perm.manova(mod.bcd.07,mod.2007$treatment, nperm=1000) #ungrazed communities are same, both differ from grazed
mod_isa_07 = multipatt(cover.relrow.07, mod.2007$treatment, control=how(nperm=999))
summary(mod_isa_07)

mod.2008<-subset(mod.data.early, year==2008)
cover.2008<-mod.2008 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.08 <- rowSums(cover.2008 [1:156])
cover.relrow.08 <- data.frame(cover.2008/cover.rowsums.08)
mod.bcd.08 <- vegdist(cover.relrow.08)
permanova08<-adonis(cover.relrow.08~mod.2008$treatment, perm=1000, method="bray")
permanova08
pairwise.perm.manova(mod.bcd.08,mod.2008$treatment, nperm=1000) #ungrazed communities are same, both differ from grazed
mod_isa_08 = multipatt(cover.relrow.08, mod.2008$treatment, control=how(nperm=999))
summary(mod_isa_08)

#####################
#successional vectors on summarized MODERATE data for burn (2005-2008)
####################
mod_yr_burn<-all.dat %>%filter(transect!="PGEM1", transect!="PGEM2", transect!="PGEM3",transect!="PGEM4")%>%dplyr::group_by(thermal, burn, graze, year, spname) %>% filter(thermal == "moderate", year>2003 & year<2009) %>% 
  summarize(mean=mean(cover))%>% arrange(burn)%>%  arrange(graze)%>% arrange(year)%>%
  spread(spname, mean) %>% dplyr::select(-Unknown)
mod_yr_burn <- mod_yr_burn %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
mod_yr_burn[is.na(mod_yr_burn)] <- 0 

#combine Festuca myuros and Festuca bromoides
mod_yr_burn <- mod_yr_burn %>% group_by(year, treatment) %>% mutate('Festuca sp.'=sum(`Festuca bromoides`, `Festuca myuros`)) %>%
  dplyr::select(-`Festuca bromoides`,-`Festuca myuros`) %>% ungroup()

cover.yr <- mod_yr_burn %>% ungroup %>% dplyr::select(-thermal,-burn,-graze, -year, -treatment)
cover.rowsums <- rowSums(cover.yr [1:154])
cover.relrow <- data.frame(cover.yr /cover.rowsums)

#merge env to match full data
#yr.env<-merge(env,mod_yr) %>% dplyr::select(NH4, NO3, totalN, ppt, temp)

#make bray-curtis dissimilarity matrix
vec.bcd <- vegdist(cover.relrow)
dimcheckMDS(cover.relrow)#check for optimal dimensions

#NMDS 
vec.mds<-metaMDS(cover.relrow, distance="bray", trace = TRUE, autotransform = T, noshare=0.02, trymax=1000, k=4) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
vec.mds #solution converged after 20 tries
summary(vec.mds)

#quick plot of results
stressplot(vec.mds, vec.bcd) #stressplot to show fit
ordiplot(vec.mds)

#store scores in new dataframe
spscores1.vec<-scores(vec.mds,display="sites",choices=1)
spscores2.vec<-scores(vec.mds,display="sites",choices=2)
year<-mod_yr_burn$year
burn<-mod_yr_burn$burn
spscoresall.vec<-data.frame(burn,year,spscores1.vec,spscores2.vec)


#to plot indicator species (from above) on plot
species.e<-as.data.frame(vec.mds$species)
species.e$name<-row.names(species.e)
spc.e<- species.e %>% filter(name == "Agoseris.heterophylla"| name=="Microseris.douglasii"| name == "Chlorogalum.pomeridianum" |name=="Epilobium.sp."| name=="Muilla.maritima"|name=="Lasthenia.californica"| name=="Calandrinia.ciliata" | name=="Trifolium.depauperatum" | name=="Gilia.tricolor" | name=="Plantago.erecta"| name=="Lepidium.nitidum"|name=="Aphanes.occidentalis"|name=="Castilleja.densiflora"|name=="Brodiaea.spp."|name=="Hemizonia.congesta"|name=="Acmispon.wrangelianus"| name=="Festuca.sp."| name =="Hordeum.murinum.ssp..leporinum" | name == "Festuca.perennis" |name == "Avena.sp.")
spc.e <- spc.e %>% mutate(fontface = "bold", fontface = ifelse(name=="Festuca.sp."| name =="Hordeum.murinum.ssp..leporinum" | name == "Festuca.perennis" |name == "Avena.sp.", "italic",  fontface)) #fontface based native/non-native
#spc.e<- species.e %>% filter(name == "Trifolium.depauperatum"| name == "Rigiopappus.leptoclodus" |name=="Poa.secunda.ssp..secunda"| name=="Festuca.bromoides"|name=="Koeleria.macrantha"| name=="Galium.aparine"| name=="Festuca.myuros"| name=="Layia.gaillardiodes"| name=="Athysanus.pusilus"|name=="Sisyrinchium.bellum"|name=="Epilobium.sp."| name=="Chlorogalum.pomeridianum"|name=="Sanicula.bipinnatifida"|name=="Lessingia.micradenia.glabratai"|name=="Triteleia.laxa"| name=="Allium.serra"|name=="Plantago.erecta"|name=="Lasthenia.californica"|name=="Aphanes.occidentalis"|name=="Erodium.cicutarium"|name=="Gilia.tricolor"|name=="Lepidium.nitidum"| name=="Hemizonia.congesta"| name=="Castilleja.densiflora"| name=="Microseris.douglasii"|name=="Calandrinia.ciliata"| name=="Agoseris.heterophylla"|name=="Bromus.madritensis"|name=="Brodiaea.spp."|name=="Hordeum.murinum ssp..leporinum"|name=="Muilla.maritima"|name=="Festuca.perennis")
#spc.e<- species.e %>% filter(name == "Trifolium depauperatum"| name == "Rigiopappus leptoclodus" |name=="Poa secunda ssp. secunda"| name=="Festuca bromoides"|name=="Koeleria macrantha"| name=="Galium aparine"| name=="Festuca myuros"| name=="Layia gaillardiodes"| name=="Athysanus pusilus"|name=="Sisyrinchium bellum"|name=="Epilobium sp."| name=="Chlorogalum pomeridianum"|name=="Sanicula bipinnatifida"|name=="Lessingia micradenia glabratai"|name=="Triteleia laxa"| name=="Allium serra"|name=="Plantago erecta"|name=="Lasthenia californica"|name=="Aphanes occidentalis"|name=="Erodium cicutarium"|name=="Gilia tricolor"|name=="Lepidium nitidum"| name=="Hemizonia congesta"| name=="Castilleja densiflora"| name=="Microseris douglasii"|name=="Calandrinia ciliata"|name=="Agoseris heterophylla"|name=="Bromus.madritensis"|name=="Brodiaea spp."|name=="Hordeum murinum ssp. leporinum"|name=="Muilla maritima"|name=="Festuca perennis")
#move species around a little so don't overlap on plot
spc.e$MDS2[spc.e$name == "Festuca.sp."] <- -0.2
spc.e$MDS1[spc.e$name == "Agoseris.heterophylla"] <- -0.7
spc.e$MDS1[spc.e$name == "Lepidium.nitidum"] <- -0.35
spc.e$MDS2[spc.e$name == "Plantago.erecta."] <- 0.06
spc.e$MDS1[spc.e$name == "Calandrinia.ciliata"] <- -0.3
spc.e$MDS1[spc.e$name == "Lasthenia.californica"] <- -0.78
spc.e$MDS1[spc.e$name == "Castilleja.densiflora"] <- -0.65
spc.e$MDS2[spc.e$name == "Aphanes.occidentalis"] <- 0.35
spc.e$MDS2[spc.e$name == "Acmispon.wrangelianus"] <- 0.29
spc.e$MDS2[spc.e$name == "Hemizonia.congesta"] <- 0.55
spc.e$MDS1[spc.e$name == "Chlorogalum.pomeridianum"] <- 0.25
spc.e$MDS2[spc.e$name == "Microseris.douglasii"] <- 0.7

vec1<-ggplot(spscoresall.vec, aes(x=NMDS1, y=NMDS2))+
  xlim(-1,1)+
  ggtitle("a) Change in community composition after fire from 2005 - 2008")+
  geom_path(arrow=arrow(), aes(col=mod_yr_burn$treatment))+
  scale_color_manual(values=c("grey0", "grey60","grey85"))+
  geom_point(cex=6, aes(fill=mod_yr_burn$treatment), pch=21)+
  scale_fill_manual(values=c("grey0", "grey85","grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+

  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
#theme(legend.position=c(0.8,0.8), legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
vec1

vec1<-vec1+geom_text(data=spc.e, mapping=aes(x=MDS1, y=MDS2, label=name, fontface=fontface), cex=3)
vec1

#create layout for panel
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,3,3),
             c(2,2,3,3),
             c(2,2,3,3),
             c(4,4,5,5),
             c(4,4,5,5),
             c(4,4,5,5))
grid.arrange(vec1, fig1b, fig1c, fig1d, fig1e, layout_matrix = lay) #put panel together
#save as 900W x 1200H
#save PDF as 8x14" portrait


################
##Test effects of grazing reintroduction
###MODERATE#####
###late years####
#create wide data, first filter so year is 2008 to 2012
mod.data.late<- all.dat %>% dplyr::select(-X1.x, -X1.y, -gopher, -bare, -rock, -litter, -cowpie, -spcode) %>% filter(year>2007 & year < 2013) %>% filter(thermal=="moderate") %>% spread(spname, cover)
mod.data.late[is.na(mod.data.late)] <- 0 #replace NAs with 0 (species not counted in plots have NAs when wide dataset created)

#remove transects that do not span all years
mod.data.late <- mod.data.late %>% filter(transect!="PGEM1", transect!="PGEM2", transect!="PGEM3",transect!="PGEM4")

#combine Festuca myuros and Festuca bromoides
mod.data.late <- mod.data.late %>% group_by(year, quadratNew) %>% mutate('Festuca sp.'=sum(`Festuca bromoides`, `Festuca myuros`)) %>%
  dplyr::select(-`Festuca bromoides`,-`Festuca myuros`) %>% ungroup()

#check count of graze and burn factor
mod.data.late %>% 
  group_by(graze, burn, transect) %>%
  summarise(no_rows = length(graze))

plotnames<-mod.data.late[,1]
cover.mod.late<- mod.data.late %>% dplyr::select(-c(1:6)) #wide data with ID columns removed, only species/cover for NMDS
rownames(cover.mod.late)<-plotnames

#now relativize by row
cover.rowsums.ml <- rowSums(cover.mod.late [1:155])
cover.relrow.ml <- data.frame(cover.mod.late /cover.rowsums.ml)

#check for empty rows
cover.Biodrop.mod.late<-cover.mod.late[rowSums(cover.mod.late[, (1:155)]) ==0, ] #if no empty rows, next step not needed
#cover.Biodrop.gb<-cover.gb[rowSums(cover.gb[, (1:157)])  >0 ]#remove empty rows

######################
#NMDS
######################
#make bray-curtis dissimilarity matrix
mod.bcd.late <- vegdist(cover.relrow.ml)

#check for optimal number of dimensions (optimum is when scree plot flattens out/break in slope, we'd also like the stress to be under 10% or 0.10)
dimcheckMDS(cover.relrow.ml, distance = "bray", k = 8, trymax = 20, autotransform = F)

mod.mds.late<-metaMDS(cover.relrow.ml, trace = TRUE, autotransform=F, noshare=0.01, trymax=1000, k=5) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, k=4 for 4 dimensions
mod.mds.late #no solution reached
summary(mod.mds.late)

#quick plot of results
stressplot(mod.mds.late, mod.bcd.late) #stressplot to show fit
ordiplot(mod.mds.late)
orditorp(mod.mds.late,display="sites",cex=0.8,air=0.01)

#store scores in new dataframe
mod.trt.late <- mod.data.late %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
spscores1.mod.l<-scores(mod.mds.late,display="sites",choices=1)
spscores2.mod.l<-scores(mod.mds.late,display="sites",choices=2)
treatment<-mod.trt.late$treatment
year<-mod.data.late$year
spscoresall.mod.l<-data.frame(year,treatment,spscores1.mod.l,spscores2.mod.l)

############################
#Indicator species for late MODERATE#
############################
mod_trt_isa_late = multipatt(cover.relrow.ml, mod.trt.late$treatment, control=how(nperm=999))
summary(mod_trt_isa_late)

#create plot in ggplot 
fig2b<-ggplot(subset(spscoresall.mod.l, year==2008), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("b) 2008: grazing reintroduced")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="none")
fig2b

fig2c<-ggplot(subset(spscoresall.mod.l, year==2009), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("c) 2009: 1 year post-grazing")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="none")
fig2c
#fig1c+geom_text(subset(spscoresall.mod.e, year==2009), mapping=aes(x=NMDS1, y=NMDS2, label=mod.2009$transect), cex=4)

fig2d<-ggplot(subset(spscoresall.mod.l, year==2010), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("d) 2010: 2 years post-grazing")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
fig2d

fig2e<-ggplot(subset(spscoresall.mod.l, year==2011), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("e) 2011: 3 years post-grazing")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
#theme(legend.position="bottom", legend.title=element_text(size=11), legend.text=element_text(size=10), axis.text=element_text(size=8), axis.title=element_text(size=11))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig2e

fig2f<-ggplot(subset(spscoresall.mod.l, year==2012), aes(x=NMDS1, y=NMDS2, fill=treatment, shape=as.factor(year)))+
  geom_point(cex=4.5)+
  ggtitle("f) 2012: 4 years post-grazing")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_shape_manual(values=c(21),guide = guide_legend(title = "Year"))+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))
#theme(legend.position="bottom", legend.title=element_text(size=11), legend.text=element_text(size=10), axis.text=element_text(size=8), axis.title=element_text(size=11))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig2f

fig2g<-ggplot(spscoresall.mod.l, aes(x=NMDS1, y=NMDS2, fill=spscoresall.mod.l$treatment, shape=as.factor(mod.data.late$year)))+
  geom_point(cex=4.5,pch = 21)+
  ggtitle("g)")+
  xlim(-1,1)+
  ylim(-1,1)+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment:"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  #theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  theme(legend.position="left", legend.title=element_text(size=20), legend.text=element_text(size=17), axis.text=element_text(size=6), axis.title=element_text(size=11))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
fig2g

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend2<-g_legend(fig2g)

##test for differences in treatment by year#######
mod.2008<-subset(mod.data.late, year==2008) %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
cover.2008<-mod.2008 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.08 <- rowSums(cover.2008 [1:155])
cover.relrow.08 <- data.frame(cover.2008/cover.rowsums.08)
mod.bcd.08 <- vegdist(cover.relrow.08)
permanova08<-adonis(cover.relrow.08~mod.2008$treatment,perm=1000, method="bray")
permanova08
pairwise.perm.manova(mod.bcd.08,mod.2008$treatment, nperm=1000) #ungrazed communities are the same
mod_isa_08 = multipatt(cover.relrow.08, mod.2008$treatment, control=how(nperm=999))
summary(mod_isa_08)

mod.2009<-subset(mod.data.late, year==2009) %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
cover.2009<-mod.2009 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.09 <- rowSums(cover.2009 [1:155])
cover.relrow.09 <- data.frame(cover.2009/cover.rowsums.09)
mod.bcd.09 <- vegdist(cover.relrow.09)
permanova09<-adonis(cover.relrow.09~mod.2009$treatment, perm=1000, method="bray")
permanova09
pairwise.perm.manova(mod.bcd.09,mod.2009$treatment, nperm=1000) #all communities differ
mod_isa_09 = multipatt(cover.relrow.09, mod.2009$treatment, control=how(nperm=999))
summary(mod_isa_09)

mod.2010<-subset(mod.data.late, year==2010) %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
cover.2010<-mod.2010 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.10 <- rowSums(cover.2010 [1:155])
cover.relrow.10 <- data.frame(cover.2010/cover.rowsums.10)
mod.bcd.10 <- vegdist(cover.relrow.10)
permanova10<-adonis(cover.relrow.10~mod.2010$treatment, perm=1000, method="bray")
permanova10
pairwise.perm.manova(mod.bcd.10,mod.2010$treatment, nperm=1000) #all communities differ
mod_isa_10 = multipatt(cover.relrow.10, mod.2010$treatment, control=how(nperm=999))
summary(mod_isa_10)

mod.2011<-subset(mod.data.late, year==2011) %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
cover.2011<-mod.2011 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.11 <- rowSums(cover.2011 [1:155])
cover.relrow.11 <- data.frame(cover.2011/cover.rowsums.11)
mod.bcd.11 <- vegdist(cover.relrow.11)
permanova11<-adonis(cover.relrow.11~mod.2011$treatment, perm=1000, method="bray")
permanova11
pairwise.perm.manova(mod.bcd.11,mod.2011$treatment, nperm=1000) #all communities differ
mod_isa_11 = multipatt(cover.relrow.11, mod.2011$treatment, control=how(nperm=999))
summary(mod_isa_11)

mod.2012<-subset(mod.data.late, year==2012) %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")
cover.2012<-mod.2012 %>% dplyr::select(-quadratNew,-treatment,-thermal,-burn,-graze, -year, -transect)
cover.rowsums.12 <- rowSums(cover.2012 [1:155])
cover.relrow.12 <- data.frame(cover.2012/cover.rowsums.12)
mod.bcd.12 <- vegdist(cover.relrow.12) 
permanova12<-adonis(cover.relrow.12~mod.2012$treatment, perm=1000, method="bray")
permanova12
pairwise.perm.manova(mod.bcd.12,mod.2012$treatment, nperm=1000) #all communities differ
mod_isa_12 = multipatt(cover.relrow.12, mod.2012$treatment, control=how(nperm=999))
summary(mod_isa_12)

#####################
#successional vectors on summarized MODERATE data for graze (2008-2012)
####################
mod_yr_graze<-all.dat %>% filter(transect!="PGEM1", transect!="PGEM2", transect!="PGEM3",transect!="PGEM4")%>%dplyr::group_by(thermal, burn, graze, year, spname) %>% filter(thermal == "moderate", year>2007 & year<2013) %>% dplyr::select(-X1.x, -X1.y, -gopher, -bare, -rock, -litter, -cowpie)%>%
  summarize(mean=mean(cover))%>% arrange(burn)%>%  arrange(graze)%>% arrange(year)%>%
  spread(spname, mean) 
mod_yr_graze[is.na(mod_yr_graze)] <- 0 

#combine Festuca myuros and Festuca bromoides
mod_yr_graze <- mod_yr_graze %>% group_by(year, burn, graze) %>% mutate('Festuca sp.'=sum(`Festuca bromoides`, `Festuca myuros`)) %>%
  dplyr::select(-`Festuca bromoides`,-`Festuca myuros`) %>% ungroup()

cover.yr <- mod_yr_graze %>% ungroup %>% dplyr::select(-thermal,-burn,-graze, -year)
cover.rowsums <- rowSums(cover.yr [1:155])
cover.relrow <- data.frame(cover.yr /cover.rowsums)

#make bray-curtis dissimilarity matrix
vec.bcd <- vegdist(cover.relrow)
dimcheckMDS(cover.relrow)#check for optimal dimensions
#starts to flatten out around 5 or 6

#NMDS 
vec.mds<-metaMDS(cover.relrow, distance="bray", trace = TRUE, autotransform = T, noshare=0.02, trymax=1000, k=4) #runs several with different starting configurations
#trace= TRUE will give output for step by step what its doing
#default is 2 dimensions, can put k=4 for 4 dimensions
vec.mds #solution converged after 20 tries
summary(vec.mds)

#quick plot of results
stressplot(vec.mds, vec.bcd) #stressplot to show fit
ordiplot(vec.mds)

#store scores in new dataframe
spscores1.vec<-scores(vec.mds,display="sites",choices=1)
spscores2.vec<-scores(vec.mds,display="sites",choices=2)
year<-mod_yr_graze$year
burn<-mod_yr_graze$burn
spscoresall.vec<-data.frame(burn,year,spscores1.vec,spscores2.vec)

#make plot to show successional vectors
#first, set colors and shapes
mod_yr_graze <- mod_yr_graze %>%
  unite(treatment, c(burn, graze), remove=FALSE, sep = " ")

############################
###Make vector plot again in ggplot
###########################

#to plot indicator species (from above) on plot
species.e<-as.data.frame(vec.mds$species)
species.e$name<-row.names(species.e)
spc.e<- species.e %>% filter(name == "Agoseris.heterophylla"| name=="Microseris.douglasii"| name == "Chlorogalum.pomeridianum" |name=="Epilobium.sp."| name=="Muilla.maritima"|name=="Lasthenia.californica"| name=="Calandrinia.ciliata" | name=="Trifolium.depauperatum" | name=="Gilia.tricolor" | name=="Plantago.erecta"| name=="Lepidium.nitidum"|name=="Aphanes.occidentalis"|name=="Castilleja.densiflora"|name=="Brodiaea.spp."|name=="Hemizonia.congesta"|name=="Acmispon.wrangelianus"| name=="Festuca.sp."| name =="Hordeum.murinum.ssp..leporinum" | name == "Festuca.perennis" |name == "Avena.sp.")
spc.e <- spc.e %>% mutate(fontface = "bold", fontface = ifelse(name=="Festuca.sp."| name =="Hordeum.murinum.ssp..leporinum" | name == "Festuca.perennis" |name == "Avena.sp.", "italic",  fontface)) #fontface based native/non-native
#move label over so not outside of plot margins
spc.e$MDS1[spc.e$name == "Agoseris.heterophylla"] <- -0.7
spc.e$MDS1[spc.e$name == "Chlorogalum.pomeridianum"] <- 0.4

vec2<-ggplot(spscoresall.vec, aes(x=NMDS1, y=NMDS2))+
  geom_path(arrow=arrow(), aes(col=mod_yr_graze$treatment))+
  scale_color_manual(values = c("grey0", "grey60", "grey85")) +
  geom_point(cex=6, pch = 21, aes(fill=mod_yr_graze$treatment))+
  ggtitle("a) Change in community composition from 2008-2012 in response to grazing")+
  scale_fill_manual(values=c("grey0", "grey85", "grey100"), guide = guide_legend(title = "Treatment"), #change legend title
                     labels=c("Burned & Grazed", "Burned & Ungrazed", "Unburned & Ungrazed"))+ #change labels in the legend)+
  theme(legend.position="none")+
  theme(plot.title = element_text(color="black", face="bold.italic"))+
  xlim(-0.85,0.6)
#theme(legend.position=c(0.8,0.8), legend.title=element_text(size=14), legend.text=element_text(size=12), axis.text=element_text(size=16), axis.title=element_text(size=16))+
#theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
vec2

vec2<-vec2+geom_text(data=spc.e, mapping=aes(x=MDS1, y=MDS2, label=name, fontface=fontface), cex=4)
vec2

#create layout for panel
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,3,3),
             c(2,2,3,3),
             c(2,2,3,3),
             c(4,4,5,5),
             c(4,4,5,5),
             c(4,4,5,5),
             c(6,6,NA,NA),
             c(6,6,NA,NA),
             c(6,6,NA,NA))
grid.arrange(vec2, fig2b, fig2c, fig2d, fig2e, fig2f, layout_matrix = lay) #put panel together
#save as 900W X 1600H png
#save as 8x18.7 PDF

#save as 800w x 1200l