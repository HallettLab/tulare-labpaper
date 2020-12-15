#set up the working directory
datpath_clean <- "~/University of Oregon/O365.hallett-lab - Documents/Tulare-LabPaper/Lab day/"

#load packages
library(tidyverse)
library(nlme) #linear mixed models
library(MuMIn)
library(readr)

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

########
#litter#
########
#load environmental data
envdat <- read_csv(paste(datpath_clean, "/envdat.csv", sep = "")) %>%
  select(-1) 
#load "joindat" from grazing_recovery.R
joindat <- left_join(alldat, envdat)
#load "lit" from grazing_recovery.R
lit <- joindat %>%
  mutate(trt=paste(graze, burn)) %>%
  ungroup() %>%
  dplyr::select(quadratNew, year, trt, transect, burn, graze, litter)
#separate transect and quadrat
litlit<-lit%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 
#set the intercept to ungrazed burned for the third pair comparison
litlit2<-litlit%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))%>%
  filter(!is.na(litter))
litlit2$trt <- factor(litlit2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))
litlit2$year<- factor(litlit2$year, ordered=TRUE)
levels(litlit2$year)
litlit3<-litlit2 %>% group_by(transect2, quadrat,  year, trt) %>% summarize(lit=mean(litter))
#unique quadrat name
litlit3$quadUnique <- paste(litlit3$transect2, litlit3$quadrat, sep = "_")

#MODELS:

#We run repeated measures to native and non-native richness & cover and litter cover.
#Models for 2005-2008 will be separate from models for 2009-2012 because the grazing reintroduction happened in 2008.

#I. No transect or quadrat
#II. Transect as fixed
#III. Transect as random
#III. Quadrat  as random (each quadrat has a unique identifier)

##############################
#nonnative richness 2005-2008#
##############################
rich_IG2<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2004&year<2009))
summary(rich_IG2)
rich_IG2_aov<-anova(rich_IG2)
rich_IG2_aov #treatment significant
AIC(rich_IG2)

#nonnative richness 2005-2008 (transect as fixed)
rich_IG2_fixed<-lm(richness~trt*year+transect, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2004&year<2009))
summary(rich_IG2_fixed)
AIC(rich_IG2_fixed)
rich_IG2_fixed_aov<-anova(rich_IG2_fixed)
rich_IG2_fixed_aov #treatment significant, transect significant
emmeans(rich_IG2_fixed,  ~ transect) #THM2 stands out (low nonnative richness)

#nonnative richness 2005-2008 (transect as random)
rich_IG2_random<-lme(richness~trt*year, random=~1|transect, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2004&year<2009))
summary(rich_IG2_random)
AIC(rich_IG2_random)
rich_IG2_random_aov<-anova(rich_IG2_random)
rich_IG2_random_aov #now treament is not significant
r.squaredGLMM(rich_IG2_random) #adding transect explains an additional ~38% of the variation

#nonnative richness 2005-2008 (quadrat as random)
rich_IG2_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2004&year<2009))
AIC(rich_IG2_random)
rich_IG2_random_aov<-anova(rich_IG2_random)
rich_IG2_random_aov # treatment not significant now
r.squaredGLMM(rich_IG2_random) #similar to transect as a random effect, quad as random explains an additional ~40% of the variation

##############################
#nonnative richness 2009-2012#
##############################
rich_IG3<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2008&year<2013))
summary(rich_IG3)
AIC(rich_IG3)
rich_IG3_aov<-anova(rich_IG3)
rich_IG3_aov #treatment and year are significant

#nonnative richness 2009-2012 (transect as fixed)
rich_IG3_fixed<-lm(richness~trt*year+transect, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2008&year<2013))
summary(rich_IG3_fixed)
AIC(rich_IG3_fixed)
rich_IG3_fixed_aov<-anova(rich_IG3_fixed)
rich_IG3_fixed_aov #treatment, year, treatment*year, and transect are significant
emmeans(rich_IG3_fixed, ~transect) #again, THM2 stands out as low nonnative richness

#nonnative richness 2009-2012 (transect as random)
rich_IG3_random<-lme(richness~trt*year, random=~1|transect, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2008&year<2013))
AIC(rich_IG3_random)
rich_IG3_random_aov<-anova(rich_IG3_random)
rich_IG3_random_aov #year and interaction of trt*year are significant
r.squaredGLMM(rich_IG3_random) #additional 27% of variation explained by transect 

#nonnative richness 2009-2012 (quadrat as random)
rich_IG3_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2008&year<2013))
AIC(rich_IG3_random)
rich_IG3_random_aov<-anova(rich_IG3_random)
rich_IG3_random_aov #year and interaction of trt*year are significant
r.squaredGLMM(rich_IG3_random) #additional 49% of variation explained by quadrat

###########################
#native richness 2005-2008#
###########################
rich_NF2<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2004&year<2009))
summary(rich_NF2)
AIC(rich_NF2)
rich_NF2_aov<-anova(rich_NF2)
rich_NF2_aov #treatment is significant

#native richness 2005-2008 (transect as fixed)
rich_NF2_fixed<-lm(richness~trt*year+transect, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2004&year<2009))
summary(rich_NF2_fixed)
AIC(rich_NF2_fixed)
rich_NF2_fixed_aov<-anova(rich_NF2_fixed)
rich_NF2_fixed_aov #all terms are significant
emmeans(rich_NF2_fixed, ~transect) #transect THM2 stands out as very high native richness

#native richness 2005-2008 (transect as random)
rich_NF2_random<-lme(richness~trt*year, random=~1|transect, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2004&year<2009))
AIC(rich_NF2_random)
rich_NF2_random_aov<-anova(rich_NF2_random)
rich_NF2_random_aov #year and interaction of trt*year are significant
r.squaredGLMM(rich_NF2_random) #explains an additional ~30% of the variation

#native richness 2005-2008 (quadrat as random)
rich_NF2_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2004&year<2009))
AIC(rich_NF2_random)
rich_NF2_random_aov<-anova(rich_NF2_random)
rich_NF2_random_aov
r.squaredGLMM(rich_NF2_random)

###########################
#native richness 2009-2012#
###########################
rich_NF3<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2008&year<2013))
summary(rich_NF3)
AIC(rich_NF3)
rich_NF3_aov<-anova(rich_NF3)
rich_NF3_aov #trt, year are significant

#native richness 2009-2012 (transect as fixed)
rich_NF3_fixed<-lm(richness~trt*year+transect, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2008&year<2013))
summary(rich_NF3_fixed)
AIC(rich_NF3_fixed)
rich_NF3_fixed_aov <- anova(rich_NF3_fixed)
rich_NF3_fixed_aov #trt, year, transect are significant
emmeans(rich_NF3_fixed, ~transect) #THM2 stands out


#native richness 2009-2012 (transect as random)
rich_NF3_random<-lme(richness~trt*year, random=~1|transect, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2008&year<2013))
AIC(rich_NF3_random)
rich_NF3_random_aov<-anova(rich_NF3_random)
rich_NF3_random_aov #year is significant
r.squaredGLMM(rich_NF3_random)

#native richness 2009-2012 (quadrat as random)
rich_NF3_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2008&year<2013))
AIC(rich_NF3_random)
rich_NF3_random_aov<-anova(rich_NF3_random)
rich_NF3_random_aov 
r.squaredGLMM(rich_NF3_random)

########################
#native cover 2005-2008#
########################
cov_NF2<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2004&year<2009))
summary(cov_NF2)
AIC(cov_NF2)
cov_NF2_aov<-anova(cov_NF2)
cov_NF2_aov #trt, year, trt*year

#native cover 2005-2008 (transect as fixed)
cov_NF2_fixed<-lm(relcov~trt*year+transect, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2004&year<2009))
summary(cov_NF2_fixed)
AIC(cov_NF2_fixed)
cov_NF2_fixed_aov<-anova(cov_NF2_fixed)
cov_NF2_fixed_aov #all terms significant
emmeans(cov_NF2_fixed, ~transect) #THM2 on the high end of native cover

#native cover 2005-2008 (transect as random)
cov_NF2_random<-lme(relcov~trt*year, random=~1|transect, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2004&year<2009))
cov_NF2_random_aov<-anova(cov_NF2_random)
AIC(cov_NF2_random)
cov_NF2_random_aov #all terms significant
r.squaredGLMM(cov_NF2_random)#additional 25% of variation explained

#nonnative cover 2005-2008 (quadrat as random)
cov_NF2_random<-lme(relcov~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2004&year<2009))
cov_NF2_random_aov<-anova(cov_NF2_random)
AIC(cov_NF2_random)
cov_NF2_random_aov #all terms significant
r.squaredGLMM(cov_NF2_random) #additional 30% of variation explained by quad

########################
#native cover 2009-2012#
########################
cov_NF3<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2008&year<2013))
summary(cov_NF3)
AIC(cov_NF3)
cov_NF3_aov<-anova(cov_NF3)
cov_NF3_aov

#native cover 2009-2012 (transect as fixed)
cov_NF3_fixed<-lm(relcov~trt*year+transect, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2008&year<2013))
summary(cov_NF3_fixed)
AIC(cov_NF3_fixed)
cov_NF3_fixed_aov<-anova(cov_NF3_fixed)
cov_NF3_fixed_aov #trt, transect
emmeans(cov_NF3_fixed, ~transect)

#native cover 2009-2012 (transect as random)
cov_NF3_random<-lme(relcov~trt*year, random=~1|transect, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2008&year<2013))
cov_NF3_random_aov<-anova(cov_NF3_random)
AIC(cov_NF3_random)
cov_NF3_random_aov #nothing significant
r.squaredGLMM(cov_NF3_random) #improves by ~36%

#nonnative cover 2009-2012 (quadrat as random)
cov_NF3_random<-lme(relcov~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2008&year<2013))
cov_NF3_random_aov<-anova(cov_NF3_random)
AIC(cov_NF3_random)
cov_NF3_random_aov #trt, trt*year
r.squaredGLMM(cov_NF3_random) #45% improvement

###########################
#nonnative cover 2005-2008#
###########################
cov_IG2<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2004&year<2009))
summary(cov_IG2)
AIC(cov_IG2)
cov_IG2_aov<-anova(cov_IG2)
cov_IG2_aov 

#nonnative cover 2005-2008 (transect as fixed)
cov_IG2_fixed<-lm(relcov~trt*year+transect, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2004&year<2009))
summary(cov_IG2_fixed)
AIC(cov_IG2_fixed)
cov_IG2_fixed_aov<-anova(cov_IG2_fixed)
cov_IG2_fixed_aov #all terms significant
emmeans(cov_IG2_fixed, ~transect) #THM1 and THM2 on lower end of nonnative cover

#native cover 2005-2008 (transect as random)
cov_IG2_random<-lme(relcov~trt*year, random=~1|transect, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2004&year<2009))
cov_IG2_random_aov<-anova(cov_IG2_random)
AIC(cov_IG2_random)
cov_IG2_random_aov #year, trt*year
r.squaredGLMM(cov_IG2_random)#~25% improvement

#native cover 2005-2008 (quadrat as random)
cov_IG2_random<-lme(relcov~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2004&year<2009))
cov_IG2_random_aov<-anova(cov_IG2_random)
AIC(cov_IG2_random)
cov_IG2_random_aov #all terms significant
r.squaredGLMM(cov_IG2_random)#~25% improvement

###########################
#nonnative cover 2009-2012#
###########################
cov_IG3<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2008&year<2013))
summary(cov_IG3)
AIC(cov_IG3)
cov_IG3_aov<-anova(cov_IG3)
cov_IG3_aov #trt significant

#nonnative cover 2009-2012 (transect as fixed)
cov_IG3_fixed<-lm(relcov~trt*year+transect, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2008&year<2013))
summary(cov_IG3_fixed)
AIC(cov_IG3_fixed)
cov_IG3_fixed_aov<-anova(cov_IG3_fixed)
cov_IG3_fixed_aov #trt, trt*year, transect
emmeans(cov_IG3_fixed, ~transect) #THM2 on lower end of nonnative cover

#native cover 2009-2012 (transect as random)
cov_IG3_random<-lme(relcov~trt*year, random=~1|transect, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2008&year<2013))
cov_IG3_random_aov<-anova(cov_IG3_random)
AIC(cov_IG3_random)
cov_IG3_random_aov #trt*year
r.squaredGLMM(cov_IG3_random)

#native cover 2009-2012 (quadrat as random)
cov_IG3_random<-lme(relcov~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2008&year<2013))
cov_IG3_random_aov<-anova(cov_IG3_random)
AIC(cov_IG3_random)
cov_IG3_random_aov #trt, trt*year
r.squaredGLMM(cov_IG3_random)

###########################
# litter  years  2006-2008#
###########################
lit2<-lm(lit~trt*year, na.action=na.omit, data = subset(litlit3,year>2005&year<2009))
summary(lit2)
AIC(lit2)
lit2_aov<-anova(lit2)
lit2_aov #all terms

# litter years 2006-2008 (transect as fixed)
lit2_fixed<-lm(lit~trt*year+transect2, na.action=na.omit, data = subset(litlit3,year>2005&year<2009))
summary(lit2_fixed)
AIC(lit2_fixed)
lit_fixed_ao2<-anova(lit2_fixed)
lit_fixed_ao2 #all terms

# litter years 2006-2008 (transect as random)
lit2_random<-lme(lit~trt*year, random=~1|transect2, na.action=na.omit, data = subset(litlit3,year>2005&year<2009))
AIC(lit2_random)
r.squaredGLMM(lit2_random)
lit2_random_ao2<-anova(lit2_random)
lit2_random_ao2 #all terms

#litter 2005-2008 (quadrat as random)
lit2_random<-lme(lit~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(litlit3,year>2005&year<2009))
AIC(lit2_random)
r.squaredGLMM(lit2_random)
lit2_random_ao2<-anova(lit2_random)
lit2_random_ao2 #all terms

###########################
# litter  years  2009-2012#
###########################
lit3<-lm(lit~trt*year, na.action=na.omit, data = subset(litlit3,year>2008&year<2013))
summary(lit3)
AIC(lit3)
lit3_aov<-anova(lit3)
lit3_aov

# litter years 2009-2012 (transect as fixed)
lit3_fixed<-lm(lit~trt*year+transect2, na.action=na.omit, data = subset(litlit3,year>2008&year<2013))
summary(lit3_fixed)
AIC(lit3_fixed)
lit_ao3_fixed<-anova(lit3_fixed)
lit_ao3_fixed

# litter years 2009-2012 (transect as random)
lit3_random<-lme(lit~trt*year, random=~1|transect2, na.action=na.omit, data = subset(litlit3,year>2008&year<2013))
AIC(lit3_random)
r.squaredGLMM(lit3_random)
lit3_random_aov<<-anova(lit3_random)
lit3_random_aov

#litter 2009-2012 (quadrat as random)
lit3_random<-lme(lit~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(litlit3,year>2008&year<2013))
AIC(lit3_random)
r.squaredGLMM(lit3_random)
lit3_random_aov<<-anova(lit3_random)
lit3_random_aov
