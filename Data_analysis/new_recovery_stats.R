#load updated master data, "alldat", in grazing_recovery.R
#load "rich" from grazing_recovery.R
# richness by year
richrich<-rich%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 
#set the intercept to ungrazed burned for the third pair comparison
richrich2<-richrich%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))%>%
  filter(!is.na(richness))
richrich2$trt <- factor(richrich2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))
richrich2$quadUnique <- paste(richrich2$transect, richrich2$quadrat, sep = "_")

#load "cov" from grazing_recovery.R
covcov<-cov%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 
#set the intercept to ungrazed burned for the third pair comparison
covcov2<-covcov%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))
covcov2$trt <- factor(covcov2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))
covcov2$quadUnique <- paste(covcov2$transect, covcov2$quadrat, sep = "_")

#litter
#pull in object "joindat" from grazingrecovery
lit <- joindat %>%
  mutate(trt=paste(graze, burn)) %>%
  ungroup() %>%
  dplyr::select(quadratNew, year, trt, transect, burn, graze, litter)
#load "lit"
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
litlit3$quadUnique <- paste(litlit3$transect, litlit3$quadrat, sep = "_")



#nonnative richness 2005-2008 (quadrat as random)
rich_IG2_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2004&year<2009))
AIC(rich_IG2_random)
r.squaredGLMM(rich_IG2_random)

#nonnative richness 2009-2012
rich_IG3_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2008&year<2013))
AIC(rich_IG3_random)
r.squaredGLMM(rich_IG3_random)

#native richness 2005-2008
rich_NF2_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2004&year<2009))
AIC(rich_NF2_random)
r.squaredGLMM(rich_NF2_random)

#native richness 2009-2012
rich_NF3_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2008&year<2013))
AIC(rich_NF3_random)
r.squaredGLMM(rich_NF3_random)

#nonnative cover 2005-2008
cov_NF2_random<-lme(relcov~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2004&year<2009))
AIC(cov_NF2_random)
r.squaredGLMM(cov_NF2_random)

#nonnative cover 2009-2012
cov_NF3_random<-lme(relcov~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2008&year<2013))
AIC(cov_NF3_random)
r.squaredGLMM(cov_NF3_random)

#native cover 2005-2008
cov_IG2_random<-lme(relcov~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2004&year<2009))
AIC(cov_IG2_random)
r.squaredGLMM(cov_IG2_random)

#native cover 2009-2012
cov_IG3_random<-lme(relcov~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2008&year<2013))
AIC(cov_IG3_random)
r.squaredGLMM(cov_IG3_random)

#litter 2005-2008
lit2_random<-lme(lit~trt*year, random=(~1|quadUnique), na.action=na.omit, data = subset(litlit3,year>2005&year<2009))
AIC(lit2_random)
r.squaredGLMM(lit2_random)

#litter 2009-2012
lit3_random<-lme(lit~trt*year, random=(~1|quadUnique), na.action=na.omit, data = subset(litlit3,year>2008&year<2013))
AIC(lit3_random)
r.squaredGLMM(lit3_random)
