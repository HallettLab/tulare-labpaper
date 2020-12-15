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
litlit2$quadUnique <- paste(litlit2$transect, litlit2$quadrat, sep = "_")



#nonnative richness 2005-2008 (quadrat as random)
rich_IG2_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2004&year<2009))
AIC(rich_IG2_random)
r.squaredGLMM(rich_IG2_random)

#nonnative richness 2009-2012
rich_IG3_random<-lme(richness~trt*year, random=~1|quadUnique, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2008&year<2013))
AIC(rich_IG3_random)
r.squaredGLMM(rich_IG3_random)
#native richness 2005-2008
#native richness 2009-2012

#nonnative cover 2005-2008
#nonnative cover 2009-2012
#native cover 2005-2008
#native cover 2009-2012

#litter 2005-2008
#litter 2009-2012
