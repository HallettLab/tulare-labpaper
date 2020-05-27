library(tidyverse)
library(readr)
library(sjstats) #standardized effect size cohen's f
library(nlme) #linear mixed models
library(multcomp) #tukey

#load updated master data, "alldat", in grazing_recovery.R

########
#one-way ANOVA analyze within each year 
#richness, cover, evenness, litter as dependent variables
#focus on two functional groups: native forbs and non-native grasses
########

#Richness
#load "rich" from grazing_recovery.R

TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2005, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2006, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2007, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2008, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2009, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2010, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2011, func == "forb native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2012, func == "forb native")))

TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2005, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2006, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2007, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2008, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2009, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2010, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2011, func == "grass non-native")))
TukeyHSD(aov(richness~trt, data = rich%>%filter(year == 2012, func == "grass non-native")))

rich2005 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2005)
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2005%>%filter(func == "grass native"))))

rich2006 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2006)
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2006%>%filter(func == "grass native"))))

rich2007 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2007)
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2007%>%filter(func == "grass native"))))

rich2008 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year == 2008)
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich2008%>%filter(func == "grass native"))))

rich2009 <- rich %>%
  filter(year == 2009)
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2009%>%filter(func == "grass native"))))

rich2010 <- rich %>%
  filter(year == 2010)
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2010%>%filter(func == "grass native"))))

rich2011 <- rich %>%
  filter(year == 2011)
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2011%>%filter(func == "grass native"))))

rich2012 <- rich %>%
  filter(year == 2012)
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich2012%>%filter(func == "grass native"))))

#cover
#load "cov" from grazing_recovery.R

TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2005, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2006, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2007, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2008, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2009, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2010, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2011, func == "forb native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2012, func == "forb native")))

TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2005, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2006, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2007, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2008, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2009, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2010, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2011, func == "grass non-native")))
TukeyHSD(aov(relcov~trt, data = cov%>%filter(year == 2012, func == "grass non-native")))
cov2005 <- cov %>%
  filter(year == 2005) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(sumcov~graze, data = cov2005%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2005%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2005%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2005%>%filter(func == "grass native"))))

cov2006 <- cov %>%
  filter(year == 2006) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(sumcov~graze, data = cov2006%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2006%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2006%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2006%>%filter(func == "grass native"))))

cov2007 <- cov %>%
  filter(year == 2007) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(sumcov~graze, data = cov2007%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2007%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2007%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2007%>%filter(func == "grass native"))))

cov2008 <- cov %>%
  filter(year == 2008) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(sumcov~graze, data = cov2008%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2008%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2008%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov2008%>%filter(func == "grass native"))))

cov2009 <- cov %>%
  filter(year == 2009)
anova_stats(anova(lm(sumcov~trt, data = cov2009%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2009%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2009%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2009%>%filter(func == "grass native"))))

cov2010 <- cov %>%
  filter(year == 2010)
anova_stats(anova(lm(sumcov~trt, data = cov2010%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2010%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2010%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2010%>%filter(func == "grass native"))))

cov2011 <- cov %>%
  filter(year == 2011)
anova_stats(anova(lm(sumcov~trt, data = cov2011%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2011%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2011%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2011%>%filter(func == "grass native"))))

cov2012 <- cov %>%
  filter(year == 2012)
anova_stats(anova(lm(sumcov~trt, data = cov2012%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2012%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2012%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov2012%>%filter(func == "grass native"))))

#evenness
#load object "shan2" from grazing_recovery.R
shan3 <- shan2 %>%
  mutate(trt1 = trt) %>%
  separate(trt, into=c("grazed", "burn"), sep=" ")

TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2005, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2006, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2007, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2008, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2009, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2010, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2011, func == "forb native")))
TukeyHSD(aov(Shannon~trt, data = shan2%>%filter(year == 2012, func == "forb native")))

shan2005 <- shan3%>%
  filter(year == 2005) %>%
  filter(trt1 != "ungrazed unburned")
anova_stats(anova(lm(Shannon~grazed, data = shan2005%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2005%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2005%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2005%>%filter(func == "grass native"))))

shan2006 <- shan3%>%
  filter(year == 2006) %>%
  filter(trt1 != "ungrazed unburned")
anova_stats(anova(lm(Shannon~grazed, data = shan2006%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2006%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2006%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2006%>%filter(func == "grass native"))))

shan2007 <- shan3%>%
  filter(year == 2007) %>%
  filter(trt1 != "ungrazed unburned")
anova_stats(anova(lm(Shannon~grazed, data = shan2007%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2007%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2007%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2007%>%filter(func == "grass native"))))

shan2008 <- shan3%>%
  filter(year == 2008) %>%
  filter(trt1 != "ungrazed unburned")
anova_stats(anova(lm(Shannon~grazed, data = shan2008%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2008%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2008%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan2008%>%filter(func == "grass native"))))

shan2009 <- shan3%>%
  filter(year == 2009)
anova_stats(anova(lm(Shannon~trt1, data = shan2009%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2009%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2009%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2009%>%filter(func == "grass native"))))

shan2010 <- shan3%>%
  filter(year == 2010)
anova_stats(anova(lm(Shannon~trt1, data = shan2010%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2010%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2010%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2010%>%filter(func == "grass native"))))

shan2011 <- shan3%>%
  filter(year == 2011)
anova_stats(anova(lm(Shannon~trt1, data = shan2011%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2011%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2011%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2011%>%filter(func == "grass native"))))

shan2012 <- shan3%>%
  filter(year == 2012)
anova_stats(anova(lm(Shannon~trt1, data = shan2012%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2012%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2012%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan2012%>%filter(func == "grass native"))))

#litter
#pull in object "joindat" from grazingrecovery
lit <- joindat %>%
  mutate(trt=paste(graze, burn)) %>%
  ungroup() %>%
  dplyr::select(quadratNew, year, trt, transect, burn, graze, litter)

TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2005)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2006)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2007)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2008)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2009)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2010)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2011)))
TukeyHSD(aov(litter~trt, data = lit%>%filter(year == 2012)))

lit2005 <- lit%>%
  filter(year == 2005) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(litter~graze, data = lit2005)))

lit2006 <- lit%>%
  filter(year == 2006) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(litter~graze, data = lit2006)))

lit2007 <- lit%>%
  filter(year == 2007) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(litter~graze, data = lit2007)))

lit2008 <- lit%>%
  filter(year == 2008) %>%
  filter(trt != "ungrazed unburned")
anova_stats(anova(lm(litter~graze, data = lit2008)))

lit2009 <- lit%>%
  filter(year == 2009) 
anova_stats(anova(lm(litter~trt, data = lit2009)))

lit2010 <- lit%>%
  filter(year == 2010) 
anova_stats(anova(lm(litter~trt, data = lit2010)))

lit2011 <- lit%>%
  filter(year == 2011) 
anova_stats(anova(lm(litter~trt, data = lit2011)))

lit2012 <- lit%>%
  filter(year == 2012) 
anova_stats(anova(lm(litter~trt, data = lit2012)))

########
#one-way ANOVA analyze two blocks of years: 2005-2008 and 2008-2012 
########
rich3 <- rich %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year%in%c(2005, 2006, 2007, 2008))

anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~graze, data = rich3%>%filter(func == "grass native"))))

rich4 <- rich %>%
  filter(year%in%c(2008:2012))

anova_stats(anova(lm(richness~trt, data = rich4%>%filter(func == "forb native"))))
anova_stats(anova(lm(richness~trt, data = rich4%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(richness~trt, data = rich4%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(richness~trt, data = rich4%>%filter(func == "grass native"))))

cov_pre <- cov %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year%in%c(2005:2008))

anova_stats(anova(lm(sumcov~graze, data = cov_pre%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~graze, data = cov_pre%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov_pre%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~graze, data = cov_pre%>%filter(func == "grass native"))))

cov_post <- cov %>%
  filter(year%in%c(2008:2012))
anova_stats(anova(lm(sumcov~trt, data = cov_post%>%filter(func == "forb native"))))
anova_stats(anova(lm(sumcov~trt, data = cov_post%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov_post%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(sumcov~trt, data = cov_post%>%filter(func == "grass native"))))

shan_pre <- shan3 %>%
  filter(trt1 != "ungrazed unburned") %>%
  filter(year%in%c(2005:2008))
anova_stats(anova(lm(Shannon~grazed, data = shan_pre%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan_pre%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan_pre%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~grazed, data = shan_pre%>%filter(func == "grass native"))))

shan_post <- shan3 %>%
  filter(year%in%c(2008:2012))
anova_stats(anova(lm(Shannon~trt1, data = shan_post%>%filter(func == "forb native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan_post%>%filter(func == "grass non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan_post%>%filter(func == "forb non-native"))))
anova_stats(anova(lm(Shannon~trt1, data = shan_post%>%filter(func == "grass native"))))

lit_pre <- lit %>%
  filter(trt != "ungrazed unburned") %>%
  filter(year%in%c(2005:2008))
anova_stats(anova(lm(litter~graze, data = lit_pre)))

lit_post <- lit %>%
  filter(year%in%c(2008:2012))
anova_stats(anova(lm(litter~trt, data = lit_post)))






########
## LINEAR MIXED MODELS
#grazed and burned as fixed effects, quadrats nested within transects as random
########

# richness by year
richrich<-rich%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 
#set the intercept to ungrazed burned for the third pair comparison
richrich2<-richrich%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))%>%
  filter(!is.na(richness))
richrich2$trt <- factor(richrich2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))
#load "cov" from grazing_recovery.R
covcov<-cov%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 
#set the intercept to ungrazed burned for the third pair comparison
covcov2<-covcov%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))
covcov2$trt <- factor(covcov2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))
#load "lit"
litlit<-lit%>%
  separate(quadratNew, into=c("transect2", "quadrat"), sep="-") 
#set the intercept to ungrazed burned for the third pair comparison
litlit2<-litlit%>%
  ungroup()%>%
  mutate(trt=as.factor(trt))%>%
  filter(!is.na(litter))
litlit2$trt <- factor(litlit2$trt, levels = c("ungrazed burned", "ungrazed unburned", "grazed burned"))

library(MuMIn)
library(lme4)
##### MM:RICHNESS #####
#native 2005-2008
rich_NF2005.2008<-lme(richness~trt+year, random = ~1|transect2/quadrat, data = subset(richrich2, func=="forb native"&year<2009))
summary(rich_NF2005.2008)
anova(rich_NF2005.2008)
summary(glht(rich_NF2005.2008, linfct=mcp(trt="Tukey")))
r.squaredGLMM(rich_NF2005.2008)

#native 2009-2012
rich_NF2009.2012<-lme(richness~trt+year, random = ~1|transect2/quadrat, data = subset(richrich2, func=="forb native"&year>2008))
summary(rich_NF2009.2012)
anova(rich_NF2009.2012)
summary(glht(rich_NF2009.2012, linfct=mcp(trt="Tukey")))
r.squaredGLMM(rich_NF2009.2012)


#nonnative 2005-2008
rich_IG2005.2008<-lme(richness~trt+year, random = ~1|transect2/quadrat, data = subset(richrich2, func=="grass non-native"&year<2009))
summary(rich_IG2005.2008)
summary(glht(rich_IG2005.2008, linfct=mcp(trt="Tukey")))

#nonnative 2009-2012
rich_IG2009.2012<-lme(richness~trt+year, random = ~1|transect2/quadrat, data = subset(richrich2, func=="grass non-native"&year>2008))
summary(rich_IG2009.2012)
summary(glht(rich_IG2009.2012, linfct=mcp(trt="Tukey")))



##### MM:COVER #####

#native 2005-2008
cov_NF0508<-lme(relcov~trt+year, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year<2009))
summary(cov_NF0508)
summary(glht(cov_NF0508, linfct=mcp(trt="Tukey")))

#native 2009-2012
cov_NF0912<-lme(relcov~trt+year, random = ~1|transect2/quadrat, data = subset(covcov2, func=="forb native"&year>2008))
summary(cov_NF0912)
summary(glht(cov_NF0912, linfct=mcp(trt="Tukey")))

#nonnative 2005-2008
cov_IG0508<-lme(relcov~trt+year, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year<2009))
summary(cov_IG0508)
summary(glht(cov_IG0508, linfct=mcp(trt="Tukey")))

#non-native 2009-2012
cov_IG0912<-lme(relcov~trt+year, random = ~1|transect2/quadrat, data = subset(covcov2, func=="grass non-native"&year>2008))
summary(cov_IG0912)
summary(glht(cov_IG0912, linfct=mcp(trt="Tukey")))

##### MM:LITTER #####

#lme of litter each year from 2006-2012
lit_0508<-lme(litter~trt+year, random = ~1|transect2/quadrat, data = subset(litlit2,year<2009))
summary(lit_0508)
summary(glht(lit_0508, linfct=mcp(trt="Tukey")))

lit_0912<-lme(litter~trt+year, random = ~1|transect2/quadrat, data = subset(litlit2, year>2008))
summary(lit_0912)
summary(glht(lit_0912, linfct=mcp(trt="Tukey")))






##plot lme stats

lmestats<-read_csv(paste(datpath_clean, "/Grazing_recovery_stats_lme.csv", sep=""))%>%
  select(-Test)%>%
  rename("p"="p-value", "p1"="p-value_1", "t"="t-value", "t1"="t-value_1", "ugub"="ungrazed unburned", "ugb"="ungrazed burned")

lmestats<-lmestats%>%
  mutate(p=ifelse(p<.05, .05, ifelse(p<.1, .1, NA)))%>%
  mutate(p1=ifelse(p1<.05, .05, ifelse(p1<.1, .1, NA)))
  
lme1<-ggplot(subset(lmestats, func=="native forb"&response=="cover"), aes(year, ugb*100)) +geom_line(color="blue") +geom_point(color="blue") +
  geom_line(aes(year, ugub*100), color="red") +geom_point(aes(year, ugub*100), color="red") +
  geom_point(aes(year-.1, y=-p*100),shape=8, color="blue")+geom_point(aes(year+.1,y=-p1*100),shape=8, color="red")+
  geom_hline(yintercept=0) +
  annotate("text", x= 2004.75, y = .5, label = "fire", size = 3) +
  annotate("text", x= 2009, y = .5, label = "grazing", size = 3)+
  geom_vline(xintercept=2008.5, color = "grey66", lty =2)+
  geom_vline(xintercept=2004.5, color="grey66", lty = 2)+
  xlab("")+ylab("% Cover Effect Size vs Burned-Grazed")

lme2<-ggplot(subset(lmestats, func=="non-native grass"&response=="cover"), aes(year, ugb*100)) +geom_line(color="blue") +geom_point(color="blue") +
  geom_line(aes(year, ugub*100), color="red") +geom_point(aes(year, ugub*100), color="red") +
  geom_point(aes(year-.1, y=p*100),shape=8, color="blue")+geom_point(aes(year+.1,y=p1*100),shape=8, color="red")+
  geom_hline(yintercept=0) +
  annotate("text", x= 2004.75, y = .5, label = "fire", size = 3) +
  annotate("text", x= 2009, y = .5, label = "grazing", size = 3)+
  geom_vline(xintercept=2008.5, color = "grey66", lty =2)+
  geom_vline(xintercept=2004.5, color="grey66", lty = 2)+
  xlab("")+ylab("")


lme3<-ggplot(subset(lmestats, func=="native forb"&response=="richness"), aes(year, ugb)) +geom_line(color="blue") +geom_point(color="blue") +
  geom_line(aes(year, ugub), color="red") +geom_point(aes(year, ugub), color="red") +
  geom_point(aes(year-.1, y=p*10),shape=8, color="blue")+geom_point(aes(year+.1,y=p1*10),shape=8, color="red")+
  geom_hline(yintercept=0) +
  annotate("text", x= 2004.75, y = .5, label = "fire", size = 3) +
  annotate("text", x= 2009, y = .5, label = "grazing", size = 3)+
  geom_vline(xintercept=2008.5, color = "grey66", lty =2)+
  geom_vline(xintercept=2004.5, color="grey66", lty = 2)+
  xlab("")+ylab("Richness Effect Size vs Burned-Grazed")


lme4<-ggplot(subset(lmestats, func=="non-native grass"&response=="richness"), aes(year, ugb)) +geom_line(color="blue") +geom_point(color="blue") +
  geom_line(aes(year, ugub), color="red") +geom_point(aes(year, ugub), color="red") +
  geom_point(aes(year-.1, y=p),shape=8, color="blue")+geom_point(aes(year+.1,y=p1),shape=8, color="red")+
  geom_hline(yintercept=0) +
  annotate("text", x= 2004.75, y = .5, label = "fire", size = 3) +
  annotate("text", x= 2009, y = .5, label = "grazing", size = 3)+
  geom_vline(xintercept=2008.5, color = "grey66", lty =2)+
  geom_vline(xintercept=2004.5, color="grey66", lty = 2)+
  xlab("")+ylab("")

library(ggpubr)
ggarrange(lme1, lme2, lme3, lme4,  ncol = 2, nrow = 2, 
          labels = c("a) Native forb cover", "b) Non-native grass cover",
                     "c) Native forb richness", "d) Non-native grass richness"),
          common.legend = TRUE, legend = "right", 
          font.label = list(size = 10),
          hjust = c(-0.5, -0.35, -0.5, -0.35))
  
#########################
#repeated measures
########################
library(lsmeans)
richrich2$year<- factor(richrich2$year, ordered=TRUE)
levels(richrich2$year)

#nonnative richness all years
rich_IG<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2004&year<2013))
summary(rich_IG)
rich_IG_aov<-anova(rich_IG)
rich_IG_aov
LS.IG1<-lsmeans(rich_IG, ~year*trt)
LS.IG1.cont<-contrast(LS.IG1, "pairwise", by="year")
LSig1<-summary(LS.IG1.cont)

#nonnative richness 2005-2008
rich_IG2<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2004&year<2009))
summary(rich_IG2)
rich_IG2_aov<-anova(rich_IG2)
rich_IG2_aov
LS.IG2<-lsmeans(rich_IG2, ~year*trt)
LS.IG2.cont<-contrast(LS.IG2, "pairwise", by="year")
LSig2<-summary(LS.IG2.cont)

#nonnative richness 2009-2012
rich_IG3<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="grass non-native"&year>2008&year<2013))
summary(rich_IG3)
rich_IG3_aov<-anova(rich_IG3)
rich_IG3_aov
LS.IG3<-lsmeans(rich_IG3, ~year*trt)
LS.IG3.cont<-contrast(LS.IG3, "pairwise", by="year")
LSig3<-summary(LS.IG3.cont)

#native richness all years
rich_NF<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2004&year<2013))
summary(rich_NF)
rich_NF_aov<-anova(rich_NF)
rich_NF_aov
LS.NF1<-lsmeans(rich_NF, ~year*trt)
LS.NF1.cont<-contrast(LS.NF1, "pairwise", by="year")
LSnf1<-summary(LS.NF1.cont)

#native richness 2005-2008
rich_NF2<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2004&year<2009))
summary(rich_NF2)
rich_NF2_aov<-anova(rich_NF2)
rich_NF2_aov
LS.NF2<-lsmeans(rich_NF2, ~year*trt)
LS.NF2.cont<-contrast(LS.NF2, "pairwise", by="year")
LSnf2<-summary(LS.NF2.cont)

#native richness 2009-2012
rich_NF3<-lm(richness~trt*year, na.action=na.omit, data = subset(richrich2, func=="forb native"&year>2008&year<2013))
summary(rich_NF3)
rich_NF3_aov<-anova(rich_NF3)
rich_NF3_aov
LS.NF3<-lsmeans(rich_NF3, ~year*trt)
LS.NF3.cont<-contrast(LS.NF3, "pairwise", by="year")
LSnf3<-summary(LS.NF3.cont)


covcov2$year<- factor(covcov2$year, ordered=TRUE)
levels(covcov2$year)

#native cover all years
cov_NF<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2004&year<2013))
summary(cov_NF)
cov_NF_aov<-anova(cov_NF)
cov_NF_aov
LS.cNF1<-lsmeans(cov_NF, ~year*trt)
LS.cNF1.cont<-contrast(LS.cNF1, "pairwise", by="year")
LSnfc1<-summary(LS.cNF1.cont)

#native cover 2005-2008
cov_NF2<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2004&year<2009))
summary(cov_NF2)
cov_NF2_aov<-anova(cov_NF2)
cov_NF2_aov
LS.cNF2<-lsmeans(cov_NF2, ~year*trt)
LS.cNF2.cont<-contrast(LS.cNF2, "pairwise", by="year")
LSnfc2<-summary(LS.cNF2.cont)

#native cover 2009-2012
cov_NF3<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="forb native"&year>2008&year<2013))
summary(cov_NF3)
cov_NF3_aov<-anova(cov_NF3)
cov_NF3_aov
LS.cNF3<-lsmeans(cov_NF3, ~year*trt)
LS.cNF3.cont<-contrast(LS.cNF3, "pairwise", by="year")
LSnfc3<-summary(LS.cNF3.cont)

#nonnative cover all years
cov_IG<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2004&year<2013))
summary(cov_IG)
cov_IG_aov<-anova(cov_IG)
cov_IG_aov
LS.cIG1<-lsmeans(cov_IG, ~year*trt)
LS.cIG1.cont<-contrast(LS.cIG1, "pairwise", by="year")
LSigc1<-summary(LS.cIG1.cont)

#nonnative cover 2005-2008
cov_IG2<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2004&year<2009))
summary(cov_IG2)
cov_IG2_aov<-anova(cov_IG2)
cov_IG2_aov
LS.cIG2<-lsmeans(cov_IG2, ~year*trt)
LS.cIG2.cont<-contrast(LS.cIG2, "pairwise", by="year")
LSigc2<-summary(LS.cIG2.cont)

#nonnative cover 2009-2012
cov_IG3<-lm(relcov~trt*year, na.action=na.omit, data = subset(covcov2, func=="grass non-native"&year>2008&year<2013))
summary(cov_IG3)
cov_IG3_aov<-anova(cov_IG3)
cov_IG3_aov
LS.cIG3<-lsmeans(cov_IG3, ~year*trt)
LS.cIG3.cont<-contrast(LS.cIG3, "pairwise", by="year")
LSigc3<-summary(LS.cIG3.cont)


litlit2$year<- factor(litlit2$year, ordered=TRUE)
levels(litlit2$year)
litlit3<-litlit2 %>% group_by(transect2, quadrat,  year, trt) %>% summarize(lit=mean(litter))


# litter all years from 2006-2012
lit.m<-lm(lit~trt*year, na.action=na.omit, data = subset(litlit3,year>2006&year<2013))
summary(lit.m)
lit_aov<-anova(lit.m)
lit_aov
LS.lit<-lsmeans(lit, ~year*trt)
LS.lit.cont<-contrast(LS.lit, "pairwise", by="year")
LSlit<-summary(LS.lit.cont)

# litter  years  2006-2008
lit2<-lm(lit~trt*year, na.action=na.omit, data = subset(litlit3,year>2006&year<2009))
summary(lit2)
lit_ao2<-anova(lit2)
lit_ao2
LS.lit2<-lsmeans(lit2, ~year*trt)
LS.lit2.cont<-contrast(LS.lit2, "pairwise", by="year")
LSlit2<-summary(LS.lit2.cont)

# litter  years  2009-2012
lit3<-lm(lit~trt*year, na.action=na.omit, data = subset(litlit3,year>2008&year<2013))
summary(lit3)
lit_ao3<-anova(lit3)
lit_ao3
LS.lit3<-lsmeans(lit3, ~year*trt)
LS.lit3.cont<-contrast(LS.lit3, "pairwise", by="year")
LSlit3<-summary(LS.lit3.cont)


                                                     