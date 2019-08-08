## Set your datpath!! (file in Data_cleaning)

library(tidyverse)
library(readr)
library(ggplot2)

##FN for Calculating SE
calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#load data
dat <- read_csv(paste(datpath_clean, "/bugrdat.csv", sep=""))
SC <- read_csv(paste(datpath_clean, "/SpeciesCodes.csv", sep=""))
clim <- read_csv(paste(datpath_clean, "/san_jose_clim.csv", sep =""), skip = 9) %>%
  mutate(stand_ppt = scale(PRCP, center = T, scale = T)) %>% #standardize to z-scores
  mutate(stand_temp = scale(TAVG, center = T, scale = T)) #standardize to z-scores

View(clim)

#plot timeseries of richness 
dat2<-left_join(dat, SC)
richness0 <- dat2 %>%
  mutate(func=paste(func, status))%>%
  filter(cover != 0, spname != c("Unknown", "Moss")) %>%
  group_by(year, quadratNew, func, thermal)%>%
  summarize(richness = length(unique(spname)))
richness<-richness0%>%
  group_by(year,func, thermal) %>%
  filter(thermal!="?")%>%
  summarize(mean_rich = mean(richness), se_rich=calcSE(richness))

ggplot(richness, aes(year, mean_rich)) +
  geom_line(aes(color=as.factor(func)))+facet_wrap(~thermal) #+ geom_errorbar(aes((color=as.factor(func)), ymin=mean_rich-se_rich, ymax=mean_rich+se_rich), width=.2) +
  geom_line(aes(color=as.factor(func))) + ylab("mean species richness +/-se")

#plot time series of shannon diversity
shandiv<- community_diversity(dat1, time.var = "year", abundance.var="count", replicate.var="site.trt", metric = c("Shannon")))%>%
  group_by(year, site, alltrt, Shannon)%>%
  summarize(richness=length(unique(as.factor(code))))%>%
  group_by(year, alltrt)%>%
  summarize(meanrich=mean(richness), meanshan=mean(Shannon))

ggplot(richsum)+geom_line(aes(as.factor(year), meanshan, group=alltrt, color=alltrt))

#join data and species key
tog<-left_join(dat, SC)
View(tog)

#plot timeseries of cover by thermal
functog<-tog%>%
  group_by(quadratNew, year, status, func, thermal)%>%
  summarize(sumcov=sum(cover))%>%
  filter(!is.na(status), !is.na(func))
functogagg<-functog%>%
  group_by(year, status, func)%>%
  summarize(meancov=mean(sumcov), se_cov=calcSE(sumcov))

View(functoagg)
View(functog)

thermal_trend <- functog %>%
  group_by(year, status, func, thermal) %>%
  summarize(meancov = mean(sumcov), se_cov=calcSE(sumcov)) %>%
  filter(thermal != "?")

ggplot(functog, aes(as.factor(year), sumcov))+geom_boxplot()+facet_grid(status~func)
ggplot(functogagg, aes((year), meancov))+geom_line(aes(color=interaction(status, func)))+geom_point(aes(color=interaction(status, func)))+
  geom_errorbar(aes(ymin=meancov-se_cov, ymax=meancov+se_cov, color=interaction(status, func)), width=.2)

#cover by thermal transposed on annual precip
ggplot(thermal_trend, aes((year), meancov)) + facet_grid(status~func) +
  geom_bar(data = clim, aes(x = DATE, y = PRCP*3), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color= thermal))+ geom_point(aes(color = thermal))  +
  geom_errorbar(aes(ymin=meancov-se_cov, ymax=meancov+se_cov, color=thermal), width=.2)+
  scale_y_continuous(sec.axis = sec_axis(~./3, name = "Annual Precipitation in inches"))

#timeseries of ppt
ggplot(clim, aes(DATE,PRCP)) + geom_line() + geom_point()

#cover by thermal transpoesd on annual mean temp
ggplot(thermal_trend, aes((year), meancov)) + facet_grid(status~func) +
  geom_bar(data = clim, aes(x = DATE, y = stand_temp*30), stat = "identity", fill = "lightgrey") +
  geom_line(aes(color= thermal))+ geom_point(aes(color = thermal))  +
  geom_errorbar(aes(ymin=meancov-se_cov, ymax=meancov+se_cov, color=thermal), width=.2)+
  scale_y_continuous(sec.axis = sec_axis(~./30, name = "Annual Mean Temp Deviation (z-scores)"))

#join clim and thermal_trend
clim1 <- clim %>%
  rename(year = DATE)
joined_dat <- left_join(thermal_trend, clim1)

#cover by temp
ggplot(joined_dat, aes((PRCP), meancov)) + facet_grid(status~func) +
   geom_point(aes(color = thermal)) + geom_smooth(aes(color = thermal), method = "lm", se = F)

dat1<-dat%>%
  select(-1)%>%
  group_by(year, spname, quadratNew, thermal)%>%
  summarize(cover=sum(cover))%>%
  ungroup()%>%
  mutate(sitetrt=as.factor(paste(quadratNew, thermal, sep="_")))%>%
  filter(thermal!="?")

#year by year species turnover
turnover1<-turnover(dat1,
                    time.var="year",
                    species.var="spname",
                    abundance.var="cover",
                    replicate.var="sitetrt")%>%
  separate(sitetrt, into=c("site", "trt"), sep="_")%>%
  group_by( trt, year)%>%
  summarize(mean.turnover=mean(as.numeric(total)), se.turnover=calcSE(as.numeric(total)))


ggplot(turnover1)+
  geom_point(aes((year), mean.turnover, color=trt))+
  geom_line(aes((year), mean.turnover, color=trt))+
  labs(x="Year", y="Annual Turnover")+
  geom_errorbar(aes(color=trt, (year), ymin=mean.turnover-se.turnover, ymax=mean.turnover+se.turnover, ), width=.2)

# year by year species rank abundance shift
rankshift <- rank_shift(dat1,
                        time.var="year",
                        species.var="spname",
                        abundance.var="cover",
                        replicate.var="sitetrt")%>%
  separate(sitetrt, into=c("site", "trt"), sep="_")%>%
  group_by( trt, year_pair)%>%
  summarize(mean.MRS=mean(MRS), se.MRS=calcSE(MRS))%>%
  mutate(year=substr(year_pair, 6,9))

ggplot(rankshift)+
  geom_line(aes(as.factor(year), mean.MRS, group=trt, color=trt))+
  geom_errorbar(aes(color=trt, as.factor(year), ymin=mean.MRS-se.MRS, ymax=mean.MRS+se.MRS, ), width=.1)+
  labs(x="Year", y="Mean Rank Shift")