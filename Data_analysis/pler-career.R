library(tidyverse)

## Set ggplot2 theme
theme_set(theme_classic())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 18),
              strip.text= element_text(size = 18), 
              axis.text = element_text(size = 18))


calcSE<-function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}


alldat<-read_csv(paste(datpath_clean, "/alldatsptrt.csv", sep="")) %>%
  select(-1)%>%
  filter(transect%in%c("THBUGM1", "THBUGM2", "THM1", "THM2", "THM3", "THM4", "THUBUGM1", "THUBUGM2"))%>%
  filter(!quadratNew%in%c("THM1-1", "THM3-3", "THM1-10"))%>%
  filter(thermal=="moderate")%>%
  group_by(year, spname, spcode, quadratNew, status, type, transect, burn, graze)%>%
  summarize(cover=sum(cover))

pler <- alldat %>%
  filter(spcode == "PLA ERE") %>%
  filter(spcode == "PLA ERE") %>%
  group_by(burn, year, graze, transect)  %>%
  mutate(trt=paste(graze, burn))%>%
  group_by(year, trt)  %>%
  summarize(meancover = mean(cover), secover = calcSE(cover)) %>%
  filter(trt != "ungrazed burned") 


ggplot(subset(pler, year > 2006), aes(x = year, y = meancover, color = trt)) + geom_line() + geom_point(size = 3) +
  geom_errorbar(aes(ymin=meancover-secover, ymax=meancover+secover ), width = .2) + theme(legend.position = "none") + 
  labs(x="Year", y = "Percent cover") + 
  annotate(geom = 'text', x= 2008.5, y = 5, 
           label = "atop(Cattle, reintroduced)", 
           parse = TRUE, size=5) +
  geom_segment(aes(x=2008.5, y=2.5, xend=2008.5, yend=1), arrow=arrow(length = unit(0.03, "npc")),  color = "grey0") +
  scale_color_manual(values = c("#625377", "#AFA1C4"))
ggsave("pler-career.pdf", width = 5, height = 4)
#855C75,#D9AF6B,#AF6458,#736F4C,#526A83,#625377,#68855C,#9C9C5E,#A06177,#8C785D,#467378,#7C7C7C


ggplot(pler, aes(x=year, y = meancover, color = interaction(burn, graze))) + geom_line() + 
  ggplot(pler, aes(x=year, y = meancover, color = trt)) + geom_line() + 
  geom_point() + facet_grid(burn~graze)

geom_point(aes(fill = trt), pch = 21) +
  geom_errorbar(aes(ymin=meancover-secover, ymax=meancover+secover, color=as.factor(trt)), width=.2) +
  labs(x = "Year", y = "Absolute Cover (%)") +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(2000, 2018), breaks = c(2000, 2005, 2010, 2015)) +
  scale_fill_manual(values = c("grey0", "grey85", "grey100"))+
  scale_color_manual(values= c("grey0", "grey36", "grey65"), guide = FALSE)  +
  annotate("text", x= 2004.5, y = 5, label = "fire", size = 3) +
  geom_segment(aes(x=2004.5, y=2.5, xend=2004.5, yend=.5), arrow=arrow(length = unit(0.03, "npc")), color = "grey0") + 
  annotate(geom = 'text', x= 2008.5, y = 30, 
           label = "atop(cattle, reintroduced)", 
           parse = TRUE, size=3) +
  geom_segment(aes(x=2008.5, y=25, xend=2008.5, yend=20), arrow=arrow(length = unit(0.03, "npc")),  color = "grey0")