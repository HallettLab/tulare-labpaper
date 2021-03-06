library(tidyverse)
library(readxl)
library(stringr)

## GOPHER, COWPIE, BARE, LITTER
dat <- read_excel(datpath, 
                  sheet = "Plants2008", skip = 2) 

datnames <- read_excel(datpath, 
                       sheet = "Plants2008")[0,]

names(dat) = names(datnames)
names(dat)[1:4] = c("year", "quadrat", "transect", "site")


dat2008env <- dat %>%
  select(year, quadrat, transect, site, GOPHER, BARE, ROCK, LITTER, COWPIE)
 

dat2008 <- dat %>%
  select(-c(GOPHER, BARE, ROCK, LITTER, COWPIE)) %>%
  gather(species, cover, "X__5":"X__28") %>%
  mutate(dummyspp = substr(species, 1,3)) %>%
  filter(dummyspp != "X__") %>%
  select(-dummyspp) %>%
  #mutate(site = substr(quadrat, 1,2)) %>%
  filter(site == "TH" | site =="TH-BUG") %>%
  mutate(year = 2008)


key2008 <- dat2008 %>%
  select(transect, quadrat, site, year) %>%
  unique()

rm(dat, datnames)
