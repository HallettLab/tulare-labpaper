## to make manuscript plots with ellipses, run BUGR_nmds_msFigs.R script first

# function for ellipsess 
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

e2005<-subset(spscoresall.mod.e, year==2005)
e2005$treatment<-as.factor(e2005$treatment)
levels(e2005$treatment)
#data for ellipse, in this case using the management factor
df_ell.fire <- data.frame() #sets up a data frame before running the function.
for(g in levels(e2005$treatment)){
  df_ell.fire  <- rbind(df_ell.fire , cbind(as.data.frame(with(e2005 [e2005$treatment==g,],
      veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2005=aggregate(e2005[ ,c("NMDS1", "NMDS2")], 
                         list(group = e2005$treatment), mean)

fig1b.e<-fig1b+ 
  geom_path(data = df_ell.fire, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse, separate ones by treatment 


e2006<-subset(spscoresall.mod.e, year==2006)
e2006$treatment<-as.factor(e2006$treatment)
levels(e2006$treatment)
#data for ellipse, in this case using the management factor
df_ell.fire.06 <- data.frame() #sets up a data frame before running the function.
for(g in levels(e2006$treatment)){
  df_ell.fire.06  <- rbind(df_ell.fire.06 , cbind(as.data.frame(with(e2006 [e2006$treatment==g,],
                                  veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2006=aggregate(e2006[ ,c("NMDS1", "NMDS2")], 
                         list(group = e2006$treatment), mean)

fig1c.e<-fig1c+ 
  geom_path(data = df_ell.fire.06, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse

e2007<-subset(spscoresall.mod.e, year==2007)
e2007$treatment<-as.factor(e2007$treatment)
levels(e2007$treatment)
#data for ellipse, in this case using the management factor
df_ell.fire.07 <- data.frame() #sets up a data frame before running the function.
for(g in levels(e2007$treatment)){
  df_ell.fire.07  <- rbind(df_ell.fire.07 , cbind(as.data.frame(with(e2007 [e2007$treatment==g,],
                             veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2007=aggregate(e2007[ ,c("NMDS1", "NMDS2")], 
                         list(group = e2007$treatment), mean)

fig1d.e<-fig1d+ 
  geom_path(data = df_ell.fire.07, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse, sepa

e2008<-subset(spscoresall.mod.e, year==2008)
e2008$treatment<-as.factor(e2008$treatment)
levels(e2008$treatment)
#data for ellipse, in this case using the management factor
df_ell.fire.08 <- data.frame() #sets up a data frame before running the function.
for(g in levels(e2008$treatment)){
  df_ell.fire.08  <- rbind(df_ell.fire.08 , cbind(as.data.frame(with(e2008 [e2008$treatment==g,],
                              veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2008=aggregate(e2008[ ,c("NMDS1", "NMDS2")], 
                         list(group = e2008$treatment), mean)

fig1e.e<-fig1e+ 
  geom_path(data = df_ell.fire.08, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse, sepa

#put new fig1 panel together
grid.arrange(vec1, fig1b.e, fig1c.e, fig1d.e, fig1e.e,  layout_matrix = lay) #put panel together
#save as 1200W x 1800L


###new figure 2, post-grazing
l2008<-subset(spscoresall.mod.l, year==2008)
l2008$treatment<-as.factor(l2008$treatment)
levels(l2008$treatment)
#data for ellipse, in this case using the management factor
df_ell.graze.08 <- data.frame() #sets up a data frame before running the function.
for(g in levels(l2008$treatment)){
  df_ell.graze.08  <- rbind(df_ell.graze.08 , cbind(as.data.frame(with(l2008 [l2008$treatment==g,],
                                                               veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2008=aggregate(l2008[ ,c("NMDS1", "NMDS2")], 
                         list(group = l2008$treatment), mean)

fig2b.e<-fig2b+ 
  geom_path(data = df_ell.graze.08, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse, separate ones by treatment 
fig2b.e

l2009<-subset(spscoresall.mod.l, year==2009)
l2009$treatment<-as.factor(l2009$treatment)
levels(l2009$treatment)
#data for ellipse, in this case using the management factor
df_ell.graze.09 <- data.frame() #sets up a data frame before running the function.
for(g in levels(l2009$treatment)){
  df_ell.graze.09  <- rbind(df_ell.graze.09 , cbind(as.data.frame(with(l2009 [l2009$treatment==g,],
                                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2009=aggregate(l2009[ ,c("NMDS1", "NMDS2")], 
                         list(group = l2009$treatment), mean)

fig2c.e<-fig2c+ 
  geom_path(data = df_ell.graze.09, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse
fig2c.e

l2010<-subset(spscoresall.mod.l, year==2010)
l2010$treatment<-as.factor(l2010$treatment)
levels(l2010$treatment)
#data for ellipse, in this case using the management factor
df_ell.graze.10 <- data.frame() #sets up a data frame before running the function.
for(g in levels(l2010$treatment)){
  df_ell.graze.10  <- rbind(df_ell.graze.10 , cbind(as.data.frame(with(l2010 [l2010$treatment==g,],
                                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2010=aggregate(l2010[ ,c("NMDS1", "NMDS2")], 
                         list(group = l2010$treatment), mean)

fig2d.e<-fig2d+ 
  geom_path(data = df_ell.graze.10, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse, sepa
fig2d.e

l2011<-subset(spscoresall.mod.l, year==2011)
l2011$treatment<-as.factor(l2011$treatment)
levels(l2011$treatment)
#data for ellipse, in this case using the management factor
df_ell.graze.11 <- data.frame() #sets up a data frame before running the function.
for(g in levels(l2011$treatment)){
  df_ell.graze.11  <- rbind(df_ell.graze.11 , cbind(as.data.frame(with(l2011 [l2011$treatment==g,],
                                                                     veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2011=aggregate(l2011[ ,c("NMDS1", "NMDS2")], 
                         list(group = l2011$treatment), mean)

fig2e.e<-fig2e+ 
  geom_path(data = df_ell.graze.11, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse, sepa
fig2e.e

l2012<-subset(spscoresall.mod.l, year==2012)
l2012$treatment<-as.factor(l2012$treatment)
levels(l2012$treatment)
#data for ellipse, in this case using the management factor
df_ell.graze.12 <- data.frame() #sets up a data frame before running the function.
for(g in levels(l2012$treatment)){
  df_ell.graze.12  <- rbind(df_ell.graze.12, cbind(as.data.frame(with(l2012 [l2012$treatment==g,],
                                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,treatment=g))
}

# data for labelling the ellipse
NMDS.mean.2012=aggregate(l2012[ ,c("NMDS1", "NMDS2")], 
                         list(group = l2012$treatment), mean)

fig2f.e<-fig2f+ 
  geom_path(data = df_ell.graze.12, aes(x = NMDS1, y = NMDS2, group = treatment)) #this is the ellipse, sepa
fig2f.e

#put new fig1 panel together
grid.arrange(vec2, fig2b.e, fig2c.e, fig2d.e, fig2e.e, fig2f.e, layout_matrix = lay2) #put panel together
#save as 1200W x 2350L
