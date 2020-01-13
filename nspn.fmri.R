# Analysis code for NSPN manuscript "Conservative and disruptive modes of adolescent change in brain functional connectivity"
# by František Váša, Rafael Romero-Garcia, Manfred G. Kitzbichler, Jakob Seidlitz, Kirstie J. Whitaker, Matilde M. Vaghi, Prantik Kundu, 
# Ameera X. Patel, Peter Fonagy, Raymond J. Dolan, Peter B. Jones, Ian M. Goodyer, the NSPN Consortium, Petra E. Vértes* and Edward T. Bullmore*.
#
# František Váša, January 2020
# fdv247@gmail.com

# library dependencies (all must be installed)

# lme functions
library(nlme)         # linear mixed-effect models
library(MuMIn)        # function "r.squaredGLMM" (Pseudo-R-Squared For Generalized Mixed-Effect Models)

# plotting functions
library(ggplot2)      # lme plotting
library(AICcmodavg)   # function "predictSE.lme" for plotting of SE's on lme
library(seqinr)       # function "col2alpha" for color transparency
library(vioplot)      # violin plot
library(RColorBrewer) # color palettes
library(viridis)      # color palettes
library(hexbin)       # hexagonal binning of scatterplots

# other
library(moments)      # skewness and kurtosis (of FC distributions)

# code dependencies
code.path = '~/Desktop' # modify depending on where the dependencies were downloaded to
source(paste(code.path,'/dependencies/rp.main.R',sep=''))           # formats correlations and p-values for plot titles
source(paste(code.path,'/dependencies/write.fMRI.subset.R',sep='')) # writes out text files for plotting using KJW's visualisation code (see link below)
source(paste(code.path,'/dependencies/auto.lim.R',sep=''))          # automatic plot limits
source(paste(code.path,'/dependencies/perm.sphere.p.R',sep=''))     # "P_spin" p-value from permutation of the sphere
source(paste(code.path,'/dependencies/subc.plot.R',sep=''))         # subcortical region plot < >

# Cortical plots in the manuscript were made using Pysurfer-based cortical visualisation code written by Kirstie J. Whitaker.
# For details, see https://github.com/KirstieJane/DESCRIBING_DATA/wiki/Making-Surface-Plots-with-Pysurfer

##### Released data:
# Data are divided into (I) "specific" files containing FC matrices and participant demographics (specific to denoising/analysis stream: (i) main, (ii) low-motion, (iii) GSR) 
#                   and (II) a "general" file containing important auxiliary variables such as region ID's, names, etc

### (I)   Specific
#   (i)   nspn.fmri.main.RData
#   (ii)  nspn.fmri.lowmot.RData
#   (iii) nspn.fmri.gsr.RData

## Values
# age           subject ages
# fc            functional connectivity matrices [N(ROI) x N(ROI) x N(subj)]
# fd            average frame-wise displacement (in "main" data, this was regressed from FC, edge-wise)
# id            subject ID's
# male          sex (male / female)
# site          scanner site

### (II) General
# nspn.fmri.general.vars.RData

## Data
# cort.map      cortical maps (of cortical structure, PET, gene expression etc) used for relationships with maturational index (for detailed description of variables, see relevant section of code below)
# nm            region names (all 376 regions)
# perm.id       regional permutation ID's for spherical permutation (spin) test [N(cortical-ROI) x N(perm)] [here, N(perm) = 10'000]
# pos           regional centroid coordinates [N(ROI) x 3 (= x,y,z)]

## Values
# cort.map.nm   names of external cortical maps (from variable cort.map; for detailed description of variables, see relevant section of code below)
# hcp.346.cort  which of 346 retained regions are cortical
# hcp.346.subc  which of 346 retained regions are subcortical
# hcp.drop.id   which regions (in 1:376 = all cortical + all subcortical) are dropouts
# hcp.keep.id   which regions (in 1:376 = all cortical + all subcortical) are retained
# nm.simple     region names - "simple" format
# nm.subc       subcortical names (single hemisphere, for subcortical plots)
# nroi          number of regions of interest (after dropout exclusion)
# ve.id.subc    id's of von economo cytoarchitectonic classes, with subcortex (all 376 regions)
# yeo.id.subc   id's of yeo networks, with subcortex (all 376 regions)

# Data are available for download from [insert figshare link here]
data.path = '~/Desktop' # modify depending on where the data was downloaded to

# Choose analysis stream
data.stream = 'main'
#data.stream = 'lowmot'
#data.stream = 'gsr'

# load data
load(paste(data.path,'/nspn.fmri.',data.stream,'.RData',sep='')) # stream-specific variables
load(paste(data.path,'/nspn.fmri.general.vars.RData',sep=''))    # general variables

### rename stream-specific variables prior to running analysis script, and clear original variables
if (data.stream=='main') { # main data stream
  age = age.main
  fc = fc.main
  id = id.main
  male = male.main
  site = site.main
  fd = fd.main
  rm(age.main,fc.main,id.main,male.main,site.main,fd.main)
} else if (data.stream=='lowmot') { # low-motion data stream
  age = age.lowmot
  fc = fc.lowmot
  id = id.lowmot
  male = male.lowmot
  site = site.lowmot
  fd = fd.lowmot
  rm(age.lowmot,fc.lowmot,id.lowmot,male.lowmot,site.lowmot,fd.lowmot)
} else if (data.stream=='gsr') { # gsr data stream
  age = age.gsr
  fc = fc.gsr
  id = id.gsr
  male = male.gsr
  site = site.gsr
  fd = fd.gsr
  rm(age.gsr,fc.gsr,id.gsr,male.gsr,site.gsr,fd.gsr)
}

# define path for plot output (and create necessary folders)
plot.root = '~/Desktop' # path to folder
plot.path = paste(plot.root,'/nspn_fmri_',data.stream,sep='')
# master folder (nspn_fmri_[stream]) + folder for text files for brain-surface plots (surf_txt) including individual subcortical regions (/surf_txt/subc_ind/)
if (!dir.exists(plot.path)) dir.create(paste(plot.path,'/surf_txt/subc_ind/',sep=''),showWarning=T, recursive=T)  

### additional code dependency for spherical permutation "spin" test
# 1) either, use function "perm.sphere.p.R" with pre-generated permutations "perm.id$perm.all"
use.pre.perm.id = T
# 2) or, (re)generate your own permutations with "rotate.parcellation.R" (and "sphere_HCP.txt" data)
# functions available from https://github.com/frantisekvasa/rotate_parcellation
if (!use.pre.perm.id) {
  library(matrixStats)              # rowMins function, for "rotate.parcellation"
  
  # read-in spherical coordinates
  hcp.centroids = as.matrix(read.table('sphere_HCP.txt')) # available from https://github.com/frantisekvasa/rotate_parcellation
  
  # obtain coordinates of regions on L and R hemispheres
  coord.l = hcp.centroids[1:180,];              # coordinates of L hemisphere
  coord.r = hcp.centroids[181:360,];            # coordinates of R hemisphere
  keep.l = hcp.360.keep[hcp.360.keep<=180];     # exclude coordinates of "dropout" regions - L
  keep.r = hcp.360.keep[hcp.360.keep>180]-180;  # exclude coordinates of "dropout" regions - R
  # recreate coords and sizes
  coord.l = coord.l[keep.l,];
  coord.r = coord.r[keep.r,];
  
  # run spherical permutation script
  perm.id = rotate.parcellation(coord.l,coord.r,nrot=10000)
}

### Additional set-up and variable definition

# subset of retained von economo classes and yeo networks (with subcortex as the eighth class / network)
ve.id.subc.keep = ve.id.subc[hcp.keep.id]     # von economo ID's
nve = max(ve.id.subc.keep)                    # number of von economo classes
col.ve = c('darkmagenta','blue1','forestgreen','orange1','gold1','cyan','magenta','grey25') # von economo colors

yeo.id.subc.keep = yeo.id.subc[hcp.keep.id]   # yeo ID's
nyeo = max(yeo.id.subc.keep)                  # number of yeo networks
col.yeo = c('darkmagenta','cadetblue','forestgreen','hotpink2','cornsilk1','orange','salmon','grey25') # yeo colors

ns = length(age)                                        # number of subjects                          
triup = upper.tri(matrix(nrow=nroi,ncol=nroi))          # indices for upper triangular part of matrix
age.pred = seq(from=min(age),to=max(age),length.out=ns) # predicted age (for plotting of linear models)

# additional ID's (for writing of text files, for cortical surface plotting)
nroi.subc.tot = 16  # total subcortical regions
nroi.cort.tot = 360 # total subcortical regions
hcp.360.keep = hcp.keep.id[(nroi.subc.tot+1):nroi]-nroi.subc.tot # which of 1:360 cortical regions are kept?

# blue-red colobar 
x=matrix(rnorm(100),nrow=10)*100; xmin=0; xmax=100; x[x<xmin]=xmin; x[x>xmax]=xmax;
collist <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0","#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
ColorRamp<-colorRampPalette(collist)(10000); ColorLevels<-seq(from=xmin, to=xmax, length=10000)
bluered.col <- ColorRamp[round(1+(min(x)-xmin)*10000/(xmax-xmin)) : round( (max(x)-xmin)*10000/(xmax-xmin) )]
rm(x,xmin,xmax,ColorRamp,ColorLevels) # clear set-up variables

# decimal floor and ceiling functions for automatic plot limits
# https://stackoverflow.com/questions/35807523/r-decimal-ceiling
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level) # decimal floor
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level) # decimal ceiling

##### Analyses

### Global maturation - moments of correlation distribution 

# mean, variance, skewness and kurtosis
fc.m = fc.var = fc.skew = fc.kurt = vector(length=ns)
for (n in 1:ns) {
  if (n%%30==0) print(paste('participant ',toString(n),' of ',toString(ns),sep='')) # track progress
  fc.m[n] = mean(fc[,,n][triup])          # mean
  fc.var[n] = var(fc[,,n][triup])         # variance
  fc.skew[n] = skewness(fc[,,n][triup])   # skewness
  fc.kurt[n] = kurtosis(fc[,,n][triup])   # kurtosis
}

## head motion QC

# mean subject motion (mean FD) as a function of age (main = Fig. S2A)
df = data.frame(x1 = age, x2 = male, x3 = site, y = fd, id = id)
l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
newdat = expand.grid(x3 = 'wbic', x2='male',x1=range(age))
pred.m = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='male',x1=range(age)), se.fit=T, level=0)
pred.f = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='female',x1=range(age)), se.fit=T, level=0)
pred = list(); pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
pdf(paste(plot.path,'/age_vs_head_motion_lme.pdf',sep=''),width=6,height=5)
ggplot(df, aes(x=x1, y=y)) + #, colour=x2   
  labs(x = 'age (years)', y = expression(paste(mu,' FD (mm)',sep='')), title = bquote(r^2*' = '*.(toString(signif(r.squaredGLMM(l)[1,1],2)))*', '*p[age]*' = '*.(toString(signif(summary(l)$tTable[2,5],2))))) +
  geom_point(size=1, colour = 'grey') +
  geom_path(aes(group = id), colour = 'grey') + # spaghetti plot
  geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
  geom_line(data=newdat, aes(y=pred$fit),size=1) + #, size="Population")) +
  theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #, axis.line = element_line(colour = "black")
dev.off()

# mean FC as a function of mean subject motion (mean FD) (main = Fig. S2F)
# (! uncorrected imperfection spotted during consolidation of code: in "main" stream, to exactly reproduce statistics in Fig. S2F, remove "x3='wbic'" from lines 214-216 below)
df = data.frame(x1 = fd, x2 = male, x3 = site, y = fc.m, id = id)
l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
newdat = expand.grid(x3='wbic',x2='male',x1=range(fd))
pred.m = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='male',x1=range(fd)), se.fit=T, level=0)
pred.f = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='female',x1=range(fd)), se.fit=T, level=0)
pred = list();pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
pdf(paste(plot.path,'/head_motion_vs_mean_fc_lme.pdf',sep=''),width=6,height=5)
ggplot(df, aes(x=x1, y=y)) + #, colour=x2
  labs(x = expression(paste(mu, ' FD (mm)',sep='')), y = expression(paste(mu, ' FC',sep='')), title = bquote(r^2*' = '*.(toString(signif(r.squaredGLMM(l)[1,1],2)))*', '*p[FD]*' = '*.(toString(signif(summary(l)$tTable[2,5],2))))) +
  geom_point(size=1, colour = 'grey') +
  geom_path(aes(group = id), colour = 'grey') + # spaghetti plot
  geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
  geom_line(data=newdat, aes(y=pred$fit),size=1) + #, size="Population")) +
  theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

## plot moments of correlation distribution (main = Fig. S3A)
# for the models below, relevant statistics (embedded in plots in the paper) can be obtained using:
# t_age (age = "x1"): summary(l)$tTable[2,4]
# p_age (age = "x1"): summary(l)$tTable[2,5]

# mean
df = data.frame(x1 = age, x2 = male, x3 = site, y = fc.m, id = id)
l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
newdat = expand.grid(x3='wbic',x2='male',x1=range(age))
pred.m = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='male',x1=range(age)), se.fit=T, level=0)
pred.f = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='female',x1=range(age)), se.fit=T, level=0)
pred = list(); pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
pdf(paste(plot.path,'/fc_mean_lme.pdf',sep=''),width=6,height=5)
ggplot(df, aes(x=x1, y=y)) + 
  labs(x = 'age (years)', y = 'mean', title = '') + 
  geom_point(size=1, colour = 'grey') +
  geom_path(aes(group = id), colour = 'grey') + 
  geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
  geom_line(data=newdat, aes(y=pred$fit),size=1) + 
  theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

# variance
df = data.frame(x1 = age, x2 = male, x3 = site, y = fc.var, id = id)
l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
newdat = expand.grid(x2='male',x1=range(age))
pred.m = predictSE.lme(l, newdata=expand.grid(x3='wbic', x2='male',x1=range(age)), se.fit=T, level=0)
pred.f = predictSE.lme(l, newdata=expand.grid(x3='wbic', x2='female',x1=range(age)), se.fit=T, level=0)
pred = list(); pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
pdf(paste(plot.path,'/fc_variance_lme.pdf',sep=''),width=6,height=5)
ggplot(df, aes(x=x1, y=y)) + 
  labs(x = 'age (years)', y = 'variance', title = '') + 
  geom_point(size=1, colour = 'grey') +
  geom_path(aes(group = id), colour = 'grey') + 
  geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
  geom_line(data=newdat, aes(y=pred$fit),size=1) + 
  theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# skewness
df = data.frame(x1 = age, x2 = male, x3 = site, y = fc.skew, id = id)
l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
newdat = expand.grid(x3='wbic',x2='male',x1=range(age))
pred.m = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='male',x1=range(age)), se.fit=T, level=0)
pred.f = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='female',x1=range(age)), se.fit=T, level=0)
pred = list(); pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
pdf(paste(plot.path,'/fc_skewness_lme.pdf',sep=''),width=6,height=5)
ggplot(df, aes(x=x1, y=y)) + 
  labs(x = 'age (years)', y = 'skewness', title = '') + 
  geom_point(size=1, colour = 'grey') +
  geom_path(aes(group = id), colour = 'grey') + 
  geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
  geom_line(data=newdat, aes(y=pred$fit),size=1) + 
  theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# kurtosis
df = data.frame(x1 = age, x2 = male, x3 = site, y = fc.kurt, id = id)
l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
newdat = expand.grid(x3='wbic',x2='male',x1=range(age))
pred.m = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='male',x1=range(age)), se.fit=T, level=0)
pred.f = predictSE.lme(l, newdata=expand.grid(x3='wbic',x2='female',x1=range(age)), se.fit=T, level=0)
pred = list(); pred$fit = (pred.m$fit+pred.f$fit)/2; pred$se.fit = (pred.m$se.fit+pred.f$se.fit)/2;
pdf(paste(plot.path,'/fc_kurtosis_lme.pdf',sep=''),width=6,height=5)
ggplot(df, aes(x=x1, y=y)) + 
  labs(x = 'age (years)', y = 'kurtosis', title = '') + 
  geom_point(size=1, colour = 'grey') +
  geom_path(aes(group = id), colour = 'grey') + 
  geom_ribbon(data=newdat, aes(x=x1, ymin=pred$fit-2*pred$se.fit, ymax=pred$fit+2*pred$se.fit), alpha=0.5,fill="grey", inherit.aes=F) +
  geom_line(data=newdat, aes(y=pred$fit),size=1) + 
  theme_bw(base_size=22) + theme(plot.title = element_text(size=24, hjust = 0.5)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

### percentiles of correlation distribution (main = Fig. S3B)
prct.seq = seq(from=0,to=1,by=0.01) # percentiles to evaluate
nprct = length(prct.seq)            # number of percentiles
prct.col = viridis(101)             # color palette

fc.prct = array(NA,dim=c(nprct,ns));                            # intialise FC at each percentile
for (i in 1:ns) fc.prct[,i] = quantile(fc[,,i][triup],prct.seq) # calculate

# fit mixed-effect linear model and extract parameters
fc.prct.int = fc.prct.age = vector(length=nprct) # initialise variables
for (i in 1:nprct) {
  if (i%%10 == 0) print(paste('percentile ',toString(i),' of ',toString(nprct),sep='')) # track progress
  
  df = data.frame(x1 = age, x2 = male, x3 = site, y = fc.prct[i,], id = id) # data frame for linear mixed-effects model (lme)
  l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df) # lme fitting
  fc.prct.int[i] = l$coefficients$fixed[1]            # intercept
  fc.prct.age[i] = l$coefficients$fixed[2]            # beta coefficient for effect of age
}

# plot *with centering*
pdf(paste(plot.path,'/fc_prct_emp_oneplot_centered_prct_col.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
plot(0, type='n', yaxt='n', xlab='age (years)', ylab=expression(paste('lm(corr. prct ~ age)',sep='')), main='', xlim=range(age.pred), ylim=c(-0.1,0.1))
axis(side=2,at=seq(-0.1,0.1,by=0.1)); axis(side=2,at=seq(-0.05,0.05,by=0.1)) 
for (i in 1:nprct) abline(a = fc.prct.int[i]-c( fc.prct.int[i] + fc.prct.age[i]*mean(range(age)) ), b = fc.prct.age[i], col=prct.col[i], lwd=3)
rect(xleft = 13, ybottom = -0.13, xright = min(age), ytop = 0.13, col = 'white', border = NA) # white rectangles so that lines span only age range (left)
rect(xleft = 27, ybottom = -0.13, xright = max(age), ytop = 0.13, col = 'white', border = NA) # white rectangles so that lines span only age range (right)
box(which = "plot", lty = "solid") # replot plot box
dev.off()

# # plot *without centering*
# pdf(paste(plot.path,'/fc_prct_emp_oneplot_prct_col.pdf',sep=''),width=6,height=5)
# par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
# plot(0, type='n', xlab='age (years)', ylab=expression(paste('lm(corr. prct ~ age)',sep='')), main='', xlim=range(age.pred), ylim=c(-0.8,2.5))
# for (i in 1:nprct) abline(a = fc.prct.int[i], b = fc.prct.age[i], col=prct.col[i], lwd=3)
# rect(xleft = 13, ybottom = -1, xright = min(age), ytop = 3, col = 'white', border = NA) # white rectangles so that lines span only age range (left)
# rect(xleft = 27, ybottom = -1, xright = max(age), ytop = 3, col = 'white', border = NA) # white rectangles so that lines span only age range (right)
# box(which = "plot", lty = "solid") # replot plot box
# dev.off()

### Regional analyses (main = Fig. 1)

# split each region into cortical and subcortical node strength
str.cort = str.subc = array(NA,dim=c(nroi,ns))
for (n in 1:ns) {
  str.cort[,n] = rowMeans(fc[,,n][,hcp.346.cort],na.rm=T) # cortical
  str.subc[,n] = rowMeans(fc[,,n][,hcp.346.subc],na.rm=T) # subcortical
}

# strength for individual subcortical regions (averaged between left and right)
nsubc.hemi = length(hcp.346.subc)/2 # number of subcortical regions (in each hemisphere)
str.indsubc = array(NA,dim=c(nroi,ns,nsubc.hemi))
for (i in 1:ns) {
  for (n in 1:nsubc.hemi) str.indsubc[,i,n] = rowMeans(fc[,,i][,c(n,n+nsubc.hemi)],na.rm=T)
}

# nodal trajectories
# cortex-all + subcortex-all
str.cort.14 = str.cort.age.p = str.cort.age.t = vector(length = nroi) # cortical node strength
str.subc.14 = str.subc.age.p = str.subc.age.t = vector(length = nroi) # subcortical node strength
cort.conv.error = subc.conv.error = c()                               # catch errors for any node where lme does not converge  (esp. in smaller (eg: low-motion) subsets of data)
for (n in 1:nroi) {
  if (n%%20 == 0) print(paste('region ',toString(n),' of ',toString(nroi),sep='')) # track progress
  
  # str.cort
  tryCatch(
    expr = {
      # fit lme
      df = data.frame(x1 = age, x2 = male, x3 = site, y = str.cort[n,], id = id)  # data frame for lme
      l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)                       # fit lme
      str.cort.age.t[n] = summary(l)$tTable[2,4]                                  # t-statistic for effect of age (dFC_14-26)
      str.cort.age.p[n] = summary(l)$tTable[2,5]                                  # p-value for effect of age
      # predicted value at 14 (FC_14; includes "average" effect of sex (0.5), and average effect of each of three scanner sites (1/3))
      str.cort.14[n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*0.5 + l$coefficients$fixed[4]*(1/3) + l$coefficients$fixed[5]*(1/3)
    },error=function(e){
      # in case of non-convergence, store region that did not converge (+ replace values by "NA" below)
      cort.conv.error <<- c(cort.conv.error,n) # <<- force assign (https://stackoverflow.com/questions/21956031/trycatch-does-not-seem-to-return-my-variable)
      cat('cortical n =',n,':',conditionMessage(e), '\n')
    })
  
  # str.subc
  tryCatch(
    expr = {
      # fit lme
      df = data.frame(x1 = age, x2 = male , x3 = site, y = str.subc[n,], id = id) # data frame for lme
      l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)                       # fit lme
      str.subc.age.t[n] = summary(l)$tTable[2,4]                                  # t-statistic for effect of age (dFC_14-26)
      str.subc.age.p[n] = summary(l)$tTable[2,5]                                  # p-value for effect of age
      # predicted value at 14 (FC_14; includes "average" effect of sex (0.5), and average effect of each of three scanner sites (1/3))
      str.subc.14[n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*0.5 + l$coefficients$fixed[4]*(1/3) + l$coefficients$fixed[5]*(1/3)
    },error=function(e){
      # in case of non-convergence, store region that did not converge (+ replace values by "NA" below)
      subc.conv.error <<- c(subc.conv.error,n) # <<- force assign
      cat('subcortical n =',n,':',conditionMessage(e), '\n')
    })
}
# replace any non-converged values by "NA" (P=1)
# cort
if (!is.null(cort.conv.error)) {
  str.cort.age.t[cort.conv.error] = NA
  str.cort.age.p[cort.conv.error] = 1
  str.cort.14[cort.conv.error] = NA
}
# subc
if (!is.null(subc.conv.error)) {
  str.subc.age.t[subc.conv.error] = NA
  str.subc.age.p[subc.conv.error] = 1
  str.subc.14[subc.conv.error] = NA
}

# as above for connectivity of individual subcortical regions to the rest of the brain
str.indsubc.14 = str.indsubc.age.p = str.indsubc.age.t = array(NA,dim=c(nroi,nsubc.hemi))
indsubc.conv.error.n = indsubc.conv.error.i = c() # catch errors for any node where lme does not converge  (esp. in smaller (eg: low-motion) subsets of data)
for (n in 1:nroi) {
  if (n%%10 == 0) print(paste('region ',toString(n),' of ',toString(nroi),sep='')) # track progress
  for (i in 1:nsubc.hemi) {   # loop over subcortical regions
  
    tryCatch(
      expr = {
        # fit lme
        df = data.frame(x1 = age, x2 = male, x3 = site, y = str.indsubc[n,,i], id = id)   # data frame for lme
        l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)                             # fit lme
        str.indsubc.age.t[n,i] = summary(l)$tTable[2,4]                                   # t-statistic for effect of age (dFC_14-26)
        str.indsubc.age.p[n,i] = summary(l)$tTable[2,5]                                   # p-value for effect of age
        # predicted value at 14 (FC_14; includes "average" effect of sex (0.5), and average effect of each of three scanner sites (1/3))
        str.indsubc.14[n,i] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*0.5 + l$coefficients$fixed[4]*(1/3) + l$coefficients$fixed[5]*(1/3)
      },error=function(e){
        # in case of non-convergence, store id's for regions that did not converge (+ replace values by "NA" below)
        indsubc.conv.error.n <<- c(indsubc.conv.error.n,n) # <<- force assign
        indsubc.conv.error.i <<- c(indsubc.conv.error.i,i)
        cat('individual subcortical n =',n,', i =',i,':',conditionMessage(e), '\n')
      })
  }
}
# replace any non-converged values by "NA" (P=1)
if (!is.null(indsubc.conv.error.n)) {
  for (j in length(indsubc.conv.error.n)) {
    str.indsubc.age.t[indsubc.conv.error.n[j],indsubc.conv.error.i[j]] = NA
    str.indsubc.age.p[indsubc.conv.error.n[j],indsubc.conv.error.i[j]] = 1
    str.indsubc.14[indsubc.conv.error.n[j],indsubc.conv.error.i[j]] = NA
  }
}

# correct p-values (cort and subc) for multiple comparisons (FDR)
str.cort.age.p.fdr = p.adjust(str.cort.age.p, method = 'fdr')   # cortex-all
str.subc.age.p.fdr = p.adjust(str.subc.age.p, method = 'fdr')   # subcortex-all
str.indsubc.age.p.fdr = array(NA,dim=c(nroi,nsubc.hemi))        # individual subcortical regions
for (i in 1:nsubc.hemi) str.indsubc.age.p.fdr[,i] = p.adjust(str.indsubc.age.p[,i], method = 'fdr')

### text files for cortical surface plots
# For details, see https://github.com/KirstieJane/DESCRIBING_DATA/wiki/Making-Surface-Plots-with-Pysurfer
# Following installation, surface plots can be generated with the following example *bash* command:
# pysurfer_plot_parcellation_surface_values.py --annot_name HCP -sd [required_data_directory_released_with_pysurfer_code] --cmap plasma_r -s inflated ~/Desktop/nspn_fmri/surf_txt/str_cort_14.txt -l -0.2 -u 0.6
# Plot files are saved in a "PNGS" folder, in directory where plotting text files are located (in this case, '~/Desktop/nspn_fmri/surf_txt/PNGS')

# FC_14 (pysurfer cmap: plasma_r)
write.fMRI.subset(str.cort.14[hcp.346.cort], hcp.360.keep, nroi.cort.tot, paste(plot.path,'/surf_txt/str_cort_14.txt',sep=''))
write.fMRI.subset(str.subc.14[hcp.346.cort], hcp.360.keep, nroi.cort.tot, paste(plot.path,'/surf_txt/str_subc_14.txt',sep=''))

# dFC_14-26 (effect of age, t-statistic) (pysurfer cmap: seismic)
write.fMRI.subset(str.cort.age.t[hcp.346.cort], hcp.360.keep, nroi.cort.tot, paste(plot.path,'/surf_txt/str_cort_age_t.txt',sep='')) # cortex-all
write.fMRI.subset(str.subc.age.t[hcp.346.cort], hcp.360.keep, nroi.cort.tot, paste(plot.path,'/surf_txt/str_subc_age_t.txt',sep='')) # subcortex-all

# subcortex ind.
for (i in 1:nsubc.hemi) {
  # FC_14
  write.fMRI.subset(str.indsubc.14[hcp.346.cort,i], hcp.360.keep, nroi.cort.tot, paste(plot.path,'/surf_txt/subc_ind/str_subc_',nm.subc[i],'_14.txt',sep=''))
  # dFC_14-26 (effect of age, t-statistic)
  write.fMRI.subset(str.indsubc.age.t[hcp.346.cort,i], hcp.360.keep, nroi.cort.tot, paste(plot.path,'/surf_txt/subc_ind/str_subc_',nm.subc[i],'_age_t.txt',sep=''))
}

### Average subcortical plots (Fig. S4)

# order subcortical regions according to average rate of change as a function of age
subc.ord = sort(colMeans(str.indsubc.age.t,na.rm=T),decreasing=T,index.return=T)$ix
#subc.ord = c(6,5,8,3,2,7,4,1) # order used throughout manuscript, based on "main" stream analysis

# "col.l" and "col.u" parameters (setting colorbar range and x-axis limits) can be manually adjusted for nicer plots (default values use decimal floor and ceiling functions)
# legend placement also needs to be adjusted

# FC_14 (main = Fig. S4)

# subc-cort
pdf(paste(plot.path,'/str_14_cort_subc.pdf',sep=''),width=4,height=3)
par(mar=c(3, 5, 1, 2) + 0.1, cex.lab = 2, cex.axis = 1.3, cex.main = 2, font.main = 1, bg='white')
subc.plot(str.cort.14[hcp.346.subc],subc.ord,nm.subc,colBar = rev(plasma(256)),col.l = min(floor_dec(str.cort.14[hcp.346.subc],1)),col.u = max(ceiling_dec(str.cort.14[hcp.346.subc],1)))
#legend('topright',title='hemis.',legend=c('left','right'),col='grey60',text.col='grey30',pch=c(60,62),cex=1.2,pt.cex=1.35,bty='n')
dev.off()

# subc-subc
pdf(paste(plot.path,'/str_14_subc_subc.pdf',sep=''),width=4,height=3)
par(mar=c(3, 5, 1, 2) + 0.1, cex.lab = 2, cex.axis = 1.3, cex.main = 2, font.main = 1, bg='white')
subc.plot(str.subc.14[hcp.346.subc],subc.ord,nm.subc,colBar = rev(plasma(256)),col.l = min(floor_dec(str.subc.14[hcp.346.subc],1)),col.u = max(ceiling_dec(str.subc.14[hcp.346.subc],1)))
#legend('bottomleft',title='hemis.',legend=c('left','right'),col='grey60',text.col='grey30',pch=c(60,62),cex=1.2,pt.cex=1.35,bty='n')
dev.off()

# dFC_14-26 (effect of age, t-statistic; main = Fig. S4)

# subc-cort
pdf(paste(plot.path,'/str_tstat_cort_subc.pdf',sep=''),width=4,height=3)
par(mar=c(3, 5, 1, 2) + 0.1, cex.lab = 2, cex.axis = 1.3, cex.main = 2, font.main = 1, bg='white')
subc.plot(str.cort.age.t[hcp.346.subc],subc.ord,nm.subc,colBar = colorRampPalette(c('blue','white','red'))(100),col.l = min(floor_dec(str.cort.age.t[hcp.346.subc],1)),col.u = max(ceiling_dec(str.cort.age.t[hcp.346.subc],1)))
#legend('topleft',title='hemis.',legend=c('left','right'),col='grey60',text.col='grey30',pch=c(60,62),cex=1.2,pt.cex=1.35,bty='n')
dev.off()

# subc-subc
pdf(paste(plot.path,'/str_tstat_subc_subc.pdf',sep=''),width=4,height=3)
par(mar=c(3, 5, 1, 2) + 0.1, cex.lab = 2, cex.axis = 1.3, cex.main = 2, font.main = 1, bg='white')
subc.plot(str.subc.age.t[hcp.346.subc],subc.ord,nm.subc,colBar = colorRampPalette(c('blue','white','red'))(100),col.l = min(floor_dec(str.subc.age.t[hcp.346.subc],1)),col.u = max(ceiling_dec(str.subc.age.t[hcp.346.subc],1)))
#legend('bottomright',title='hemis.',legend=c('left','right'),col='grey60',text.col='grey30',pch=c(60,62),cex=1.2,pt.cex=1.35,bty='n')
dev.off()

##### Fit lme at level of individual FC edges (takes ~1h+ to run)
fc.age.t = fc.age.p = fc.14 = array(0,dim=c(nroi,nroi))
fc.conv.error.i = fc.conv.error.j = c() # catch errors for any edge where lme does not converge (esp. in smaller (eg: low-motion) subsets of data) 
for (i in 1:(nroi-1)) {     # only loop over upper triangular (as matrices are symmetric), and then fill lower triangular with the transpose
  print(paste('region ',toString(i),' of ',toString(nroi),sep='')) # track progress
  for (j in (i+1):nroi) {   # as above - loop over upper trinagular only
  
    tryCatch(
      expr = {
        # fit lme
        df = data.frame(x1 = age, x2 = male, x3 = site, y = fc[i,j,], id = id)  # data frame for model fitting
        l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)                   # fit model
        fc.age.t[i,j] = summary(l)$tTable[2,4]                                  # t-statistic of effect of age
        fc.age.p[i,j] = summary(l)$tTable[2,5]                                  # p-value of effect of age
        # predicted value at 14 (includes "average" effect of sex (0.5), and average effect of each of three scanner sites (1/3))
        fc.14[i,j] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*0.5 + l$coefficients$fixed[4]*(1/3) + l$coefficients$fixed[5]*(1/3)  # y = mx+b (+ average effect of gender)
      },error=function(e){
        # in case of non-convergence, store id's for edges that did not converge (+ replace values by "NA" below)
        fc.conv.error.i <<- c(fc.conv.error.i,i) # <<- force assign
        fc.conv.error.j <<- c(fc.conv.error.j,j)
        cat('fc edge i =',i,', j =',j,':',conditionMessage(e), '\n')
      })
  }
}
# replace non-converged regions by NA (P=1), and reflect matrices to fill empty lower triangular parts
if (!is.null(fc.conv.error.i)) {
  for (j in length(fc.conv.error.i)) {
    fc.age.t[fc.conv.error.i[i],fc.conv.error.j[j]] = NA
    fc.age.p[fc.conv.error.i[i],fc.conv.error.j[j]] = 1
    fc.14[fc.conv.error.i[i],fc.conv.error.j[j]] = NA
  }
}
# fill empty lower triangulars of matrices of lme parameters
if (all(t(fc.age.t)[triup]==0)) { # reflection condition
  fc.age.t = fc.age.t + t(fc.age.t)
  fc.age.p = fc.age.p + t(fc.age.p)
  fc.14 = fc.14 + t(fc.14)
}

##### Maturational Index (MI): Relationship between edge-wise FC at 14 and change as a function of age at each node

# regional relationships (evaluated region is removed)
mi.rho = mi.p = vector(length=nroi)
for (i in 1:nroi) {
  mi.rho[i] = cor.test(fc.14[i,-i],fc.age.t[i,-i],method='spearman')$estimate   # spearman rho
  mi.p[i] = cor.test(fc.14[i,-i],fc.age.t[i,-i],method='spearman')$p.value      # p-value
}
mi.p.fdr = p.adjust(mi.p, method = 'fdr') # correct p-values for multiple comparisons

# spherical permutation ("spin") p-values (takes ~20-30 mins to run; evaluated region is left in, for computational tractability)
mi.pspin = vector(length=nroi)
for (i in 1:nroi) {
  if (i%%10==0) print(paste('region ',toString(i),' of ',toString(nroi),sep='')) # track progress
  mi.pspin[i] = perm.sphere.p(fc.14[i,][hcp.346.cort],fc.age.t[i,][hcp.346.cort],perm.id$perm.all,corr.type='spearman') 
}
mi.pspin.fdr = p.adjust(mi.pspin, method = 'fdr') # correct p-values for multiple comparisons

# plots of example regions (main = Fig. 2B)

# map maturational index to colorbar
colBar = colorRampPalette(c('blue','white','red'))(100)
mi.col.l = -0.8; mi.col.u = 0.8                           # lower and upper colorbar limits (might require adjustment, based on range(mi.rho))
temp = ceiling(100*((mi.rho-mi.col.l)/(mi.col.u-mi.col.l)));
col.slt_14 = matrix(colBar[temp])

# positive
n = which(nm.simple[hcp.keep.id]=='L 1') # left motor -> positive /  n = 67
int.pred = seq(from=0,to=1.5,length.out=nroi)
pdf(paste(plot.path,'/str_mat_index_example_L1.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
x = fc.14[n,-n]; l = lm(fc.age.t[n,-n]~x);
plot(0, type='n', xlab='FC at 14', ylab=expression(paste(Delta,FC,sep='')), xlim=c(range(int.pred)), ylim=c(-2.5,7)) # ylim might require adjustment
points(fc.14[n,-n],fc.age.t[n,-n],col='black',pch=19,cex=0.75)
pred = predict(l,newdata=data.frame('x'=int.pred),se.fit=T)
polygon(c(rev(int.pred), int.pred), c(rev(pred$fit+2*pred$se.fit), pred$fit-2*pred$se.fit), col = col2alpha('grey',alpha=0.5), border = NA)
lines(int.pred,pred$fit,col=col.slt_14[n],lwd=3);
dev.off()

# negative
n = which(nm.simple[hcp.keep.id]=='L RSC') # left cingulate -> negative /  n = 30
int.pred = seq(from=0,to=1.65,length.out=nroi)
pdf(paste(plot.path,'/str_mat_index_example_LRSC.pdf',sep=''),width=6,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
x = fc.14[n,-n]; l = lm(fc.age.t[n,-n]~x);
plot(0, type='n', xlab='FC at 14', ylab=expression(paste(Delta,FC,sep='')), xlim=c(range(int.pred)), ylim=c(-4.5,4.5)) # ylim might require adjustment
points(fc.14[n,-n],fc.age.t[n,-n],col='black',pch=19,cex=0.75)
pred = predict(l,newdata=data.frame('x'=int.pred),se.fit=T)
polygon(c(rev(int.pred), int.pred), c(rev(pred$fit+2*pred$se.fit), pred$fit-2*pred$se.fit), col = col2alpha('grey',alpha=0.5), border = NA)
lines(int.pred,pred$fit,col=col.slt_14[n],lwd=3);
dev.off()

# write out text files for surface plotting (main = Fig. 2C, and Fig. S5) (1 = all, 2 = P < 0.05, 3 = P_spin < 0.05)
write.fMRI.subset(mi.rho[hcp.346.cort], hcp.360.keep, nroi.cort.tot, paste(plot.path,'/surf_txt/mat_index.txt',sep=''))
write.fMRI.subset(mi.rho[hcp.346.cort][which(mi.p.fdr[hcp.346.cort]<0.05)], hcp.360.keep[which(mi.p.fdr[hcp.346.cort]<0.05)], nroi.cort.tot, paste(plot.path,'/surf_txt/mat_index_p_fdr.txt',sep=''))
write.fMRI.subset(mi.rho[hcp.346.cort][which(mi.pspin.fdr[hcp.346.cort]<0.05)], hcp.360.keep[which(mi.pspin.fdr[hcp.346.cort]<0.05)], nroi.cort.tot, paste(plot.path,'/surf_txt/mat_index_pspin_fdr.txt',sep=''))

# plots of individual subcortical regions (main = Fig. 2D, and Fig. S5)
# (! uncorrected imperfection spotted during consolidation of code: in "gsr" stream, subcortical regions in Fig. S31D are ordered according to their default order (different from analogous figures from other analysis streams))
# all
pdf(paste(plot.path,'/mat_index_subc.pdf',sep=''),width=5,height=3)
par(mar=c(3, 5, 1, 2) + 0.1, cex.lab = 2, cex.axis = 1.3, cex.main = 2, font.main = 1, bg='white')
subc.plot(mi.rho[hcp.346.subc],subc.ord,nm.subc,colBar = colorRampPalette(c('blue','white','red'))(100),col.l = mi.col.l,col.u = mi.col.u)
legend('bottomright',title='hemis.',legend=c('left','right'),col='grey60',text.col='grey30',pch=c(60,62),cex=1.2,pt.cex=1.35,bty='n')
dev.off()
# P < 0.05
pdf(paste(plot.path,'/mat_index_subc_p_fdr.pdf',sep=''),width=5,height=3)
par(mar=c(3, 5, 1, 2) + 0.1, cex.lab = 2, cex.axis = 1.3, cex.main = 2, font.main = 1, bg='white')
subc.plot.sig(mi.rho[hcp.346.subc],mi.p.fdr[hcp.346.subc],0.05,subc.ord,nm.subc,colBar = colorRampPalette(c('blue','white','red'))(100),col.l = mi.col.l,col.u = mi.col.u)
legend('bottomright',title='hemis.',legend=c('left','right'),col='grey60',text.col='grey30',pch=c(60,62),cex=1.2,pt.cex=1.35,bty='n')
dev.off()
# P_spin < 0.05
pdf(paste(plot.path,'/mat_index_subc_pspin_fdr.pdf',sep=''),width=5,height=3)
par(mar=c(3, 5, 1, 2) + 0.1, cex.lab = 2, cex.axis = 1.3, cex.main = 2, font.main = 1, bg='white')
subc.plot.sig(mi.rho[hcp.346.subc],mi.pspin.fdr[hcp.346.subc],0.05,subc.ord,nm.subc,colBar = colorRampPalette(c('blue','white','red'))(100),col.l = mi.col.l,col.u = mi.col.u)
legend('bottomright',title='hemis.',legend=c('left','right'),col='grey60',text.col='grey30',pch=c(60,62),cex=1.2,pt.cex=1.35,bty='n')
dev.off()

### MI in context

# violin plots, coloured by average maturational index

# function to plot colorbar
color.bar <- function(lut, alpha, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  lut = sapply(lut, col2alpha, alpha)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, at=ticks, labels=ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# von Economo w/ subcortex (main = Fig. 3A)

# calculate mean for each class, to color violin plots
ve.slt_14 = vector(length=nve)
for (i in 1:nve) ve.slt_14[i] = mean(c(mi.rho[ve.id.subc.keep==i]))

# map values to colorbar
colBar = colorRampPalette(c('blue','white','red'))(100)
col.l = -0.5; col.u = 0.5 # lower and upper range of colorbar (might require manual adjustment, depending on analysis stream)
col.slt_14 = matrix(colBar[ceiling(100*((ve.slt_14-col.l)/(col.u-col.l)))])

pdf(paste(plot.path,'/mat_ind_von_economo.pdf',sep=''),width=5,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0.5,nve+0.5), ylim=c(-0.8,0.8),xlab='vE class',ylab='Maturational Index (MI)')
for (i in 1:nve) axis(side=1,at=i,col.axis=col.ve[i])
axis(side=2,at=c(col.l,0,col.u))
abline(h=c(col.l,0,col.u),lty=2,col='grey')
for (i in 1:nve) vioplot(c(mi.rho[ve.id.subc.keep==i]),at=i,col=col.slt_14[i],add=T,wex=0.75) # violin plot
dev.off()

# colorbar plot
pdf(paste(plot.path,'/mat_ind_von_economo_cbar.pdf',sep=''),width=1.5,height=6)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
color.bar(colBar, alpha = 0.8, min = col.l, max =  col.u, nticks = 3)
dev.off()

# Yeo w/ subcortex (main = Fig. 3B)

# calculate mean for each Yeo network, to color violin plots
yeo.slt_14 = vector(length=nve)
for (i in 1:nve) yeo.slt_14[i] = mean(c(mi.rho[yeo.id.subc.keep==i]))

# map values to colorbar
colBar = colorRampPalette(c('blue','white','red'))(100)
col.l = -0.5; col.u = 0.5 # lower and upper range of colorbar (might require manual adjustment, depending on analysis stream)
col.slt_14 = matrix(colBar[ceiling(100*((yeo.slt_14-col.l)/(col.u-col.l)))])

pdf(paste(plot.path,'/mat_ind_yeo.pdf',sep=''),width=5,height=5)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0.5,nyeo+0.5), ylim=c(-0.8,0.8),xlab='Yeo network',ylab='Maturational Index (MI)')
for (i in 1:nyeo) axis(side=1,at=i,col.axis=col.yeo[i])
axis(side=2,at=c(col.l,0,col.u))
abline(h=c(col.l,0,col.u),lty=2,col='grey')
for (i in 1:nyeo) vioplot(c(mi.rho[yeo.id.subc.keep==i]),at=i,col=col.slt_14[i],add=T,wex=0.75) # violin plot
dev.off()

# colorbar plot
pdf(paste(plot.path,'/mat_ind_yeo_cbar.pdf',sep=''),width=1.5,height=6)
par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
color.bar(colBar, alpha = 0.8, min = col.l, max =  col.u, nticks = 3)
dev.off()

##### Relationship of edge-wise FC slopes to Euclidean distance (main = Fig. S7)
# i.e. test of "distance-dependence" of maturation

# Euclidean distance between regional centroids (using centroid variable "pos")
dist = matrix(nrow=nroi,ncol=nroi)
for (i in 1:nroi) {; for (j in 1:nroi) {; dist[i,j] = sqrt( sum( (pos[i,]-pos[j,])^2 ) ) ; }; }

# hexagonal 2D binning
bin = hexbin(dist[triup],fc.age.t[triup], xbins=55)
my_colors = colorRampPalette(rev(brewer.pal(11,'Spectral')))

# plot
pdf(paste(plot.path,'/fc_dist_vs_age_t_hexbin.pdf',sep=''),width=6,height=5)
par(mar=c(1, 1, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 2, font.main = 1, bg='white')
x = dist[triup]; l = lm(fc.age.t[triup]~x);
spear = cor.test(dist[triup],fc.age.t[triup],method='spearman'); sp.rho = spear$estimate; sp.p = spear$p.value
pl = plot(bin, colramp=my_colors , legend=F,xlab='distance (mm)', ylab=expression(paste(Delta,FC,sep=' ')),main=rp.main.sp(sp.rho,sp.p,2)) 
hexVP.abline(pl$plot.vp,l,col = 'white', lwd = 5)
hexVP.abline(pl$plot.vp,l,col = gray(.6), lwd = 2)
dev.off()

# plot colorbar
pdf(paste(plot.path,'/fc_dist_vs_age_t_hexbin_cbar.pdf',sep=''),width=1.5,height=6)
color.bar(rev(brewer.pal(11,'Spectral')), alpha = 0.8, min = min(bin@count), max =  max(bin@count), nticks = 6, ticks = rev(seq(min(bin@count), max(bin@count), len=3)))
dev.off()

### Comparison of maturational index to external cortical maps (main = Fig. 4, Fig. S9)
# For a detailed description of maps, and full references, see SI section "Contextualising findings" -> "Independent cortical maps" & "Gene expression analyses"

     ### structural development (from NSPN data; Whitaker*, Vértes* et al. PNAS 113(32), 2016 and Váša et al. Cereb. Cortex 28, 2018)
# 1) dCT              age-related changes (d) in cortical thickness (CT)
# 2) dMT_70           age-related changes (d) in magnetization transfer (MT) at 70% fractional depth into cortex (70)
# 3) dSC_CT           age-related changes (d) in node strength of structural covariance (SC) of CT across subjects
     ### PET measures (Vaishnavi et al PNAS 107(41) 2010)
# 4) glyc_ind         glycolytic index (measure of aerobic glycolysis)
# 5) CMRO2            metabolic rate of oxygen
# 6) CMRGlu           metabolic rate of glucose
# 7) CBV              cerebral blood volume
# 8) CBF              cerebral blood flow
     ### cortical expansion / scaling (Hill et al. PNAS 107(29), 2010 and Reardon*, Seidlitz* et al. Science 360(6394), 2018)
# 9) evol_exp         evolutionary expansion
# 10) postnat_exp     post-natal expansion
# 11) areal_scaling   scaling of cortical surface area
      ### expression of aerobic-glycolysis (AG) related genes (Goyal et al. Cell Metab. 19(2), 2014)
# 12) AG_gene_expr    mean regional expression genes related to aerobic glycolysis (AG) (left hemisphere only; right hemisphere = NA)

# (! uncorrected imperfection spotted during consolidation of code: in the "lowmot" stream, some statistics differ very slightly [O(0.01)] between this code and Fig. S28...)
# (...this is due to a very small numerical difference between the "lowmot MI" reported in the text, and the "lowmot MI" generated by this consolidated code (these two lowmot MI's correlate at r = 0.999994)

ncort.map = length(cort.map.nm)

# surface plots
for (n in 1:ncort.map) write.fMRI.subset(cort.map[,n][hcp.360.keep], hcp.360.keep, nroi.cort.tot, paste(plot.path,'/surf_txt/cort_map_',cort.map.nm[n],'_drop.txt',sep=''))

# i) perform statistical comparisons and store all P-values (for subsequent multiple comparisons correction)
rho.cort.map = p.cort.map = p.spin.cort.map = vector(length=ncort.map) # initialise P-value vector
for (n in 1:ncort.map) {
  print(paste(cort.map.nm[n],' (',toString(n),' of ',toString(ncort.map),')',sep='')) # track progress
  
  # for all maps except AG gene expression, use data from both hemispheres
  if (cort.map.nm[n]!='AG_gene_expr') {
    
    spear = cor.test(cort.map[,n],mi.rho[hcp.346.cort],method='spearman')
    rho.cort.map[n] = spear$estimate  # spearman rho
    p.cort.map[n] = spear$p.value     # spearman P
    p.spin.cort.map[n] = perm.sphere.p(cort.map[,n],mi.rho[hcp.346.cort],perm.id$perm.all,corr.type='spearman')  # P_spin

  # for AG gene expression, only use data from left hemisphere
  } else {
    
    lh.ind = which(!is.na(cort.map[,n])) # index of LH regions
    spear = cor.test(cort.map[,n][lh.ind],mi.rho[hcp.346.cort][lh.ind],method='spearman')
    rho.cort.map[n] = spear$estimate  # spearman rho
    p.cort.map[n] = spear$p.value     # spearman P                                                                             
    p.spin.cort.map[n] = perm.sphere.p(cort.map[,n][lh.ind],mi.rho[hcp.346.cort][lh.ind],perm.id$perm.all[lh.ind,],corr.type='spearman')   # P_spin
    
  }
}

# ii) FDR correction
p.cort.map.fdr = p.adjust(p.cort.map, method = 'fdr')
p.spin.cort.map.fdr = p.adjust(p.spin.cort.map, method = 'fdr')

# iii) plotting (main = Fig. 4 and Fig. S9)
for (n in 1:ncort.map) {
  
  # for all maps except AG gene expression, use data from both hemispheres
  if (cort.map.nm[n]!='AG_gene_expr') {
    
    pdf(paste(plot.path,'/mat_index_vs_cort_map_',cort.map.nm[n],'.pdf',sep=''),width=6,height=5)
    par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 1.85, font.main = 1, bg='white')
    x = mi.rho[hcp.346.cort]; l = lm(cort.map[,n]~x);
    plot(0, type='n', ylab=cort.map.nm[n], xlab='Maturational Index (MI)', main=rp.main.sp.spin(rho.cort.map[n],p.cort.map.fdr[n],p.spin.cort.map.fdr[n],2), xlim=c(-0.9,0.9), ylim=auto.lim(cort.map[,n]))
    points(mi.rho[hcp.346.cort],cort.map[,n],col='black',pch=19,cex=0.75)
    pred = predict(l,newdata=data.frame('x'=seq(-0.9,0.9,length.out=100)),se.fit=T)
    polygon(c(rev(seq(-0.9,0.9,length.out=100)), seq(-0.9,0.9,length.out=100)), c(rev(pred$fit+2*pred$se.fit), pred$fit-2*pred$se.fit), col = col2alpha('grey',alpha=0.5), border = NA)
    lines(seq(-0.9,0.9,length.out=100),pred$fit,col='black',lwd=3);
    dev.off()
    
  # for AG gene expression, only use data from left hemisphere
  } else {
    
    lh.ind = which(!is.na(cort.map[,n])) # index of LH regions
    pdf(paste(plot.path,'/mat_index_vs_cort_map_',cort.map.nm[n],'.pdf',sep=''),width=6,height=5)
    par(mar=c(5, 5, 3, 1) + 0.1, cex.lab = 2, cex.axis = 1.5, cex.main = 1.85, font.main = 1, bg='white')
    x = mi.rho[hcp.346.cort][lh.ind]; l = lm(cort.map[,n][lh.ind]~x);
    plot(0, type='n', ylab=cort.map.nm[n], xlab='Maturational Index (MI)', main=rp.main.sp.spin(rho.cort.map[n],p.cort.map.fdr[n],p.spin.cort.map.fdr[n],2), xlim=c(-0.9,0.9), ylim=auto.lim(cort.map[,n][lh.ind]))
    points(mi.rho[hcp.346.cort][lh.ind],cort.map[,n][lh.ind],col='black',pch=19,cex=0.75)
    pred = predict(l,newdata=data.frame('x'=seq(-0.9,0.9,length.out=100)),se.fit=T)
    polygon(c(rev(seq(-0.9,0.9,length.out=100)), seq(-0.9,0.9,length.out=100)), c(rev(pred$fit+2*pred$se.fit), pred$fit-2*pred$se.fit), col = col2alpha('grey',alpha=0.5), border = NA)
    lines(seq(-0.9,0.9,length.out=100),pred$fit,col='black',lwd=3);
    dev.off()
    
  }
}
