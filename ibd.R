setwd('your/working/directory')

library(tidyverse)
library(data.table)
library(hierfstat)
library(adegenet)

# =============================== IBD ANALYSIS =================================
# Read in data for IBD analysis
pop.gps <- fread('../pops_gps.txt') %>% arrange(pop)
samps <- fread('sampleinfo_filtered2.txt')
dat.bed <- read.VCF('filtered4_thin5k.recode.vcf', convert.chr = FALSE)

# Males and females ------------------------------------------------------------
dat.df <- as.matrix(dat.bed) %>% as.data.frame()
rownames(dat.df) <- NULL

pops <- unlist(strsplit(dat.bed@ped$famid,'_'))[seq(1,346,2)]
pops[pops %in% c('WN5', 'WN6')] <- 'WN'

dat <- mutate(dat.df, pop = as.factor(pops))
dat <- dat[,c(12376,1:12375)]

popdist <- as.matrix(dist(cbind(pop.gps$lat, pop.gps$lon)))

uq.pops <- unique(pops)

# Calculate pairwise FST between each population. This runs slowly
pwf <- pairwise.WCfst(dat)

pwf.df <- data.frame(
  pop1 = c(rep('EL', 7), rep('EM', 6), rep('EU', 5), rep('WL', 4), rep('WM', 3), 
           'WN', 'WN', 'WS'),
  pop2 = c(uq.pops[2:8], uq.pops[3:8], uq.pops[4:8], uq.pops[5:8],
           uq.pops[6:8], uq.pops[7:8], uq.pops[8]),
  fst = c(pwf[2:8,1], pwf[3:8,2], pwf[4:8,3], pwf[5:8,4], pwf[6:8,5], 
          pwf[7:8,6], pwf[8,7]),
  dist = c(popdist[2:8,1], popdist[3:8,2], popdist[4:8,3], popdist[5:8,4], 
           popdist[6:8,5], popdist[7:8,6], popdist[8,7]),
  bank = c(rep('same', 2), rep('opposite', 5), 'same', rep('opposite', 10),
           rep('same', 10))
)

dist.lm <- lm(fst ~ dist, pwf.df)

pwf.df <- pwf.df %>%
  mutate(residual = unname(residuals(dist.lm)))

# Run the first time, load if coming back to this to avoid  needing to 
# rerun pairwise FST:
# write.csv(pwf.df, 'pwf.df.csv', row.names = FALSE)

# Mantel test ------------------------------------------------------------------
mant <- mantel.rtest(as.dist(pwf), as.dist(popdist), nrepet = 1000)

# Randomization test -----------------------------------------------------------
set.seed(888)
diff <- c()
for(i in 1:10000){
  same <- sample(1:28, 13, replace = FALSE)
  opposite <- c(1:28)[which(!c(1:28) %in% same)]
  nuld.s <- mean(lm(fst ~ dist, pwf.df)$residuals[same])
  nuld.o <- mean(lm(fst ~ dist, pwf.df)$residuals[opposite])
  diff[i] <- nuld.o - nuld.s
}

# Difference in residuals
mean(dist.lm$residuals[pwf.df$bank == 'opposite']) - 
  mean(dist.lm$residuals[pwf.df$bank == 'same'])
# Confidence interval
quantile(diff, c(0.025,0.975))
quantile(diff, 0.998075)
# P-value
1 - 0.998075

# Males only -----------------------------------------------------------------
males <- samps[samps$Sex == 'M',]$Sample

datm.bed <- dat.bed[dat.bed@ped$famid %in% males]
datm.df <- as.matrix(datm.bed) %>% as.data.frame()
rownames(datm.df) <- NULL

mpops <- unlist(strsplit(datm.bed@ped$famid,'_'))[seq(1,270,2)]
mpops[mpops %in% c('WN5', 'WN6')] <- 'WN'

mdat <- mutate(datm.df, pop = as.factor(mpops))
mdat <- mdat[,c(12376,1:12375)]

# Calculate pairwise FST between each population, runs slowly
mpwf <- pairwise.WCfst(mdat)

mpwf.df <- data.frame(
  pop1 = c(rep('EL', 7), rep('EM', 6), rep('EU', 5), rep('WL', 4), rep('WM', 3), 
           'WN', 'WN', 'WS'),
  pop2 = c(uq.pops[2:8], uq.pops[3:8], uq.pops[4:8], uq.pops[5:8],
           uq.pops[6:8], uq.pops[7:8], uq.pops[8]),
  fst = c(mpwf[2:8,1], mpwf[3:8,2], mpwf[4:8,3], mpwf[5:8,4], mpwf[6:8,5], 
          mpwf[7:8,6], mpwf[8,7]),
  dist = c(popdist[2:8,1], popdist[3:8,2], popdist[4:8,3], popdist[5:8,4], 
           popdist[6:8,5], popdist[7:8,6], popdist[8,7]),
  bank = c(rep('same', 2), rep('opposite', 5), 'same', rep('opposite', 10),
           rep('same', 10))
)

mdist.lm <- lm(fst ~ dist, mpwf.df)

mpwf.df <- mpwf.df %>%
  mutate(residual = unname(residuals(mdist.lm)))

# Run the first time, load if coming back to this to avoid  needing to 
# rerun pairwise FST:
# write.csv(mpwf.df, 'mpwf.df.csv', row.names = FALSE)

# Mantel test ------------------------------------------------------------------
mmant <- mantel.rtest(as.dist(mpwf), as.dist(popdist), nrepet = 1000)

# Randomization test -----------------------------------------------------------
set.seed(988)
diff <- c()
for(i in 1:10000){
  same <- sample(1:28, 13, replace = FALSE)
  opposite <- c(1:28)[which(!c(1:28) %in% same)]
  nuld.s <- mean(lm(fst ~ dist, mpwf.df)$residuals[same])
  nuld.o <- mean(lm(fst ~ dist, mpwf.df)$residuals[opposite])
  diff[i] <- nuld.o - nuld.s
}

# Difference in residuals
mean(mdist.lm$residuals[mpwf.df$bank == 'opposite']) - 
  mean(mdist.lm$residuals[mpwf.df$bank == 'same'])
# Confidence interval
quantile(diff, c(0.025,0.975))
quantile(diff, 0.999806)
# P-value
1 - 0.999806
