setwd('your/working/directory')

library(adegenet)
library(poppr)
library(tidyverse)

# ============================== AMOVA =========================================
# Read in data for AMOVA
dat <- read.PLINK('plink.raw', parallel = FALSE)
samps <- fread('sampleinfo_filtered2.txt')

# Males and females ------------------------------------------------------------

oldnames <- dat$ind.names
pops <- c(rep('EL', 18), rep('EM', 20), rep('EU', 25), rep('WL', 19),
          rep('WM', 19), rep('WN', 31), rep('WS', 19), rep('WU', 22))
newnames <- paste0(pops, '_', oldnames)
dat$ind.names <- newnames

dat$other$bank <- as.factor(c(rep('E', 63), rep('W', 110)))
dat$other$pop <- as.factor(pops)
strata(dat) <- data.frame(other(dat))
hier(dat) <- ~bank/pop

dat.teribe <- dat[dat$other$pop %in% c('WL', 'WM', 'WU', 'WN', 'WN', 'WS')]
dat.teribe$other$bank <- as.factor(c(rep('N', 19), rep('S', 19), 
                                     rep('N', 31), rep('S', 41)))
strata(dat.teribe) <- data.frame(other(dat.teribe))
hier(dat.teribe) <- ~bank/pop

set.seed(90)
# No sig bank effect
dat.amova <- poppr.amova(dat, ~bank/pop, within = FALSE, filter = TRUE)
dat.amova.test <- randtest(dat.amova, nrepet = 9999)

# Sig bank effect
dat.nowu.amova <- poppr.amova(dat[c(1:151),], ~bank/pop, 
                              within = FALSE, filter = TRUE)
dat.nowu.amova.test <- randtest(dat.nowu.amova, nrepet = 9999)

set.seed(64)
# Teribe bank effect test, no sig effect
dat.teribe.amova <- poppr.amova(dat.teribe, ~bank/pop, 
                                within = FALSE, filter = TRUE)
dat.teribe.amova.test <- randtest(dat.teribe.amova, nrepet = 9999)


# Males only -------------------------------------------------------------------
datm <- dat
males <- samps[samps$Sex == 'M',]$Sample
pops2 <- c(rep('EL', 18), rep('EM', 20), rep('EU', 25), rep('WL', 19),
          rep('WM', 19), rep('WN5', 10), rep('WN6', 21), rep('WS', 19), 
          rep('WU', 22))
newnames2 <- paste0(pops2, '_', oldnames)
datm$ind.names <- newnames2
datm <- datm[newnames2 %in% males,]

set.seed(19)
# Sig bank effect
datm.amova <- poppr.amova(datm, ~bank/pop, within = FALSE, filter = TRUE)
datm.amova.test <- randtest(datm.amova, nrepet = 9999)




