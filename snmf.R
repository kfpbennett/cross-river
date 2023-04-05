library(LEA)
library(tidyverse)
library(data.table)

setwd('your/working/directory')

# ========================== Admixture =========================================

# Read in the data
vcf2geno('filtered4_thin5k.recode.vcf', 'filtered4_thin5k.geno')
geno2lfmm('filtered4_thin5k.geno', 'f4t5k.lfmm')

vcf2geno('f4t5k_males.recode.vcf', 'f4t5k_males.geno')
geno2lfmm('f4t5k_males.geno', 'f4t5k_males.lfmm')

# Males and females ------------------------------------------------------------
set.seed(323)
proj45k <- snmf('f4t5k.lfmm', K = 1:8, entropy = TRUE, project = 'new')

# K = 2 is best
plot(proj45k)

qs <- as.data.frame(Q(proj45k, K = 2, run = 1))
qs <- qs[c(102:132, 64:82, 133:151, 83:101, 152:173, 1:63),]

mycolors <- c('salmon','lavender')
barplot(t(-qs[,1:2]), border = NA, space = 0, col = mycolors)

# Males only -------------------------------------------------------------------
set.seed(670)
proj45km <- snmf('f4t5k_males.lfmm', K = 1:8, entropy = TRUE, project = 'new')

plot(proj45km)

mqs <- as.data.frame(Q(proj45km, K = 2, run = 1))
mqs <- mqs[c(82:104, 48:66, 105:120, 67:81, 121:135, 1:47),]
barplot(t(mqs[,1:2]), border = NA, space = 0, col = mycolors)
