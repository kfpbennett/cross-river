source('your/working/directory/pca-source.R')

# ================================ PCA =========================================
# Read in data for PCA
dat <- read.PLINK('full_filtered4_thin5k/plink.raw', parallel = FALSE)
samps <- fread('sampleinfo_filtered2.txt')

# Males and females ------------------------------------------------------------

# PCA --------------------------------------------------------------------------
oldnames <- dat$ind.names
pops <- c(rep('EL', 18), rep('EM', 20), rep('EU', 25), rep('WL', 19),
          rep('WM', 19), rep('WN', 31), rep('WS', 19), rep('WU', 22))
newnames <- paste0(pops, '_', oldnames)
dat$ind.names <- newnames

pca <- glPcaFast(dat, nf = 50)

pca$eig[1]/sum(pca$eig)
pca$eig[2]/sum(pca$eig)

pca_cols <- c('#800404','#C30E0E','#FF6B6B',
                       '#077408','#15AD17',
                       '#1AC391', '#0EA87A', 
                       '#66ED58')
                       
plot(-1*pca$scores[,1],pca$scores[,2], 
     col = alpha(pca_cols[as.factor(pops)], 0.5), pch = 16)

# DAPC -------------------------------------------------------------------------
dat$pop <- as.factor(pops)
x <- 10:80
set.seed(807)
res <- xvalDapc.custom(tab(dat, NA.method = 'mean'), pop(dat), 
                n.pca = x, n.rep = 100)
resdf <- data.frame(x = x, y = res, xsq = x^2)
mod <- lm(y ~ x + xsq, data = resdf)
xvals <- seq(10, 80, 0.2)
pred <- predict(mod, list(x = xvals, xsq = xvals^2))
plot(x, res)
lines(xvals, pred)

da <- dapc(dat, grp = as.factor(pops), glPca = pca, n.pca = 50)

scatter(da, scree.da = FALSE)

# Males only -------------------------------------------------------------------

# PCA --------------------------------------------------------------------------
datm <- dat
males <- samps[samps$Sex == 'M',]$Sample
pops2 <- c(rep('EL', 18), rep('EM', 20), rep('EU', 25), rep('WL', 19),
           rep('WM', 19), rep('WN5', 10), rep('WN6', 21), rep('WS', 19), 
           rep('WU', 22))
newnames2 <- paste0(pops2, '_', oldnames)
datm$ind.names <- newnames2

datm <- datm[newnames2 %in% males,]

pcam <- glPcaFast(datm, nf = 50)

pcam$eig[1]/sum(pcam$eig)
pcam$eig[2]/sum(pcam$eig)

pops3 <- as.factor(substr(males, 1, 2))
                       
plot(-1*pcam$scores[,1],-1 * pcam$scores[,2], 
     col = alpha(pca_cols[pops3], 0.5), pch = 16,
     xlab = 'PC1 (2.33% variance explained)',
     ylab = 'PC2 (1.47% variance explained')

# DAPC -------------------------------------------------------------------------
dam <- dapc(datm, grp = pops3, glPca = pcam)
scatter(dam, scree.da = FALSE)

