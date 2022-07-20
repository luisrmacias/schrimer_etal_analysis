##proof of principle analysis of the study by Schirmer et al.

library(data.table)
library(openxlsx)
library(mixOmics)
set.seed(3498)
setwd("bioinfo/GUIDE-IBD/")
data <- read.xlsx("pediatric_UC-mmc7.xlsx", startRow = 2)
d.data <- setDT(data)

ind.out52w <- dcast(d.data, SubjectID ~ ., value.var = "CSFREE_REMISSION_WK52", fun.aggregate = function(x) any(x=="Remission"))
names(ind.out52w)[2] <- "Remission"
table(ind.out52w$Remission, exclude = NULL)
initial <- data[data$collectionWeek==0,]

##leave out variables different from microbiome
X <- initial[, -(1:24)]
##leave out OTUs with all 0 values
X <- subset(X, select = colSums(X)>0)

explainedVariance <- tune.pca(X, ncomp = 20, center = TRUE, scale = FALSE) 
plot(explainedVariance)


res.pca <- pca(X, ncomp = 6, center = TRUE, scale = TRUE)
plotIndiv(res.pca, group = initial$INITIAL_TRT_C4, legend = TRUE)
plotIndiv(res.pca, group = initial$sampleType, legend = TRUE)
plotIndiv(res.pca, group = initial$sampleType, legend = TRUE, comp = c(2,3))
plotIndiv(res.pca, group = initial$CSFREE_REMISSION_WK52, legend = TRUE)

## plotVar(res.pca)

## Sparse partial least squares discriminant analysis tuning

Y <- initial$CSFREE_REMISSION_WK52
selection <-complete.cases(Y)
uc.splsda <- splsda(X[selection,], Y[selection], ncomp = 10)

plotIndiv(uc.splsda , comp = 1:2, 
          group = initial$CSFREE_REMISSION_WK52[selection], ind.names = initial$SubjectID[selection],  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')


# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(uc.splsda, comp.predicted=2, dist = "max.dist")

plotIndiv(uc.splsda, comp = 1:2,
          group = initial$CSFREE_REMISSION_WK52[selection], ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")


##number of component selection for splsda by classification error
perf.uc.srbct <- perf(uc.splsda, validation = "Mfold", 
                      folds = 10, nrepeat = 50, # use repeated cross-validation
                      auc = TRUE, progressBar = TRUE) # include AUC values

plot(perf.uc.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.uc.srbct$choice.ncomp
##first component seem to be best outcome classifier (with no other component)

##variable number selection for splsda
list.keepX <- c(1:10,  seq(20, 1000, 10))

# undergo the tuning process to determine the optimal number of variables
tune.splsda.uc <- tune.splsda(X[selection,], Y[selection], ncomp = 3, # calculate for first 2 components
                              validation = 'Mfold',
                              folds = 10, nrepeat = 50, # use repeated cross-validation
                              dist = 'max.dist', # use max.dist measure
                              measure = "BER", # use balanced error rate of dist measure
                              test.keepX = list.keepX,
                              cpus = 2, progressBar = TRUE) # allow for parallelization to decrease runtime

plot(tune.splsda.uc, col = color.jet(3))

tune.splsda.uc$choice.ncomp$ncomp
## Result is 1 but splsda can´t be performed on 1 component only

tune.splsda.uc$choice.keepX 
## what is the optimal amount of variables according to tune.splsda()
## comp1: 520, comp2: 220, comp3: 3

##optimal.ncomp <- tune.splsda.uc$choice.ncomp$ncomp
##normally we would input number of components from tune procedure but in this case we set to 2
optimal.ncomp <- 2
optimal.keepX <- tune.splsda.uc$choice.keepX[1:optimal.ncomp]


final.splsda <- splsda(X[selection,], Y[selection], 
                       ncomp = 2, 
                       keepX = optimal.keepX)

plotIndiv(final.splsda, comp = c(1,2), # plot samples from final model
          group = initial$CSFREE_REMISSION_WK52[selection], ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = ' (a) sPLS-DA on SRBCT, comp 1 & 2')

vip <- vip(final.splsda)
hist(vip[,1])

highest_vip <- sort(vip[vip[,1]>3,1], decreasing = TRUE)