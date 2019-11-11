setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output")
library(dplyr)
library(tidyr)
library(stringr)
library(kernlab)
library(pROC)

# test data
xtest <- rbind(matrix(rnorm(120),,2),matrix(rnorm(120,mean=3),,2))
# test kernel matrix i.e. inner/kernel product of test data with
# Support Vectors
Ktest <- as.kernelMatrix(crossprod(t(xtest),t(x[SVindex(svp2), ])))
test <- predict(svp2, Ktest)
plot(test,data=xtest)     

#exclude empty labels
noLabel = which(block0$DX == "")
block0 = block0[-noLabel, ]
block0$DX = factor(block0$DX)

Ks = Ks[-noLabel, -noLabel]
K = K[-noLabel, -noLabel]


#LOOCV code
test <- as.vector(block0$DX)
testS <- as.vector(block0$DX)
testP <- as.vector(block0$DX)
testSP <- as.vector(block0$DX)
loocv_tmp <- matrix(NA, nrow = nrow(block0), 1)
for (k in 1:nrow(block0)) {
  #k = 1
  train <- ksvm(K[-k, -k], block0$DX[-k], type="C-svc", kernel='matrix', prob.model = TRUE)
  Ktest <- as.kernelMatrix(K[k, -k, drop=F][,SVindex(train), drop=F])
  testP[k] <- predict(train, Ktest, type = "probabilities")
  test[k] <- predict(train, Ktest)
  loocv_tmp[k, ] <- as.numeric(as.character(test[k]))-1==(as.numeric(block0$DX[k])+1)%%2
}
loocv <- colMeans(loocv_tmp)

loocvS_tmp <- matrix(NA, nrow = nrow(block0), 1)
for (k in 1:nrow(block0)) {
  #k = 1
  train <- ksvm(Ks[-k, -k], block0$DX[-k], type="C-svc", kernel='matrix', prob.model = TRUE)
  Ktest <- as.kernelMatrix(Ks[k, -k, drop=F][,SVindex(train), drop=F])
  testSP[k] <- predict(train, Ktest, type = "probabilities")
  testS[k] <- predict(train, Ktest)  
  loocvS_tmp[k, ] <- as.numeric(as.character(testS[k]))-1==(as.numeric(block0$DX[k])+1)%%2
}
loocvS <- colMeans(loocvS_tmp)

roc1 = roc((as.numeric(block0$DX))%%2, as.numeric(testP))
pdf('/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/clusteringResults/RAWoriginal.pdf')
plot.roc(roc1,print.auc = TRUE,print.auc.cex=2,
         auc.polygon = TRUE,
         grid=c(0.1, 0.2),
         grid.col = c("green", "red"),
         max.auc.polygon = TRUE,
         auc.polygon.col = "skyblue",
         #print.thres = TRUE, 
         cex.lab=1.6, cex.axis=1.6, cex.main=1.6, cex.sub=1.6)
dev.off()

roc1 = roc((as.numeric(block0$DX))%%2, as.numeric(testSP))
pdf('/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/clusteringResults/RAWregressed.pdf')
plot.roc(roc1,print.auc = TRUE,print.auc.cex=2,
         auc.polygon = TRUE,
         grid=c(0.1, 0.2),
         grid.col = c("green", "red"),
         max.auc.polygon = TRUE,
         auc.polygon.col = "skyblue",
         #print.thres = TRUE, 
         cex.lab=1.6, cex.axis=1.6, cex.main=1.6, cex.sub=1.6)
dev.off()


library(caret)
# define training control
train_control <- trainControl(method="LOOCV",savePredictions = T)
# set custom models
lpSVM <- list(type = "Classification", library = "kernlab", loop = NULL)
# model parameters
prm <- data.frame(parameter = c("C", "sigma"),
                  class = rep("numeric", 2),
                  label = c("Cost", "Sigma"))
lpSVM$parameters <- prm
# parameter tuning grid
svmGrid <- function(x, y, len = NULL, search = "grid") {
  library(kernlab)
  ## This produces low, middle and high values for sigma 
  ## (i.e. a vector with 3 elements). 
  sigmas <- sigest(as.matrix(x), na.action = na.omit, scaled = TRUE)
  expand.grid(sigma = mean(sigmas[-2]),
              C = 2 ^((1:len) - 3))
}
lpSVM$grid <- svmGrid
# the fit element
svmFit <- function(x, y, wts, lev, last, weights, classProbs, ...) {
  ksvm(as.kernelMatrix(x), y = y,
       prob.model = classProbs,
       ...)
}
lpSVM$fit <- svmFit
# the predict element
svmPred <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
lpSVM$predict <- svmPred

#The prob Element
svmProb <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type="probabilities")
lpSVM$prob <- svmProb


svmSort <- function(x) x[order(x$C), ]
lpSVM$sort <- svmSort
lpSVM$levels <- function(x) lev(x)


# train the model
#set.seed(998)
#inTraining <- createDataPartition(iris$Species, p = .75, list = FALSE)
#training <- iris[ inTraining,]
#testing  <- iris[-inTraining,]
block1$tag = block0$DX
set.seed(825)
Laplacian <- caret::train(x = Ks@.Data, y = block0$DX, 
                   method = lpSVM,`
                   preProc = c("center", "scale"),
                   tuneLength = 8,
                   trControl = train_control)
print(Laplacian, digits = 3)
ggplot(Laplacian)
# summarize results
library(plyr)
library(pROC)
plot(roc(predictor = as.numeric(Laplacian$pred$pred), response = Laplacian$pred$obs))