library(dplyr)
library(tidyr)
library(stringr)
library(kernlab)
library(igraph)
library(Rcpp)
library(inline)
library(mixKernel)
#library(mixOmics)

cppFunction('double xfun2(const NumericMatrix x, const NumericMatrix y, double sigma) {
         int n = x.nrow(), m = y.nrow();
         double temp = 0;
         for (int i = 0; i < n; i++) {
           for (int j = 0; j < m ; j++) {
             temp += exp(-pow(pow(x(i,0)-y(j,0),2)+pow(x(i,1)-y(j,1),2),0.5)/8/sigma)
                      -exp(-pow(pow(x(i,0)-y(j,1),2)+pow(x(i,1)-y(j,0),2),0.5)/8/sigma);
           }
         }
         return temp/8/sigma/3.14159265;
      }')
start.time0 <- Sys.time()
setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/thickness_PH_normed_max/PH0_thinTOthick/")
#setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/gene_PH_normed_max/PH0_test/")
temp = list.files(pattern="*.csv")
PH0_thickTOthin = lapply(temp, read.csv)
Kmatrix <- diag(length(PH0_thickTOthin))
start.time <- Sys.time()
call.time = 0
for (i in 1:(length(PH0_thickTOthin)-1)){
  for (j in (i):length(PH0_thickTOthin)){
    x <- data.matrix(PH0_thickTOthin[[i]][,1:2])
    y <- data.matrix(PH0_thickTOthin[[j]][,1:2])
    x <- x[is.finite(rowSums(x)),]
    y <- y[is.finite(rowSums(y)),]
    call.time <- call.time + system.time(Kmatrix[i,j] <- xfun2(x, y, 1))[3]
  }
}
end.time <- Sys.time()
Kmatrix[lower.tri(Kmatrix)]  <- t(Kmatrix)[lower.tri(Kmatrix)]
end.time - start.time
call.time/3600
is.positive.semi.definite(data.matrix(test$mat), tol=1e-8)
write.table(Kmatrix, file='/N/u/yan30/Karst/R/thickPH0thinKernel.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)

setwd("/N/u/yan30/Karst/R/")
block1 = read.csv("thickPH0thickKernel_new.tsv", header=FALSE, sep='\t');
block2 = read.csv("thickPH1thickKernel_new.tsv", header=FALSE, sep='\t');
block3 = read.csv("thickPH0thinKernel_new.tsv", header=FALSE, sep='\t');
block4 = read.csv("thickPH1thinKernel_new.tsv", header=FALSE, sep='\t');
#block4 = read.csv("thickPH1thin2thickKernel.tsv", header=FALSE, sep='\t');
#plot(K)
# 
# f <- function(x,y) sqrt(sum((x-y)^2))
# Kentry <- outer( 
#   1:nrow(PH0_thickTOthin[[1]]), 1:nrow(PH0_thickTOthin[[2]]), 
#   Vectorize( function(i,j) f(PH0_thickTOthin[[1]][i,1:2], PH0_thickTOthin[[2]][j,1:2]) )
# )

#read the subject social demographic features
Kthick1 <- compute.kernel(data.matrix(block1), kernel.func = "kidentity");
Kthick2 <- compute.kernel(data.matrix(block2), kernel.func = "kidentity");
Kthick3 <- compute.kernel(data.matrix(block3), kernel.func = "kidentity");
Kthick4 <- compute.kernel(data.matrix(block4), kernel.func = "kidentity");
setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/")
block0 = read.csv("subFeaturesAll.csv", header=TRUE, sep=",");

#delete emtpy subjects, 159-4=155 remaining subs
noLabel = which(block0$DX == "")
block0 = block0[-noLabel, ]
block0$DX = factor(block0$DX)

feature <- block0$AGE
dummy <- mean(feature, na.rm=TRUE)
feature[is.na(feature)] <- dummy
Kage <- compute.kernel(data.frame(feature), kernel.func="gaussian.radial.basis", sigma = 100)
feature <- block0$EDU
dummy <- mean(feature, na.rm=TRUE)
feature[is.na(feature)] <- dummy
Kedu <- compute.kernel(data.frame(feature), kernel.func="gaussian.radial.basis", sigma = 100)
feature <- block0$GENDER
#feature[is.na(feature)] <- "E"
feature[feature==""] <- "E"
feature <- as.numeric(feature)
dummy <- mean(feature, na.rm=TRUE)
Kgender <- compute.kernel(data.frame(feature), kernel.func="gaussian.radial.basis", sigma = 100)

#kernel comparison and merging
pdf('/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/clusteringResults/kernelSimilarityAll_new.pdf')
cim.kernel(Thick1=Kthick1, Thick2=Kthick2, Thick3=Kthick3, Thick4=Kthick4, Age=Kage, Edu=Kedu, Gender=Kgender)
dev.off()
meta.kernel0 <- combine.kernels(Thick1=Kthick1, Thick2=Kthick2, Thick3=Kthick3, Thick4=Kthick4, method = "STATIS-UMKL")
meta.kernel1 <- combine.kernels(Age=Kage, Edu=Kedu, Gender=Kgender,method = "STATIS-UMKL")
meta.kernel <- combine.kernels(Thick1=Kthick1, Thick2=Kthick2, Thick3=Kthick3, Thick4=Kthick4, Age=Kage, Edu=Kedu, Gender=Kgender,method = "sparse-UMKL")
write.table(meta.kernel$kernel, file='/N/u/yan30/Karst/R/metaKernel.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)

library(data.table)
setwd("/gpfs/projects/UITS/IUNI/IMAGENE/AliceCSVs_12122018/")
block2 = read.csv("ThicknessValues.csv", header=TRUE, sep=",");
block3 = read.csv("ThicknessValues_AgeResidualized.csv", header=TRUE, sep=",");
block2 <- transpose(block2)
block3 <- transpose(block3)

#mixkernel return kernels with NaN, use kernelab instead
dt <- as.matrix(block2)
## initialize kernel function
rbf <- rbfdot(sigma = 0.0001)
rbf
## calculate kernel matrix
Korig = kernelMatrix(rbf, dt)
dt <- as.matrix(block3)
## initialize kernel function
rbf <- rbfdot(sigma = 0.0001)
rbf
## calculate kernel matrix
temp = kernelMatrix(rbf, dt)


#block2 <- block2[, colSums(abs(block2))!= 0, drop = F]
Korig <- compute.kernel(data.matrix(Korig), kernel.func = "kidentity");
Kregress <- compute.kernel((temp+t(temp))/2, kernel.func = "kidentity");

meta.kernel0 <- combine.kernels(Thick1=Korig, Age=Kage, Edu=Kedu, Gender=Kgender,method = "sparse-UMKL")
meta.kernel1 <- combine.kernels(Thick1=Kregress, Age=Kage, Edu=Kedu, Gender=Kgender,method = "full-UMKL")

K <- as.kernelMatrix(as.matrix(meta.kernel0$kernel));
Ks <- as.kernelMatrix(as.matrix(meta.kernel1$kernel));

#kernel.pca.result <- kernel.pca(Kthick,  ncomp = 4) #problematic, use kernel lab instead
##visalizations 
# G <- graph_from_adjacency_matrix(Kmatrix, weighted=TRUE, mode="undirected");
# plot(G)
# tkid <- tkplot(G) #tkid is the id of the tkplot that will open
# l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot
# tk_close(tkid, window.close = T)
# plot(G, layout=l)

#kernel for kernel lab
K <- as.kernelMatrix(as.matrix(meta.kernel$kernel));
K0 <- as.kernelMatrix(as.matrix(meta.kernel0$kernel));
K1 <- as.kernelMatrix(as.matrix(meta.kernel1$kernel));

Ks <- as.kernelMatrix(as.matrix(meta.kernel$kernel));
K0s <- as.kernelMatrix(as.matrix(meta.kernel0$kernel));
K1s <- as.kernelMatrix(as.matrix(meta.kernel1$kernel));

sc <- kkmeans(K, centers = 3);
sc
centers(sc)
size(sc)
withinss(sc)
plot(K, col=sc)

sc(K, centers = 3);
sc
centers(sc)
size(sc)
withinss(sc)
plot(K, col=sc)

pca_K <- kpca(K, features = 3, th = 1e-4)
pcv(pca_K)

##t-sne plots
library(tsne)
ecb = function(x,y){ plot(x,t='n'); text(x,col=sc)}
pdf('/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/clusteringResults/test.pdf')
tsne_K = tsne(K, epoch_callback = ecb, perplexity=5)
dev.off()

output = kernel.pca.result$x
output <- cbind(output, feature)

write.table(output, file='/N/u/yan30/Karst/R/thickPH1(2)KPCA.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)
write.table(sc@.Data, file='/N/u/yan30/Karst/R/thickPH1(2)label.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)
write.table(temp, file='/N/u/yan30/Karst/R/Covarates.tsv', quote=FALSE, sep='\t', col.names = T, row.names = F)


# colors = rainbow(length(unique(iris$Species)))
# names(colors) = unique(iris$Species)
# ecb = function(x,y){ plot(x,t='n'); text(x,labels=iris$Species, col=colors[iris$Species]) }
# tsne_iris = tsne(iris[,1:4], epoch_callback = ecb, perplexity=50)
# 
# pca_iris = princomp(iris[,1:4])$scores[,1:2]
# plot(pca_iris, t='n')
# text(pca_iris, labels=iris$Species,col=colors[iris$Species])
# 
# data(iris)
# sc <- kkmeans(as.matrix(iris[,-5]), centers=3)
# sc
# centers(sc)
# size(sc)
# withinss(sc)
