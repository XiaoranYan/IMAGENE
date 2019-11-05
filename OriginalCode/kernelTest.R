setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output")
library(dplyr)
library(tidyr)
library(stringr)
library(kernlab)
library(igraph)
library(Rcpp)
library(inline)
library(mixKernel)

#For gene distance matrices based on betti curves
setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/")
block1 = read.csv("Betticurves_differences_genes_PH0.csv", header=FALSE, sep=',');
block1 <- block1[-c(1),]
block1 <- block1[,-c(1)]
#simple RBF kernel
blockK <- exp(data.matrix(block1)*(-0.001))


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
#setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/thickness_PH_normed_max/PH0_thinTOthick/")
setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/gene_PH_normed_max/PH0_test/")
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
write.table(Kmatrix, file='/N/u/yan30/Karst/R/genePH0thinKernel.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)


Kthick1 <- compute.kernel(data.matrix(blockK), kernel.func = "kidentity");
pdf('/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/clusteringResults/kernelSimilarity_betti.pdf')
cim.kernel(Thick1=Kthick1, Thick2=Kthick2)
dev.off()

meta.kernel0 <- combine.kernels(Thick1=Kthick1, Thick2=Kthick2, method = "STATIS-UMKL")

#kernel for kernel lab
K <- as.kernelMatrix(as.matrix(Kthick1$kernel));

sc <- kkmeans(K, centers = 3);
sc
centers(sc)
size(sc)
withinss(sc)
plot(K, col=sc)

sc <- specc(K, centers = 3);
sc
centers(sc)
size(sc)
withinss(sc)
plot(K, col=sc)

pca_K <- kpca(K, features = 6, th = 1e-4)
pcv(pca_K)
rotated(pca_K)

##t-sne plots
library(tsne)
ecb = function(x,y){ plot(x,t='n'); text(x,col=sc)}
pdf('/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/clusteringResults/test.pdf')
tsne_K = tsne(K, epoch_callback = ecb, perplexity=5)
dev.off()

output = kernel.pca.result$x
output <- cbind(output, feature)

write.table(pcv(pca_K), file='/N/u/yan30/Karst/R/genePH0KPCA.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)
write.table(sc@.Data, file='/N/u/yan30/Karst/R/genePH0label.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)
write.table(temp, file='/N/u/yan30/Karst/R/Covarates.tsv', quote=FALSE, sep='\t', col.names = T, row.names = F)


