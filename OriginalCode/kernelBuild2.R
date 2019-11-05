library(dplyr)
library(tidyr)
library(stringr)
library(kernlab)
library(igraph)
library(Rcpp)
library(inline)
#library(mixKernel)
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

setwd("/gpfs/projects/UITS/IUNI/IMAGENE/Phom_output/thickness_PH_zscores/PH0_thinTOthick/")
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
write.table(Kmatrix, file='/N/u/yan30/Karst/R/thickPH0thinKernel_new.tsv', quote=FALSE, sep='\t', col.names = F, row.names = F)
