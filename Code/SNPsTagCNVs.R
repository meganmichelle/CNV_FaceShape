#####
# Calculate the LD between SNPs and 
# BL2 CNV windows
# Chromosome  20
#####


## Kinship Matrix
kin <- read.table('Chr20kinship.kindump')
kin <- as.matrix(kin)

## Format the matrix for the calculation
vkin <- chol2inv(chol(kin))
rm(kin)

## Load in the  CNV data:
cnv <- read.table('Chr20_CNV.txt', 
                 header = TRUE)

## Load  inthe SNP data
snp <- read.table('Chr20_SNP.txt', header  = TRUE)

### Info needed  for the calculation
n <- nrow(snp) # Number of SNPs
n1 <- matrix(1, nrow = n, ncol = 1) # column with n 1's
n1_mat <- matrix(1, nrow = n, ncol = n) # n by n matrix of 1's

### calculation used for each CNV/SNP pair:
mid <- (n1_mat %*% vkin) / as.numeric(t(n1) %*% vkin %*% n1) 

###  Function to calculate LD
structured_corr <- function(xm,xl,mid){
  snp <- cbind(xl,xm)
  out <- snp - (mid %*% snp)
  lm <- t(out) %*% vkin %*% out
  as.numeric(lm[1,2])^2 / (as.numeric(lm[1,1]) * 
                             as.numeric(lm[2,2])) 
}

### xl will be  the column of SNP information (0/1/2)
### xm will be  the column of CNV information (0/1 for Same  direction model)

#### Note that this calculation is tedious
#### Recommend to  only  calculate LD between CNVs and SNPs
###### within a window (1 Mb was used in the CNV paper)

data  <- as.data.frame(NA, nrow = ncol(snps), ncol = ncol(cnv))

Sys.time()
for(j in c(1:ncol(cnv))){
    for( k in 1:ncol(snps)){
      temp[k,1] <- structured_corr(snps[,k], df[,j], mid = mid)
    }
    data[,j] <- temp
    print(j)
}
Sys.time()




