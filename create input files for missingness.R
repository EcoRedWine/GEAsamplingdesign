rm(list=ls())

setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/ALLsims/")
dir()

data <- read.csv(dir()[1])    # read in the first data file
data

### Sample 100 individuals 

samp <- read.csv("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/sample100.csv", header=F)  # 500 or 100 sampled individuals
samp <- samp[,1]                              # convert to a vector
sampdat <- data[samp,]                        # subsample the full dataset

# Reduce the SNP data set to only one column per locus (from two)
snps <- sampdat[,8:207] 
loci <- seq(from=1, to=200, by=2)
newsnps <- snps[,loci]
newsnps

snps <- newsnps; colnames(snps) <- seq(0, 99, by=1)
View(snps)

snp.mat <- as.matrix(snps)
head(snp.mat)
dim(snp.mat)
class(snp.mat)

### Missing random 
n <- c(0.05, 0.1, 0.2, 0.5) ## percentages of missing data
missing.loci <- list()

for (i in 1:length(n)) {
  snp.mat2 <- snp.mat
  snp.mat2[,1][sample(1:length(snp.mat2[,1]), round(n[i]*length(snp.mat2[,1])), replace=FALSE)] <- 'NA'
  snp.mat2[,2:100][sample(1:length(snp.mat2[,2:100]), round(n[i]*length(snp.mat2[,2:100])), replace=FALSE)] <- 'NA'
  missing.loci[[i]] <- which(snp.mat2 == 'NA')
  setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/")
  write.table(missing.loci[[i]], file = paste0("random_missing_loci_",n[i],"_100.csv"),row.names=FALSE, col.names=FALSE)
}


# This generates 4 input files: a list of indexes for 5, 10, 20 and 50% random missing data to subset full datasets 

length(missing.loci[[1]]) # 2500  (5% * 50,000)  OR  500 (5% * 10,000)
length(missing.loci[[2]]) # 5000  (10% * 50,000) OR  1000 (10% * 10,000)
length(missing.loci[[3]]) # 10000 (20% * 50,000) OR  2000 (20% * 10,000)
length(missing.loci[[4]]) # 25000 (50% * 50,000) OR  5000 (50% * 10,000)

# Check proportions of neutral and adaptive loci that are missing on one of the actual datasets

missing10per <- read.csv("random_missing_loci_0.1_100.csv",header=FALSE)

missing10per <- missing10per[,1]
missing10per
length(missing10per)

snps.missing.10 <- as.matrix(snps)
snps.missing.10[missing10per] <- 'NA'

length(which(snps.missing.10[,1] == 'NA'))
length(which(snps.missing.10[,2:100] == 'NA'))
length(which(snps.missing.10 == 'NA'))

### Missing adaptive (locus0) 

dim(snps)  # 100 x 500
dim(snp.mat)

n <- c(0.05, 0.1, 0.2, 0.5) ## percentages of missing adaptive locus
missing.adapt <- list()

for (i in 1:length(n)) {
  snp.mat3 <- snp.mat
  snp.mat3[sample(1:length(snp.mat3[,1]), round(n[i]*length(snp.mat3[,1])), replace=FALSE)] <- 'NA'
  missing.adapt[[i]] <- which(snp.mat3 == 'NA')
  setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/")
  write.table(missing.adapt[[i]], file = paste0("adaptive_missing_loci_",n[i],"_100.csv"),row.names=FALSE, col.names=FALSE)
}

## This generates 4 input files: a list of indexes for 5, 10, 20 and 50% of missing data in the adaptive locus (locus 0) only

length(snp.mat3[,1])       # 500 individuals (rows) OR 100 individuals
length(missing.adapt[[1]]) # 25 (5% * 500)    OR  5 (5% * 100)
length(missing.adapt[[2]]) # 50 (10% * 500)   OR  10 (10% * 100)
length(missing.adapt[[3]]) # 100 (20% * 500)  OR  20 (20% * 100)
length(missing.adapt[[4]]) # 250 (50% * 500)  OR  50 (50% * 100)

which(snp.mat3[,2:100] == 'NA')  # No NAs in columns 2:99


### Missing neutral (locus1 - locus99) 

n <- c(0.05, 0.1, 0.2, 0.5) ## percentages of missing adaptive locus
missing.neutral <- list()

for (i in 1:length(n)) {
  snp.mat4 <- snp.mat
  snp.mat.neut <- snp.mat[,2:100]
  snp.mat.neut[sample(1:length(snp.mat.neut), round(n[i]*length(snp.mat.neut)), replace=FALSE)] <- 'NA'
  snp.mat4[,2:100] <- snp.mat.neut
  missing.neutral[[i]] <- which(snp.mat4 == 'NA')
  setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/")
  write.table(missing.neutral[[i]], file = paste0("neutral_missing_loci_",n[i],"_100.csv"),row.names=FALSE, col.names=FALSE)
}

## This generates 4 input files: a list of indexes for 5, 10, 20 and 50% of missing data in the neutral loci (locus1 - locus99) only

length(snp.mat4[,2:100])       # 49500 neutral snps   OR 9900 neutral snps
length(missing.neutral[[1]]) # 2475 (5% * 49500)      OR 495 (5% * 9900)
length(missing.neutral[[2]]) # 4950 (10% * 49500)     OR 990 (10% * 9900)
length(missing.neutral[[3]]) # 9900 (20% * 49500)     OR 1980 (20% * 9900)
length(missing.neutral[[4]]) # 24740 (50% * 49500)    OR 4950 (50% * 9900)

which(snp.mat4[,1] == 'NA')  # No NAs in column 1


