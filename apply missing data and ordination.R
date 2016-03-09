rm(list=ls())

library(raster)
library(vegan)

setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/AllH9sims/")   # CHANGE THE WORKING DIRECTORY 
dir()
list.dir <- dir()
length(list.dir)

data <- lapply(list.dir, read.csv) # this takes some time 
length(data)
class(data)
dim(data[[5]])

coord <- data[[1]][,2:3]

L10.files <- dir("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Surfaces_Sample/", pattern = glob2rx("L10H9*_aa.asc"))   #CHANGE THE HABITAT CONFIGURATION
L10.files

qrule = list()
for (q in 1:10) {
  setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Surfaces_Sample/")
  qrule[[q]] <- raster(L10.files[q])
}

qrule.rep <- rep(qrule, 24)
plot(qrule.rep[[12]]) # just for checking 

habitat = list()
for (h in 1:length(qrule)) {
  habitat[[h]] <- extract(qrule[[h]],coord)
}

## now there is a list of 10 habitat files 

habitat.all <- rep(habitat,24)
habitat.all[[1]] == habitat.all[[2]]   ## These should not be the same 
habitat.all[[11]] == habitat.all[[1]]  ## These should be the same


## Now we have a list of 240 habitat values ordered by replicate 0-9 and repeated; list of datasets with coordinates, age, sex, and snps (207 columns - need to remove snps and combine to 1 allele per locus). 

head(coord)

xvar <- coord[,1]
yvar <- coord[,2]

fulldata <- lapply(seq_along(data), function(x) cbind(xvar, yvar, habitat.all[[x]], data[[x]][,8:207]))

length(fulldata)
fulldata[[1]] == fulldata[[2]]  # should be false 

## fulldata is a list of all datasets containing x and y coordinates, habitat value, and snps 

## sample to 100 individuals CHANGE THE SAMPLE FILE 
samp <- read.csv("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/sample500.csv", header=F)
samp <- samp[,1]  

sampdata <- llply(fulldata, function(x) subset(x[samp,]))
View(sampdata[[186]])
length(sampdata)
dim(sampdata[[1]])
head(sampdata[[1]])
sampdata[[1]] == sampdata[[2]]  # should be false

# extract the environmental data:
env <- llply(sampdata, function(x) subset(x[,1:3])) 
env <- lapply(env, function(x) scale(x,center=T, scale=T)) # scale and center env vars so they are comparable
env <- llply(env, function(x) data.matrix(x))
class(env[[1]])

for (i in 1:240) {
  colnames(env[[i]]) <- c("xvar","yvar","habitat")
}

# extract the genetic data: 

snps <- llply(sampdata, function(x) subset(x[,4:203]))
snps <- llply(snps, function(x) subset(x[,seq(1, 200, 2)]))
dim(snps[[1]])
snps[[1]] == snps[[2]]  # should be false
View(snps[[186]])

# convert to matrix
snps_mat <- llply(snps, function(x) as.matrix(x))
class(snps_mat[[1]])
length(snps_mat)

for (j in 1:240) {
  colnames(snps[[j]]) <- seq(0,99,by=1)
  colnames(snps_mat[[j]]) <- seq(0,99,by=1)
}

# generate the missing data. Here, for each matrix in the list of 240 datasets, I think we should create the NAs using the first input file for 5% missing data and run the RDA/dbRDA analysis, then go back to the top and re run for 10% and so on. 

## remove NA individuals
## change data to NAs 
## impute the missing data (mean)
## scale 
## vegdist (bray)


## TEST ##
# table_missing <- list()
# new_list <- list()

# for (i in 1:length(snps_mat)) {
#   allsnps <- snps_mat[[i]]
#   for (j in 1:length(missing_files)) {
#     missing <- read.csv(missing_files[j],header=FALSE)
#     missing <- missing[,1]  # change to a vector 
#     new_snp_mat <- allsnps
#     new_snp_mat[missing] <- NA
#     table_missing[[j]] <- table(new_snp_mat[,1], useNA = "always")
#   }
#   new_list[[i]] <- ldply(table_missing, data.frame)
# }

## REAL DATA ##

setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/Input Files /")

#Read in the missing data input files (4 at a time)	CHANGE INPUT FILE 
missing_files <- list.files("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/Input Files /", pattern = glob2rx("random_*_500.csv"), full.names = F, recursive = F) 
missing_files

for (i in 1:length(snps_mat)) {
  allsnps <- snps_mat[[i]]
  missing <- read.csv(missing_files[1],header=FALSE)  # just the first missing data input file 
  missing <- missing[,1]  # change to a vector 
  new_snp_mat <- allsnps
  new_snp_mat[missing] <- NA 
  snps_mat[[i]] <- new_snp_mat
}

num_missing <- ldply(snps_mat, function(x) table(x,useNA='always'))
num_missing

# remove NA individuals 

setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/AllH9sims/")
tempnames <- list()
snps_test <- snps_mat
View(snps_test[[186]])

for (j in 1:length(snps_test)) {
  subsetNA <- apply(snps_test[[j]],1,function(x) length(unique(x))>1) ## applied to rows (NAs)
  snps_test[[j]] <- snps_test[[j]][subsetNA,]
  env[[j]] <- env[[j]][subsetNA,]
  tempnames[[j]] <- subsetNA
}

# We know that dataset 186 in H9 sims has a missing individual 
tempnames[[186]][213]  # == FALSE
dim(env[[186]])
dim(snps_test[[186]])

# check final # individuals:
dir()
list.dir <- dir()
messy <- list()

for (l in 1:length(snps_test)) {
  dims <- dim(snps_test[[l]])
  messy[[l]] <- cbind(list.dir[l], dims[1], dims[2], names(which(tempnames[[l]]==FALSE)))
}

NA_inds <- ldply(messy, data.frame)
colnames(NA_inds) <- c("sim","#indiv","#loci","NA individual")
write.csv(NA_inds, "/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Test_Results/H9_individual_NAs.csv")

dim(snps_test[[186]])  # for H9S50D03-rep05, should be 499 rows

dim(snps_test[[186]])
dim(snps_mat[[186]])

snps_mat <- snps_test
dim(snps_mat[[186]])

# impute missing data
# replace NAs with colMeans 

colMeans(snps_mat[[2]], na.rm = TRUE)
mean(snps_mat[[1]][,1], na.rm=TRUE)

for (i in 1:length(snps_mat)) {
  cM <- colMeans(snps_mat[[i]], na.rm =TRUE)
  indx <- which(is.na(snps_mat[[i]]), arr.ind=TRUE)
  snps_mat[[i]][indx] <- cM[indx[,2]]
  snps_mat[[i]]
}

snps_mat[[1]][,1]

#################
## ordinations ##
#################

length(snps_mat)

for (s in 1:length(snps_mat)) {
  
  snps.scale <- scale(snps_mat[[s]], center = TRUE, scale = TRUE)     ## scale and center for RDA 
  snps.bray <- vegdist(snps_mat[[s]], method = 'bray')                ## bray-curtis distance for dbRDA
  
  ##  RDA  ##
  snp.rda <- rda(snps.scale, env[[s]], scale=F)
  snp.rda.sum <- summary(snp.rda)
  snpload.rda <- snp.rda.sum$species[,1:3]
  
  ##  dbRDA  ##
  snp.dbrda <- capscale(snps.bray ~ env[[s]], comm=snps_mat[[s]])  # RDA on a PCoA
  snp.dbrda.sum <- summary(snp.dbrda)
  snpload.dbrda <- snp.dbrda.sum$species[,1:3]   # just like in RDA
  rownames(snpload.dbrda) <- colnames(snps_mat[[s]])
  
  #### DETECT OUTLIERS 
  
  tests <- list(snpload.rda, snpload.dbrda)
  
  test.names <- c("rda","dbrda")
  
  for (i in 1:length(tests)) {
    
    for (j in 1:3) {                                    
      x <- tests[i]
      x <- matrix(unlist(x), ncol = 3)
      rownames(x) <- colnames(snps_mat[[s]])
      x <- x[,j]
      lims <- mean(x) + c(-1, 1) * 3 * sd(x)         
      out <- x[x < lims[1] | x > lims[2]]
      
      if (length(out) > 0) {
        names <- as.numeric(names(out))                    
        axis <- rep(j, length(out))
        outdata <- t(rbind(axis, names, as.numeric(out)))
        
        snpcors <- matrix(NA, nrow=length(names), ncol=6)
        
        for (k in 1:length(names)) {
          outlocus <- names[k] + 1   # to access correct column in snp dataframe
          
          if (ncol(snps_mat[[s]]) == 100) { 
            outsnp <- snps_mat[[s]][,outlocus]
          } else if (ncol(snps_mat[[s]]) == 99) {
            outsnp <- outsnp[!is.na(outsnp)] # remove NA individuals, if any
          }  
          
          corr.x <- cor.test(outsnp, env[[s]][,1])
          snpcors[k,1] <- corr.x$estimate
          snpcors[k,2] <- corr.x$p.value
          
          corr.y <- cor.test(outsnp, env[[s]][,2])
          snpcors[k,3] <- corr.y$estimate
          snpcors[k,4] <- corr.y$p.value
          
          corr.h <- cor.test(outsnp, env[[s]][,3])
          snpcors[k,5] <- corr.h$estimate
          snpcors[k,6] <- corr.h$p.value
        }
        
        outdata <- cbind(outdata, snpcors)
        fname <- paste(test.names[i],j, sep="")
        assign(fname, outdata)
      }
      else if (length(out) == 0) {
        fname <- paste(test.names[i],j, sep="")
        assign(fname, NA)
      }
      
    }
    
    
    
  }    
  out.rda <- rbind(rda1, rda2, rda3)
  out.rda <- as.data.frame(out.rda)
  label0 <- rep("RDA", nrow(out.rda))
  label1 <- rep(dir()[s], nrow(out.rda))
  out.rda <- cbind(label0,label1,out.rda)
  out.rda <- out.rda[complete.cases(out.rda),]
  
  out.dbrda <- rbind(dbrda1, dbrda2, dbrda3)
  out.dbrda <- as.data.frame(out.dbrda)
  label0 <- rep("dbRDA", nrow(out.dbrda))
  label1 <- rep(dir()[s], nrow(out.dbrda))
  out.dbrda <- cbind(label0,label1,out.dbrda)
  out.dbrda <- out.dbrda[complete.cases(out.dbrda),]
  
  #outs <- out.rda
  outs <- rbind(out.rda, out.dbrda)
  
  if (nrow(outs) > 0) {   
    colnames(outs) <- c("ord","sim","axis","locus","loading", "x-corr", "x-pval", "y-corr", "y-pval", "h-corr", "h-pval")
    fname <- paste("out_", dir()[s], sep="")
    assign(fname,outs)  
  }
  
}



## save output!
save <- ls()[grep("out_", ls())]

bar = NA

for (l in 1:length(save)) {
  foo <- get(save[l])
  bar <- rbind(foo, bar)
}

bar <- bar[-nrow(bar),]
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Test_Results/", "H9", "_RawData.csv", sep="")
write.csv(bar, file=fname)