## This code will loop through all replicates for a set of 24 simulations (Gsims, H1sims, H5sims, H9sims) and 6 simulations for Nsims, with 1 missing data input file.  Keep the same file structure as when you download the files from Dryad (i.e. one folder per sim, and 10 replicate files in each folder named "R0.csv" to "R9.csv"). 


rm(list=ls())

library(raster)
library(vegan)

setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims")
sim <- dir()
repl <- c(0:9)

# Begin looping through all sims in the working directory 

for (s in 1:length(sim)) {
  for (r in 1:length(repl)) {
    
    # read in the habitat files; CHANGE HERE FOR H1, H5, OR H9
    # R1 goes with sim rep0
    fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Surfaces_Sample/L10H9R",r,"_aa.asc", sep="")
    qrule <- raster(fname)                                
    
    fname1 <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims/",sim[s],"/R",repl[r],".csv", sep="")
    data <- read.csv(fname1)
    
    coord <- data[,2:3]
    habitat <- extract(qrule,coord)
    
    xvar <- (data[,2] - 477895.2)/1000         # convert UTM - subtract corner coordinate and /1000 to convert from m to km
    yvar <- (data[,3] - 4905702)/1000 
    
    fulldata <- cbind(xvar, yvar, habitat, data[,8:207])
    
    # CHANGE samp FOR 100 OR 500 INDIVIDUAL SAMPLE SIZES 
    samp <- read.csv("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/sample500.csv", header=F)
    samp <- samp[,1]  
    
    sampdata <- fulldata[samp,]
    
    env <- sampdata[,1:3] 
    env <- scale(env, center=T, scale=T) # scale and center env vars so they are comparable
    env <- data.matrix(env)
    colnames(env)<-c("xvar","yvar","habitat")
    
    ## fix the snp dataset
    snps <- sampdata[,4:203]
    snps <- snps[,seq(1, ncol(snps), 2)]
    colnames(snps) <- seq(0, 99, by=1)
    snps <- data.matrix(snps)
    
    ## enter NAs; CHANGE MISSING DATA TYPE BY CHANGING random TO adaptive OR neutral, AND 500 to 100 
    missing_files <- list.files("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/Input Files /", pattern = glob2rx("random_*_500.csv"), full.names = F, recursive = F) 
    
    # CHANGE missing_files 1-4  
    missing <- read.csv(paste0("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/Input Files /", missing_files[1]),header=FALSE)  # just the first missing data input file; change [1] for 1-4 
    missing <- missing[,1]  # change to a vector 
    snps[missing] <- NA
    
    # this removes NA rows (individuals) from both snp and env datasets; there are a few replicates where NA individuals so there will be 499 instead of 500 individuals 
    subsetNA <- apply(snps,1,function(x) length(unique(x))>1)
    snps <- snps[subsetNA,]                                      ## applied to rows (NAs)
    env <- env[subsetNA,] 
    
    # remove monomorphic snps
    subsetMM <- apply(snps,2,function(x) length(unique(x))>1)
    snps <- snps[,subsetMM]                                       ## applied to columns (monomorphic snps)
    
    # check which loci have been removed & final # individuals
    tempname <- names(subsetMM)
    dims <- dim(snps)
    ind_missing <- cbind(sim[s], repl[r], dims[1], dims[2], tempname[subsetMM==F])
    fname <- paste("NAinds_", sim[s], repl[r], sep="")
    assign(fname, ind_missing)
    
    snps_missing <- cbind(sim[s],repl[r],length(which(is.na(snps)==TRUE)))
    fname <- paste("NAsnps_", sim[s], repl[r], sep="")
    assign(fname, snps_missing)
    
    # impute missing data
    # replace NAs with colMeans 
    
    colMeans(snps, na.rm = TRUE)
    mean(snps[,1], na.rm=TRUE)
    
    cM <- colMeans(snps, na.rm =TRUE)
    indx <- which(is.na(snps), arr.ind=TRUE)
    snps[indx] <- cM[indx[,2]]
    snps[,1]
    
    ## format the snp data
    snps.scale <- scale(snps, center=T, scale=T)  ## scale and center for PCA/RDA
    snps.bray <- vegdist(snps, method="bray")  ## bray-curtis distance for PCoA and dbRDA
    
    
    ##  RDA  ##
    snp.rda <- rda(snps.scale, env, scale=F)
    snp.rda.sum <- summary(snp.rda)
    snpload.rda <- snp.rda.sum$species[,1:3]
    
    ##  dbRDA  ##
    snp.dbrda <- capscale(snps.bray ~ env, comm=snps)  # RDA on a PCoA
    snp.dbrda.sum <- summary(snp.dbrda)
    snpload.dbrda <- snp.dbrda.sum$species[,1:3]   # just like in RDA
    rownames(snpload.dbrda) <- colnames(snps)
    
    
    ## detect outliers ##
    
    tests <- list(snpload.rda, snpload.dbrda)
    
    test.names <- c("rda","dbrda")
    
    for (i in 1:length(tests)) {
      
      for (j in 1:3) {                                    
        x <- tests[i]
        x <- matrix(unlist(x), ncol = 3)
        rownames(x) <- colnames(snps)
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
            
            if (ncol(snps) == 100) { 
              outsnp <- snps[,outlocus]
            } else if (ncol(snps) == 99) {
              outsnp <- newsnps[,outlocus] # use snps matrix without monomorphic removed
              outsnp <- outsnp[!is.na(outsnp)] # remove NA individuals, if any
            }  
            
            corr.x <- cor.test(outsnp, env[,1])
            snpcors[k,1] <- corr.x$estimate
            snpcors[k,2] <- corr.x$p.value
            
            corr.y <- cor.test(outsnp, env[,2])
            snpcors[k,3] <- corr.y$estimate
            snpcors[k,4] <- corr.y$p.value
            
            corr.h <- cor.test(outsnp, env[,3])
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
    label1 <- rep(sim[s], nrow(out.rda))
    label2 <- rep(repl[r], nrow(out.rda))
    out.rda <- cbind(label0,label1,label2,out.rda)
    out.rda <- out.rda[complete.cases(out.rda),]
    
    out.dbrda <- rbind(dbrda1, dbrda2, dbrda3)
    out.dbrda <- as.data.frame(out.dbrda)
    label0 <- rep("dbRDA", nrow(out.dbrda))
    label1 <- rep(sim[s], nrow(out.dbrda)) 
    label2 <- rep(repl[r], nrow(out.dbrda))
    out.dbrda <- cbind(label0,label1,label2,out.dbrda)
    out.dbrda <- out.dbrda[complete.cases(out.dbrda),]
    
    outs <- rbind(out.rda, out.dbrda)
    
    if (nrow(outs) > 0) {   
      colnames(outs) <- c("ord","sim","rep", "axis","locus","loading", "x-corr", "x-pval", "y-corr", "y-pval", "h-corr", "h-pval")
      fname <- paste("out_", sim[s],"_R",repl[r], sep="")
      assign(fname,outs)  
    }
    
  }
}


## save output; there will be 5 output files in total
# 1. Raw Data
# 2. True Positives 
# 3. Summary file 
# 4. Number of individuals after removing NA individuals 
# 5. Number of snps converted to NAs 

save <- ls()[grep("out_", ls())]

bar = NA

for (l in 1:length(save)) {
  foo <- get(save[l])
  bar <- rbind(foo, bar)
}

bar <- bar[-nrow(bar),]
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/", sim[s], "_random5percent_500ind_RawData.csv", sep="")
write.csv(bar, file=fname)


## true positives:
tp <- bar[bar$locus == 0,]
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/", sim[s], "_random5percent_500ind_TruePos.csv", sep="")
write.csv(tp, file=fname)

# summary of true and false positives
fp <- bar[bar$locus != 0,]
ord <- c("RDA","dbRDA")
summary <- matrix(NA, nrow=length(sim)*2, ncol=5)
colnames(summary) <- c("sim", "ord", "tp", "fp", "fp.sd")

summary[,1] <- as.vector(sapply(sim, function(x) rep(x,2)))
summary[,2] <- as.vector(rep(ord,length(sim)))


for (i in 1:length(sim)) {
  foo <- tp[tp$sim == sim[i],]
  baz <- fp[fp$sim == sim[i],]
  
  for (j in 1:length(ord)) {
    
    # true positives    
    bar <- foo[foo$ord == ord[j],]
    bar1 <- bar[!duplicated(bar$rep),]
    rowindex <- (i-1)*2
    summary[j+rowindex,3] <- nrow(bar1)
    
    #false positives
    qux <- baz[baz$ord == ord[j],]
    
    temp <- vector(mode="integer", length=10)
    
    for (k in 1:length(repl)) {
      rux <- qux[qux$rep == repl[k],]
      rux1 <- rux[!duplicated(rux$locus),] # removes any loci detected on >1 axis in same replicate
      temp[k] <- nrow(rux1)
    }
    
    summary[j+rowindex,4] <- sum(temp)
    summary[j+rowindex,5] <- sd(temp)
  }
}

fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/", sim[s], "_random5percent_500ind_Summary.csv", sep="")
write.csv(summary, file=fname)


## save file telling how many NA individuals removed 

library(gtools)

save <- ls()[grep("NAinds_", ls())]

bar = data.frame()

for (l in 1:length(save)) {
  foo <- get(save[l])
  foo <- as.data.frame(foo)
  bar <- smartbind(foo, bar)
}

bar <- bar[-nrow(bar),]
colnames(bar) <- c("sim","repl","num inds", "loci")
rownames(bar) <- seq(1,nrow(bar),1)
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/", sim[s], "_random5percent_500ind_NAinds.csv", sep="")
write.csv(bar, file=fname)

## save file telling how many snps are NA

save <- ls()[grep("NAsnps_", ls())]

bar = data.frame()

for (l in 1:length(save)) {
  foo <- get(save[l])
  foo <- as.data.frame(foo)
  bar <- smartbind(foo, bar)
}

bar <- bar[-nrow(bar),]
colnames(bar) <- c("sim","repl","num NA snps")
rownames(bar) <- seq(1, nrow(bar), 1)
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/", sim[s], "_random5percent_500ind_NAsnps.csv", sep="")
write.csv(bar, file=fname)
