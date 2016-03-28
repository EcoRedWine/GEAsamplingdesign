rm(list=ls())


library(raster)
library(vegan)

#############
# DATA PREP #
#############

sim <- c("H9S10D50")

##### To try lfmm without looping bells and whistles ####
##### Use lappy sections for looping through scenarios ####

setwd(dir = "Code/GEA/")
dir()
## list.dir <- dir()
## length(list.dir)  


## data <- lapply(list.dir, read.csv) # this takes some time 
data <- read.csv("H5S05D15-R5.csv")
length(data)
class(data)
# dim(data[[5]])

qrule <- raster("Qrule-H5-R5.asc")  

coord <- data[,2:3]
habitat <- extract(qrule,coord)

## L10.files <- dir("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Surfaces_Sample/", pattern = glob2rx("L10H9*_aa.asc"))   #CHANGE THE HABITAT CONFIGURATION
## L10.files

## qrule = list()
## for (q in 1:10) {
##  setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Surfaces_Sample/")
##  qrule[[q]] <- raster(L10.files[q])
## }

## qrule.rep <- rep(qrule, 24)
## plot(qrule.rep[[12]]) # just for checking 

## habitat = list()
## for (h in 1:length(qrule)) {
##  habitat[[h]] <- extract(qrule[[h]],coord)
# }




## now there is a list of 10 habitat files 

## habitat.all <- rep(habitat,24)
## habitat.all[[1]] == habitat.all[[2]]   ## These should not be the same 
## habitat.all[[11]] == habitat.all[[1]]  ## These should be the same


## Now we have a list of 240 habitat values ordered by replicate 0-9 and repeated; list of datasets with coordinates, age, sex, and snps (207 columns - need to remove snps and combine to 1 allele per locus). 

head(coord)
xvar <- (data[,2] - 477895.2)/1000         # convert UTM - subtract corner coordinate and /1000 to convert from m to km
yvar <- (data[,3] - 4905702)/1000 

## fulldata <- lapply(seq_along(data), function(x) cbind(xvar, yvar, habitat.all[[x]], data[[x]][,8:207]))
## length(fulldata)

fulldata <- cbind(xvar, yvar, habitat, data[,8:207])
fulldata[[1]] == fulldata[[2]]  # should be false 

## CHANGE THE SAMPLE FILE FOR DIFFERENT SAMPLING SCHEMES
samp <- read.csv("sample500.csv", header=F)
samp <- samp[,1]  # convert to a vector
sampdata <- fulldata[samp,]

## sampdata <- llply(fulldata, function(x) subset(x[samp,]))
## View(sampdata[[186]])
## length(sampdata)
## dim(sampdata[[1]])
## head(sampdata[[1]])
## sampdata[[1]] == sampdata[[2]]  # should be false

# extract the environmental data:
env <- sampdata[,1:3] 
env <- scale(env, center=T, scale=T) # scale and center env vars so they are comparable
env <- data.matrix(env)

### Use lappy method to loop through data sets
# env <- lapply(samp, function(x) subset(x[,1:3])) 
# env <- lapply(env, function(x) scale(x,center=T, scale=T)) # scale and center env vars so they are comparable
# env <- llply(env, function(x) data.matrix(x))
# class(env[[1]])

# for (i in 1:240) {
#  colnames(env[[i]]) <- c("xvar","yvar","habitat")
# }

# extract the genetic data: 

snps <- sampdata[,4:203] 
snps <- snps[,seq(1, ncol(snps), 2)]
colnames(snps) <- seq(0, 99, by=1)

# Use lapply method to loop through genetic data sets
# snps <- llply(sampdata, function(x) subset(x[,4:203]))
# snps <- llply(snps, function(x) subset(x[,seq(1, 200, 2)]))
# dim(snps[[1]])
# snps[[1]] == snps[[2]]  # should be false
# View(snps[[186]])

# this removes NA rows (individuals) from both snp and env datasets:
subset <- apply(snps,1,function(x) length(unique(x))>1)
snps <- snps[subset,]                                      ## applied to rows (N/As)
env <- env[subset,] 

#################
## CALCULATE K ##
#################    
dat <- t(snps)

library(RMTstat)

# estimate posterior allele frequencies             
p <- vector("numeric",nrow(dat))
for (i in 1:nrow(dat)) p[i] <- (1+sum(dat[i,],na.rm=TRUE))/(2+2*sum(!is.na(dat[i,])))

# center
mu <- apply(dat,1,mean,na.rm=TRUE)
dat <- dat - mu
dat[is.na(dat)] <- 0

# normalize
for (i in 1:nrow(dat)) dat[i,] <- dat[i,]/sqrt(p[i]*(1-p[i]))

# perform eigendecomposition of the covariance matrix
a <- eigen(cov(dat))
m <- length(a$values)

## PATTERSON K:
count_patt <- 0

for (j in 1:(m-1)) {
  L1 <- sum(a$values[j:m])
  L2 <- sum(a$values[j:m]^2)
  S2 <- (m-j)^2*L2/L1^2                   #Patterson equation 16
  nhat <- (m-j)*((m-j)+2)/(S2-(m-j))           #Patterson equation 13
  lambda <- a$values[j]*(m-j)/L1                 #Patterson equation 9
  mu <- (sqrt(nhat-1)+sqrt(m-j))^2/nhat                  #Patterson equation 5
  sigma <- (sqrt(nhat-1)+sqrt(m-j))/nhat*(1/sqrt(nhat-1)+1/sqrt(m-j))^(1/3)                         #Patterson equation 6
  twstat <- (lambda-mu)/sigma                    #Patterson equation 8
  twpvalue <- ptw(twstat,lower.tail=FALSE)
  
  if (is.na(twpvalue) == TRUE) {
    count_patt <- count_patt }
  else if (twpvalue < 0.0500) {
    count_patt <- count_patt + 1 }
  else if (twpvalue <= 0.0500) {
    count_patt <- count_patt }
}                         

if (count_patt == 0) {        ### K cannot be 0, so if it is 0, make it 1
  count_patt <- 1}

print(count_patt)

## MAP K:
count_map <- 0

for (j in 1:(m-1)) {
  L1 <- sum(a$values[j:m])
  L2 <- sum(a$values[j:m]^2)
  nhat <- L1^2/L2
  lambda <- a$values[j]*(m-j)/L1                 #Patterson equation 9
  mu <- (sqrt(nhat-1)+sqrt(m-j))^2/nhat                  #Patterson equation 5
  sigma <- (sqrt(nhat-1)+sqrt(m-j))/nhat*(1/sqrt(nhat-1)+1/sqrt(m-j))^(1/3)                         #Patterson equation 6
  twstat <- (lambda-mu)/sigma                    #Patterson equation 8
  twpvalue <- ptw(twstat,lower.tail=FALSE)
  
  if (is.na(twpvalue) == TRUE) {
    count_map <- count_map }
  else if (twpvalue < 0.0500) {
    count_map <- count_map + 1 }
  else if (twpvalue <= 0.0500) {
    count_map <- count_map }
}

if (count_map == 0) {
  count_map <- 1}

print(count_map)

# Write out files
fname <- paste("LFMMinputs/", sim,".pattK", sep="")
write.table(count_patt, file=fname, row.names=F, col.names=F)

fname <- paste("LFMMinputs/", sim,".mapK", sep="")
write.table(count_map, file=fname, row.names=F, col.names=F)

fname <- paste("LFMMinputs/", sim,".lfmm", sep="")
write.table(snps, file=fname, row.names=F, col.names=F)

# Incorporate the following line for multiple replicates
# fname <- paste("LFMMinputs/", sim[s], "-R", repl[r],".env", sep="")
fname <- paste("LFMMinputs/", sim, ".env", sep="")
write.table(env, file=fname, row.names=F, col.names=F)

## Keep the Ks for later:
outKs <- cbind(count_map, count_patt)
outks <- as.data.frame(outKs)
#For multiple sims, reps
#outKs <- cbind(sim[s], repl[r], outKs)
outKs <- cbind(sim, outKs)

# For multiple sims, reps
# colnames(outKs) <- c("sim", "rep", "mapK", "pattK")

colnames(outKs) <- c("sim", "mapK", "pattK")
fname <- paste("outKs_", sim, sep="")
# fname <- paste("outKs_", sim[s], repl[r], sep="")
assign(fname, outKs)

## save output!
save <- ls()[grep("outKs_", ls())]

bar = NA

for (l in 1:length(save)) {
  foo <- get(save[l])
  bar <- rbind(foo, bar)
}

bar <- bar[-nrow(bar),]
#fname <- paste("F:/7-Simulation Project/2-LFMM/LFMMresults/", sim[s], "_Ks.csv", sep="")
fname <- paste("LFMMresults/", sim, "_Ks.csv", sep="")
write.csv(bar, file=fname)


