rm(list=ls())

library(LEA)

setwd(dir = "LFMMinputs/")

sim <- c("H9S10D50")

################
# runs of lfmm #
################

# Runs with K = 9 and 5 repetitions.
# The runs are composed of 6000 iterations including 3000 iterations
# for burnin.
# around 30 seconds per run.
project = NULL
project = lfmm("H9S10D50.lfmm", "H9S10D50.env", K = 6, repetitions = 5, 
               project = "new")

# get the zscores of each run for K = 6
z1 = z.scores(project, K = 6, d = 1)
z2 = z.scores(project, K = 6, d = 2)
z3 = z.scores(project, K = 6, d = 3)
z4 = z.scores(project, K = 6, d = 4)

# Combine the z-scores using the Stouffer method
zs.stouffer = apply(zs, MARGIN = 1, median)

# calculate the inflation factor
lambda = median(zs.stouffer^2)/.456

# calculate adjusted p-values
cp.values = pchisq(zs.stouffer^2/lambda, df = 1, lower = FALSE)

for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("expected FDR:", alpha))
  L = length(cp.values)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp.values) < alpha * (1:L) / L)
  candidates = order(cp.values)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("FDR:", estimated.FDR, "True Positive:", estimated.TP))
}

################
# Post-treatment #
################

# show the project
show(project)

# summary of the project
summary(project)

# get the p-values for the 2nd run for K = 6
p1 = p.values(project, K = 6, run = 2, d = 1) 
p2 = p.values(project, K = 6, run = 2, d = 2) 
p3 = p.values(project, K = 6, run = 2, d = 3) 
p <- c(p1,p2,p3) # create a matrix with headings

# get the -log10(p-values) for the 2nd run for K = 6
mp = mlog10p.values(project, K = 6, run = 2, d =1)