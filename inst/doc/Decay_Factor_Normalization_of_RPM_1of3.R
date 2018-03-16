## ------------------------------------------------------------------------
# you can see that in this fabricated data the RNAs have a half-life of ~30 min
RPM_geneX = data.frame(T0 = c(150, 135, 148), T30 = c(72, 76, 69), T60 = c(35, 35, 30), 
                     row.names = paste0("rep", 1:3))
RPM_geneX
RPM_geneX_T0norm = RPM_geneX/mean(RPM_geneX$T0)
RPM_geneX_T0norm
colMeans(RPM_geneX_T0norm)

## ------------------------------------------------------------------------
library(RNAdecay)
RPMs[1:2, c(1:11, 128)] # built-in example data
# NOTE: not showing all columns here

## ------------------------------------------------------------------------
cols(RPMs, "WT", "00") #gives the column indices that contain both tags in the column names

colnames(RPMs)[cols(RPMs, "WT", "00")]

# NOTE: this is based on grep pattern matching so tags MUST be unique and non-overlapping (i.e. one entire label should not be part of another).

## ------------------------------------------------------------------------
# create a directory for results output
wrdir = "~/DecayAnalysis/Example analysis results 1"
if(!file.exists(wrdir)) dir.create(wrdir, recursive = TRUE)

# specify sample tags from colnames of RPMs
treatment = c("WT", "sov", "vcs", "vs") # NOTE: I did not use 'vcs sov' for the name of the double mutant, or cols(RPMs, "vcs", "00") would have indexed all T0 samples of both "vcs" and "vcs" sov"; instead, I used 'vs'.
reps = c("r1", "r2", "r3", "r4")
tdecay = c("00", "07", "15", "30", "60", "120", "240", "480") #again NO nesting of one label in another hence "00" instead of "0".

## ------------------------------------------------------------------------
mean_RPM = data.frame("geneID" = rownames(RPMs))
SE_RPM = data.frame("geneID" = rownames(RPMs))
for(g in treatment){
  for(t in tdecay){
    mean_RPM = cbind(mean_RPM, rowMeans(RPMs[, cols(RPMs, g, t)]))
    names(mean_RPM)[length(mean_RPM)] = paste0(g, "_", t)
    SE_RPM = cbind(SE_RPM, apply(X = RPMs[, cols(RPMs, g, t)], MARGIN = 1, FUN = stats::sd)/sqrt(length(reps)))
    names(SE_RPM)[length(SE_RPM)] = paste0(g, "_", t)
  }}
mean_RPM[1:2, ]
SE_RPM[1:2, ]
# write output to file
write.table(x = mean_RPM, paste0(wrdir, "/RPM_mean.txt"), sep = "\t")
write.table(x = SE_RPM,   paste0(wrdir, "/RPM_SE.txt"), sep = "\t")

## ------------------------------------------------------------------------
filt1 = rep(TRUE, 118) # we will not filter in this example, so the filter value for each gene is set to TRUE.

## ------------------------------------------------------------------------
mT0norm = data.frame(row.names = rownames(RPMs)[filt1])
for(g in treatment){
  mean_T0reps = rowMeans(RPMs[filt1, cols(RPMs, g, "00")])
  for(r in reps){
    df = RPMs[filt1, colnames(RPMs)[cols(RPMs, g, r)]]
    df = df[, 1:length(tdecay)]/mean_T0reps
    mT0norm = cbind(mT0norm, df)
  }}
write.table(x = mT0norm, file = paste0(wrdir, "/mean T0 normalized.txt"),  sep = "\t")


## ------------------------------------------------------------------------
mean_mT0norm = data.frame(row.names = rownames(mT0norm))
for(g in treatment){
  for(t in tdecay){
    mean_mT0norm = cbind(mean_mT0norm, rowMeans(mT0norm[, cols(mT0norm, g, t)]))
    names(mean_mT0norm)[length(names(mean_mT0norm))] = paste0(g, "_", t)
  }}

SE_mT0norm = data.frame(row.names = rownames(mT0norm))
for(g in treatment){
  for(t in tdecay){
    SE_mT0norm = cbind(SE_mT0norm, apply(X = mT0norm[, cols(mT0norm, g, t)], MARGIN = 1, FUN = function(x) stats::sd(x)/sqrt(length(reps))))
    names(SE_mT0norm)[length(names(SE_mT0norm))] = paste0(g, "_", t)
  }}

# write output to file
write.table(x = mean_mT0norm, file = paste0(wrdir, "/mean T0 normalized_Mean.txt"),  sep = "\t")
write.table(x = SE_mT0norm, file = paste0(wrdir, "/mean T0 normalized_SE.txt"),  sep = "\t")

## ------------------------------------------------------------------------
stablegenes = c( "ATCG00490", "ATCG00680", "ATMG00280", "ATCG00580", "ATCG00140", "AT4G38970", "AT2G07671", "ATCG00570", "ATMG00730", "AT2G07727", "AT2G07687", "ATMG00160" ,"AT3G11630", "ATMG00060", "ATCG00600", "ATMG00220", "ATMG01170", "ATMG00410", "AT1G78900", "AT3G55440", "ATMG01320", "AT2G21170" ,"AT5G08670", "AT5G53300", "ATMG00070", "AT1G26630", "AT5G48300", "AT2G33040", "AT5G08690", "AT1G57720")

## ------------------------------------------------------------------------
stabletable = mean_mT0norm[stablegenes, ]
normFactors = colMeans(stabletable)
write.table(x = normFactors, paste0(wrdir, "/Normalziation Factors.txt"), sep = "\t")
normFactors_mean = matrix(normFactors, nrow = length(tdecay))
normFactors_SE = matrix(apply(X = stabletable, MARGIN = 2, function(x) stats::sd(x)/sqrt(length(stablegenes))), nrow = length(tdecay))
t.decay = c(0, 7.5, 15, 30, 60, 120, 240, 480)
rownames(normFactors_mean) = t.decay
rownames(normFactors_SE) = t.decay
colnames(normFactors_mean) = treatment
colnames(normFactors_SE) = treatment
list(normalizationFactors = normFactors_mean, SE = normFactors_SE)

## ------------------------------------------------------------------------
nF = vector()
ind = sapply(names(normFactors), function(x) grep(x, colnames(mT0norm)))
for(i in 1:ncol(ind)){
nF[ind[, i]] = colnames(ind)[i]
}
normFactorsM = t(matrix(rep(normFactors[nF], nrow(mT0norm)), ncol = nrow(mT0norm)))
rm(nF, ind)

## ------------------------------------------------------------------------
mT0norm_2 = data.frame(mT0norm/normFactorsM, 
                     row.names = rownames(mT0norm))

write.table(mT0norm_2, paste0(wrdir, "/mean T0 normalized and decay factor corrected.txt"), sep = "\t")

## ------------------------------------------------------------------------
mT0norm_2.1 = reshape2::melt(as.matrix(mT0norm_2), varnames = c("geneID", "variable"))
mT0norm_2.1 = cbind(mT0norm_2.1, reshape2::colsplit(mT0norm_2.1$variable, "_", names = c("treatment", "t.decay", "rep")))

mT0norm_2.1 = mT0norm_2.1[, colnames(mT0norm_2.1) !=  "variable"]
mT0norm_2.1 = mT0norm_2.1[, c(1, 3, 4, 5, 2)]
colnames(mT0norm_2.1) = c("geneID", "treatment", "t.decay", "rep", "value")
mT0norm_2.1$rep = gsub("r", "rep", mT0norm_2.1$rep)
mT0norm_2.1$t.decay = as.numeric(gsub("7", "7.5", as.numeric(mT0norm_2.1$t.decay)))
mT0norm_2.1$treatment = factor(mT0norm_2.1$treatment,levels = c("WT","sov","vcs","vs"))
mT0norm_2.1$rep = factor(mT0norm_2.1$rep, levels = paste0("rep",1:4))
write.table(x = mT0norm_2.1,  file = paste0(wrdir, 
"/ExampleDecayData+stableGenes.txt"), sep = "\t")

## ------------------------------------------------------------------------
table(RNAdecay::decaydata[,1:4] == mT0norm_2.1[,1:4]) # should be all TRUE no FALSE

