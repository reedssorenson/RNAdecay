## ------------------------------------------------------------------------
# you can see that in this fabricated data the RNAs have a half-life of ~30 min
RPM_geneX <- data.frame(T0 = c(150, 135, 148), T30 = c(72, 76, 69), T60 = c(35, 35, 30), 
                     row.names = paste0("rep", 1:3))
RPM_geneX
RPM_geneX_T0norm <- RPM_geneX/mean(RPM_geneX$T0)
RPM_geneX_T0norm
colMeans(RPM_geneX_T0norm)

## ------------------------------------------------------------------------
library(RNAdecay)
RPMs <- RNAdecay::RPMs # built-in example data
RPMs[1:2, c(1:11, 128)] # NOTE: not showing all columns here

## ------------------------------------------------------------------------
cols(RPMs, "WT", "00") #gives the column indices that contain both tags in the column names

colnames(RPMs)[cols(RPMs, "WT", "00")]

# NOTE: this is based on grep pattern matching so tags MUST be unique and non-nested (i.e. one entire label should not be part of another).

## ------------------------------------------------------------------------
# create a directory for results output
wrdir <- "~/DecayAnalysis/Example analysis results 1"
if(!file.exists(wrdir)) dir.create(wrdir, recursive = TRUE)

# specify sample tags from colnames of RPMs
treatment <- c("WT", "sov", "vcs", "vs") # NOTE: I did not use 'vcs sov' for the name of the double mutant, or cols(RPMs, "vcs", "00") would have indexed all T0 samples of both "vcs" and "vcs" sov"; instead, I used 'vs'.
reps <- c("r1", "r2", "r3", "r4")
tdecay <- c("00", "07", "15", "30", "60", "120", "240", "480") #again NO nesting of one label in another hence "00" instead of "0".

## ------------------------------------------------------------------------
mean_RPM <- data.frame("geneID" = rownames(RPMs))
SE_RPM <- data.frame("geneID" = rownames(RPMs))
for(g in treatment){
  for(t in tdecay){
    mean_RPM <- cbind(mean_RPM, rowMeans(RPMs[, cols(RPMs, g, t)]))
    names(mean_RPM)[length(mean_RPM)] <- paste0(g, "_", t)
    SE_RPM <- cbind(SE_RPM, apply(X = RPMs[, cols(RPMs, g, t)], MARGIN = 1, FUN = stats::sd)/sqrt(length(reps)))
    names(SE_RPM)[length(SE_RPM)] <- paste0(g, "_", t)
  }}
mean_RPM[1:2, ]
SE_RPM[1:2, ]
# write output to file
write.table(x = mean_RPM, paste0(wrdir, "/RPM_mean.txt"), sep = "\t")
write.table(x = SE_RPM,   paste0(wrdir, "/RPM_SE.txt"), sep = "\t")

## ------------------------------------------------------------------------
filt1 <- rep(TRUE, 118) # we will not filter in this example, so the filter value for each gene is set to TRUE.

## ------------------------------------------------------------------------
mT0norm <- data.frame(row.names = rownames(RPMs)[filt1])
for(g in treatment){
  mean_T0reps <- rowMeans(RPMs[filt1, cols(RPMs, g, "00")])
  for(r in reps){
    df <- RPMs[filt1, colnames(RPMs)[cols(RPMs, g, r)]]
    df <- df[, 1:length(tdecay)]/mean_T0reps
    mT0norm <- cbind(mT0norm, df)
  }}
write.table(x = mT0norm, file = paste0(wrdir, "/T0 normalized.txt"),  sep = "\t")


## ------------------------------------------------------------------------
mean_mT0norm <- data.frame(row.names = rownames(mT0norm))
for(g in treatment){
  for(t in tdecay){
    mean_mT0norm <- cbind(mean_mT0norm, rowMeans(mT0norm[, cols(mT0norm, g, t)]))
    names(mean_mT0norm)[length(names(mean_mT0norm))] <- paste0(g, "_", t)
  }}

SE_mT0norm <- data.frame(row.names = rownames(mT0norm))
for(g in treatment){
  for(t in tdecay){
    SE_mT0norm <- cbind(SE_mT0norm, apply(X = mT0norm[, cols(mT0norm, g, t)], MARGIN = 1, FUN = function(x) stats::sd(x)/sqrt(length(reps))))
    names(SE_mT0norm)[length(names(SE_mT0norm))] <- paste0(g, "_", t)
  }}

# write output to file
write.table(x = mean_mT0norm, file = paste0(wrdir, "/T0 normalized_Mean.txt"),  sep = "\t")
write.table(x = SE_mT0norm, file = paste0(wrdir, "/T0 normalized_SE.txt"),  sep = "\t")

## ------------------------------------------------------------------------
stablegenes <- c( "ATCG00490", "ATCG00680", "ATMG00280", "ATCG00580", "ATCG00140", "AT4G38970", "AT2G07671", "ATCG00570", "ATMG00730", "AT2G07727", "AT2G07687", "ATMG00160" ,"AT3G11630", "ATMG00060", "ATCG00600", "ATMG00220", "ATMG01170", "ATMG00410", "AT1G78900", "AT3G55440", "ATMG01320", "AT2G21170" ,"AT5G08670", "AT5G53300", "ATMG00070", "AT1G26630", "AT5G48300", "AT2G33040", "AT5G08690", "AT1G57720")

## ------------------------------------------------------------------------
stabletable <- mean_mT0norm[stablegenes, ]
normFactors <- colMeans(stabletable)
write.table(x <- normFactors, paste0(wrdir, "/Normalziation Decay Factors.txt"), sep = "\t")
normFactors_mean <- matrix(normFactors, nrow = length(tdecay))
normFactors_SE <- matrix(apply(X = stabletable, MARGIN = 2, function(x) stats::sd(x)/sqrt(length(stablegenes))), nrow = length(tdecay))
t.decay <- c(0, 7.5, 15, 30, 60, 120, 240, 480)
rownames(normFactors_mean) <- t.decay
rownames(normFactors_SE) <- t.decay
colnames(normFactors_mean) <- treatment
colnames(normFactors_SE) <- treatment
list(normalizationFactors = normFactors_mean, SE = normFactors_SE)

## ------------------------------------------------------------------------
nF <- vector()
ind <- sapply(names(normFactors), function(x) grep(x, colnames(mT0norm)))
for(i in 1:ncol(ind)){
nF[ind[, i]] <- colnames(ind)[i]
}
normFactorsM <- t(matrix(rep(normFactors[nF], nrow(mT0norm)), ncol = nrow(mT0norm)))
rm(nF, ind)

## ------------------------------------------------------------------------
mT0norm_2 <- data.frame(mT0norm/normFactorsM, 
                     row.names = rownames(mT0norm))

write.table(mT0norm_2, paste0(wrdir, "/T0 normalized and decay factor corrected.txt"), sep = "\t")

## ------------------------------------------------------------------------
mT0norm_2.1 <- reshape2::melt(as.matrix(mT0norm_2), varnames = c("geneID", "variable"))
mT0norm_2.1 <- cbind(mT0norm_2.1, reshape2::colsplit(mT0norm_2.1$variable, "_", names = c("treatment", "t.decay", "rep")))

mT0norm_2.1 <- mT0norm_2.1[, colnames(mT0norm_2.1) !=  "variable"]
mT0norm_2.1 <- mT0norm_2.1[, c(1, 3, 4, 5, 2)]
colnames(mT0norm_2.1) <- c("geneID", "treatment", "t.decay", "rep", "value")
mT0norm_2.1$rep <- gsub("r", "rep", mT0norm_2.1$rep)
mT0norm_2.1$t.decay <- as.numeric(gsub("7", "7.5", as.numeric(mT0norm_2.1$t.decay)))
mT0norm_2.1$treatment <- factor(mT0norm_2.1$treatment,levels = c("WT","sov","vcs","vs"))
mT0norm_2.1$rep <- factor(mT0norm_2.1$rep, levels = paste0("rep",1:4))
write.table(x = mT0norm_2.1,  file = paste0(wrdir, 
"/ExampleDecayData+stableGenes.txt"), sep = "\t")

## ------------------------------------------------------------------------
table(RNAdecay::decay_data[,1:4] == mT0norm_2.1[,1:4]) # should be all TRUE no FALSE

## ------------------------------------------------------------------------
library(RNAdecay)
decay_data <- RNAdecay::decay_data # built-in example data

# make new directory for results
wrdir <- "~/DecayAnalysis/Example analysis results 2"
if(!file.exists(wrdir)) dir.create(wrdir, recursive = TRUE)

## ------------------------------------------------------------------------
levels(decay_data$treatment) # This can not be longer than 4 elements for modeling.

# reorder factor levels if necessary (skipped here):
# decay_data$treatment <- factor(decay_data$treatment, levels = levels(decay_data$treatment)[c(4, 1, 2, 3)]) # numbers here refer to the current order position in the factor levels() - the order of the numbers designates the new position

levels(decay_data$treatment)[4] <- "vcs.sov"
levels(decay_data$treatment)

## ------------------------------------------------------------------------
decay_data <- decay_data[order(decay_data$t.decay), ]
decay_data <- decay_data[order(decay_data$rep), ]
decay_data <- decay_data[order(as.numeric(decay_data$treatment)), ] 
decay_data <- decay_data[order(decay_data$geneID), ]
ids <- as.character(unique(decay_data$geneID)) # 118 in example set
decay_data[1:10, ]

## ------------------------------------------------------------------------
# (no changes needed here - just execute the code)
nEquivGrp <- if (length(unique(decay_data$treatment)) == 2) {2} else
  if (length(unique(decay_data$treatment)) == 3) {5} else 
    if (length(unique(decay_data$treatment)) == 4) {15}
genoSet <- 1:(length(unique(decay_data$rep)) * length(unique(decay_data$t.decay)))
nTreat <- length(unique(decay_data$treatment))
nSet <- length(genoSet)*nTreat
groups <- groupings(decay_data)
mods <- data.frame(
  "a" = as.vector(sapply(1:nEquivGrp, function(x) {rep(x, nEquivGrp + 1)})),
  "b" = rep(1:(nEquivGrp + 1), nEquivGrp),
  row.names = paste0("mod", 1:(nEquivGrp*(nEquivGrp+1)))
)

## ------------------------------------------------------------------------
group_map(decaydata = decay_data, path = paste0(wrdir, "/Model grouping colormap.pdf"), nEquivGrp = nEquivGrp, groups = groups, mods = mods)

## ------------------------------------------------------------------------
a_bounds <- c(a_low(max(decay_data$t.decay)),
             a_high(min(unique(decay_data$t.decay)[unique(decay_data$t.decay)>0])))
b_bounds <- c(b_low(max(decay_data$t.decay)), 0.075)
a_bounds;b_bounds

## ------------------------------------------------------------------------
mod_optimization("AT2G18150", decay_data, group = groups, mod = mods, a_bounds, b_bounds, c("mod1", "mod2", "mod16", "mod239", "mod240"), file_only = FALSE) 
mod_optimization("AT4G09680", decay_data, group = groups, mod = mods, a_bounds, b_bounds, c("mod1", "mod2", "mod16", "mod239", "mod240"), file_only = FALSE) 

## ------------------------------------------------------------------------
####### To run all genes in parallel use:
# parallel::mclapply(ids, FUN = mod_optimization, 
#                    data = decay_data, group = groups, mod = mods, 
#                    alpha_bounds = a_bounds, beta_bounds = b_bounds,
#                    models = paste0("mod", 1:240), 
#                    path = paste0(wrdir, "/modeling_results"),
#   mc.cores = getOption("mc.cores",  15L), # set the number of compute cores to use here (e.g., 9L = 9 cores, 11L = 11 cores)
#   mc.preschedule = TRUE,
#   mc.set.seed = TRUE,
#   mc.silent = FALSE,
#   mc.cleanup = TRUE,
#   mc.allow.recursive = TRUE)

## ------------------------------------------------------------------------
test_ids <- sample(ids, 1) # NOTE: that everytime this line is run it generates a different random sampling, therefore the gene modeled below will be different each time this code is run. To test the exact gene shown in the vignette make a new character vector of the gene id reported below and pass it to the gene argument of mod_optimization using lapply instead of passing 'test_ids' as we do here.  
test_ids

a <- proc.time()[3]
models <- lapply(X = test_ids, # alternatively use `ids` here to complete the entire sample data set, but be prepared to wait ~ 10 h. These gene IDs will get passed one at time to the "gene" argument of mod_optimization() and return a list of the results data frame.
                FUN = mod_optimization,
                data = decay_data, 
                group = groups, 
                mod = mods,
                alpha_bounds = a_bounds, 
                beta_bounds = b_bounds,
                models = rownames(mods), 
                path = paste0(wrdir, "/modeling_results"),
                file_only = FALSE) 
names(models) <- test_ids
b <- proc.time()[3]
(b-a)/60/length(test_ids) # gives you average min per gene 

## ------------------------------------------------------------------------
models <- lapply( paste0( wrdir, "/modeling_results/", test_ids, "_results.txt"), read.delim, header = TRUE )   
names(models) <- test_ids

## ------------------------------------------------------------------------
models <- RNAdecay::models # built-in example data, comment this line out to continue with your own modelling output
results <- t(sapply(models, function(x) x[x[, "AICc"] == min(x[, "AICc"]), ]))
results <- as.data.frame(results)
results[, 1:2] <- sapply(as.data.frame(results[, 1:2]), function(x) as.character(unlist(x)))
results[, -c(1,2)] <- sapply(results[, -c(1,2)], unlist)
write.table(results, file = paste0(wrdir,"/best model results.txt"), sep = "\t")
results <- read.delim(paste0(wrdir,"/best model results.txt"))

## ------------------------------------------------------------------------
library(ggplot2)
pdf(paste0(wrdir,"/distributions of stats.pdf"))
p <- ggplot(results)
print(p+geom_histogram(aes(x = sigma2), bins = 300)+
        coord_cartesian(xlim = c(0,0.5))+
        geom_vline(xintercept = 0.0625,color = "red2")+
        plain_theme())
print(p+stat_bin(aes(x = sigma2), breaks = c(seq(0,0.25,0.25/50),1), geom = "bar")+coord_cartesian(xlim = c(0,0.25)))
print(p+stat_ecdf(aes(sigma2), geom = "line")+
        coord_cartesian(xlim = c(0,0.5)))
print(p+stat_ecdf(aes(sqrt(sigma2)), geom = "line")+
        coord_cartesian(xlim = c(0,sqrt(0.5))))
print(p+geom_histogram(aes(x = range.LL), bins = 60))
print(p+geom_histogram(aes(x = nUnique.LL), bins = 60))
dev.off()

pdf(paste0(wrdir,"/lowest AICc model counts.pdf"), height = 8, width = 32)
p <- ggplot(data = data.frame(
  model = as.integer(gsub("mod","",names(table(results$mod)))),
  counts = as.integer(table(results$mod))))+
  geom_bar(aes(x = model, y = counts), stat = "identity")+
  scale_x_continuous(limits = c(0,nrow(mods)), breaks = seq(0,nrow(mods),5))+
  ggtitle("model distribution of absolute lowest AICs")
print(p+plain_theme(25))
dev.off()

## ------------------------------------------------------------------------
min_mods <- sapply(models, function(x) which (x[, "AICc"] < (2+min(x[, "AICc"])))) 
min_alpha_mods <- lapply(min_mods, function(x) unique(mods[x, "a"]))

pdf(paste0(wrdir,"/number of models that performed similar to the one selected.pdf"))
barplot(height = table(sapply(min_mods, length)), xlab = "No. models in the lowest AICc group (not more than 2 different from lowest)",
        ylab = "No. genes")
barplot(height = table(sapply(min_alpha_mods, length)), xlab = "No. alpha groups in the lowest AICc group (not more than 2 different from lowest)",
        ylab = "No. genes")
dev.off()

## ------------------------------------------------------------------------
results <- read.delim(paste0(wrdir,"/best model results.txt"))
results$alpha_grp <- mods[as.character(results$mod), "a"]
results$beta_grp <- mods[as.character(results$mod), "b"]
results$mod <- as.numeric(gsub("mod", "", as.character(results$mod)))

results$alphaPattern <- sapply(rownames(results), function(x) {
  paste0(gsub("alpha_", "", colnames(results)[3:(2+nTreat)][order(round(results[x, 3:(2+nTreat)], 4))]), collapse = "<=")
  })
results$alphaPattern <- paste0(results$alpha_grp, "_", results$alphaPattern)
results$betaPattern <- sapply(rownames(results), function(x){
  paste0(gsub("beta_", "", colnames(results)[(3+nTreat):(2+2*nTreat)][order(round(results[x, (3+nTreat):(2+2*nTreat)], 4))]), collapse = "<=")
  })
results$betaPattern <- paste0(results$beta_grp, "_", results$betaPattern)

results <- results[order(rownames(results)), ]
results <- results[order(results$beta_grp), ]
results <- results[order(results$alphaPattern), ]
results <- results[order(results$alpha_grp), ]

results$alphaPattern <- factor(results$alphaPattern, levels = as.character(unique(results$alphaPattern)))

results <- data.frame(results[, 3:(2*nTreat+3), 2], results[, c("AICc", "alpha_grp", "beta_grp", "alphaPattern", "betaPattern")])
results$nEqMods <- sapply(min_mods[rownames(results)], length)
results$nEqAgp <- sapply(min_alpha_mods[rownames(results)], length)

# Customize: add columns of relative alphas and betas as desired, e.g.:
results$rA_sov.WT    <- results$alpha_sov      / results$alpha_WT
results$rA_vcs.WT    <- results$alpha_vcs      / results$alpha_WT
results$rA_vcssov.WT <- results$alpha_vcs.sov  / results$alpha_WT

write.table(results, paste0(wrdir,"/alphas+betas+mods+grps+patterns+relABs.txt"), sep = "\t")

results <- read.delim(paste0(wrdir,"/alphas+betas+mods+grps+patterns+relABs.txt"), header = TRUE, colClasses =  c(NA, rep("numeric", 2+2*nTreat), rep("integer", 2), rep("character", 2), rep("integer", 2), rep("numeric", 3)))
# results$alpha_subgroup <- factor(results$alpha_subgroup, levels = unique(results$alpha_subgroup))
results$alphaPattern <- factor(results$alphaPattern, levels = unique(results$alphaPattern))
results$betaPattern <- factor(results$betaPattern, levels = unique(results$betaPattern))
results[1:3, ]

## ------------------------------------------------------------------------
library(RNAdecay)
decay_data <- RNAdecay::decay_data # built-in package example data; comment this line out to use your own decay_data
decay_data$treatment <- factor(decay_data$treatment, levels = c("WT", "sov", "vcs", "vs")) # you must type them identically in the new order  
levels(decay_data$treatment)[4] <- "vcs sov"
decay_data[1:4,]

## ------------------------------------------------------------------------
results <- RNAdecay::results # built-in package example data; comment this line out to use your own results
# For example:
# results <- read.delim("~/DecayAnalysis/Example analysis results 2/alphas+betas+mods+grps+patterns+relABs.txt", header = TRUE, colClasses =  c(NA, rep("numeric", 9), rep("integer", 3), rep("character", 3), rep("numeric", 6)))
# results$alpha_subgroup <- factor(results$alpha_subgroup, levels = unique(results$alpha_subgroup))
results[1:3, ]

## ------------------------------------------------------------------------
# gene_desc <- read.delim("~/gene_description_20140101.txt", quote = NULL, comment = '', header = FALSE)
# gene_desc[, 1] <- substr(gene_desc[, 1], 1, 9)
# gene_desc <- data.frame(gene_desc[!duplicated(gene_desc[, 1]), ], row.names = 1)
# colnames(gene_desc) <- c("type", "short description", "long description", "computational description")
# descriptions <- gene_desc[, "long description"]
# names(descriptions) <- rownames(gene_desc)

descriptions <- c(
  "Encodes a ubiquitin E3 ligase with RING and SPX domains that is involved in mediating immune responses and mediates degradation of PHT1s at plasma membranes.  Targeted by MIR827. Ubiquitinates PHT1;3, PHT1;2, PHT1;1/AtPT1 and PHT1;4/AtPT2.",
  "",
  "Related to Cys2/His2-type zinc-finger proteins found in higher plants. Compensated for a subset of calcineurin deficiency in yeast. Salt tolerance produced by ZAT10 appeared to be partially dependent on ENA1/PMR2, a P-type ATPase required for Li+ and Na+ efflux in yeast. The protein is localized to the nucleus, acts as a transcriptional repressor and is responsive to chitin oligomers. Also involved in response to photooxidative stress.",
  "Encodes a stress enhanced protein that localizes to the thylakoid membrane and whose mRNA is upregulated in response to high light intensity.  It may be involved in chlorophyll binding."
  )
names(descriptions) <- c("AT1G02860","AT5G54730", "AT1G27730", "AT4G34190")

## ---- fig.show = 'hold'--------------------------------------------------
p <- decay_plot(geneID = "AT1G02860", 
              xlim = c(0, 350), 
              ylim = c(0, 1.25), 
              xticks = 1:5*100, 
              yticks = 0:5/4,
              alphaSZ  =  12, 
              what  =  c("meanSE", "reps", "models"),
              treatments  =  c("WT", "vcs sov"), 
              colors = c("darkblue", "red3"),
              DATA = decay_data, 
              mod.results = results, 
              gdesc = descriptions)
print(p+plain_theme(8, 1, leg.pos = c(0.7, 0.8))) #this will print the plot to your current graphics device (dev.cur() tells you what that is), if you do not have a graphics device open (e.g., "null device") initiate one (e.g., use quartz(),pdf(), or windows(); dev.off() closes the device and writes the file).

## ------------------------------------------------------------------------
# make new directory for results
wrdir <- "~/DecayAnalysis/Example analysis results 3"
if(!file.exists(wrdir)) dir.create(wrdir, recursive = TRUE)

pdf(paste0(wrdir, "/decay plot example.pdf"), width = 10, height = 8)
p <- decay_plot(geneID = "AT1G02860", 
              xlim = c(0, 500), 
              ylim = c(0, 1.25), 
              xticks = 1:5*100, 
              yticks = 0:5/4,
              alphaSZ = 10, 
              what = c("Desc","meanSE", "reps", "models","alphas&betas"),
              treatments = c("WT", "vcs sov"),
              colors = c("darkblue", "red3"),
              DATA = decay_data,
              mod.results = results,
              gdesc = descriptions)
print(p+plain_theme(32, 1, leg.pos = 'right'))
dev.off()

# plot multiple genes
ids <- c("AT5G54730", "AT1G27730", "AT4G34190")
# or e.g.,
# ids <- rownames(results[results$alpha_WT > 0.1, ]) 

pdf(paste0(wrdir, "/multiple RNA decay plots.pdf"), width = 10, height = 7)
for(i in ids){
  p <- decay_plot(i, 
                what = c("meanSE", "models", "Desc"),
                mod.results = results, 
                gdesc = descriptions, 
                DATA = decay_data)
  print(p+plain_theme(13.8, 1, leg.pos = 'right'))
  cat(i, " plotted.\n")
}
dev.off()

