## ------------------------------------------------------------------------
library(RNAdecay)

# make new directory for results
wrdir = "~/DecayAnalysis/Example analysis results 2"
if(!file.exists(wrdir)) dir.create(wrdir, recursive = TRUE)

## ------------------------------------------------------------------------
levels(decaydata$treatment) # This can not be longer than 4 elements

# decaydata$treatment = factor(decaydata$treatment, levels = levels(decaydata$treatment)[c(4, 1, 2, 3)]) # numbers here refer to the current order position in the factor levels() - the order of the numbers designates the new position

levels(decaydata$treatment)[4] = "vcs.sov"
levels(decaydata$treatment)

## ------------------------------------------------------------------------
decaydata = decaydata[order(decaydata$t.decay), ]
decaydata = decaydata[order(decaydata$rep), ]
decaydata = decaydata[order(as.numeric(decaydata$treatment)), ] 
decaydata = decaydata[order(decaydata$geneID), ]
ids = unique(decaydata$geneID) # 118 in example set
decaydata[1:10, ]

## ------------------------------------------------------------------------
# (no changes needed here - just execute the code)
nEquivGrp = if (length(unique(decaydata$treatment)) == 2) {2} else
  if (length(unique(decaydata$treatment)) == 3) {5} else 
    if (length(unique(decaydata$treatment)) == 4) {15}
genoSet = 1:(length(unique(decaydata$rep)) * length(unique(decaydata$t.decay)))
nTreat = length(unique(decaydata$treatment))
nSet = length(genoSet)*nTreat
groups = groupings(decaydata)
mods = data.frame(
  "a" = as.vector(sapply(1:nEquivGrp, function(x) {rep(x, nEquivGrp + 1)})),
  "b" = rep(1:(nEquivGrp + 1), nEquivGrp),
  row.names = paste0("mod", 1:(nEquivGrp*(nEquivGrp+1)))
)

## ------------------------------------------------------------------------
a.bounds = c(a.low(max(decaydata$t.decay)), a.high(min(unique(decaydata$t.decay)[unique(decaydata$t.decay)>0])))
b.bounds = c(b.low(max(decaydata$t.decay)), b.high())########################  NEED UPPER BOUND FUNCTION FOR BETA ###
a.bounds;b.bounds

## ------------------------------------------------------------------------
groupMap(decaydata = decaydata, path = paste0(wrdir, "/Model grouping colormap.pdf"), nEquivGrp = nEquivGrp, groups = groups, mods = mods)

## ------------------------------------------------------------------------
# only run once after installing RNAdecay pacakage
src_files = dir(paste0(find.package("RNAdecay"), "/src"), full.names = T)
for (i in src_files) {TMB::compile(i)}

## ------------------------------------------------------------------------
libs = gsub(".cpp", "", dir(paste0(find.package("RNAdecay"), "/src"), pattern = ".cpp", full.names = T))
for (i in libs) {dyn.load(TMB::dynlib(i))}

## ------------------------------------------------------------------------
modOptimization("AT2G18150", decaydata, group = groups, mod = mods, a.bounds, b.bounds, c("mod1", "mod2", "mod16", "mod239", "mod240"), file.only = FALSE) 
modOptimization("AT4G09680", decaydata, group = groups, mod = mods, a.bounds, b.bounds, c("mod1", "mod2", "mod16", "mod239", "mod240"), file.only = FALSE) 

## ------------------------------------------------------------------------
####### To run all genes in parallel use:
# parallel::mclapply(ids, data = decaydata, group = groups, mod = mods, 
#                    alpha.bounds = a.bounds, beta.bounds = b.bounds,
#                    models = paste0("mod", 1:240), 
#                    path = paste0(wrdir, "/modeling_results"),
#   mc.cores = getOption("mc.cores",  15L), # set the number of compute cores to use here (e.g., 9L = 9 cores, 11L = 11 cores)
#   mc.preschedule = TRUE,
#   mc.set.seed = TRUE,
#   mc.silent = FALSE,
#   mc.cleanup = TRUE,
#   mc.allow.recursive = TRUE)

## ------------------------------------------------------------------------
test.ids = sample(ids, 3) # NOTE: that everytime this line is run it generates a different random sampling, therefore the genes modeled below will be different each time this code is run. To test the exact set of genes shown in the vignette make a new character vector of the gene ids reported below and pass it to the gene argument using lapply instead of passing 'test.ids' as we do here.  
test.ids

a = proc.time()[3]
models = lapply(X = test.ids, # alternatively use `ids` here to complete the entire sample data set, but be prepared to wait 10 h. These gene IDs will get passed one at time to the "gene" argument of modOptimization() and return a list of the results data.frame.
                FUN = modOptimization,
                data = decaydata, 
                group = groups, 
                mod = mods,
                alpha.bounds = a.bounds, 
                beta.bounds = b.bounds,
                models = rownames(mods), 
                path = paste0(wrdir, "/modeling_results"),
                file.only = FALSE) 
names(models) = test.ids
b = proc.time()[3]
(b-a)/60/length(test.ids) # gives you average min per gene 

## ------------------------------------------------------------------------
models = lapply( paste0( wrdir, "/modeling_results/", test.ids, "_results.txt"), read.delim, header = TRUE )   
names(models) = test.ids

## ------------------------------------------------------------------------
models = RNAdecay::models # built-in data, comment this line out to continue with your own modelling output
results = t(sapply(models, function(x) x[x[, "AICc"]==min(x[, "AICc"]), ]))
results = as.data.frame(results)
results[, 1:2] = sapply(as.data.frame(results[, 1:2]), function(x) as.character(unlist(x)))
results[, -c(1,2)] = sapply(results[, -c(1,2)], unlist)
write.table(results, file = paste0(wrdir,"/best model results.txt"), sep = "\t")
results = read.delim(paste0(wrdir,"/best model results.txt"))

## ------------------------------------------------------------------------
library(ggplot2)
pdf(paste0(wrdir,"/distributions of stats.pdf"))
p = ggplot(results)
print(p+geom_histogram(aes(x = sigma2), bins = 300)+
        xlim(0,0.5)+
        geom_vline(xintercept = 0.0625)+
        myTheme())
print(p+stat_bin(aes(x = sigma2), breaks = c(seq(0,0.25,0.25/50),Inf), geom = "bar"))
print(p+stat_ecdf(aes(sigma2), geom = "line")+
        xlim(c(0,sqrt(0.5))))
print(p+stat_ecdf(aes(sqrt(sigma2)), geom = "line")+
        xlim(c(0,0.5)))
print(p+geom_histogram(aes(x = range.LL), bins = 60))
print(p+geom_histogram(aes(x = nUnique.LL), bins = 60))
dev.off()

pdf(paste0(wrdir,"/lowest AICc model counts.pdf"), height = 8, width = 32)
p = ggplot(data = data.frame(
  model = as.integer(gsub("mod","",names(table(results$mod)))),
  counts = as.integer(table(results$mod))))+
  geom_bar(aes(x = model, y = counts), stat = "identity")+
  scale_x_continuous(limits = c(0,nrow(mods)), breaks = seq(0,nrow(mods),5))+
  ggtitle("Model distribution of absolute lowest AICs")
print(p+myTheme(25))
dev.off()

## ------------------------------------------------------------------------
minMods = sapply(models, function(x) which (x[, "AICc"] < (2+min(x[, "AICc"])))) 
minAMods = lapply(minMods, function(x) unique(mods[x, "a"]))

pdf(paste0(wrdir,"/number of models that performed similar to the one selected.pdf"))
barplot(height = table(sapply(minMods, length)), xlab = "No. models in the lowest AICc group (not more than 2 different from lowest)",
        ylab = "No. genes")
barplot(height = table(sapply(minAMods, length)), xlab = "No. alpha groups in the lowest AICc group (not more than 2 different from lowest)",
        ylab = "No. genes")
dev.off()

## ------------------------------------------------------------------------
results = read.delim(paste0(wrdir,"/best model results.txt"))
results$alpha_grp = mods[as.character(results$mod), "a"]
results$beta_grp = mods[as.character(results$mod), "b"]
results$mod = as.numeric(gsub("mod", "", as.character(results$mod)))

results$alphaPattern = sapply(rownames(results), function(x) {
  paste0(gsub("alpha_", "", colnames(results)[3:(2+nTreat)][order(round(results[x, 3:(2+nTreat)], 4))]), collapse = "<=")
  })
results$alphaPattern = paste0(results$alpha_grp, "_", results$alphaPattern)
results$betaPattern = sapply(rownames(results), function(x){
  paste0(gsub("beta_", "", colnames(results)[(3+nTreat):(2+2*nTreat)][order(round(results[x, (3+nTreat):(2+2*nTreat)], 4))]), collapse = "<=")
  })
results$betaPattern = paste0(results$beta_grp, "_", results$betaPattern)

results = results[order(rownames(results)), ]
results = results[order(results$beta_grp), ]
results = results[order(results$alphaPattern), ]
results = results[order(results$alpha_grp), ]

results$alphaPattern = factor(results$alphaPattern, levels = as.character(unique(results$alphaPattern)))

results = data.frame(results[, 3:(2*nTreat+3), 2], results[, c("AICc", "alpha_grp", "beta_grp", "alphaPattern", "betaPattern")])
results$nEqMods = sapply(minMods[rownames(results)], length)
results$nEqAgp = sapply(minAMods[rownames(results)], length)

# Customize: add columns of relative alphas and betas as desired, e.g.:
results$rA_sov.WT    = results$alpha_sov      / results$alpha_WT
results$rA_vcs.WT    = results$alpha_vcs      / results$alpha_WT
results$rA_vcssov.WT = results$alpha_vcs.sov  / results$alpha_WT

write.table(results, paste0(wrdir,"/alphas+betas+mods+grps+patterns+relABs.txt"), sep = "\t")

results = read.delim(paste0(wrdir,"/alphas+betas+mods+grps+patterns+relABs.txt"), header = TRUE, colClasses =  c(NA, rep("numeric", 10), rep("integer", 2), rep("character", 2), rep("integer", 2), rep("numeric", 3)))
# results$alpha_subgroup = factor(results$alpha_subgroup, levels = unique(results$alpha_subgroup))
results$alphaPattern = factor(results$alphaPattern, levels = unique(results$alphaPattern))
results$betaPattern = factor(results$betaPattern, levels = unique(results$betaPattern))
results[1:3, ]

