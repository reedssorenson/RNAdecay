## ------------------------------------------------------------------------
library(RNAdecay)
decaydata = RNAdecay::decaydata # built-in package example data; comment this line out to use your own decaydata
decaydata$treatment = factor(decaydata$treatment, levels = c("WT", "sov", "vcs", "vs")) # you must type them identically in the new order  
levels(decaydata$treatment)[4] = "vcs sov"
decaydata[1:4,]

## ------------------------------------------------------------------------
results=RNAdecay::results # built-in package example data; comment this line out to use your own results
# For example:
# results = read.delim("~/DecayAnalysis/Example analysis results 2/alphas+betas+mods+grps+patterns+relABs.txt", header = TRUE, colClasses =  c(NA, rep("numeric", 9), rep("integer", 3), rep("character", 3), rep("numeric", 6)))
# results$alpha_subgroup = factor(results$alpha_subgroup, levels = unique(results$alpha_subgroup))
results[1:3, ]

## ------------------------------------------------------------------------
# gene_desc = read.delim("~/gene_description_20140101.txt", quote = NULL, comment = '', header = FALSE)
# gene_desc[, 1] = substr(gene_desc[, 1], 1, 9)
# gene_desc = data.frame(gene_desc[!duplicated(gene_desc[, 1]), ], row.names = 1)
# colnames(gene_desc) = c("type", "short description", "long description", "computational description")
# descriptions = gene_desc[, "long description"]
# names(descriptions) = rownames(gene_desc)

descriptions=c("Encodes a ubiquitin E3 ligase with RING and SPX domains that is involved in mediating immune responses and mediates degradation of PHT1s at plasma membranes.  Targeted by MIR827. Ubiquitinates PHT1;3, PHT1;2, PHT1;1/AtPT1 and PHT1;4/AtPT2.","","Related to Cys2/His2-type zinc-finger proteins found in higher plants. Compensated for a subset of calcineurin deficiency in yeast. Salt tolerance produced by ZAT10 appeared to be partially dependent on ENA1/PMR2, a P-type ATPase required for Li+ and Na+ efflux in yeast. The protein is localized to the nucleus, acts as a transcriptional repressor and is responsive to chitin oligomers. Also involved in response to photooxidative stress.","Encodes a stress enhanced protein that localizes to the thylakoid membrane and whose mRNA is upregulated in response to high light intensity.  It may be involved in chlorophyll binding.")
names(descriptions) = c("AT1G02860","AT5G54730", "AT1G27730", "AT4G34190")

## ---- fig.show = 'hold'--------------------------------------------------
p = DecayPlot(geneID = "AT1G02860", 
              xlim = c(0, 350), 
              ylim = c(0, 1.25), 
              xticks = 1:5*100, 
              yticks = 0:5/4,
              alphaSZ  =  12, 
              what  =  c("meanSE", "reps", "models"),
              treatments  =  c("WT", "vcs sov"), 
              colors = c("darkblue", "red3"),
              DATA = decaydata, 
              mod.results = results, 
              gdesc = descriptions)
print(p+myTheme(8, 1, leg.pos = c(0.7, 0.8))) #this will print the plot to your current graphics device (dev.cur() tells you what that is), if you do not have a graphics device open (e.g., "null device") initiate one (e.g., use quartz(),pdf(), or windows(); dev.off() closes the device and writes the file).

## ------------------------------------------------------------------------
# make new directory for results
wrdir = "~/DecayAnalysis/Example analysis results 3"
if(!file.exists(wrdir)) dir.create(wrdir, recursive = TRUE)

pdf(paste0(wrdir, "/decay plot example.pdf"), width = 10, height = 8)
p = DecayPlot(geneID = "AT1G02860", 
              xlim = c(0, 500), 
              ylim = c(0, 1.25), 
              xticks = 1:5*100, 
              yticks = 0:5/4,
              alphaSZ = 10, 
              what = c("Desc","meanSE", "reps", "models","alphas&betas"),
              treatments = c("WT", "vcs sov"),
              colors = c("darkblue", "red3"),
              DATA = decaydata,
              mod.results = results,
              gdesc = descriptions)
print(p+myTheme(32, 1, leg.pos = 'right'))
dev.off()

# plot multiple genes
ids = c("AT5G54730", "AT1G27730", "AT4G34190")
# or e.g.,
# ids = rownames(results[results$alpha_WT > 0.1, ]) 

pdf(paste0(wrdir, "/multiple gene decay plots.pdf"), width = 10, height = 7)
for(i in ids){
  p = DecayPlot(i, 
                what = c("meanSE", "models", "Desc"),
                mod.results = results, 
                gdesc = descriptions, 
                DATA = decaydata)
  print(p+myTheme(13.8, 1, leg.pos = 'right'))
  cat(i, " plotted.\n")
}
dev.off()

