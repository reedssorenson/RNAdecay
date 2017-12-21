
#' Normalized RNA abundance RNA decay timecourse
#'
#' A long form dataset of RNA abundance of 118 genes in four Arabidopsis thaliana genotypes (WT, sov, vcs, vcs sov). Four biological replicates were collectred 0, 7.5, 15, 30, 60, 120, 240, 480 min after blocking transcription. RNA was extracted, subjected to ribodepletion, and sequenced by RNA-seq (Illumina 50 nt single end reads). RPM values were normalized to mean T0 abundance and corrected by a decay factor.
#'
#' @format a data frame with 5 columns and 15104 rows.
#'  \describe{
#'     \item{geneID}{gene identifier; AGI}
#'     \item{treatment}{Arabidopsis genotype}
#'     \item{t.decay}{time of decay, in minutes}
#'     \item{rep}{replicate number}
#'     \item{value}{RPM value normalized to the replicate samples' mean T0 abundance and decay factor corrected}
#'     }
#'
#' @source Sorenson et al. (2017) Submitted; \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86361}
#'
"decaydata"
