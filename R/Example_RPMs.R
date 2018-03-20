
#' RNA abundance reads per million over RNA decay timecourse
#'
#' A dataset of RNA abundance of 118 genes in four Arabidopsis thaliana genotypes
#'     (WT, sov, vcs, vcs sov). Four biological replicates were collectred 0,
#'     7.5, 15, 30, 60, 120, 240, 480 min after blocking transcription. RNA was
#'     extracted, subjected to ribodepletion, and sequenced by RNA-seq (Illumina
#'     50 nt single end reads).
#'
#' @format a data frame with 118 rows and 128 columns; data are all RNA
#'     abundance values presented as reads per million. Column names indicate
#'     genotype, time point, and replicate number separated by underscores.
#'
#' @source Sorenson et al. (2018);
#'     \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86361}
#'
"RPMs"
