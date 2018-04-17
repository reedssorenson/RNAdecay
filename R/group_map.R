
#' model color map
#'
<<<<<<< HEAD:R/group_map.R
#' group_map makes a color map of alpha and beta equivalence groups by model. Similar colors in a row indicate constrained parameter equivalence between treatements. Gray indicates values of 0.
=======
#' GroupMap makes a color map of alpha and beta equivalence groups by model.
#'     Similar colors in a row indicate constrained parameter equivalence
#'     between treatements. Gray indicates values of 0.
>>>>>>> 6348625a4176d22804ebafedf857ed5938765d20:R/groupMap.R
#'
#' @param decaydata 5 column data.frame with colnames "geneID","treatment",
#'     "t.decay","rep","value"
#' @param path write path and file name, must end in ".pdf"
#' @param nEquivGrp number of equivalence groups based on number of treatments
#' @param groups equivalence group matrix
#' @param mods alpha beta equivalence group usage index (matrix)
#'
#' @return creates a model colormap and writes it to a pdf file named \code{path}
#'
#' @export
#' @examples
<<<<<<< HEAD:R/group_map.R
#' group_map(decaydata=data.frame(geneID=paste0("gene",1:4),
#'                     treatment=as.factor(rep(paste0("treat",1:2),2)),
=======
#' groupMap(decaydata=data.frame(geneID=paste0("gene",seq_len(4)),
#'                     treatment=as.factor(rep(paste0("treat",seq_len(2)),2)),
>>>>>>> 6348625a4176d22804ebafedf857ed5938765d20:R/groupMap.R
#'                     t.decay=0:3,
#'                     rep=rep("rep1"),
#'                     value=c(1,0.5,0.25,0.12)),
#'          path=paste0(getwd(),"/parameter equivalence colormap.pdf"),
#'          nEquivGrp = 2,
#'          groups = t(matrix(c(1,2,1,1,NA,NA),nrow=2,
#'                     dimnames=list(c("treat1","treat2"),
#'                     c("grp1","grp2","grp3")))),
#'          mods = t(matrix(c(1,1,1,2,1,3,2,1,2,2,2,3),nrow=2,
#'                          dimnames=list(c("a","b"),
#'                          paste0("mod",seq_len(6))))))


group_map <- function(decaydata,
                    path,
                    nEquivGrp = nEquivGrp,
                    groups = groups,
                    mods = mods) {
  nTreat <- length(unique(decaydata$treatment))
  groupingsA <- t(matrix(rep(groups[1,], nEquivGrp + 1), nrow = nTreat))
  for (i in 2:nEquivGrp) {
    groupingsA <- rbind(groupingsA, t(matrix(rep(
      groups[i,], nEquivGrp + 1
    ), nrow = nTreat)))
  }
<<<<<<< HEAD:R/group_map.R
  colnames(groupingsA) <- paste0("alpha_", unique(decaydata$treatment))
  groupingsB <- rbind(groups[1:nEquivGrp,] + nTreat, rep(NA, nTreat))
  rownames(groupingsB)[nEquivGrp + 1] <- paste0("grp", nEquivGrp + 1)
  groupingsB <- t(matrix(rep(t(groupingsB), nEquivGrp), nrow = nTreat))
  colnames(groupingsB) <- paste0("beta_", unique(decaydata$treatment))
  groupings <- cbind(groupingsA, groupingsB)
=======
  colnames(groupingsA) = paste0("alpha_", unique(decaydata$treatment))
  groupingsB = rbind(groups[seq_len(nEquivGrp),] + nTreat, rep(NA, nTreat))
  rownames(groupingsB)[nEquivGrp + 1] = paste0("grp", nEquivGrp + 1)
  groupingsB = t(matrix(rep(t(groupingsB), nEquivGrp), nrow = nTreat))
  colnames(groupingsB) = paste0("beta_", unique(decaydata$treatment))
  groupings = cbind(groupingsA, groupingsB)
>>>>>>> 6348625a4176d22804ebafedf857ed5938765d20:R/groupMap.R
  rm(groupingsA, groupingsB)
  rownames(groupings) <- rownames(mods)
  grDevices::pdf(path, width = 8, height = 10)
  heats <- c("darkblue",
            "green",
            "orange",
            "yellow",
            "magenta",
            "cyan",
            "red",
            "darkgray")[seq_len(nTreat * 2)]
  gplots::heatmap.2(
    groupings,
    trace = "none",
    breaks = seq_len(nTreat * 2 + 1) - 0.5,
    col = heats,
    Colv = FALSE,
    Rowv = FALSE,
    dendrogram = 'none',
    #labCol=T,labRow=T,
    main = "Model groupings - similar colors within\n each model have
    constrained equivalence",
    margins = c(10, 5),
    lwid = c(0.25, 0.75),
    lhei = c(0.20, 0.8),
    # relative width and height of the box with the key vs the box with the
    # heatmap
    cexRow = 0.25,
    key = FALSE
  )
  grDevices::dev.off()
}
