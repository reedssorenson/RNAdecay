
#' decay_plot() function
#'
#' Plots RNA decay data and/or decay models using the ggplot2 package.
#'
#' @param geneID single gene ID from data set (e.g. "AT1G00100") for which to plot data/model
#' @param xlim,ylim vector of length 2 defineing the limits of the plot (zooms in on data)
#' @param xticks,yticks vectors specifyng tick marks for the x and y axes
#' @param alphaSZ text size of alpha and beta parameter labels if plotted
#' @param what character vector specifying what to plot; any or all (default) of "Desc","models","reps","meanSE","alphas&betas"
#'    "Desc"         - plots gene descriptions behind data
#'    "models"       - plots the selected fit model
#'    "reps"         - plots individual replicate data as distinct shapes
#'    "meanSE        - plots the replicate means and standard errors
#'    "alphas&betas" - plots the values of the alphas and betas for each model below the model at the greatest x position
#' @param treatments what treatments/genotypes to plot from the supplied data
#' @param colors vector of R recognized colors (e.g. "red","darkblue")
#' @param DATA (required) normalized abundance decay data with column names: "geneID", "treatment","t.decay", "rep","value"
#' @param mod.results (optional; required for plotting models) data.frame of the model results as output from the modeling (e.g. "alphas+betas+mods+grps+patterns+relABs.txt")
#' @param gdesc (optional; required for plotting gene descriptions) gene descriptions (geneID-named vector of gene descriptions geneID must match those of data)
#' @param desc.width width of gene descriptions (in number of characters) before word wrap
#'
#' @return returns a ggplot to be used with print; could also be modified using the syntax of ggplot2 e.g.'+geom_XXXX(...)'
#'
#' @export
#'
#' @examples
# p<-RNAdecay::decay_plot(
#   geneID = "GOI_id",
#   mod.results = data.frame(alpha_WT = 0.0830195, beta_WT = 0.04998945,
#                            model = 1, alpha_grp = 1, beta_grp = 1, alpha_subgroup = 1.1,
#                            row.names = "GOI_id"),
#   what = c("meanSE","alphas&betas","models"),
#   treatments = "WT",
#   colors = "black",
#   DATA = data.frame(geneID=rep("GOI_id",15),
#                     treatment=rep("WT",15),
#                     t.decay=rep(c(0,7.5,15,30,60),3),
#                     rep=paste0("rep",c(rep(1,5),rep(2,5),rep(3,5))),
#                     value= c(0.9173587, 0.4798672, 0.3327807, 0.1990708, 0.1656554,
#                              0.9407511, 0.7062988, 0.3450886, 0.3176824, 0.2749946,
#                              1.1026497, 0.6156978, 0.4563346, 0.2865779, 0.1680075)),
#   xlim = c(0, 65),
#   alphaSZ = 10)
# print(p)



decay_plot=
  function (
    geneID,
    xlim = c(0, 500),
    ylim = c(0, 1.25),
    xticks = NA,
    yticks = 0:5/4,
    alphaSZ = 8,
    what = c("Desc", "models", "reps", "meanSE", "alphas&betas"),
    DATA,
    treatments = NA,
    colors = NA,
    mod.results = NA,
    gdesc = NA,
    desc.width = 55) {
    if(any(! geneID %in% rownames(mod.results))) {stop(paste0("geneID:",geneID," not found in the dataset."),call. = F)}
    if(any(!treatments %in% gsub("alpha_","",colnames(mod.results)[grep("alpha_",colnames(mod.results)[1:4])]))) {
      stop(paste0("Supplied 'treatments' are not found in 'mod.results'."))
    }
    if(any(!treatments %in% unique(DATA$treatment))) {
      stop(paste0("'treatments' indicated are not found in 'DATA'."))
    }
    dExp <- function(t, par) {
      a <- par[1]
      b <- par[2]
      exp(-(a/b) * (1 - exp(-b * t)))
    }
    fun_exp <- function(t, a) {
      exp(-a * t)
    }
    wrapper <- function(x, ...) {
      paste(strwrap(x, ...), collapse = "\n")
    }
    if (any(what %in% "models"))     mod <- mod.results[geneID, "mod"]  else mod <- NA
    if (any(what %in% "models"))     A_grp <- mod.results[geneID, "alpha_subgroup"]   else ""
    if (is.na(treatments[1]))     treatments <- unique(DATA$treatment)  else treatments <- treatments[treatments %in% unique(DATA$treatment)]
    if (is.na(colors[1]))     colors <- grDevices::rainbow(length(treatments),alpha = 1)
    if (any(is.na(xticks)))     xticks <- c(0, 1:5 * diff(xlim)/5 + xlim[1])
    names(colors) <- treatments
    p <- ggplot2::ggplot(data = DATA[DATA$geneID == geneID & DATA$treatment %in% treatments, ])
    if (any(what %in% "Desc")) {
      p <- p + ggplot2::annotate(geom = "text", x = xlim[2]/2, y = ylim[2],
                                 label = wrapper(paste0(geneID,
                                                        if (any(what %in% "models" & !is.na(mod))) {
                                                          paste0(" - model ",mod)
                                                        } else { "" }, paste0(" - ", gdesc[geneID])), desc.width),
                                 vjust = 1,
                                 hjust = 0.5,
                                 color = grDevices::gray(0.4),
                                 size = alphaSZ *0.4,
                                 family = c("mono"), fontface = "plain",
                                 angle = 0)
    } else {
      p <- p + ggplot2::ggtitle(if (any(what %in% "models") & !is.na(mod)) {
        paste0(geneID," - model ", mod)
      } else {geneID}
      )
    }
    if (any(what %in% "alphas&betas")) {
      p <- p + ggplot2::geom_text(parse = TRUE,
                                  data = data.frame(
                                    treatment = treatments,
                                    text = sapply(treatments,function(g) {
                                      paste0(
                                        "alpha==", as.character(round(mod.results[geneID, paste0("alpha_", g)], 4)),
                                        "~beta==", as.character(round(mod.results[geneID, paste0("beta_", g)], 4))
                                      )
                                    }),
                                    x = rep(xlim[2], length(treatments)),
                                    y = if (any(mod.results[geneID, (length(unique(DATA$treatment)) + 1):(length(unique(DATA$treatment)) * 2)] == 0)) {
                                      sapply(treatments,function(g) {
                                        c(max(fun_exp(xlim[2], unlist(mod.results[geneID, paste0("alpha_", g)])) - 0.2, 0))})
                                    } else {
                                      sapply(treatments,function(g) {
                                        c(max(dExp(xlim[2], unlist(mod.results[geneID, paste0(c("alpha_", "beta_"), g)])) - 0.2, 0))})
                                    }
                                  ),
                                  mapping = ggplot2::aes(label = text, x = x, y = y, color = treatment),
                                  # color = "black",
                                  vjust = -0.5,
                                  hjust = 1,
                                  size = (alphaSZ/1.3) * 1.5,
                                  alpha = 0.75)
    }
    if (any(what %in% "meanSE")) {
      p <- p + ggplot2::stat_summary(ggplot2::aes(x = t.decay,
                                                  y = value, color = treatment), fun = mean, geom = "line", #fun was fun.y
                                     size = if (any(what %in% "models")) {
                                       0.35} else {1}, alpha = if (any(what %in% "models")) {
                                         0.6}  else {1}) + ggplot2::stat_summary(ggplot2::aes(x = t.decay,
                                                                                              y = value, color = treatment), fun.data = ggplot2::mean_se,
                                                                                 geom = "errorbar", size = 0.35, alpha = 0.6)
    }
    if (any(what %in% "reps")) {
      p <- p + ggplot2::geom_point(ggplot2::aes(x = t.decay,
                                                y = value, color = treatment, shape = rep), size = 1.2,
                                   alpha = 0.5)
    }
    if (any(what %in% "models")) {
      if (any(mod.results[geneID, (length(unique(DATA$treatment)) +
                                   1):(length(unique(DATA$treatment)) * 2)] == 0)) {
        for (g in treatments) {
          p <- p + ggplot2::geom_line(data = data.frame(x = xlim[1]:xlim[2],
                                                        y = fun_exp(xlim[1]:xlim[2], a = unlist(mod.results[geneID,
                                                                                                            paste0("alpha_", gsub(" ", ".", g))]))),
                                      ggplot2::aes(x = x, y = y), color = colors[g],
                                      linewidth = 0.5, alpha = 1)
        }
      } else {
        for (g in treatments) {
          p <- p + ggplot2::geom_line(data = data.frame(x = xlim[1]:xlim[2],
                                                        y = dExp(xlim[1]:xlim[2], par = unlist(mod.results[geneID,
                                                                                                           paste0(c("alpha_", "beta_"), gsub(" ", ".",
                                                                                                                                             g))]))), ggplot2::aes(x = x, y = y), color = colors[g],
                                      linewidth = 0.5, alpha = 1)
        }
      }
    }
    p <- p + ggplot2::scale_color_manual("", breaks = names(colors),
                                         values = colors, labels = gsub("\\.", " ", treatments)) + ggplot2::coord_cartesian(xlim = xlim,
                                                                                                                            ylim = ylim) + ggplot2::scale_y_continuous(breaks = yticks) +
      ggplot2::scale_x_continuous(breaks = xticks) + ggplot2::ylab("relative abundance") +
      ggplot2::xlab("time (min)") + ggplot2::scale_shape(guide = "none")
    return(p)
  }
