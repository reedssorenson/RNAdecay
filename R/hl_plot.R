#' hl_plot() function
#'
#' Plots RNA half-life distribution with select half-lives of select RNAs as large arrows colored by treatment using the ggplot2 package.
#'



#' @param geneID single gene ID from data set (e.g. "AT3G14100") for which to plot data/model
#' @param gene_symbol (optional) pasted to gene ID in plot label (e.g., "AT3G14100/UBP1C)
#' @param df_decay_rates data.frame of modeling results with decay rate columns labeled as alpha_<treatment>
#' @param hl_dist_treatment name of the treatment for which the background distribution will be plotted
#' @param hl_treatment names of the treatments for which arrows indicating half-life will be plotted
#' @param arrow_colors (optional) character vector of R colors; named with corresponding treatments
#' @param arrow_lab_loc label arrows on plot ("plot") or in a key ("key")
#' @param x_limits x-axis (half-life) limits in min; default is log(2)/c(0.25,4.5e-4)
#' @param x_breaks x-axis (half-life) breaks/tick marks in min defaults to c(5,1:12*10,180,240,300,360,420,480,720,1080,1440)
#' @param x_tick_labels x-axis (half-life) break labels, defaults to c("5","10","","30","","","60","","","","","","2h","","4h","","","","8h","12h","","24h")+
#'


#'
#' @return returns a ggplot to be used with print; could also be modified using the syntax of ggplot2 e.g.'+geom_XXXX(...)'
#'
#' @export
#'
#' @examples
#'
#' p <- hl_plot(
#' geneID = rownames(RNAdecay::results)[4],
#' df_decay_rates = RNAdecay::results,
#' hl_treatment = c("WT","sov","vcs","vcs.sov"),
#' hl_dist_treatment = "WT",
#' arrow_colors = c(WT = "#88CCEE", sov = "#CC6677", vcs = "#117733",vcs.sov = "#882255"),
#' arrow_lab_loc = "key",
#' gene_symbol = ""
#' )
#'
#'print(p)
#'
#' p <- hl_plot(
#' geneID = rownames(RNAdecay::results)[4],
#' gene_symbol = "",
#' df_decay_rates = RNAdecay::results,
#' hl_dist_treatment = "WT",
#' hl_treatment = c("WT","sov","vcs","vcs.sov"),
#' arrow_colors = c(WT = "#88CCEE", sov = "#CC6677", vcs = "#117733",vcs.sov = "#882255"),
#' arrow_lab_loc = "plot"
#' )
#'
#'print(p)



hl_plot <-
  function(
    geneID,
    gene_symbol = "",
    df_decay_rates,
    hl_dist_treatment,
    hl_treatment,
    arrow_colors = NA,
    arrow_lab_loc = c("key"),
    x_limits = log(2)/c(0.25,4.5e-4),
    x_breaks = c(5,1:12*10,180,240,300,360,420,480,720,1080,1440),
    x_tick_labels = c("5","10","","30","","","60","","","","","","2h","","4h","","","","8h","12h","","24h")

  ) {

    treatment.cols <- paste0("alpha_",hl_treatment)
    names(treatment.cols) <- hl_treatment
    hls = round(log(2)/ df_decay_rates[geneID, treatment.cols],1)

    if(all(is.na(arrow_colors))) {

      arrow_colors <- c("#88CCEE","#CC6677","#DDCC77","#117733",
                        "#332288","#AA4499","#44AA99","999933",
                        "#882255","#661100","#6699CC")[seq_along(hl_treatment)]
      names(arrow_colors) = hl_treatment

    }

    if(!hl_dist_treatment %in% hl_treatment) {stop("'hl_dist_treatment' must be one of 'hl_treatment'.")}

    if(length(geneID)>1) {

      warning("'length(geneID)' > 1, only the first will be used.")
      geneID <- geneID[1]

      }

    if(all(is.na(hls[geneID,]))){

      stop("'geneID' is not found in 'df_decay_rates'.")

      }

    if(length(hl_treatment) != length(arrow_colors)) {

      stop("'length(hl_treatment)' does not equal 'length(arrow_colors'. These must be the same.")

      }

    if(is.null(names(arrow_colors)) | !all(names(arrow_colors) %in% hl_treatment)) {

      warning("Not all colors' names match hL_treatment or names are missing. Colors will be arbitrarily assigned.")
      names(arrow_colors) <- hl_treatment

      }

    labs <- sapply(hl_treatment, function(x){

      hl <- hls[geneID,treatment.cols[x]]

      if(hl > max(x_breaks)){

        paste0("paste(\"",x," \",italic(t)[\"1/2\"] > \"24 h\")")

      } else{

        paste0("paste(\"",x," \",italic(t)[\"1/2\"] == \"",

               if(hl < 120)    {paste0(as.character(round(hl))," min")}    else    {paste0(round(hl/60,1)," h")},

               "\")")}

    })

    hl_dist <- df_decay_rates[df_decay_rates[,treatment.cols[hl_dist_treatment]] > 4.8e-5,
                              treatment.cols[hl_dist_treatment]]
    quants <- round(stats::quantile(log(2)/hl_dist))
    peak <- max(table(cut(log10(log(2)/hl_dist),include.lowest = TRUE,
                          breaks = c(seq(from =log10(x_limits[1]),
                                         to = x.max<-if(max(x_breaks) < max(log(2)/hl_dist)) {
                                           log10(max(x_limits))}    else    {max(log10(log(2)/hl_dist))},
                                         by = 0.02)))))

    arrow_lengths <- peak*0.9^seq_along(hl_treatment)
    names(arrow_lengths) <- hl_treatment
    arrow_label_heights <- arrow_lengths + arrow_lengths*0.1
    names(arrow_label_heights) <- hl_treatment
    arrow_thickness <- 3*0.85^seq_along(hl_treatment)
    names(arrow_thickness) <- hl_treatment
    arrowhead_length_in <- grDevices::dev.size(units = "in")[2]*0.04*0.85^seq_along(hl_treatment)
    names(arrowhead_length_in) <- hl_treatment

    p <- ggplot2::ggplot(df_decay_rates[df_decay_rates[,treatment.cols[hl_dist_treatment]] > 4.8e-5,])

    # add histogram and x scale
    p <- p + ggplot2::geom_histogram(data = data.frame(hl_dist = hl_dist),
                                     mapping = ggplot2::aes(x = log(2)/hl_dist),
                                     binwidth=0.02,
                                     fill = grDevices::gray(0.8))+
      ggplot2::geom_vline(xintercept=quants[2],color = "white",linetype = "dashed",size = 0.25)+
      ggplot2::geom_vline(xintercept=quants[3],color = "white",linetype = "dashed",size = 0.25)+
      ggplot2::geom_vline(xintercept=quants[4],color = "white",linetype = "dashed",size = 0.25)+
      ggplot2::scale_x_log10(limits=x_limits,
                             breaks = x_breaks,
                             labels = x_tick_labels)+
      ggplot2::annotate(geom = "text",
                        x = max(quants[3]),
                        y = 0,
                        label = paste0(hl_dist_treatment," distribution"),
                        hjust = 0.575,
                        vjust = -0.2,
                        color = grDevices::gray(0.65),
                        size = 72.27*grDevices::dev.size(units = "in")[2]/60)

    if(arrow_lab_loc == "plot"){

      for(i in hl_treatment){

        # add arrow
        hl <- log(2)/df_decay_rates[geneID,treatment.cols[i]]
        p <- p+ggplot2::geom_segment(data = data.frame(x=if(hl > x_limits[2]) {x_limits[2]} else {hl},
                                                       y=arrow_lengths[i]),#arrow height
                                     mapping = ggplot2::aes(x = x, y=y,xend = x, yend = 0),
                                     arrow = ggplot2::arrow(angle = 18,
                                                         length = grid::unit(arrowhead_length_in,"in"),
                                                         type = "closed"),#arrow head
                                     color = arrow_colors[i],
                                     size = arrow_thickness[i],#arrow width in pts
                                     lineend = "butt",
                                     linejoin = "mitre")

        # add arrow label
        p <- p+ggplot2::annotate(parse = TRUE,geom = "text",
                                 x = if(hl>max(x_limits[2])){max(x_limits[2])}else{hl},
                                 y = arrow_label_heights[i], # arrow lebel height
                                 label = labs[i],
                                 hjust = if(hl>103.4){0.9} else {0.15},
                                 vjust = 0.5,
                                 color = arrow_colors[i],
                                 size = 72.27*grDevices::dev.size(units = "in")[2]/75)
      }
    }

    if(arrow_lab_loc == "key"){

      hl <- log(2)/df_decay_rates[geneID,treatment.cols]

      p <- p + ggplot2::geom_segment(data = data.frame(x=unlist(hl),
                                                       y=arrow_lengths[hl_treatment], #arrow height
                                                       labs = hl_treatment),
                                     mapping = ggplot2::aes(x = x, y=y,xend = x, yend = 0,
                                                            color = labs,
                                                            size = labs),
                                     arrow = ggplot2::arrow(angle = 18,
                                                         length = grid::unit(arrowhead_length_in,"in"),
                                                         type = "closed"), #arrow head
                                     lineend = "butt",
                                     linejoin = "mitre")+
        ggplot2::scale_color_manual(NULL,values = arrow_colors,breaks = hl_treatment,labels=scales::parse_format()(labs[hl_treatment]))+
        ggplot2::scale_size_manual(values = arrow_thickness)+
        ggplot2::guides(size=FALSE)
    }

    p <- p + ggplot2::annotate(geom="text",
                               label = if(gene_symbol != "") {paste0(geneID,"/",gene_symbol)} else {geneID},
                               x= 10^x.max,
                               y = peak*1.25,
                               hjust = 1,
                               vjust = 1,
                               color = "black",
                               size = 72.27*grDevices::dev.size(units = "in")[2]/50)+
      ggplot2::xlab(expression(italic(t)["1/2"]))+
      ggplot2::coord_cartesian(xlim = x_limits,
                               ylim = c(0,peak*1.25))+
      ggplot2::theme_classic()+
      ggplot2::theme(title = ggplot2::element_text(color = "darkblue"),
                     axis.text = ggplot2::element_text(color = grDevices::gray(0)),
                     axis.title = ggplot2::element_text(color = grDevices::gray(0))
      )
  }
