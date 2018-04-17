
#' a custom ggplot2 theme
#'
#' A custom ggplot2 theme generating function for ggplot2 plots; can be further
#'     manipulated using standard ggplot2 syntax.
#'
#' @param bigFont larger font size of axis labels in points (used for plot
#'     title, axis titles, facet titles)
#' @param smFont fractional mulitiplier of \code{bigFont} (used for axis text)
#' @param x.ang x-axis label angle
#' @param leg.pos legend position on plot as relative coordinates c(x,y) (i.e.,
#'     range is [0,1]) or 'right', 'left', 'above', 'below'
#'
#' @return returns a ggplot2 theme of class "theme" "gg"
#'
#' @export
#'
#' @examples plain_theme(10)


plain_theme <- function(bigFont = 30,
                   smFont = 0.85,
                   x.ang = 0,
                   leg.pos = c(0.85, 0.85)) {
  ggplot2::theme(
    #text=element_text(family=c("Courier")),
    panel.background = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(size = bigFont, vjust = 0.5),
    #element_rect(color="white",fill="white",size = 3,linetype="solid"),
    panel.grid.minor = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(
      colour = 'black',
      size = 1,
      linetype = 'solid'
    ),
    axis.line.y = ggplot2::element_line(
      colour = 'black',
      size = 1,
      linetype = 'solid'
    ),
    plot.margin = grid::unit(c(1, 1, 1, 1), "lines"),
    axis.title.x = ggplot2::element_text(size = bigFont),
    axis.text.x = ggplot2::element_text(
      colour = "black",
      size = bigFont * smFont,
      angle = x.ang,
      vjust = 0.5
    ),
    axis.title.y = ggplot2::element_text(size = bigFont, angle = 90),
    axis.text.y = ggplot2::element_text(colour = "black", size = bigFont *
                                          smFont),
    legend.box.background = ggplot2::element_blank(),
    #colour="white",fill = 'none'),element_blank(),
    legend.key = ggplot2::element_rect(colour = "white", fill = 'white'),
    legend.position = leg.pos,
    #c(0.85,0.85),
    legend.text = ggplot2::element_text(size = bigFont, face = c("italic")),
    legend.title = ggplot2::element_text(size = bigFont),
    #face=c("italic")),
    legend.key.height = grid::unit(bigFont / 12, "lines"),
    #,
    legend.title.align = 0.1,
    strip.text.x = ggplot2::element_text(
      size = bigFont,
      colour = "black",
      face = 'bold'
    ),
    strip.text.y = ggplot2::element_text(
      size = bigFont,
      colour = "black",
      face = 'bold',
      angle = 270
    )
    #aspect.ratio=1
  )
}
