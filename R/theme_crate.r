#'
#' A publication-ready theme for ggplot2
#' 
#' @import ggplot2
#' @import ggthemes
#' 
#' @param base_size Base font size (default is 10)
#' 
#' @return A ggplot2 theme
#' @export
theme_create <- function(base_size=10) {
    return(
        ggthemes::theme_foundation(
            base_size = base_size,
            base_family = 'Helvetica'
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = 1.2*base_size,
                color = "#282A36",
                face = 'bold',
                hjust = 0,
                vjust = 0.5
            ),
            plot.background = ggplot2::element_blank(),
            plot.tag = ggplot2::element_text(
                size = 1.6*base_size,
                face = "bold",
                color = "#282A36"
            ),
            plot.margin = ggplot2::margin(0.04, 0.04, 0.04, 0.04, unit = "in"),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(color = "#282A36"),
            axis.line.x = NULL,
            axis.line.y = NULL,
            axis.text = ggplot2::element_text(size = ggplot2::rel(0.95), color = "#282A36"),
            axis.text.x = ggplot2::element_text(
                margin = ggplot2::margin(t = 0.8 * base_size/4),
                vjust = 1, color = "#282A36"
            ),
            axis.text.x.top = ggplot2::element_text(
                margin = ggplot2::margin(b = 0.8 * base_size/4),
                vjust = 0, color = "#282A36"
            ),
            axis.text.y = ggplot2::element_text(
                margin = ggplot2::margin(r = 0.5 * base_size/4),
                hjust = 1, color = "#282A36"
            ),
            axis.text.y.right = ggplot2::element_text(
                margin = ggplot2::margin(l = 0.5 * base_size/4),
                hjust = 0, color = "#282A36"
            ),
            axis.ticks = ggplot2::element_line(), 
            axis.ticks.length = ggplot2::unit(base_size/2.5, "pt"),
            axis.ticks.length.x = NULL,
            axis.ticks.length.x.top = NULL,
            axis.ticks.length.x.bottom = NULL,
            axis.ticks.length.y = NULL,
            axis.ticks.length.y.left = NULL,
            axis.ticks.length.y.right = NULL,
            strip.text = ggplot2::element_text(
                size = base_size, face = "bold", color = "#282A36"
            ),
            strip.background = ggplot2::element_blank(),
            strip.placement = "outside",
            legend.key.size= ggplot2::unit(base_size/40, "in"),
            legend.spacing = ggplot2::unit(0, "in"),
            legend.title = ggplot2::element_text(face="italic", color = "#282A36"),
            legend.text = ggplot2::element_text(face = "bold", color = "#282A36"),
            legend.justification = c(0.1, 0.75),
            legend.box.just = "right",
            legend.margin = ggplot2::margin(6, 6, 6, 6),
            legend.box.spacing = ggplot2::unit(-base_size/500, "in"),
            legend.position = "right",
            legend.background = ggplot2::element_blank()
        )
    )
}
