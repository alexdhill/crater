#'
#' A publication-ready theme for ggplot2
#'
#' @param base_size Base font size (default is 10)
#'
#' @return A ggplot2 theme
#'
#' @export
theme_crate <- function(base_size = 6) {
    return(
        ggthemes::theme_foundation(
            base_size = base_size,
            base_family = 'Helvetica'
        ) +
        ggplot2::theme(
            plot.tag = ggplot2::element_text(
                size = 1.5*base_size,
                face = "bold",
                color = "#282A36"
            ),
            plot.title = ggplot2::element_text(
                size = 1.25*base_size,
                color = "#282A36",
                face = 'bold',
                hjust = 0,
                vjust = 0.5
            ),
            plot.background = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(1, 1, 1, 1),
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(color = "#282A36"),
            axis.line.x = NULL,
            axis.line.y = NULL,
            axis.text = ggplot2::element_text(
                size = base_size,
                color = "#282A36"
            ),
            axis.text.x = ggplot2::element_text(
                margin = ggplot2::margin(t = 1),
                vjust = 1, color = "#282A36"
            ),
            axis.text.x.top = ggplot2::element_text(
                margin = ggplot2::margin(b = 1),
                vjust = 0, color = "#282A36"
            ),
            axis.text.y = ggplot2::element_text(
                margin = ggplot2::margin(r = 1),
                hjust = 1, color = "#282A36"
            ),
            axis.text.y.right = ggplot2::element_text(
                margin = ggplot2::margin(l = 1),
                hjust = 0, color = "#282A36"
            ),
            axis.ticks = ggplot2::element_line(),
            axis.ticks.length = ggplot2::unit(base_size/2, "pt"),
            axis.ticks.length.x = NULL,
            axis.ticks.length.x.top = NULL,
            axis.ticks.length.x.bottom = NULL,
            axis.ticks.length.y = NULL,
            axis.ticks.length.y.left = NULL,
            axis.ticks.length.y.right = NULL,
            strip.text = ggplot2::element_text(
                size = base_size,
                face = "bold",
                color = "#282A36"
            ),
            strip.background = ggplot2::element_blank(),
            strip.placement = "outside",
            legend.key.size = ggplot2::unit(base_size, "pt"),
            legend.spacing = ggplot2::unit(0, "cm"),
            legend.title = ggplot2::element_text(
                size = base_size,
                face = "bold",
                color = "#282A36"
            ),
            legend.text = ggplot2::element_text(
                size = base_size,
                face = "italic",
                color = "#282A36"
            ),
            legend.background = ggplot2::element_blank(),
            legend.justification = c(0, 0.75),
            legend.box.just = "right",
            legend.margin = ggplot2::margin(1, 1, 1, 1),
            legend.box.spacing = ggplot2::unit(0, "cm"),
            legend.position = "right"
        )
    )
}
