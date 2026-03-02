#'
#' Refits a already generated plot to a new size
#'
#' @param plot The plot object to be refitted
#' @param owid The original width of the plot (in)
#' @param ohgt The original height of the plot (in)
#' @param scale The new width of the plot (in)
#'
#' @return The refitted plot object
#' @export
scale_plot <- function(plot, owid, ohgt, scale) {
    nwid = owid / scale
    nheight = nwid * ratio
    message("Refitting plot from ", owid, "x", ohgt, " to ", nwid, "x", nheight)
    base_size = 10 / scale**1.33
    message("New base size: ", base_size)
    return(plot & theme_create(base_size))
}