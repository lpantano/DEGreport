StatCor <- ggproto("StatCor", Stat,
                   compute_group = function(data, scales,
                                            xpos = NULL, ypos = NULL,
                                            method = "spearman") {
                       
                       data = data[!is.na(data$x) & !is.na(data$y),]
                       cr = cor.test(data$x, data$y, method = method)
                       x = xpos
                       if (is.null(x))
                           x = min(data$x) + abs((max(data$x) - min(data$x)) / 2)
                       y = ypos
                       if (is.null(y))
                           y = scales$y$range$range[2] + 0.1 * scales$y$range$range[2]
                       data.frame(x = x,
                                  y = y,
                                  label = paste0("Cor (", method, "): ",
                                                 round(cr[["estimate"]], digits = 2),
                                                 ", pval: ",
                                                 round(cr[["p.value"]], digits = 2)))
                       
                   },
                   required_aes = c("x", "y")
)

GeomCor <- ggproto("GeomText", Geom,
                   required_aes = c("x", "y"),
                   
                   default_aes = aes(
                       colour = "black", size = 3.88, angle = 0, hjust = 0.5,
                       vjust = 0.5, alpha = 1, family = "", fontface = 1, lineheight = 1.2
                   ),
                   
                   
                   draw_panel = function(data, panel_params, coord, parse = FALSE,
                                         na.rm = FALSE, check_overlap = FALSE) {
                       
                       lab <- data$label
                       if (parse) {
                           lab <- parse(text = as.character(lab))
                       }
                       data <- coord$transform(data, panel_params)
                       if (is.character(data$vjust)) {
                           data$vjust <- data$y
                       }
                       if (is.character(data$hjust)) {
                           data$hjust <- data$x
                       }

                       textGrob(
                           lab,
                           data$x, data$y, default.units = "native",
                           hjust = data$hjust, vjust = data$vjust,
                           gp = grid::gpar(
                               col = data$colour,
                               fontsize = data$size * 3,
                               fontfamily = data$family,
                               fontface = data$fontface,
                               lineheight = data$lineheight
                           )
                       )
                   },
                   
                   draw_key = draw_key_text
)

#' Add correlation and p-value to a [ggplot2] plot
#' 
#' `geom_cor` will add the correlatin, method and p-value to the plot
#'   automatically guessing the position if nothing else specidfied.
#'   family font, size and colour can be used to change the format.
#' 
#' @seealso [ggplot2::layer()]
#' @param mapping Set of aesthetic mappings created by [aes()] or
#'   [aes_()]. If specified and `inherit.aes = TRUE` (the
#'   default), it is combined with the default mapping at the top level of the
#'   plot. You must supply `mapping` if there is no plot mapping.
#' @param data The data to be displayed in this layer. There are three
#'    options:
#'
#'    If `NULL`, the default, the data is inherited from the plot
#'    data as specified in the call to [ggplot()].
#'
#'    A `data.frame`, or other object, will override the plot
#'    data. All objects will be fortified to produce a data frame. See
#'    [fortify()] for which variables will be created.
#'
#'    A `function` will be called with a single argument,
#'    the plot data. The return value must be a `data.frame.`, and
#'    will be used as the layer data.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics,
#'   rather than combining with them. This is most useful for helper functions
#'   that define both data and aesthetics and shouldn't inherit behaviour from
#'   the default plot specification, e.g. [borders()].
#' @param method Method to calculate the correlation. Values are
#'   passed to [cor.test()]. (Spearman, Pearson, Kendall).
#' @param ypos Locate text at that position on the y axis.
#' @param xpos Locate text at that position on the x axis.
#' @param ... other arguments passed on to [layer()]. These are
#'   often aesthetics, used to set an aesthetic to a fixed value, like
#'   `color = "red"` or `size = 3`. They may also be parameters
#'   to the paired geom/stat.
#' @details It was integrated after reading this tutorial to extend
#' ggplot2 [layers](https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html)
#' @examples
#' data(humanGender)
#' library(SummarizedExperiment)
#' library(ggplot2)
#' ggplot(as.data.frame(assay(humanGender)[1:1000,]),
#'        aes(x = NA20502, y = NA20504)) +
#'   geom_point() +
#'   ylim(0,1.1e5) +
#'   geom_cor(method = "kendall", ypos = 1e5) 
#' @export
geom_cor <- function(mapping = NULL, data = NULL,
                     method = "spearman",
                     inherit.aes = TRUE, ...) {
    layer(stat = StatCor,
                   data = data, mapping = mapping, geom = GeomCor,
                   position = "identity",
                   inherit.aes = inherit.aes,
                   params = list(na.rm = FALSE, method = method, ...)
    )
}
