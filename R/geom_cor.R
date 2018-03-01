StatCor <- ggproto("StatCor", Stat,
                   compute_group = function(data, scales, method = "spearman") {
                       data = data[!is.na(data$x) & !is.na(data$y),]
                       cr = cor.test(data$x, data$y, method = method)
                       
                       data.frame(x = (max(data$x, na.rm = T) - min(data$x, na.rm = T)) / 2,
                                  y = max(data$y, na.rm = T),
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
#' @param method Method to calculate the correlation. Values are
#'   passed to [cor.test()]. (Spearman, Pearson, Kendall).
#' @examples
#' data(humanGender)
#' library(ggplot2)
#' ggplot(as.data.frame(assay(humanGender)[1:1000,]),
#'        aes(x = NA20502, y = NA20504)) +
#'   geom_point() +
#'   geom_cor(method = "kendall") 
#' @export
geom_cor <- function(mapping = NULL, data = NULL, geom = GeomCor,
                     method = "spearman",
                     position = "identity", na.rm = FALSE, show.legend = NA, 
                     inherit.aes = TRUE, ...) {
    layer(stat = StatCor,
                   data = data, mapping = mapping, geom = geom, 
                   position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                   params = list(na.rm = na.rm, method = method, ...)
    )
}
