#' Plot EDR
#' 
#' @description
#' A plot method for an object of class `edr` which is returned
#' by `edr()`. 
#' 
#' @param x `[numeric]` - An object of class "edr", 
#' created using the `edr()` function.
#' @param include_cases `[logical(1)]` Should case counts be plotted alongside the EDR 
#' estimates? The default is `TRUE`.
#' @param ... Unused.
#' 
#' @returns
#' A `ggplot` of the EDR trend, optionally alongside case counts.
#' 
#' @seealso [edr()]
#' 
#' @importFrom ggplot2 .data
#' 
#' @exportS3Method base::plot
#' @export plot.edr
plot.edr <- function(x, include_cases = TRUE, ...){
  check_valid_edr(x)
  
  plot_data <- x
  
  indigo <- "#0077B6"
  
  max_cases <- max(plot_data[["cases"]], na.rm = TRUE)
  max_edr <- max(plot_data[["edr"]][is.finite(plot_data[["edr"]])], na.rm = TRUE)
  min_edr <- min(plot_data[["edr"]], na.rm = TRUE)
  scale_factor <- max_cases / max_edr
  
  # median_cases <- stats::median(plot_data[["cases"]], na.rm = TRUE)
  # median_edr <- stats::median(plot_data[["edr"]][is.finite(plot_data[["edr"]])], na.rm = TRUE)
  # scale_factor <- median_cases / median_edr
  
  # mean_cases <- mean(plot_data[["cases"]], na.rm = TRUE)
  # mean_edr <- mean(plot_data[["edr"]][is.finite(plot_data[["edr"]])], na.rm = TRUE)
  # scale_factor <- mean_cases / mean_edr
  
  edr_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["time"]], 
                                                      y = .data[["edr"]])) +
    ggplot2::geom_line() +
    # ggplot2::coord_cartesian(ylim = c(0, max_edr)) +
    ggplot2::labs(x = "", y = "EDR")
  if (!is.null(plot_data[["lower"]])){
    edr_plot <- edr_plot +
      ggplot2::geom_smooth(ggplot2::aes(ymin = .data[["lower"]], ymax = .data[["upper"]]), 
                           linetype = 0, stat = "identity")
  }
  edr_plot <- edr_plot +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(
      breaks = scales::extended_breaks(n = 10),
      # transform = if (!include_cases) "log2" else "identity",
      sec.axis = if (include_cases){
        ggplot2::sec_axis(function(x) x * scale_factor,
                          name = "Cases",
                          breaks = scales::extended_breaks(n = 10))
      } else {
        ggplot2::waiver()
      }
    )
  if (include_cases){
    
    # Smallest difference between time points
    # This will determine geom_col bar width
    
    t <- as.double(plot_data[["time"]])
    d <- diff2(t)
    bar_width <- min(d[which(d > sqrt(.Machine$double.eps))])
    if (length(bar_width) == 0 || bar_width == 0 || is.infinite(bar_width) || is.na(bar_width)){
      bar_width <- 1
    }
    
    edr_plot <- edr_plot +
      ggplot2::geom_col(ggplot2::aes(y = .data[["cases"]] / scale_factor), 
                        # width = 1, 
                        width = bar_width,
                        alpha = 0.4, fill = indigo) +
      ggplot2::theme(
        axis.line.y.right = ggplot2::element_line(color = indigo),
        axis.ticks.y.right = ggplot2::element_line(color = indigo),
        axis.text.y.right = ggplot2::element_text(color = indigo),
        axis.title.y.right = ggplot2::element_text(color = indigo)
      )
  }
  
  edr_plot
}

