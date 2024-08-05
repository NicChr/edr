#' Estimated Dissemination Ratio
#'
#' @param x `[numeric]` - A vector of case numbers.
#' @param window `[integer(1)]` - An integer denoting the window size. \cr
#' Default is `1`.
#' @param na_rm `[logical(1)]` - Should `NA` values be ignored? \cr
#' Default is `FALSE`.
#' @param simulations `[numeric(1)]` - Number of Poisson simulations for
#' both the numerator and denominator. \cr 
#' If simulations is 0 (the default) then the result is
#' a numeric vector, otherwise a data frame of the estimate, 
#' lower confidence interval and upper confidence interval 
#' is returned.
#' @param alpha `[numeric(1)]` - Alpha significance level. \cr
#' The default is a 95% (1 - 0.05) confidence interval.
#'
#' @returns
#' A numeric vector of edr estimates when `simulations` is 0, otherwise
#' a 3-col `data.frame` of estimates, lower and upper confidence intervals.
#'
#' @details
#'
#' ### Ratio
#' EDR is a simple ratio calculated as the number of cases in a time interval, 
#' divided by the number of cases in a preceding time interval 
#' of the same width.
#' Therefore it is easy to produce, both conceptually and computationally, 
#' making it a useful metric in analysing epidemic curves.
#' 
#' For example, to calculate the 7-day EDR, the first ratio would be 
#' the number of cases in week 2 (days 8 to 14), 
#' divided by the number of cases in week 1 (days 1 to 7). 
#' You then apply this on a rolling basis, e.g. the number of cases between 
#' days 9 to 15 divided by the number of cases between days 2 to 8, and so on.
#' 
#' A ratio \bold{equal to} 1 signifies the cases have remained stable in 
#' the current time period relative to the previous one. \cr
#' A ratio \bold{above} 1 signifies the cases have increased in 
#' the current time period relative to the previous one.\cr
#' A ratio \bold{below} 1 signifies the cases have decreased in 
#' the current time period relative to the previous one. 
#' 
#' ### Confidence intervals
#' Confidence intervals are calculated by assuming the cases across both 
#' intervals are Poisson distributed random variables. 
#' Random deviates are then generated for each and ratios are calculated
#' across these two sets. Percentile confidence intervals are produced 
#' by taking the `alpha/2` and `1 - (alpha/2)` percentiles of the ratios.
#'
#' @references
#' Perez-Reche FJ, Taylor N, McGuigan C, Conaglen P, Forbes K, Strachan N, 
#' Honhold N (2021) Estimated Dissemination Ratio --- 
#' A practical alternative to the reproduction number for infectious diseases. 
#' Frontiers in Public Health 9. DOI:  10.3389/fpubh.2021.675065.
#' 
#' @examples
#' library(edr)
#' library(outbreaks)
#' library(data.table)
#' 
#' ebola <- as.data.table(ebola_sim_clean$linelist)
#' cases <- ebola[, .(n = .N), keyby = date_of_onset]
#' 
#' cases[, edr := edr(n, window = 7)]
#' 
#' cases <- cases[-(1:13)] # EDR starts at time = 2 * window
#' 
#' plot(cases$date_of_onset, cases$edr,
#'      xlab = "Date of onset",
#'      ylab = "EDR")
#' abline(h = 1, lty = 2, lwd = 2, col = "purple")
#' lines(cases$date_of_onset,
#'       fitted(loess(cases$edr ~ as.numeric(cases$date_of_onset))),
#'       lwd = 4, col = "blue")
#' title(main = "Cases above purple line are increasing
#' from week to week.\nCases below are decreasing.")
#' @export
edr <- function(x, window = 1, na_rm = FALSE, simulations = 0, alpha = 0.05){
  
  N <- length(x)
  
  # EDR
  
  top <- data.table::frollsum(x, n = window, align = "right", na.rm = na_rm)
  bottom <- data.table::frollsum(
    data.table::shift(x, n = window, type = "lag"), 
    n = window, align = "right", na.rm = na_rm
  )
  edr_est <- top / bottom
  
  # Start of confint loop
  
  start <- as.integer(window) * 2L
  end <- N
  confint_length <- max(N - start + 1L, 0L)
  
  if (simulations <= 0){
    return(edr_est)
  }
  
  check_alpha(alpha)
  lower_prob <- alpha / 2
  upper_prob <- 1 - lower_prob
  probs <- c(lower_prob, upper_prob)
  
  # Confidence interval calculation
  
  if (simulations > 0 && N >= start){

    lcl <- rep_len(NA_real_, N)
    ucl <- lcl

    # All simulations at once
    if ( (confint_length * simulations) <= 1e08 ){
      iss <- start:end
      g <- rep.int(seq_len(confint_length), simulations)
      edr_sim <- stats::rpois(simulations * confint_length, top[iss]) /
        stats::rpois(simulations * confint_length, bottom[iss])
      lcl[iss] <- collapse::fnth(edr_sim, lower_prob, g = g, use.g.names = FALSE)
      ucl[iss] <- collapse::fnth(edr_sim, upper_prob, g = g, use.g.names = FALSE)
    } else {
      for (i in start:end){
        edr_sim <- stats::rpois(simulations, top[i]) /
          stats::rpois(simulations, bottom[i])
        quantiles_sim <- collapse::fquantile(edr_sim, probs,
                                             type = 7L, names = FALSE)
        lcl[i] <- quantiles_sim[1L]
        ucl[i] <- quantiles_sim[2L]
      }
    }
  }
  data.frame(est = edr_est, lower = lcl, upper = ucl)
}

check_alpha <- function(x){
  stopifnot(is.numeric(x) && length(x) == 1 && x >= 0 && x <= 1)
}

