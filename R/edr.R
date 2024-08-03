#' Estimated Dissemination Ratio
#'
#' @param x `[numeric]` - A vector of case numbers.
#' @param window `[integer(1)]` - An integer denoting the window size.
#' @param na_rm `[logical(1)]` - Should `NA` values be ignored? 
#' Default is `FALSE`.
#' @param simulations `[numeric(1)]` - Number of Poisson simulations for
#' both the numerator and denominator. \cr 
#' If simulations is 0 then the result is
#' a numeric vector, otherwise a matrix of the estimate, 
#' lower confidence interval and upper confidence interval 
#' are returned.
#' @param alpha `[numeric(1)]` - Alpha significance level.
#' The default is a 95% (1 - 0.05) confidence interval.
#'
#' @returns
#' A numeric vector of edr estimates when `simulations` is 0, otherwise
#' a 3-col matrix of estimates, lower and upper confidence intervals.
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
#' A ratio above 1 signifies the cases have increased in 
#' the current time period relative to the previous one. 
#' A ratio below 1 signifies the cases have decreased in 
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
  
  start <- window * 2
  end <- N
  confint_length <- max(N - start + 1, 0)
  
  if (simulations <= 0){
    return(edr_est)
  }
  
  check_alpha(alpha)
  lower_prob <- alpha / 2
  upper_prob <- 1 - lower_prob
  probs <- c(lower_prob, upper_prob)
  
  out <- matrix(c(edr_est, rep_len(NA_real_, N * 2)), ncol = 3,
                byrow = FALSE)
  
  # Confidence interval calculation
  
  if (simulations > 0 && N >= start){
    
    # If we have relatively little data compared to sims we loop through x
    # Otherwise we loop through sims
    
    if (N <= (simulations * 10)){
      for (i in start:end){
        edr_sim <- stats::rpois(simulations, top[i]) / 
          stats::rpois(simulations, bottom[i])
        quantiles_sim <- collapse::fquantile(edr_sim, probs, 
                                             type = 7L, names = FALSE)
        out[i, 2] <- quantiles_sim[1]
        out[i, 3] <- quantiles_sim[2]
      }
    } else {
      sims <- matrix(numeric(confint_length * simulations), ncol = simulations)
      tops <- top[start:end]
      bottoms <- bottom[start:end]
      
      for (i in seq_len(simulations)){
        sims[, i] <- stats::rpois(confint_length, tops) / 
          stats::rpois(confint_length, bottoms)
      }
      j <- 1
      for (i in start:end){
        quantile_sim <- collapse::fquantile(sims[j, ], probs, 
                                            type = 7L, names = FALSE)
        out[i, 2] <- quantile_sim[1]
        out[i, 3] <- quantile_sim[2]
        j <- j + 1L
      }
    }
    
  }
  colnames(out) <- c("est", "lower", "upper")
  out
}

# edr2 <- function(x, window = 1, na_rm = FALSE, simulations = 0, alpha = 0.05){
# 
#   N <- length(x)
#   top <- data.table::frollsum(x, n = window, align = "right", na.rm = na_rm)
#   bottom <- data.table::frollsum(data.table::shift(x, n = window, type = "lag"),
#                                  n = window, align = "right", na.rm = na_rm)
#   edr_est <- top / bottom
# 
#   # Start of confint loop
#   start <- window * 2
#   end <- N
#   confint_length <- max(N - start + 1, 0)
# 
#   # edr is valid only on complete windows
#   edr_est[seq_len(max(start - 1L, 0L))] <- NA
# 
#   if (simulations <= 0){
#     return(edr_est)
#   }
# 
#   check_alpha(alpha)
#   lower_prob <- alpha / 2
#   upper_prob <- 1 - lower_prob
#   probs <- c(lower_prob, upper_prob)
# 
#   out <- matrix(c(edr_est, rep_len(NA_real_, N * 2)), ncol = 3,
#                 byrow = FALSE)
# 
#   if (simulations > 0 && N >= start){
# 
#     # top <- as.integer(top)
#     # bottom <- as.integer(bottom)
# 
#     top_sim <- stats::rpois(simulations, top[start])
#     bottom_sim <- stats::rpois(simulations, bottom[start])
# 
#     quantiles_sim <- collapse::fquantile(top_sim / bottom_sim, probs, type = 7L,
#                                          names = FALSE)
#     out[start, 2L] <- quantiles_sim[1L]
#     out[start, 3L] <- quantiles_sim[2L]
# 
#     if (N > start){
#       for (i in (start + 1L):end){
#         
#         # Instead of generating simulated variables for every 
#         # numerator/denominator pair, we simply do it once for the first pair
#         # And then update those values so the mean and variance match 
#         # each new respective numerator/denominator pair
#         
#         top_sim <- ( (top_sim - top[i - 1L]) * sqrt(top[i] /  top[i - 1L]) ) + top[i]
#         bottom_sim <- ( (bottom_sim - bottom[i - 1L]) * sqrt(bottom[i] /  bottom[i - 1L]) ) + bottom[i]
#         
# 
#         edr_sim <- top_sim / bottom_sim
#         quantiles_sim <- collapse::fquantile(edr_sim, probs, type = 7L,
#                                              names = FALSE)
#         out[i, 2L] <- quantiles_sim[1L]
#         out[i, 3L] <- quantiles_sim[2L]
#       }
#     }
# 
#   }
#   colnames(out) <- c("est", "lower", "upper")
#   out
# }

check_alpha <- function(x){
  stopifnot(is.numeric(x) && length(x) == 1 && x >= 0 && x <= 1)
}

