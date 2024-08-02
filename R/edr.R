
#' Estimated Dissemination Ratio
#'
#' @param x `[numeric]` - A vector of case numbers.
#' @param window `[integer(1)]` - An integer denoting the window size.
#' @param na_rm `[logical(1)]` - Should `NA` values be ignored? Default is `FALSE`.
#' @param simulations `[numeric(1)]` - Number of poisson simulations for
#' both the numerator and denominator. If simulations is 0 then the result is
#' a numeric vector, otherwise a matrix of the estimate, lower confidence interval
#' and upper confidence interval are returnd.
#' @param alpha `[numeric(1)]` - Alpha significance level.
#' The default is a 95% (1 - 0.05) confidence interval.
#'
#' @returns
#' A numeric vector of edr estimates when `simulations` is 0, otherwise
#' a 3-col matrix of estimates, lower and upper confidence intervals.
#'
#' @details
#'
#' ### Confidence intervals
#' To generate confidence intervals, simply specify
#' a `simulations` value greater than 0.
#'
#' @references
#' Perez-Reche FJ, Taylor N, McGuigan C, Conaglen P, Forbes K, Strachan N, Honhold N (2021) Estimated Dissemination Ratio --- A practical alternative to the reproduction number for infectious diseases. Frontiers in Public Health 9. DOI:  10.3389/fpubh.2021.675065.
#'
#'
#' @export
edr <- function(x, window = 1, na_rm = FALSE, simulations = 0, alpha = 0.05){

  N <- length(x)
  top <- data.table::frollsum(x, n = window, align = "right", na.rm = na_rm)
  bottom <- data.table::frollsum(data.table::shift(x, n = window, type = "lag"),
                                 n = window, align = "right", na.rm = na_rm)
  edr_est <- top / bottom

  # Start of confint loop
  start <- window * 2
  end <- N
  confint_length <- max(N - start + 1, 0)

  # edr is valid only on complete windows
  edr_est[seq_len(max(start - 1L, 0L))] <- NA

  if (simulations <= 0){
    return(edr_est)
  }

  check_alpha(alpha)
  lower_prob <- alpha / 2
  upper_prob <- 1 - lower_prob
  probs <- c(lower_prob, upper_prob)

  out <- matrix(c(edr_est, rep_len(NA_real_, N * 2)), ncol = 3,
                byrow = FALSE)

  if (simulations > 0 && N >= start){
    # If we have relatively little data compared to sims we loop through x
    # Otherwise we loop through sims
    if (N <= (simulations * 10)){
      for (i in start:end){
        edr_sim <- stats::rpois(simulations, top[i]) / stats::rpois(simulations, bottom[i])
        quantiles_sim <- collapse::fquantile(edr_sim, probs, type = 7L,
                                             names = FALSE)
        out[i, 2] <- quantiles_sim[1]
        out[i, 3] <- quantiles_sim[2]
      }
    } else {
      sims <- matrix(numeric(confint_length * simulations), ncol = simulations)
      tops <- top[start:end]
      bottoms <- bottom[start:end]

      for (i in seq_len(simulations)){
        sims[, i] <- stats::rpois(confint_length, tops) / stats::rpois(confint_length, bottoms)
      }
      j <- 1
      for (i in start:end){
        quantile_sim <- collapse::fquantile(sims[j, ], probs, type = 7L, names = FALSE)
        out[i, 2] <- quantile_sim[1]
        out[i, 3] <- quantile_sim[2]
        j <- j + 1L
      }
    }

  }
  colnames(out) <- c("est", "lower", "upper")
  out
}

check_alpha <- function(x){
  stopifnot(is.numeric(x) && length(x) == 1 && x >= 0 && x <= 1)
}

