
# edr

<!-- badges: start -->

[![R-CMD-check](https://github.com/NicChr/edr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/NicChr/edr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Estimated Dissemination Ratio

This package is a fast implementation of the Estimated Dissemination
Ratio, described here

[Estimated Dissemination Ratio—A Practical Alternative to the
Reproduction Number for Infectious
Diseases](https://dx.doi.org/10.3389/fpubh.2021.675065)

## Installation

``` r
# install.packages("devtools")
devtools::install_github("NicChr/edr")
```

### Libraries

``` r
library(edr)
library(data.table)
library(ggplot2)
```

### Recreating Figure 1

To recreate figure 1 from the paper, we use data from the John Hopkins
COVID-19 data repository.  
Source:
[CSSEGISandData/COVID-19](https://github.com/CSSEGISandData/COVID-19)

### Data

``` r
data(uk_covid_cases)
setDT(uk_covid_cases)
```

### EDR calculation

Here we calculate the 7-day EDR with `edr()` along with 99% percentile
confidence intervals

``` r
edr_seven_day <- edr(uk_covid_cases$new, window = 7, simulations = 1e04, alpha = 0.01)

# Join estimates and confint to data
uk_covid_cases <- cbind(uk_covid_cases, edr_seven_day)
```

We also calculate the 7-day rolling average of new confirmed cases

``` r
uk_covid_cases[, ma7 := frollmean(new, n = 7, align = "right")]
```

Finally plotting everything

``` r
uk_covid_cases <- uk_covid_cases[
  reporting_date >= as.Date("2020-03-10") & 
    reporting_date < as.Date("2020-05-01")]

scale_factor <- 500
# You could also calculate it more generally as below 
# scale_factor <- max(uk_covid_cases$ma7) / max(uk_covid_cases$est)
uk_lockdown <- as.Date("2020-03-24")

edr_plot <- uk_covid_cases |> 
  ggplot(aes(x = reporting_date, y = est)) + 
  geom_col(aes(y = ma7 / scale_factor), width = 1, alpha = 0.6, fill = "lightblue") +
  geom_smooth(aes(y = est, ymin = lower, ymax = upper), 
              stat = "identity", 
              linewidth = 1.25, 
              col = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(aes(x = uk_lockdown, y = 6, xend = uk_lockdown, yend = 3.25),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = uk_lockdown, y = 6, 
           label = "full UK lockdown \nimplemented",
           vjust = -0.1,
           fontface = "bold") +
  labs(x = "Reporting Date", y = "UK EDR") +
  scale_y_continuous(breaks = seq(0, 10, 1),
                     sec.axis = sec_axis(\(x) x * scale_factor,
                                         breaks = seq(0, 5000, 500),
                                         name = "UK new case 7-day rolling average")) +
  theme_bw() + 
  theme(axis.line.y.right = element_line(color = "lightblue"), 
        axis.ticks.y.right = element_line(color = "lightblue"),
        axis.text.y.right = element_text(color = "lightblue"), 
        axis.title.y.right = element_text(color = "lightblue")
  )

edr_plot
#> Warning in geom_segment(aes(x = uk_lockdown, y = 6, xend = uk_lockdown, : All aesthetics have length 1, but the data has 52 rows.
#> ℹ Please consider using `annotate()` or provide this layer with data containing
#>   a single row.
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />
