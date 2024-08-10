
# library(tidyverse)
# cases <- read_csv("https://raw.github.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
# 
# uk_covid_cases <- cases |>
#   filter(`Country/Region` == "United Kingdom") |>
#   select(-(1:4)) |>
#   pivot_longer(everything(), values_to = "n") |>
#   mutate(date = mdy(name)) |>
#   summarise(n = sum(n), .by = date) |>
#   arrange(date) |>
#   mutate(new = n - lag(n, default = 0)) |>
#   mutate(new = as.integer(round(new))) |>
#   filter(date <= dmy(01042021)) |>
#   select(reporting_date = date, new) |>
#   as.data.frame()
# 
# usethis::use_data(uk_covid_cases, overwrite = TRUE)
