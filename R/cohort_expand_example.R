library(tidyverse)
library(lubridate)
library(popEpi)
BL <- list(fot = seq(0, 5, by = 1/12))
x <- lexpand(sire, breaks = BL, 
             birth = bi_date, entry = dg_date, exit = ex_date,
             pophaz=popmort)



d <- tibble(sire)

d <- mutate(d, across(where(is.Date), decimal_date),
            id = row_number())

fotd <- tibble(fot = BL$fot)
perd <- tibble(per = decimal_date(ydm(BL$per)))

## by = character() to do full cross join
fotd <- full_join(d, fotd, by = character())
perd <- full_join(d, perd, by = character())


initper <- d |> mutate(age = dg_age,
                       per = dg_date,
                       fot = 0)

finper <- d |> mutate(fot = ex_date - dg_date,
                      per = ex_date,
                      age = dg_age + fot)

fotdp <- mutate(fotd, age = fot + dg_age,
       per = dg_date + fot) |>
  filter(per < ex_date, per > dg_date)

perdp <- mutate(perd, fot = per - dg_date,
                age = fot + dg_age) |>
  filter(per < ex_date & per > dg_date)


dd <- bind_rows(initper, perdp, fotdp, finper)


dd <- arrange(dd, id, per)
                
