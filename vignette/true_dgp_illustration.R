library(tidyverse)
library(splines)
set.seed(42)
t <- seq(from = 0, to =  1, length.out = 5000)
x <-  ns(t, df = 15) %*% (rnorm(15))
y <-  ns(t, df = 7) %*% rexp(7)


dat <- tibble(t, x, y) |>
    mutate(
        roughx = x + rnorm(1000, 0, 0.05),
        roughy = y + rnorm(1000, 0, 0.05))

## Motivating plot: suppose that the "truth" that generates real-life data is a complicated mathematical structure (not even a function necessarily).
## in the motivating figure we think of the truth being the black curve and the noisy path around representing the "data".


ggplot(dat, aes(x,y)) +
    geom_path() +
    geom_path(aes(roughx, roughy), color = 'gray', alpha = 0.4) +
    geom_point(aes(roughx, roughy), size = 0.01, alpha = 0.1) + theme_void()


py <- predict(lm(y ~ ns(t, 8), data = dat))
px <- predict(lm(x ~ ns(t, 8), data = dat))


mutate(dat, py = py, px = px) |>
    ggplot( aes(x,y)) +
    geom_path() +
    geom_path(data = dat, aes(px, py), color = 'red', alpha = 0.4) +
    geom_point(aes(roughx, roughy), size = 0.01, alpha = 0.1) + theme_void()

