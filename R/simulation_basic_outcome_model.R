set.seed(34)
source('outcome_models.R')

library(tidyverse)
library(spatstat)
## swedishpines <- murchison$gold

## swedishpines <- shift(swedishpines, -c(352782.9, 6699742))


## it appears there's some sort of interpolation bug due to window scale going on that results in covariates not scaling properly.
## may also have something to do with distance calculations
myexposure <- rescale(swedishpines, c(10))
## myexposure <- redwood
## myexposure <- rescale(redwood, 0.67)

## myexposure <- swedishpines
exp_im <- (as.matrix(as.im(myexposure)))
kernx <- seq(from = -12, to = 12, length = 128 * 2)
kern <- fft2shift(outer(kernx, kernx, function(x,y) exp(- x ^2 - y^2)))

exp_imb <- matrix(0, nrow = 256, ncol = 256)
exp_imb[1:128, 1:128] <- exp_im


exp_im <- Re((fft(fft(exp_imb) * fft(kern), inverse = TRUE) / 256^2))[1:128, 1:128]
## intens_surface <- as.im(exp_im / 100, W = myexposure$window)
intens_surface <- as.im(exp_im * 10, W = myexposure$window)



## dd <- dd |> mutate(mytreatment = treatcov)
make_summary_plot <- function(model, intens, newpoints){

    
    plot_data <- mutate(newpoints,
                        pred = predict(model, newpoints, type = 'response'),
                        truth = interp.im(intens, newpoints),
                        res = truth - pred
                        ) 

    plot_data |>
        pivot_longer(c(truth, res, pred), names_to = 'plot', values_to = 'value') |>
        ggplot(aes(x, y, fill = value)) + geom_raster() + facet_wrap(~plot) + scale_fill_viridis_c()
}




mytreatment <- (intens_surface > median(intens_surface))



treated_intens_surface <- exp(log(intens_surface))
myoutcome <- (rpoispp(lambda = treated_intens_surface))
nd <- 200
q <- quadscheme(myoutcome, nd = nd)
mydummy <- as_tibble(coords(q$dummy))
mydd <- mydummy |>
    slice(-((nrow(mydummy) - 3):nrow(mydummy)))



fit2 <- ppmod(myoutcome, q, ppcov = list(myexposure = myexposure), k = 40)
mydd2 <- mutate(mydd, mytreatment = mytreatment[list(x = x, y = y)])
make_summary_plot(fit2, treated_intens_surface, mydd2)


show visual outcome model in basic model, including other covariates
# outcome model estimation




