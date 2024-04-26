library(RandomFields)
library(tidyverse)
library(sf)
library(spatstat)
source('covariate_types.R')
source('Fcov.R')
source('Pcov.R')
source('outcome_models.R')
source('exposure_models_continuous.R')
source('coord.R')
source('treatment_effects.R')

#####################
## BASIC DATA PREP ##
#####################


## for all simulations we will use the basic map of allegheny county
allegheny <- read_rds('simulation_utilities/allegheny_region_tracts.rds') |>
    mutate(region = id, id = row_number()) |>
    st_make_valid() |>
    st_transform('+proj=aeqd +lon_0=-79.3212891 +lat_0=40.2041883 +datum=WGS84 +units=m +no_def')

covering_grid <- st_make_grid(allegheny, n = c(20, 20))
alltess <-as.function(tess(tiles = map(allegheny$geometry, as.owin)))
window <- as.owin(allegheny)

set.seed(12)
## generate point process valued outcome with confounding by covariates for basic binary exposure.
## for simulation 1 we consider a simple, spatially located binary treatment with a non-confounding covariate

treat_point <- runifpoint(1, window)
treatment1 <- function(x,y) distfun(treat_point)(x,y) < 12000

covar_values <- rnorm(nrow(allegheny))
covarfun <-function(x,y) covar_values[as.numeric(alltess(x,y))]


dat <- data.frame(treatment = Fcov(treatment1), tess = Fcov(covarfun))


## the intercept term is specified as the expected number of points with treatment and covariate set to 0
## to get a reasonable and roughly controllable process it is best to set this term manually as log(N /area(w))
## where N is the expected number of points for the "neutral" process and area is the area of the observation window


#######################
## BINARY TREATMENTS ##
#######################

intercept <- log(100 / area(window))
intens <- function(x,y) exp(evaluate(with(dat, treatment * 2  + intercept), coord(x,y)))


intens_im <- as.im(intens, window)
pp <- rpoispp(intens_im)

Q <- quadscheme(pp)
Qc <- coords(Q)
Qc <- coord(Qc[,1], Qc[,2])

covariates <- map_dfc(dat, ~setup(., Qc))
outcome <- prep_outcome(pp, Q)
for_mod <- bind_cols(outcome, covariates)

fit <- gam(outcome ~ treatment + tess , weights = w,  family =poisson(), data = for_mod)
muf <- outcome_model(fit)

exposure_fit <- gam(treatment ~ tess + s(x,y), data = for_mod, family = binomial())
ef <- exposure_model(exposure_fit, numerator_model = NULL)



ppm(pp ~ treatment1 + covarfun)
ppm(pp ~ treatment + s(x,y, k = 40, bs = 'gp'), use.gam = TRUE)
glmdata <- getglmdata(ppm(pp ~ alltessf + treatment))

ff <- glm(.mpl.Y ~ treatment + alltessf, weights = .mpl.W, data = glmdata, family = poisson())

treatment



## simulation 2
## pointwise binary exposure, including covariates, but not confounding
## pointwise on the fine 128 by 128 grid used for computation
## technically this is still blockwise, but for small blocks

treatment2 <- rbinom(128 * 128, size = 1, prob = 0.5)
treatment2 <- matrix(treatment2, 128, 128)
treatment2 <- as.im(treatment2, window)

intens <- exp(intercept + treatment2 + covar)
pp <- rpoispp(intens)
ppm(pp ~ treatment2 + covar)



## simulation3
## binary random field exposure
## generated from matern random field and then clipped on the basis of being greater or lower than zero
gridx <- seq(from = 0, to = 1, length.out = 128)
treatment3 <- as.im(as.matrix(RFsimulate(RMmatern(var = 2, scale = 0.3, nu = 1.3), gridx, gridx)) > 0, window)


intens <- exp(intercept + treatment3 + 0.1 * covar)

plot(intens)
pp <- rpoispp(intens)
fit <- ppm(pp ~ treatment3 + covar, use.gam = TRUE)
fit <- ppm(pp ~ treatment3 + covar + s(x,y, bs = 'gp', k = 40),  use.gam = TRUE)

## simulation 4, binary exposure with lower degree of spatial autocorrelation
treatment4 <- as.im(as.matrix(RFsimulate(RMmatern(var = 2, scale = 0.05, nu = 1.3), gridx, gridx)) > 0, window)


intens <- exp(intercept + treatment4 + 0.1 * covar)
pp <- rpoispp(intens)
fit <- ppm(pp ~ treatment4 + covar, use.gam = TRUE)
fit <- ppm(pp ~ treatment4 + covar + s(x,y, bs = 'gp', k = 40),  use.gam = TRUE)
dd <- getglmdata(fit)
fit <- gam(.mpl.Y ~ treatment4 + covar + te(x,y, bs = 'gp', k = 10), weights = .mpl.W, data = dd, family = poisson())


treatment5 <- as.im(as.matrix(RFsimulate(RMmatern(var = 2, scale = 0.001, nu = 1.3), gridx, gridx)) > 0, window)

intens <- exp(intercept + treatment5 + 0.1 * covar)
pp <- rpoispp(intens)
fit <- ppm(pp ~ treatment5 + covar, use.gam = TRUE)
fit <- ppm(pp ~ treatment5 + covar + s(x,y, bs = 'gp', k = 40),  use.gam = TRUE)
dd <- getglmdata(fit)
fit <- gam(.mpl.Y ~ treatment5 + covar + te(x,y, bs = 'gp', k = 10), weights = .mpl.W, data = dd, family = poisson())

## simulation structure



## next round of simulations
## consider scenarios where covariates are themselves spatially autocorrelated, so there is residual spatially determined confounding.
spatial_covar <- as.im(as.matrix(RFsimulate(RMmatern(var = 2, scale = 0.1, nu = 1.3), gridx, gridx)), window)
intens <- exp(intercept + treatment3 + spatial_covar)

pp <- rpoispp(intens)
fit <- ppm(pp ~ treatment3 + spatial_covar, use.gam = TRUE)
fit <- ppm(pp ~ treatment3)
fit <- ppm(pp ~ treatment3 +  s(x,y, bs = 'gp', k = 40),  use.gam = TRUE)


## explicitly confounded treatment on the basis of spatially autocorrelated covariate
treatment6 <- as.im(matrix(rbinom(128^2, 1, plogis(c(as.matrix(spatial_covar)) * 0.3 - 1)), 128, 128), window)

intens <- exp(intercept + treatment6 + spatial_covar)
pp <- rpoispp(intens)

fit <- ppm(pp ~ treatment6,  use.gam = TRUE)
fit <- ppm(pp ~ treatment6 +  spatial_covar,  use.gam = TRUE)
fit <- ppm(pp ~ treatment6 +  s(x,y, bs = 'gp', k = 40),  use.gam = TRUE)



## treatment 7 is a spatially autocorrelated treatment resulting counfounded by a spatially autocorrelated covariate,
## but that is realized as a direct sum of the covariate and a secondary random process.
## as a result, this treatment has additional spatial autocorrelation and is 'clumpy'
residual_process <- RFsimulate(RMmatern(1.2, var = 0.8, scale = 0.7), gridx, gridx)
treatment7t <- as.im(as.matrix(residual_process), window) + spatial_covar
treatment7t <- as.function(treatment7t)
treatment7 <- function(x,y) treatment7t(x,y) >0
treatment7 <- as.im(treatment7, window)


intens <- exp(intercept + treatment7 + spatial_covar)
pp <- rpoispp(intens)

fit <- ppm(pp ~ treatment7,  use.gam = TRUE)
fit <- ppm(pp ~ treatment7 +  spatial_covar,  use.gam = TRUE)
fit <- ppm(pp ~ treatment7 +  s(x,y, bs = 'gp', k = 40),  use.gam = TRUE)



treatment8
## generate treatment as set of circles where circle centers are generated from a poisson process on the basis of the confounding covariate
set.seed(0)
pp_base <- spatial_covar - min(spatial_covar) + 1

## generate approximately 10 exposure centers
pp_base <- 10 * pp_base / integral(pp_base)

centers <- rpoispp(pp_base)

treat_surf <- distfun(centers)

treatment8 <- as.im(as.function(function(x,y) treat_surf(x,y) < 5000), window)

intens <- exp(intercept + treatment8 + spatial_covar)
pp <- rpoispp(intens)

fit <- ppm(pp ~ treatment8,  use.gam = TRUE)
fit <- ppm(pp ~ treatment8 +  spatial_covar,  use.gam = TRUE)
fit <- ppm(pp ~ treatment8 +  s(x,y, bs = 'gp', k = 40),  use.gam = TRUE)

###########################
## CONTINUOUS TREATMENTS ##
###########################


## treatment is spatial covar plus a residual surface
treatment9 <- as.im(treatment7t, window)

intens <- exp(intercept + treatment9 + spatial_covar)
pp <- rpoispp(intens)

fit <- ppm(pp ~ treatment9,  use.gam = TRUE)
fit <- ppm(pp ~ treatment9 +  spatial_covar,  use.gam = TRUE)
fit <- ppm(pp ~ treatment9 +  s(x,y, bs = 'gp', k = 40),  use.gam = TRUE)


## realistic and potentially advantageous scenario
## we have a continuous treatment generated as a function of space, but spatial confounders have discontinuities due to geographical features
sim_surface <- function(){
    as.im(as.matrix(RFsimulate(RMmatern(nu = 1, var = 2, scale = 1), gridx, gridx)), window)
}
s <- map(1:4, ~sim_surface())
## determine data generation mechanism
## pick optimal model from candidate model set according to specified criteria
## calculate truth value of specified estimand under dgp
## generate data
## calculate estimator values on data

## aggregate estimator values
## calculate summary statistics on estimators
## combine summary statistics across simulations for comparisons

broken_covar <- function(x,y){s[[allegheny$region[as.numeric(alltessf(x,y))]]](x,y)}
covar <- as.im(broken_covar, window)    
