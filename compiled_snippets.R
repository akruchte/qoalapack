## basic AIPW 

library(spatstat)

win <- owin(c(-1, 1), c(-1, 1))

beta <-t( t(rnorm(6)))
treatment_link <- function(x,y) cbind(1, x, y, x^2, y^2, x * y) %*% beta

treatment_mu <- function(x,y) plogis(treatment_link(x,y))

treatment <- function(x,y) rbinom(length(x), 1, prob = treatment_mu(x,y))

beta2 <- t(t(rnorm(6)))
confound <- function(x,y) cbind(1, x, y, x^2, y^2, x * y) %*% beta2

treatment <- as.function(as.im(treatment, win))
confound <- as.function(as.im(confound, win))

p <- process <- rpoispp(function(x,y) exp(2 + 2 * treatment(x,y) + confound(x,y)), win = win)


tt <- as.im(treatment)

for_reg <- mutate(expand_grid(x = tt$xcol, y = tt$yrow), v = c(tt$v))

prop_fit <- glm(v ~ poly(x,y, degree = 2), data = for_reg, family = binomial())

prop <- predict(prop_fit, type = 'response')
vv <- tt$v

vv[] <- prop
propensity_score <- tt
propensity_score$v <- vv
outcome <- ppm(p ~ treatment +poly(x,y, degree = 2))


tt_new = as.im(1, win)
new <- list(treatment = tt_new)
mu <- predict(outcome, type = 'intensity', covariates = new)
new <- list(treatment = tt_new - 1)
mu2 <- predict(outcome, type = 'intensity', covariates = new)

Y <- as.im(p)
a <- c(tt$v)

## looks correct
## first term is integration with respect to intensity measure, i.e. the expectaiton.
second term is integration with respect to residual measure, i.e. 
y <- integral(mu) + sum(tt * Y / propensity_score) - integral(tt * mu / propensity_score)
y2 <- integral(mu2) + sum((1-tt) * Y / (1-propensity_score)) - integral((1-tt) * mu2 / (1-propensity_score))





## basic composition

library(tidyverse)
library(tidycensus)
library(sf)

d <- read_csv('~/Desktop/decennial_race_data/DECENNIALPL2020.P1_data_with_overlays_2021-12-20T114314.csv', skip = 1)
dr <- d %>% select(1:2, contains('Population of one race')) 

dr <- dr %>% select(!3)


library(compositions)

dr %>% pack(comp = c(3:8)) %>%
    mutate(comp = ilr(comp))
bg <- tigris::block_groups(state = 'PA')

bg <- bg %>% transmute(id = GEOID)
dr <- dr %>% mutate(id = str_extract(id, '(?<=US).+'))

dr <- dr %>% select(1, 2, 3, 4, 6, 8) %>% mutate(across(!c(1:2), ~ . + ))

dr <- dr %>% pack(comp = (3:6)) %>%
    mutate(comp = ilr(comp))
dr <- dr %>% mutate(comp = as.data.frame(comp)) %>%
    unpack(cols = comp)
dr <- dr %>% select(-2)

data <- inner_join(dr, bg)
data <- st_as_sf(data)

dd <- dist(select(as_tibble(data), V1, V2, V3))

dd <- as.matrix(dd)
example <- 1

datat <- data %>% mutate(dist_from = c(dd[example,]))

point <- st_centroid(data[example,])

ggplot() + scale_fill_viridis_c() +
    geom_sf(aes(fill = dist_from), data = datat) +
    geom_sf(data = point, size = 2)



## basic raster fits

library(tidyverse)
library(stars)
library(partykit)
d <- read_ncdf('~/Downloads/V5GL02.HybridPM25.NorthAmerica.201907-201907.nc')

d <- filter(d, lon > -91, lon < -90, lat < 41, lat > 40)

r <- d$GWRPM25
y <- c(r)
co <- st_coordinates(d)
dd <- bind_cols(co, y = y)
library(mgcv)
fit <- gam(y ~ s(lon, lat, bs = 'gp', m = 5), data = dd)
r2y <- predict(fit)
r2 <- r
r2[!is.nan(r2)] <- r2y

par(mfrow = c(3,1))
image(r)
image(r2)
image(r - r2)

fit <- lmtree(y ~ lon + lat | lon + lat, data = dd)
p <- predict(fit, newdata = dd)

r3 <- r2
r3[] <- p
image(r3)

resid <- r - r3

fit <- gam(y ~ s(lon, lat, bs = 'gp', m = 5), data = mutate(dd, y = c(resid)))

p2 <- predict(fit, newdata = dd)


library(rSPDE)
library(INLA)

## loc_2d_mesh = co
## mesh_2d = inla.mesh.2d(
##   loc=loc_2d_mesh,
##   cutoff=0.05,
##   max.edge=c(0.1,0.5) )


## Abar <- rspde.make.A(mesh = mesh_2d, loc = as.matrix(loc_2d_mesh))
## mesh.index <- rspde.make.index(name = 'field', mesh = mesh_2d)

## st.dat <- inla.stack(
##   data = list(y = as.vector(r)),
##   A = Abar,
##   effects = mesh.index)


## rspde_model <- rspde.matern(mesh = mesh_2d, nu_upper_bound = 1)


## f = y ~ -1 + f(field, model = rspde_model)
## ## inla currently not working
## ## rspde_fit <- inla(f, data = inla.stack.data(st.dat),
##                   ## family = 'gaussian',
##                   ## control.predictor = list(A = inla.stack.A(st.dat)),
## ## inla.mode = 'experimental', verbose = TRUE)

mlik <- function(theta, Y, A, op) {
  sigma = exp(theta[1])
  kappa = exp(theta[2])
  nu = exp(theta[3])
  return(-rSPDE::rSPDE.matern.loglike(object = op,
                                      Y = Y,
                                      A = A,
                                      user_kappa = kappa,
                                      user_sigma = sigma,
                                      user_nu = nu, sigma. = exp(theta[4])))
}

loc_2d_mesh = co[!is.nan(y), ]
y <- y[!is.nan(y)]
mesh_2d = inla.mesh.2d(
  loc=loc_2d_mesh,
  cutoff=0.05,
  max.edge=c(0.1,0.5) )


library(optimParallel)
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)

setDefaultCluster(cl = cl)
sigma <- 1
range <- 0.2
nu <- 0.5
kappa <- sqrt(8*nu)/range
op <- matern.operators(mesh=mesh_2d,nu=nu,
                                   kappa=kappa,sigma=sigma,m=2)

theta0 = c(get.inital.values.rSPDE(mesh = mesh_2d),
           log(0.1 * sqrt(var(as.vector(y), na.rm = TRUE))))

A <- inla.spde.make.A(
  mesh=mesh_2d,
  loc=as.matrix(loc_2d_mesh))

sigma.e <- 0.1
pars <- optimParallel(theta0, mlik, Y = y, A = A, op = op)
results <- data.frame(sigma = exp(pars$par[1]), 
                      kappa =  exp(pars$par[2]),
                      nu = exp(pars$par[3]),
                      sigma.e = exp(pars$par[4]),
                      row.names = "Estimate")


predgrid <- inla.mesh.projector(mesh_2d, xlim = c(0,1), 
ylim = c(0,1))
A.prd2 <- predgrid$proj$A


sigma.e.est <- exp(pars$par[4])
op.prd <- update(op, user_sigma = exp(pars$par[1]),
                 user_kappa = exp(pars$par[2]),
                 user_nu = exp(pars$par[3]))

pred.rspde <- predict(op.prd, A = A, Aprd = A, Y = y, 
                      sigma.e = sigma.e.est,
                      compute.variances = TRUE)
rp <- r
rp[!is.nan(r)] <- pred.rspde$mean
