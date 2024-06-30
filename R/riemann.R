#' @export
base_bump <- function(
                 a=1, b=1, c=1,
                 d = 1, p = 2,
                 trans = identity,
                 tropical_op = min,
                 exponential_map = exp,
                 bump_parity = -1
                 )
{
     exponential_map(bump_parity * (a/(b - c * tropical_op(d, trans(x)^2))))
}

bump_location_scale <-function(base_bump, location, scale) {
    function (x){base_bump(x - location) *scale}
}


## take a boundary region and deform to a circular space
## find be easiest to do this via simplex structure, open simplex structure to euclidean two dimensional space and then encode as a complex value
#' @export
boundary_to_circle <- function(){}

#' @export
circle_to_segment <- function(){}

#' @export
boundary <- function(sf_rectifiable_or_spatstat_window, win = sf_retifiable_or_spatstat_window)
{
    opening <- window_internal
}

## zeta prime is the computable zeros equivalent of the zeta function
## General idea is we can use zeta "renormalization" via compact bumps truncations to zero  and actual zeta function applications.
zetap <- function(x){}


## point process and associated persistent homology with softcore constraints?



