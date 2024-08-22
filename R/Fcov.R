
##Fcov
## convert to Fcov
#' @export
as_Fcov <- function(object, ...){
  UseMethod('as_Fcov')
}

#' @export
as_Fcov.Fcov <- function(object, ...){
  object
}

#' @export
as.im.Fcov <- function(object, W) {
  as.im(function(x,y) evaluate(object, coord(x,y)), W = W)
}

#' @export
as_Fcov.Pcov <- function(object, kernel, ...){
  force(kernel)
  f <- function(x,y) evaluate(object, locations = coord(x,y), kernel = kernel)
  Fcov(f)
}

#' @export
as_Fcov.Lcov <- function(object, kernel, ...){
  force(kernel)
  f <- function(x,y) evaluate(object, locations = coord(x,y), kernel = kernel)
  Fcov(f)
}




## Fcov is similar in spirit to Rcov, but handles data via a functional representation
## this may be under some circumstances more amenable to handling uncertainty due to measurement error
## Fcov requires a basis or model for the representation of the data

#' @export
Fcov <- function(..., formula, method, family){
  object <- list2(...)
  fitted <- lapply(object, function(ob) Fcov_prepare(ob, ...))

  vctrs::new_rcrd(list(fitted = fitted), class = c('Fcov', 'spatial_covariate'))
}
#' @export
format.Fcov <- function(ob, ...){
  rep('Fcov', vec_size(ob))
}

#' @export
`+.Fcov` <- function(ob, ob2) {
  if (inherits(ob2, 'numeric')) {
    result <- map2(field(ob, 'fitted'),
                   ob2,
                   function(f1, f2) function(x,y) f1(x,y) + f2)
    return(Fcov(!!!result))
  }
  if(inherits(ob2, 'Fcov')){

    result <- map2(field(ob, 'fitted'),
                   field(ob2, 'fitted'),
                   function(f1, f2) function(x,y) f1(x,y) + f2(x,y))
    return(Fcov(!!!result))
  }
  stop('Not implemented')
}

#' @export
`*.Fcov` <- function(ob, ob2) {
  if (inherits(ob2, 'numeric')) {
    result <- map2(field(ob, 'fitted'),
                   ob2,
                   function(f1, f2) function(x,y) f1(x,y) * f2)
    return(Fcov(!!!result))
  }
  if(inherits(ob2, 'Fcov')){

    result <- map2(field(ob, 'fitted'),
                   field(ob2, 'fitted'),
                   function(f1, f2) function(x,y) f1(x,y) * f2(x,y))
    return(Fcov(!!!result))
  }
  stop('Not implemented')
}
#' @export
`-.Fcov` <- function(ob, ob2) {
  if (inherits(ob2, 'numeric')) {
    result <- map2(field(ob, 'fitted'),
                   ob2,
                   function(f1, f2) function(x,y) f1(x,y) - f2)
    return(Fcov(!!!result))
  }
  if(inherits(ob2, 'Fcov')){

    result <- map2(field(ob, 'fitted'),
                   field(ob2, 'fitted'),
                   function(f1, f2) function(x,y) f1(x,y) - f2(x,y))
    return(Fcov(!!!result))
  }
  stop('Not implemented')
}
#' @export
`/.Fcov` <- function(ob, ob2) {
  if (inherits(ob2, 'numeric')) {
    result <- map2(field(ob, 'fitted'),
                   ob2,
                   function(f1, f2) function(x,y) f1(x,y) / f2)
    return(Fcov(!!!result))
  }
  if(inherits(ob2, 'Fcov')){

    result <- map2(field(ob, 'fitted'),
                   field(ob2, 'fitted'),
                   function(f1, f2) function(x,y) f1(x,y) / f2(x,y))
    return(Fcov(!!!result))
  }
  stop('Not implemented')
}



#' @export
evaluate.Fcov <- function(ob, locations){
  x <- coordx(locations)
  y <- coordy(locations)
  preds <- map(field(ob, 'fitted'), function(f) f(x,y))
  reduce(preds, c)
}

setup.Fcov <- evaluate.Fcov

#' @export
Fcov_prepare <- function(object, ...){
  UseMethod('Fcov_prepare')
}
#' @export
Fcov_prepare.function <- function(object, ...) {
  stopifnot(names(formals(object)) == c('x', 'y'))
  object
}
## basis_select should, given a particular basis type, be able to pick the relevant basis elements and perform any sort of low-rank computations required
## basis select needs to be able to consider selecting between various types of basis, (cubic spline, thin plate, wavelet, fpca, gaussian process)
## it also needs to be able to make basis constructions either jointly across all covariate values at the same time (pooling the data)
## it should also be able to make pooling decisions using random subsets of the pooled data (for representativeness while also reducing computation costs).
## Alternatively it should be able to select bases for each model individually
## A final alternative is to use a specified portion of the data to make a single global basis selection
#' @export
Fcov_prepare.data.frame <- function(data, formula, method, family){
  method(formula, family = family, data = data)
}
#' @export
basis_select <- function(...){

}




## basis_project should perform model fitting to get basis representations of each of the functional covariates
## It should be able to fit each portion entirely seperately with its own smoothing parameter selection
## Alternatively, it should be additionally possible to pool smoothing parameter selection
## Finally, it should also be possible to fit models with a single totally pooled model across all covariate data (such as a large spatiotemporal regression)



## evaluate should be able to evaluate the fitted functional components at a set of specified coordinates

## Arguments to Fcov could be reasonably provide as either a list of data frames corresponding to the component regressions,
## or a single data frame corresponding to the total regression

## there could potentially be different basis selections for different components, e.g. if data is highly non-stationary over time

## new_rcrd(fields = list(coefs = coefs,
##                        ## the actual basis representation, coupled with the ability to interpret the basis at new points is critical
##                        which_basis,
##                        ## could possibly need influence function for down the line work
##                        ## this may be better computed via autodiff or explicit derivatives
##                        influence = influence,
##                        ## keeping the data as is would be way too memory intensive
##                        ## may need to refer to data in order to be able to calculate confidence intervals down the line
##                        data = data))
## ## the attributes of the record should contain the basis, or bases if different bases are used at varying times

