## Ecov is for representing exposure based covariates.
## These are often used in general in ways that are quite similar to Fcovs, but may be more appropriate for modeling as the outcome


#' @export
Ecov <- function(..., outcome) {
    object <- rlang::list2(...)
    fitted <- lapply(object, function(ob) Ecov_prepare(ob, ...))

    vctrs::new_rcrd(list(outcome), class = c('Ecov', 'spatial_covariate'))
}

#' @export
format.Ecov <- function(ob, ...){
    rep('Ecov', vctrs::vec_size(ob))
}


