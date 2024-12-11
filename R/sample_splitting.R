#' Generate Sample Splitting Folds
#'
#' Sample splitting folds ensure independence  for various inferential purposes.
#' In common use sample splitting and cross-validation are two names for the same general idea.
#' However, it the author's view of the literature as of the current writing (September 2024), cross validation is more commonly used
#' in a general predictive context where its primary purpose is for assesment of model fit.
#' We use the term sample splitting to refer to a general inferential purpose instead.
#' In particular sample splitting type estimators use sample splitting strategies to ensure independence of estimates of influence functions from the estimated nuisance parameters.
#'
#' While we emphasize the two ideas are essentially the same, we use this distinction to convey intent of the design of the library.
#'
#' Consequently for many functions there will be a "ss" class which is sample splitting and a "cv" class for cross validiation.
#' Depending on class default behavior of methods will generally differ. 
random_sample_split <- function (df, n_folds, proportion_exclude = 0.20, seed = NULL) {
    folds <- random_folds(df, n_folds = n_folds, proportion_exclude = proportion_exclued, ...)

        structure(df,
              cv_folds = folds,
              class = c("random_ss", "ss", class(df)))
}

#' Spatial Sample Splitting
#'
#' @export
spatial_sample_split <- function(df,
                                 geo_regions,
                                 sampling_contexts = geo_regions,
                                 buffer = NULL,
                                 markov = !is.null(buffer),
                                 independent = is.null(buffer) & !markov)
{
    
    
}




#' Random Folds
#'
#' Generates random folds for use in sample splitting or cross validation
random_folds <- function(df, n_folds, proportion_exclude = 0.20, seed = NULL) {
    if (!is.null(seed)) stop("random_cv currently not implemented to deal with the seed argument")
    ## TODO check if df is grouped
    n_indices <- nrow(df)
    indices <- 1:n_indices

    n_out <- ceiling(proportion_exclude * n_indices)
    folds <- array(1L, dim = c(n_indices, n_folds))
    
    for (fold in 1:n_folds) {
        exclude <- sample(indices,  n_out, replace = FALSE)
        folds[exclude, fold] <- 0L
    }
    folds
}


## spatio-temporal cross-validation framework
random_cv <- function (df, n_folds, proportion_exclude, ...) { 
    folds <- random_folds(df = df, n_folds = n_folds, proportion_exclude = proportion_exclude, ...)
    
    folds
}

print.cv <- function(object) {
    NextMethod()
    cat("With CV folds\n")
}

print.ss <- function(object) {
    cat("Sample Splitting Environment\n")
    NextMethod()
}


fit_with_ss <- function ( ) {
    
}


fit_with_cv <- function ( ) {

}

