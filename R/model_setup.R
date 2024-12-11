
convolve_basis <- function(fft_raster, basis_fun, window) {
    ## zero_pad basis_fun
    resolution <- dim(fft_raster)
    n_rows <- resolution[[1]]
    n_cols <- resolution[[2]]
    
    buf <- matrix(0, n_rows, n_cols)
    buf[1:(n_rows %/% 2), 1:(n_cols %/% 2)] <- basis_fun
    
    conv_result <- Re(fft(fft(buf) * fft_raster, inverse = TRUE) / prod(resolution))
    outx <- n_rows %/% 4 + 1:(n_rows %/% 2)
    outy <- n_cols %/% 4 + 1:(n_cols %/% 2)
    conv_out <- conv_result[outx, outy]
        
    outrast <- spatstat.geom::as.im(conv_out, W = window)
    outrast
}

rasterize <- function(ppp, resolution) {
    as.im(ppp, dimyx = resolution)
}

#' @export
convolutional_basis_setup <- function (covariate,
                                       quadrature_points, 
                                       n_basis, knots, 
                                       resolution,
                                       order,
                                       max_distance, distance_scale,
                                       keep_interpolation = FALSE,
                                       force_window = NULL) 
{
  if (missing(max_distance)) 
  {
    stopifnot(missing(distance_scale))
    distance_prop <- 0.50
    message("Using default distance proportion of 0.50\n")
  }
  else {
    distance_prop <- max_distance / distance_scale
  }
  
  quadrature_coords <- as_matrix(quadrature_points)

    if(!missing(knots)) stop("Only default knots currently supported")
    if(!missing(order)) stop("Currently only implemented for order = 4")
    if (missing(resolution))
    {
      message ("Using default 128 x 128 resolution.\n")
      resolution <- c(128, 128)
    }
  

    n_rows <- resolution[1]
    n_cols <- resolution[2]


    ## setup point process covariate
    if (inherits(covariate, "ppp")) {
      win <- covariate$window
      covariate <- rasterize(covariate, resolution = resolution)
      covariate <- list(covariate)
      window_list <- list(win)
      ncovar <- 1L
    }
    else if (inherits(covariate, "im")) {
      if(!is.null(covariate$window)){
          win <- covariate$window
        }
        else {
          win <- owin(covariate$xrange, covariate$yrange)
        }

      covariate <- covariate
      covariate <- list(covariate)
      window_list <- list(win)
      ncovar <- 1L
    }
    else if (is.list(covariate)){
      stopifnot(all(sapply(covariate, \(ob) inherits(ob, "ppp"))))
      
      ncovar <- length(covariate)
      window_list <- vector("list", ncovar)
      
      for (i in 1:ncovar){
        if(!is.null(covariate[[i]]$window)){
          win <- covariate[[i]]$window
        }
        else {
          win <- owin(covariate[[i]]$xrange, covariate[[i]]$yrange)
        }
        window_list[[i]] <- win
        covariate[[i]] <- rasterize(covariate[[i]], resolution = resolution)
      }
    }
    
    if (!is.null(force_window)) {
      for (i in 1:length(window_list)){
        window_list[[i]] <- union.owin(window_list[[i]], force_window)
      }
    }
    
    else {
      stop("Unknown covariate type")
    }
    ## setup basis expansion
    xseq <- seq(from = -1L, to = 1L, length.out = n_rows)
    yseq <- seq(from = -1L, to = 1L, length.out = n_cols)

    distrast <- outer(xseq, yseq, function(x,y) sqrt(x^2 + y^2))
    
    basis <- bspline_basis(c(distrast), n_basis, order = 4, max_knot = distance_prop)
    penalty <- wiggliness_penalty(basis)

    interpolation_basis_single <- vector("list", ncol(basis))
    interpolation_basis <- vector("list", ncovar)
    convolved_basis_list <- vector("list", ncovar)

    
    for (covari in 1:ncovar) {
      covariate_buffer <- matrix(0, n_rows * 2, n_cols * 2)
      vtmp <- covariate[[covari]]$v
      vtmp <- ifelse(is.na(vtmp), 0, vtmp)
      covariate_buffer[1:n_rows, 1:n_cols] <- vtmp
      covariate_buffer <- fft(covariate_buffer)
      
      convolved_basis <- matrix(0, nrow = nrow(quadrature_coords), ncol = ncol(basis))
      win <- window_list[[covari]]
      
      for (basis_column in 1:ncol(basis)) {
          local_conv <- convolve_basis(covariate_buffer, basis[,basis_column], win)
          interpolation_basis_single[[basis_column]] <- local_conv
          convolved_basis[ ,basis_column] <- interp.im(local_conv, quadrature_coords)
      }
      if(keep_interpolation) {
        interpolation_basis[[covari]] <- interpolation_basis_single
      }
      convolved_basis_list[[covari]] <- convolved_basis
    }
    
    convolved_basis <- do.call(rbind, convolved_basis_list)

    if(!keep_interpolation){
      interpolation_basis = NULL
    }
    out <- list(
        X = convolved_basis,
        interpolation_basis = interpolation_basis,
        distance_basis = basis,
        penalty = penalty
    )
    class(out) <- "convolutional_spline_basis"
    out
}

#' @exportS3Method
print.convolutional_spline_basis <- function ( object, ...){
    cat("Convolutional spline basis \n")
}
    




#' @export
Predict.matrix.Convspline.smooth <- function(object, data) {
    
    coords <- extract_coords(data[[object$term]])
    
    interp_basis <- object$interpolation_basis

    ncoord <- length(coords)
    
    nr <- ncoord * length(interp_basis)
    ## assumes all bases have same dimension which they certainly should
    nc <- length(interp_basis[[1]])
    Xmat <- matrix(0, nrow = nr, ncol = nc)

    for (i in 1:length(interp_basis)){

        locinterp <- interp_basis[[i]]

        for (j in 1:nc){
            ## e <<- environment()
            row_range <- (((i - 1) * ncoord) + 1) : (((i) * ncoord) )
            interp <- evaluate(locinterp[[j]], coords)
            Xmat[row_range, j] <- interp
        }
    }

    Xmat
}



## TODO
#' @export
update_exposure <- function(model, new_exposure) {

    newdata <- model$gam_data
    ## data checks
    stopifnot(is.list(new_exposure))
    if(any(is.null(names(new_exposure)))) stop('New exposures must have names')
    stopifnot(all(names(new_exposure) %in% model$ppcovs))


    which_covs <- model$ppcovs %in% names(new_exposure)
    covs <- model$ppcovs[which_covs]


    for (cov in covs) {
        newdata[[cov]] <- covariate_placeholder(
            Pcov(new_exposure[[cov]], model$Y, model$dim),
            extract_coords(model$gam_data[[cov]])
        )
    }
    smooths <- fit$smooth

    for(i in seq_along(smooths)){
        if(is(smooths[[i]], 'Convspline.smooth')){
            object <- smooths[[i]]
            cov <- object$term
            conv_data <- extract_data(newdata[[cov]])

            covlen <- length(conv_data)

            

            basis <- list( c(field(conv_data, 'pcov')[[1]]$distances))
            names(basis) <- object$term
            
            new_basis <- Predict.matrix(object$internal_basis, basis)

            
            object$interpolation_basis <- apply(
                new_basis, 2,
                function(basis) convolve_basis(basis,
                                               field(conv_data, 'pcov')[[1]]$covariate,
                                               field(conv_data, 'pcov')[[1]]$dims,
                                               field(conv_data, 'pcov')[[1]]$window,
                                               coords)
            )

            smooths[[i]] <- object
        }
    }
    model$smooth <- smooths
    model$counterfactual_covs <- covs

    if(!is(model, 'counterfactual')) class(model) <- c('counterfactual', class(model))

    model
}

