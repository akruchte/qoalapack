#' Parameter Objects
#' @name param
#' @param sym A symbolic (unquoted) R expression denoting the name of a parameter
#' @param par A parameter object created with `param`
#' @export
param <- function(sym){
    if (missing(sym)) {
        if(!exists("...GlobalParCounter...")) {
            ...GlobalParCounter... <<- 0L
        }
        ...GlobalParCounter... <<- ...GlobalParCounter... + 1L
        sym <- str2lang(paste0("...anon_par", ...GlobalParCounter...))
    }
    structure(list(substitute(sym)),
              shape = NULL,
              value = NULL,
              fixed = FALSE,
              class = "param")
}

#' @rdname param
#' @exportS3Method
print.param <- function (object, ... ) {
    cat(paste0("param: ", object, "\n"))
}

#' @export
`bind<-` <- function(par, value ) {
    if (attr(par, "fixed")) stop("Attempted to re-assign fixed param")
    attr(par, "value") <- value
    par
}

#' @rdname param
#' @export
`fix` <- function(par, value) {
    attr(par, "fixed") <- TRUE
    par
}

#' @export
`unfix` <- function (par, value) {
    attr(par, "fixed") <- FALSE
    par
}


#' @rdname param
#' @export
is_fixed <- function(par) {
    attr(par, "fixed") 
}


#' @rdname param
#' @export
value_of <- function(par, value) {
    attr(par, "value")
}

#' @rdname param
#' @export
name_of <- function(par) {
    as.character(par)
}
