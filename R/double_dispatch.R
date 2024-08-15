generic <- function(x, y, ...){
    UseMethod('generic')
}

generic.default <- function(x,y, ...) {
    cat('Did not make it')
}

generic.test <- function(x, y, ...) {
    UseMethod('generic.test', y)
}

generic.test.test <- function(y, x, ...){
    cat('Double dispatched on clone')
}

generic.test.test2 <- function(y, x, ...) {
    cat('Double dispatched')
}

x <- structure(1, class = 'test')
y <- structure(2, class = 'test2')
