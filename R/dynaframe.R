dynadouble  <- function(init_values) {
    if (!is.double(init_values)) abort("init_values must be of type (double)")

    new_vctr(init_values, class = 'dynadouble')
}





`+.fun`<- function(x,y){function(z) x(z) + y(z)}


