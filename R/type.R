## ... possible configuration options
Define <- function(sym, ...){
    structure(
        list(sym = substitute(sym)),
        class = "Definition_declaration"
    )
}

`[.Definition_declaration` <- function(fdec, ...){

    structure(
        list = (fdec = fdec),
        class = "Type_declaration")

    function(body){}
}


f <- function(body) {
    substitute(body)
}

function(){}  

Type <- function (typename){
    structure(list(typename = substitute(typename)), class = "Type")
}


`*.Type` <- function (x, y )
{
    structure(
        list(
            typename = paste0(typename(x), ", ", typename(y)),
            t0 = x,
            t1 = y
        ),
        class = c("ProductType", "Type")
    )
}

print.ProductType <- function(ob, ...) {
    cat(paste0("( ", typename(ob$t0), ", ", typename(ob$t1), " )"))
    cat("\n")
}


FunctionType <- function(Dom, Im) {
    structure(
        typename = paste0(typename(Dom), " -> ", typename(Im)),
        list(
            Dom = Dom,
            Im = Im
        ),
        class = c("FunctionType", "Type")
    )
}

print.FunctionType <- function(ob, ...) {
    cat(paste0(typename(ob$Dom), " -> ",  typename(ob$Im), "\n"))
}



typename <- function(type) {
    type$typename
}


print.Type <- function (ob, ...) {
    cat(ob$typename)
    cat("\n")
}

Real <- Type(Real)








