## If we want to make an object revisable/loggable
## we can turn that function into an active binding that returns the original function, captures the results of that call,
## logs those results,and does prior initialization as required.

## MAKE_REVISABLE(gam) should add to the global namespace a new variable gam._1
## or alternatively given by a specified name as well as an associated logger
## .log._1 or .log.name

MAKE_REVISABLE <- function(name, logger){

    env <- rlang::quo_get_env(fun)
    sym <- rlang::quo_get_expr(fun)

    if (logger != NULL) {
        logger <- default_logger_initialize(fun)
    }
    
    act_binding <- function ( ) {
        
    }

    makeActiveBinding(sym, act_binding ,env)
}








