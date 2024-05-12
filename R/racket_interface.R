out <- NULL
set_out <- function(str) out <<- sys::as_text(str)


f <- sys::exec_wait('ls', std_out = set_out)
