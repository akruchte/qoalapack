



aipw_term <- mu + (1/ pi) * (y - mu)
logaipw_term <- log(exp(log(mu)) + exp(log(y-mu) - log(pi)))



y <- 0
pi <- 0.05
mu <- 1



terms <- c(log(mu), log(y-mu + 0i) - log(pi))

exp(terms[2] + log(sum(exp(terms - terms[2]))))


