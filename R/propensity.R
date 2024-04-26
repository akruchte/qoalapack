library(spatstat)
fit <- ppm(swedishpines ~ x + y, interaction = Strauss(5))

ss <- simulate(fit)

sample <- ss[[1]]
model <- fit


log_density <- function(sample, model)
{
    intens <- predict(model, locations = sample, type = 'trend')

    intensity_potential <- sum(log(intens))

    interaction_potential <- 0

    if (! is.null(model$interaction)) {
        d <- pairdist(sample)
        d <- d[upper.tri(d)]
        pot <- 2 * model$interaction$pot(d, model$interaction$par)
        interaction_potential <- sum(coef(model)[['Interaction']] * pot)
    }


    potential <- interaction_potential + intensity_potential
    potential
}

bridge <- function(n1, n2, sample, model) {
    s1 <- n1 / (n1 + n2)
    s2 <- n2 / (n1 + n2)

    density(sample, model) * s1 * density(sample,
}

bridge_sample <- function(proposal, sample, density, model) {
    numerator <- mean(bridge(proposal) * normalization_constant(proposal) * density(proposal, model)
    denominator <- mean(bridge(sample) * density(sample, model))
}
