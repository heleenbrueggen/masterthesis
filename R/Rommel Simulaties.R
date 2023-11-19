##########
# Rommel #
##########
#####
# Simulation data with only one x and z variable
#####
simdata <- lapply(1:1000, function(i) {
  data <- tibble(
    id = 1:(ngroup * groupsize),
    group = rep(1:ngroup, each = groupsize),
    eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
    x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
    z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u0 = rep(rnorm(n = ngroup,
                   mean = 0,
                   sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, eij, icc) {

                     daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var(eij))))

                     return(daticc)
                   }, interval = c(0, 100),
                   tol = .0001,
                   extendInt = 'yes',
                   maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
    beta0j = g00 + g01 * z + u0,
    beta1j = g10 + g11 * z + u1,
    y = beta0j + beta1j * x1 + eij)
}
)
#####
# Simulation like in thesis proposal
#####
simdata <- lapply(1:1000, function(i) {
  data <- tibble(
    id = 1:(ngroup * groupsize),
    group = rep(1:ngroup, each = groupsize),
    eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
    x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
    x2 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3),
    z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u2 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u0 = rep(rnorm(n = ngroup,
                   mean = 0,
                   sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, g20, u2, eij, icc) {

                     daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var((g20 + u2)*x2) + var(eij))))

                     return(daticc)
                   }, interval = c(0, 100),
                   tol = .0001,
                   extendInt = 'yes',
                   maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, g20 = g20, u2 = u2, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
    beta0j = g00 + g01 * z + u0,
    beta1j = g10 + g11 * z + u1,
    y = beta0j + beta1j * x1 + eij)
}
)

simdata <- replicate(n = 1000,
                     expr = tibble(
                       id = 1:(ngroup * groupsize),
                       group = rep(1:ngroup, each = groupsize),
                       eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
                       x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
                       x2 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3),
                       z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u2 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u0 = rep(rnorm(n = ngroup,
                                      mean = 0,
                                      sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, g20, u2, eij, icc) {

                                        daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var((g20 + u2)*x2) + var(eij))))

                                        return(daticc)
                                      }, interval = c(0, 100),
                                      tol = .0001,
                                      extendInt = 'yes',
                                      maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, g20 = g20, u2 = u2, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
                       beta0j = g00 + g01 * z + u0,
                       beta1j = g10 + g11 * z + u1,
                       beta2j = g20 + u2,
                       y = beta0j + beta1j * x1 + beta2j * x2 + eij),
                     simplify = FALSE)
#####
# Multivariate data simulation
#####
simdata <- replicate(n = 1000,
                     expr = rmvnorm(n = (ngroup*groupsize),
                                    mean = c(0,0),
                                    sigma = matrix(c(6.25, 0,
                                                     0, 9), ncol = 2)) %>%
                       as_tibble() %>%
                       rename(x1 = V1, eij = V2) %>%
                       mutate(id = 1:(ngroup*groupsize),
                              group = rep(1:ngroup, each = groupsize)) %>%
                       cbind(., rep(rmvnorm(ngroup,
                                            mean = c(0, 0, 0),
                                            sigma = matrix(data = c(1, 0, 0,
                                                                    0, .49, 0,
                                                                    0, 0, 20),
                                                           ncol = 3)), each = groupsize) %>% matrix(., ncol = 3)) %>%
                       mutate(z = `1`,
                              u1 = `2`,
                              u0 = `3`,
                              beta0j = g00 + g01*z + u0,
                              beta1j = g10 + g11*z + u1,
                              y = beta0j + beta1j*x1 + rnorm(n = ngroup * groupsize, mean = 0, sd = var(beta0j)*(1/icc-1)-var(beta1j*x1))) %>%
                       select(-`1`, -`2`, -`3`),
                     simplify = FALSE)


#####
# Using a constraint for data generation to make sure the ICC is .5. Process takes very long
#####
library(foreach)
library(doParallel)

cores <- detectCores()
cl <- makeCluster(cores - 1)  # one less than total cores
registerDoParallel(cl)
simdata <- foreach(icount(10)) %dopar% {
  library(tidyverse)
  repeat {
    # Simulate data
    data <- rnorm(n = (ngroup * groupsize), mean = 0, sd = 2.5) %>%
      as_tibble() %>%
      rename(x1 = value) %>%
      mutate(
        id = 1:(ngroup * groupsize),
        group = rep(1:ngroup, each = groupsize),
        eij = rnorm((ngroup * groupsize), mean = 0, sd = 3),
        z = rep(rnorm(ngroup, mean = 0, sd = 1), each = groupsize),
        u1 = rep(rnorm(ngroup, mean = 0, sd = 0.7), each = groupsize),
        u0 = rep(rnorm(ngroup, mean = 0, sd = sqrt((-icc*(g01^2*var(z) + var(x1)*(g11^2*var(z) + var(u1)) + var(eij)) + g01^2*var(z)) / (icc-1))), each = groupsize),
        beta0j = g00 + g01 * z + u0,
        beta1j = g10 + g11 * z + u1,
        y = beta0j + beta1j * x1 + eij
      ) %>%
      select(-eij, -u1)

    # Calculate ICC
    iccdat <- var(data$beta0j) / var(data$y)

    # Check if ICC is close to the desired value
    if (abs(iccdat - icc) < 0.01) {
      break  # Exit the loop if ICC is close enough to the desired value
    } else {
      # Adjust the variance of u0 to get closer to the desired ICC
      data$u0 <- data$u0 * sqrt(icc / iccdat)
    }
  }
  data <- data %>% select(-u0)
  return(data)
}
stopCluster(cl)  # Stop the parallel processing

#####
# Non multivariate generation of the multilevel data set including manual equation for the variance of u0 to scale to the right ICC
#####
simdata <- replicate(n = 1000,
                     expr = rnorm(n = (ngroup*groupsize),
                                  mean = 0,
                                  sd = 2.5) %>%
                       as_tibble() %>%
                       rename(x1 = value) %>%
                       mutate(id = 1:(ngroup*groupsize),
                              group = rep(1:ngroup, each = groupsize),
                              eij = rnorm((ngroup*groupsize), mean = 0, sd = 3),
                              z = rep(rnorm(ngroup, mean = 0, sd = 1), each = groupsize),
                              u1 = rep(rnorm(ngroup, mean = 0, sd = .7), each = groupsize),
                              u0 = rep(rnorm(ngroup, mean = 0, sd = 1.79*sqrt((-icc*(g01^2*var(z) + var(x1)*(g11^2*var(z) + var(u1)) + var(eij)) + g01^2*var(z)) / (icc-1))), each = groupsize),
                              beta0j = g00 + g01*z + u0,
                              beta1j = g10 + g11*z + u1,
                              y = beta0j + beta1j*x1 + eij) %>%
                       select(-eij, -u1, -u0),
                     simplify = FALSE)
#####
# ICC with formula eij
#####
simdata <- lapply(1:1000, function(i) {
  data <- tibble(
    id = 1:(ngroup * groupsize),
    group = rep(1:ngroup, each = groupsize),
    x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
    z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u0 = rep(rnorm(n = ngroup, mean = 0, sd = 5), each = groupsize),
    beta0j = g00 + g01 * z + u0,
    beta1j = g10 + g11 * z + u1,
    y = beta0j + beta1j * x1 + rnorm(n = ngroup * groupsize, mean = 0, sd = sqrt(var(beta0j)*(1/icc-1)-var(beta1j*x1))))
}
)
#####
# ICC with uniroot eij
#####
simdata <- lapply(1:1000, function(i) {
  data <- tibble(
    id = 1:(ngroup * groupsize),
    group = rep(1:ngroup, each = groupsize),
    x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
    z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u0 = rep(rnorm(n = ngroup, mean = 0, sd = 20), each = groupsize),
    beta0j = g00 + g01 * z + u0,
    beta1j = g10 + g11 * z + u1,
    sdy = uniroot(function(vareij, beta0j, beta1j, x1, icc) {
      daticc <- (icc - (var(beta0j) / (var(beta0j)+var(beta1j*x1)+vareij)))
      return(daticc)
    }, interval = c(0, 100),
    tol = .0001,
    extendInt = 'yes',
    maxiter = 1000, beta0j = beta0j, beta1j = beta1j, x1 = x1,icc = icc)$root %>% as.numeric() %>% sqrt(),
    y = beta0j + beta1j * x1 + rnorm(n = ngroup * groupsize,
                                     mean = 0,
                                     sd = sdy))
}
)

#####
# Simulation data with only one x and z variable
#####
simdata <- lapply(1:1000, function(i) {
  data <- tibble(
    id = 1:(ngroup * groupsize),
    group = rep(1:ngroup, each = groupsize),
    eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
    x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
    z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u0 = rep(rnorm(n = ngroup,
                   mean = 0,
                   sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, eij, icc) {

                     daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var(eij))))

                     return(daticc)
                   }, interval = c(0, 100),
                   tol = .0001,
                   extendInt = 'yes',
                   maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
    beta0j = g00 + g01 * z + u0,
    beta1j = g10 + g11 * z + u1,
    y = beta0j + beta1j * x1 + eij)
}
)
#####
# Simulation like in thesis proposal
#####
simdata <- lapply(1:1000, function(i) {
  data <- tibble(
    id = 1:(ngroup * groupsize),
    group = rep(1:ngroup, each = groupsize),
    eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
    x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
    x2 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3),
    z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u2 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u0 = rep(rnorm(n = ngroup,
                   mean = 0,
                   sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, g20, u2, eij, icc) {

                     daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var((g20 + u2)*x2) + var(eij))))

                     return(daticc)
                   }, interval = c(0, 100),
                   tol = .0001,
                   extendInt = 'yes',
                   maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, g20 = g20, u2 = u2, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
    beta0j = g00 + g01 * z + u0,
    beta1j = g10 + g11 * z + u1,
    y = beta0j + beta1j * x1 + eij)
}
)

simdata <- replicate(n = 1000,
                     expr = tibble(
                       id = 1:(ngroup * groupsize),
                       group = rep(1:ngroup, each = groupsize),
                       eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
                       x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
                       x2 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3),
                       z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u2 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u0 = rep(rnorm(n = ngroup,
                                      mean = 0,
                                      sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, g20, u2, eij, icc) {

                                        daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var((g20 + u2)*x2) + var(eij))))

                                        return(daticc)
                                      }, interval = c(0, 100),
                                      tol = .0001,
                                      extendInt = 'yes',
                                      maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, g20 = g20, u2 = u2, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
                       beta0j = g00 + g01 * z + u0,
                       beta1j = g10 + g11 * z + u1,
                       beta2j = g20 + u2,
                       y = beta0j + beta1j * x1 + beta2j * x2 + eij),
                     simplify = FALSE)
