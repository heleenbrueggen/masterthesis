##########
# Rommel #
##########
uniroot(varu0,
        interval = c(0, 100),
        tol = .0001,
        extendInt = 'yes',
        maxiter = 1000,
        g00=10, g01=.5, g11=.35, g02=.5, g21=.35, g32=.35,
        g10=.2, g20=.2, g30=.2, g40=.2, g50=.2, g60=.2, g70=.2,
        varz1=3, varz2=2,
        varx1=6, varx2=12, varx3=6.5, varx4=3.3, varx5=8, varx6=15, varx7=20,
        varu1=1, varu2=1, varu3=1, varu4=1, varu5=1, varu6=1,
        icc=.5)
uniroot(varu0, # Original
        interval = c(0, 100),
        tol = .0001,
        extendInt = 'yes',
        maxiter = 1000,
        g00=g00, g01=g01, g10=g10, g11=g11, g02=g02, g20=g20, g21=g21, g30=g30, g32=g32, g40=g40, g50=g50, g60=g60, g70=g70,
        varz1=varz1, varz2=varz2,
        varx1=varx1, varx2=varx2, varx3=varx3, varx4=varx4, varx5=varx5, varx6=varx6, varx7=varx7,
        varu1=varu1, varu2=varu2, varu3=varu3, varu4=varu4, varu5=varu5, varu6=varu6,
        vare=vare,
        icc=icc)$root %>% sqrt()

varu0 <- function(varu0,
                  g00, g10, g01, g11, g02, g20, g21, g30, g32, g40, g50, g60, g70,
                  z1, z2,
                  x1, x2, x3, x4, x5, x6, x7,
                  u1, u2, u3, u4, u5, u6,
                  eij,
                  icc) {
  daticc <- icc - ((g00^2 + g01^2 * var(z1) + g02^2 * var(z2) + varu0) / (g00^2 + g01^2 * var(z1) + g02^2 * var(z2) + varu0 +
                                                                             (g10^2 + g11^2 * var(z1) + var(u1))*var(x1) +
                                                                             (g20^2 + g21^2 * var(z1) + var(u2))*var(x2) +
                                                                             (g30^2 + g32^2 * var(z2) + var(u3))*var(x3) +
                                                                             (g40^2 + var(u4))*var(x4) +
                                                                             (g50^2 + var(u5))*var(x5) +
                                                                             (g60^2 + var(u6))*var(x6) +
                                                                             g70^2*var(x7) + var(eij)))
  return(daticc)
}
varu0 <- function(varu0, # Original
                  g00, g10, g01, g11, g02, g20, g21, g30, g32, g40, g50, g60, g70,
                  z1, z2,
                  x1, x2, x3, x4, x5, x6, x7,
                  u1, u2, u3, u4, u5, u6,
                  eij,
                  icc) {
  daticc <- icc - ((var(g00 + g01 * z1 + g02 * z2) + varu0) / (var(g00 + g01 * z1 + g02 * z2) + varu0 +
                                                                 var((g10 + g11 * z1 + u1)*x1) +
                                                                 var((g20 + g21 * z1 + u2)*x2) +
                                                                 var((g30 + g32 * z2 + u3)*x3) +
                                                                 var((g40 + u4)*x4) +
                                                                 var((g50 + u5)*x5) +
                                                                 var((g60 + u6)*x6) +
                                                                 var(g70*x7) + var(eij)))
  return(daticc)
}
varu0e <- function(varu0,
                   g00, g01, g11, g02, g21, g32,
                   g10, g20, g30, g40, g50, g60, g70,
                   varz1, varz2,
                   varx1, varx2, varx3, varx4, varx5, varx6, varx7,
                   varu1, varu2, varu3, varu4, varu5, varu6,
                   vare,
                   icc) {
  daticc <- icc - ((g00^2 + g01^2 * varz1 + g02^2 * varz2 + varu0) / (g00^2 + g01^2 * varz1 + g02^2 * varz2 + varu0 +
                                                                        (g10^2 + g11^2 * varz1 + varu1)*varx1 +
                                                                        (g20^2 + g21^2 * varz1 + varu2)*varx2 +
                                                                        (g30^2 + g32^2 * varz2 + varu3)*varx3 +
                                                                        (g40^2 + varu4)*varx4 +
                                                                        (g50^2 + varu5)*varx5 +
                                                                        (g60^2 + varu6)*varx6 +
                                                                        g70^2*varx7 + vare))
  return(daticc)
}


varu0e <- function(parameters, g00, g01, g11, g02, g21, g32, g10, g20, g30, g40, g50, g60, g70, varz1, varz2, varx1, varx2, varx3, varx4, varx5, varx6, varx7, varu1, varu2, varu3, varu4, varu5, varu6, icc) {
  varu0 <- parameters[1]
  vare <- parameters[2]

  daticc <- icc - ((g00^2 + g01^2 * varz1 + g02^2 * varz2 + varu0) / (g00^2 + g01^2 * varz1 + g02^2 * varz2 + varu0 +
                                                                        (g10^2 + g11^2 * varz1 + varu1)*varx1 +
                                                                        (g20^2 + g21^2 * varz1 + varu2)*varx2 +
                                                                        (g30^2 + g32^2 * varz2 + varu3)*varx3 +
                                                                        (g40^2 + varu4)*varx4 +
                                                                        (g50^2 + varu5)*varx5 +
                                                                        (g60^2 + varu6)*varx6 +
                                                                        g70^2*varx7 + vare))
  return(abs(daticc))  # Return the absolute value of the result for minimization
}

optim(c(varu0 = 10, vare = 25), varu0e,
      g00=10, g01=.5, g11=.35, g02=.5, g21=.35, g32=.35,
      g10=.8, g20=.8, g30=.8, g40=.8, g50=.8, g60=.8, g70=.8,
      varz1=3, varz2=2,
      varx1=6, varx2=12, varx3=6.5, varx4=3.3, varx5=8, varx6=15, varx7=20,
      varu1=1, varu2=1, varu3=1, varu4=1, varu5=1, varu6=1, # pass other parameters as needed
      method = "L-BFGS-B",
      lower = c(0,25), upper = c(Inf, Inf),
      icc = .5)$par

simdata <- simulationdata(nsim = 1, ngroup = 50, groupsize = 50, icc = .5,
                          g00 = 10, # Overall intercept
                          g10 = .8, g20 = .8, g30 = .8, g40 = .8, g50 = .8, g60 = .8, g70 = .8, # Individual effects
                          g01 = .5, g02 = .5, # Group level effects
                          g11 = .35, g21 = .35, g32 = .35, # Cross level interactions
                          vare = 25, # Residual variance
                          varz1 = 3, varz2 = 2, # Variance of group level variables
                          varx1 = 6, varx2 = 12, varx3 = 6.5, varx4 = 3.3, varx5 = 8, varx6 = 15, varx7 = 20) # Variance of individual level variables
lmer(y ~ 1 + (1|group), REML = FALSE, data = simdata[[1]]) %>% summ()
#####
# ICC function
#####
# Generating function for ICC calculation
iccfunction <- function (data) {
  icc <- var(data$beta0j)/var(data$y)
  return(icc)
}
# Checking ICC over all simulated data sets
iccvalues <- rep(0, 288)
for(i in 1:288) {
  for (j in 1:10) {
    iccvalues[i] <- mean(iccfunction(simdatasets[[i]][[j]]))
  }
}
mean(iccvalue)
#####
# Simulated data set with function simulationdata
#####
simdata <- simulationdata(nsim = 10,
                   ngroup = 50,
                   groupsize = 50,
                   icc = .5,
                   g00 = 10,
                   g01 = .5, g02 = .5,
                   g11 = .35, g21 = .35, g32 = .35,
                   g10 = .8, g20 = .8, g30 = .8, g40 = .8, g50 = .8, g60 = .8, g70 = .8,
                   varz1 = 3, varz2 = 2,
                   varx1 = 6, varx2 = 12, varx3 = 6.5, varx4 = 3.3, varx5 = 8, varx6 = 15, varx7 = 20,
                   vare = 25)
lmer(y ~ 1 + (1|group), REML = FALSE, data = simdata[[1]]) %>% summ()
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
# Working Simulation with ICC in u0
#####
replicate(n = 1000,
          expr = tibble(
            # individual id
            id = 1:(ngroup * groupsize),
            # group id
            group = rep(1:ngroup, each = groupsize),
            # residual variance
            eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
            # level 1 variables
            x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
            x2 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3),
            x3 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2),
            x4 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3.4),
            x5 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2),
            x6 = rnorm(n = ngroup * groupsize, mean = 0, sd = 1.5),
            x7 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.3),
            # level 2 variables
            z1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
            z2 = rep(rnorm(n = ngroup, mean = 0, sd = 1.6), each = groupsize),
            # random slopes
            u6 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
            u5 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
            u4 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
            u3 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
            u2 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
            u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
            u0 = rep(rnorm(n = ngroup,
                           mean = 0,
                           sd = uniroot(function(varu0, g00, g01, g11, g20, g21, g30, g32, g40, g50, g60, g70, z1, z2, x1, x2, x3, x4, x5, x6, x7, eij, icc) {

                             daticc <- (icc - ((var(g00 + g01 * z1) + varu0) / (var(g00 + g01 * z1) + varu0 +
                                                                                  var((g10 + g11 * z1 + u1)*x1) +
                                                                                  var((g20 + g21 * z1 + u2)*x2) +
                                                                                  var((g30 + g32 * z2 + u3)*x3) +
                                                                                  var((g40 + u4)*x4) +
                                                                                  var((g50 + u5)*x5) +
                                                                                  var((g60 + u6)*x6) +
                                                                                  var(g70*x7) + var(eij))))

                             return(daticc)
                           }, interval = c(0, 100),
                           tol = .0001,
                           extendInt = 'yes',
                           maxiter = 1000, g00=g00, g01=g01, g11=g11, g20=g20, g21=g21, g30=g30, g32=g32, g40=g40, g50=g50, g60=g60, g70=g70, z1=z1, z2=z2, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, x6=x6, x7=x7, eij=eij, icc=icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
            # coefficient generation (including random slopes and cross-level interactions)
            beta0j = g00 + g01 * z1 + u0,
            beta1j = g10 + g11 * z1 + u1,
            beta2j = g20 + g21 * z1 + u2,
            beta3j = g30 + g32 * z2 + u3,
            beta4j = g40 + u4,
            beta5j = g50 + u5,
            beta6j = g60 + u6,
            beta7j = g70,
            # generation of dependent variable y
            y = beta0j + beta1j * x1 + beta2j * x2 + beta3j * x3 + beta4j * x4 + beta5j * x5 + beta6j * x6 * beta7j * x7 + eij) %>%
            # taking out terms that are only used for model generation
            select(-u0, -u1, -u2, -u3, -u4, -u5, -u6, -eij, -beta0j, -beta1j, -beta2j, -beta3j, -beta4j, -beta5j, -beta6j, -beta7j),
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
