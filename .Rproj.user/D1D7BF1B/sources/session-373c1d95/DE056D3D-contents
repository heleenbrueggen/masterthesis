#Stan4bart

# Installing package
# devtools::install_github("vdorie/stan4bart")

# Libraries
library(stan4bart)
library(BART)
library(lme4)
library(dbarts)


(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
summary(fm1)
stan4bart(y ~ 1 + x1 + z + x1*z + (x1|group), data = simdata[[1]])


fit <- rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2, group.by = group, data = simdata[[1]])
summary(fit)
