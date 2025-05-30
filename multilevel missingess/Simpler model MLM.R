#################
# Load packages #
#################
library(lme4)
library(lmerTest)
library(tidyverse)
library(magrittr)
library(mlmpower)
library(mvtnorm)
library(varTestnlme)
library(mice)


############
# Set seed #
############
set.seed(123)


###########################
# Defining design factors #
###########################
ngroups <- c(30, 50) # Number of groups
groupsizes <- c(15, 50) # Group sizes
iccs <- c(.5) # Intraclass correlation coefficient
mar_mcar <- c("mar", "mcar") # Missing data mechanism
miss <- c(0, 50) # Percentage of missing data
g <- c(.5) # Within-group effect size
combinations <- expand.grid(
    ngroup = ngroups,
    groupsize = groupsizes,
    icc = iccs,
    mar_mcar = mar_mcar,
    miss = miss,
    g = g
)


################
# Define model #
################
ranslope2 <- (
    effect_size(
        icc          = c(.50),
        within       = 0.2307692,
        between      = 0.2307692,
        product      = .01,
        random_slope = .03
    )
    + outcome("y", mean = 10, sd = 5)
        + within_predictor("x1", weight = .8, icc = 0)
        + within_predictor("x2", weight = .8, icc = 0)
        + within_predictor("x3", weight = .8)
        + within_predictor("x4", weight = .2)
        + between_predictor("z1", weight = .8)
        + between_predictor("z2", weight = .8)
        + random_slope("x1", weight = .5)
        + random_slope("x2", weight = .5)
        + product("x1", "z1", weight = .5)
        + correlations(
            within = fixed(.3),
            between = fixed(.3),
            randeff = 0
        )
)


#################
# Generate data #
#################
combinations <- rbind(read_rds(paste(path, "data/nomissing/combinations.rds", sep = "")), combinations)

simdata <- list()
names <- list()
for (i in seq_len(nrow(combinations))) {
    # Generate data
    simdata[[i]] <- generate(
        ranslope2,
        n_within = combinations$groupsize[i],
        n_between = combinations$ngroup[i],
        ndata = 100
    )
    simdata[[i]] <- map(simdata[[i]], \(x) setNames(x, c("group", "y", "x1", "x2", "x3", "x4", "z1", "z2")))

    # Generate names for each dataset
    names[[i]] <- paste("mlmsimdata2",
        colnames(combinations)[1], combinations[i, 1],
        colnames(combinations)[2], combinations[i, 2],
        colnames(combinations)[3], combinations[i, 3],
        colnames(combinations)[4], combinations[i, 4],
        colnames(combinations)[5], combinations[i, 5],
        colnames(combinations)[6], combinations[i, 6],
        sep = "_"
    )
}


#########################
# Creating missing data #
#########################
combinations <- combinations %>%
    filter(miss != 0)

# Defining patterns for missing mechanism
patterns <- expand.grid(c(1, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1)) %>%
    filter(rowSums(.) == 6 | rowSums(.) == 5 | rowSums(.) == 4 | rowSums(.) == 3) %>%
    as.matrix()
colnames(patterns) <- c("group", "y", "x1", "x2", "x3", "x4", "z1", "z2")
patterns <- patterns[sample(nrow(patterns), 10, replace = FALSE), ] # Sample only 10 patterns for all possibilities

# Determining the frequency of each pattern
freq <- ampute.default.freq(patterns)

# Determining the weights of each pattern
weights <- ampute.default.weights(patterns, "MAR")
colnames(weights) <- c("group", "y", "x1", "x2", "x3", "x4", "z1", "z2")
# Increasing weights for z2 and x4
weights[, "z2"] <- weights[, "z2"] * 1.5
weights[, "x3"] <- weights[, "x4"] * 2

# Creating missing data
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")

    if (combinations[i, "mar_mcar"] == "mcar") {
        simdata_miss[[i]] <-
            simdata[[i]] %>%
            map(
                ~ ampute(
                    data = .x,
                    prop = (combinations[i, "miss"] * .01),
                    mech = "MCAR",
                    patterns = patterns,
                    freq = freq,
                )$amp,
                .progress = TRUE
            )
    } else {
        simdata_miss[[i]] <-
            simdata[[i]] %>%
            map(
                ~ ampute(
                    data = .x,
                    prop = (combinations[i, "miss"] * .01),
                    mech = "MAR",
                    patterns = patterns,
                    freq = freq,
                    weights = weights
                )$amp,
                .progress = TRUE
            )
    }
}


############
# Analyses #
############


############################
# Plan parallel processing #
############################
cl <- makeForkCluster(5)


#######################
# Multilevel analysis #
#######################
# Define model
lmer.model <- function(x) {
    model <- x %>% lme4::lmer(y ~ x1 + x2 + x3 + x4 + z1 + z2 + x1 * z1 + (1 + x1 + x2 | group),
        REML = TRUE,
        data = .
    )
    results <- broom.mixed::tidy(model, conf.int = TRUE)

    return(list(results = results, model = model))
}
results_ldlc <- list()
for (i in seq_len(nrow(combinations))) {
    # Perform analyses
    results_ldlc[[i]] <- pblapply(simdata_miss[[i]], lmer.model.ld, cl = cl)
}
results_complete <- list()
for (i in seq_len(nrow(combinations))) {
    # Perform analyses
    results_complete[[i]] <- pblapply(simdata[[i]], lmer.model.complete, cl = cl)
}


############################
# Stop parallel processing #
############################
stopCluster(cl)


#####################################
# Creating functions for evaluation #
#####################################
# Bias function
bias.mlm <- function(estimated) {
    # Extracting estimates
    estimates <- estimated %>%
        map(function(x) {
            fixed <- x %>%
                filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "z1" | term == "z2" | term == "x1:z1") %>%
                pull(estimate, term)
            random <- x %>%
                filter(term == "sd__(Intercept)" | term == "sd__x1" | term == "sd__x2" | term == "sd__Observation") %>%
                pull(estimate, term) %>%
                sapply(function(x) x^2)
            c(fixed, random)
        })
    # Defining truth
    truth <- tibble(
        beta0j = 10,
        beta1j = 1.042100,
        beta2j = 1.042100,
        beta3j = 1.473752,
        beta4j = 0.368438,
        z1 = 1.679151,
        z2 = 1.679151,
        `x1:z1` = 0.500000,
        u0 = 1.227667,
        u1 = .375,
        u2 = .375,
        eij = 5.73077
    )
    # Calculating bias for all datasets
    bias.datasets <- map(estimates, \(x) (x - truth)) %>%
        list_rbind() %>%
        as_tibble()
    # Average bias
    bias <- colMeans(bias.datasets) %>%
        t() %>%
        as_tibble()
    # Mean of estimates
    mean_estimates <- estimates %>%
        do.call(rbind, .) %>%
        colMeans() %>%
        t() %>%
        as_tibble()
    # Calculating MCSE
    mcse <- map(estimates, \(x) (x - mean_estimates)^2) %>%
        list_rbind() %>%
        colSums() %>%
        map_vec(\(x) sqrt(x / (length(estimates) * (length(estimates) - 1)))) %>%
        t() %>%
        as_tibble()
    colnames(mcse) <- colnames(bias)

    return(list(bias.datasets = bias.datasets, bias = bias, bias.mcse = mcse))
}
# Coverage function
coverage.mlm <- function(estimated) {
    # Extracting estimates
    estimates <- estimated %>%
        map(~ .x %>%
            filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "z1" | term == "z2" | term == "x1:z1") %>%
            select(conf.low, conf.high))
    # Defining truth
    truth <- tibble(
        beta0j = 10,
        beta1j = 1.042100,
        beta2j = 1.042100,
        beta3j = 1.473752,
        beta4j = 0.368438,
        z1 = 1.679151,
        z2 = 1.679151,
        `x1:z1` = 0.500000
    ) %>%
        t() %>%
        as_tibble() %>%
        rename(value = V1)
    # Combining estimates and truth
    combined <- map(estimates, \(x) cbind(x, truth))
    # Coverage of all data sets
    coverage.datasets <- map(combined, \(x) x %>%
        mutate(coverage = value > conf.low & value < conf.high) %>%
        select(coverage) %>%
        mutate(coverage = as.numeric(coverage)) %>%
        t() %>%
        as_tibble() %>%
        rename(beta0j = V1, beta1j = V2, beta2j = V3, beta3j = V4, beta4j = V5, z1 = V6, z2 = V7, `x1:z1` = V8)) %>%
        list_rbind()
    # Average coverage
    coverage <- colMeans(coverage.datasets) %>%
        t() %>%
        as_tibble()
    # Calculating MCSE of coverage
    mcse <- sqrt(coverage * (1 - coverage) / length(estimates)) %>%
        as_tibble()

    return(list(coverage.datasets = coverage.datasets, coverage = coverage, coverage.mcse = mcse))
}
# Confidence interval width function
ciw.mlm <- function(estimated) {
    # Extracting estimates
    estimates <- estimated %>%
        map(~ .x %>%
            filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "z1" | term == "z2" | term == "x1:z1") %>%
            select(conf.low, conf.high))
    # Calculating CIW for all datasets
    ciw.datasets <- map(estimates, ~ .x %>%
        mutate(ciw = conf.high - conf.low) %>%
        select(ciw) %>%
        t() %>%
        as_tibble() %>%
        rename(beta0j = V1, beta1j = V2, beta2j = V3, beta3j = V4, beta4j = V5, z1 = V6, z2 = V7, `x1:z1` = V8)) %>%
        list_rbind()
    # Average CIW
    ciw <- colMeans(ciw.datasets) %>%
        t() %>%
        as_tibble()

    return(list(ciw.datasets = ciw.datasets, ciw = ciw))
}


###########################
# Formatting bias results #
###########################
bias.datasets_ldlc_mlm <- list()
bias.datasets_complete <- list()
for (i in seq_len(nrow(combinations))) {
    # Bias
    bias.datasets_ldlc_mlm[[i]] <- bias.mlm(analyseslc_mlm[[i]])
    bias.datasets_complete[[i]] <- bias.mlm(analyses_complete[[i]])
}

format.bias <- function(x, y, method, name = c("Bias", "MCSE")) {
    data <- cbind(x, y) %>%
        pivot_longer(cols = beta0j:eij, names_to = "term", values_to = name) %>%
        mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g), method = method)
    colnames(data) <- c("Number of groups", "Group size", "ICC", "Missingness mechanism", "Percentage of missing data", "Gamma", "Term", name, "Method")
    data <- data %>%
        mutate(`Missingness mechanism` = ifelse(`Missingness mechanism` == "mar", "MAR", "MCAR"))
    return(data)
}

bias_ldlc_mlm <- map(bias.datasets_ldlc_mlm, \(x) x$bias) %>%
    list_rbind() %>%
    as_tibble() %>%
    format.bias(combinations, method = "ldlc_mlm", name = "Bias") %>%
    cbind(map(bias.datasets_ldlc_mlm, \(x) x$bias.mcse) %>%
        list_rbind() %>%
        as_tibble() %>%
        format.bias(combinations, method = "ldlc_mlm", name = "MCSE") %>%
        dplyr::select(MCSE))

bias_complete <- map(bias.datasets_complete, \(x) x$bias) %>%
    list_rbind() %>%
    as_tibble() %>%
    format.bias(combinations, method = "complete", name = "Bias") %>%
    cbind(map(bias.datasets_complete, \(x) x$bias.mcse) %>%
        list_rbind() %>%
        as_tibble() %>%
        format.bias(combinations, method = "complete", name = "MCSE") %>%
        dplyr::select(MCSE))

bias_ld <- bind_rows(bias_ld, bias_ldlc_mlm, bias_complete)


##############
# Plots bias #
##############
# Define color palette
cbbPalette <- c("#000000", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#5D478B", "#0072B2", "#F0E442")

bias_ld %>%
    filter(Term == "eij" | Term == "u0") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .25) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # scale_x_continuous(n.breaks = 15, minor_breaks = seq(-11, 11, .5), limits = c(-10, 10)) +
    scale_y_discrete(limits = c("eij", "u0"), labels = c(expression(paste(epsilon, "ij")), expression(paste(upsilon, "0")))) +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 1.5, x = -2.5, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 1.5, x = 2.5, linetype = "dashed", color = "gray40") +
    geom_segment(y = 1.5, yend = 3, x = -8.474903, linetype = "dashed", color = "gray40") +
    geom_segment(y = 1.5, yend = 3, x = 8.474903, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )

bias_ld %>%
    filter(Term == "beta0j") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .3)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .3), width = .08) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.4, .4, .025)) +
    scale_y_discrete(limits = c("beta0j"), labels = c(expression(paste(gamma, "00")))) +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 2, x = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 2, x = .1, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )

bias_ld %>%
    filter(Term == "beta1j" | Term == "beta2j" | Term == "beta3j" | Term == "beta4j" | Term == "beta5j" | Term == "beta6j" | Term == "beta7j") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .3) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.08, .08, .005)) +
    scale_y_discrete(limits = c("beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j"), labels = c(expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")))) +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_vline(xintercept = -.05, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .05, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )

bias_ld %>%
    filter(Term == "z1" | Term == "z2") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .25) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.4, .4, .025)) +
    scale_y_discrete(limits = c("z1", "z2"), labels = c(expression(paste(gamma, "01")), expression(paste(gamma, "02")))) +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 3, x = -.05, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 3, x = .05, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )

bias_ld %>%
    filter(Term == "x1:z1" | Term == "x2:z1" | Term == "x3:z2") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .25) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.07, .07, .005)) +
    scale_y_discrete(limits = c("x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 4, x = -.035, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 4, x = .035, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )

bias_ld %>%
    filter(Term == "u1" | Term == "u2" | Term == "u3") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .25) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.1, .3, .01)) +
    scale_y_discrete(limits = c("u1", "u2", "u3"), labels = c(expression(paste(upsilon, "1")), expression(paste(upsilon, "2")), expression(paste(upsilon, "3")))) +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 4, x = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 4, x = .1, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
