# ---------------------------------------
# Effect of Pseudoreplication on a t-test
# ---------------------------------------

set.seed(123)

# number of points
N <- 20
# data set with N points without and N points with a treatment, coded 0 and 1
# the effect is 0.5 and random numbers or summed with mean zero and stdev 1
data <- data.frame(treatment = rep(c(0, 1), each = N),
                   effect = rep(c(0, 0.5), each = N) + rnorm(n = N, mean = 0, sd = 1))

boxplot(effect ~ treatment, data, xlab = "Treatment", ylab = "Effect")

# unpaired t-test between the non-treated and treated subjects
tt1 <- t.test(x = data$effect[data$treatment == 0], y = data$effect[data$treatment == 1],
              paired = FALSE, alternative = "two.sided")

# double the data by copying the data set below itself
data.rep <- rbind(data, data)

tt2 <- t.test(x = data.rep$effect[data.rep$treatment == 0], y = data.rep$effect[data.rep$treatment == 1],
              paired = FALSE, alternative = "two.sided")

print(paste("Doubling the data changed the t-tests p-value from", 
            round(tt1$p.value, 5), "to", round(tt2$p.value, 5), "."))


# lm1 <- lm(formula = effect ~ treatment, data)
# summary(lm1)
# 
# lm2 <- lm(formula = effect ~ treatment, data.rep)
# summary(lm2)
