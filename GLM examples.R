# example of a Poisson GLM
# --------------------------

library(gtools) # for logit calculations
# kernel bandwith
# https://cran.r-project.org/web/packages/kedd/vignettes/kedd.pdf

# number of points
N = 1000

# accident data
accident.df <- data.frame(id = 1:N,
                          danger.level = runif(N, min = 1, max = 3),
                          slope = runif(N, min = 0, max = 45),
                          orientation = runif(N, min = 0, max = 360),
                          accident = NA)

# risk model
p.acc <- function{
  
}