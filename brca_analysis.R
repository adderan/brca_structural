source("../multitask/multitask.R")
source("pair-cross-validation.R")
library(ggplot2)

load("data/gray-5000.RData")

#analytic.scores <- cross.validation(X, Y, Xu, drug = "Erlotinib", drug.target = "EGFR", n.tests = 5, n.features = 100, feature.start = 100, use.cache = FALSE, analytic = TRUE)

numerical.results <- cross.validation(X, Y, Xu, lambda = 2, drug = "Erlotinib", drug.target = "EGFR", n.tests = 5, n.features = 2000, feature.start = 0, use.cache = FALSE, analytic = FALSE)


#barplot(analytic.scores)
barplot(numerical.results$scores)


