source("cross-validation.R")
library(ggplot2)

load("data/gray-5000.RData")


#pdf("gray-5000-ando-zhang.pdf")

small.featureset.size <- 100
n.targets <- 10


drug.targets <- read.table("data/gray-targets.tab", fill = TRUE)


random.problem.list <- sample(small.featureset.size, n.targets)
h <- 100
iters <- 3
lambda <- 1
for(i in 1:3) {
	drug <- colnames(Y)[[i]]
	targets.for.drug <- as.character(drug.targets[i,2])
	targets.for.drug <- unlist(strsplit(targets.for.drug, split = ","))
	print(targets.for.drug)
	title <- drug
	title <- paste(title, " - Numerical solution with 5000 features", sep = "")
	graph <- score.targets(X, Y, Xu, drug, random.problem.list, targets.for.drug, graph.title = title)
	print(graph)
	
}



#dev.off()
