source("../multitask/multitask.R")
library(glmnet)
load("data/gray-5000.RData")

predictions.as.auxiliary <- function(data, lambda) {
	X.train <- data$train
	X.test <- data$test
	labels <- data$labels
	X.unlabeled <- data$unlabeled

	labeled.glmnet.fit <- glmnet(x = t(X.train), y = labels, alpha = 0)
	extra.labels <- predict(labeled.glmnet.fit, t(X.unlabeled), s = lambda)

	base.predictions <- predict(labeled.glmnet.fit, t(X.test), s = lambda)

	aso.X <- list()
	aso.Y <- list()
	aso.X[["primary"]] <- X.train
	aso.X[["auxiliary"]] <- X.unlabeled
	aso.Y[["primary"]] <- labels
	aso.Y[["auxiliary"]] <- extra.labels

	aso.fit <- aso.train(aso.X, aso.Y)
	aso.predictions <- aso.predict(aso.fit, new.x = X.test, primary.problem = "primary")

	return(list(aso.predictions = aso.predictions, base.predictions = base.predictions))
}

split.data <- function(X, Y, Xu, drug) {
	labels <- Y[,drug]
	good.samples <- !is.na(labels)
	X <- X[,good.samples]
	labels <- labels[good.samples]

	n.samples <- length(labels)
	split <- floor(n.samples/2)
	X.test <- X[,1:split]
	X.train <- X[,(split+1):n.samples]

	labels.test <- labels[1:split]
	labels.train <- labels[(split+1):n.samples]


	return(list(train = X.train, test = X.test, labels = labels.train, answers = labels.test, unlabeled = Xu))
}
score.drug <- function(drug) {
	#load("data/gray-5000.RData")
	data <- split.data(X, Y, Xu, drug)
	predictions <- predictions.as.auxiliary(data, 1)
	aso.score <- cor(as.vector(predictions$aso.predictions), data$answers, method = "spearman")
	glm.score <- cor(as.vector(predictions$base.predictions), data$answers, method = "spearman")
	return(list(aso.score = aso.score, base.score = glm.score))
}
graph.all.drug.scores <- function() {
	library(ggplot2)
	library(reshape)
	n.drugs <- dim(Y)[[2]]
	scores <- matrix(0, n.drugs, 2)
	colnames(scores) <- c("ASO", "glmnet")
	rownames(scores) <- colnames(Y)
	for(i in 1:n.drugs) {
		drug <- colnames(Y)[[i]]
		drug.score <- score.drug(drug)
		scores[drug, 1] <- drug.score$aso.score
		scores[drug, 2] <- drug.score$base.score
	}
	scores <- melt(scores)
	return(scores)
}

		



