source("../multitask/multitask.R")
library(glmnet)
library(ggplot2)


cross.validation.aso <- function(X, Y, Xu, drug, target, h, iters, lambda, analytic) {

		
	#X <- X[feature.start:(n.features+feature.start),]
	#Xu <- Xu[feature.start:(n.features+feature.start),]

	labels <- Y[,drug]
	good.samples <- !is.na(labels)
	X.fixed <- X[, good.samples]
	labels <- labels[good.samples]

	n.samples <- length(labels)
	split <- floor(n.samples/2)

	labels.test <- labels[1:split]
	labels.train <- labels[(split + 1):n.samples]

	cat("Using ", target, " to predict ", drug, "\n")

	target.index <- which(rownames(X) == target)
		
	#Y.aux <- X[target,]
	X.reduced <- X.fixed[-target.index,]
		
	Y.aux <- Xu[target,]
	Xu.reduced <- Xu[-target.index,]


	X.test <- X.reduced[,1:split]
	X.train <- X.reduced[, (split+1):n.samples]


		
	ando.X <- list()
	ando.Y <- list()
	ando.X[[drug]] <- X.train
	ando.X[[target]] <- Xu.reduced
	ando.Y[[drug]] <- labels.train
	ando.Y[[target]] <- Y.aux


	fit <- aso.train(ando.X, ando.Y, h, iters, lambda = lambda, analytic = analytic, use.cache = FALSE)
	predictions <- aso.predict(aso.trained.model = fit, new.x = X.test, primary.problem = drug)

	score <- score.model(predictions, labels.test)
	return(score)

	#ando.X <- list()
	#ando.Y <- list()
	#X.train <- X.fixed[,1:split]
	#X.test <- X.fixed[,(split + 1):n.samples]
	#ando.X[[drug]] <- X.train
	#ando.Y[[drug]] <- labels.train
	#fit <- aso.train(ando.X, ando.Y, h = h, iters = iters, lambda = lambda, analytic = analytic, use.cache = use.cache)
	#predictions <- aso.predict(aso.trained.model = fit, new.x = X.test, primary.problem = drug)
	#scores[["None"]] <- mse(predictions, labels.test)
	#return(list(scores = scores, models = models))

}
cross.validation.glmnet <- function(X, Y, drug, lambda) {
	labels <- Y[,drug]
	good.samples <- !is.na(labels)
	X.fixed <- X[,good.samples]
	labels <- labels[good.samples]

	n.samples <- length(labels)
	split <- floor(n.samples/2)
	X.test <- X[,1:split]
	X.train <- X[,(split+1):n.samples]
	labels.test <- labels[1:split]
	labels.train <- labels[(split+1):n.samples]


	fit <- glmnet(t(X.train), labels.train, alpha = 0)
	predictions <- predict(fit, newx = t(X.test), s = lambda)

	score <- score.model(predictions, labels.test)
	return(score)
}
score.model <- function(predictions, answers) {
	mse <- mean((predictions - answers)^2)
	score <- cor(as.vector(predictions), answers, method = "spearman")
	return(score)
}

score.targets <- function(X, Y, Xu, drug, target.list, drug.targets, graph.title = "default", h, iters, lambda, analytic) {
	#target.list <- sample(dim(X)[[1]], n.targets)
	plot.data <- data.frame(gene = character(), drug = character(), score = numeric(), description = character(), stringsAsFactors = FALSE)
	print(drug)
	n.targets <- length(target.list)
	for(i in 1:n.targets) {
		target <- target.list[[i]]
		target.name <- rownames(X)[[target]]
		target.score <- cross.validation.aso(X, Y, Xu, drug, target.name, h = h, iters = iters, lambda = lambda, analytic = analytic)
		plot.data <- rbind(plot.data, data.frame(gene = target.name, drug = drug, score = target.score, description = "ASO with randomly selected auxiliary problem"))		
	}
	for(i in 1:length(drug.targets)) {
		drug.target <- drug.targets[[i]]
		if(!(drug.target %in% rownames(X))) {
			next
		}
		primary.target.score <- cross.validation.aso(X, Y, Xu, drug, drug.target, h = h, iters = iters, lambda = lambda, analytic = analytic)
		plot.data <- rbind(plot.data, data.frame(gene = drug.target, drug = drug, score = primary.target.score, description = "ASO with drug target as auxiliary problem"))
	}



	glmnet.score <- cross.validation.glmnet(X, Y, drug, lambda)

	plot.data <- rbind(plot.data, data.frame(gene = "glmnet", drug = drug, score = glmnet.score, description = "glmnet only"))

	graph <- ggplot(data = plot.data, aes(x = gene, y = score, fill = description)) + geom_bar(stat="identity") + labs(title = graph.title)
	print(graph)
}
