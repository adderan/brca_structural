source("../multitask/multitask.R")
source("cross-validation.R")
load("data/gray-5000.RData")
library(reshape)
library(ggplot2)

drug.targets <- read.table("data/gray-targets.tab", fill = TRUE)

n.drugs <- dim(Y)[[2]]

scores <- matrix(NA, n.drugs, 2)

h <- 100
lambda <- 1
iters <- 3

colnames(scores) <- c("ASO with drug target as auxiliary problem", "glmnet")
rownames(scores) <- colnames(Y)
for(i in 1:n.drugs) {
	drug <- colnames(Y)[[i]]
	drug.target <- drug.targets[drug.targets$V1 == drug,]
	drug.target <- as.character(drug.target[,2])
	drug.target <- unlist(strsplit(drug.target, split = ","))
	for(i in 1:length(drug.target)) {
		if(!is.null(drug.target[i]) && drug.target[i] %in% rownames(X)) {
			drug.target <- drug.target[i]
			break
		}
	}
	if(is.null(drug.target) || !(drug.target %in% rownames(X))) next

	scores[drug,1] <- cross.validation.aso(X, Y, Xu, drug, drug.target, h = h, iters = iters, lambda = lambda)

	scores[drug,2] <- cross.validation.glmnet(X, Y, drug, lambda)


}
scores <- na.omit(scores)
scores <- melt(scores)
colnames(scores) <- c("drug", "description", "score")

png(filename="gray-5000-aso.png", width = 900, height = 700)
graph <- ggplot(data = scores, aes(x = drug, y = score, fill = description)) + geom_bar(stat = "identity", position = position_dodge()) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(graph)
dev.off()
#ggsave(graph, file = "gray-5000-aso.png", dpi = 500)

