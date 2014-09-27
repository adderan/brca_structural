library(reshape)
#plot score as a function of h


drug.list <- c("Erlotinib", "Lapatinib", "Sorafenib")

n <- 3
drug.targets <- c("EGFR", "EGFR", "EGFR")
h.scores <- matrix(0, n, length(drug.list) + 1)
colnames(h.scores) <- c("h", drug.list)
for(i in 1:n) {
	h.scores[i,1] <- h
	for(j in 1:length(drug.list)) {
		drug <- drug.list[[j]]
		drug.target <- drug.targets[[j]] 
		score <- cross.validation.aso(X[100:200,], Y, Xu[100:200,], drug, drug.target, h = h, iters = iters, lambda = lambda)
		h.scores[i,j+1] <- score
	}
	h <- h + 10

}
print(h.scores)
h.scores <- data.frame(h.scores)
h.scores <- melt(h.scores, id = c('h'))
print(h.scores)

graph <- ggplot(h.scores, aes(x = h, y = value, fill = variable)) + geom_line()
print(graph)

