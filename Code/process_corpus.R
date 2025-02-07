process_corpus <- function(text.list, remove.words = stopwords("english"), 
	top.terms = 25, topic.level = "###", file.name = "abstracts.pdf"){

	has.vals <- which(sapply(text.list, function(x) x != ""))
	text.list <- text.list[has.vals]

	corp <- VCorpus(VectorSource(unlist(text.list)))
	corp <- tm_map(corp, content_transformer(tolower))
	corp <- tm_map(corp, removeWords, remove.words)
	corp <- tm_map(corp, stripWhitespace)

	#sometimes removing prohibited words reduces titles to basically nothing
	#e.g. "Forever 21" in a Klotho search.
	#remove anything with fewer than three characters here
	has.vals <- which(sapply(corp, function(x) nchar(x$content)) > 10)
	corp <- corp[has.vals]

	dtm <- DocumentTermMatrix(corp)
	widf <- weightTfIdf(dtm)

	#high.term.idx <- which(term.max > 1)
	#length(high.term.idx)
	
	doc.cor <- cor(t(as.matrix(widf)))
	#pheatmap(doc.cor, show_rownames = FALSE, show_colnames = FALSE)
	max.k <- min(c(50, nrow(doc.cor) - 2))
	print(max.k)

	#cluster the document correlation matrix into groups
	doc.cl.test <- test.pam.k(doc.cor, 2:max.k, plot.results = FALSE)
	#boxplot(doc.cl.test[[1]])
	#plot(sapply(doc.cl.test[[1]], mean))
	best.k <- which.max(sapply(doc.cl.test[[1]], mean))
	k <- as.numeric(names(best.k))

	doc.cl <- pam(doc.cor, k = k)

	#row.order <- order(doc.cl$clustering)
	#pheatmap(doc.cor[row.order,row.order], show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE)

	#pdf(file.name, width = 14, height = 7)

	num.docs <- sapply(1:k, function(x) length(which(doc.cl$clustering == x)))
	names(num.docs) <- paste0("Cluster", 1:k)
	#barplot(sort(num.docs), las = 2)

	word.importance.list <- vector(mode = "list", length = k)

	for(cl in 1:k){
		cat(paste(topic.level, "Cluster", cl, "\n\n"))
		#par(mfrow = c(2,2))
		par(mfrow = c(1,2))
		cl.idx <- which(doc.cl$clustering == cl)
		cl.idf <- widf[cl.idx,]
		
		par(mar = c(4,4,4,2))
		#cl.decomp <- plot.decomp(t(cl.idf), label.points = TRUE, main = paste("Cluster", cl, "\n", num.docs[cl], "Documents"))
		
		cl.word.importance <- apply(cl.idf, 2, mean) #use the mean for the cluster as word importance
		word.importance.list[[cl]] <- cl.word.importance
		#I have tried other word importance measures, and they all favor individual 
		#articles. I really want the importance overall for the cluster
		#mean gives us terms that are important acruss multiple dosuments
		#as well as some that are very highly important in individual articles
		#cl.word.importance <- apply(cl.idf, 2, function(x) median(x[which(x > 0)])) #use the median of the non-zero terms for the cluster as word importance
		#cl.word.importance <- abs(cl.decomp$u[,1]) #use first PC as word importance
		#cl.word.importance <- abs(cl.decomp$u[,2]) #use second PC as word importance; this is better than the first
		#names(cl.word.importance) <- colnames(widf)

		#plot(sort(cl.word.importance), main = "Word Importance", ylab = "Importance")
		
		top.words <- tail(sort(cl.word.importance), top.terms)
		par(mar = c(4,12,2,2))
		mean.order <- order(top.words)
		boxplot(as.matrix(cl.idf)[,names(top.words)[mean.order]], las = 2, horizontal = TRUE,
			main = paste("Cluster", cl, "\n", num.docs[cl], "Documents"), xlab = "term importance (tf-idf)")

		#barplot(top.words, horiz = TRUE, las = 2)
		#hist(word.importance[which(word.importance > 0)])
		#plot(sort(word.importance[which(word.importance > 0)]))
		#par(mar = c(4,8,2,2)); barplot(tail(sort(word.importance), 10), las = 2, horiz = TRUE)
		par(mar = c(0,0,0,0))
		wordcloud(names(top.words), top.words)
		cat("\n\n")
	}
	#dev.off()
	
	result <- list("widf" = widf, "num_docs" = num.docs, 
		"clusters" = doc.cl$clustering, 
		"word_importance" = word.importance.list)
	invisible(result)
}
