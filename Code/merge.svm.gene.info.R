#This function 
merge.svm.gene.info <- function(results.dir = ".", gene.info.table, 
	entrezgene.column = "entrezgene_id", gene.start.column = "start_position",
	gene.end.column = "end_position", gene.name.column = "external_gene_name"){
	
	module.dir.info <- get.module.dir(results.dir, dir.table = TRUE)
	module.dir <- module.dir.info$module.dir
	dir.table <- module.dir.info$dir.table

	all.results <- vector(mode = "list", length = length(module.dir))
	names(all.results) <- dir.table[,2]
	
	for(i in 1:length(module.dir)){
		
		svm.csv.file <- file.path(module.dir[i], "Candidate.Gene.SVM.Scores.csv")
		fp.csv.file <- file.path(module.dir[i], "Candidate.Gene.FP.Rates.csv")
		results.file <- file.path(module.dir[i], "Candidate.Gene.Results.csv")
		plot.file <- file.path(module.dir[i], "Candidate.Gene.SVM.Results.jpg")
		
		svm.scores <- read.csv(svm.csv.file, stringsAsFactors = FALSE)
		mean.svm <- colMeans(svm.scores)
		gene.ids <- gsub("X", "", colnames(svm.scores))
		
		fp.rates <- read.csv(fp.csv.file, stringsAsFactors = FALSE)
		mean.fp <- colMeans(fp.rates)

		common.genes <- intersect(as.numeric(gene.info.table[,entrezgene.column]), gene.ids)
		gene.locale.table <- match(common.genes, as.numeric(gene.info.table[,entrezgene.column]))
		# head(cbind(common.genes, gene.info.table[gene.locale.table,]))
		gene.locale.svm <- match(common.genes, gene.ids)
		# head(cbind(common.genes, gene.ids[gene.locale.svm]))
		
		final.table <- cbind(gene.info.table[gene.locale.table,], 
		mean.svm[gene.locale.svm], mean.fp[gene.locale.svm])
		ncol.final.table <- ncol(final.table)
		colnames(final.table)[tail(1:ncol.final.table, 2)] <- c("Mean.SVM.Score", "Mean.FP.Rate")
		
		write.table(final.table, results.file, sep = ",", quote = FALSE, row.names = FALSE)
		
		mean.gene.position <- rowMeans(apply(final.table[,c(gene.start.column, gene.end.column)], 2, as.numeric))
		
		jpeg(plot.file, height = 7, width = 10, units = "in", res = 300)
		plot.new()
		plot.window(xlim = c(min(mean.gene.position), max(mean.gene.position)), 
			ylim = c(min(as.numeric(final.table[,"Mean.SVM.Score"])), max(as.numeric(final.table[,"Mean.SVM.Score"]))))
		text(x = as.numeric(mean.gene.position), y = as.numeric(final.table[,"Mean.SVM.Score"]), 
			labels = final.table[,gene.name.column], cex = 0.7)
		axis(1);axis(2)
		mtext("Genomic Position", side = 1, line = 2.5)				
		mtext("SVM Score", side = 2, line = 2.5)
		abline(h = 0)
		dev.off()
		
		all.results[[i]] <- final.table
		}
	invisible(all.results)

}
