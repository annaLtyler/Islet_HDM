#This function finds specific motifs in the gene interaction network

find.motifs <- function(data.obj, p_or_q = 0.05, include.covar = FALSE){


	var.inf <- plot_variant_influences(cross, covar_width = 1, pheno_width = 1,
		show_marker_labels =TRUE, p_or_q = p_or_q)
		dev.off()
	just.int <- var.inf[,1:nrow(var.inf)]
	just.main <- var.inf[,(nrow(var.inf)+1):ncol(var.inf)]
	num.pheno <- ncol(just.main)
	
	if(!include.covar){
		covar.info <- get_covar(data.obj)
		covar.names <- covar.info$covar_names
		covar.locale <- match(covar.names, rownames(just.main))
		if(length(covar.locale) > 0 && !is.na(covar.locale)){
			just.int <- just.int[-covar.locale,-covar.locale]
			just.main <- just.main[-covar.locale,]
			}
		}
		
	#find all the edges
	gene.edge.locale <- which(just.int != 0, arr.ind = TRUE)
	
	if(length(gene.edge.locale) == 0){stop("There are no edges at this p value.")}
	
	#to find the complete three-node motifs, look at each
	#edge and verify that each gene has a connection to the
	#phenotype as well. Each phenotype will get its own set
	#of motifs
	all.triplets <- vector(mode = "list", length = num.pheno)
	names(all.triplets) <- colnames(just.main)
	
	all.triplet.signs <- vector(mode = "list", length = num.pheno)
	names(all.triplet.signs) <- colnames(just.main)
	
	for(ph in 1:num.pheno){
		triplet.list <- matrix(NA, ncol = 3, nrow = nrow(gene.edge.locale))
		triplet.signs <- matrix(NA, ncol = 3, nrow = nrow(gene.edge.locale))
		colnames(triplet.list) <- c("source", "target", "pheno")
		colnames(triplet.signs) <- c("source-target", "source-Ph", "target-Ph")
		for(g in 1:nrow(gene.edge.locale)){
			source.locale <- gene.edge.locale[g,1]
			target.locale <- gene.edge.locale[g,2]
			
			source.main <- just.main[source.locale,ph]
			target.main <- just.main[target.locale,ph]
			
			#motifs must have two main effects
			if(!is.na(source.main) && !is.na(target.main)){
				triplet.names <- c(colnames(just.int)[source.locale], 
				colnames(just.int)[target.locale], colnames(just.main)[ph])
				triplet.list[g,] <- triplet.names
				
				int.val <- just.int[source.locale,target.locale]
				trip.signs <- sign(c(source.main, target.main, int.val))
				triplet.signs[g,] <- trip.signs
				}
			}
		#take out the NA rows, and the non-interactions 
		#for each phenotype
		non.na.locale <- which(!is.na(triplet.list[,1]))
		triplet.list <- triplet.list[non.na.locale,,drop=FALSE]
		triplet.signs <- triplet.signs[non.na.locale,,drop=FALSE]

		all.triplets[[ph]] <- triplet.list	
		all.triplet.signs[[ph]] <- triplet.signs
	}
	
	result <- list("triplets" = all.triplets, "triplet-signs" = all.triplet.signs)
	return(result)		
	
	
}