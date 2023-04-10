#This function gets the mean phenotype effects for
#two continuous markers


get.pheno.effects.continuous <- function(data.obj, geno.obj, marker1.name, marker2.name, 
pheno.name, covar = NULL, scan.what = c("normalized.traits", "raw.traits", "eigentraits"),
n.geno.bins = 50){

	pheno <- get_pheno(data.obj, scan_what = scan.what, covar = covar)
	geno <- get_geno(data.obj, geno.obj)

	covar.info <- get_covar(data.obj)
	if(!is.null(covar.info$covar_table)){
		not.na.locale <- which(!is.na(rowSums(covar.info$covar_table)))
	}else{
		not.na.locale <- 1:nrow(pheno)	
	}
			
	split.m1 <- strsplit(marker1.name, "_")
	m1.name <- sapply(split.m1, function(x) x[1])
	m1.allele <- sapply(split.m1, function(x) x[2])

	split.m2 <- strsplit(marker2.name, "_")
	m2.name <- sapply(split.m2, function(x) x[1])
	m2.allele <- sapply(split.m2, function(x) x[2])

	m1.geno <- geno[not.na.locale,which(dimnames(geno)[[2]] == m1.allele), which(dimnames(geno)[[3]] == m1.name)]
	m2.geno <- geno[not.na.locale,which(dimnames(geno)[[2]] == m2.allele), which(dimnames(geno)[[3]] == m2.name)]
	
	pheno.locale <- which(colnames(pheno) == pheno.name)
	if(length(pheno.locale) == 0){
		stop("Can't find phenotype, check to make sure scan.what is set properly.")
	}
	pheno.vals <- pheno[,pheno.locale]
	
	#pheno.col <- colors.from.values(pheno.vals, use.pheatmap.colors = TRUE)
	#plot.with.model(m1.geno, m2.geno, col = pheno.col, pch = 16)
	#add.model <- lm(pheno.vals~m1.geno+m2.geno);summary(add.model)
	#int.model <- lm(pheno.vals~m1.geno*m2.geno);summary(int.model)
	#pred.add <- predict(add.model)
	#pred.int <- predict(int.model)
	#plot.with.model(pheno.vals, pred.add)
	#plot.with.model(pheno.vals, pred.int)
	#add.diff <- pheno.vals - pred.add
	#int.diff <- pheno.vals - pred.int


	pred.mats <- int_heat(phenoV = pheno.vals, marker1_vals = m1.geno, 
	marker2_vals = m2.geno, pheno_name = pheno.name, marker1_label = m1.name, 
	marker2_label = m2.name, nbins = n.geno.bins)

	#pheatmap(pred.mats[[1]], cluster_rows = FALSE, cluster_cols = FALSE)
	#pheatmap(pred.mats[[2]], cluster_rows = FALSE, cluster_cols = FALSE)
	#pheatmap(pred.mats[[3]], cluster_rows = FALSE, cluster_cols = FALSE)

	return(pred.mats)
	}
