#This function gets the mean phenotype effects for
#two markers and their combined effect using three
#different patterns of genotype coding
#To use a binning other than the ones offered by 
#additive, dominant, and recessive, use geno.bins
#to define custom binning


get.pheno.effects <- function(data.obj, geno.obj, marker1.name, marker2.name, 
pheno.name, covar = NULL, ref.center = TRUE,
scan.what = c("normalized.traits", "raw.traits", "eigentraits"), 
collapsed.net = FALSE, geno.coding = c("Additive", "Dominant", "Recessive"), 
geno.bins = NULL){

	pheno <- get_pheno(data.obj, scan_what = scan.what, covar = covar)
	geno <- get_geno(data.obj, geno.obj)

	covar.info <- get_covar(data.obj)
	if(!is.null(covar.info$covar_table)){
		not.na.locale <- which(!is.na(rowSums(covar.info$covar_table)))
	}else{
		not.na.locale <- 1:nrow(pheno)	
	}

	#=========================================================
	#internal functions
	#=========================================================
	get.marker.name <- function(marker.name){
		if(length(grep("Chr", marker.name)) > 0){
			marker.locale <- which(names(data.obj$linkage_blocks_full) == marker.name)
			marker.label <- data.obj$linkage_blocks_full[[marker.locale]]
			return(marker.label)
			}else{
			if(is.null(allele)){
				full.name <- marker.name
				}else{
				full.name <- paste0(marker.name, "_", allele)
				}
			return(marker.name)
			}
		}
	#=========================================================


	if(collapsed.net){
		m1.markers <- data.obj$linkage_blocks_collapsed[[which(names(data.obj$linkage_blocks_collapsed) == marker1.name)]]
		split.m1 <- strsplit(m1.markers, "_")
		m1.name <- sapply(split.m1, function(x) x[1])
		m1.allele <- sapply(split.m1, function(x) x[2])

		m2.markers <- data.obj$linkage_blocks_collapsed[[which(names(data.obj$linkage_blocks_collapsed) == marker2.name)]]
		split.m2 <- strsplit(m1.markers, "_")
		m2.name <- sapply(split.m2, function(x) x[1])
		m2.allele <- sapply(split.m2, function(x) x[2])
		
		m1.marker.ind <- match(m1.name, dimnames(geno)[[3]])
		m1.allele.ind <- match(m1.allele, dimnames(geno)[[2]])
		
		m2.marker.ind <- match(m2.name, dimnames(geno)[[3]])
		m2.allele.ind <- match(m2.allele, dimnames(geno)[[2]])

		all.int <- data.obj$full_net[m1.marker.ind, m2.marker.ind, drop=FALSE]
		max.int <- which(abs(all.int) == max(abs(all.int)), arr.ind = TRUE)

		m1.geno <- geno[not.na.locale,m1.allele.ind[max.int[1,1]], m1.marker.ind[max.int[1,1]]]
		m2.geno <- geno[not.na.locale,m2.allele.ind[max.int[1,1]], m2.marker.ind[max.int[1,2]]]

		m1.name <- 	m1.markers[max.int[1,1]]
		m2.name <- 	m2.markers[max.int[1,2]]
		}else{
			
		m1.marker <- get.marker.name(marker1.name)
		split.m1 <- strsplit(m1.marker, "_")
		m1.name <- sapply(split.m1, function(x) x[1])
		m1.allele <- sapply(split.m1, function(x) x[2])

		m2.marker <- get.marker.name(marker2.name)
		split.m2 <- strsplit(m2.marker, "_")
		m2.name <- sapply(split.m2, function(x) x[1])
		m2.allele <- sapply(split.m2, function(x) x[2])

		m1.geno <- geno[not.na.locale,which(dimnames(geno)[[2]] == m1.allele), which(dimnames(geno)[[3]] == m1.name)]
		m2.geno <- geno[not.na.locale,which(dimnames(geno)[[2]] == m2.allele), which(dimnames(geno)[[3]] == m2.name)]
		}

	
	if(is.null(geno.bins)){
		if(geno.coding == "Dominant"){
			geno.bins <- c(0, 1)
			}
			
		if(geno.coding == "Recessive"){
			geno.bins <- c(0, 1)
			}
	
		if(geno.coding == "Additive"){
			geno.bins <- c(0, 0.5, 1)
		}
	}

	binned.m1 <- bin.vector(m1.geno, geno.bins)
	binned.m2 <- bin.vector(m2.geno, geno.bins)
	#plot(jitter(binned.m1), jitter(binned.m2))

	pheno.locale <- which(colnames(pheno) == pheno.name)
	if(length(pheno.locale) == 0){
		stop("Can't find phenotype, check to make sure scan.what is set properly.")
	}
	pheno.vals <- pheno[,pheno.locale]

	geno.pairs <- pair.matrix(geno.bins, ordered = TRUE, self.pairs = TRUE)
	
	pair.locale <- apply(geno.pairs, 1, 
	function(x) intersect(which(binned.m1 == x[1]), which(binned.m2 == x[2])))
	
	pair.pheno <- sapply(pair.locale, function(x) mean(pheno.vals[x], na.rm = TRUE))

	#subtract the phenotype of the referenct genotype for both 
	#genotypes. This is phenotype for the minimum genotype at
	#both markers. 

	ref.locale <- min(which(!is.na(pair.pheno)))

	if(!is.na(pair.pheno[ref.locale]) && ref.center){
		geno.pheno <- cbind(geno.pairs, pair.pheno - pair.pheno[ref.locale])
	}else{
		geno.pheno <- cbind(geno.pairs, pair.pheno)
	}
	colnames(geno.pheno) <- c("marker1", "marker2", "avg.pheno")

	int.mat <- matrix(NA, nrow = length(geno.bins), ncol = length(geno.bins))
	rownames(int.mat) <- geno.bins
	colnames(int.mat) <- geno.bins
	for(i in 1:nrow(geno.pheno)){
		m1.locale <- which(rownames(int.mat) == geno.pheno[i,1])
		m2.locale <- which(colnames(int.mat) == geno.pheno[i,2])
		int.mat[m1.locale, m2.locale] <- geno.pheno[i,3]
	}
	#pheatmap(int.mat, cluster_rows = FALSE, cluster_cols = FALSE)

	return(list(geno.pheno, int.mat))
	}
