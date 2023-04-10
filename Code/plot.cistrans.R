#This function takes in an R/qtl2 scan1 object, a map,
#and a gene information matrix, and plots the cis/trans 
#eQTL plot

plot.cistrans <- function(scan1.result, map, gene.info, lod.threshold = 3, 
id.col = "ensembl_gene_id", chr.column = "chromosome_name", pos.column = "start_position"){
	
	#find max lod scores for pt colors
	max.lod <- max(as.vector(scan1.result))
	val.seq <- seq(0, max.lod, 0.01)
	pt.cols <- colors.from.values(val.seq, use.pheatmap.colors = TRUE)


	#make relative positions for each SNP by chromosome
	chr.max <- sapply(map, function(x) max(x)) #get the max position for making relative positions of transcripts
	rel.snp <- sapply(map, function(x) x/max(x))
	for(i in 1:length(rel.snp)){
		rel.snp[[i]] <- rel.snp[[i]] + i
	}
	snp.pos.table <- Reduce("c", rel.snp)

	rel.transcript.pos <- function(transcript.id){
		chr <- unique(gene.info[which(gene.info[,id.col] == transcript.id),chr.column])
		if(length(chr) == 0){
			return(NA)
		}else{
			chr.locale <- which(names(chr.max) == chr)
			chr.size <- chr.max[chr.locale]
			trans.loc <- as.numeric(gene.info[which(gene.info[,id.col] == transcript.id),pos.column])/1e6
			rel.loc <- (trans.loc/chr.size) + chr.locale
			return(rel.loc[1])
		}
	}

	#find where each transcript in the object is encoded on the genome
	#these values are used for the x position of each point
	transcript.pos <- sapply(colnames(scan1.result), rel.transcript.pos)
	
	#find the lod scores above the threshold
	lod.thresh <- apply(scan1.result, 2, function(x) x[which(x >= lod.threshold)])
	#find the position of lod scores above the threshold
	lod.thresh.idx <- apply(scan1.result, 2, function(x) which(x >= lod.threshold))

	#now plot each lod score where
	#x values are the position of the eQTL peak (from lod.thresh.idx)
	#y values are the position of the transcript for which the peak was found (from gene.loc)	
	

	#for each transcript, at the y value where it is encoded, plot
	#points at the positions of SNPs with high lod scores for that transcript
	#I'm still confused about x and y. Need to review this carefully

	plot.new()
	plot.window(xlim = c(0, 20), ylim = c(0, 20))

	#add chromosome boundaries and labels
	par(xpd = TRUE)
	for(i in seq(1,length(map), 2)){
		draw.rectangle(i, i+1, 1, length(map)+1, fill = rgb(189/256 ,189/256 ,189/256, alpha = 0.5),
		border = NA)
		text(x = i+0.5, y = 0.5, labels = i)
		#if(i < length(map)){text(x = i+1.5, y = 0.5, labels = i+1)}
	}


	for(i in 1:length(transcript.pos)){
		tp <- transcript.pos[[i]]
		if(length(tp) > 0 && !is.na(tp)){
			related.snps <- names(lod.thresh[[i]])
			if(length(related.snps) > 0){
				sp <- snp.pos.table[match(related.snps, names(snp.pos.table))]
				pt.size <- lod.thresh[[i]]
				pt.col <- colors.from.values(pt.size, col.scale = "gray", global.color.scale = TRUE,
				global.min = 0, global.max = max.lod, light.dark = "d")
				points(sp, rep(tp, length(sp)), pch = 16, cex = 0.3, col = pt.col)
			}
		}
	}

	
}