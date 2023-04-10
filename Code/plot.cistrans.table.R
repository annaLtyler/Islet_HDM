#This function takes in an eqtl.table with positions
# of eQTLs for each transcript with columns
#transcript.id, chromosome, position, LOD score

#and a transcript.pos.table, for the physical position
#of each transcript with columns
#transcript.id, chromosome, and position
	#positions should be in Mb

#col can be a vector of the same length as the
#number of rows in transcript.pos.table.
#if add is FALSE a new plot will be generated.
#if add is TRUE, the points will be plotted on 
#an existing plot.


plot.cistrans.table <- function(eqtl.table, transcript.pos.table, map, col = NULL, 
	add = FALSE, cex = 0.3, show.cross.chromosome.boundaries = TRUE, label.cex = 1){
	
	if(class(eqtl.table)[1] != "matrix"){stop("It looks as if eqtl.table might be a data.frame or a tibble. Please convert to a matrix.")}
	if(class(transcript.pos.table)[1] != "matrix"){stop("It looks as if transcript.pos.table might be a data.frame or a tibble. Please convert to a matrix.")}

    #we will plot one point per eQTL. Make a color vector
    #based on the entries in the eQTL table
	if(is.null(col)){col = rep("black", nrow(eqtl.table))}
	if(length(col) == 1){col <- rep(col, nrow(eqtl.table))}

	#make relative positions for each SNP by chromosome
	 #get the max position for making relative positions of transcripts
	chr.max <- sapply(map, function(x) max(x))
	chr.sum <- sum(chr.max)
	
	chr.max.coord <- cumsum(chr.max)
	chr.min.coord <- chr.max.coord - chr.max

	if(!add){
		plot.max <- sum(chr.max)
		label.pos <- plot.max*-0.02

		plot.new()
		plot.window(xlim = c(0, plot.max), ylim = c(0, plot.max))

		#add chromosome boundaries and labels
		par(xpd = TRUE)
		#draw a rectangle around the whole plot
		draw.rectangle(0, plot.max, 0, plot.max, border = "darkgray")

		for(i in 1:length(map)){
			chr.mean.x <- mean(c(chr.max.coord[i], chr.min.coord[i]))
			if(i %% 2 == 1){
				#QTL location rectangles
				draw.rectangle(chr.min.coord[i], chr.max.coord[i], 0, plot.max, 
				fill = rgb(189/256 ,189/256 ,189/256, 
				alpha = 0.5), border = NA)

				if(show.cross.chromosome.boundaries){
					#transcript position rectangles
					draw.rectangle(min.y = chr.min.coord[i], max.y = chr.max.coord[i], 
					min.x = 0, max.x = plot.max, fill = rgb(189/256 ,189/256 ,189/256, 
					alpha = 0.5), border = NA)
				}
				
				#lines for transcript position chromosome boundaries
				#segments(x0 = 0, x1 = plot.max, y0 = chr.max.coord[i], lty = 2, col = "darkgray")
				#segments(x0 = 0, x1 = plot.max, y0 = chr.min.coord[i], lty = 2, col = "darkgray")

			}
			text(x = chr.mean.x, y = label.pos, labels = names(map)[i], cex = label.cex)
			text(y = chr.mean.x, x = label.pos, labels = names(map)[i], cex = label.cex)

		}
	
		par(xpd = FALSE)
    }
    
	#for each eQTL, figure out its relative position on the x axis
    rel.pos <- function(chr, pos){
        chr.locale <- which(names(map) == chr)
        rel.loc <- pos/chr.max[chr.locale]
        chr.len <- (chr.max.coord[chr.locale] - chr.min.coord[chr.locale])
        adjust.pos <- chr.len * rel.loc
        rel.x <-  chr.min.coord[chr.locale] + adjust.pos
        return(rel.x)
    }

	#find the relative position of each eQTL
    eqtl.x <- apply(eqtl.table, 1, function(x) rel.pos(x[2], as.numeric(x[3])))
    
    #build a table to hold the transcription 
	#position for each transcript in the eQTL table
	eqtl.idx <- match(eqtl.table[,1], transcript.pos.table[,1])
	eqtl.transcripts <- transcript.pos.table[eqtl.idx,]
	
	#get relative positions for these too.
	transcript.y <- apply(eqtl.transcripts, 1, function(x) rel.pos(x[2], as.numeric(x[3])))
    
    #we won't get positions for anything on chromosomes Y or MT
	has.position <- which(sapply(transcript.y, length) > 0)
    eqtl.x <- eqtl.x[has.position]
    transcript.y <- unlist(transcript.y[has.position])
    all.col <- col[has.position]

    #plot(eqtl.x, transcript.y, col = all.col, pch = 16, cex = cex)
    points(eqtl.x, transcript.y, col = all.col, pch = 16, cex = cex)

	coord.table <- cbind(eqtl.x, transcript.y, all.col)
    colnames(coord.table) <- c("x", "y", "col")



	invisible(coord.table)
}
