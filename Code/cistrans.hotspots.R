#This function takes the same input as 
#plot.cistrans.table and plots 

#eqtl.table contains positions
# of eQTLs for each transcript with columns:
#transcript.id, chromosome, position, LOD score

#transcript.pos.table contains the physical position
#of each transcript with columns:
#transcript.id, chromosome, and position
	#positions should be in Mb

#mb buffer dictates how far from the gene 
#boundaries an eQTL is before it is considered
#distal. The default is 4 Mb.

#This function returns lists of genes that are 
#targeted by each window, split by whether they
#are local or distal targets.

cistrans.hotspots <- function(eqtl.table, transcript.pos.table, map, 
    mb.buffer = 4, window.mb = 4, window.gap.mb = 2, col = NULL, 
    plot.results = TRUE, hotspot.min = 20, ymax = NULL){

    if(class(eqtl.table)[1] != "matrix"){stop("It looks as if eqtl.table might be a data.frame or a tibble. Please convert to a matrix.")}
	if(class(transcript.pos.table)[1] != "matrix"){stop("It looks as if transcript.pos.table might be a data.frame or a tibble. Please convert to a matrix.")}
	if(is.null(col)){col = "black"}

    local.counts <- local.names <- distal.counts <- distal.names <- vector(mode = "list", length = length(map))
    names(local.counts) <- names(local.names) <- names(distal.counts) <- names(distal.names) <- names(map)

    #go through eQTL on each chromosome
    for(ch in 1:length(map)){
        chr.name <- names(map)[ch]
        eqtl.chr.idx <- which(eqtl.table[,2] == chr.name)
        chr.eqtl.table <- eqtl.table[eqtl.chr.idx,]
        
        #classify each eQTL as cis or trans
        #start with everything as cis. 
        #convert to trans any eQTL that is (a) on a different chromosome
        #or (b) more than 4 mb away from gene body.
        matched.gene.table <- transcript.pos.table[match(chr.eqtl.table[,1], transcript.pos.table[,1]),]
        cis.trans <- rep("cis", nrow(chr.eqtl.table))
        cis.trans[which(matched.gene.table[,2] != chr.name)] <- "trans"
        eqtl.dist <- as.numeric(chr.eqtl.table[,3]) - as.numeric(matched.gene.table[,3])
        cis.trans[which(eqtl.dist > mb.buffer)] <- "trans"

        #separate the cis and trans eQTL and count the eQTL in a sliding window.
        windows <- sliding.window.el(1:max(map[[ch]]), window.mb, window.gap.mb)
        cis.pos <- as.numeric(chr.eqtl.table[which(cis.trans == "cis"),3])
        cis.idx <- lapply(windows, function(x) intersect(which(cis.pos >= min(x)), which(cis.pos <= max(x))))
        cis.counts <- sapply(cis.idx, length)
        cis.names <- lapply(cis.idx, function(x) chr.eqtl.table[which(cis.trans == "cis")[x],"gene.id"])

        #plot(cis.counts, type = "l")
        names(cis.counts) <- names(cis.names) <- sapply(windows, mean)        
        local.counts[[ch]] <- cis.counts
        local.names[[ch]] <- cis.names
        
        trans.pos <- as.numeric(chr.eqtl.table[which(cis.trans == "trans"),3])
        trans.idx <- lapply(windows, function(x) intersect(which(trans.pos >= min(x)), which(trans.pos <= max(x))))
        trans.counts <- sapply(trans.idx, length)
        trans.names <- lapply(trans.idx, function(x) chr.eqtl.table[which(cis.trans == "trans")[x],"gene.id"])
        names(trans.counts) <- names(trans.names) <- sapply(windows, mean)
        #plot(trans.counts, type = "l")
        distal.counts[[ch]] <- trans.counts   
        distal.names[[ch]] <- trans.names
    }

    if(plot.results){

        plot_counts <- function(eqtl.counts, plot.label = "eQTL Counts", line.height = NA){
            layout(matrix(1:(length(map)+1), nrow = 1), widths = c(0.9, rep(1, length(map))))
            if(is.null(ymax)){
                ylim = c(0, max(unlist(eqtl.counts)))
            }else{
                ylim = c(0, ymax)
            }
            par(mar = c(4,0,4,0))
            plot.new()
            plot.window(xlim = c(0,1), ylim = ylim)
            tick.y <- axTicks(side = 2)
            segments(x0 = 0.9, y0 = 0, y1 = max(tick.y))
            segments(x0 = rep(0.9, length(tick.y)), x1 = rep(0.8, length(tick.y)), 
                y0 = tick.y)
            text(x = rep(0.75, length(tick.y)), y = tick.y, labels = tick.y, adj = 1)
            text(x = 0.15, y = mean(ylim), labels = "Count", cex = 1.5, srt = 90)
            for(ch in 1:length(map)){
                xlim <- c(min(map[[ch]]), max(map[[ch]]))
                plot.new()
                plot.window(xlim = xlim, ylim = ylim)
                if(ch%%2 == 0){
                    draw.rectangle(min(map[[ch]]), max(map[[ch]]), ylim[1], max(tick.y), 
                        fill = "lightgray", border = NA)
                }
                abline(h = tick.y, col = rgb(0.8, 0.8, 0.8)) #add grid lines
                points(as.numeric(names(eqtl.counts[[ch]])), eqtl.counts[[ch]], type = "l", 
                    col = col)
                abline(h = 0)
                abline(h = line.height, col = "red")
                par(xpd = NA)
                text(x = mean(map[[ch]]), y = ylim[2]*-0.05, labels = names(map)[ch])
                par(xpd = TRUE)
            }        
            mtext(plot.label, side = 3, line = -2.5, outer = TRUE)
        }

        #pdf("~/Desktop/eQTL_counts.pdf", width = 14, height = 3)
        plot_counts(local.counts, "Local eQTL Counts")
        plot_counts(distal.counts, "Distal eQTL Counts", line.height = hotspot.min)
        #dev.off()

    }

    results <- list("local_targets" = local.names, "distal_targets" = distal.names)
    invisible(results)

}
