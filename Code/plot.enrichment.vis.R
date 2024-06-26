#This function takes in an enrichment table returned from gProfileR
#and *plots* the results in an easier-to-read table
#the "gprofiler" option of order.by uses the default ordering from
#gprofiler
#This function has been updated to 

plot.enrichment.vis <- function(enrichment, num.terms = 10, text.size = 1, 
order.by = c("p_value", "gprofiler", "overlap_size", "term_size"), 
decreasing = FALSE, plot.label = "Enrichment", highlight.terms = NULL, 
highlight.col = "#1f78b4", mar = c(5,20,4,4), max.term.size = NULL){
	
	if(class(enrichment) == "list"){
		enrichment <- enrichment[[1]]
	}

	if(is.null(enrichment) || nrow(enrichment) == 0){
		plot.new()
		plot.window(xlim = c(0, 1), ylim = c(0, 1))
		text(x = 0.5, y = 0.5, "No Significant Enrichment")
		text(x = 0.5, y = 0.75, plot.label)
		return()
		}
		
	order.by <- order.by[1]
	
	if(!is.null(max.term.size)){
		small.terms <- which(enrichment[,"term_size"] <= max.term.size)
		enrichment <- enrichment[small.terms,,drop=FALSE]
	}

	if(order.by != "gprofiler"){
		enrichment <- enrichment[order(enrichment[,order.by], decreasing = decreasing),]
		}	
	
	#subset to only the terms we are interested in
	n.row <- min(c(num.terms, nrow(enrichment)))
	enrichment <- enrichment[1:n.row,]
	
	log.p <- -log(enrichment[,"p_value"], base = 10)
	
	col <- rep("gray", nrow(enrichment))
	
	if(!is.null(highlight.terms)){
		highlight.idx <- unique(unlist(sapply(highlight.terms, function(x) grep(x, enrichment[,"term_name"], ignore.case = TRUE))))
		if(length(highlight.idx) > 0){
			col[highlight.idx] <- highlight.col
			}
		}

	par(mar = mar)
	barplot(rev(log.p), horiz = TRUE, xlab = "-log10 p value", 
	names = rev(enrichment[,"term_name"]), las = 2, main = plot.label, col = col)
	
	term.table <- cbind(enrichment, col)
	invisible(term.table)
	
	}	

