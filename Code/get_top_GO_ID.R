get_top_GO_ID <- function(enrichment.results, max_pval = 0.05, max_term_size = NULL,
sort.by = "p_value", num.terms = 10){

    result.table <- enrichment.results$result
    sig.vals <- which(result.table[,"p_value"] <= max_pval)
    if(length(sig.vals) == 0){
        return("no significant enrichments")
    }
    
    sig.table <- result.table[sig.vals,]
    if(sort.by == "p_value"){
        sig.table <- sig.table[order(sig.table[,"p_value"], decreasing = FALSE),]
    }

    if(!is.null(max_term_size)){
        to.keep <- which(sig.table[,"term_size"] <= max_term_size)
        if(length(to.keep) == 0){return("No terms smaller than maximum size.")}
        sig.table <- sig.table[to.keep,,drop=FALSE]
    }

    to.return <- min(c(num.terms), nrow(sig.table))

    top.terms <- sig.table[1:to.return,c("term_id", "p_value")]
    return(top.terms)
}
