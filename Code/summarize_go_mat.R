#go_term_mat is from plot.enrichment.group
#go_term_sim is from get_sim_mat in High_Dimentional_Mediation.Rmd
#if use.pc.mat is TRUE, we do an SVD of the similarity matrix
#before clustering.
summarize_go_mat <- function(go_term_mat, go_term_sim, min.cl = 2, 
    max.cl = 20, top.n.terms = 10, use.pc.mat = FALSE, pc = 2, plot.results = FALSE){

    #for all domains
    cl_by_domain <- domain_term_prob <- vector(mode = "list", length = length(go_term_sim))
    names(cl_by_domain) <- names(domain_term_prob) <- names(go_term_sim)

    for(d in 1:length(go_term_sim)){
    
        term_mat <- go_term_sim[[d]]
        common.terms <- intersect(names(rownames(go_term_mat)), names(rownames(term_mat)))

        term.idx <- match(common.terms, names(rownames(term_mat)))    
        sub_go_sim <- term_mat[term.idx, term.idx]

        term.idx <- match(common.terms, names(rownames(go_term_mat)))
        sub_go_term <- go_term_mat[term.idx,]
    
        #cluster the similarity matrix
        #use k means clustering, check for the best separation
        #and use that.
        sim_decomp <- plot.decomp(sub_go_sim, pc = pc, plot.results = FALSE)
        #pdf("~/Desktop/test.pdf", width = 10, height = 10); plot.decomp(sub_go_sim, pc = pc, label.points = TRUE); dev.off()

        if(use.pc.mat){
            k.test <- test.pam.k(sim_decomp$u, min.cl:max.cl, plot.results = FALSE)
        }else{
            k.test <- test.pam.k(sub_go_sim, min.cl:max.cl, plot.results = FALSE)
        }

        mean.cl <- sapply(k.test[[1]], mean)
        best_clustering <- which.max(mean.cl)
        best_cl <- k.test[[2]][,best_clustering]
        u_cl <- unique(best_cl)    
        #cluster membership
        cl_mem <- lapply(u_cl, function(x) which(best_cl == x))
        
        #plot PC with clusters colored
        if(plot.results){
            pairs(sim_decomp$u, col = best_cl, pch = 16)
        }

        #merge the rows in the go_term_mat based on the identified clusters
        cl_term <- lapply(cl_mem, function(x) go_term_mat[x,,drop=FALSE])
        
        #term_sig <- lapply(cl_term, function(x) sort(apply(x, 1, max), decreasing = TRUE))
        #pdf("~/Desktop/test.pdf", width = 10, height = 20)
        #par(mar = c(4,10,4,2))
        #barplot(sort(term_sig[[1]]), las = 2, horiz = TRUE, cex.names = 0.7)
        #pheatmap(cl_term[[1]][,mean.order], show_colnames = FALSE, cluster_cols = FALSE)
        #dev.off()
        
        #individual-level mean -log10(p) across whole term cluster
        ind_mean <- t(sapply(cl_term, function(x) apply(x, 2, function(y) mean(y))))
        
        #term-level mean -log10(p) across individuals
        term_mean <- lapply(cl_term, function(x) apply(x, 1, function(y) mean(y)))
        #i = 3
        #wordcloud(names(term_mean[[i]]), freq = term_mean[[i]])

        #make summary descriptions for each cluster 
        num_terms <- lapply(term_mean, function(x) min(c(length(x), top.n.terms)))
        cl_names <- sapply(1:length(term_mean), function(x) paste(names(sort(term_mean[[x]], decreasing = TRUE))[1:num_terms[[x]]], collapse = "; "))
        rownames(ind_mean) <- cl_names

        cl_by_domain[[d]] <- ind_mean
    }

    final_result <- Reduce("rbind", cl_by_domain)
    
    #pdf("~/Desktop/test.pdf", width = 40, height = 6)
    #pheatmap(final_result, show_colnames = FALSE, cluster_rows = FALSE)
    #dev.off()

    return(final_result)
}
