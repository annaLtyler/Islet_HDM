#perform CCA over a grid of penalties
#calculate R2 between latent traits
#and raw traits

opt_r2 <- function(tissue.mats, 
num.components = 2, pen_seq = seq(0.1, 0.9, 0.1),
search_grid = FALSE, filename){

    if(!file.exists(filename)){
        if(search_grid){
            trans.pen.seq <- pair.matrix(pen_seq, ordered = TRUE, self.pairs = TRUE)
        }else{
            trans.pen.seq <- cbind(pen_seq, pen_seq)
        }
        colnames(trans.pen.seq) <- c("X", "Z")
        all.r2 <- matrix(NA, nrow = nrow(trans.pen.seq), ncol = ncol(tissue.mats$Z))
        
        for(i in 1:nrow(trans.pen.seq)){
            if(is.interactive){
                report.progress(i, nrow(trans.pen.seq))
            }
            tissue.results <- CCA(tissue.mats$X, tissue.mats$Z, 
            typex = "standard", typez = "standard", K = num.components, 
            penaltyx = trans.pen.seq[i,1], penaltyz = trans.pen.seq[i,2], 
            niter = 100, trace = FALSE)
            lt.r2 <- get_latent_var(tissue.mats, tissue.results)
            all.r2[i,] <- lt.r2[[3]]
         }
        result <- list("Penalties" = trans.pen.seq, 
        "All_R2" = all.r2)
        saveRDS(result, filename)
    }else{
        result <- readRDS(filename)
    }
    return(result)
}