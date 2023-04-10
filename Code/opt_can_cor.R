#This function searches a grid, or series of penalties,
#and calculates R2 across all latent variable pairs
#The default number of components to fit is the number
#of traits. 
opt_can_cor <- function(tissue.mats, num.components = ncol(scaled.pheno), 
pen_seq = seq(0.1, 0.9, 0.1), search_grid = FALSE, filename){

    if(!file.exists(filename)){
        if(search_grid){
            trans.pen.seq <- pair.matrix(pen_seq, ordered = TRUE, self.pairs = TRUE)
        }else{
            trans.pen.seq <- cbind(pen_seq, pen_seq)
        }
        colnames(trans.pen.seq) <- c("X", "Z")
        all.can.cor <- matrix(NA, nrow = nrow(trans.pen.seq), ncol = ncol(tissue.mats$Z))
        rownames(all.can.cor) <- trans.pen.seq[,1]
        for(i in 1:nrow(trans.pen.seq)){
            if(is.interactive){
                report.progress(i, nrow(trans.pen.seq))
            }
            tissue.results <- CCA(tissue.mats$X, tissue.mats$Z, 
            typex = "standard", typez = "standard", K = num.components, 
            penaltyx = trans.pen.seq[i,1], penaltyz = trans.pen.seq[i,2], 
            niter = 100, trace = FALSE)
            all.can.cor[i,] <- tissue.results$cor
         }
        result <- list("Penalties" = trans.pen.seq, 
        "All_Can_Cor" = all.can.cor)
        saveRDS(result, filename)
    }else{
        result <- readRDS(filename)
    }
    return(result)
}