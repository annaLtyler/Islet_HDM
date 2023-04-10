mouse_to_human_ensembl  <- function(mouse_ensembl = NULL, human_ensembl = NULL, hum.mus.ortho){
    
    if(!is.null(mouse_ensembl)){
        id.locale <- match(mouse_ensembl, hum.mus.ortho[,"Mouse.Ortholog.Ensembl"])
    }else{
        id.locale <- match(human_ensembl, hum.mus.ortho[,"Human.Ensembl"])
    }
    alt.ensembl <- hum.mus.ortho[id.locale,]
    
    return(alt.ensembl)
}