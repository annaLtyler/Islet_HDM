#use any ID type to get back a mouse-level ID
#if the input type is known, specify input type.
#this will speed up the translation.

get_mouse_info <- function(input.ID, manifest, input.type = NULL){
    if(is.null(input.type)){
        id.col <- which.max(apply(manifest, 2, function(x) length(which(x == input.ID))))
    }else{
        id.col <- which(colnames(manifest) == input.type)
    }
    id.locale <- which(manifest[,id.col] == input.ID)[1]
    mouse.info <- manifest[id.locale,c("Barcode.1", "CLIMB.ID", "Sex", "Diet", "Treatment", "User.Defined.Strain", "Timepoint")]
    mouse.info[which(mouse.info == "10% fat + fiber")] <- "LFD"
    mouse.info[which(mouse.info == "44% fat + fiber")] <- "HFD"
  
    return(mouse.info)
}