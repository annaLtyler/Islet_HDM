#This function takes in a gene set, from 

set_to_list <- function(gene.set, set.name = NULL){
  all.sets <- gene.set[[3]]
  all.genes <- gene.set[[5]]
  if(!is.null(set.name)){
    set.locale <- grep(set.name, all.sets)
  }else{
    set.locale <- 1:length(all.sets)
  }
  set.types <- unique(all.sets[set.locale])
  set.genes <- vector(mode = "list", length = length(set.types))
  names(set.genes) <- set.types
  for(i in 1:length(set.genes)){
    set.locale <- which(all.sets == set.types[i])  
    set.genes[[i]] <- all.genes[set.locale]
  }
  return(set.genes)
}