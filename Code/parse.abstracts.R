#This function attempts to pull the abstract out of a list
#returned by fetch_pubmed_data


parse.abstracts <- function(abstract.list, keyword){

    section.boundaries <- which(abstract.list == "")
    separate.sections <- vector(mode = "list", length = length(section.boundaries))

    start.idx <- 1
    for(i in 1:length(section.boundaries)){
        end.idx <- section.boundaries[i]
        separate.sections[[i]] <- abstract.list[start.idx:end.idx]
        start.idx <- end.idx + 1
    }

    keyword.locale <- grep(keyword, separate.sections, ignore.case = TRUE)
    relevant.sections <- separate.sections[keyword.locale]
    return(relevant.sections)
}