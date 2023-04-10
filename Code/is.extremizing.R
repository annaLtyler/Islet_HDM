#This function determines whether a motif is extremizing or not
#based on predicted additive effects

is.extremizing  <- function(add.mat, act.mat){

    add.main1 <- colMeans(add.mat, na.rm = TRUE)
    #barplot(add.main1)
    add.main2 <- rowMeans(add.mat, na.rm = TRUE)
    #barplot(add.main2)

    main1.dir <- instant.slope(1:length(add.main1), add.main1)
    main2.dir <- instant.slope(1:length(add.main2), add.main2)

    if(all(main1.dir) > 0){

    }
    

    diff.mat <- add.mat - act.mat

}