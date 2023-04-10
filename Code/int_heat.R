#' Plot phenotypic effects for two markers as a heat map
#' 
#' This internal function is called by 
#' \code{\link{plot_effects}} to generate a 
#' heat map showing the effects of genotype on
#' phenotype. This function fits linear models
#' to the markers and traits. It then uses
#' these models to predict trait values at 
#' different genotype combinations in a 2D 
#' grid. It plots these predicted values in
#' a heat map.
#' 
#' @param phenoV A vector of trait values 
#' @param marker1_vals A vector of genotype values 
#' for marker1
#' @param marker2_vals A vector of genotype values
#' for marker2.
#' @param pheno_name A string indicating the name of
#' the trait being plotted.
#' @param marker1_label A string indicating the name
#' of marker1
#' @param marker2_label A string indicating the name
#' of marker2
#' @param nbins The number of bins for the markers over 
#' which to predict values of the trait. This can also 
#' be a vector specifying those bins.
#' @param grid.min A minimum for the grid over which to predict
#' phenotype values. If left null, the grid.min is the minimum 
#' value found in marker1_vals and marker2_vals
#' @param grid.max A maximum for the grid over which to predict
#' phenotype values. If left null, the grid.max is the maximum
#' value found in marker1_vals and marker2_vals
#'
#' @return None
#' 
#' @importFrom stats as.formula lm predict
#' @keywords internal

int_heat <- function(phenoV, marker1_vals, marker2_vals, pheno_name = NULL, 
                        marker1_label = NULL, marker2_label = NULL, nbins = 50,
						grid.min = NULL, grid.max = NULL){


	has_dash1 <- grep("-", marker1_label)
	if(length(has_dash1) > 0){
		marker1_label <- gsub("-", "_", marker1_label)
	}
	has_dash2 <- grep("-", marker2_label)
	if(length(has_dash2) > 0){
		marker2_label <- gsub("-", "_", marker2_label)
	}
	
    if(length(nbins) == 1){
		if(is.null(grid.min)){
			min_marker1 <- min(signif(marker1_vals, 2), na.rm = TRUE)
		}else{
			min_marker1 = grid.min
		}
		if(is.null(grid.max)){
			max_marker1 <- max(signif(marker1_vals, 2), na.rm = TRUE)
		}else{
			max_marker1 = grid.max
		}
        marker_grid1 <- segment_region(min_marker1, max_marker1, nbins, alignment = "ends")
    }else{
    	marker_grid1 <- nbins
    }
    
    if(length(nbins) == 1){
		if(is.null(grid.min)){
			min_marker2 <- min(signif(marker2_vals, 2), na.rm = TRUE)
		}else{
			min_marker2 = grid.min
		}
		if(is.null(grid.max)){
			max_marker2 <- max(signif(marker2_vals, 2), na.rm = TRUE)
		}else{
			max_marker2 = grid.max
		}
        marker_grid2 <- segment_region(min_marker2, max_marker2, nbins, alignment = "ends")
    }else{
       	marker_grid2 <- nbins
   	}
        
    marker1_bins <- bin_vector(signif(marker1_vals, 2), bins = marker_grid1)
    marker2_bins <- bin_vector(signif(marker2_vals, 2), bins = marker_grid2)
	#quartz();par(mfrow = c(2,2))
	#plot(marker1_vals, marker1_bins)
	#plot(marker2_vals, marker2_bins)
	#plot(jitter(marker1_bins), jitter(marker2_bins), xlim = c(0, 1), ylim = c(0,1))
	#pheno.col <- colors.from.values(phenoV, use.pheatmap.colors = TRUE)
	#plot.with.model(marker1_vals, marker2_vals, xlim = c(0, 1), ylim = c(0,1),col = pheno.col)

	#as best as we can, get the actual trait values for all bin pairs
	marker.bins <- segment.region(min(marker1_bins), max(marker1_bins), nbins, "ends")
	bin_pairs <- pair.matrix(marker.bins, self.pairs = TRUE, ordered = TRUE)
	pheno_mean <- pheno_count <- matrix(NA, nrow = nbins, ncol = nbins)
	rownames(pheno_mean) <- rownames(pheno_count) <- colnames(pheno_mean) <- colnames(pheno_count) <- marker.bins
	for(i in 1:nrow(bin_pairs)){
		row_locale <- which(marker.bins == bin_pairs[i,1])
		col_locale <- which(marker.bins == bin_pairs[i,2])
		pheno_locale <- intersect(which(marker1_bins == bin_pairs[i,1]), which(marker2_bins == bin_pairs[i,2]))
		if(length(pheno_locale) > 0){
			pheno_mean[row_locale, col_locale] <- mean(phenoV[pheno_locale])
			pheno_count[row_locale, col_locale] <- length(pheno_locale)
		}
	}
	#pheatmap(pheno_mean, cluster_rows = FALSE, cluster_cols = FALSE)
	#pheatmap(pheno_count, cluster_rows = FALSE, cluster_cols = FALSE)
	
	add.fmla <- paste(pheno_name, "~", marker1_label, "+", marker2_label)
	int.fmla <- paste(pheno_name, "~", marker1_label, "*", marker2_label)    
	fit_df <- data.frame(cbind(phenoV, marker1_vals, marker2_vals))
	colnames(fit_df) <- c(pheno_name, marker1_label, marker2_label)

	model.add <- lm(as.formula(add.fmla), data = fit_df)
	model.int  <- lm(as.formula(int.fmla), data = fit_df)
	#summary(model.add)
	#summary(model.int)
	
	predict_grid <- cbind(rep(marker_grid1, length(marker_grid2)), rep(marker_grid2, 
	each = length(marker_grid1)))
	colnames(predict_grid) <- c(marker1_label, marker2_label)
	predict_df <- data.frame(predict_grid)

	pred_data_add <- predict(model.add, newdata = predict_df)
	pred_data_int <- predict(model.int, newdata = predict_df)
	#plot(pred_data_add, pred_data_int)

	pred_mat.add <- matrix(pred_data_add, nrow = length(marker_grid2), ncol = length(marker_grid1), 
	byrow = FALSE)
	pred_mat.int <- matrix(pred_data_int, nrow = length(marker_grid2), ncol = length(marker_grid1), 
	byrow = FALSE)
	rownames(pred_mat.add) <- rownames(pred_mat.int) <- marker_grid1
	colnames(pred_mat.add) <- colnames(pred_mat.int) <- marker_grid2
	
	#pheatmap(pred_mat.add, cluster_rows = FALSE, cluster_cols = FALSE)
	#pheatmap(pred_mat.int, cluster_rows = FALSE, cluster_cols = FALSE)
	#pheatmap(pred_mat.add - pheno_mean, cluster_rows = FALSE, cluster_cols = FALSE)
	#pheatmap(pred_mat.int - pheno_mean, cluster_rows = FALSE, cluster_cols = FALSE)
	#par(mfrow = c(1,2));plot.with.model(as.vector(pheno_mean), as.vector(pred_mat.add))
	#plot.with.model(as.vector(pheno_mean), as.vector(pred_mat.int))
	add.error <- pred_mat.add - pheno_mean
	add.mse <- mean(add.error^2, na.rm = TRUE)
	int.error <- pred_mat.int - pheno_mean
	int.mse <- mean(int.error^2, na.rm = TRUE)
	
	return(list("Additive_Prediction" = pred_mat.add, 
	"Interactive_Prediction" = pred_mat.int, "Actual" = pheno_mean, 
	"Additive.MSE" = add.mse, "Interactive.MSE" = int.mse))
}
