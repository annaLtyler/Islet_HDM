#glucose in DO

test <- read.csv("~/Documents/Data/DO/Phenotype_Data/svenson850_phenotypes_v2.csv")


compare_dist <- function(colname = "glucose1", comp.col = "diet", thresh = 250, 
	cols = c("black", "#3182bd")){
	
	max.glu1 <- max(test[, colname], na.rm = TRUE)

	groups <- unique(test[,comp.col])

	g1.idx <- which(test[,comp.col] == groups[1])
	g2.idx <- which(test[,comp.col] == groups[2])

	no.g1.na <- which(!is.na(test[g1.idx, colname]))
	dens.g1 <- density(test[g1.idx[no.g1.na], colname])

	no.g2.na <- which(!is.na(test[g2.idx, colname]))
	dens.g2 <- density(test[g2.idx[no.g2.na], colname])

	maxy <- max(c(dens.g1$y, dens.g2$y))
	maxx <- max(c(dens.g1$x, dens.g2$x))
	plot(density(test[g1.idx[no.g1.na], colname]), xlim = c(0, maxx), ylim = c(0, maxy), 
		xlab = colname, lwd = 3, main = paste(colname, " levels in the DO"), col = cols[1])
	points(density(test[g2.idx[no.g2.na], colname]), col = cols[2], type = "l", lwd = 3)
	abline(v = thresh, lwd = 2, lty = 2)
	legend("topright", legend = c(groups), col = cols, lty = 1, lwd = 3)
	
	above.thresh.g1 <- length(which(test[g1.idx[no.g1.na],colname] >= thresh))
	below.thresh.g1 <- length(which(test[g1.idx[no.g1.na],colname] < thresh))
	prop.g1 <- above.thresh.g1/length(no.g1.na)

	above.thresh.g2 <- length(which(test[g2.idx[no.g2.na],colname] >= thresh))
	below.thresh.g2 <- length(which(test[g2.idx[no.g2.na],colname] < thresh))
	prop.g2 <- above.thresh.g2/length(no.g2.na)
	
	result <- rbind(c(below.thresh.g1, above.thresh.g1, prop.g1), 
		c(below.thresh.g2, above.thresh.g2, prop.g2))	
	colnames(result) <- c("Num Below", "Num Above", "Proportion Above")	
	rownames(result) <- groups
	return(result)
}

par(mfrow = c(1,2))
compare_dist(colname = "glucose1", comp.col = "diet", thresh = 250)
compare_dist("glucose2", comp.col = "diet", thresh = 250)

par(mfrow = c(1,2))
compare_dist(colname = "glucose1", comp.col = "sex", thresh = 250)
compare_dist("glucose2", comp.col = "sex", thresh = 250)

boxplot(test[,"glucose1"]~test[,"diet"]*test[,"sex"], xlab = "Group", ylab = "Glucose (mg/dL)")
abline(h = 250)

boxplot(test[,"glucose2"]~test[,"diet"]*test[,"sex"], xlab = "Group", ylab = "Glucose (mg/dL)")
abline(h = 250)
