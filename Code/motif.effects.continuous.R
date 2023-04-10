#This function looks at the main effects in different motifs

motif.effects.continuous <- function(data.obj, geno.obj, 
	scan.what = c("normalized", "eigentraits", "raw"), 
	interaction = c("enhancing", "suppressing"), 
	main = c("coherent", "incoherent"),
	source.sign = c("pos", "neg", "any"), p_or_q = 0.05,
	n.geno.bins = 50, separate.windows = FALSE){
	

	main.dir <- function(mean.effects){
		if(mean.effects[1] < tail(mean.effects, 1)){
			return("pos")
		}
		if(mean.effects[1] > tail(mean.effects, 1)){
			return("neg")
		}
		if(mean.effects[1] == tail(mean.effects, 1)){
			return("none")
		}
	}

	# data.obj <- get.network(data.obj, p.or.q = p.or.q, collapse.linked.markers = TRUE)
	motif.obj <- find.motifs(data.obj, p_or_q = p_or_q)
	
	#collect the unique motifs for all phenotypes
	motifs <- Reduce("rbind", motif.obj[[1]])
	motif.dir <- Reduce("rbind", motif.obj[[2]])
	motif.table <- cbind(motifs, motif.dir)

	u_pheno <- unique(motifs[,3])
	motif.pheno.effects <- vector(mode = "list", length = length(u_pheno))
	names(motif.pheno.effects) <- u_pheno
	
	for(ph in 1:length(u_pheno)){
		pheno.locale <- which(motifs[,3] == u_pheno[ph])
		pheno.motif.table <- motif.table[pheno.locale,]
		sign.locale <- 1:nrow(pheno.motif.table)

		if(nrow(pheno.motif.table) > 0){		
			
			motif.effects <- lapply(1:nrow(pheno.motif.table), 
				function(m) get.pheno.effects.continuous(data.obj, geno.obj, 
				marker1.name = pheno.motif.table[m,1], marker2.name = pheno.motif.table[m,2], 
				pheno.name = u_pheno[ph], covar = data.obj$p_covar, 
				scan.what = scan.what, n.geno.bins = n.geno.bins))

			add.effects <- lapply(motif.effects, function(x) x$Additive_Prediction - mean(x$Additive_Prediction))
			int.effects <- lapply(motif.effects, function(x) x$Interactive_Prediction - mean(x$Interactive_Prediction))
			actual.effects <- lapply(motif.effects, function(x) x$Actual - mean(x$Actual, na.rm = TRUE))
			geno.text <- c("0", rep("", (ncol(actual.effects[[1]])-2)), "1")

			source.main <- lapply(add.effects, function(x) colMeans(x, na.rm = TRUE))
			target.main <- lapply(add.effects, function(x) rowMeans(x, na.rm = TRUE))
			source.dir <- sapply(source.main, main.dir)
			target.dir <- sapply(target.main, main.dir)
			dir.table <- cbind(source.dir, target.dir)

			#detertmine which motifs to calculate effects for
			if(main == "coherent"){
				main.locale <- which(apply(dir.table, 1, function(x) x[1] == x[2]))
			}else{
				main.locale <- which(apply(dir.table, 1, function(x) x[1] != x[2]))
			}

			if(interaction == "enhancing"){
				int.locale <- which(pheno.motif.table[,"source-target"] == 1)
			}else{
				int.locale <- which(pheno.motif.table[,"source-target"] == -1)
			}

			if(source.sign != "any"){
				sign.locale <- which(dir.table[,1] == source.sign)
			}
			
			final.locale <- Reduce("intersect", list(main.locale, int.locale, sign.locale))
			motif.pheno.effects[[ph]] <- pheno.motif.table[final.locale,,drop=FALSE]

			if(length(final.locale) > 0){
				sub.add.effects <- add.effects[final.locale]
				sub.int.effects <- int.effects[final.locale]
				sub.actual.effects <- actual.effects[final.locale]
				sub.dir <- dir.table[final.locale,]

				#pheatmap(sub.add.effects[[1]], cluster_rows = FALSE, cluster_cols = FALSE)
				#pheatmap(sub.int.effects[[1]], cluster_rows = FALSE, cluster_cols = FALSE)
				#pheatmap(sub.actual.effects[[1]], cluster_rows = FALSE, cluster_cols = FALSE)

				#if(separate.windows){quartz(width = 9, height = 6)}
				#par(mfrow = c(2,3))
				#boxplot(sub.add.effects[[m]], main = "Main Effect 1 Additive Prediction")
				#boxplot(sub.int.effects[[m]], main = "Main Effect 1 Interactive Prediction")
				#boxplot(sub.actual.effects[[m]], main = "Main Effect 1 Actual Values")
				#boxplot(t(sub.add.effects[[m]]), main = "Main Effect 1 Additive Prediction")
				#boxplot(t(sub.int.effects[[m]]), main = "Main Effect 1 Interactive Prediction")
				#boxplot(t(sub.actual.effects[[m]]), main = "Main Effect 1 Actual Values")
				
				add.mse <- sapply(motif.effects[final.locale], function(x) x$Additive.MSE)
				int.mse <- sapply(motif.effects[final.locale], function(x) x$Interactive.MSE)
				#boxplot(list(add.mse, int.mse), names = c("Additive", "Interactive"))
				#plot(add.mse, int.mse);abline(0,1)
								
				addV <- unlist(sub.add.effects)	
				intV <- unlist(sub.int.effects)	
				actV <- unlist(sub.actual.effects)	

				if(separate.windows){quartz(width = 8, height = 4)}
				par(mfrow = c(1,2))
				plot.with.model(actV, addV, main = "Additive Model", xlab = "Trait Value",
				ylab = "Predicted Trait Value")
				plot.with.model(actV, intV, main = "Interactive Model", xlab = "Trait Value",
				ylab = "Predicted Trait Value")
				mtext(paste(u_pheno[ph], main, interaction, "with", source.sign, 
				"main effects"), side = 3, outer = TRUE, line = -1)

				cross.section.add <- lapply(1:length(sub.add.effects[[1]]), 
				function(x) sapply(sub.add.effects, function(y) y[x]))
				cross.section.int <- lapply(1:length(sub.int.effects[[1]]), 
				function(x) sapply(sub.int.effects, function(y) y[x]))
				cross.section.actual <- lapply(1:length(sub.actual.effects[[1]]), 
				function(x) sapply(sub.actual.effects, function(y) y[x]))

				#find the deviation of the actual value from the predicted 
				#value for each combination
				add.dev <- lapply(1:length(sub.add.effects), 
				function(x) sub.add.effects[[x]] - sub.actual.effects[[x]])
				int.dev <- lapply(1:length(sub.int.effects), 
				function(x) sub.int.effects[[x]] - sub.actual.effects[[x]])

				#depending on the motif, we expect the 1,1 trait to be
				#either more extreme or less extreme than the additive
				#model predicts. Coherent-enhancing, and Incoherent-suppressing
				#interactions are likely to have more extreme effects than
				#predicted, while coherent-suppressing, and incoherent-enhancing
				#are likely to have less extreme effects

				add.dev <- matrix(sapply(1:length(cross.section.add), 
					function(x) mean(abs(cross.section.add[[x]] - cross.section.actual[[x]]), 
					na.rm = TRUE)), nrow = nrow(sub.add.effects[[1]]), 
					ncol = ncol(sub.add.effects[[1]]), byrow = FALSE)
				int.dev <- matrix(sapply(1:length(cross.section.int), 
					function(x) mean(abs(cross.section.int[[x]] - cross.section.actual[[x]]), 
					na.rm = TRUE)), nrow = nrow(sub.int.effects[[1]]), 
					ncol = ncol(sub.int.effects[[1]]), byrow = FALSE)
				dimnames(add.dev) <- dimnames(int.dev) <- list(geno.text, geno.text)

				#par(mfrow = c(1,2))
				#imageWithText(rotate.mat(rotate.mat(rotate.mat(add.dev))), 
				#use.pheatmap.colors = TRUE, show.text = FALSE)
				#imageWithText(rotate.mat(rotate.mat(rotate.mat(int.dev))), 
				#use.pheatmap.colors = TRUE, show.text = FALSE)

				if(separate.windows){quartz(width = 8, height = 4)}
				par(mfrow = c(1,2))
				plot(as.vector(add.dev), as.vector(int.dev),
				xlab = "Deviation from Additive Model", pch = 16,
				ylab = "Deviation from Interactive Model")
				abline(0,1)

				add.better <- which(abs(add.dev) < abs(int.dev)) 
				int.better <- which(abs(add.dev) > abs(int.dev)) 
				err.mat <- add.dev - int.dev
				#class.mat <- add.dev
				#class.mat[add.better] <- 1
				#class.mat[int.better] <- 2
				#dimnames(err.mat) <- dimnames(add.dev)
				rot.err <- rotate.mat(rotate.mat(rotate.mat(err.mat)))
				#rot.class <- rotate.mat(rotate.mat(rotate.mat(class.mat)))
				#imageWithText(signif(rot.err, 2), class.mat = rot.class, 
				#col.scale = c("blue", "brown"), col.text.rotation = 0, light.dark = "f",
				#grad.dir = "ends")

				imageWithText(signif(rot.err, 2), split.at.vals = TRUE, 
				col.scale = c("blue", "brown"), col.text.rotation = 0, light.dark = "f",
				grad.dir = "ends")

				mtext(paste(u_pheno[ph], main, interaction, "with", source.sign, 
				"main effects"), side = 3, outer = TRUE, line = -2)
		}else{
			plot.text("No interactions for these parameters.")	
			} #end case for if there are phenotypes matching this description
		} #end case for if there are motifs for this phenotype
	} #end looping through phenotypes
	invisible(motif.pheno.effects)
}