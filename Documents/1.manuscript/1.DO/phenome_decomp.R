test <- bigcor(scaled.expr[[1]])


low.cor <- vector(mode = "list", length = nrow(test))
names(low.cor) <- colnames(scaled.expr[[1]])
for(i in 1:10){
    report.progress(i, ncol(scaled.expr[[1]]))
    low.cor[[i]]  <- which(abs(test[i,]) <= 0.14)
}

library(glmnet)
alpha = 1 #for lasso, alpha = 0 for ridge regression

#bin individuals into 10 groups
ind.bins <- bin.sequence(sample(1:nrow(scaled.expr[[1]])), 10)

for(i in 1:ncol(scaled.expr[[1]])){
    y <- scaled.expr[[1]][,i,drop=FALSE] #one transcript
    x <- scaled.expr[[1]][,low.cor[[i]]] #transcripts that are uncorrelated with this transcript

    #perform 10-fold cv 
    pred.mat <- matrix(NA, nrow(scaled.expr[[1]]), ncol = length(ind.bins))    
    for(cv in 1:length(ind.bins)){
        print(cv)
        leave.out <- rownames(scaled.expr[[1]])[ind.bins[[cv]]]
        leave.in <- setdiff(rownames(scaled.expr[[1]]), leave.out)
        cv_model <- cv.glmnet(x[leave.in,], y[leave.in,,drop=FALSE], alpha = alpha)
        best_lambda <- cv_model$lambda.min
        #plot(cv_model)
        best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
        pred.mat[,cv] <- predict(best_model, newx = x)
    }
    pred.trans <- rowMeans(pred.mat)
    names(pred.trans) <- rownames(scaled.expr[[1]])
    plot(pred.trans, scaled.expr[[1]][,i])    
    
    orig.scan <- scan1(genoprobs, scaled.expr[[1]][,i])
    pred.scan <- scan1(genoprobs, pred.trans)
    maxy <- max(c(pred.scan, orig.scan))
    plot(orig.scan, map = map, ylim = c(0, maxy*1.05))
    plot(pred.scan, map = map, add = TRUE, col = "lightblue")

    plot(orig.scan[,1], pred.scan[,1])

    orig.herit <- as.numeric(est_herit(pheno = scaled.expr[[1]][,i,drop=FALSE], kinship = Kt[[1]]))
    pred.herit <- as.numeric(est_herit(pheno = pred.trans, kinship = Kt[[1]]))
}


x = matrix(rnorm(100 * 20), 100, 20)
y = rnorm(100)
fit1 = glmnet(x, y)
print(fit1)
coef(fit1, s = 0.01)  # extract coefficients at a single value of lambda
predict(fit1, newx = x[1:10, ], s = c(0.01, 0.005))  # make predictions
