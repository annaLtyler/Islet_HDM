# estimated SE of estimated heritability
#
# written by Karl Broman: https://github.com/rqtl/qtl2/issues/193
#
# Visscher & Goddard (2015), https://doi.org/10.1534/genetics.114.171017
#
# lambda_i = eigenvalues of K
#
# a = sum[ (lam_i-1)^2 / (1 + h^2(lam_i - 1))^2 ]
#
# b = sum[ (lam_i - 1) / (1 + h^2(lam_i - 1)) ]
#
# var(h_hat^2) = 2 / (a - b^2/N)
#
# N = sample size

est_herit_error <-
    function(pheno, kinship)
{

    h2 <- as.numeric(est_herit(pheno, kinship))

    d <- decomp_kinship(kinship)

    lam <- d$values
    a <- sapply(h2, function(hsq) sum( (lam-1)^2 / (1 + hsq*(lam-1))^2 ) )
    b <- sapply(h2, function(hsq) sum( (lam-1) / (1 + hsq*(lam-1)) ) )
    n <- length(lam)

    setNames( sqrt(2 / (a - b^2/n)), colnames(pheno) )
}