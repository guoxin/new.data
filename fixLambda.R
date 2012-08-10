load("readpan20120114 md5.cc4a88a33be3296333121e9c4d524a8f")
source("internals.R")
lambda.set <- exp(seq(-18,-8))
LooNewReport <- function(peptideStrings, alleleStrings, 
                         peptideList, alleleList,
                         K1, rValues,
                         lambda.set){
    #Warning: assume "diag" works well below.
    thisReport <- list()
    beta.allele <- 0.06
    beta.peptide <- 0.11387

    peptide.K3Hat <- stringToK3(string.set = peptideStrings, K1 = K1, beta = beta.peptide)
    allele.K3Hat <- stringToK3(string.set = alleleStrings, K1 = K1, beta = beta.allele)
    K3Hat <- peptide.K3Hat[peptideList, peptideList] * allele.K3Hat[alleleList, alleleList]; gc()
    gpgx.time.stamp("K3Hat obtained")

    rmse.trn <- rep(0, length(lambda.set))
    names(rmse.trn) <- log(lambda.set)

    for(lambda.id in 1:length(lambda.set)){
        diag(K3Hat) <- 1 + 33930 * lambda.set[lambda.id]
        G.inv <- solve(K3Hat); gc();  gpgx.time.stamp("G.inv obtained")
        coeff <- G.inv %*% rValues; G.inv <- diag(G.inv); gc(); gpgx.time.stamp("coeff obtained")
        rmse.trn[lambda.id] <- sqrt(mean((coeff / (G.inv))^2)); gpgx.time.stamp("rmse.trn found")
    }
    return(rmse.trn)
}
reportsPan <- LooNewReport(
    peptideStrings = readpan20120114$peptide.set, 
    alleleStrings = readpan20120114$allele.set,
    peptideList = readpan20120114$peptide.list,
    alleleList = readpan20120114$allele.list,
    K1 = readpan20120114$K1,
    rValues = readpan20120114$rValue,
    lambda.set = lambda.set
)
save(reportsPan, file = "reportsPan")
