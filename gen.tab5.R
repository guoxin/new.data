load("new.Rdata")
TheLambda <- exp(-13)
load("readpan20120114 md5.cc4a88a33be3296333121e9c4d524a8f")
source("internals.R")
beta.allele <- 0.06
beta.peptide <- 0.11387

#Nielsen's codes:http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan

total.peptides <- unique(c(newData$peptide.set, readpan20120114$peptide.set))
names(total.peptides) <- total.peptides
total.alleles <- unique(c(newData$allele.set, readpan20120114$allele.set))
temp.al <- c(newData$allele.set, readpan20120114$allele.set)
temp.nm <- c()
for(i in 1:length(total.alleles))temp.nm[i] <- names(temp.al)[which(total.alleles[i] == temp.al)[1]]
names(total.alleles) <- temp.nm

peptide.K3Hat <- stringToK3(string.set = total.peptides, K1 = readpan20120114$K1, beta = beta.peptide)
allele.K3Hat <- stringToK3(string.set = total.alleles, K1 = readpan20120114$K1, beta = beta.allele)

K3Training <- peptide.K3Hat[readpan20120114$peptide.list, readpan20120114$peptide.list] *
    allele.K3Hat[readpan20120114$allele.list, readpan20120114$allele.list]; gc()

diag(K3Training) <- 1 + TheLambda * 33931
coeff <- solve(K3Training, readpan20120114$rValue)
rm(K3Training); gc()

K3Testing <- peptide.K3Hat[newData$peptide.list, readpan20120114$peptide.list] *
    allele.K3Hat[names(newData$allele.list), readpan20120114$allele.list]; gc()

predicted <- K3Testing %*% coeff

rm(K3Testing); gc()
report <- list(coeff=coeff, predicted = predicted, peptide = newData$peptide.list, allele = names(newData$allele.list))
save(report, file = "report.Rdata")
