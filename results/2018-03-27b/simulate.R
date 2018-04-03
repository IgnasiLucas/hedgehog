ObserveGenotypes <- function(x){
   # x is a 2 * n matrix of numbers of reads of each of two alleles in
   # n individuals.
   observed <- c(0, 0, 0)
   for (i in 1:dim(x)[2]) {
      if ((x[1,i] >  0) & (x[2,i] == 0)) observed[1] <- observed[1] + 1
      if ((x[1,i] >  0) & (x[2,i] >  0)) observed[2] <- observed[2] + 1
      if ((x[1,i] == 0) & (x[2,i] >  0)) observed[3] <- observed[3] + 1
   }
   return(observed)
}

EstimateWell <- function(x, MaxDist=0.0001, MaxRep=100){
   # x is a 2 * n matrix of numbers of reads of each of two alleles in n
   # individuals. This is a cumbersome way to estimate allele frequency.

   AA <- (x[1,] >  0) & (x[2,] == 0)
   Aa <- (x[1,] >  0) & (x[2,] >  0)
   aa <- (x[1,] == 0) & (x[2,] >  0)
   ObsGenoNum <- c(sum(AA), sum(Aa), sum(aa))
   ExpGenoNum <- numeric(length = 3)
   GenoFreq <- ObsGenoNum / dim(x)[2]
   distance <- 1
   NumRep <- 0
   while ((distance > MaxDist) && (NumRep < MaxRep)) {
      # I came up with the formulaes below. The idea is that the expected number of
      # homozygous individuals is less than observed, because a portion of those
      # apparently homozygous must be undersampled heterozygous. The expected number
      # of undersampled heterozygous genotypes contributing to the homozygous counts
      # is a summation across the posterior probability of observed homozygous being
      # in fact undersampled heterozygous. Because we need estimates of the real
      # genotype frequencies, I run an iteration similar to an expectation-maximization
      # algorithm: I first estimate the real number of each genotype, then I improve
      # with the genotype frequencies used in the previous step, and so on.

      ExpGenoNum[1] <- ObsGenoNum[1] - GenoFreq[2] * sum(1 / (GenoFreq[2] + GenoFreq[1] * 2^x[1,AA]))
      ExpGenoNum[2] <- ObsGenoNum[2] + GenoFreq[2] * sum(1 / (GenoFreq[2] + GenoFreq[1] * 2^x[1,AA])) + 
                                       GenoFreq[2] * sum(1 / (GenoFreq[2] + GenoFreq[3] * 2^x[2,aa]))
      ExpGenoNum[3] <- ObsGenoNum[3] - GenoFreq[2] * sum(1 / (GenoFreq[2] + GenoFreq[3] * 2^x[2,aa]))
      NewGenoFreq <- ExpGenoNum / dim(x)[2]
      distance <- sqrt(sum((NewGenoFreq - GenoFreq)^2))
      GenoFreq <- NewGenoFreq
      NumRep <- NumRep + 1
   }
   return(GenoFreq[1] + GenoFreq[2]/2)
}

NumSim <- 100
Sizes <- 5:12
Depths <- 1:10
Frequencies <- seq(from = 0.05, to = 0.50, by = 0.05)
TotalSize <- NumSim * length(Sizes) * length(Depths) * length(Frequencies)
results <- data.frame(NumInd = numeric(length = TotalSize),
                      NumRead = numeric(length = TotalSize),
                      MAF = numeric(length = TotalSize),
                      Estim1 = numeric(length = TotalSize),
                      Estim2 = numeric(length = TotalSize))
index <- 0
for (SampleSize in Sizes) {
   for (Depth in Depths) {
      for (Freq in Frequencies) {
         genotypes <- rmultinom(NumSim, SampleSize, prob=c(Freq^2, 2*Freq*(1-Freq), (1-Freq)^2))
         for (i in 1:NumSim) {
            index <- index + 1
            reads <-cbind(rmultinom(genotypes[1,i], Depth, prob=c(1,0)),
                          rmultinom(genotypes[2,i], Depth, prob=c(0.5, 0.5)),
                          rmultinom(genotypes[3,i], Depth, prob=c(0,1)))
            observed <- ObserveGenotypes(reads)
            results$Estim1[index] <- (2 * observed[1] + observed[2]) / (2 * SampleSize)
            results$Estim2[index] <- EstimateWell(reads, MaxDist=0.001, MaxRep=50)
            results$NumInd[index] <- SampleSize
            results$NumRead[index] <- Depth
            results$MAF[index]    <- Freq
         }
      }
   }
}

rm(NumSim, TotalSize, index)

attach(results)
variances <- data.frame(NumInd = numeric(length = length(Sizes) * length(Depths) * length(Frequencies)),
                       NumRead = numeric(length = length(Sizes) * length(Depths) * length(Frequencies)),
                           MAF = numeric(length = length(Sizes) * length(Depths) * length(Frequencies)),
                        ExpVar = numeric(length = length(Sizes) * length(Depths) * length(Frequencies)),
                        MaxVar = numeric(length = length(Sizes) * length(Depths) * length(Frequencies)),
                       VarEst1 = numeric(length = length(Sizes) * length(Depths) * length(Frequencies)),
                       VarEst2 = numeric(length = length(Sizes) * length(Depths) * length(Frequencies)))
m <- 1
for (i in Sizes) {
   for (j in Depths) {
      for (k in Frequencies) {
         f <- (NumInd == i) & (NumRead == j) & (MAF == k)
         variances$NumInd[m]  <- i
         variances$NumRead[m] <- j
         variances$MAF[m]     <- k
         variances$MaxVar[m]  <- k * (1 - k) / i         # Expected variance if we sampled only one chromosome per individual.
         variances$ExpVar[m]  <- k * (1 - k) / (2 * i)   # Sampling 2 chromosomes per individual.
         variances$VarEst1[m] <- var(Estim1[f])
         variances$VarEst2[m] <- var(Estim2[f])
         m <- m + 1
      }
   }
}

detach(results)
attach(variances)
f <- (MAF > 0.3) & (NumRead > 1) & (NumRead < 6)
png(filename = "estimation/variances.png")
plot(VarEst1[f], VarEst2[f], xlab = "Variance of conventional allele frequency estimate",
                             ylab = "Variance of new allele frequency estimate")
abline(0, 1, col='red')
dev.off()
detach(variances)
