import argparse
import numpy
import math

###################################################################
#                      LOAD OPTIONS AND FILES                     #
###################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--individual", help="Name of the sample to be analysed.")
parser.add_argument("-s", "--samples", type=argparse.FileType('r'), help="Text file with names of all samples in order.")
parser.add_argument("-k", "--groups", type=int, default=4, help="Number of ancestry groups. Default: 4.")
parser.add_argument("-p", "--frequencies", type=argparse.FileType('r'), help="Admixture's output file with allele frequencies in every position (rows) and ancestry groups (columns).")
parser.add_argument("-q", "--contributions", type=argparse.FileType('r'), help="Admixture's output file with proportion of ancestry contributed by each ancestry (columns) to each individual (columns).")
parser.add_argument("-g", "--genotypes", type=argparse.FileType('r'), help="Matrix of genotypes, coded as number of major alleles (0, 1, 2; -1 for missing data), with positions in columns and samples in rows. First column specifies the sample's index.")
parser.add_argument("-m", "--positions", type=argparse.FileType('r'), help="File with two columns: chromosome or contig name and position within, for each locus.")
parser.add_argument("-o", "--output", type=argparse.FileType('w'), default="z1.out", help="Output file.")

args = parser.parse_args()

Samples = args.samples.read().splitlines()
K = args.groups

# Note that the following loads large files in memory. The first column of
# Genotypes is just an index of the samples, from 0.
P = numpy.loadtxt(args.frequencies)
Q = numpy.loadtxt(args.contributions)
Genotypes = numpy.loadtxt(args.genotypes, dtype=int)
Chromosome = {}
Position = {}
Locus = 0
for line in args.positions:
   Chromosome[Locus], Position[Locus] = line.split()
   Locus += 1

assert (P.shape[0] == Genotypes.shape[1] - 1), "Different number of loci in the P and the genotypes matrices."
assert (P.shape[1] == K), "The number of ancestry groups does not match with the number of columns in the P matrix."
assert (Q.shape[0] == len(Samples)), "The number of samples does not match the number of rows in the Q matrix."
assert (Q.shape[1] == K), "The number of columns in the Q matrix is different from the number of ancestry groups."
assert (Genotypes.shape[0] == len(Samples)), "The number of rows in the genotypes matrix does not match the number of samples."
assert (args.individual in Samples), "The individual to be analysed is not in the list."
assert (len(Chromosome) == len(Position) == P.shape[0]), "The numer of loci in the P matrix does not match the number of rows in the chromosomic positions file."

NumLoci = P.shape[0]
TargetIndividualIndex = Samples.index(args.individual)

##################################################################
#                            Functions                           #
##################################################################

def GetGenotype(AltAlleleCount, index):
   if AltAlleleCount[index] == -1:
      return -1
   if AltAlleleCount[index] == 1:
      return 1
   SamplesWithData = (AltAlleleCount >= 0)
   if sum(SamplesWithData) > 0:
      AltAlleleFreq = sum(AltAlleleCount[SamplesWithData]) / (2.0 * sum(SamplesWithData))
      if AltAlleleFreq > 0.5:
         return AltAlleleCount[index]
      if AltAlleleFreq < 0.5:
         return 2 - AltAlleleCount[index] # Here, it should not be -1.
      if AltAlleleFreq == 0.5:
         return 3 # This means either 0 or 2.

def CalculatePosterior(genotype, allele, freqArray, priorArray, ancestry):
   '''Calculates the posterior probability of an ancestry for a particular allele of a genotype'''
   if genotype == -1:
      return priorArray[ancestry]   # without data, the posterior is equal to the prior
   if genotype == 0:
      return (1.0 - freqArray[ancestry]) * priorArray[ancestry] / sum(priorArray * (1.0 - freqArray))
   if genotype == 1:
      if allele == 0:
         return (1.0 - freqArray[ancestry]) * priorArray[ancestry] / sum(priorArray * (1.0 - freqArray))
      if allele == 1:
         return freqArray[ancestry] * priorArray[ancestry] / sum(priorArray * freqArray)
   if genotype == 2:
      return freqArray[ancestry] * priorArray[ancestry] / sum(priorArray * freqArray)
   if genotype == 3:
      posterior0 = (1.0 - freqArray[ancestry]) * priorArray[ancestry] / sum(priorArray * (1.0 - freqArray))
      posterior2 = freqArray[ancestry] * priorArray[ancestry] / sum(priorArray * freqArray)
      return (posterior0, posterior2)


###################################################################
#                           Main                                  #
###################################################################

for locus in range(NumLoci):
   genotype = GetGenotype(Genotypes[:, locus + 1], TargetIndividualIndex)
   PostAncProb = {0: {}, 1: {}}     # Posterior Ancestry Probabilities: allele -> ancestry -> probability
   for ancestry in range(K):
      PostAncProb[0][ancestry] = CalculatePosterior(genotype, 0, P[locus, :], Q[TargetIndividualIndex, :], ancestry)
      PostAncProb[1][ancestry] = CalculatePosterior(genotype, 1, P[locus, :], Q[TargetIndividualIndex, :], ancestry)
   if genotype == 3:   # In this case, GetGenotype returned a tuple, and I need to choose one of the two values.
      MLAncestry = Q[TargetIndividualIndex, :].argmax()
      ExpCopies0 = [PostAncProb[0][x][0] + PostAncProb[1][x][0] for x in range(K)]
      MostProb0 = ExpCopies0.index(max(ExpCopies0))
      ExpCopies2 = [PostAncProb[0][x][1] + PostAncProb[1][x][1] for x in range(K)]
      MostProb2 = ExpCopies2.index(max(ExpCopies2))
      if MLAncestry == MostProb0 and MostProb0 != MostProb2:
         for ancestry in range(K):
            PostAncProb[0][ancestry] = PostAncProb[0][ancestry][0]
            PostAncProb[1][ancestry] = PostAncProb[1][ancestry][0]
      elif MLAncestry == MostProb2 and MostProb0 != MostProb2:
         for ancestry in range(K):
            PostAncProb[0][ancestry] = PostAncProb[0][ancestry][1]
            PostAncProb[1][ancestry] = PostAncProb[1][ancestry][1]
      else:
         continue
   AtLeastOne = [ PostAncProb[0][x] * (1.0 - PostAncProb[1][x]) + 
                  (1.0 - PostAncProb[0][x]) * PostAncProb[1][x] + 
                  PostAncProb[0][x] * PostAncProb[1][x] for x in range(K) ]
   try:
      LogOdds = [ math.log(AtLeastOne[x]) - math.log(1.0 - AtLeastOne[x]) for x in range(K) ]
   except ValueError:
      LogOdds = []
      for i in range(K):
         if AtLeastOne[i] >= 1.0:
            AtLeastOne[i] = 0.9999999999999999
         LogOdds.append(math.log(AtLeastOne[i]) - math.log(1.0 - AtLeastOne[i]))

   args.output.write('{}\t{}\t{:f}\t{:f}\t{:f}\t{:f}\n'.format(Chromosome[locus], Position[locus], LogOdds[0], LogOdds[1], LogOdds[2], LogOdds[3]))

