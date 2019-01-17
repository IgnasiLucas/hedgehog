import argparse

# vcfpy is developed by Manuel Holtgrewe (Berlin Institute of Health), and released under
# the MIT license. Similar to pyvcf, but newer, and less tested.
import vcfpy

def parse_popmap(popmap):
   '''Creates a dictionary with sample names as keys and population names as values.'''
   table = {}
   for line in popmap:
      if not line.startswith('#'):
         table[line.split()[0]] = line.split()[1]
   return table

def make_poplist(popmap):
   '''Creates a dictionary with population name as key and a list of sample names as value.'''
   table = {}
   for line in popmap:
      if not line.startswith('#'):
         (sample, population) = line.split()
         try:
            table[population].append(sample)
         except KeyError:
            table[population] = [sample]
   return table

def main():
   parser = argparse.ArgumentParser(description = 'Takes a vcf file and a population file and produces the input files required for Bayesian Genomic Clines estimation (bgc), namely the two parental populations data files and the admixed population(s) data file. The genetic map is not implemented yet.')
   parser.add_argument('-0', '--pop0', type=str, default='roumanicus', help='Name of parental population 0, as in the populations mapping file. Default: roumanicus.')
   parser.add_argument('-1', '--pop1', type=str, default='europaeus', help='Name of parental population 1, as in the population mapping file. Default: europaeus.')
   parser.add_argument('-a', '--admixed', type=str, action='append', default='admixed', help='Name of an admixed population, as in the population mapping file. Repeat to add more than one. Default: admixed.')
   parser.add_argument('-p', '--prefix', type=str, default="", help='Prefix to be added to the predefined output names, which are based on the names of the populations. For example, to save output in a particular folder.')
   parser.add_argument('-u', '--uncertain', action='store_true', help='Genotypes are uncertain, and an estimate of sequencing errors must be included in output. Requires specification of the error through -e. Default: false.')
   parser.add_argument('-e', '--errors', type=float, help='Common sequencing error value for all sites. If zero, site-specific error values will be estimated from vcf statistics, if available, provided -u. If -u is not used, -e is ignored. Default: 0.0.')
   parser.add_argument('-v', '--vcf', type=argparse.FileType('r'), default='-', help='Input vcf file. Default: standard input.')
   parser.add_argument('popmap', type=argparse.FileType('r'), help='Input populations mapping file, tab or space separated: sample name and population.' )
   args = parser.parse_args()
#   PopulationDict = parse_popmap(args.popmap)
   PopulationList = make_poplist(args.popmap)
   assert args.pop0 in PopulationDict.values()
   assert args.pop1 in PopulationDict.values()
   assert args.admixed in PopulationDict.values()
   pop0 = open(args.prefix + 'pop0.txt', 'w')
   pop1 = open(args.prefix + 'pop1.txt', 'w')
   admixed = open(args.prefix + 'admixed.txt', 'w')
   reader = vcfpy.Reader.from_stream(vcf)
   if args.uncertain:
      if args.errors == 0.0:

      elif args.errors > 0.0:

      else:
   else:
      for record in reader:
         for sample in PopulationList[args.pop0]:
            genotype = record.call_for_sample[sample].data.get('GT')
