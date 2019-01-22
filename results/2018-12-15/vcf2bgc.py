import argparse

# vcfpy is developed by Manuel Holtgrewe (Berlin Institute of Health), and released under
# the MIT license. Similar to pyvcf, but newer, and less tested.
import vcfpy

def MakePopDict(popmap):
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
   parser.add_argument('-a', '--admixed', action='append', help='Name of an admixed population, as in the population mapping file. Repeat to add more than one. Default: admixed.')
   parser.add_argument('-p', '--prefix', type=str, default="", help='Prefix to be added to the predefined output names, which are based on the names of the populations. For example, to save output in a particular folder.')
   parser.add_argument('-u', '--uncertain', action='store_true', help='Genotypes are uncertain, and an estimate of sequencing errors must be included in output. Requires specification of the error through -e. Default: false.')
   parser.add_argument('-e', '--errors', type=float, help='Common sequencing error value for all sites. If zero, site-specific error values will be estimated from vcf statistics, if available, provided -u. If -u is not used, -e is ignored. Default: 0.0.')
   parser.add_argument('-v', '--vcf', type=argparse.FileType('r'), default='-', help='Input vcf file. Default: standard input.')
   parser.add_argument('popmap', type=argparse.FileType('r'), help='Input populations mapping file, tab or space separated: sample name and population.' )
   args = parser.parse_args()
   SamplesInPop = MakePopDict(args.popmap)
   assert args.pop0 in SamplesInPop.keys()
   assert args.pop1 in SamplesInPop.keys()
   for pop in args.admixed:
      assert pop in SamplesInPop.keys()
   pop0 = open(args.prefix + 'pop0.txt', 'w')
   pop1 = open(args.prefix + 'pop1.txt', 'w')
   admixed = open(args.prefix + 'admixed.txt', 'w')
   reader = vcfpy.Reader.from_stream(args.vcf)
   # I want to remove from the populations' lists of samples the samples that are not present in the vcf.
   for pop in SamplesInPop.keys():
      tmp = [sample for sample in SamplesInPop[pop] if sample in reader.header.samples.names]
      SamplesInPop[pop] = tmp
   if args.uncertain:
      # Site-specific errors
      if args.errors == 0.0:
         pass
      # Common error rate in all sites
      else:
         for record in reader:
            if './.' in [cal.data.get('GT') for call in record.calls]:
               continue
            pop0.write("l_{}:{}\n".format(record.CHROM, record.POS))
            pop1.write("l_{}:{}\n".format(record.CHROM, record.POS))
            for sample in SamplesInPop[args.pop0]:
               pop0.write("{} {}\n".format(record.call_for_sample[sample].data.get('RO'), record.call_for_sample[sample].data.get('AO')[0]))
            for sample in SamplesInPop[args.pop1]:
               pop1.write("{} {}\n".format(record.call_for_sample[sample].data.get('R0'), record.call_for_sample[sample].data.get('AO')[0]))
            for pop in args.admixed:
               admixed.write("p_{}\n".format(pop))
               for sample in SamplesInPop[pop]:
                  admixed.write("{} {}\n".format(record.call_for_sample[sample].data.get('RO'), record.call_for_sample[sample].data.get('AO')[0]))
   else:
      locus = 0
      for record in reader:
         # I need to skip incomplete records, to limit the amount of loci
         if './.' in [call.data.get('GT') for call in record.calls]:
            continue
         P0RefCount = 0
         P0AltCount = 0
         P1RefCount = 0
         P1AltCount = 0
         for sample in SamplesInPop[args.pop0]:
            P0RefCount += record.call_for_sample[sample].data.get('GT').count('0')
            P0AltCount += record.call_for_sample[sample].data.get('GT').count('1')
         for sample in SamplesInPop[args.pop1]:
            P1RefCount += record.call_for_sample[sample].data.get('GT').count('0')
            P1AltCount += record.call_for_sample[sample].data.get('GT').count('1')
         admixed.write("l_{}:{}\n".format(record.CHROM, record.POS))
#        admixed.write("locus_{}\n".format(locus))
         for pop in args.admixed:
            # bgc.h requires that population names start with 'p', and loci name, with 'l'.
            admixed.write("p_{}\n".format(pop))
            # Here I assume that samples are iterated over always in the same order.
            for sample in SamplesInPop[pop]:
               AdRefCount = record.call_for_sample[sample].data.get('GT').count('0')
               AdAltCount = record.call_for_sample[sample].data.get('GT').count('1')
               if AdRefCount + AdAltCount != 2:
                  AdRefCount = -9
                  AdAltCount = -9
               admixed.write("{} {}\n".format(AdRefCount, AdAltCount))
         pop0.write("l_{}:{}\n{} {}\n".format(record.CHROM, record.POS, P0RefCount, P0AltCount))
         pop1.write("l_{}:{}\n{} {}\n".format(record.CHROM, record.POS, P1RefCount, P1AltCount))
#        pop0.write("locus_{}\n{} {}\n".format(locus, P0RefCount, P0AltCount))
#        pop1.write("locus_{}\n{} {}\n".format(locus, P1RefCount, P1AltCount))
         locus += 1
   pop0.close()
   pop1.close()
   admixed.close()

if __name__ == '__main__':
   main()
