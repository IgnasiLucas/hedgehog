#!/bin/bash

ipyrad -p params-hhog.txt -r | gawk 'BEGIN{
   POP["Er73_SNG1"] = "Atelerix";  POP["Er38_LB1"] = "concolor";     POP["Er49_GR38"] = "concolor";   POP["Er53_ASR7"] = "concolor";  POP["Er64_IS1"] = "concolor"
   POP["Er75_TRC2A"] = "concolor"; POP["Er29_122"] = "europeaus";    POP["Er30_79"] = "europeaus";    POP["Er31_453"] = "europeaus";  POP["Er33_211"] = "europeaus"
   POP["Er34_197"] = "europeaus";  POP["Er36_SK24"] = "europeaus";   POP["Er51_436"] = "europeaus";   POP["Er52_451"] = "europeaus";  POP["Er54_AU1"] = "europeaus"
   POP["Er55_AU7"] = "europeaus";  POP["Er56_AZ5"] = "europeaus";    POP["Er57_COR4"] = "europeaus";  POP["Er58_FI7"] = "europeaus";  POP["Er59_FR1"] = "europeaus"
   POP["Er63_IR6"] = "europeaus";  POP["Er67_IT5"] = "europeaus";    POP["Er68_PRT1B"] = "europeaus"; POP["Er71_SAR2"] = "europeaus"; POP["Er72_SIE1B"] = "europeaus"
   POP["Er74_SP16"] = "europeaus"; POP["Er65_IS25"] = "Hemiechinus"; POP["Er37_SK27"] = "hybrid";     POP["Er26_JUG4"] = "romanicus"; POP["Er27_SK32"] = "romanicus"
   POP["Er28_112"] = "romanicus";  POP["Er32_183"] = "romanicus";    POP["Er35_209"] = "romanicus";   POP["Er39_PL1"] = "romanicus";  POP["Er40_M2"] = "romanicus"
   POP["Er41_GR36"] = "romanicus"; POP["Er42_GR35"] = "romanicus";   POP["Er43_SL7"] = "romanicus";   POP["Er44_VOJ1"] = "romanicus"; POP["Er45_BLG3"] = "romanicus"
   POP["Er46_RMN7"] = "romanicus"; POP["Er47_CR4"] = "romanicus";    POP["Er48_BH16"] = "romanicus";  POP["Er50_R3"] = "romanicus";   POP["Er60_GR5"] = "romanicus"
   POP["Er61_GR87"] = "romanicus"; POP["Er62_GR95"] = "romanicus";   POP["Er66_IT3"] = "romanicus";   POP["Er69_R2"] = "romanicus";   POP["Er70_RMN42"] = "romanicus"

   ORDER[1] = "Er29_122";   ORDER[2] = "Er30_79";     ORDER[3] = "Er31_453";    ORDER[4] = "Er33_211";    ORDER[5] = "Er34_197"
   ORDER[6] = "Er36_SK24";  ORDER[7] = "Er51_436";    ORDER[8] = "Er52_451";    ORDER[9] = "Er54_AU1";    ORDER[10] = "Er55_AU7"
   ORDER[11] = "Er56_AZ5";  ORDER[12] = "Er57_COR4";  ORDER[13] = "Er58_FI7";   ORDER[14] = "Er59_FR1";   ORDER[15] = "Er63_IR6"
   ORDER[16] = "Er67_IT5";  ORDER[17] = "Er68_PRT1B"; ORDER[18] = "Er71_SAR2";  ORDER[19] = "Er72_SIE1B"; ORDER[20] = "Er74_SP16"
   ORDER[21] = "Er37_SK27"; ORDER[22] = "Er26_JUG4";  ORDER[23] = "Er27_SK32";  ORDER[24] = "Er28_112";   ORDER[25] = "Er32_183"
   ORDER[26] = "Er35_209";  ORDER[27] = "Er39_PL1";   ORDER[28] = "Er40_M2";    ORDER[29] = "Er41_GR36";  ORDER[30] = "Er42_GR35"
   ORDER[31] = "Er43_SL7";  ORDER[32] = "Er44_VOJ1";  ORDER[33] = "Er45_BLG3";  ORDER[34] = "Er46_RMN7";  ORDER[35] = "Er47_CR4"
   ORDER[36] = "Er48_BH16"; ORDER[37] = "Er50_R3";    ORDER[38] = "Er60_GR5";   ORDER[39] = "Er61_GR87";  ORDER[40] = "Er62_GR95"
   ORDER[41] = "Er66_IT3";  ORDER[42] = "Er69_R2";    ORDER[43] = "Er70_RMN42"; ORDER[44] = "Er38_LB1";   ORDER[45] = "Er49_GR38"
   ORDER[46] = "Er53_ASR7"; ORDER[47] = "Er64_IS1";   ORDER[48] = "Er75_TRC2A"; ORDER[49] = "Er73_SNG1";  ORDER[50] = "Er65_IS25"

   OFFSET = -1
   HOWMUCH = 0
}(/^      /){
   OFFSET += HOWMUCH
   HOWMUCH = NF - 1
   for (i = 1; i < NF; i++) VARNAME[OFFSET + i] = $i
}(/^Er/){
   SAMPLE[$1] = 1
   for (i = 1; i < NF; i++) M[$1, OFFSET + i] = $(i + 1)
}END{
   HEADER = "#\tPopulation"
   for (i = 1; i <= OFFSET + HOWMUCH; i++) {
      HEADER = HEADER "\t" VARNAME[i]
   }
   print HEADER
   for (i = 1; i <= 50; i++) {
      STRING = ORDER[i] "\t" POP[ORDER[i]]
      for (j = 1; j <= OFFSET + HOWMUCH; j++) {
         STRING = STRING "\t" M[ORDER[i], j]
      }
      print STRING
   }
}' > z.summary.txt

gnuplot < mapping_success.gnp
