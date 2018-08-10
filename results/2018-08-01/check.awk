BEGIN{
   POP["Er29_122"]   = "europaeus";  POP["Er30_79"]   = "europaeus";  POP["Er31_453"]  = "europaeus";  POP["Er33_211"]   = "europaeus"
   POP["Er34_197"]   = "europaeus";  POP["Er51_436"]  = "europaeus";  POP["Er52_451"]  = "europaeus";  POP["Er54_AU1"]   = "europaeus"
   POP["Er56_AZ5"]   = "europaeus";  POP["Er57_COR4"] = "europaeus";  POP["Er63_IR6"]  = "europaeus";  POP["Er71_SAR2"]  = "europaeus"
   POP["Er72_SIE1B"] = "europaeus";  POP["Er74_SP16"] = "europaeus";  POP["Er26_JUG4"] = "roumanicus"; POP["Er27_SK32"]  = "roumanicus"
   POP["Er28_112"]   = "roumanicus"; POP["Er32_183"]  = "roumanicus"; POP["Er35_209"]  = "roumanicus"; POP["Er37_SK27"]  = "roumanicus"
   POP["Er39_PL1"]   = "roumanicus"; POP["Er40_M2"]   = "roumanicus"; POP["Er41_GR36"] = "roumanicus"; POP["Er42_GR35"]  = "roumanicus"
   POP["Er43_SL7"]   = "roumanicus"; POP["Er44_VOJ1"] = "roumanicus"; POP["Er45_BLG3"] = "roumanicus"; POP["Er46_RMN7"]  = "roumanicus"
   POP["Er47_CR4"]   = "roumanicus"; POP["Er48_BH16"] = "roumanicus"; POP["Er50_R3"]   = "roumanicus"; POP["Er55_AU7"]   = "roumanicus"
   POP["Er60_GR5"]   = "roumanicus"; POP["Er66_IT3"]  = "roumanicus"; POP["Er69_R2"]   = "roumanicus"; POP["Er70_RMN42"] = "roumanicus"
   POP["Er38_LB1"]   = "concolor";   POP["Er49_GR38"] = "concolor";   POP["Er53_ASR7"] = "concolor";   POP["Er64_IS1"]   = "concolor"
   POP["Er75_TRC2A"] = "concolor"
}(FILENAME == "LOD3_ordered.012"){
   for (pos = 1; pos <= NF - 1; pos ++) {
      if ($(pos+1) >= 0) {
          ALTFREQ[POP[$1]][pos]  += $(pos+1)
         WITHDATA[POP[$1]][pos] += 1
      } else {
          ALTFREQ[POP[$1]][pos]  += 0
         WITHDATA[POP[$1]][pos] += 0
      }
      if ($1 == "Er26_JUG4") GENO_JUG4[pos] = $(pos+1)
   }
}
(FILENAME == "LOD3_P.txt"){
   if (WITHDATA["concolor"][FNR] > 0) {
      ALTCON = ALTFREQ["concolor"][FNR] / (2 * WITHDATA["concolor"][FNR])
   } else {
      ALTCON = "NA"
   }
   if (WITHDATA["roumanicus"][FNR] > 0) {
      ALTROU = ALTFREQ["roumanicus"][FNR] / (2 * WITHDATA["roumanicus"][FNR])
   } else {
      ALTROU = "NA"
   }
   if (WITHDATA["europaeus"][FNR] > 0) {
      ALTEUR = ALTFREQ["europaeus"][FNR] / (2 * WITHDATA["europaeus"][FNR])
   } else {
      ALTEUR = "NA"
   }
   if (WITHDATA["concolor"][FNR] + WITHDATA["roumanicus"][FNR] + WITHDATA["europaeus"][FNR] > 0) {
      GLOBAL_ALT_FREQ = (ALTFREQ["concolor"][FNR] + ALTFREQ["roumanicus"][FNR] + ALTFREQ["europaeus"][FNR]) / (2 * (WITHDATA["concolor"][FNR] + WITHDATA["roumanicus"][FNR] + WITHDATA["europaeus"][FNR]))
   } else {
      GLOBAL_ALT_FREQ = "NA"
   }
   if ((ALTCON != "NA") && (ALTROU != "NA") && (ALTEUR != "NA")) {
      MEAN_ALT_FREQ = (ALTCON + ALTROU + ALTEUR) / 3.0
   } else {
      MEAN_ALT_FREQ = "NA"
   }
   if (GLOBAL_ALT_FREQ > 0.5) {
      MAJOR_GLOB_FREQ = GLOBAL_ALT_FREQ
      MAJOR_CON_FREQ = ALTCON
      MAJOR_ROU_FREQ = ALTROU
      MAJOR_EUR_FREQ = ALTEUR
      MAJOR_MEAN_FREQ = MEAN_ALT_FREQ
   } else {
      MAJOR_GLOB_FREQ = 1.0 - GLOBAL_ALT_FREQ
      MAJOR_CON_FREQ = 1.0 - ALTCON
      MAJOR_ROU_FREQ = 1.0 - ALTROU
      MAJOR_EUR_FREQ = 1.0 - ALTEUR
      MAJOR_MEAN_FREQ = 1.0 - MEAN_ALT_FREQ
      if (GENO_JUG4[FNR] >= 0) GENO_JUG4[FNR] = 2 - GENO_JUG4[FNR]
   }
   LINE = $0   sprintf("\t%.6f\t%.6f\t%.6f", MAJOR_EUR_FREQ, MAJOR_ROU_FREQ, MAJOR_CON_FREQ)
   LINE = LINE "\t" GENO_JUG4[FNR]
   LINE = LINE sprintf("\t%3f", ($3 - MAJOR_EUR_FREQ)**2 + ($4 - MAJOR_ROU_FREQ)**2 + ($5 - MAJOR_CON_FREQ)**2)
   LINE = LINE sprintf("\t%.4f\t%.4f", MAJOR_GLOB_FREQ, MAJOR_MEAN_FREQ)
   LINE = LINE sprintf("\t%i\t%i\t%i", 14 - WITHDATA["europaeus"][FNR], 22 - WITHDATA["roumanicus"][FNR], 5 - WITHDATA["concolor"][FNR])
   print LINE
}
