#((/^#/) || (($6 >= 40) && ($8 ~ /TYPE=snp;|TYPE=complex;|TYPE=ins;|TYPE=del;|TYPE=mnp;|TYPE=complex,snp;|TYPE=snp,complex;|TYPE=snp,snp;|TYPE=mnp,snp;|TYPE=snp,mnp;|TYPE=del,snp;|TYPE=complex,complex;|TYPE=ins,ins;|TYPE=ins,snp;|TYPE=snp,del;|TYPE=del,del;|TYPE=snp,ins;|TYPE=del,ins;|TYPE=del,complex;|TYPE=complex,del;|TYPE=complex,mnp;|TYPE=ins,del;|TYPE=ins,complex;|TYPE=complex,ins;|TYPE=mnp,mnp;/))){
((/^#/) || (($6 >= 40) && ($8 ~ /TYPE=snp;|TYPE=complex;|TYPE=mnp;|TYPE=snp,complex;|TYPE=snp,snp;|TYPE=snp,mnp;|TYPE=complex,snp;|TYPE=mnp,snp;|TYPE=complex,complex;|TYPE=snp,complex,snp;|TYPE=complex,snp,snp;|TYPE=snp,complex,complex;|TYPE=snp,snp,complex;|TYPE=mnp,mnp;|TYPE=complex,snp,complex;|TYPE=complex,mnp;|TYPE=snp,mnp,snp;|TYPE=mnp,complex;|TYPE=complex,complex,snp;|TYPE=snp,complex,mnp;|TYPE=mnp,snp,snp;|TYPE=snp,snp,mnp;|TYPE=snp,mnp,complex;|TYPE=snp,mnp,mnp;|TYPE=snp,snp,snp;|TYPE=mnp,snp,complex;|TYPE=mnp,snp,mnp;|TYPE=complex,complex,complex;|TYPE=complex,snp,mnp;|TYPE=mnp,complex,snp;|TYPE=complex,mnp,snp;|TYPE=mnp,mnp,snp;/))){
   print
}

# snp                2332548     0.924759
# complex              58484     0.947945
# mnp                  48885     0.967326
# snp,complex          24729     0.97713
# snp,snp              15998     0.983472
# snp,mnp              13579     0.988856
# complex,snp          13535     0.994222
# mnp,snp               7280     0.997108
# complex,complex       1329     0.997635
# snp,complex,snp        918     0.997999
# complex,snp,snp        535     0.998211
# snp,complex,complex    509     0.998413
# snp,snp,complex        399     0.998571
# mnp,mnp                373     0.998719
# complex,snp,complex    343     0.998855
# complex,mnp            340     0.99899
# snp,mnp,snp            305     0.999111
# mnp,complex            284     0.999223
# complex,complex,snp    250     0.999322
# snp,complex,mnp        224     0.999411
# mnp,snp,snp            161     0.999475
# snp,snp,mnp            143     0.999532
# snp,mnp,complex        135     0.999585
# snp,mnp,mnp            110     0.999629
# snp,snp,snp             96     0.999667
# mnp,snp,complex         95     0.999705
# mnp,snp,mnp             70     0.999732
# complex,complex,complex 62     0.999757
# complex,snp,mnp         61     0.999781
# mnp,complex,snp         59     0.999805
# complex,mnp,snp         56     0.999827
# mnp,mnp,snp             50     0.999847
