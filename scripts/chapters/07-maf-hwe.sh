#!/usr/bin/env bash
# Commands from the chapter: Allele frequency and Hardy-Weinberg equilibrium
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle --cow --freq --out freq
head -3 freq.frq
plink --bfile cattle --cow --freq counts --out freqc
head -3 freqc.frq.counts
plink --bfile cattle --cow --freqx --out freqx
head -3 freqx.frqx
plink --bfile cattle --cow --freq --within breed.clst --out freq_breed
grep -E "snp00001|snp00706" freq_breed.frq.strat
plink --bfile cattle --cow --hardy --out hwe
for b in ANG HOL NDA; do
  awk -v b=$b '$3==b {print $1, $2}' breed.clst > keep_$b.txt
  plink --bfile cattle --cow --keep keep_$b.txt --hardy --out hwe_$b
done
awk '$3=="ALL" && $9 < 1e-6 {print $2}' hwe_ANG.hwe hwe_HOL.hwe hwe_NDA.hwe | sort -u > hwe_fail.txt
plink --bfile cattle --cow --exclude hwe_fail.txt --make-bed --out cattle_hwe
plink --bfile cattle --cow --hwe 1e-6 --make-bed --out h1                     # controls (and founders) only
plink --bfile cattle --cow --hwe 1e-6 midp --make-bed --out h2                # mid-p exact test
plink --bfile cattle --cow --hwe 1e-6 include-nonctrl --make-bed --out h3     # cases too
plink2 --bfile cattle --cow --keep-if "PHENO1 == 1" --hardy --out hwe_ctrl   # 217 controls
awk 'NR > 1 && $NF < 1e-6 {print $2}' hwe_ctrl.hardy > hwe_fail.txt
plink2 --bfile cattle --cow --exclude hwe_fail.txt --make-pgen --out cattle_hwe
