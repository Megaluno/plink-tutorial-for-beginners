#!/usr/bin/env bash
# Commands from the chapter: Allele frequency and Hardy-Weinberg equilibrium
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle --cow --freq --out freq
plink --bfile cattle --cow --hardy --out hwe
for b in ANG HOL NDA; do
  awk -v b=$b '$3==b {print $1, $2}' breed.clst > keep_$b.txt
  plink --bfile cattle --cow --keep keep_$b.txt --hardy --out hwe_$b
done
awk '$3=="ALL" && $9 < 1e-6 {print $2}' hwe_ANG.hwe hwe_HOL.hwe hwe_NDA.hwe | sort -u > hwe_fail.txt
plink --bfile cattle --cow --exclude hwe_fail.txt --make-bed --out cattle_hwe
