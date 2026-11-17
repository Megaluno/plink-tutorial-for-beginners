#!/usr/bin/env bash
# Commands from the chapter: Haplotypes, proxy association and imputation
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle_qc --cow --chr 5 --from-mb 60 --to-mb 64 --make-bed --out brd_region   # 59 SNPs
plink107 --noweb --bfile brd_region --cow --hap-window 3 --hap-freq --out hfreq
plink107 --noweb --bfile brd_region --cow --hap-snps snp00705-snp00707 --hap-assoc --out hsnps
plink107 --noweb --bfile brd_region --cow --hap-window 3 --hap-logistic --covar pcs.txt --out hlogpc
# 1. is there any haplotype effect at all?
plink107 --noweb --bfile brd_region --cow --chap --hap-snps snp00705-snp00707 --covar pcs.txt --out chap1
# 2. is anything left once we allow for snp00706?
plink107 --noweb --bfile brd_region --cow --chap --hap-snps snp00705-snp00707 --control snp00706 --covar pcs.txt --out chap2
# 3. does snp00706 have an effect beyond its neighbours?
plink107 --noweb --bfile brd_region --cow --chap --hap-snps snp00705-snp00707 --independent-effect snp00706 --covar pcs.txt --out chap3
plink107 --noweb --bfile brd_region --cow --proxy-assoc snp00706 --proxy-verbose --out proxy1
plink107 --noweb --file masked --cow --proxy-impute all --proxy-replace --make-bed --out imputed
