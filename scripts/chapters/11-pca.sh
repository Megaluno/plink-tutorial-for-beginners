#!/usr/bin/env bash
# Commands from the chapter: Population structure and principal component analysis
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle --cow --maf 0.05 --autosome --indep-pairwise 50 5 0.2 --out prune
plink --bfile cattle --cow --extract prune.prune.in --pca 10 header tabs var-wts --out pca
# 1. PCs from the 301 founders only, saving the SNP loadings
plink2 --bfile cattle --cow --extract prune.prune.in --keep founders.txt \
       --maf 0.05 --freq --pca 10 allele-wts --out pcaF
# 2. project all 361 animals onto those PCs
plink2 --bfile cattle --cow --extract prune.prune.in --read-freq pcaF.afreq \
       --score pcaF.eigenvec.allele 2 6 header-read no-mean-imputation variance-standardize \
       --score-col-nums 7-16 --out proj
plink --bfile cattle --cow --extract prune.prune.in --cluster --mds-plot 4 --out mds
