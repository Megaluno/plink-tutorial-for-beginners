#!/usr/bin/env bash
# Commands from the chapter: Missing genotypes: per SNP and per animal
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle --cow --geno 0.10 --make-bed --out qc1   # round 1, SNPs
plink --bfile qc1    --cow --mind 0.10 --make-bed --out qc2   # round 1, animals
plink --bfile qc2    --cow --mind 0.03 --make-bed --out qc3   # round 2, animals
plink --bfile qc3    --cow --geno 0.01 --make-bed --out qc4   # round 2, SNPs
plink --bfile cattle --cow --missing --out miss
plink --bfile cattle --cow --geno 0.10 --make-bed --out qc1
plink --bfile qc1    --cow --mind 0.10 --make-bed --out qc2
plink --bfile qc2    --cow --mind 0.03 --make-bed --out qc3
plink --bfile qc3    --cow --geno 0.01 --make-bed --out qc4
plink --bfile cattle --cow --allow-no-sex --test-missing --out tmiss
