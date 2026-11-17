#!/usr/bin/env bash
# Commands from the chapter: Beyond single SNPs: gene sets, epistasis, environments and meta-analysis
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle_qc --cow --make-set genes.txt --write-set --out sets
# --make-set: 8 sets defined.
plink --bfile cattle_qc --cow --allow-no-sex --pheno pheno.txt --pheno-name WWT \
      --make-set genes.txt --linear mperm=10000 set-test --covar cov_all.txt --hide-covar --seed 2 --out settest
plink --bfile cattle_qc --cow --allow-no-sex --autosome --pheno pke.txt --fast-epistasis --epi1 1e-5 --out pke_fe
plink --bfile cattle_qc --cow --allow-no-sex --autosome --pheno pke.txt --fast-epistasis boost --epi1 1e-5 --out pke_feb
plink --bfile cattle_qc --cow --allow-no-sex --autosome --pheno pke.txt --epistasis --epi1 1e-5 --out pke_ep
plink --bfile cattle_qc --cow --allow-no-sex --pheno pke.txt --twolocus snp01294 snp03388 --out pke_two
# AGE_C = (AGE - 205) / 10
plink ... --linear --covar cov_c.txt --covar-name PC1,PC2,SEX_M,FEEDLOT,AGE_C \
          --interaction --parameters 1-6,11 --out int_agec
for b in ANG HOL NDA; do
  plink --bfile cattle_qc --cow --allow-no-sex --pheno pheno.txt --pheno-name WWT --keep-allele-order \
        --keep k_$b.txt --linear --covar cov_all.txt --covar-name SEX_M,FEEDLOT,AGE --hide-covar --ci 0.95 --out mk_$b
done
plink --cow --meta-analysis mk_ANG.assoc.linear mk_HOL.assoc.linear mk_NDA.assoc.linear + qt --out meta_fixed
