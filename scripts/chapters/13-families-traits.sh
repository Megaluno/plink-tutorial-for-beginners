#!/usr/bin/env bash
# Commands from the chapter: Families, quantitative traits, interaction and permutation
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle_qc --cow --allow-no-sex --pheno pheno.txt --pheno-name WWT --assoc --adjust --out q_adj
plink --bfile cattle_qc --cow --allow-no-sex --pheno pheno.txt --pheno-name WWT \
      --linear --covar cov_all.txt --hide-covar --adjust --out lin
plink --bfile cattle_qc --cow --allow-no-sex --pheno pheno.txt --pheno-name WWT \
      --linear --covar cov_all.txt --interaction --parameters 1-6,10 --out lin_int
plink --bfile cattle_qc --cow --allow-no-sex --pheno pheno.txt --pheno-name WWT --gxe --covar env.txt --out gxe
plink ... --linear --covar cov_all.txt --hide-covar --chr 14 --from-mb 23 --to-mb 27 \
          --condition snp02145 --out lin_chr14_cond
plink2 --bfile cattle_qc --cow --autosome --maf 0.05 --make-grm-bin --out grm   # grm.grm.bin, .grm.N.bin, .grm.id
gcta64 --mlma-loco --bfile cattle_qc --autosome-num 29 --pheno wwt.phen \
       --qcovar qcovar.txt --covar covar.txt --out wwt_mlma          # GCTA
gemma -bfile cattle_qc -gk 1 -o grm                                  # GEMMA: 1. relationship matrix
gemma -bfile cattle_qc -k output/grm.cXX.txt -lmm 4 -c covar_gemma.txt -o wwt   # 2. association
plink --bfile cattle_qc --cow --autosome --recode A --out geno
python3 scripts/lmm_gwas.py geno.raw cattle_qc.bim pheno.txt WWT cov_all.txt SEX_M,FEEDLOT,AGE lmm_noPC.txt
python3 scripts/lmm_gwas.py geno.raw cattle_qc.bim pheno.txt WWT cov_all.txt PC1,PC2,SEX_M,FEEDLOT,AGE lmm_PC.txt
plink --bfile cattle_qc --cow --allow-no-sex --tdt --out tdt
plink ... --logistic --covar pcs.txt --hide-covar --mperm 1000 --seed 1 --out lg_mperm
plink --bfile cattle_qc --cow --allow-no-sex --assoc --mperm 1000 --seed 1 --within breed.clst --out mperm_within
plink --fam cattle.fam --cow --allow-no-sex --dosage chr14.dosage format=1 --map chr14.dosage.map \
      --pheno pheno.txt --pheno-name WWT --covar cov_all.txt --out dos
