#!/usr/bin/env bash
# Commands from the chapter: Real data: cattle and sheep
# Run from a folder that holds the input files named in the commands.

plink --file cattle_real --cow --make-bed --out cattle_real
plink --file Kijas2016_Sheep_subset --sheep --make-bed --out sheep
plink --bfile cattle_real --cow --missing --out c_miss
plink --bfile cattle_real --cow --geno 0.05 --mind 0.05 --maf 0.01 --hwe 1e-6 --make-bed --out cattle_real_qc
plink --bfile sheep --sheep --geno 0.05 --mind 0.05 --maf 0.01 --make-bed --out sheep_qc
plink --bfile cattle_real_qc --cow --indep-pairwise 50 5 0.2 --out c_prune     # 2,984 of 7,211 SNPs kept
plink --bfile cattle_real_qc --cow --extract c_prune.prune.in --het --out c_het
plink --bfile cattle_real_qc --cow --extract c_prune.prune.in --pca 10 header tabs --out c_pca
plink --bfile cattle_real_qc --cow --homozyg --homozyg-density 1000 --homozyg-gap 3000 \
      --homozyg-snp 30 --homozyg-kb 2000 --homozyg-window-snp 30 --out c_roh
plink --bfile cattle_real_qc --cow --extract c_prune.prune.in --genome --out c_ibd
plink --bfile cattle_real_qc --cow --make-rel square --out c_grm
plink --bfile sheep_qc --sheep --indep-pairwise 50 5 0.2 --out s_prune         # 843 of 4,686 SNPs kept
plink --bfile sheep_qc --sheep --extract s_prune.prune.in --pca 10 header tabs --out s_pca
plink --bfile sheep_qc --sheep --fst --within sheep_breed.clst --out s_fst
plink --bfile cattle_real_qc --cow --r2 --ld-window-kb 5000 --ld-window 99999 --ld-window-r2 0 --out c_ld
