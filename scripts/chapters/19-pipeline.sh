#!/usr/bin/env bash
# Commands from the chapter: A complete QC pipeline and the answer key
# Run from a folder that holds the input files named in the commands.

# 0. Convert and look (nothing is removed yet)
plink --file cattle --cow --make-bed --out cattle
plink --bfile cattle --cow --missing --out miss                                # Chapter 6
plink --bfile cattle --cow --freq --out freq                                   # Chapter 7
for b in ANG HOL NDA; do                                                       # Chapter 7: HWE within breed
  plink --bfile cattle --cow --keep k_$b.txt --hardy --out hwe_$b
done
plink --bfile cattle --cow --allow-no-sex --test-missing --out tmiss           # Chapter 6
plink --bfile cattle --cow --autosome --indep-pairwise 50 5 0.2 --out prune
plink --bfile cattle --cow --extract prune.prune.in --het --out het            # Chapter 8
plink --bfile cattle --cow --check-sex --out sexcheck                          # Chapter 8
plink --bfile cattle --cow --mendel --out mendel                               # Chapter 8
plink --bfile cattle --cow --extract prune.prune.in --genome --min 0.2 --out ibd   # Chapter 10
plink --bfile cattle --cow --extract prune.prune.in --maf 0.05 --pca 10 header tabs --out pca  # Chapter 11

# 1. Build the lists from those results (by reading the outputs and plots)
#    remove_animals.txt : 6 poor-DNA animals, contaminated ANG_030, duplicate HOL_050_rep,
#                         calf with the wrong sire HOL_F03_C2, crossbred HOL_077
#    sex_update.txt     : 4 animals set to their genetic sex
#    snp_exclude.txt    : 12 SNPs failing HWE (p < 1e-6) in any breed
#                         + 8 SNPs with --test-missing p < 1e-4
#    breed_fixed.clst   : ANG_066 relabelled as HOL

# 2. Apply them: animals first, then SNPs
plink --bfile cattle --cow --remove remove_animals.txt --update-sex sex_update.txt --make-bed --out step1
plink --bfile step1  --cow --exclude snp_exclude.txt --geno 0.05 --make-bed --out step2
plink --bfile step2  --cow --mind 0.05 --make-bed --out step3
plink --bfile step3  --cow --maf 0.01 --make-bed --out cattle_qc

# 3. Structure for the association models
plink --bfile cattle_qc --cow --autosome --indep-pairwise 50 5 0.2 --out qc_prune
plink --bfile cattle_qc --cow --extract qc_prune.prune.in --pca 10 header tabs --out qc_pca
