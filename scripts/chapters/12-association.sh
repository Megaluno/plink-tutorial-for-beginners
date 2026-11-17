#!/usr/bin/env bash
# Commands from the chapter: Case/control association
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle_qc --cow --allow-no-sex --assoc --adjust --out a_basic
plink --bfile cattle_qc --cow --allow-no-sex --logistic --covar pcs.txt --hide-covar --adjust --out lg_pc
plink --bfile cattle_qc --cow --allow-no-sex --assoc fisher --out a_fisher
plink --bfile cattle_qc --cow --allow-no-sex --model --out m_all
plink --bfile cattle_qc --cow --allow-no-sex --logistic --covar pcs.txt --genotypic --out lg_geno
plink --bfile cattle_qc --cow --allow-no-sex --mh --within breed.clst --out mh
plink --bfile cattle_qc --cow --allow-no-sex --bd --within breed.clst --out bd
plink --bfile cattle_qc --cow --allow-no-sex --homog --within breed.clst --out homog
# one breed, rare SNPs or a small study: exact test
plink --bfile cattle_qc --cow --allow-no-sex --assoc fisher --out a_fisher          # p = 1.4e-08
# one breed, common SNPs: allelic and genotype tests, and check lambda
plink --bfile cattle_qc --cow --allow-no-sex --model --out m_all                    # TREND p = 9.8e-09
plink --bfile cattle_qc --cow --allow-no-sex --assoc --adjust --out a_basic         # lambda = 10.19 here!
# several breeds: model the structure
plink --bfile cattle_qc --cow --allow-no-sex --logistic --covar pcs.txt --hide-covar --adjust --out lg_pc   # p = 2.3e-08, lambda 1.04
plink --bfile cattle_qc --cow --allow-no-sex --mh --within breed.clst --out mh                              # p = 9.9e-09
# then: same effect in every breed? additive?
plink --bfile cattle_qc --cow --allow-no-sex --bd --within breed.clst --out bd                              # Breslow-Day p = 0.59
plink --bfile cattle_qc --cow --allow-no-sex --logistic --covar pcs.txt --genotypic --out lg_geno          # DOMDEV p = 0.92
plink2 --bfile cattle_qc --cow --glm hide-covar firth-fallback --covar pcs.txt --adjust --out glm_brd
