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
plink2 --bfile cattle_qc --cow --glm hide-covar firth-fallback --covar pcs.txt --adjust --out glm_brd
