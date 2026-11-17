#!/usr/bin/env bash
# Commands from the chapter: Sample checks: sex, pedigree and heterozygosity
# Run from a folder that holds the input files named in the commands.

plink --bfile cattle --cow --autosome --indep-pairwise 50 5 0.2 --out prune
plink --bfile cattle --cow --extract prune.prune.in --het --out het
plink --bfile cattle --cow --autosome --ibc --out ibc
grep -E "IID|ANG_001|ANG_030|NDA_040" ibc.ibc
plink --bfile cattle --cow --check-sex --out sexcheck
plink --bfile cattle --cow --check-sex 0.3 0.7 --out sex_37        # female if F < 0.3, male if F > 0.7
plink --bfile cattle --cow --check-sex ycount --out sex_y          # also count Y-chromosome calls
plink --bfile cattle --cow --impute-sex --make-bed --out cattle_sexed
plink --bfile cattle --cow --mendel --out mendel
# remove families with > 5% errors and SNPs with > 10% errors
plink --bfile cattle --cow --me 0.05 0.1 --make-bed --out me_filtered
# or keep everything but set the offending genotypes to missing
plink --bfile cattle --cow --set-me-missing --make-bed --out me_missing
